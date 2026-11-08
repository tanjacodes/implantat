"""
Implant shape optimisation: one-factor-at-a-time sensitivity sweeps around the
baseline 3-D strut plate, followed by evaluation of a combined lightweight design.

Design load: left (fracture-side) molar clench x LOAD_FACTOR (= ~240 N bite force).
Constraints at the design load:
  plate von Mises (max)          <= 240 MPa  (fatigue limit of CP-Ti grade 4, safety ~2 on yield)
  interfragmentary motion (max)  <= 150 um   (secondary bone healing limit)
  peri-screw bone strain (p95)   <= 3000 ue  (overload / resorption threshold)
Objective: minimise implant volume (profile & bulk), maximise stability margin.
"""
import json, sys, os, time, gc
import numpy as np
from geometry import build
from fe import Model, CASES

LOAD_FACTOR = 2.0
MESH = dict(h_implant=0.55, h_near=1.3, h_far=3.8)
LIM = dict(vm=240.0, ifm=150.0, e1=3000.0)


def evaluate(name, P, cases=("left_molar",), mesh=MESH, keep=False):
    p = {**P, **mesh}
    msh = f"opt_{name}.msh"
    meta = build(p, msh)
    m = Model(msh, p, verbose=False); m.assemble()
    res = dict(name=name, params=P, volume_mm3=meta["implant_volume_mm3"], cases={})
    for c in cases:
        u, o = m.solve_static(c); m.free_solver()
        r, F, I = m.summary(u, o)
        res["cases"][c] = r
    lm = res["cases"]["left_molar"]
    res["design"] = dict(vm=lm["plate_vm_max"] * LOAD_FACTOR, ifm=lm["ifm_max_um"] * LOAD_FACTOR,
                         e1=lm["peri_screw_e1_p95_ue"] * LOAD_FACTOR)
    res["feasible"] = all(res["design"][k] <= LIM[k] for k in LIM)
    if not keep: os.remove(msh)
    del m; gc.collect()
    d = res["design"]
    print(f"{name:22s} V={res['volume_mm3']:6.1f} mm3 | vM {d['vm']:6.1f} | IFM {d['ifm']:6.1f} um | "
          f"e1p95 {d['e1']:6.0f} ue | {'OK' if res['feasible'] else 'FAIL'}", flush=True)
    return res


if __name__ == "__main__":
    sweeps = [("baseline", {})]
    sweeps += [(f"t={t}", dict(tp=t)) for t in (0.6, 0.8, 1.3)]
    sweeps += [(f"bar={w}", dict(wb=w)) for w in (1.2, 2.0, 2.4)]
    sweeps += [("struts=outer", dict(struts="outer")), ("struts=none", dict(struts="none")),
               ("struts=all+X", dict(diagonals=True))]
    sweeps += [("4-hole square", dict(holes_p=[-5.5, 5.5]))]
    sweeps += [(f"height z={z}", dict(center_z=z)) for z in (12.5, 15.5)]
    sweeps += [("bars 13 mm apart", dict(Q=6.5))]
    out = []
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    for name, P in sweeps:
        if only and name not in only: continue
        try:
            out.append(evaluate(name, P))
        except Exception as e:
            print(name, "failed:", e, flush=True)
        json.dump(out, open("results/optim_sweeps.json", "w"), indent=1)
