"""Static analysis of the baseline 3-D strut plate: three physiological clenching tasks."""
import json, sys, os
import numpy as np
from fe import Model, CASES, MUSCLES, MAT
from vis import grid

msh = sys.argv[1] if len(sys.argv) > 1 else "base.msh"
tag = sys.argv[2] if len(sys.argv) > 2 else "base"
os.makedirs("results", exist_ok=True); os.makedirs("vtu", exist_ok=True)
m = Model(msh)
m.assemble()
allres = {}
U = {}
for case in CASES:
    u, out = m.solve_static(case)
    m.free_solver()
    res, F, I = m.summary(u, out)
    allres[case] = res
    U[case] = u
    np.savez_compressed(f"results/ifm_{tag}_{case}.npz", pts=I["_pts"], mag=I["_mag"])
    if tag == "base":
        g = grid(m, u, cell={"von Mises (MPa)": F["vm_max"], "eps1 (ue)": F["e1_max"] * 1e6,
                             "eps3 (ue)": F["e3_min"] * 1e6})
        g.save(f"vtu/{tag}_{case}.vtu")
np.savez_compressed(f"results/u_{tag}.npz", **U)
json.dump(allres, open(f"results/static_{tag}.json", "w"), indent=1)
for c, r in allres.items():
    print(f"{c:12s} bite {r['bite_force_N']:6.1f} N | plate vM max {r['plate_vm_max']:6.1f} p99 {r['plate_vm_p99']:6.1f} | "
          f"screw {r['screw_vm_max']:5.1f} | IFM max {r['ifm_max_um']:5.1f} um | peri-screw e1 {r['peri_screw_e1_max_ue']:6.0f} ue | umax {r['umax_mm']*1e3:5.1f} um")
