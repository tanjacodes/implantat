import json, sys
from geometry import build
from fe import Model
levels = {"coarse": dict(h_implant=0.7, h_near=1.6, h_far=4.5),
          "medium": dict(h_implant=0.55, h_near=1.3, h_far=3.8),
          "base":   dict(),
          }
out = {}
for name, p in levels.items():
    meta = build(p, f"conv_{name}.msh")
    m = Model(f"conv_{name}.msh", p, verbose=False); m.assemble()
    u, o = m.solve_static("left_molar"); m.free_solver()
    r, F, I = m.summary(u, o)
    r.update(n_tet=int(m.mesh.nelements), ndof=int(m.basis.N))
    out[name] = r; json.dump(r, open(f"results/conv_{name}.json", "w"), indent=1)
    print(name, r["ndof"], round(r["plate_vm_max"],1), round(r["plate_vm_p99"],1), round(r["ifm_max_um"],2), round(r["peri_screw_e1_max_ue"]), round(r["bite_force_N"],1), flush=True)
    del m
json.dump(out, open("results/convergence.json", "w"), indent=1)
