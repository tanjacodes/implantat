"""Evaluate the optimised plate on the base mesh for all load cases; save fields for figures."""
import json, numpy as np
from geometry import build
from fe import Model, CASES
from vis import grid
OPT = dict(holes_p=[-5.5, 5.5], wb=2.4)
meta = build(OPT, "opt.msh")
m = Model("opt.msh", OPT); m.assemble()
res, U = {}, {}
for c in CASES:
    u, o = m.solve_static(c); m.free_solver()
    r, F, I = m.summary(u, o); res[c] = r; U[c] = u
    np.savez_compressed(f"results/ifm_opt_{c}.npz", pts=I["_pts"], mag=I["_mag"])
    grid(m, u, cell={"von Mises (MPa)": F["vm_max"], "eps1 (ue)": F["e1_max"]*1e6, "eps3 (ue)": F["e3_min"]*1e6}).save(f"vtu/opt_{c}.vtu")
res["volume_mm3"] = meta["implant_volume_mm3"]; res["params"] = OPT
np.savez_compressed("results/u_opt.npz", **U)
json.dump(res, open("results/static_opt.json", "w"), indent=1)
for c in CASES:
    r = res[c]; print(c, round(r["bite_force_N"],1), round(r["plate_vm_max"],1), round(r["plate_vm_p99"],1), round(r["screw_vm_max"],1), round(r["ifm_max_um"],2), round(r["peri_screw_e1_max_ue"]), round(r["peri_screw_e1_p95_ue"]))
print("V", meta["implant_volume_mm3"])
