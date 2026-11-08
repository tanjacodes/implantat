"""
Dynamic analyses of the baseline model:
  1. modal analysis (eigenfrequencies / mode shapes) for two support conditions
  2. transient chewing cycle (left molar power stroke, 0.8 s, Newmark-beta)
  3. sudden hard bite (1 ms rise, step load) -> dynamic amplification factor
Rayleigh damping, zeta = 5 % at the 1st and 4th eigenfrequency of the occluded system.
"""
import json, time, gc
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, LinearOperator
from fe import Model, CASES
from skfem import Basis

MESH_P = dict(h_implant=0.55, h_near=1.3, h_far=3.8)   # "medium" mesh (static results within ~4 % of base)
m = Model("conv_medium.msh", MESH_P)
m.assemble(mass=True)
K, M = m.K, m.M
N = m.basis.N
case = "left_molar"
f0 = m.load_vector(case)
Dcase = m.dirichlet(case)
Dcond = m.D_cond.copy(); Dbite = m.D_bite.copy()
out = {}

# ------------------------------------------------------------------ 1. modal
def modes(D, k=8, label=""):
    I = m.basis.complement_dofs(D)
    t0 = time.time()
    m.factor_solve(I, np.zeros(len(I)))
    n = len(I)
    op = LinearOperator((n, n), matvec=lambda x: m.resolve(np.asarray(x).ravel()), dtype=float)
    Mii = M[I][:, I].tocsr()
    lam, V = eigsh(m._Aii, k=k, M=Mii, sigma=0.0, OPinv=op, which="LM")  # Aii (triu) only used via OPinv
    del Mii
    m.free_solver(); gc.collect()
    fr = np.sqrt(np.abs(lam)) / (2 * np.pi)
    Vf = np.zeros((N, k)); Vf[I] = V
    print(f"modes {label}: {np.round(fr, 0)} Hz  ({time.time()-t0:.0f}s)", flush=True)
    return fr, Vf

import sys
IMPACT_ONLY = "impact" in sys.argv
if IMPACT_ONLY:
    out = json.load(open("results/dynamic.json"))
    fr_occ = np.array(out["modal"]["occluded_Hz"])
else:
    fr_free, V_free = modes(Dcond, 6, "condyles fixed, teeth apart")
    fr_occ, V_occ = modes(Dcase, 6, "condyles fixed + left molar contact")
    np.savez_compressed("results/modes.npz", fr_free=fr_free, V_free=V_free, fr_occ=fr_occ, V_occ=V_occ)
    out["modal"] = dict(free_Hz=fr_free.tolist(), occluded_Hz=fr_occ.tolist())

# Rayleigh damping
zeta = 0.05
w1, w2 = 2 * np.pi * fr_occ[0], 2 * np.pi * fr_occ[3]
alpha = 2 * zeta * w1 * w2 / (w1 + w2); beta = 2 * zeta / (w1 + w2)
out["rayleigh"] = dict(alpha=alpha, beta=beta, zeta=zeta)
print("Rayleigh", alpha, beta)

# monitors
plate_el = np.where((m.tag == 2) | (m.tag == 3))[0]
pb = m.basis.with_elements(plate_el)
E, nu = 105000.0, 0.34
lamT, muT = E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))

def plate_vm(u):
    G = pb.interpolate(u).grad
    e = 0.5 * (G + np.swapaxes(G, 0, 1)); tr = np.trace(e)
    s = 2 * muT * e + lamT * tr * np.eye(3)[:, :, None, None]
    vm = np.sqrt(0.5 * ((s[0, 0] - s[1, 1]) ** 2 + (s[1, 1] - s[2, 2]) ** 2 + (s[2, 2] - s[0, 0]) ** 2)
                 + 3 * (s[0, 1] ** 2 + s[1, 2] ** 2 + s[2, 0] ** 2))
    v = vm.max(1)
    return float(v.max()), float(np.percentile(v, 99)), v

# gap node pairs for IFM
from geometry import fracture_geom
from scipy.spatial import cKDTree
A, B, d, n = fracture_geom(m.P)
Cf = 0.5 * (A + B)
gv = np.unique(m.mesh.t[:, m.tag == 4])
L = m.to_local(m.mesh.p[:, gv], "L")
sd = (L[0] - Cf[0]) * n[0] + (L[2] - Cf[1]) * n[1]
fa = gv[sd < -m.P["gap"] / 2 + 0.05]; fb = gv[sd > m.P["gap"] / 2 - 0.05]
dd, j = cKDTree(m.mesh.p[:, fb].T).query(m.mesh.p[:, fa].T)
pa, pbn = fa[dd < 1.5], fb[j[dd < 1.5]]
ND = m.basis.nodal_dofs

def ifm(u):
    return float(np.linalg.norm(u[ND][:, pbn] - u[ND][:, pa], axis=0).max() * 1e3)

Kb = K[Dbite]; Mb = M[Dbite]

# static reference
I = m.basis.complement_dofs(Dcase)
us = np.zeros(N); us[I] = m.factor_solve(I, f0[I]); m.free_solver(); gc.collect()
ref = dict(plate_vm_max=plate_vm(us)[0], ifm_um=ifm(us), bite=float(-(K[Dbite] @ us - f0[Dbite]).sum()))
out["static_reference"] = ref
print("static ref", ref, flush=True)


def newmark(g, dt, nsteps, save_every, label):
    """average-acceleration Newmark, f(t) = g(t) f0, zero initial conditions"""
    gam, bet = 0.5, 0.25
    a0 = 1 / (bet * dt ** 2); a1 = gam / (bet * dt); a2 = 1 / (bet * dt); a3 = 1 / (2 * bet) - 1
    a4 = gam / bet - 1; a5 = dt / 2 * (gam / bet - 2)
    C = alpha * M + beta * K
    Keff = (K + a0 * M + a1 * C).tocsr()
    t0 = time.time()
    u = np.zeros(N); v = np.zeros(N); a = np.zeros(N)
    m.factor_solve(I, np.zeros(len(I)), A=Keff)
    hist = dict(t=[], g=[], plate_vm=[], plate_p99=[], ifm=[], bite=[])
    frames = []
    for k in range(1, nsteps + 1):
        t = k * dt
        rhs = g(t) * f0 + M @ (a0 * u + a2 * v + a3 * a) + C @ (a1 * u + a4 * v + a5 * a)
        un = np.zeros(N); un[I] = m.resolve(rhs[I])
        an = a0 * (un - u) - a2 * v - a3 * a
        vn = v + dt * ((1 - gam) * a + gam * an)
        u, v, a = un, vn, an
        if k % save_every == 0:
            vmx, p99, _ = plate_vm(u)
            r = Kb @ u + Mb @ a + (alpha * Mb + beta * Kb) @ v - g(t) * f0[Dbite]
            hist["t"].append(t); hist["g"].append(g(t)); hist["plate_vm"].append(vmx)
            hist["plate_p99"].append(p99); hist["ifm"].append(ifm(u)); hist["bite"].append(float(-r.sum()))
            frames.append(u.astype(np.float32))
    m.free_solver(); gc.collect()
    print(f"{label}: {nsteps} steps in {time.time()-t0:.0f}s", flush=True)
    return hist, np.array(frames)

# ------------------------------------------------------------------ 2. chewing cycle
T = 0.8
def g_chew(t):
    tt = t % T
    if tt < 0.30:      # jaw closing, light muscle activity (10 %)
        return 0.10 * np.sin(np.pi * tt / 0.30) ** 2
    if tt < 0.55:      # power stroke on the bolus: activation rises to 100 % and back
        return np.sin(np.pi * (tt - 0.30) / 0.25) ** 2
    return 0.0         # opening
dt = 2e-3
if not IMPACT_ONLY:
    h_chew, fr_chew = newmark(g_chew, dt, int(T / dt), 5, "chewing")
    np.savez_compressed("results/chew_frames.npz", frames=fr_chew, t=np.array(h_chew["t"]))
    out["chewing"] = h_chew
h_chew = out["chewing"]

# ------------------------------------------------------------------ 3. sudden hard bite
# rise times: 0.25 ms (impact-like, e.g. biting unexpectedly on a hard object) and 2 ms
for tr, key in ((0.25e-3, "impact_fast"), (2.0e-3, "impact_slow")):
    g_imp = (lambda trr: (lambda t: min(t / trr, 1.0)))(tr)
    dt2 = 2.5e-5
    h_imp, fr_imp = newmark(g_imp, dt2, int(8e-3 / dt2), 4, f"sudden bite tr={tr*1e3} ms")
    if key == "impact_fast":
        np.savez_compressed("results/impact_frames.npz", frames=fr_imp[:60], t=np.array(h_imp["t"])[:60])
    h_imp["rise_time_s"] = tr
    h_imp["DAF_plate"] = max(h_imp["plate_vm"]) / ref["plate_vm_max"]
    h_imp["DAF_ifm"] = max(h_imp["ifm"]) / ref["ifm_um"]
    out[key] = h_imp
    del fr_imp; gc.collect()
out.pop("impact", None); out.pop("impact_DAF_plate", None); out.pop("impact_DAF_ifm", None)
out["chew_peak_ratio_plate"] = max(h_chew["plate_vm"]) / ref["plate_vm_max"]
json.dump(out, open("results/dynamic.json", "w"), indent=1)
print({k: (out[k]["DAF_plate"], out[k]["DAF_ifm"]) for k in ("impact_fast", "impact_slow")}, out["chew_peak_ratio_plate"])
