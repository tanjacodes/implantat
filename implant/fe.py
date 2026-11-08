"""
Linear-elastic 3-D finite element solver (scikit-fem, quadratic tetrahedra)
for the fractured mandible + 3-D strut plate model built by geometry.py.

Units: mm, N, MPa, s, tonne (density in t/mm^3) -> consistent set.
"""
import math, time
import numpy as np
import meshio
import scipy.sparse as sp
from scipy.spatial import cKDTree
import pypardiso
from skfem import (MeshTet, ElementTetP2, ElementVector, Basis, FacetBasis,
                   BilinearForm, LinearForm, asm, condense)
from skfem.helpers import dot, ddot, sym_grad, trace, eye

from geometry import DEFAULT, frames, fracture_geom

# ----------------------------------------------------------------- materials
MAT = {
    "cortical":   dict(E=13700.0, nu=0.30, rho=1.90e-9),
    "cancellous": dict(E=1370.0,  nu=0.30, rho=0.70e-9),
    "titanium":   dict(E=105000.0, nu=0.34, rho=4.51e-9),   # CP-Ti grade 4 (ASTM F67)
    "gap":        dict(E=3.0,     nu=0.40, rho=1.05e-9),    # granulation tissue / haematoma
}
TI_YIELD = 483.0       # MPa, ASTM F67 grade 4 minimum 0.2% proof stress
TI_FATIGUE = 240.0     # MPa, conservative endurance estimate used for design
CORTEX_T = 2.0         # mm cortical thickness

# muscle maximal forces per side (N) and directions (lateral, anterior, superior)
MUSCLES = {
    "masseter":          dict(F=136.0, dir=(0.21, 0.42, 0.89)),
    "medial_pterygoid":  dict(F=76.8,  dir=(-0.55, 0.36, 0.76)),
    "temporalis":        dict(F=176.6, dir=(0.22, -0.50, 0.83)),
}
# activation factors (working side W / balancing side B)
CASES = {
    "incisal":   dict(bite="I",  act={"L": (0.40, 0.75, 0.15), "R": (0.40, 0.75, 0.15)}),
    "left_molar":  dict(bite="L", act={"L": (1.0, 1.0, 1.0), "R": (0.6, 0.6, 0.6)}),
    "right_molar": dict(bite="R", act={"L": (0.6, 0.6, 0.6), "R": (1.0, 1.0, 1.0)}),
}
MOLAR_S = 38.0   # first molar position along the body (local s, mm)


def lame(E, nu):
    return E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))


class Model:
    def __init__(self, msh, P=None, verbose=True):
        self.P = {**DEFAULT, **(P or {})}
        self.v = verbose
        t0 = time.time()
        m = meshio.read(msh)
        T, tags = [], []
        for cb, tg in zip(m.cells, m.cell_data["gmsh:physical"]):
            if cb.type == "tetra10":
                T.append(cb.data[:, :4]); tags.append(tg)
        T = np.vstack(T); tags = np.concatenate(tags)
        used, inv = np.unique(T, return_inverse=True)
        self.mesh = MeshTet(m.points[used].T.copy(), inv.reshape(T.shape).T.copy())
        self.tag = tags  # 1 bone 2 plate 3 screws 4 gap
        self.e = ElementVector(ElementTetP2())
        self.basis = Basis(self.mesh, self.e, intorder=2)
        self._assign_bone()
        self._log(f"mesh: {self.mesh.nelements} P2 tets, {self.basis.N} dofs  ({time.time()-t0:.1f}s)")

    def _log(self, s):
        if self.v: print(s, flush=True)

    # ---------------------------------------------------- local coordinates
    def to_local(self, X, side):
        X = np.array(X, dtype=float).copy()
        if side == "R":
            X[0] = -X[0]
        return frames(self.P).T @ X

    def _assign_bone(self):
        """cortical shell = bone within CORTEX_T of the outer bone surface"""
        mesh = self.mesh
        f2t = mesh.f2t
        bone = self.tag == 1
        tg = np.zeros((2, f2t.shape[1]), dtype=int)
        for k in (0, 1):
            ok = f2t[k] >= 0
            tg[k, ok] = self.tag[f2t[k, ok]]
        a, b = tg
        outer = ((a == 1) & np.isin(b, [0, 2])) | ((b == 1) & np.isin(a, [0, 2]))
        F = mesh.facets[:, outer]
        p = mesh.p
        bary = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1/3, 1/3, 1/3],
                         [.5, .5, 0], [0, .5, .5], [.5, 0, .5],
                         [2/3, 1/6, 1/6], [1/6, 2/3, 1/6], [1/6, 1/6, 2/3]])
        pts = np.einsum("kj,ijn->nki", bary, p[:, F]).reshape(-1, 3)
        tree = cKDTree(pts)
        cen = p[:, mesh.t].mean(axis=1).T
        d, _ = tree.query(cen)
        self.cortical = bone & (d < CORTEX_T)
        self.cancellous = bone & ~self.cortical
        self.dist_surf = d
        self.outer_facets = np.where(outer)[0]

    def material_groups(self):
        return {"cortical": np.where(self.cortical)[0],
                "cancellous": np.where(self.cancellous)[0],
                "titanium": np.where((self.tag == 2) | (self.tag == 3))[0],
                "gap": np.where(self.tag == 4)[0]}

    # ---------------------------------------------------- matrices
    def assemble(self, mass=False):
        t0 = time.time()
        K = None; M = None
        for name, el in self.material_groups().items():
            if len(el) == 0: continue
            lam, mu = lame(MAT[name]["E"], MAT[name]["nu"])
            b = self.basis.with_elements(el)

            @BilinearForm
            def a(u, v, w):
                e = sym_grad(u)
                return ddot(2 * mu * e + lam * eye(trace(e), 3), sym_grad(v))
            Ki = asm(a, b)
            K = Ki if K is None else K + Ki
            if mass:
                rho = MAT[name]["rho"]
                @BilinearForm
                def m(u, v, w):
                    return rho * dot(u, v)
                Mi = asm(m, b)
                M = Mi if M is None else M + Mi
        self.K = K.tocsr(); self.M = M.tocsr() if mass else None
        self._log(f"assembly {time.time()-t0:.1f}s")

    # ---------------------------------------------------- boundary regions
    def facets_where(self, fn, bone_only=True):
        mesh = self.mesh
        bf = mesh.boundary_facets()
        if bone_only:
            el = mesh.f2t[0, bf]
            bf = bf[self.tag[el] == 1]
        c = mesh.p[:, mesh.facets[:, bf]].mean(axis=1)
        return bf[fn(c)]

    def region_condyle(self, side):
        P = self.P
        sc = P["Lb"] + P["ramus_lean"] - 6.0; zc = P["ramus_h"] + 1.0
        def fn(c):
            L = self.to_local(c, side)
            r = np.hypot(L[0] - sc, L[2] - zc)
            return (r > P["cond_r"] - 0.3) & (L[2] > zc + 1.0) & (np.abs(L[1]) < P["cond_len"] / 2 - 0.5)
        return self.facets_where(fn)

    def region_bite(self, which, radius=3.0):
        P = self.P
        if which == "I":
            top = self.facets_where(lambda c: c[2] > P["H"] - 0.05)
            cc = self.mesh.p[:, self.mesh.facets[:, top]].mean(axis=1)
            ymax = cc[1][np.abs(cc[0]) < 2].max()
            target = np.array([0.0, ymax - 3.0, P["H"]])
        else:
            target = (frames(P) @ np.array([MOLAR_S, 0.0, P["H"]]))
            if which == "R": target[0] *= -1
        fn = lambda c: (c[2] > P["H"] - 0.05) & (np.hypot(c[0] - target[0], c[1] - target[1]) < radius)
        return self.facets_where(fn), target

    def region_muscle(self, name, side):
        P = self.P; Lb, T, H = P["Lb"], P["T"], P["H"]
        s_top0 = Lb - P["ramus_w"] + P["ramus_lean"]
        if name == "masseter":
            f = lambda L: (np.abs(L[1] + T / 2) < 0.15) & (L[0] > Lb - 26) & (L[2] < 30)
        elif name == "medial_pterygoid":
            f = lambda L: (np.abs(L[1] - T / 2) < 0.15) & (L[0] > Lb - 24) & (L[2] < 28)
        else:
            f = lambda L: (L[2] > P["ramus_h"] - 1.5) & (L[0] < s_top0 + 11) & (L[0] > s_top0 - 4)
        return self.facets_where(lambda c: f(self.to_local(c, side)))

    # ---------------------------------------------------- load vector
    def load_vector(self, case, scale=1.0):
        f = np.zeros(self.basis.N)
        self.muscle_info = []
        for side in ("L", "R"):
            for k, (mname, mu) in enumerate(MUSCLES.items()):
                F = mu["F"] * CASES[case]["act"][side][k] * scale
                if F == 0: continue
                lat, ant, sup = mu["dir"]
                d = np.array([-lat if side == "L" else lat, ant, sup]); d /= np.linalg.norm(d)
                fac = self.region_muscle(mname, side)
                fb = FacetBasis(self.mesh, self.e, facets=fac, intorder=3)
                area = asm(LinearForm(lambda v, w: v), FacetBasis(self.mesh, ElementTetP2(), facets=fac)).sum()
                tr = F * d / area
                f += asm(LinearForm(lambda v, w: tr[0] * v[0] + tr[1] * v[1] + tr[2] * v[2]), fb)
                self.muscle_info.append((side, mname, F, area))
        return f

    def dirichlet(self, case):
        D = [self.basis.get_dofs(facets=self.region_condyle(s)).all() for s in ("L", "R")]
        bf, self.bite_point = self.region_bite(CASES[case]["bite"])
        Db = self.basis.get_dofs(facets=bf).all("u^3")
        self.D_cond = np.unique(np.concatenate(D)); self.D_bite = np.unique(Db)
        return np.unique(np.concatenate(D + [Db]))

    # ---------------------------------------------------- linear solver
    def factor_solve(self, I, rhs, A=None):
        """symmetric positive-definite PARDISO (mtype=2) on the free dofs"""
        import gc
        A = self.K if A is None else A
        Aii = sp.triu(A[I][:, I]).tocsr(); gc.collect()
        s = pypardiso.PyPardisoSolver(mtype=2)
        s.set_iparm(1, 1); s.set_iparm(2, 3); s.set_iparm(24, 1)
        s.factorize(Aii)
        x = s.solve(Aii, rhs)
        self._solver, self._Aii = s, Aii
        return x

    def resolve(self, rhs):
        return self._solver.solve(self._Aii, rhs)

    def free_solver(self):
        try:
            self._solver.free_memory(everything=True)
        except Exception:
            pass
        self._solver = None; self._Aii = None

    # ---------------------------------------------------- static solve
    def solve_static(self, case, scale=1.0):
        t0 = time.time()
        f = self.load_vector(case, scale)
        D = self.dirichlet(case)
        u = np.zeros(self.basis.N)
        I = self.basis.complement_dofs(D)
        u[I] = self.factor_solve(I, f[I])
        r = self.K @ u - f
        self.u = u
        out = dict(case=case,
                   bite_force_N=float(-r[self.D_bite].sum()),
                   muscle_total_N=float(sum(x[2] for x in self.muscle_info)))
        # condyle reactions per side (vector)
        for s in ("L", "R"):
            Dc = self.basis.get_dofs(facets=self.region_condyle(s))
            vec = [float(r[Dc.all(f"u^{i}")].sum()) for i in (1, 2, 3)]
            out[f"condyle_{s}_reaction_N"] = vec
        self._log(f"  {case}: solve {time.time()-t0:.1f}s, bite reaction {out['bite_force_N']:.1f} N")
        return u, out

    # ---------------------------------------------------- post-processing
    def fields(self, u):
        """element-wise (max over quadrature points) von Mises and principal strains"""
        ui = self.basis.interpolate(u)
        G = ui.grad  # (3,3,nel,nqp)
        eps = 0.5 * (G + np.swapaxes(G, 0, 1))
        E = np.empty(self.mesh.nelements); nu = np.empty(self.mesh.nelements)
        for name, el in self.material_groups().items():
            E[el] = MAT[name]["E"]; nu[el] = MAT[name]["nu"]
        lam, mu = lame(E, nu)
        tr = np.trace(eps)
        sig = 2 * mu[None, None, :, None] * eps + (lam[:, None] * tr)[None, None] * np.eye(3)[:, :, None, None]
        s = sig
        vm = np.sqrt(0.5 * ((s[0, 0] - s[1, 1]) ** 2 + (s[1, 1] - s[2, 2]) ** 2 + (s[2, 2] - s[0, 0]) ** 2)
                     + 3 * (s[0, 1] ** 2 + s[1, 2] ** 2 + s[2, 0] ** 2))
        ep = np.linalg.eigvalsh(np.moveaxis(eps, (0, 1), (-2, -1)))  # (nel,nqp,3) ascending
        e1 = ep[..., 2]; e3 = ep[..., 0]
        return dict(vm_max=vm.max(1), vm_mean=vm.mean(1), e1_max=e1.max(1), e3_min=e3.min(1),
                    e1_mean=e1.mean(1), e3_mean=e3.mean(1))

    def ifm(self, u):
        """interfragmentary motion across the fracture gap (vertex based)"""
        P = self.P
        A, B, d, n = fracture_geom(P)
        Cf = 0.5 * (A + B)
        gap_el = np.where(self.tag == 4)[0]
        gv = np.unique(self.mesh.t[:, gap_el])
        L = self.to_local(self.mesh.p[:, gv], "L")
        sd = (L[0] - Cf[0]) * n[0] + (L[2] - Cf[1]) * n[1]
        g = P["gap"]
        fa = gv[sd < -g / 2 + 0.05]; fb = gv[sd > g / 2 - 0.05]
        U = u[self.basis.nodal_dofs]  # (3, nverts)
        tree = cKDTree(self.mesh.p[:, fb].T)
        dd, j = tree.query(self.mesh.p[:, fa].T)
        ok = dd < 2.5 * g
        rel = U[:, fb[j[ok]]] - U[:, fa[ok]]
        R = frames(P)
        nglob = R @ np.array([n[0], 0, n[1]])
        rn = nglob @ rel
        rt = np.linalg.norm(rel - np.outer(nglob, rn), axis=0)
        mag = np.linalg.norm(rel, axis=0)
        Lf = L[:, np.isin(gv, fa[ok])] if False else self.to_local(self.mesh.p[:, fa[ok]], "L")
        return dict(max_um=float(mag.max() * 1e3), p95_um=float(np.percentile(mag, 95) * 1e3),
                    mean_um=float(mag.mean() * 1e3),
                    normal_max_open_um=float(rn.max() * 1e3), normal_max_close_um=float(-rn.min() * 1e3),
                    shear_max_um=float(rt.max() * 1e3),
                    _pts=Lf, _mag=mag)

    def summary(self, u, out):
        F = self.fields(u)
        imp_p = self.tag == 2; imp_s = self.tag == 3
        cort = self.cortical; bone = self.tag == 1
        # peri-screw bone: bone elements touching screw elements' vertices
        sv = np.unique(self.mesh.t[:, imp_s])
        peri = bone & np.isin(self.mesh.t, sv).any(axis=0)
        U = u[self.basis.nodal_dofs]
        res = dict(out)
        res.update(
            umax_mm=float(np.linalg.norm(U, axis=0).max()),
            plate_vm_max=float(F["vm_max"][imp_p].max()),
            plate_vm_p99=float(np.percentile(F["vm_max"][imp_p], 99.0)),
            screw_vm_max=float(F["vm_max"][imp_s].max()),
            screw_vm_p99=float(np.percentile(F["vm_max"][imp_s], 99.0)),
            bone_e1_p99_ue=float(np.percentile(F["e1_max"][bone], 99.0) * 1e6),
            bone_e3_p01_ue=float(np.percentile(F["e3_min"][bone], 1.0) * 1e6),
            peri_screw_e1_max_ue=float(F["e1_max"][peri].max() * 1e6),
            peri_screw_e1_p95_ue=float(np.percentile(F["e1_max"][peri], 95) * 1e6),
            peri_screw_e3_p05_ue=float(np.percentile(F["e3_min"][peri], 5) * 1e6),
            cortical_vm_p99=float(np.percentile(F["vm_max"][cort], 99.0)),
        )
        I = self.ifm(u)
        res.update({f"ifm_{k}": v for k, v in I.items() if not k.startswith("_")})
        return res, F, I
