"""
Parametric geometry + mesh of a simplified human mandible with a left
mandibular-angle fracture fixed by a 3-D (strut / "box") titanium plate.

Units: mm.  Global axes: x = right(+), y = anterior(+), z = superior(+).
The left hemimandible is built in a local frame (s, w, z):
  s : along the body, from symphysis (s=0) backwards to the gonion (s=Lb)
  w : medial (+) ; the lateral (buccal) face is at w = -T/2
  z : superior
and mapped to global coordinates with a rigid rotation about z.
The right hemimandible is the mirror image (x -> -x).

Volumes (physical groups):
  1 = bone (cortical/cancellous assigned later by distance to the surface)
  2 = plate (titanium)
  3 = screws (titanium, bonded to plate and bone)
  4 = fracture-gap tissue (soft granulation tissue / haematoma)
"""
import math, json, sys
import numpy as np
import gmsh

DEFAULT = dict(
    # mandible
    theta_deg=35.0,   # angle of each hemibody to the mid-sagittal plane
    Lb=72.0,          # body length symphysis -> gonion
    H=28.0,           # body height
    T=11.0,           # body / ramus thickness (bucco-lingual)
    ramus_w=32.0,     # antero-posterior width of the ramus at its base
    ramus_h=58.0,     # ramus height
    ramus_lean=20.0,  # backwards lean of the ramus (gonial angle ~110 deg)
    cond_r=5.0, cond_len=20.0,
    # fracture (in local s-z plane): from alveolar crest (Lb-30, H) to inferior border (Lb-12, 0)
    frac_top_s=-30.0, frac_bot_s=-12.0, gap=0.6,
    # 3-D strut plate (plate frame: p along bars, q across)
    tp=1.0,           # plate thickness
    wb=1.6,           # bar width
    ws=1.4,           # strut width
    ring_r=2.0,       # outer radius of the screw rings
    holes_p=[-11.0, -5.0, 5.0, 11.0],
    Q=5.0,            # half distance between the two bars
    struts="all",     # "all" | "outer" | "inner" | "none"
    diagonals=False,  # add an X-brace between the inner holes
    alpha_deg=None,   # bar angle; None = perpendicular to the fracture
    center_z=14.0,    # plate centre height on the fracture line
    screw_r=1.0, screw_len=6.0,   # 2.0 mm monocortical screws, 6 mm in bone
    # mesh
    h_implant=0.45, h_near=1.0, h_far=3.2, d_near=3.0, d_far=22.0,
)


def frames(P):
    th = math.radians(P["theta_deg"])
    # columns: s-dir, w-dir (medial), z-dir   (det = +1)
    R = np.array([[-math.sin(th), math.cos(th), 0.0],
                  [-math.cos(th), -math.sin(th), 0.0],
                  [0.0, 0.0, 1.0]])
    return R


def local_to_global(P, slz, side="L"):
    R = frames(P)
    X = np.atleast_2d(slz) @ R.T
    if side == "R":
        X = X * np.array([-1, 1, 1])
    return X


def fracture_geom(P):
    Lb, H = P["Lb"], P["H"]
    A = np.array([Lb + P["frac_top_s"], H])     # (s,z) top point
    B = np.array([Lb + P["frac_bot_s"], 0.0])   # bottom point
    d = (A - B) / np.linalg.norm(A - B)         # direction along fracture (upwards)
    n = np.array([d[1], -d[0]])                 # normal (points posterior/down)
    if n[0] < 0:
        n = -n
    return A, B, d, n


def plate_layout(P):
    """Return plate centre (s,z), bar angle alpha (rad), and hole centres (s,z)."""
    A, B, d, n = fracture_geom(P)
    zc = P["center_z"]
    t = (zc - B[1]) / d[1]
    C = B + t * d
    alpha = math.atan2(n[1], n[0]) if P["alpha_deg"] is None else math.radians(P["alpha_deg"])
    ca, sa = math.cos(alpha), math.sin(alpha)
    holes = []
    for p in P["holes_p"]:
        for q in (-P["Q"], P["Q"]):
            holes.append((C[0] + p * ca - q * sa, C[1] + p * sa + q * ca))
    return C, alpha, np.array(holes)


def build(P, msh_path, verbose=False):
    P = {**DEFAULT, **P}
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
    gmsh.option.setNumber("General.NumThreads", 2)
    gmsh.model.add("mandible")
    occ = gmsh.model.occ
    Lb, H, T = P["Lb"], P["H"], P["T"]

    # ---------------- left hemimandible in local coordinates -------------
    body = occ.addBox(-6.0, -T / 2, 0.0, Lb + 6.0, T, H)
    rect = occ.addRectangle(Lb - P["ramus_w"], -T / 2, 0.0, P["ramus_w"], T)
    ram = occ.extrude([(2, rect)], P["ramus_lean"], 0.0, P["ramus_h"])
    ramus = [dt for dt in ram if dt[0] == 3][0][1]
    zt = P["ramus_h"]
    s_top0 = Lb - P["ramus_w"] + P["ramus_lean"]          # anterior end of ramus top
    # condyle (cylinder, axis medio-lateral) at the posterior top
    sc = Lb + P["ramus_lean"] - 6.0
    cond = occ.addCylinder(sc, -P["cond_len"] / 2, zt + 1.0, 0, P["cond_len"], 0, P["cond_r"])
    # coronoid process (temporalis insertion) at the anterior top
    rc = occ.addRectangle(s_top0, -T / 2 + 1.0, zt - 1.0, 10.0, T - 2.0)
    cor = occ.extrude([(2, rc)], -2.0, 0.0, 9.0)
    coronoid = [dt for dt in cor if dt[0] == 3][0][1]
    left, _ = occ.fuse([(3, body)], [(3, ramus), (3, cond), (3, coronoid)])
    occ.synchronize()

    # fracture slab (local), plate, screws (local)
    A, B, d, n = fracture_geom(P)
    Cf = 0.5 * (A + B)
    slab = occ.addBox(-25.0, -T, -P["gap"] / 2, 50.0, 2 * T, P["gap"])
    phi = math.atan2(-d[1], d[0])  # rotate local x-axis onto d (rotation about y)
    occ.rotate([(3, slab)], 0, 0, 0, 0, 1, 0, phi)
    occ.translate([(3, slab)], Cf[0], 0, Cf[1])

    C, alpha, holes = plate_layout(P)
    tp, wb, ws, Q, rr = P["tp"], P["wb"], P["ws"], P["Q"], P["ring_r"]
    y0 = -T / 2 - tp
    hp = P["holes_p"]
    pmin, pmax = min(hp), max(hp)
    parts = []
    for q in (-Q, Q):   # bars
        parts.append(occ.addBox(pmin, y0, q - wb / 2, pmax - pmin, tp, wb))
    for p in hp:        # rings
        for q in (-Q, Q):
            parts.append(occ.addCylinder(p, y0, q, 0, tp, 0, rr))
    sel = {"all": hp, "outer": [pmin, pmax], "inner": [p for p in hp if abs(p) < max(abs(pmin), pmax)],
           "none": []}[P["struts"]]
    for p in sel:
        parts.append(occ.addBox(p - ws / 2, y0, -Q, ws, tp, 2 * Q))
    if P["diagonals"]:
        pin = sorted([p for p in hp if abs(p) == min(abs(x) for x in hp)])
        L = math.hypot(pin[1] - pin[0], 2 * Q)
        for sgn in (1, -1):
            b = occ.addBox(-L / 2, y0, -ws / 2, L, tp, ws)
            occ.rotate([(3, b)], 0, 0, 0, 0, 1, 0, -sgn * math.atan2(2 * Q, pin[1] - pin[0]))
            parts.append(b)
    plate, _ = occ.fuse([(3, parts[0])], [(3, t) for t in parts[1:]])
    screws = []
    for p in hp:
        for q in (-Q, Q):
            screws.append(occ.addCylinder(p, y0, q, 0, tp + P["screw_len"], 0, P["screw_r"]))
    occ.rotate(plate + [(3, s) for s in screws], 0, 0, 0, 0, 1, 0, -alpha)
    occ.translate(plate + [(3, s) for s in screws], C[0], 0, C[1])

    # ---------------- map to global, mirror, fuse ------------------------
    R = frames(P)
    M = [R[0, 0], R[0, 1], R[0, 2], 0, R[1, 0], R[1, 1], R[1, 2], 0, R[2, 0], R[2, 1], R[2, 2], 0]
    every = left + [(3, slab)] + plate + [(3, s) for s in screws]
    occ.affineTransform(every, M)
    right = occ.copy(left)
    occ.mirror(right, 1, 0, 0, 0)
    mand, _ = occ.fuse(left, right)
    occ.synchronize()

    # gap tissue = slab ∩ bone ; bone = bone − slab
    gapv, _ = occ.intersect(mand, [(3, slab)], removeObject=False, removeTool=False)
    bone, _ = occ.cut(mand, [(3, slab)], removeObject=True, removeTool=True)
    occ.synchronize()
    implant_vol = sum(occ.getMass(3, t) for _, t in plate)

    objs = bone + gapv + plate + [(3, s) for s in screws]
    out, omap = occ.fragment(objs, [])
    occ.synchronize()
    nb, ng, npl = len(bone), len(gapv), len(plate)
    vb = sorted({t for m in omap[:nb] for _, t in m})
    vg = sorted({t for m in omap[nb:nb + ng] for _, t in m})
    vp = sorted({t for m in omap[nb + ng:nb + ng + npl] for _, t in m})
    vs = sorted({t for m in omap[nb + ng + npl:] for _, t in m})
    # fragment may assign overlapping regions to several inputs: screws win over bone/plate
    vb = [t for t in vb if t not in vs and t not in vp and t not in vg]
    vp = [t for t in vp if t not in vs]
    vg = [t for t in vg if t not in vs]
    gmsh.model.addPhysicalGroup(3, vb, 1, "bone")
    gmsh.model.addPhysicalGroup(3, vp, 2, "plate")
    gmsh.model.addPhysicalGroup(3, vs, 3, "screws")
    gmsh.model.addPhysicalGroup(3, vg, 4, "gap")

    # ---------------- mesh size fields -----------------------------------
    def bnd(vols):
        return [t for _, t in gmsh.model.getBoundary([(3, v) for v in vols], oriented=False, combined=False)]
    s_impl = sorted(set(bnd(vp + vs)))
    s_gap = sorted(set(bnd(vg)))
    f = gmsh.model.mesh.field
    f1 = f.add("Distance"); f.setNumbers(f1, "SurfacesList", s_impl + s_gap); f.setNumber(f1, "Sampling", 40)
    f2 = f.add("Threshold"); f.setNumber(f2, "InField", f1)
    f.setNumber(f2, "SizeMin", P["h_near"]); f.setNumber(f2, "SizeMax", P["h_far"])
    f.setNumber(f2, "DistMin", P["d_near"]); f.setNumber(f2, "DistMax", P["d_far"])
    f3 = f.add("Constant"); f.setNumber(f3, "VIn", P["h_implant"]); f.setNumber(f3, "VOut", 1e9)
    f.setNumbers(f3, "VolumesList", vp + vs)
    f4 = f.add("Constant"); f.setNumber(f4, "VIn", min(P["h_near"], 0.6)); f.setNumber(f4, "VOut", 1e9)
    f.setNumbers(f4, "VolumesList", vg)
    f5 = f.add("Min"); f.setNumbers(f5, "FieldsList", [f2, f3, f4])
    f.setAsBackgroundMesh(f5)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)   # HXT (parallel)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)  # straight-sided P2
    gmsh.model.mesh.setOrder(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(msh_path)
    ntet = sum(len(e[0]) for e in [gmsh.model.mesh.getElementsByType(11)])
    gmsh.finalize()
    meta = dict(params=P, implant_volume_mm3=implant_vol, n_tet10=int(ntet),
                holes_local=plate_layout(P)[2].tolist())
    return meta


if __name__ == "__main__":
    P = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
    m = build(P, sys.argv[1] if len(sys.argv) > 1 else "mandible.msh", verbose=True)
    print(json.dumps({k: v for k, v in m.items() if k != "params"}, indent=1))
