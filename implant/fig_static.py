"""Figures for the baseline static analysis."""
import json, os
import numpy as np
import pyvista as pv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fe import Model, CASES, MUSCLES
from geometry import frames, plate_layout, fracture_geom
from vis import grid, sub, view_lateral_left, view_oblique, BONE, TI

os.makedirs("figures", exist_ok=True)
m = Model("base.msh", verbose=False)
U = np.load("results/u_base.npz")
R = json.load(open("results/static_base.json"))
W, Hh = 1400, 1000
NICE = {"incisal": "Incisal clench", "left_molar": "Left (fractured-side) molar clench",
        "right_molar": "Right (contralateral) molar clench"}


def facet_poly(fac):
    tri = m.mesh.facets[:, fac].T
    faces = np.hstack([np.full((len(tri), 1), 3), tri]).ravel()
    return pv.PolyData(m.mesh.p.T, faces)


# ---------------- Fig 1: model, materials and boundary conditions -------------
g = grid(m)
p = pv.Plotter(off_screen=True, window_size=(W, Hh))
p.set_background("white")
p.add_mesh(sub(g, [1]).extract_surface(), color=BONE, smooth_shading=False)
p.add_mesh(sub(g, [2, 3]).extract_surface(), color=TI, specular=0.6)
cols = {"masseter": "#e07b39", "medial_pterygoid": "#3a9d5d", "temporalis": "#7b5ea7"}
for s in "LR":
    p.add_mesh(facet_poly(m.region_condyle(s)), color="#1f5fbf")
    for k in MUSCLES:
        poly = facet_poly(m.region_muscle(k, s))
        p.add_mesh(poly, color=cols[k], opacity=0.85)
        c = poly.extract_surface().clean().points.mean(0) if poly.n_points else None
        lat, ant, sup = MUSCLES[k]["dir"]
        d = np.array([-lat if s == "L" else lat, ant, sup]); d /= np.linalg.norm(d)
        cc = m.mesh.p[:, m.mesh.facets[:, m.region_muscle(k, s)]].mean(axis=(1, 2))
        p.add_mesh(pv.Arrow(cc, d, scale=22, tip_length=0.3, shaft_radius=0.04), color=cols[k])
for b in "IL":
    bf, tgt = m.region_bite(b)
    p.add_mesh(pv.Arrow(tgt + np.array([0, 0, 16]), (0, 0, -1), scale=15, tip_length=0.35, shaft_radius=0.05), color="#c62828")
view_oblique(p, m, 1.25)
p.add_text("Bone (beige) · Ti plate+screws (grey) · condyles fixed (blue)\n"
           "Masseter (orange) · Medial pterygoid (green) · Temporalis (purple) · Bite supports (red)",
           font_size=11, position="upper_left")
p.screenshot("figures/fig01_model_bcs.png"); p.close()

# ---------------- Fig 2: mesh close-up at the fracture --------------------------
p = pv.Plotter(off_screen=True, window_size=(W, Hh)); p.set_background("white")
gb = sub(g, [1]).extract_surface()
p.add_mesh(gb, scalars=None, color=BONE, show_edges=True, edge_color="#8a7d66", line_width=0.4)
p.add_mesh(sub(g, [2, 3]).extract_surface(), color=TI, show_edges=True, edge_color="#333333", line_width=0.4)
p.add_mesh(sub(g, [4]).extract_surface(), color="#d84a4a")
view_lateral_left(p, m, zoom=2.6)
p.add_text("Quadratic tetrahedral mesh (local refinement at plate, screws and 0.6 mm fracture gap)", font_size=11)
p.screenshot("figures/fig02_mesh_closeup.png"); p.close()

# cross-section through the bone showing cortical shell
p = pv.Plotter(off_screen=True, window_size=(W, Hh)); p.set_background("white")
Rm = frames(m.P)
C, a, holes = plate_layout(m.P)
origin = Rm @ np.array([holes[2][0], 0, holes[2][1]])
normal = Rm @ np.array([np.cos(a), 0, np.sin(a)])
cut = g.slice(normal=normal, origin=origin)
b = 18.0
cut = cut.clip_box([origin[0]-b, origin[0]+b, origin[1]-b, origin[1]+b, origin[2]-b, origin[2]+b], invert=False)
mat_names = np.array(["cancellous", "cortical", "plate", "screw", "gap"])
lut = pv.LookupTable(cmap=["#d9b98c", "#f3ead8", "#8d99a6", "#56606b", "#d84a4a"], n_values=5)
lut.scalar_range = (0, 4)
from matplotlib.colors import ListedColormap
p.add_mesh(cut, scalars="mat", cmap=ListedColormap(["#c9a46c", "#f6efe0", "#8d99a6", "#3c444d", "#d84a4a"]), clim=(-0.5, 4.5),
           n_colors=5, preference="cell", lighting=False, show_edges=True, edge_color="#666666", line_width=0.3, show_scalar_bar=False)
p.camera.focal_point = origin
p.camera.position = origin - 150 * normal + np.array([0, 0, 30])
p.camera.up = (0, 0, 1); p.reset_camera(); p.camera.zoom(1.3)
p.add_text("Section through the posterior screw column:\n cortical shell 2 mm (light), cancellous core (tan), "
           "screws (dark grey), plate (grey)", font_size=11)
p.screenshot("figures/fig03_section_materials.png"); p.close()

# ---------------- Fig 4: von Mises in the implant for three cases ----------------
vmax = max(R[c]["plate_vm_max"] for c in CASES)
p = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(2100, 760), border=False)
p.set_background("white")
for i, c in enumerate(CASES):
    F = m.fields(U[c])
    gc = grid(m, U[c], cell={"von Mises (MPa)": F["vm_max"]})
    p.subplot(0, i)
    p.add_mesh(sub(gc, [1]).extract_surface(), color=BONE, opacity=0.35)
    p.add_mesh(sub(gc, [2, 3]).extract_surface(), scalars="von Mises (MPa)", cmap="turbo", clim=(0, 120),
               scalar_bar_args=dict(title="von Mises stress (MPa)", vertical=False, position_x=0.15, width=0.7,
                                    title_font_size=14, label_font_size=12, fmt="%.0f"))
    view_lateral_left(p, m, zoom=2.3)
    p.add_text(f"{NICE[c]}\nmax {R[c]['plate_vm_max']:.0f} MPa (99th pct {R[c]['plate_vm_p99']:.0f})", font_size=10)
p.screenshot("figures/fig04_plate_vonmises.png"); p.close()

# ---------------- Fig 5: displacement of the whole mandible ----------------------
c = "left_molar"
gc = grid(m, U[c], scale=400.0)
p = pv.Plotter(off_screen=True, window_size=(W, Hh)); p.set_background("white")
p.add_mesh(grid(m).extract_surface(), style="wireframe", color="#bbbbbb", opacity=0.15)
p.add_mesh(gc.extract_surface(), scalars="|u| (um)", cmap="viridis",
           scalar_bar_args=dict(title="|u| (micrometre)", vertical=True, position_x=0.86, position_y=0.2, height=0.6))
view_oblique(p, m, 1.25)
p.add_text(f"Displacement magnitude, left molar clench (deformation ×400, grey = undeformed)\n"
           f"max {R[c]['umax_mm']*1e3:.1f} micrometre", font_size=11)
p.screenshot("figures/fig05_displacement.png"); p.close()

# ---------------- Fig 6: bone strain around the fracture ------------------------
F = m.fields(U[c])
gc = grid(m, U[c], cell={"eps1": F["e1_max"] * 1e6, "eps3": F["e3_min"] * 1e6})
p = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(2000, 860), border=False); p.set_background("white")
for i, (k, lab, cm, clim) in enumerate([("eps1", "Max principal strain (microstrain)", "Reds", (0, 400)),
                                         ("eps3", "Min principal strain (microstrain)", "Blues_r", (-500, 0))]):
    p.subplot(0, i)
    p.add_mesh(sub(gc, [1]).extract_surface(), scalars=k, cmap=cm, clim=clim,
               scalar_bar_args=dict(title=lab, vertical=False, position_x=0.15, width=0.7))
    p.add_mesh(sub(gc, [2, 3]).extract_surface(), color=TI, opacity=0.25)
    view_lateral_left(p, m, zoom=2.4)
p.add_text("Bone strain, left molar clench (plate shown transparent)", font_size=11)
p.screenshot("figures/fig06_bone_strain.png"); p.close()

# ---------------- Fig 7: interfragmentary motion along the fracture -------------
A, B, d, n = fracture_geom(m.P)
fig, ax = plt.subplots(figsize=(7.5, 4.6), dpi=150)
colors = {"incisal": "#4c72b0", "left_molar": "#c44e52", "right_molar": "#55a868"}
for c in CASES:
    z = np.load(f"results/ifm_base_{c}.npz")
    L, mag = z["pts"], z["mag"] * 1e3
    pos = (L[0] - B[0]) * d[0] + (L[2] - B[1]) * d[1]      # distance along fracture from inferior border
    bins = np.linspace(0, np.linalg.norm(A - B), 18)
    idx = np.digitize(pos, bins)
    xb = [0.5 * (bins[i - 1] + bins[i]) for i in range(1, len(bins)) if np.any(idx == i)]
    yb = [mag[idx == i].max() for i in range(1, len(bins)) if np.any(idx == i)]
    ax.plot(xb, yb, "-o", ms=3, color=colors[c], label=NICE[c])
ax.axhline(150, color="k", ls="--", lw=1); ax.text(1, 140, "150 µm stability limit", va="top", fontsize=8)
ax.set_xlabel("Position along fracture line from inferior border (mm)")
ax.set_ylabel("Interfragmentary motion (µm)")
ax.set_yscale("log"); ax.set_ylim(1, 300)
ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8, frameon=False)
ax.set_title("Interfragmentary motion across the angle fracture (baseline plate)", fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig07_ifm.png"); plt.close(fig)
print("done")
