"""Figures for the plate optimisation."""
import json
import numpy as np
import pyvista as pv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fe import Model
from vis import grid, sub, view_lateral_left, BONE

S = json.load(open("results/optim_sweeps.json")) + json.load(open("results/optim_combos.json"))
LF = 2.0
fig, ax = plt.subplots(figsize=(8.2, 5.4), dpi=150)
ax.axhspan(240, 320, color="#c44e52", alpha=0.08); ax.axhline(240, color="#c44e52", lw=1, ls="--")
ax.text(84, 243, "fatigue limit 240 MPa (design load = 2x physiological)", fontsize=7.5, color="#c44e52")
grp = {"t=": "#4c72b0", "bar=": "#55a868", "struts": "#8172b2", "4-hole": "#dd8452", "height": "#937860",
       "bars 13": "#937860", "baseline": "k"}
for r in S:
    n = r["name"]
    col = next((c for k, c in grp.items() if n.startswith(k)), "#c44e52")
    combo = "," in n
    mk = "*" if combo else ("s" if n == "baseline" else "o")
    ax.plot(r["volume_mm3"], r["design"]["vm"], mk, color=col if not combo else "#c44e52", ms=13 if combo else 7,
            mec="k" if combo else col, mew=0.6)
    ax.annotate(n, (r["volume_mm3"], r["design"]["vm"]), xytext=(5, 3), textcoords="offset points", fontsize=6.5)
best = [r for r in S if r["name"] == "4-hole, bar 2.4"][0]
ax.annotate("selected design", (best["volume_mm3"], best["design"]["vm"]), xytext=(25, -38), textcoords="offset points",
            fontsize=8, color="#c44e52", arrowprops=dict(arrowstyle="->", color="#c44e52"))
ax.set_xlabel("Implant volume (plate + screws), mm³"); ax.set_ylabel("Peak plate von Mises at design load (MPa)")
ax.set_title("Design space: one-factor sweeps (dots) and combined candidates (stars)\n"
             "left molar clench, all designs: IFM ≤ 24 µm, peri-screw strain ≤ 440 µε (limits 150 µm / 3000 µε)", fontsize=9)
ax.set_xlim(75, 230); ax.set_ylim(140, 310); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("figures/fig11_design_space.png"); plt.close(fig)

# sensitivity bars (relative to baseline)
base = [r for r in S if r["name"] == "baseline"][0]
names = [r["name"] for r in S if r["name"] != "baseline" and "," not in r["name"]]
dv = [100 * (r["volume_mm3"] / base["volume_mm3"] - 1) for r in S if r["name"] in names]
ds = [100 * (r["design"]["vm"] / base["design"]["vm"] - 1) for r in S if r["name"] in names]
fig, ax = plt.subplots(figsize=(8.2, 4.6), dpi=150)
y = np.arange(len(names))
ax.barh(y - 0.2, dv, 0.4, color="#8d99a6", label="Δ volume (%)")
ax.barh(y + 0.2, ds, 0.4, color="#c44e52", label="Δ peak plate stress (%)")
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8); ax.axvline(0, color="k", lw=0.8)
ax.invert_yaxis(); ax.legend(fontsize=8, frameon=False); ax.grid(axis="x", alpha=0.3)
ax.set_title("Sensitivity of the baseline 8-hole plate (t 1.0, bar 1.6 mm) to single changes", fontsize=9)
fig.tight_layout(); fig.savefig("figures/fig12_sensitivity.png"); plt.close(fig)

# baseline vs optimised plate stress, left molar
p = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1800, 820), border=False); p.set_background("white")
for i, (msh, P, ufile, lab) in enumerate([("base.msh", {}, "results/u_base.npz", "Baseline: 8-hole 3-D plate, bar 1.6 mm (168 mm³)"),
                                          ("opt.msh", dict(holes_p=[-5.5, 5.5], wb=2.4), "results/u_opt.npz",
                                           "Optimised: 4-hole 3-D plate, bar 2.4 mm (103 mm³)")]):
    m = Model(msh, P, verbose=False)
    u = np.load(ufile)["left_molar"]
    F = m.fields(u)
    g = grid(m, u, cell={"von Mises (MPa)": F["vm_max"]})
    p.subplot(0, i)
    p.add_mesh(sub(g, [1]).extract_surface(algorithm="dataset_surface"), color=BONE, opacity=0.35)
    p.add_mesh(sub(g, [2, 3]).extract_surface(algorithm="dataset_surface"), scalars="von Mises (MPa)", cmap="turbo",
               clim=(0, 110), scalar_bar_args=dict(title="von Mises (MPa)", vertical=False, position_x=0.2, width=0.6))
    view_lateral_left(p, m, zoom=2.3)
    p.add_text(f"{lab}\npeak {F['vm_max'][(m.tag==2)|(m.tag==3)].max():.0f} MPa, left molar clench", font_size=11)
p.screenshot("figures/fig13_baseline_vs_optimised.png"); p.close()
print("done")
