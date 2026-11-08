"""Figures and animations for the dynamic analyses."""
import json, os
import numpy as np
import pyvista as pv
import imageio.v2 as imageio
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fe import Model
from vis import grid, sub, view_lateral_left, view_oblique, BONE, TI

os.makedirs("figures", exist_ok=True); os.makedirs("animations", exist_ok=True)
m = Model("conv_medium.msh", dict(h_implant=0.55, h_near=1.3, h_far=3.8), verbose=False)
D = json.load(open("results/dynamic.json"))
plate_el = np.where((m.tag == 2) | (m.tag == 3))[0]
pb = m.basis.with_elements(plate_el)
E, nu = 105000.0, 0.34
lamT, muT = E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))


def vm_all(u):
    G = pb.interpolate(u).grad
    e = 0.5 * (G + np.swapaxes(G, 0, 1)); tr = np.trace(e)
    s = 2 * muT * e + lamT * tr * np.eye(3)[:, :, None, None]
    vm = np.sqrt(0.5 * ((s[0, 0] - s[1, 1]) ** 2 + (s[1, 1] - s[2, 2]) ** 2 + (s[2, 2] - s[0, 0]) ** 2)
                 + 3 * (s[0, 1] ** 2 + s[1, 2] ** 2 + s[2, 0] ** 2)).max(1)
    out = np.zeros(m.mesh.nelements); out[plate_el] = vm
    return out

import sys
ONLY_ANIM = "anim" in sys.argv
if not ONLY_ANIM:
    # ---------------------------------------------------------------- mode shapes
    Z = np.load("results/modes.npz")
    fr, V = Z["fr_occ"], Z["V_occ"]
    p = pv.Plotter(off_screen=True, shape=(2, 2), window_size=(1600, 1200), border=False); p.set_background("white")
    base = grid(m).extract_surface(algorithm="dataset_surface")
    for k in range(4):
        phi = V[:, k] / np.abs(V[:, k]).max()
        g = grid(m, phi, scale=8.0)
        p.subplot(k // 2, k % 2)
        p.add_mesh(base, style="wireframe", color="#cccccc", opacity=0.12)
        p.add_mesh(g.extract_surface(algorithm="dataset_surface"), scalars="|u| (um)", cmap="plasma", show_scalar_bar=False)
        view_oblique(p, m, 1.15)
        p.add_text(f"Mode {k+1}: {fr[k]:.0f} Hz", font_size=13)
    p.screenshot("figures/fig08_modes.png"); p.close()

    # mode-1 animation
    g0 = grid(m, V[:, 0] / np.abs(V[:, 0]).max(), scale=0.0)
    P0 = g0.points.copy(); U0 = g0.point_data["u"]
    p = pv.Plotter(off_screen=True, window_size=(900, 700)); p.set_background("white")
    surf = g0.extract_surface(algorithm="dataset_surface")
    ids = surf.point_data["vtkOriginalPointIds"]
    actor = p.add_mesh(surf, scalars="|u| (um)", cmap="plasma", show_scalar_bar=False)
    view_oblique(p, m, 1.15); p.add_text(f"Mode 1 ({fr[0]:.0f} Hz), occluded, amplified", font_size=11)
    w = imageio.get_writer("animations/mode1.gif", fps=20, loop=0)
    for ph in np.linspace(0, 2 * np.pi, 30, endpoint=False):
        surf.points = P0[ids] + 8.0 * np.sin(ph) * U0[ids]
        p.render(); w.append_data(p.screenshot(return_img=True))
    w.close(); p.close()

    # ---------------------------------------------------------------- chewing history
    h = D["chewing"]; ref = D["static_reference"]
    t = np.array(h["t"]) * 1e3
    fig, ax = plt.subplots(3, 1, figsize=(7.5, 7), dpi=150, sharex=True)
    ax[0].plot(t, np.array(h["g"]) * 100, color="#555"); ax[0].set_ylabel("Muscle activation (%)")
    ax[0].axvspan(0, 300, color="#4c72b0", alpha=0.08); ax[0].axvspan(300, 550, color="#c44e52", alpha=0.08)
    ax[0].text(150, 104, "closing", ha="center", fontsize=8); ax[0].text(425, 104, "power stroke", ha="center", fontsize=8)
    ax[0].text(675, 104, "opening", ha="center", fontsize=8); ax[0].set_ylim(0, 115)
    ax[1].plot(t, h["bite"], color="#c44e52"); ax[1].set_ylabel("Bite force (N)")
    ax1b = ax[1].twinx(); ax1b.plot(t, h["ifm"], color="#4c72b0", ls="--"); ax1b.set_ylabel("IFM (µm)", color="#4c72b0")
    ax[2].plot(t, h["plate_vm"], color="#8172b2", label="max"); ax[2].plot(t, h["plate_p99"], color="#8172b2", ls=":", label="99th percentile")
    ax[2].axhline(ref["plate_vm_max"], color="k", lw=0.8, ls="--"); ax[2].text(5, ref["plate_vm_max"] + 3, "static solution", fontsize=8)
    ax[2].set_ylabel("Plate von Mises (MPa)"); ax[2].set_xlabel("Time (ms)"); ax[2].legend(fontsize=8, frameon=False)
    for a in ax: a.grid(alpha=0.3)
    fig.suptitle("Transient chewing cycle, left (fractured-side) molar", fontsize=10)
    fig.tight_layout(); fig.savefig("figures/fig09_chewing_history.png"); plt.close(fig)

    # ---------------------------------------------------------------- impact history
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9), dpi=150)
    T1 = 1 / D["modal"]["occluded_Hz"][0] * 1e3
    for key, col in (("impact_fast", "#c44e52"), ("impact_slow", "#4c72b0")):
        h = D[key]; t = np.array(h["t"]) * 1e3; tr = h["rise_time_s"] * 1e3
        ax[0].plot(t, np.array(h["plate_vm"]) / ref["plate_vm_max"], color=col,
                   label=f"rise {tr:g} ms: DAF plate {h['DAF_plate']:.2f}, IFM {h['DAF_ifm']:.2f}")
        ax[0].plot(t, h["g"], color=col, lw=0.7, ls=":")
    ax[0].set_xlabel("Time (ms)"); ax[0].set_ylabel("Plate peak von Mises / static value")
    ax[0].legend(fontsize=7.5, frameon=False, loc="lower right"); ax[0].grid(alpha=0.3)
    ax[0].set_title("Sudden hard bite (dotted = applied load ramp)", fontsize=10)
    trr = np.linspace(0.02, 5, 300)
    ax[1].plot(trr, 1 + np.abs(np.sin(np.pi * trr)) / (np.pi * trr), color="k", label="undamped SDOF theory")
    for key, col in (("impact_fast", "#c44e52"), ("impact_slow", "#4c72b0")):
        x = D[key]["rise_time_s"] * 1e3 / T1
        ax[1].plot([x], [D[key]["DAF_plate"]], "o", color=col, ms=7, label=f"FE, t_r/T1 = {x:.2f}")
    ax[1].set_xlabel("Load rise time / first natural period  t_r / T1"); ax[1].set_ylabel("Dynamic amplification factor")
    ax[1].grid(alpha=0.3); ax[1].legend(fontsize=8, frameon=False); ax[1].set_title("Chewing (t_r ~ 100 ms, t_r/T1 ~ 100) is quasi-static", fontsize=10)
    fig.tight_layout(); fig.savefig("figures/fig10_impact.png"); plt.close(fig)

# ---------------------------------------------------------------- chewing animation
C = np.load("results/chew_frames.npz"); frames, tf = C["frames"], C["t"]
gC = grid(m, frames[0].astype(float), scale=0.0)
P0 = gC.points.copy()
def surf_of(tags):
    sg = sub(gC, tags)
    sf = sg.extract_surface(algorithm="dataset_surface")
    pid = np.asarray(sg.point_data["vtkOriginalPointIds"])[np.asarray(sf.point_data["vtkOriginalPointIds"])]
    cid = np.asarray(sg.cell_data["vtkOriginalCellIds"])[np.asarray(sf.cell_data["vtkOriginalCellIds"])]
    return sf, pid, cid
bone, bid, _ = surf_of([1])
impl, iid, icell = surf_of([2, 3])
ND = m.basis.nodal_dofs
p = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1600, 720), border=False); p.set_background("white")
p.subplot(0, 0)
bone.point_data["|u|"] = np.zeros(bone.n_points)
p.add_mesh(bone, scalars="|u|", cmap="viridis", clim=(0, 45),
           scalar_bar_args=dict(title="|u| (micrometre)", vertical=False, position_x=0.2, width=0.6))
p.add_mesh(impl.copy(), color=TI)
view_oblique(p, m, 1.2)
txt = p.add_text("", font_size=11)
p.subplot(0, 1)
impl.cell_data["vM"] = np.zeros(impl.n_cells)
p.add_mesh(impl, scalars="vM", cmap="turbo", clim=(0, 110),
           scalar_bar_args=dict(title="von Mises (MPa)", vertical=False, position_x=0.2, width=0.6))
p.add_mesh(bone.copy(), color=BONE, opacity=0.3)
view_lateral_left(p, m, zoom=2.3)
SC = 300.0
w = imageio.get_writer("animations/chewing_cycle.mp4", fps=12, codec="libx264", quality=8, macro_block_size=8)
gw = imageio.get_writer("animations/chewing_cycle.gif", fps=12, loop=0)
for k in range(len(tf)):
    u = frames[k].astype(float)
    U = u[ND].T
    mag = np.linalg.norm(U, axis=1) * 1e3
    bone.points = P0[bid] + SC * U[bid]; bone.point_data["|u|"] = mag[bid]
    impl.points = P0[iid] + SC * U[iid]; impl.cell_data["vM"] = vm_all(u)[icell]
    p.subplot(0, 0)
    txt.SetText(2, f"Chewing cycle, left molar  t = {tf[k]*1e3:5.0f} ms   (deformation x{SC:.0f})")
    p.render(); img = p.screenshot(return_img=True)
    w.append_data(img); gw.append_data(img[::2, ::2])
w.close(); gw.close(); p.close()

# ---------------------------------------------------------------- static load-up animation (base mesh)
print("done")
