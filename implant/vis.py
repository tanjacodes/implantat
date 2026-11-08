"""Visualisation helpers (pyvista off-screen) for the mandible / implant model."""
import numpy as np
import pyvista as pv
pv.OFF_SCREEN = True
pv.global_theme.font.family = "arial"
pv.global_theme.font.color = "black"

BONE = "#e9dcc4"; TI = "#8d99a6"


def grid(model, u=None, cell=None, scale=0.0):
    m = model.mesh
    cells = np.hstack([np.full((m.nelements, 1), 4), m.t.T]).ravel()
    pts = m.p.T.copy()
    g = pv.UnstructuredGrid(cells, np.full(m.nelements, pv.CellType.TETRA, dtype=np.uint8), pts)
    g.cell_data["tag"] = model.tag
    g.cell_data["mat"] = np.where(model.tag == 1, np.where(model.cortical, 1, 0), model.tag)
    if u is not None:
        U = u[model.basis.nodal_dofs].T
        g.point_data["u"] = U
        g.point_data["|u| (um)"] = np.linalg.norm(U, axis=1) * 1e3
        if scale:
            g.points = pts + scale * U
    for k, v in (cell or {}).items():
        g.cell_data[k] = v
    return g


def sub(g, tags):
    return g.extract_cells(np.where(np.isin(g.cell_data["tag"], tags))[0])


def view_lateral_left(p, model, zoom=1.0, focus=None):
    """camera looking at the left (lateral) side of the angle"""
    P = model.P
    from geometry import frames, plate_layout
    R = frames(P)
    C, a, _ = plate_layout(P)
    c = R @ np.array([C[0], 0, C[1]]) if focus is None else np.array(focus)
    n = R @ np.array([0, -1.0, 0])         # lateral outward normal (left)
    p.camera.focal_point = c
    p.camera.position = c + 120 * n + np.array([0, 25, 25])
    p.camera.up = (0, 0, 1)
    p.camera.zoom(zoom)


def view_oblique(p, model, zoom=1.0):
    p.camera.focal_point = (0, -35, 28)
    p.camera.position = (-210, 90, 150)
    p.camera.up = (0, 0, 1)
    p.camera.zoom(zoom)
