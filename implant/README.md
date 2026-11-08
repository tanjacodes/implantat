# Project 30 – 3-D mechanical analysis of a mandibular-angle fracture implant

Simulation models, results and COMSOL files that accompany the written report.

## Contents

| Folder / file | What it is |
|---|---|
| `geometry.py` | Parametric mandible + fracture + 3-D strut plate + screws, meshed with Gmsh (quadratic tets). All design variables are in `DEFAULT`. |
| `fe.py` | Linear-elastic FE solver (scikit-fem, P2 tets, PARDISO): materials, muscle loads, supports, stress/strain/IFM post-processing. |
| `run_static.py` | Three static clenching tasks for the baseline plate → `results/static_base.json`, `vtu/base_*.vtu`. |
| `conv.py` | Mesh-convergence study → `results/conv_*.json`. |
| `run_dynamic.py` | Modal analysis, transient chewing cycle, sudden-bite cases → `results/dynamic.json` (`python3 run_dynamic.py impact` reruns only the sudden-bite cases). |
| `run_optim.py`, `run_combo.py` | Optimisation sweeps and combined candidates → `results/optim_*.json`. |
| `run_final.py` | Optimised plate, all tasks, base mesh → `results/static_opt.json`, `vtu/opt_*.vtu`. |
| `fig_*.py`, `vis.py` | Figures and animations (PyVista, Matplotlib). |
| `base.msh`, `opt.msh` | Meshes of the baseline and optimised designs (Gmsh 4.1). |
| `vtu/` | Full-field results for ParaView: displacement `u`, element von Mises, ε1, ε3, material id. Use *Warp By Vector* on `u` (×300–500) to see deformation. |
| `figures/`, `animations/` | All report figures; chewing-cycle MP4/GIF and mode-1 GIF. |
| `comsol/MandibleImplant.java` | COMSOL 6.x Java API model (geometry, materials, physics, mesh, 5 studies, plots). |
| `comsol/COMSOL_build_guide.md` | Click-by-click GUI build of the same model, with all values. |

## Reproduce

```
pip install gmsh scikit-fem meshio pypardiso pyvista matplotlib imageio imageio-ffmpeg
python3 geometry.py base.msh            # baseline mesh (≈20 s)
python3 run_static.py base.msh base     # ≈2 min, needs ~6 GB RAM
python3 conv.py                         # meshes conv_*.msh; dynamics use conv_medium.msh
python3 run_dynamic.py                  # ≈8 min
python3 run_optim.py && python3 run_combo.py && python3 run_final.py
python3 fig_static.py && python3 fig_dynamic.py && python3 fig_optim.py
```

Change a design: `python3 geometry.py my.msh '{"wb":2.4,"holes_p":[-5.5,5.5]}'`.

## COMSOL

```
comsolcompile comsol/MandibleImplant.java
comsol batch -inputfile comsol/MandibleImplant.class -outputfile mandible_implant.mph
```
then open the .mph and run studies std1–std5. The Java file was syntax-checked but not executed against a COMSOL licence; if any API call differs in your version, follow the GUI guide (same names and order).

Units throughout: mm, N, MPa, s.
