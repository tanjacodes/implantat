# Building the mandible + 3-D plate model in COMSOL (GUI guide)

This guide rebuilds, click by click, the model that `MandibleImplant.java` creates through the API.
Names (`wp1`, `ext3`, `csel2` …) match the script, so you can mix both routes.
Tested target: COMSOL Multiphysics 6.x + Structural Mechanics Module. Units: mm, N, MPa.

## 0. Coordinate convention

The fractured **left** hemimandible is built axis-aligned:
x runs along the body from the symphysis (x = 0) back to the gonion (x = Lb), y is bucco-lingual (+y medial, buccal face at y = −T/2), z is superior.
The **right** hemimandible is a mirror copy about the mid-sagittal plane through the origin with normal (−sin θ, cos θ, 0).

## 1. New model

Model Wizard → 3D → Structural Mechanics → **Solid Mechanics (solid)** → Stationary → Done.
Geometry 1 → Length unit **mm**.

## 2. Parameters (Global Definitions → Parameters → Load from file)

Copy the `model.param().set(...)` lines from the script into a text file (`name value description`) and load it, or type them. Key values:

| Name | Value | Meaning |
|---|---|---|
| theta | 35[deg] | hemibody angle to mid-sagittal plane |
| Lb, H, T | 72, 28, 11 [mm] | body length, height, thickness |
| rw, rh, lean | 32, 58, 20 [mm] | ramus width, height, backward lean |
| rc, rl | 5, 20 [mm] | condyle radius, length |
| c | 2[mm] | cortical thickness |
| fAs, fBs | Lb−30, Lb−12 [mm] | fracture at alveolar crest / inferior border |
| g | 0.6[mm] | fracture gap |
| Lf | sqrt((fAs−fBs)^2+H^2) | fracture length |
| dgap | 0.5·g·Lf/H | horizontal half-width of gap strip |
| alpha | atan2(−(fAs−fBs),H) | bar angle ⟂ fracture (≈ 32.7°) |
| tp, wb, ws, rr | 1.0, 1.6, 1.4, 2.0 [mm] | plate thickness, bar width, strut width, ring radius |
| Q, p1, p2 | 5, 5, 11 [mm] | half bar spacing, inner/outer hole offsets |
| zc, Cx | 14[mm], fBs+zc/H·(fAs−fBs) | plate centre on fracture line |
| sr, sl | 1.0, 6 [mm] | screw radius, length in bone |
| FM, FP, FT | 136, 76.8, 176.6 [N] | max muscle forces per side |
| lc | 2 | load case 1/2/3 = incisal / left molar / right molar |
| aL, aR | if(lc==1,1,if(lc==2,1,0.6)), if(lc==1,1,if(lc==2,0.6,1)) | side scaling |
| kM, kP, kT | if(lc==1,0.40,1), if(lc==1,0.75,1), if(lc==1,0.15,1) | task activation |
| dyn | 0 | 1 = chewing time profile |

**Optimised plate:** set `wb = 2.4[mm]`, delete the two outer columns (holes at ±p2, struts at ±p2) and set the bar length to `2*p1` with `p1 = 5.5[mm]` (4-hole square plate).

## 3. Geometry

Under Geometry 1 create three **Cumulative Selections** (right-click Selections): `csel1` Cancellous core, `csel2` Fracture gap, `csel3` Titanium implant.

1. **Work Plane wp1**: plane type Quick, plane **xz**, y-coordinate **T/2**. In its Plane Geometry add four **Polygons** (Data source: Table, coordinates in (x, z)):
   * pol1 (anterior body): (−6,0) (Lb−26,0) (Lb−26,H) (−6,H)
   * pol2 (angle region): (Lb−26,0) (Lb,0) (Lb+lean·30/rh,30) (Lb−rw+lean·30/rh,30) (Lb−rw+lean·H/rh,H) (Lb−26,H)
   * pol3 (upper ramus): (Lb−rw+lean·30/rh,30) (Lb+lean·30/rh,30) (Lb+lean,rh) (Lb−rw+lean,rh)
   * pol4 (coronoid): (Lb−rw+lean,rh) (Lb−rw+lean+10,rh) (Lb−rw+lean+8,rh+8) (Lb−rw+lean−2,rh+8)
   * **Union** uni1 of all four, *Keep interior boundaries* ON (these splits create the muscle insertion faces).
2. **Extrude ext1** from wp1, distance **T** (the xz work plane normal points to −y, so the solid spans −T/2 … T/2; tick *Reverse direction* if your preview shows otherwise).
3. **Work Plane wp2** (xz, y = **T/2−c**) with one polygon (the profile offset by c):
   (0,c) (Lb−1.06c+lean·c/rh, c) (Lb+lean−1.06c−lean·c/rh, rh−c) (Lb−rw+lean+1.06c−lean·c/rh, rh−c) (Lb−rw+lean·(H−c)/rh+1.06c, H−c) (0, H−c).
   **Extrude ext2**, distance **T−2c**, *Contribute to* `csel1`.
4. **Cylinder cyl1**: r = rc, h = rl, axis type **y**, position (Lb+lean−6, −rl/2, rh+1).
5. **Union uni1** (3D) of ext1, ext2, cyl1, *Keep interior boundaries* ON.
6. **Mirror mir1**: input uni1, *Keep input objects* ON, point (0,0,0), normal (−sin(theta), cos(theta), 0).
7. **Work Plane wp3** (xz, y = T/2), polygon (fBs−dgap,0) (fBs+dgap,0) (fAs+dgap,H) (fAs−dgap,H); **Extrude ext3** distance T, contribute to `csel2`.
8. **Work Plane wp4** (xz, y = **−T/2**) – the plate outline (all in the work-plane (x,z) coordinates):
   * 2 Rectangles (bars): size (2·p2, wb), base Center, position (Cx∓Q·sin(alpha), zc±Q·cos(alpha)), rotation **alpha**.
   * 4 Rectangles (struts): size (ws, 2Q), base Center, position (Cx+p·cos(alpha), zc+p·sin(alpha)) for p = −p2, −p1, p1, p2, rotation alpha.
   * 8 Circles (rings): radius rr at (Cx+p·cos(alpha)−q·sin(alpha), zc+p·sin(alpha)+q·cos(alpha)) for p ∈ {−p2,−p1,p1,p2}, q ∈ {−Q,Q}.
   * Union of all (interior boundaries OFF).
   **Extrude ext4** distance **tp** (outwards, to y = −T/2−tp), contribute to `csel3`.
9. **Work Plane wp5** (xz, y = **−T/2+sl**) with 8 Circles of radius sr at the same hole centres; **Extrude ext5** distance **sl+tp**, contribute to `csel3`.
10. Three **Points**: left molar (38, 0, 28); right molar = mirror of it (38·cos2θ, 38·sin2θ, 28) ≈ (13.0, 35.7, 28); incisal (−cos θ·4.6, −sin θ·4.6, 28) ≈ (−3.77, −2.64, 28).
11. **Form Union** (default) → Build All. All interfaces are bonded (conformal), as in the Python model.

Check: Domains in `csel2` = 1 thin strip (left only). The mirrored right side has no gap.

## 4. Selections (Definitions → Selections → Ball, entity level as stated, condition *Intersects* unless noted)

| Tag | Level | Centre | r | Purpose |
|---|---|---|---|---|
| condL | Boundary | (Lb+lean−6, 0, rh+1+rc) = (86, 0, 64) | 0.5 | left condyle top faces |
| condR | Boundary | mirror of condL ≈ (29.4, 80.8, 64) | 0.5 | right condyle |
| masL | Boundary | (69, −5.5, 3) and (56, −5.5, 1) (two Balls + Union) | 0.3 | masseter, buccal angle faces both sides of the fracture |
| mptL | Boundary | (69, 5.5, 3) and (56, 5.5, 1) | 0.3 | medial pterygoid, lingual faces |
| temL | Boundary, condition **Inside** | (64, 0, 62) | 9 | coronoid faces (temporalis) |
| masR | Boundary | (18.43, 66.72, 3) and (13.98, 54.50, 1) | 0.3 | right masseter |
| mptR | Boundary | (28.77, 62.96, 3) and (24.32, 50.74, 1) | 0.3 | right medial pterygoid |
| temR | Boundary, Inside | (21.89, 60.14, 62) | 9 | right coronoid |
| biteL, biteR, biteI | Point | the three Points from step 10 | 0.2 | occlusal contacts |

The script computes all mirrored coordinates for you (`mir()` helper).

## 5. Materials (order matters – later ones override earlier ones)

| Material | Selection | E | ν | ρ |
|---|---|---|---|---|
| Cortical bone | All domains | 13.7 GPa | 0.30 | 1900 kg/m³ |
| Cancellous bone | csel1 | 1.37 GPa | 0.30 | 700 kg/m³ |
| Fracture gap tissue | csel2 | 3 MPa | 0.40 | 1050 kg/m³ |
| Titanium CP grade 4 | csel3 | 105 GPa | 0.34 | 4510 kg/m³ |

Use property group *Young's modulus and Poisson's ratio*. If the right-side core did not inherit `csel1` through the mirror, add its core domains to the Cancellous material by hand.

## 6. Physics (Solid Mechanics, Linear Elastic Material)

* **Fixed Constraint** fix1 = condL, fix2 = condR (TMJ, pinned).
* **Prescribed Displacement (point)** disp1 = biteI, disp2 = biteL, disp3 = biteR: only *z* prescribed, value 0.
* **Boundary Load**, Load type **Total force**, one per muscle and side (6 features). Force vector = magnitude × activation × unit direction in COMSOL axes:

| Muscle | Left side Ftot (x, y, z) | Right side Ftot (x, y, z) |
|---|---|---|
| Masseter | FM·kM·aL·k(t)·(−0.2222, −0.4103, 0.8844) | FM·kM·aR·k(t)·(−0.4616, −0.0685, 0.8844) |
| Medial pterygoid | FP·kP·aL·k(t)·(−0.6074, 0.2429, 0.7563) | FP·kP·aR·k(t)·(0.0205, −0.6539, 0.7563) |
| Temporalis | FT·kT·aL·k(t)·(0.5392, 0.1073, 0.8353) | FT·kT·aR·k(t)·(0.2852, 0.4700, 0.8353) |

with k(t) = `(dyn*gchew(t)+(1-dyn))`. (Directions in anatomical axes (lateral, anterior, superior): masseter (0.21, 0.42, 0.89), medial pterygoid (−0.55, 0.36, 0.76), temporalis (0.22, −0.50, 0.83), normalised and rotated into the COMSOL frame.)

* **Damping** (under Linear Elastic Material): Rayleigh, α_dM = 427 1/s, β_dK = 5.3e−6 s (ζ = 5 % at 1041 Hz and 1956 Hz).
* **Analytic function** gchew(t), argument t [s]: `if(t<0.3,0.1*sin(pi*t/0.3)^2,if(t<0.55,sin(pi*(t-0.3)/0.25)^2,0))`.

## 7. Mesh

Free Tetrahedral, element order quadratic (default). Size: maximum 3.2 mm, minimum 0.2 mm, growth 1.3.
Add Size nodes: `csel3` (implant) max **0.45 mm**, `csel2` (gap) max **0.6 mm**. Expect ~100 000 elements / ~450 000 DOF.

## 8. Studies

| Study | Step | Settings |
|---|---|---|
| std1 Incisal | Stationary | Auxiliary sweep lc = 1; *Modify model configuration*: disable disp2, disp3 |
| std2 Left molar | Stationary | lc = 2; disable disp1, disp3 |
| std3 Right molar | Stationary | lc = 3; disable disp1, disp2 |
| std4 Eigenfrequency | Eigenfrequency | 6 modes, search around 0 Hz; disable disp1, disp3 |
| std5 Chewing | Time Dependent | times range(0,0.005,0.8) s; sweep lc = 2, dyn = 1; disable disp1, disp3 |

For the design load (2× physiological) multiply FM, FP, FT by 2 or simply scale the results (linear model).
For the sudden-bite case use a second analytic function `min(t/0.25e-3,1)` (t in s) with times range(0,2.5e-5,8e-3).

## 9. Results to reproduce

* Volume plot `solid.mises` [MPa] on `csel3` (implant) – expected peak ≈ 110 MPa (left molar, baseline); 87 MPa (optimised).
* Surface plot `solid.disp` [µm] with Deformation ×400 – expected max ≈ 41 µm.
* Surface plot `solid.ep1*1e6` on bone (principal strain, µε).
* Interfragmentary motion: Derived Values → Point Evaluation or a **Cut Line** across the gap, evaluate `solid.disp` difference between the two gap faces (expected ≤ 11 µm).
* Reaction force at the bite point: Derived Values → Point Integration of `solid.RFz` on biteL (≈ −117 N).
* Animation: Results → Export → Animation (Player) on the std5 dataset, 'Sweep type: Solution', for the von Mises plot.

Small differences from the Python results (a few %) are expected from different meshes and from COMSOL's
geometric cortical shell (the Python model assigns cortex by distance to the surface).
