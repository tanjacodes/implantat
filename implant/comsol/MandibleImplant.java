/*
 * MandibleImplant.java  -  COMSOL Multiphysics model of a fractured mandible
 * fixed with a 3-D (strut / box) titanium plate.  Project 30.
 *
 * Builds the same simplified, parametric model that was solved with the
 * open-source FE code (geometry.py / fe.py): geometry, materials, Solid
 * Mechanics physics, mesh, three static clenching studies, an eigenfrequency
 * study, a transient chewing-cycle study and result plots.
 *
 * Requires: COMSOL Multiphysics 6.x with the Structural Mechanics Module
 * (only basic Solid Mechanics features are used).
 *
 * Compile and run (Windows: use comsolcompile.exe / comsolbatch.exe):
 *     comsolcompile MandibleImplant.java
 *     comsol batch -inputfile MandibleImplant.class -outputfile mandible_implant.mph
 * or open the .class from the COMSOL Desktop: File > Open > MandibleImplant.class
 *
 * COORDINATE SYSTEM USED IN COMSOL
 *   The fractured LEFT hemimandible is built axis-aligned:
 *     x = along the body, from the symphysis (x=0) backwards to the gonion (x=Lb)
 *     y = bucco-lingual, +y medial; the buccal (lateral) face is y = -T/2
 *     z = superior
 *   The RIGHT hemimandible is a mirror copy about the mid-sagittal plane,
 *   which passes through the origin with normal (-sin(theta), cos(theta), 0).
 *   (The Python model uses a rotated "anatomical" frame; results are identical.)
 *
 * Units: mm, N, MPa.
 *
 * NOTE: this file was written against the COMSOL 6.x Java API but could not be
 * executed where it was written (no COMSOL licence).  If a call fails in your
 * version, the companion guide COMSOL_build_guide.md reproduces every step in
 * the GUI, in the same order and with the same names.
 */
import com.comsol.model.*;
import com.comsol.model.util.*;

public class MandibleImplant {

    // ---------- numeric copies of the key parameters (used for point/ball selections)
    static final double TH = Math.toRadians(35.0), LB = 72, H = 28, T = 11, RW = 32, RH = 58,
            LEAN = 20, RC = 5, RL = 20, MOLAR = 38, YINC = 4.6;

    /** mirror a local point (x,y,z) about the mid-sagittal plane */
    static double[] mir(double x, double y, double z) {
        double nx = -Math.sin(TH), ny = Math.cos(TH);
        double k = x * nx + y * ny;
        return new double[]{x - 2 * k * nx, y - 2 * k * ny, z};
    }

    static String[] s(double[] p) {
        return new String[]{String.valueOf(p[0]), String.valueOf(p[1]), String.valueOf(p[2])};
    }

    /** ball selection of boundaries touching a point */
    static void ball(Model model, String tag, String label, double[] p, double r, int dim, String cond) {
        model.component("comp1").selection().create(tag, "Ball");
        model.component("comp1").selection(tag).label(label);
        model.component("comp1").selection(tag).set("entitydim", dim);
        model.component("comp1").selection(tag).set("posx", String.valueOf(p[0]));
        model.component("comp1").selection(tag).set("posy", String.valueOf(p[1]));
        model.component("comp1").selection(tag).set("posz", String.valueOf(p[2]));
        model.component("comp1").selection(tag).set("r", String.valueOf(r));
        model.component("comp1").selection(tag).set("condition", cond);
    }

    /** muscle direction (lat, ant, sup) in anatomical axes -> COMSOL local axes */
    static double[] dirLocal(double lat, double ant, double sup, boolean left) {
        double n = Math.sqrt(lat * lat + ant * ant + sup * sup);
        double gx = (left ? -lat : lat) / n, gy = ant / n, gz = sup / n;
        // global = R * local,  R = [[-sin, cos],[-cos, -sin]]  ->  local = R^T * global
        double lx = -Math.sin(TH) * gx - Math.cos(TH) * gy;
        double ly = Math.cos(TH) * gx - Math.sin(TH) * gy;
        return new double[]{lx, ly, gz};   // same coordinate frame for both sides
    }

    public static Model run() {
        Model model = ModelUtil.create("Model");
        model.label("mandible_angle_fracture_3D_plate.mph");

        // =====================================================================
        // 1. PARAMETERS
        // =====================================================================
        model.param().set("theta", "35[deg]", "Angle of each hemibody to the mid-sagittal plane");
        model.param().set("Lb", "72[mm]", "Body length symphysis-gonion");
        model.param().set("H", "28[mm]", "Body height");
        model.param().set("T", "11[mm]", "Body/ramus thickness");
        model.param().set("rw", "32[mm]", "Ramus width at base");
        model.param().set("rh", "58[mm]", "Ramus height");
        model.param().set("lean", "20[mm]", "Backward lean of ramus");
        model.param().set("rc", "5[mm]", "Condyle radius");
        model.param().set("rl", "20[mm]", "Condyle medio-lateral length");
        model.param().set("c", "2[mm]", "Cortical thickness");
        model.param().set("fAs", "Lb-30[mm]", "Fracture: x at alveolar crest");
        model.param().set("fBs", "Lb-12[mm]", "Fracture: x at inferior border");
        model.param().set("g", "0.6[mm]", "Fracture gap width");
        model.param().set("Lf", "sqrt((fAs-fBs)^2+H^2)", "Fracture length");
        model.param().set("dgap", "0.5*g*Lf/H", "Horizontal half-width of gap strip");
        model.param().set("alpha", "atan2(-(fAs-fBs),H)", "Plate bar angle (perpendicular to fracture)");
        // plate (optimised design: tp=1.0, wb=2.4, 4-hole -> change hole list below; baseline wb=1.6, 8-hole)
        model.param().set("tp", "1.0[mm]", "Plate thickness");
        model.param().set("wb", "1.6[mm]", "Bar width");
        model.param().set("ws", "1.4[mm]", "Strut width");
        model.param().set("rr", "2.0[mm]", "Ring outer radius");
        model.param().set("Q", "5[mm]", "Half distance between bars");
        model.param().set("p1", "5[mm]", "Inner hole distance from fracture");
        model.param().set("p2", "11[mm]", "Outer hole distance from fracture");
        model.param().set("zc", "14[mm]", "Plate centre height");
        model.param().set("Cx", "fBs+zc/H*(fAs-fBs)", "Plate centre x");
        model.param().set("sr", "1.0[mm]", "Screw radius (2.0 mm screws)");
        model.param().set("sl", "6[mm]", "Screw length in bone (monocortical)");
        // loads
        model.param().set("FM", "136[N]", "Masseter max force (per side)");
        model.param().set("FP", "76.8[N]", "Medial pterygoid max force");
        model.param().set("FT", "176.6[N]", "Temporalis max force");
        model.param().set("lc", "2", "Load case: 1 incisal, 2 left molar, 3 right molar");
        model.param().set("aL", "if(lc==1,1,if(lc==2,1,0.6))", "Activation scale left side");
        model.param().set("aR", "if(lc==1,1,if(lc==2,0.6,1))", "Activation scale right side");
        model.param().set("kM", "if(lc==1,0.40,1)", "Masseter activation factor");
        model.param().set("kP", "if(lc==1,0.75,1)", "Medial pterygoid activation factor");
        model.param().set("kT", "if(lc==1,0.15,1)", "Temporalis activation factor");
        model.param().set("dyn", "0", "1 = time-dependent chewing load");

        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 3);
        GeomSequence geom = model.component("comp1").geom("geom1");
        geom.lengthUnit("mm");
        geom.selection().create("csel1", "CumulativeSelection"); geom.selection("csel1").label("Cancellous core");
        geom.selection().create("csel2", "CumulativeSelection"); geom.selection("csel2").label("Fracture gap");
        geom.selection().create("csel3", "CumulativeSelection"); geom.selection("csel3").label("Titanium implant");

        // =====================================================================
        // 2. GEOMETRY
        // =====================================================================
        // 2a. lateral profile of the hemimandible (work plane xz at y=+T/2, extruded T towards -y)
        geom.create("wp1", "WorkPlane");
        geom.feature("wp1").set("quickplane", "xz");
        geom.feature("wp1").set("quicky", "T/2");
        GeomSequence w1 = geom.feature("wp1").geom();
        // P1 anterior body, P2 angle region (masseter / medial pterygoid insertion), P3 upper ramus, P4 coronoid
        w1.create("pol1", "Polygon");
        w1.feature("pol1").set("source", "table");
        w1.feature("pol1").set("table", new String[][]{{"-6", "0"}, {"Lb-26", "0"}, {"Lb-26", "H"}, {"-6", "H"}});
        w1.create("pol2", "Polygon");
        w1.feature("pol2").set("source", "table");
        w1.feature("pol2").set("table", new String[][]{{"Lb-26", "0"}, {"Lb", "0"}, {"Lb+lean*30/rh", "30"},
                {"Lb-rw+lean*30/rh", "30"}, {"Lb-rw+lean*H/rh", "H"}, {"Lb-26", "H"}});
        w1.create("pol3", "Polygon");
        w1.feature("pol3").set("source", "table");
        w1.feature("pol3").set("table", new String[][]{{"Lb-rw+lean*30/rh", "30"}, {"Lb+lean*30/rh", "30"},
                {"Lb+lean", "rh"}, {"Lb-rw+lean", "rh"}});
        w1.create("pol4", "Polygon");
        w1.feature("pol4").set("source", "table");
        w1.feature("pol4").set("table", new String[][]{{"Lb-rw+lean", "rh"}, {"Lb-rw+lean+10", "rh"},
                {"Lb-rw+lean+8", "rh+8"}, {"Lb-rw+lean-2", "rh+8"}});
        w1.create("uni1", "Union");
        w1.feature("uni1").selection("input").set("pol1", "pol2", "pol3", "pol4");
        w1.feature("uni1").set("intbnd", true);
        geom.create("ext1", "Extrude");
        geom.feature("ext1").set("workplane", "wp1");
        geom.feature("ext1").selection("input").set("wp1");
        geom.feature("ext1").setIndex("distance", "T", 0);

        // 2b. cancellous core (profile offset by the cortical thickness c)
        geom.create("wp2", "WorkPlane");
        geom.feature("wp2").set("quickplane", "xz");
        geom.feature("wp2").set("quicky", "T/2-c");
        GeomSequence w2 = geom.feature("wp2").geom();
        w2.create("pol1", "Polygon");
        w2.feature("pol1").set("source", "table");
        w2.feature("pol1").set("table", new String[][]{{"0", "c"}, {"Lb-1.06*c+lean*c/rh", "c"},
                {"Lb+lean-1.06*c-lean*c/rh", "rh-c"}, {"Lb-rw+lean+1.06*c-lean*c/rh", "rh-c"},
                {"Lb-rw+lean*(H-c)/rh+1.06*c", "H-c"}, {"0", "H-c"}});
        geom.create("ext2", "Extrude");
        geom.feature("ext2").set("workplane", "wp2");
        geom.feature("ext2").selection("input").set("wp2");
        geom.feature("ext2").setIndex("distance", "T-2*c", 0);
        geom.feature("ext2").set("contributeto", "csel1");

        // 2c. condyle
        geom.create("cyl1", "Cylinder");
        geom.feature("cyl1").set("r", "rc");
        geom.feature("cyl1").set("h", "rl");
        geom.feature("cyl1").set("axistype", "y");
        geom.feature("cyl1").set("pos", new String[]{"Lb+lean-6", "-rl/2", "rh+1"});

        geom.create("uni1", "Union");
        geom.feature("uni1").selection("input").set("ext1", "ext2", "cyl1");
        geom.feature("uni1").set("intbnd", true);
        geom.feature("uni1").label("Left hemimandible");

        // 2d. right hemimandible = mirror copy about the mid-sagittal plane
        geom.create("mir1", "Mirror");
        geom.feature("mir1").selection("input").set("uni1");
        geom.feature("mir1").set("keep", true);
        geom.feature("mir1").set("pos", new String[]{"0", "0", "0"});
        geom.feature("mir1").set("axis", new String[]{"-sin(theta)", "cos(theta)", "0"});
        geom.feature("mir1").label("Right hemimandible");

        // 2e. fracture gap (left angle) - strip from alveolar crest to inferior border
        geom.create("wp3", "WorkPlane");
        geom.feature("wp3").set("quickplane", "xz");
        geom.feature("wp3").set("quicky", "T/2");
        GeomSequence w3 = geom.feature("wp3").geom();
        w3.create("pol1", "Polygon");
        w3.feature("pol1").set("source", "table");
        w3.feature("pol1").set("table", new String[][]{{"fBs-dgap", "0"}, {"fBs+dgap", "0"},
                {"fAs+dgap", "H"}, {"fAs-dgap", "H"}});
        geom.create("ext3", "Extrude");
        geom.feature("ext3").set("workplane", "wp3");
        geom.feature("ext3").selection("input").set("wp3");
        geom.feature("ext3").setIndex("distance", "T", 0);
        geom.feature("ext3").set("contributeto", "csel2");
        geom.feature("ext3").label("Fracture gap tissue");

        // 2f. 3-D strut plate on the buccal face (work plane at y=-T/2, extruded tp outwards)
        geom.create("wp4", "WorkPlane");
        geom.feature("wp4").set("quickplane", "xz");
        geom.feature("wp4").set("quicky", "-T/2");
        GeomSequence w4 = geom.feature("wp4").geom();
        String[] qs = {"-Q", "Q"};
        String[] ps = {"-p2", "-p1", "p1", "p2"};
        int k = 0;
        for (String q : qs) {             // two bars
            k++;
            w4.create("bar" + k, "Rectangle");
            w4.feature("bar" + k).set("base", "center");
            w4.feature("bar" + k).set("size", new String[]{"2*p2", "wb"});
            w4.feature("bar" + k).set("pos", new String[]{"Cx-(" + q + ")*sin(alpha)", "zc+(" + q + ")*cos(alpha)"});
            w4.feature("bar" + k).set("rot", "alpha");
        }
        k = 0;
        for (String p : ps) {             // struts
            k++;
            w4.create("str" + k, "Rectangle");
            w4.feature("str" + k).set("base", "center");
            w4.feature("str" + k).set("size", new String[]{"ws", "2*Q"});
            w4.feature("str" + k).set("pos", new String[]{"Cx+(" + p + ")*cos(alpha)", "zc+(" + p + ")*sin(alpha)"});
            w4.feature("str" + k).set("rot", "alpha");
        }
        k = 0;
        String[] holeX = new String[8], holeZ = new String[8];
        for (String p : ps) for (String q : qs) {   // rings
            holeX[k] = "Cx+(" + p + ")*cos(alpha)-(" + q + ")*sin(alpha)";
            holeZ[k] = "zc+(" + p + ")*sin(alpha)+(" + q + ")*cos(alpha)";
            k++;
            w4.create("ring" + k, "Circle");
            w4.feature("ring" + k).set("r", "rr");
            w4.feature("ring" + k).set("pos", new String[]{holeX[k - 1], holeZ[k - 1]});
        }
        w4.create("uni1", "Union");
        w4.feature("uni1").selection("input").set("bar1", "bar2", "str1", "str2", "str3", "str4",
                "ring1", "ring2", "ring3", "ring4", "ring5", "ring6", "ring7", "ring8");
        w4.feature("uni1").set("intbnd", false);
        geom.create("ext4", "Extrude");
        geom.feature("ext4").set("workplane", "wp4");
        geom.feature("ext4").selection("input").set("wp4");
        geom.feature("ext4").setIndex("distance", "tp", 0);
        geom.feature("ext4").set("contributeto", "csel3");
        geom.feature("ext4").label("3-D strut plate");

        // 2g. screws (cylinders from the plate surface sl mm into the bone)
        geom.create("wp5", "WorkPlane");
        geom.feature("wp5").set("quickplane", "xz");
        geom.feature("wp5").set("quicky", "-T/2+sl");
        GeomSequence w5 = geom.feature("wp5").geom();
        for (int i = 0; i < 8; i++) {
            w5.create("sc" + (i + 1), "Circle");
            w5.feature("sc" + (i + 1)).set("r", "sr");
            w5.feature("sc" + (i + 1)).set("pos", new String[]{holeX[i], holeZ[i]});
        }
        geom.create("ext5", "Extrude");
        geom.feature("ext5").set("workplane", "wp5");
        geom.feature("ext5").selection("input").set("wp5");
        geom.feature("ext5").setIndex("distance", "sl+tp", 0);
        geom.feature("ext5").set("contributeto", "csel3");
        geom.feature("ext5").label("Screws");

        // 2h. occlusal contact points (embedded vertices)
        double[] pL = {MOLAR, 0, H};
        double[] pR = mir(MOLAR, 0, H);
        double[] pI = {-Math.cos(TH) * YINC, -Math.sin(TH) * YINC, H};
        geom.create("ptL", "Point"); geom.feature("ptL").setIndex("p", String.valueOf(pL[0]), 0, 0);
        geom.feature("ptL").setIndex("p", String.valueOf(pL[1]), 1, 0); geom.feature("ptL").setIndex("p", String.valueOf(pL[2]), 2, 0);
        geom.create("ptR", "Point"); geom.feature("ptR").setIndex("p", String.valueOf(pR[0]), 0, 0);
        geom.feature("ptR").setIndex("p", String.valueOf(pR[1]), 1, 0); geom.feature("ptR").setIndex("p", String.valueOf(pR[2]), 2, 0);
        geom.create("ptI", "Point"); geom.feature("ptI").setIndex("p", String.valueOf(pI[0]), 0, 0);
        geom.feature("ptI").setIndex("p", String.valueOf(pI[1]), 1, 0); geom.feature("ptI").setIndex("p", String.valueOf(pI[2]), 2, 0);

        geom.feature("fin").set("action", "union");   // Form Union: bonded (conformal) interfaces everywhere
        geom.run();

        // =====================================================================
        // 3. SELECTIONS (ball selections at known points)
        // =====================================================================
        double sc = LB + LEAN - 6, zcd = RH + 1 + RC;
        ball(model, "condL", "Condyle L (top)", new double[]{sc, 0, zcd}, 0.5, 2, "intersects");
        ball(model, "condR", "Condyle R (top)", mir(sc, 0, zcd), 0.5, 2, "intersects");
        // masseter: buccal face of the angle region, both sides of the fracture
        ball(model, "masL1", "Masseter L post", new double[]{LB - 3, -T / 2, 3}, 0.3, 2, "intersects");
        ball(model, "masL2", "Masseter L ant", new double[]{LB - 16, -T / 2, 1}, 0.3, 2, "intersects");
        ball(model, "masR1", "Masseter R post", mir(LB - 3, -T / 2, 3), 0.3, 2, "intersects");
        ball(model, "masR2", "Masseter R ant", mir(LB - 16, -T / 2, 1), 0.3, 2, "intersects");
        ball(model, "mptL1", "Med. pterygoid L post", new double[]{LB - 3, T / 2, 3}, 0.3, 2, "intersects");
        ball(model, "mptL2", "Med. pterygoid L ant", new double[]{LB - 16, T / 2, 1}, 0.3, 2, "intersects");
        ball(model, "mptR1", "Med. pterygoid R post", mir(LB - 3, T / 2, 3), 0.3, 2, "intersects");
        ball(model, "mptR2", "Med. pterygoid R ant", mir(LB - 16, T / 2, 1), 0.3, 2, "intersects");
        double cx = LB - RW + LEAN + 4, cz = RH + 4;
        ball(model, "temL", "Temporalis L (coronoid)", new double[]{cx, 0, cz}, 9.0, 2, "inside");
        ball(model, "temR", "Temporalis R (coronoid)", mir(cx, 0, cz), 9.0, 2, "inside");
        ball(model, "biteL", "Bite point L molar", pL, 0.2, 0, "intersects");
        ball(model, "biteR", "Bite point R molar", pR, 0.2, 0, "intersects");
        ball(model, "biteI", "Bite point incisal", pI, 0.2, 0, "intersects");
        for (String[] u : new String[][]{{"masL", "masL1", "masL2"}, {"masR", "masR1", "masR2"},
                {"mptL", "mptL1", "mptL2"}, {"mptR", "mptR1", "mptR2"}}) {
            model.component("comp1").selection().create(u[0], "Union");
            model.component("comp1").selection(u[0]).set("entitydim", 2);
            model.component("comp1").selection(u[0]).set("input", new String[]{u[1], u[2]});
        }

        // =====================================================================
        // 4. MATERIALS (later materials override earlier ones)
        // =====================================================================
        String[][] mats = {
                // tag, label, E, nu, rho, selection
                {"mat1", "Cortical bone", "13.7[GPa]", "0.30", "1900[kg/m^3]", "all"},
                {"mat2", "Cancellous bone", "1.37[GPa]", "0.30", "700[kg/m^3]", "geom1_csel1_dom"},
                {"mat3", "Fracture gap tissue", "3[MPa]", "0.40", "1050[kg/m^3]", "geom1_csel2_dom"},
                {"mat4", "Titanium CP grade 4", "105[GPa]", "0.34", "4510[kg/m^3]", "geom1_csel3_dom"}};
        for (String[] m : mats) {
            Material mt = model.component("comp1").material().create(m[0], "Common");
            mt.label(m[1]);
            mt.propertyGroup().create("Enu", "Young's modulus and Poisson's ratio");
            mt.propertyGroup("Enu").set("E", m[2]);
            mt.propertyGroup("Enu").set("nu", m[3]);
            mt.propertyGroup("def").set("density", m[4]);
            if (m[5].equals("all")) mt.selection().all(); else mt.selection().named(m[5]);
        }

        // =====================================================================
        // 5. PHYSICS - Solid Mechanics, linear elastic
        // =====================================================================
        Physics solid = model.component("comp1").physics().create("solid", "SolidMechanics", "geom1");
        solid.create("fix1", "Fixed", 2);
        solid.feature("fix1").label("Condyles (TMJ) pinned");
        solid.feature("fix1").selection().set(new int[0]);
        solid.feature("fix1").selection().named("condL");
        solid.create("fix2", "Fixed", 2);
        solid.feature("fix2").selection().named("condR");
        // occlusal contacts: vertical constraint at the bite point (one active per study)
        String[][] bites = {{"disp1", "biteI", "Incisal contact"}, {"disp2", "biteL", "Left molar contact"},
                {"disp3", "biteR", "Right molar contact"}};
        for (String[] b : bites) {
            solid.create(b[0], "Displacement0", 0);
            solid.feature(b[0]).label(b[2]);
            solid.feature(b[0]).selection().named(b[1]);
            solid.feature(b[0]).setIndex("Direction", "prescribed", 2);
            solid.feature(b[0]).setIndex("U0", "0", 2);
        }
        // muscle forces: total force over the insertion area (dyn=1 -> chewing time profile)
        String gfun = "(dyn*gchew(t)+(1-dyn))";
        Object[][] muscles = {
                {"Masseter", "FM", "kM", new double[]{0.21, 0.42, 0.89}, "masL", "masR"},
                {"Medial pterygoid", "FP", "kP", new double[]{-0.55, 0.36, 0.76}, "mptL", "mptR"},
                {"Temporalis", "FT", "kT", new double[]{0.22, -0.50, 0.83}, "temL", "temR"}};
        k = 0;
        for (Object[] mu : muscles) {
            double[] d = (double[]) mu[3];
            for (int side = 0; side < 2; side++) {
                k++;
                boolean left = side == 0;
                double[] v = dirLocal(d[0], d[1], d[2], left);
                String amp = mu[1] + "*" + mu[2] + "*" + (left ? "aL" : "aR") + "*" + gfun;
                String tag = "bndl" + k;
                solid.create(tag, "BoundaryLoad", 2);
                solid.feature(tag).label(mu[0] + (left ? " L" : " R"));
                solid.feature(tag).selection().named((String) mu[left ? 4 : 5]);
                solid.feature(tag).set("LoadType", "TotalForce");
                solid.feature(tag).set("Ftot", new String[]{
                        amp + "*(" + v[0] + ")", amp + "*(" + v[1] + ")", amp + "*(" + v[2] + ")"});
            }
        }
        // Rayleigh damping (5 % at 1st and 4th eigenfrequency of the occluded system)
        solid.feature("lemm1").create("dmp1", "Damping", 3);
        solid.feature("lemm1").feature("dmp1").set("DampingType", "RayleighDamping");
        solid.feature("lemm1").feature("dmp1").set("alpha_dM", "427[1/s]");
        solid.feature("lemm1").feature("dmp1").set("beta_dK", "5.3e-6[s]");

        // chewing-cycle activation (0.8 s: closing 10 %, power stroke 100 %, opening 0)
        model.component("comp1").func().create("an1", "Analytic");
        model.component("comp1").func("an1").set("funcname", "gchew");
        model.component("comp1").func("an1").set("args", new String[]{"t"});
        model.component("comp1").func("an1").set("expr",
                "if(t<0.3,0.1*sin(pi*t/0.3)^2,if(t<0.55,sin(pi*(t-0.3)/0.25)^2,0))");
        model.component("comp1").func("an1").set("argunit", new String[]{"s"});

        // =====================================================================
        // 6. MESH - quadratic tetrahedra, refined at implant and fracture
        // =====================================================================
        Mesh mesh = model.component("comp1").mesh().create("mesh1");
        mesh.create("ftet1", "FreeTet");
        mesh.feature("ftet1").create("size1", "Size");
        mesh.feature("ftet1").feature("size1").selection().geom("geom1", 3);
        mesh.feature("ftet1").feature("size1").selection().named("geom1_csel3_dom");
        mesh.feature("ftet1").feature("size1").set("custom", true);
        mesh.feature("ftet1").feature("size1").set("hmaxactive", true);
        mesh.feature("ftet1").feature("size1").set("hmax", "0.45");
        mesh.feature("ftet1").create("size2", "Size");
        mesh.feature("ftet1").feature("size2").selection().geom("geom1", 3);
        mesh.feature("ftet1").feature("size2").selection().named("geom1_csel2_dom");
        mesh.feature("ftet1").feature("size2").set("custom", true);
        mesh.feature("ftet1").feature("size2").set("hmaxactive", true);
        mesh.feature("ftet1").feature("size2").set("hmax", "0.6");
        mesh.feature("size").set("custom", true);
        mesh.feature("size").set("hmax", "3.2");
        mesh.feature("size").set("hmin", "0.2");
        mesh.feature("size").set("hgrad", "1.3");
        mesh.run();

        // =====================================================================
        // 7. STUDIES
        // =====================================================================
        String[][] st = {{"std1", "Static - incisal clench", "1", "solid/disp2", "solid/disp3"},
                {"std2", "Static - left (fractured side) molar clench", "2", "solid/disp1", "solid/disp3"},
                {"std3", "Static - right molar clench", "3", "solid/disp1", "solid/disp2"}};
        for (String[] x : st) {
            Study sd = model.study().create(x[0]);
            sd.label(x[1]);
            sd.create("stat", "Stationary");
            sd.feature("stat").set("useparam", true);
            sd.feature("stat").set("pname", new String[]{"lc"});
            sd.feature("stat").set("plistarr", new String[]{x[2]});
            sd.feature("stat").set("punit", new String[]{""});
            sd.feature("stat").set("useadvanceddisable", true);
            sd.feature("stat").set("disabledphysics", new String[]{x[3], x[4]});
        }
        // eigenfrequencies, condyles pinned + left molar contact
        Study se = model.study().create("std4");
        se.label("Eigenfrequency (left molar occlusion)");
        se.create("eig", "Eigenfrequency");
        se.feature("eig").set("neigs", 6);
        se.feature("eig").set("shift", "0");
        se.feature("eig").set("useadvanceddisable", true);
        se.feature("eig").set("disabledphysics", new String[]{"solid/disp1", "solid/disp3"});
        // transient chewing cycle, left molar
        Study stt = model.study().create("std5");
        stt.label("Time dependent - chewing cycle (left molar)");
        stt.create("time", "Transient");
        stt.feature("time").set("tlist", "range(0,0.005,0.8)");
        stt.feature("time").set("useparam", true);
        stt.feature("time").set("pname", new String[]{"lc", "dyn"});
        stt.feature("time").set("plistarr", new String[]{"2", "1"});
        stt.feature("time").set("punit", new String[]{"", ""});
        stt.feature("time").set("useadvanceddisable", true);
        stt.feature("time").set("disabledphysics", new String[]{"solid/disp1", "solid/disp3"});

        // =====================================================================
        // 8. RESULTS (plots are created; run the studies, then press Plot)
        // =====================================================================
        String[] dsets = {"dset1", "dset2", "dset3"};
        for (int i = 0; i < 3; i++) {
            String pg = "pgVM" + (i + 1);
            model.result().create(pg, "PlotGroup3D");
            model.result(pg).label("von Mises in implant - " + st[i][1]);
            model.result(pg).set("data", dsets[i]);
            model.result(pg).create("vol1", "Volume");
            model.result(pg).feature("vol1").set("expr", "solid.mises");
            model.result(pg).feature("vol1").set("unit", "MPa");
            model.result(pg).feature("vol1").set("colortable", "Turbo");
            model.result(pg).feature("vol1").create("sel1", "Selection");
            model.result(pg).feature("vol1").feature("sel1").selection().named("geom1_csel3_dom");
            String pu = "pgU" + (i + 1);
            model.result().create(pu, "PlotGroup3D");
            model.result(pu).label("Displacement - " + st[i][1]);
            model.result(pu).set("data", dsets[i]);
            model.result(pu).create("surf1", "Surface");
            model.result(pu).feature("surf1").set("expr", "solid.disp");
            model.result(pu).feature("surf1").set("unit", "um");
            model.result(pu).feature("surf1").create("def1", "Deform");
            model.result(pu).feature("surf1").feature("def1").set("scaleactive", true);
            model.result(pu).feature("surf1").feature("def1").set("scale", "400");
            String pe = "pgE" + (i + 1);
            model.result().create(pe, "PlotGroup3D");
            model.result(pe).label("1st principal strain in bone - " + st[i][1]);
            model.result(pe).set("data", dsets[i]);
            model.result(pe).create("surf1", "Surface");
            model.result(pe).feature("surf1").set("expr", "solid.ep1*1e6");
            model.result(pe).feature("surf1").set("descr", "1st principal strain (microstrain)");
            model.result(pe).feature("surf1").set("rangecoloractive", true);
            model.result(pe).feature("surf1").set("rangecolormax", "400");
            model.result(pe).feature("surf1").set("rangecolormin", "0");
        }
        // interfragmentary motion: max |u_prox - u_dist| needs two boundary probes; simplest is
        // a Surface plot of solid.disp on the gap domain boundaries, or a Line Graph along the gap.

        return model;
    }

    public static void main(String[] args) {
        Model model = run();
        // Uncomment to solve in batch:
        // model.study("std1").run(); model.study("std2").run(); model.study("std3").run();
        // model.study("std4").run(); model.study("std5").run();
    }
}
