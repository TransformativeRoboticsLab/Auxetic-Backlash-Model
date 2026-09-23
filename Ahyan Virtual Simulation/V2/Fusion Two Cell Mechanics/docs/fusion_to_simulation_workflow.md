# Fusion To Simulation Workflow

## Objective

Create a validated bridge from the real two-connected-cell Fusion model to the
browser and Python RAD simulator. The first physics target is the uniform
non-backlash pair. The simulator should reproduce the geometry and motion of two
connected cells before backlash, locks, vertical pin-hole clearance, sheet
coupling, inverse design, or neural surrogates are layered on top.

## Core Idea

Fusion should remain the source of truth for the CAD mechanism. Codex should not
guess the assembly topology from screenshots or a raw mesh. The CAD model should
be exposed as structured data:

```text
Fusion model -> Fusion API or MCP -> design graph + pose samples + geometry -> simulator
```

The structured export must preserve:

- component and body identity;
- parent-child occurrence hierarchy;
- user parameters and their units;
- joint names, types, axes, limits, and connected occurrences;
- body transforms at sampled actuation states;
- pin and hole locations where possible;
- exported separated geometry for visualization or external physics engines.

## Mathematical Target

For the two-cell, no-backlash mechanism, define generalized coordinates

```text
q = [u]
```

where `u` is the physical actuation coordinate or the CAD parameter that drives
the pair. If the two cells have independent commands, use

```text
q = [u_left, u_right].
```

The CAD-derived kinematics should produce:

```text
T_i(q)      = rigid transform of body i
p_j(q)      = tracked point or marker j
c_k(q) = 0  = exact joint/closure constraints
d_ab(q)     = clearance or contact distance between pin/hole primitives
```

The first simulator correction should fit and validate a reduced map:

```text
alpha_left(q), alpha_right(q)
theta_left(q), theta_right(q)
center_left(q), center_right(q)
pitch(q)
```

For a uniform, non-backlash two-cell pair, the simulated center spacing should be
derived from the CAD pose samples, not hand-tuned constants. Backlash should be
added only after the no-backlash baseline is reproducible.

## Recommended Workflow

1. Prepare the Fusion file.
   - Use the exact two-connected-cell mechanism.
   - Remove backlash-specific lock features for the baseline if possible.
   - Confirm the model has named components, bodies, joints, and parameters.
   - Identify the CAD parameter that corresponds to alpha, theta, servo angle,
     or linear actuation.

2. Export a design graph.
   - Preferred: query live Fusion through MCP/API if available.
   - Fallback: run `fusion_scripts/export_two_cell_design_graph.py` inside
     Fusion's Scripts and Add-Ins environment.
   - Export JSON and STEP into `exports/`.

3. Export pose samples.
   - Sweep the actuation parameter through the safe mechanical range.
   - At each sample, regenerate the model.
   - Record all occurrence transforms, tracked point coordinates, joint values,
     and any collision/interference status.
   - Use at least 21 samples across the range for the first calibration.

4. Identify mechanism variables.
   - Define left and right cell centers.
   - Define the cell angle used in the simulator.
   - Define the pitch between the connected cells.
   - Define plate/link rotations relative to the cell center.
   - Define which points should remain coincident because of joints.

5. Fit reduced kinematics.
   - Fit smooth maps from actuation to center spacing and body rotations.
   - Prefer analytic trigonometric or linkage equations where the CAD topology
     makes them clear.
   - Use interpolation only as a temporary bridge when the joint structure is
     not yet clean enough for first-principles equations.

6. Validate the pair.
   - Compare simulator and CAD positions at held-out actuation values.
   - Report RMS and maximum coordinate error.
   - Check that the pair remains symmetric when driven symmetrically.
   - Check that the shared connector geometry stays closed.
   - Check that the no-backlash pair has no artificial dead zone.

7. Integrate into the simulator.
   - Add a two-cell calibration profile.
   - Replace current hand-tuned pair expansion constants with CAD-derived
     functions.
   - Keep the older normalized model selectable as "abstract mode."
   - Make all CAD-derived parameters visibly labeled as unvalidated until bench
     retroreflector data confirms them.

8. Add later layers only after the baseline passes.
   - Backlash dead-zone operator.
   - Discrete locks and event locks.
   - Vertical pin-hole clearance and residual neighbor movement.
   - Multi-cell graph propagation.
   - Rigid-body contact in MuJoCo, Isaac Sim, Gazebo, or another external engine.
   - Neural surrogate trained on retroreflector data.

## Questions To Resolve Before Calibration

- What is the exact Fusion file for the two-connected-cell, no-backlash pair?
- What parameter or joint value is the real actuation input?
- What is the physically meaningful alpha definition in the CAD model?
- What is the physically meaningful theta definition in the CAD model?
- Which bodies form one cell, and which bodies are shared connectors?
- Are the joints declared explicitly in Fusion, or are they implied by geometry?
- What coordinate frame should become the simulator frame?
- Which points should correspond to retroreflector marker locations?
- What actuation range is safe for the real mechanism?
- Which two-cell behavior should be treated as ground truth if CAD and bench data
  disagree?
