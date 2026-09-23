# RAD Research Grounding

This note records what the current simulator can treat as paper-supported versus
what should remain a modeling assumption until measured on hardware. Page
numbers refer to the local PDFs supplied with this project.

## Source Files

- `C:\Users\ahyan\Downloads\3002882.pdf`: RLABS proceedings paper,
  "Reconfigurable Lattice of Auxetic Backlash Structures for Shape-Changing
  Surfaces".
- RADs preprint, located in `C:\Users\ahyan\Downloads`:
  `Reconfigurable_Auxetic_Devices__RADs__for_Robotic_Surface_Manipulation (2).pdf`.
  Title: "Reconfigurable Auxetic Devices (RADs) with Backlash for Soft Robotic
  Surface Manipulation".
- `C:\Users\ahyan\Downloads\programmable_mechanics_monograph.pdf`:
  project monograph, "Programmable Discontinuities and the Foundations of
  Programmable Mechanics".

## Evidence Index

- RLABS proceedings pages 1-4: configuration-space framing, conformal-map
  language, rotating-square dilation, and the `alpha(x, y)` field.
- RLABS proceedings pages 5-7: RoR/RoM definitions, die-off distance, locked and
  free cell behavior, and the bidirectional ReLU backlash model.
- RLABS proceedings page 9: airfoil/conformal molding demonstrations and
  gradient-descent correction of locking angles.
- RAD preprint page 5: RAD cell design equations, normalized backlash,
  dead-zone/ReLU coupling, theta-to-alpha relation, and rubber skin framing.
- RAD preprint page 6: two concentric parts, four joints per part, target
  Poisson ratio, normalized backlash, 35 mm side length, PLA/Prusa fabrication,
  0.1 mm hole tolerance, underside servo placement, and gravity/suspension
  description.
- RAD preprint pages 8-9: 0.1 radian servo calibration, force/contact feedback,
  surface area/conformity measurement, and reported modeled/measured area
  agreement.
- Monograph pages 13, 23, 25, 42, 52, and 55: programmable-discontinuity
  definition, open questions, structured system state, event operators, event
  algebra, validation ladder, and computational-experiment patterns.

## Paper-Supported RAD/RLABS Mechanics

### Rotating-square dilation

- The RLABS paper frames the unit cells as a rotating-square auxetic lattice.
  Its Figure 1 describes dilation across the lattice as arising from the
  rotating-square mechanism and states that this produces angle-preserving
  conformal behavior through a spatial dilation factor `alpha(x, y)`.
- The RAD preprint keeps the same rotating-square basis: the unit cells are
  "rotating square auxetic mechanisms" and local shape control is described as
  variation of the dilation factor `alpha(x, y)` across the lattice.
- The RAD preprint gives the current simulator's angle law:

```text
theta[degrees] = 70 * alpha - 60
```

This is paper-supported for the simplified cell-coupling model, not a universal
law for all future CAD geometries.

### Backlash and dead-zone coupling

- The RAD preprint defines normalized backlash as:

```text
b_norm = b / L
```

- The same preprint defines adjacent-cell coupling with a bidirectional
  dead-zone/ReLU law:

```text
f(x) = max(0, x - b) + min(x + b, 0)
```

- The RLABS paper uses the same dead-zone idea to describe rotational stiffness:
  below the backlash threshold the output is zero; outside the threshold the
  cell transmits motion. It connects this to a finite die-off distance from a
  locked or actuated cell.
- Paper-supported interpretation: backlash is a local discontinuity that
  changes when neighboring cells become mechanically coupled. The current
  simulator's superposition residual, interaction hotspot, interaction degree,
  and event-order metrics are diagnostics built on top of this idea; they are
  not separately paper-derived constitutive laws.

### Cell geometry and prototype values

- The RAD preprint states that the CAD unit cell consists of two concentric
  parts with four joints each.
- It reports a target effective Poisson ratio of `-0.4`, normalized backlash
  `b = 0.1`, and prototype cell side length `35 mm`.
- It reports prototype lattices with `8x11` and `11x15` unit cells.
- Fabrication details in the RAD preprint: Prusa MK4, PLA filament, and hole
  fabrication consistency within `0.1 mm` tolerance.
- The RAD preprint says servos are mounted underneath the lattice, preserving
  the upper surface, reducing mechanical interference, and allowing the soft
  surface to remain unobstructed.
- Prototype servo responses were calibrated at intervals of `0.1 radians`.

### Surface and membrane behavior

- The RAD preprint describes a soft rubber skin placed on the lattice surface,
  acting as a non-permeable interpolating membrane.
- The preprint frames surface adaptation as variable auxetic behavior across the
  manipulation interface via backlash between cells.
- The suspended lattice is described as deforming under its own weight, with
  gravity and boundary conditions setting a neutral curvature. That means the
  current browser membrane is a useful visualization but not a calibrated skin
  mechanics model.

### Conformal-map framing

- The RLABS paper states that conformal mapping preserves local angles while not
  preserving size or curvature.
- In the RLABS framing, the dilation field `alpha(x, y)` is the link between the
  rotating-square lattice state and conformal surface deformation.
- The proceedings paper describes the configuration space as growing with the
  lattice and uses locking/free-cell states to reconfigure the mapping function.
- Current implication: target-surface inverse design should reason over
  `alpha(x, y)`, boundary conditions, locks, and actuator commands. It should
  not assume a single fabricated conformal map.

## Programmable-Discontinuity Framework Support

The monograph is a project blueprint rather than an external mechanics paper,
but it provides a useful formal language for simulator architecture:

- Mechanical system state:

```text
S = (G, Q, q, q_dot, C, Pi, M, sigma, theta)
```

- Binary lock modes:

```text
sigma in {0, 1}^m
C_sigma = C_0 union C_lock(sigma)
```

- Discontinuity event:

```text
D_e: S^- -> S^+
```

In the simulator, this maps naturally to:

- `G`: cell-neighbor graph and inter-cell linkage graph.
- `Q, q`: cell dilation, theta, in-plane center positions, and vertical state.
- `C`: lock constraints, boundary constraints, and actuator constraints.
- `Pi`: kinematic residual or spring-hinge energy.
- `sigma`: lock/free/actuated modes.
- `theta`: geometry, backlash, pin/hole clearance, stiffness, and solver
  parameters.
- `D_e`: local actuation, lock, release, command clearing, and future contact or
  snap-through events.

The monograph's validation ladder supports the current implementation sequence:
single spring, two-bar kinematics, square frame rigidity, rotational hinge,
small rotating-unit auxetic, lock event, event sequence, and path-dependence
tests.

## Current Simulator Assumptions

These assumptions are useful for v1, but should remain explicitly labeled:

- The browser `Paper RAD cell` view is normalized and CAD-like, not an exact
  fabrication model.
- Plate thickness, boss dimensions, exact pin radius, exact hole radius, screw
  geometry, friction, and joint compliance are configurable estimates unless
  measured from the physical parts.
- Vertical residual coupling is a current model extension motivated by
  pin-hole clearance and user-observed hardware behavior. The supplied RAD text
  supports backlash and hole tolerance, but does not provide a calibrated
  vertical die-off law.
- The spring-preview/browser relaxation and Python spring-hinge solver are
  reduced-order physical approximations. They should be compared against scans,
  metrology, or motion tracking before being used as predictive hardware models.
- The MuJoCo export/run/comparison layer is an independent-engine validation
  path, but the current MJCF body and pin-hole geometry remains a coarse
  normalized proxy. The contact-parameter profile makes stiffness, damping,
  friction, and solver assumptions explicit for export, but those parameters
  remain uncalibrated until measured. The contact-parameter calibration packet
  adds fit and holdout measurement rows for those assumptions, but blank
  templates are not bench evidence. The contact-parameter bench-validation
  comparison can reject missing fields and residual mismatches in filled rows;
  it still does not prove a physical contact law. The interval-calibration layer
  can bound the chosen proxy parameters against filled fit/holdout rows, but
  those bounds are regime-specific empirical claims rather than constitutive
  mechanics. Agreement with MuJoCo is not enough for hardware accuracy without
  bench measurements, measured pin/hole dimensions, contact stiffness, friction,
  and assembly tolerance data.
- Pairwise non-additivity, interaction hotspots, interaction degree, decay fits,
  and event-order errors are numerical diagnostics for programmable mechanics.
  They help detect where operator composition is non-additive or
  noncommutative, but they are not paper-derived material laws.
- The 12-cell discrete-lock dataset summarized in
  `docs/lock_dataset_presentation_notes.md` should become the calibration
  source for lock operators. The notes distinguish nominal lock placement from
  realized lock engagement, crown angle from effective lattice theta, and
  consecutive-lock failure where a middle crown can fall out.

## Calibration Gaps To Measure

- Actual pin radius, hole radius, and clearance distribution across printed
  cells.
- Plate and linkage thicknesses, boss diameters, joint stack height, screw or
  rivet details, and assembly offsets.
- Friction and hysteresis at joints under repeated actuation.
- Real alpha-to-theta curve across the physical travel range.
- Vertical free play as a function of pin/hole clearance, load, and neighboring
  cell state.
- Rubber skin stiffness, membrane pretension, and coupling to lattice nodes.
- Gravity sag and boundary-condition sensitivity for suspended lattices.
- Servo command-to-cell displacement calibration beyond the reported
  `0.1 radian` sampling interval.

## Implementation Implications

- Preserve three geometry modes: abstract topology, paper-inspired RAD, and a
  future calibrated CAD-like model.
- Keep the paper-supported equations in shared math code and expose them in the
  UI/readouts.
- Keep inverse-design results labeled as proposals until a physical validation
  pass compares the command set against spring-hinge or measured response.
- Treat locks, backlash dead zones, actuator commands, and clearance as
  composable operators with measurable locality, reachability, and
  noncommutativity.
- Use the current interaction hotspot and degree maps to prioritize which local
  operator pairs should be studied next with physical pair-cell experiments.
