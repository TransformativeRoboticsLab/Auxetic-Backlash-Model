# Lock Dataset Notes for Final Presentation

Status: working presentation/research notes  
Dataset path: `C:\Users\ahyan\OneDrive\Desktop\logs`  
Observed files: 37 configuration CSV files plus `marker_xyz_means_all_37_files.csv`

Implemented support:

- `rad_sim.lock_dataset.load_lock_dataset(...)` loads the cleaned marker
  summary and ignores non-configuration CSV rows.
- `predict_lock_strand_shape(...)` gives a conservative nearest/RBF empirical
  surrogate from nominal locks, crown angles, realized locks, and stable-state
  index to measured marker shape.
- Known consecutive-lock realization failures from the notes are encoded,
  including the edge case where `10_11_12` realizes as only lock `12`.
- `calibrate_single_strand_from_lock_dataset(...)` compares a simulated single
  strand against measured marker geometry after remapping lab axes
  (`X/Z` ground plane, `Y` vertical) into simulator axes.
- `calibrate_repeated_row_locks_from_lock_dataset(...)` builds a repeated-row
  sheet overlay when every row has the same nominal lock mask.

Current limitation: these are calibration overlays and diagnostics. They make
the measured data available to the simulator without claiming that the current
kinematic lock solver is physically complete.

## Why This Data Matters

The current simulator lock model should not be treated as physically accurate.
The cleaned 12-cell lock dataset gives direct evidence for how discrete locks
actually reshape the lattice. This should become the calibration source for the
lock operator rather than assuming that a locked cell simply freezes an ideal
cell angle.

Presentation framing:

```text
ideal lock model -> measured lock response -> calibrated lock operator
```

The neural-network approximation can be presented as a surrogate model trained
on measured configurations, while the mathematical framework should treat the
underlying lock as an empirical discontinuity operator.

## Important Physical Nuance

Although the crowns may be labeled as 40, 30, or 20 degree locks, the effect on
the actual lattice angle theta is smaller. The nonzero width of the cell arm
means the crown angle does not transfer one-to-one into the effective rotating
lattice angle.

Modeling implication:

```text
theta_effective != crown_angle
theta_effective = f(crown_angle, arm_width, contact_location, backlash, compliance)
```

This should be shown as a correction to the naive geometric lock assumption.

## Consecutive Lock Failure Mode

For configurations where locks are placed consecutively, the middle lock often
causes no measurable difference when placed or omitted. Because crowns need to
be placed underneath the cells, the middle crown in a consecutive series
repeatedly falls out.

Presentation phrasing:

```text
Consecutive lock patterns contain a hardware-realization failure mode:
the middle lock may be nominally present in the configuration file but
physically absent during the realized stable state.
```

Modeling implication:

```text
nominal_lock_mask != realized_lock_mask
```

The simulator should eventually distinguish commanded/nominal locks from
physically realized locks.

## File Naming Convention

Ground-state examples:

```text
AOM_L0_#
```

Locked configuration examples:

```text
AOM_L3_#*#*#_303030_1
```

Interpretation:

- `AOM`: auxetic/origami/metamaterial experiment prefix.
- `L0`, `L3`, etc.: number of locks in the configuration.
- Numbers after `L#`: cell locations where locks were placed.
- `303030`: lock types/angles, here three 30 degree locks.
- Final number: stable state index.

State convention:

- For bistable configurations, state `1` begins in the lowest potential state
  and state `2` is the higher potential state.
- The same ordering convention extends to multistable configurations.

## Ground States

There are three ground states:

```text
AOM_L0_# - 3 stable states
```

These should be used as baseline shapes before comparing lock-induced
configuration changes.

## Three-Lock 30 Degree Dataset

Completed or listed configurations:

| Configuration | Stable states | Notes |
| --- | ---: | --- |
| `AOM_L0_#` | 3 | Ground states |
| `AOM_L3_1_2_3_303030_#` | 2 | Lock 2 fell |
| `AOM_L3_4_5_6_303030_#` | 2 | Lock 5 fell |
| `AOM_L3_7_8_9_303030_#` | 2 | Lock 8 fell |
| `AOM_L3_10_11_12_303030_#` | 2 | Locks 10 and 11 fell |
| `AOM_L3_1_4_7_303030_#` | 2 | Construction vibration affected measurement |
| `AOM_L3_1_5_8_303030_#` | 2 |  |
| `AOM_L3_1_6_9_303030_#` | 2 |  |
| `AOM_L3_2_4_8_303030_#` | 2 |  |
| `AOM_L3_2_5_7_303030_#` | 2 |  |
| `AOM_L3_2_6_10_303030_#` | 2 |  |
| `AOM_L3_2_9_11_303030_#` | 1 |  |
| `AOM_L3_3_4_9_303030_#` | 2 |  |
| `AOM_L3_3_5_10_303030_#` | 2 |  |
| `AOM_L3_3_6_7_303030_#` | 2 |  |
| `AOM_L3_3_8_11_303030_#` | 1 |  |
| `AOM_L3_5_9_12_303030_#` | 2 |  |
| `AOM_L3_6_8_12_303030_#` | 2 |  |
| `AOM_L3_1_2_4_303030_#` | 2 |  |

Checkpoint recommendation:

```text
Break and test whether the model works with only 30 degree locks before
continuing to 4-lock configurations.
```

## Measurement Caveats

- A first ground-plane estimate was badly tilted and about an hour of data had
  to be scrapped and restarted.
- Construction across the street introduced vibration and made measurements
  harder for at least part of the dataset.
- These caveats should be presented as experimental realities, not hidden.

## How This Should Enter the Framework

The lock should be modeled as a physical realization map:

```text
rho(lock_operator, hardware_parameters, placement_conditions)
  -> realized_lock_operator
```

This distinguishes:

- nominal lock placement,
- realized lock engagement,
- failed or fallen locks,
- effective theta constraint,
- neighbor deformation response,
- stable-state branch.

Recommended operator model:

```text
lock_data_record =
  (lock_mask,
   crown_angles,
   realized_lock_mask,
   state_index,
   marker_positions,
   effective_shape,
   notes)
```

The simulator should eventually learn or fit:

```text
F(lock_mask, crown_angles, boundary_conditions, state_index)
  -> measured marker/node shape
```

The neural network can be used as an empirical surrogate for `F`, while a more
interpretable model should be fit in parallel:

- effective theta reduction per crown angle;
- partial lock compliance;
- failed-lock correction for consecutive locks;
- local and nonlocal shape influence;
- multistable branch selection.

## Final Presentation Message

The most important point is that the experimental lock data shows why a naive
lock model fails. Locks are not ideal constraints. They are programmable
discontinuities with hardware realization errors, reduced effective theta, and
state-dependent outcomes. This supports the broader thesis that programmable
mechanics needs both an operator framework and calibration data.
