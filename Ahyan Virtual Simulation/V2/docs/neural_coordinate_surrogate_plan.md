# Neural Coordinate Surrogate Plan

Goal: learn a map from programmable lattice state to every cell's 3D
coordinate:

```text
F(lock mask, lock angles, realized locks, actuators, boundary conditions,
  removed cells, stable-state branch, hardware parameters)
    -> R^(rows x cols x 3)
```

## What Exists Now

- The cleaned 12-cell lock logs are converted into per-cell 3D coordinate
  targets.
- The browser loads `web/data/lock_dataset.js` and can apply the measured
  12-cell strand shape to a single strand or repeated identical row locks.
- `train_lock_coordinate_mlp(...)` trains a small NumPy MLP from lock features
  to 12 per-cell 3D coordinates.
- The neural model is a scaffold, not a validated mechanics model. The current
  dataset is too small for strong generalization.

## What To Collect Next

1. Keep the 12-cell single-strand setup fixed and collect more lock patterns.
2. Include 0, 1, 2, 3, and 4 lock configurations, not only 3-lock cases.
3. Repeat each configuration at least 3 times to estimate noise.
4. Record all stable branches with the same naming convention.
5. Save nominal locks and realized locks separately.
6. Measure actuator commands, if any, as continuous values.
7. Record boundary conditions: fixed cells, clamped ends, free ends, gravity
   orientation, and ground-plane reference.
8. Add sheet experiments where identical locks are repeated across rows.
9. Add non-identical sheet experiments where rows have different lock masks.
10. Keep raw marker IDs stable so the model learns consistent node/cell order.

## Model Progression

1. Nearest/RBF measured-shape surrogate for small data.
2. Small MLP for lock-only single-strand coordinate prediction.
3. Graph neural network for arbitrary row/column sheets.
4. Recurrent/sequence model only if dynamic actuation trajectories are recorded.
5. Physics-informed loss that penalizes impossible bar stretch, hinge motion, and
   contact violations.

## Acceptance Target

The model becomes useful when it predicts held-out marker/cell coordinates with
errors smaller than experimental repeatability. Until then, use it as an
empirical guide and keep the mechanics/Lean framework separate from learned
claims.
