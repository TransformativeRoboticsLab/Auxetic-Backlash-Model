# Contact And Motion Limits

## Why Limits Are Needed

The hard pin constraints can mathematically close for angles where bodies would
pass through each other. Those angles must be rejected. The actuator range should
only include angles where both of these are true:

```text
all pin residuals <= tolerance
max inter-cell penetration <= tolerance
```

## Primitive Collision Model

Each cross body is approximated by 2D primitives with vertical `Z` bands:

- horizontal arm capsule
- vertical arm capsule
- center hub disk
- east pad disk
- north pad disk
- west pad disk
- south pad disk

The contact guard compares only bodies from different cells. Bodies inside the
same cell are stacked and share a center axis, so their projected overlap is not
treated as invalid.

## Vertical Bands

```text
lower cross: z in [-4.0, 0.0]
upper cross: z in [0.0, 4.0]
```

Two primitives can collide only if their vertical bands overlap.

## Clearance Calculation

For every inter-cell primitive pair:

```text
clearance = distance_between_primitives - radius_1 - radius_2
```

If `clearance < 0`, penetration exists.

The current valid threshold is:

```text
BODY_EPSILON_MM = 0.0001
```

## Motion Range

The browser primitive samples the requested uniform drive range:

```text
full drive scan = [-80 deg, 80 deg]
sample step = 1 deg
```

An angle is allowed only if the no-backlash closure and contact guard both pass.
For the current default three-cell row, the visible allowed interval is:

```text
26 deg to 64 deg
```

The slider is clamped to the allowed interval. Animation also uses that allowed
interval instead of the full mathematical range.

## Future Improvements

The primitive contact guard is useful for preventing obvious body
interpenetration, but it is not exact physics. Future versions should replace or
augment it with:

- segmented CAD body meshes
- actual pin and hole surfaces
- rigid-body contact in MuJoCo, Isaac Sim, Gazebo, or another solver
- friction and restitution parameters
- validation against retroreflector data

