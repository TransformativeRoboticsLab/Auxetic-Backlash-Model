# Two Cell Mechanics Data Contract

This file defines what the simulator needs from the Fusion two-cell mechanism.
It is intentionally more specific than a normal CAD export because the goal is
to derive function, not just render geometry.

## Required Files

Place these in `exports/`:

- `two_cell_design_graph.json`: semantic CAD graph.
- `two_cell.step`: separated geometry for the full two-cell assembly.
- `two_cell_pose_samples.csv` or `two_cell_pose_samples.json`: sampled poses
  across actuation.
- Optional: separated `.step`, `.stl`, or `.obj` files for each moving body.

## Required Metadata

The export should include:

- `units`: length and angle units.
- `source_file`: native Fusion file name.
- `export_timestamp`: when the export was produced.
- `coordinate_frame`: definition of X, Y, Z, origin, and cell ordering.
- `baseline_state`: the actuation value considered ground or neutral.
- `actuation_range`: minimum, maximum, and step size.

## Required CAD Graph Nodes

Each node needs a persistent name or path:

- root assembly;
- cell occurrences;
- moving bodies;
- fixed bodies;
- pins;
- holes;
- plates or links;
- joints;
- joint origins;
- user parameters;
- model parameters if accessible;
- tracked marker points.

## Required CAD Graph Edges

At minimum:

- `contains`: assembly or component contains a body, joint, sketch, or feature;
- `instance_of`: occurrence is an instance of a component;
- `connects`: joint connects two occurrences or bodies;
- `drives`: parameter drives a joint, sketch dimension, or feature;
- `tracks`: named marker belongs to a body or occurrence;
- `contact_candidate`: pin/hole or link/link pair that may contact later.

## Required Pose Sample Columns

For each sampled actuation value:

- sample id;
- actuation value;
- left cell alpha;
- right cell alpha;
- left cell theta;
- right cell theta;
- left cell center x, y, z;
- right cell center x, y, z;
- pitch between cell centers;
- every tracked marker x, y, z;
- every moving body transform as a 4x4 matrix or translation/quaternion pair;
- interference/contact flag if Fusion can compute it.

## First Acceptance Criteria

The no-backlash two-cell model is ready to replace the proxy when:

- held-out CAD pose error is reported for every tracked point;
- center pitch is reproduced without arbitrary contraction constants;
- the pair is symmetric under symmetric actuation;
- the simulator can show exactly two connected cells and animate through the CAD
  sampled range;
- the old normalized cell view remains available as an abstraction mode;
- all limitations are documented before calling the model physically accurate.
