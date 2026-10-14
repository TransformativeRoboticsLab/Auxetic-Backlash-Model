# Auxetic Lattice Backlash Model
This repository contains a set of programs to model the backlash in auxetic lattice structures.

The backlash in an auxetic structure can be configured or programmed to design a shape-changing structure.

## Open the RAD Simulator

Start with the [V2 opening guide](Ahyan%20Virtual%20Simulation/V2/OPEN_SIMULATOR.md).
It covers downloading the repository, opening the simulator in Chrome without
installing development tools, basic controls, and troubleshooting.

The current simulator is on the
[`codex/web-simulator-upgrade` branch](https://github.com/TransformativeRoboticsLab/Auxetic-Backlash-Model/tree/codex/web-simulator-upgrade),
under `Ahyan Virtual Simulation/V2/web/two-cell-attachment/index.html`.
Downloading just the HTML file is not enough; keep its supporting files together.

Task List
https://docs.google.com/spreadsheets/d/1UB8SRYQrKouYJ8RqrWCk4LBVXZi8fhDDzx1xIQ8QW4M/edit?usp=sharing

## Main Questions
- How can we build valuable shapes using the backlash present in series and parallel networks of auxetic linkages?
- What metrics in the structure lead to changes in a configuration space?
- How do we model backlash

## Table of Contents
The programs and code for this work is split up into multiple files based on specific tasks.

- main_1D_dofs.py
  - functions to determine the DoF of 1D linkages
- main_in_2D.py
  - functions to explore 2D area differences in auxetic, high-backlash structures
- main_in_3D.py
  - functions to explore 3D deformation of auxetic, high-backlash structures
- main_airfoil.py
  - functions to examine 3D airfoil shape conformal placement
- main_configuration_space.py
  - functions to understand the configur


#### ReLU Function for Stiffness of Cells
To model the rotational stiffness of a rotating, auxetic cell, we use a realized linear unit equation.
![alt text](./figures/ReLU_demo_figure.png "ReLU model")


