# RAD Virtual Simulation V2

This folder is a source-only V2 snapshot of the Reconfigurable Auxetic Devices Digital Workbench. It keeps the current simulation framework together for GitHub review without generated temp files, bulky build artifacts, or runtime output folders.

Key entry points:
- Browser simulator: `web/index.html`
- Two-cell attachment simulator: `web/two-cell-attachment/index.html`
- Python package: `rad_sim/`
- Constraint solver report: `python -m rad_sim.constraint_kinematic_reports`
- Tests: `python -m unittest tests.test_rad_sim` and `python -m unittest tests.test_web_static`

The Python constraint-kinematic solver is the quantitative reference path. Browser constraint mode is a fast visualization/diagnostic mirror.
