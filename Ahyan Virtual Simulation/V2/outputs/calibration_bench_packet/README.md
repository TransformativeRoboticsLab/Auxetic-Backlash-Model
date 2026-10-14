# RAD Calibration Bench Packet

This folder contains simulator-derived planning artifacts for RAD
calibration experiments. Treat it as a bench handoff, not as
calibrated hardware validation.

## Files

- `calibration_bench_packet.json`: complete JSON packet with notebook, protocol, templates, and validation guidance.
- `calibration_bench_notebook.json`: bench notebook with instruments, phases, dataset plan, and scenario rows.
- `calibration_bench_notebook.csv`: one-row-per-scenario CSV table for a lab notebook or spreadsheet.
- `calibration_experiment_protocol.json`: machine-readable calibration protocol used by simulator comparison code.
- `calibration_fit_results_template.json`: blank fit measurement rows. Fill this before fitting a profile.
- `calibration_holdout_results_template.json`: blank holdout rows tied to the frozen profile ID.
- `hardware_profile.json`: hardware dimensions used to derive the effective normalized model.

## Bench Notes

- Keep the fit and holdout raw files separate.
- Do not inspect or tune against the holdout rows before freezing the profile.
- Set `profileFrozenAt` after the fit profile is generated and before holdout collection.
- Replace null measured fields with bench data; do not copy simulator predictions into measured rows.
- Passing residual checks supports calibration bookkeeping only, not contact, friction, stiffness, or material-law truth.
