# CAD Probe Notes

The V2 GitHub snapshot intentionally excludes the raw `cad_probe` cache because
the Autodesk derivative downloads contain very long generated filenames that fail
normal Git indexing on Windows.

For future CAD/Fusion reconstruction work, regenerate or store the probe
artifacts outside the source tree, then extract short-name summaries into this
folder. The simulator source, two-cell attachment web model, and Fusion two-cell
STEP export are still included in V2.
