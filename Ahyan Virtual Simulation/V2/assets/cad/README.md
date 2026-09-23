# RAD CAD Assets

`RADs_unit_cell.f3d` was downloaded from the public Autodesk A360 share:
https://a360.co/4bMlzip

The Fusion archive manifest identifies the asset as `RADs unit cell` and
references a linked cell labeled `RADs free cell 4mm tall 3.4mm hole`.
`RADs_unit_cell_preview.png` is the preview image extracted from the archive.
`outputs/two_cell_bench_packet/cad_rad_cell_layout.json` is a generated
dimension-consistent layout that uses the CAD bounding box, 3.4 mm nominal hole
label, eight radial pad sites, and three two-cell connector pairs. The generated
two-cell connector-contact reports use those pairs to track lateral slip,
vertical slip, clearance excess, and contact penalty while hole radius and lock
mode vary.

The `.f3d` file contains Fusion design streams and Autodesk ShapeManager BREP
blobs. The generated layout is not exact BREP extraction. For external physics
engines, export STEP/STL/OBJ from Fusion with separated moving bodies, pin axes,
hole surfaces, and material assignments.

For a quick one-cell mesh audit, export OBJ or STL from the A360/Fusion viewer
and run:

```powershell
python -m rad_sim.cad_mesh_intake assets/cad/RADs_unit_cell.obj --out outputs/cad_mesh_audit.json --pin-radius-mm 1.36 --hole-radius-mm 1.70 --boss-radius-mm 3.10 --backlash-mm 0.34
```

The output is a `rad-sim.cad-mesh-audit.v1` wrapper containing a loadable
`rad-sim.hardware-profile.v1`. Browser `Load Profile` accepts either the nested
hardware profile or the full audit JSON.

Use `assets/cad/segmented/` as the intake folder for those exports. The expected
filenames are listed in `assets/cad/segmented/README.md`; the generated
`outputs/two_cell_bench_packet/two_cell_segmented_cad_readiness.json` report
detects which assets are present.
