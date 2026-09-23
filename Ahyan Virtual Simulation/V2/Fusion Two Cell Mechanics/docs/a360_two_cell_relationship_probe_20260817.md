# A360 Two Cell Relationship Probe

Date: 2026-08-17

Source share: https://a360.co/3U08kVf

Local probe folder:

`C:\Users\ahyan\OneDrive\Documents\TRL\Ahyan Virtual Simulation\cad_probe\a360_two_cell_relationship_20260817`

## What Loaded

The A360 share loaded successfully in Chrome headless with no page errors.

Viewer title:

`TRL RADs two cells - AUTODESK FUSION`

Share metadata:

- `shareTitle`: `TRL RADs two cells`
- `fileName`: `TRL RADs two cells`
- `fileType`: `f3d`
- `mimeType`: `application/vnd.autodesk.fusion360`
- `downloadEnabled`: `true`
- `conversionStatus`: `FINISHED`
- viewer status: `ready`, `100%`

Downloaded viewer derivatives:

- `chrome_metadata.json`
- `chrome_manifest.json`
- `03_Autodesk.CloudPlatform.DesignDescription.Registered_designdescription.json`
- `02_Autodesk.CloudPlatform.PropertyDatabase_properties.db`
- `01_graphics_Design.svf`
- `properties_summary.json`
- extracted recording frames and contact sheet

## Relationship Evidence Found

The design-description derivative shows that the two-cell file is an XRef
assembly:

- root design: `TRL RADs two cells`
- referenced design: `RADs unit cell`
- referenced subdesigns inside the unit cell:
  - `RADs free cell 4mm tall 3.4mm hole`
  - `92125A136_18-8 Stainless Steel Hex Drive Flat Head Screw`

The property database exposes the object tree:

- `TRL RADs two cells v2`
  - `RADs unit cell v2:1`
    - `RADs free cell 4mm tall 3.4mm hole v3:1`
    - screw component
    - `Body1`
  - `RADs unit cell v2:2`
    - `RADs free cell 4mm tall 3.4mm hole v3:1`
    - `Body1`

This means the file-level Fusion relationship/XRef structure is present: the
two-cell model really is composed from unit-cell references rather than being
only an unrelated screenshot.

## What Is Not Proven

The public A360 viewer derivatives do not expose enough information to prove the
mechanical actuation relationship.

Missing from the public derivative evidence:

- explicit Fusion joints;
- joint axes;
- joint limits;
- motion links or motion studies;
- named actuation parameter;
- parametric dependency graph;
- sampled transforms over an actuation sweep;
- closure-constraint residuals.

So the answer is:

- File-level relationship/XRef: works.
- Static two-cell assembly visibility: works.
- Exact actuation relationship between the two cells: not proven from A360
  viewer data alone.

## Visual Read

The supplied screen recording and viewer screenshot show a plausible two-cell RAD
assembly with the two radial unit cells placed together and central connector
geometry present. The clip appears to be camera inspection of a static model, not
a measured actuation sweep. It therefore helps us understand layout, but not the
mathematical movement law.

## Native Download Attempt

The share reports `downloadEnabled: true`, but automated download from the
headless viewer did not emit a downloadable file. Direct Autodesk object-storage
download of the root `.f3d` returned HTTP 403 with the public viewer token.

Practical implication: download the `.f3d` manually from the A360 page or open it
directly in Fusion, then run the Fusion exporter script.

## Next Required Export

To verify and implement the true two-cell no-backlash mechanism, export one of
these:

1. Best: Fusion MCP/API output from the live file.
2. Good: the native `.f3d` placed in `Fusion Two Cell Mechanics\exports`.
3. Sufficient for first calibration: JSON graph plus STEP plus pose samples.

The needed pose sweep is:

- actuation parameter value;
- left and right cell centers;
- left and right alpha/theta;
- transforms of each moving body;
- tracked pin/hole or marker coordinates;
- closure error for connector joints.

Until those exist, the simulator should not claim exact two-cell mechanics.
