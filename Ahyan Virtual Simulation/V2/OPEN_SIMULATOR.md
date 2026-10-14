# Open the Reconfigurable Auxetic Devices Digital Workbench

This guide opens the current V2 pin-hole lattice simulator on your own computer.
The simplest route needs Chrome and the downloaded repository, not Python,
Node.js, Fusion, or an API key. The browser must support WebGL.

## Download and open without installing tools

1. Open the [simulator branch on GitHub](https://github.com/TransformativeRoboticsLab/Auxetic-Backlash-Model/tree/codex/web-simulator-upgrade).
   Confirm the branch selector says `codex/web-simulator-upgrade`, not `main`.
   Sign in with an account that has repository access if GitHub requests it.
2. Choose **Code > Download ZIP**. Extract the entire ZIP to a normal folder;
   do not open the application from inside the ZIP preview.
3. Inside the extracted repository, open these folders in order:
   `Ahyan Virtual Simulation` > `V2` > `web` > `two-cell-attachment`.
4. Open `index.html` with Google Chrome. On Windows, use right-click >
   **Open with > Google Chrome**, or drag the file into a Chrome window.
5. You should see **Pin-Hole RAD Attachment Lattice**, a 3D lattice, view buttons,
   and actuation controls. The browser address should begin with `file:///`.

Use this specific `index.html`, not the older `V2/web/index.html`. Viewing an HTML
file on GitHub shows its source; it does not run the simulator. Once downloaded,
the simulator's bundled scripts and Three.js library work locally without a CDN.

Keep the complete download together. In particular, the page needs:

```text
Ahyan Virtual Simulation/V2/web/
  vendor/three.min.js
  two-cell-attachment/
    index.html
    styles.css
    two_cell_attachment.js
    backlash_state.js
    pin_geometry.js
```

## First experiment

1. Set **Rows** to `1`, **Columns** to `10`, then click **Apply lattice size**.
2. Keep **Uniform no-backlash baseline** selected and move the rotational slider.
   This mode sets the pin-to-hole radius ratio to one automatically.
3. Switch to **Single-cell contact propagation**, adjust the pin radius ratio
   below one to introduce clearance, click a cell, and move the slider.
4. Select another cell to continue actuating from the existing configuration.
   Selection alone should not change the lattice configuration; the slider
   adopts that cell's current angle.
5. To anchor a cell, select it, choose **XY fixed, Z free** or **XYZ fixed**, and
   press **Lock selected position**. One positional anchor keeps the selected
   coordinates fixed without reducing the free actuation range. It does not
   lock the cell's internal rotation.
6. Use **Unlock selected position** or **Clear all locks** to remove constraints.
   Multiple anchors can restrict travel or jam the mechanism. Start with one
   anchor before trying more complicated boundary conditions.

Left-drag rotates the view. Right-drag or Shift-drag pans. The mouse wheel zooms.
Use **Top**, **Front**, **Side**, or **Iso** for standard views, and **Fit** to
frame the lattice. **Reset** resets the camera view, not the mechanical state.
**Animate drive** runs the selected actuation mode automatically.

Start with a small lattice. Dimensions up to 30 by 30 are supported, but large
single-cell solves can take substantially longer. This is an uncalibrated
kinematic research model, not validated full contact dynamics; see the
[model notes and limitations](README_V2.md).

## Optional local server

If your browser or organization blocks local HTML files, use Python 3 to serve
the files. From a terminal in the extracted or cloned repository root:

```sh
cd "Ahyan Virtual Simulation/V2/web"
python -m http.server 8765 --bind 127.0.0.1
```

On Windows, `py -3` can replace `python`; on macOS/Linux, use `python3` if needed.
Open [the local simulator](http://127.0.0.1:8765/two-cell-attachment/index.html).
Keep the terminal open while using it and press Ctrl+C to stop the server.
If port 8765 is occupied, use 8766 in both the command and browser address.
This localhost address runs on your computer; it is not a link others can use
to access your session. No server is needed for the normal file-opening method.

## Clone or update with Git

For a first checkout:

```sh
git clone --branch codex/web-simulator-upgrade --single-branch https://github.com/TransformativeRoboticsLab/Auxetic-Backlash-Model.git
```

Then follow the opening instructions above inside `Auxetic-Backlash-Model`.
For an existing checkout, preserve any local work before switching branches or
updating. From that repository root:

```sh
git fetch origin
git switch codex/web-simulator-upgrade
git pull --ff-only
```

If Git reports conflicting local changes or diverged history, stop and resolve
them with a collaborator; do not discard your work or force an update. If you
downloaded a ZIP instead, download a fresh ZIP into a separate folder for updates.

## Troubleshooting

- **Old interface or missing controls:** check the branch and exact V2 path above.
  Reload the browser after updating the files. This version is not yet merged
  into `main`.
- **Blank page or missing scripts:** extract the full ZIP and preserve the folder
  layout, especially `web/vendor`. Do not download only `index.html`.
- **Controls appear but no 3D lattice:** press **Fit**, try a small lattice, and
  check whether Chrome reports a WebGL/graphics error. Managed computers may
  require IT assistance to enable browser graphics. Do not disable browser security.
- **Disabled actuation:** clear all position locks to distinguish a constrained
  configuration from a loading problem. Two fixed ends with zero clearance can jam.
- **Still failing:** share your operating system, Chrome version, exact page path,
  lattice dimensions, actuation mode, pin ratio, lock locations/axes, and a screenshot.
  Include any red errors from Chrome's Developer Tools Console.
