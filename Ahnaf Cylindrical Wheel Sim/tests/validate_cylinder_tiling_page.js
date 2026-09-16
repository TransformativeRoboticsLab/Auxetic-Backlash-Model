const assert = require("assert");
const fs = require("fs");
const path = require("path");
const { pathToFileURL } = require("url");
const { PNG } = require("pngjs");
const { chromium } = require("playwright");

const root = path.resolve(__dirname, "..");
const pageUrl = pathToFileURL(path.join(root, "index.html")).href;

function colorCount(buffer, predicate) {
  const png = PNG.sync.read(buffer);
  let count = 0;
  for (let index = 0; index < png.data.length; index += 4) {
    const r = png.data[index];
    const g = png.data[index + 1];
    const b = png.data[index + 2];
    const a = png.data[index + 3];
    if (a > 20 && predicate(r, g, b)) count += 1;
  }
  return count;
}

async function launchBrowser() {
  try {
    return await chromium.launch({ channel: "chrome" });
  } catch (_error) {
    return chromium.launch();
  }
}

async function readMetrics(page) {
  return page.evaluate(() => ({
    alpha: document.getElementById("alpha").value,
    backlash: document.getElementById("backlash").value,
    ringCount: document.getElementById("ringCount").value,
    rowCount: document.getElementById("rowCount").value,
    alphaCommand: document.getElementById("alphaCommandMetric").textContent,
    alphaEffective: document.getElementById("alphaEffectiveMetric").textContent,
    theta: document.getElementById("thetaMetric").textContent,
    deadzone: document.getElementById("deadzoneState").textContent,
    pitch: document.getElementById("pitchMetric").textContent,
    turn: document.getElementById("turnMetric").textContent,
    diameter: document.getElementById("diameterMetric").textContent,
    radius: document.getElementById("radiusMetric").textContent,
    closure: document.getElementById("closureMetric").textContent,
    totalCells: document.getElementById("totalCellsMetric").textContent,
    height: document.getElementById("heightMetric").textContent,
    collisionState: document.getElementById("collisionState").textContent,
    penetration: document.getElementById("penetrationMetric").textContent,
    clearance: document.getElementById("clearanceMetric").textContent,
    pairsChecked: document.getElementById("pairsCheckedMetric").textContent,
    selectedStatus: document.getElementById("selectedCellStatus").textContent,
    actuatorCount: document.getElementById("actuatorCountMetric").textContent,
    lockedCount: document.getElementById("lockedCountMetric").textContent,
    selectedCellAlpha: document.getElementById("selectedCellAlphaMetric").textContent,
    minDiameter: document.getElementById("minDiameterMetric").textContent,
    maxDiameter: document.getElementById("maxDiameterMetric").textContent,
  }));
}

function mmValue(text) {
  const match = text.match(/^(-?\d+(?:\.\d+)?) ?mm$/);
  assert.ok(match, `could not parse mm metric: ${text}`);
  return Number(match[1]);
}

async function checkViewport(browser, name, viewport) {
  const page = await browser.newPage({ viewport });
  const errors = [];
  page.on("console", (message) => {
    if (message.type() === "error") errors.push(message.text());
  });
  page.on("pageerror", (error) => errors.push(error.message));
  await page.goto(pageUrl);
  await page.waitForSelector("#threeMount canvas", { state: "visible" });
  await page.waitForTimeout(500);

  // Expand every collapsible panel section so its controls/readouts are
  // reachable regardless of which ones default to collapsed.
  await page.evaluate(() => {
    document.querySelectorAll("details.panel-section").forEach((details) => {
      details.open = true;
    });
  });

  const canvasSize = await page.evaluate(() => {
    const canvas = document.querySelector("#threeMount canvas");
    return { width: canvas.width, height: canvas.height };
  });
  assert.ok(canvasSize.width > 300, `${name} canvas should have real width`);
  assert.ok(canvasSize.height > 300, `${name} canvas should have real height`);

  const overflow = await page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth);
  assert.ok(overflow <= 2, `${name} should not create horizontal overflow`);

  const defaults = await readMetrics(page);
  assert.strictEqual(defaults.ringCount, "10", `${name} default ring count should be 10`);
  assert.strictEqual(defaults.rowCount, "3", `${name} default row count should be 3`);
  assert.strictEqual(defaults.totalCells, "30", `${name} default total cells should be n*m`);
  assert.strictEqual(defaults.alpha, "1.3", `${name} default alpha should be 1.3`);
  assert.strictEqual(defaults.alphaCommand, "1.30", `${name} commanded alpha readout should mirror the slider`);
  // Dead-zone shifts the engaged effective alpha down by the gap width:
  // effective = 1.0 + reluDeadzone(alpha-1.0, backlash) = 1.0 + (0.3-0.1) = 1.2
  assert.strictEqual(defaults.alphaEffective, "1.20", `${name} at backlash=0.1 and |alpha-1|=0.3>0.1, effective alpha should be shifted by the gap width to 1.20`);
  assert.strictEqual(defaults.deadzone, "engaged", `${name} default alpha=1.3 should be outside the b=0.1 dead zone around 1.0`);
  const theta = Number(defaults.theta.replace(/ ?deg$/, ""));
  assert.ok(Math.abs(theta - (70 * 1.2 - 60)) <= 0.05, `${name} theta should follow theta = 70*effectiveAlpha - 60`);
  assert.ok(Math.abs(mmValue(defaults.closure)) <= 0.01, `${name} default ring closure residual should be ~0`);
  const defaultDiameter = mmValue(defaults.diameter);
  const defaultRadius = mmValue(defaults.radius);
  assert.ok(Math.abs(defaultDiameter - 2 * defaultRadius) <= 0.2, `${name} diameter should equal 2x radius`);
  assert.ok(defaultDiameter > 0, `${name} default diameter should be positive`);
  const turnDeg = Number(defaults.turn.replace(/ ?deg$/, ""));
  assert.ok(Math.abs(turnDeg - 36) <= 0.1, `${name} turn angle for n=10 should be 36 deg (360/10)`);
  const defaultHeight = mmValue(defaults.height);
  assert.ok(Math.abs(defaultHeight - 2 * 50) <= 0.5, `${name} default assembly height should be (m-1)*axialPitch`);

  // Backlash dead-zone: alpha exactly at the reference (1.0) with a nonzero
  // gap should read as "free (dead zone)" and pin theta at 70*1-60=10deg.
  const deadzoned = await page.evaluate(() => {
    document.getElementById("alpha").value = "1.0";
    document.getElementById("alpha").dispatchEvent(new Event("input", { bubbles: true }));
    return {
      deadzone: document.getElementById("deadzoneState").textContent,
      theta: document.getElementById("thetaMetric").textContent,
      alphaEffective: document.getElementById("alphaEffectiveMetric").textContent,
    };
  });
  assert.strictEqual(deadzoned.deadzone, "free (dead zone)", `${name} alpha=1.0 (the reference) should sit inside the backlash dead zone`);
  assert.strictEqual(deadzoned.alphaEffective, "1.00", `${name} effective alpha inside the dead zone should stay at the reference`);
  const deadzonedTheta = Number(deadzoned.theta.replace(/ ?deg$/, ""));
  assert.ok(Math.abs(deadzonedTheta - 10) <= 0.05, `${name} theta at the dead-zone reference should be 70*1-60=10deg`);

  // The dead-zone ruler under the Drive slider should visually agree with
  // the numeric state: at alpha=1.0 (still set from above) the marker
  // should fall inside the shaded band; the band's own bounds should match
  // the "alpha dead zone" readout (0.90-1.10 at the default backlash=0.1).
  const ruler = await page.evaluate(() => {
    const band = document.getElementById("deadzoneBand");
    const marker = document.getElementById("deadzoneMarker");
    return {
      bandLeft: parseFloat(band.style.left),
      bandWidth: parseFloat(band.style.width),
      markerLeft: parseFloat(marker.style.left),
      rangeText: document.getElementById("deadzoneRangeMetric").textContent,
    };
  });
  assert.strictEqual(ruler.rangeText, "0.90 - 1.10", `${name} the dead zone readout should match alpha=1 +/- backlash=0.1`);
  assert.ok(ruler.bandWidth > 0, `${name} the dead-zone band should have nonzero width`);
  assert.ok(
    ruler.markerLeft >= ruler.bandLeft && ruler.markerLeft <= ruler.bandLeft + ruler.bandWidth,
    `${name} at alpha=1.0 the marker (${ruler.markerLeft}%) should fall inside the band [${ruler.bandLeft}%, ${ruler.bandLeft + ruler.bandWidth}%]`
  );

  await page.evaluate(() => {
    document.getElementById("alpha").value = "1.3";
    document.getElementById("alpha").dispatchEvent(new Event("input", { bubbles: true }));
  });

  // Sweep alpha across its full range; the ring should stay closed at every setting.
  for (const which of ["min", "max"]) {
    const swept = await page.evaluate((minOrMax) => {
      const input = document.getElementById("alpha");
      input.value = input[minOrMax];
      input.dispatchEvent(new Event("input", { bubbles: true }));
      return document.getElementById("closureMetric").textContent;
    }, which);
    assert.ok(Math.abs(mmValue(swept)) <= 0.01, `${name} ring closure residual should stay ~0 at alpha ${which}`);
  }
  await page.evaluate((value) => {
    const input = document.getElementById("alpha");
    input.value = value;
    input.dispatchEvent(new Event("input", { bubbles: true }));
  }, defaults.alpha);

  // Material clearance guard should be reporting something sane at defaults.
  assert.notStrictEqual(defaults.collisionState, "not checked", `${name} clearance checking should be on by default`);
  assert.ok(defaults.clearance.endsWith("mm"), `${name} should report minimum clearance`);
  assert.ok(Number(defaults.pairsChecked) > 0, `${name} should report a positive number of clearance pairs checked`);

  // Resize the ring (n) and stack (m); pool should rebuild and stay closed.
  const resized = await page.evaluate(() => {
    const ringInput = document.getElementById("ringCount");
    const rowInput = document.getElementById("rowCount");
    ringInput.value = ringInput.max;
    ringInput.dispatchEvent(new Event("input", { bubbles: true }));
    rowInput.value = rowInput.max;
    rowInput.dispatchEvent(new Event("input", { bubbles: true }));
    return {
      totalCells: document.getElementById("totalCellsMetric").textContent,
      closure: document.getElementById("closureMetric").textContent,
      ringCount: ringInput.value,
      rowCount: rowInput.value,
    };
  });
  assert.strictEqual(
    resized.totalCells,
    String(Number(resized.ringCount) * Number(resized.rowCount)),
    `${name} total cells should track n*m after resizing`
  );
  assert.ok(Math.abs(mmValue(resized.closure)) <= 0.01, `${name} ring closure residual should stay ~0 after resizing to n=${resized.ringCount}`);

  await page.evaluate((defaultsValue) => {
    document.getElementById("ringCount").value = defaultsValue.ringCount;
    document.getElementById("ringCount").dispatchEvent(new Event("input", { bubbles: true }));
    document.getElementById("rowCount").value = defaultsValue.rowCount;
    document.getElementById("rowCount").dispatchEvent(new Event("input", { bubbles: true }));
  }, defaults);
  await page.waitForTimeout(50);

  // Cell selection: click near the center of the viewport (where a cell should be under the default iso view).
  const selection = await page.evaluate(() => {
    const rect = document.querySelector("#threeMount canvas").getBoundingClientRect();
    return { x: rect.left + rect.width * 0.5, y: rect.top + rect.height * 0.5 };
  });
  await page.mouse.click(selection.x, selection.y);
  await page.waitForTimeout(100);
  const afterClick = await readMetrics(page);
  assert.notStrictEqual(afterClick.selectedStatus, "no cell selected", `${name} clicking near a cell should select it`);

  const frameIsolateEnabled = await page.evaluate(() => ({
    frame: !document.getElementById("frameCellBtn").disabled,
    isolate: !document.getElementById("isolateCellBtn").disabled,
  }));
  assert.ok(frameIsolateEnabled.frame, `${name} Frame Cell should be enabled once a cell is selected`);
  assert.ok(frameIsolateEnabled.isolate, `${name} Isolate Cell should be enabled once a cell is selected`);

  // Frame Cell should not throw (camera-target state is internal to the module, not DOM-observable).
  await page.click("#frameCellBtn");
  await page.waitForTimeout(50);

  // Isolate Cell should visibly reduce rendered cell pixels (only the selected cell stays visible)
  // and flip its own label to "Show Lattice".
  const beforeIsolateBuffer = await page.locator("#threeMount").screenshot();
  const beforeIsolatePixels = colorCount(beforeIsolateBuffer, (r, g, b) => r + g + b < 720);
  await page.click("#isolateCellBtn");
  await page.waitForTimeout(150);
  const isolateLabel = await page.evaluate(() => document.getElementById("isolateCellBtn").textContent);
  assert.strictEqual(isolateLabel, "Show Lattice", `${name} Isolate Cell button should relabel to Show Lattice while active`);
  const afterIsolateBuffer = await page.locator("#threeMount").screenshot();
  const afterIsolatePixels = colorCount(afterIsolateBuffer, (r, g, b) => r + g + b < 720);
  assert.ok(afterIsolatePixels < beforeIsolatePixels * 0.5, `${name} isolating a cell should render far fewer non-background pixels`);
  await page.click("#isolateCellBtn");
  await page.waitForTimeout(50);
  const restoredLabel = await page.evaluate(() => document.getElementById("isolateCellBtn").textContent);
  assert.strictEqual(restoredLabel, "Isolate Cell", `${name} clicking Isolate Cell again should restore the full lattice`);

  // Focus mode should hide the control panel and relabel the toggle to Controls.
  await page.click("#focusToggle");
  await page.waitForTimeout(50);
  const focusState = await page.evaluate(() => ({
    bodyHasClass: document.body.classList.contains("focus-mode"),
    panelDisplay: getComputedStyle(document.querySelector(".control-panel")).display,
    label: document.getElementById("focusToggle").textContent,
  }));
  assert.ok(focusState.bodyHasClass, `${name} Focus should add the focus-mode class`);
  assert.strictEqual(focusState.panelDisplay, "none", `${name} Focus should hide the control panel`);
  assert.strictEqual(focusState.label, "Controls", `${name} Focus button should relabel to Controls while active`);
  await page.click("#focusToggle");
  await page.waitForTimeout(50);
  const unfocusState = await page.evaluate(() => ({
    bodyHasClass: document.body.classList.contains("focus-mode"),
    label: document.getElementById("focusToggle").textContent,
  }));
  assert.ok(!unfocusState.bodyHasClass, `${name} clicking Focus again should remove focus-mode`);
  assert.strictEqual(unfocusState.label, "Focus", `${name} button should relabel back to Focus`);

  // Per-cell actuation: making the selected cell an actuator at max alpha
  // should (a) show up in the actuator count, (b) break the exact ring
  // closure that only holds when every cell shares one alpha, and (c) pull
  // a neighboring free cell's alpha away from the shared baseline through
  // the backlash-gated propagation - not just the actuated cell itself.
  const beforeActuation = await readMetrics(page);
  assert.strictEqual(beforeActuation.selectedStatus === "no cell selected" ? "no" : "yes", "yes", `${name} a cell should still be selected here`);
  const preActuationAlpha = await page.evaluate(() => document.getElementById("selectedCellAlphaMetric").textContent);
  await page.evaluate(() => {
    const roleSelect = document.getElementById("cellRoleSelect");
    roleSelect.value = "actuator";
    roleSelect.dispatchEvent(new Event("change", { bubbles: true }));
    const cellAlpha = document.getElementById("cellAlpha");
    cellAlpha.value = cellAlpha.max;
    cellAlpha.dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.waitForTimeout(150);
  const afterActuation = await readMetrics(page);
  assert.strictEqual(afterActuation.actuatorCount, "1", `${name} setting a cell's role to actuator should count it`);
  assert.ok(Number(afterActuation.selectedCellAlpha) > Number(preActuationAlpha), `${name} the actuated cell's alpha should jump toward its commanded max`);
  assert.ok(Math.abs(mmValue(afterActuation.closure)) > 0.01, `${name} an actuated cell should generally break exact ring closure (non-uniform per-cell alpha)`);

  // Read the full per-cell alpha grid through the debug hook (window.__cylinderTilingDebug,
  // exposed specifically because guessing screen coordinates to click a
  // particular (row, i) cell is fragile across camera angles/viewports) and
  // confirm propagation reached the actuated cell's circumferential
  // neighbor, not just the actuated cell itself.
  const propagation = await page.evaluate(() => {
    const debugHook = window.__cylinderTilingDebug;
    const selectedCell = debugHook.getSelected();
    const alphas = debugHook.getCellAlphas();
    const roles = debugHook.getCellRoles();
    if (!selectedCell || !alphas) return null;
    const n = alphas[selectedCell.row].length;
    const neighborIndex = (selectedCell.i + 1) % n;
    return {
      selectedAlpha: alphas[selectedCell.row][selectedCell.i],
      neighborAlpha: alphas[selectedCell.row][neighborIndex],
      neighborRole: roles[selectedCell.row][neighborIndex].role,
    };
  });
  assert.ok(propagation, `${name} debug hook should report the selected cell and alpha grid`);
  assert.ok(propagation.selectedAlpha >= 1.99, `${name} the actuated cell's own alpha should be at its commanded max (~2.0)`);
  assert.strictEqual(propagation.neighborRole, "free", `${name} the actuated cell's circumferential neighbor should still be "free"`);
  assert.ok(
    Math.abs(propagation.neighborAlpha - 1.2) > 0.01,
    `${name} the actuated cell's free neighbor should shift away from the 1.20 no-actuator baseline (got ${propagation.neighborAlpha})`
  );

  // Alpha heatmap: read the actuated cell's actual material color through
  // the debug hook (precise; a screenshot pixel-color approach turned out
  // to be too fragile here, since lighting/shading shifts rendered pixels
  // away from the material's raw hex color, and several row-palette colors
  // incidentally fall inside any reasonably wide "reddish" band). With the
  // actuator still at alpha~2.0, its material color should be near the
  // heatmap's "expanded" red once the toggle is on, and back to its row
  // palette color once off.
  const colorBeforeHeatmap = await page.evaluate(() => {
    const selectedCell = window.__cylinderTilingDebug.getSelected();
    return window.__cylinderTilingDebug.getCellTopColor(selectedCell.row, selectedCell.i);
  });
  assert.notStrictEqual(colorBeforeHeatmap, "#d63a2f", `${name} the actuated cell should not already be heatmap-red before the toggle is on`);
  await page.evaluate(() => {
    const heatmapInput = document.getElementById("heatmapEnabled");
    heatmapInput.checked = true;
    heatmapInput.dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.waitForTimeout(150);
  const colorDuringHeatmap = await page.evaluate(() => {
    const selectedCell = window.__cylinderTilingDebug.getSelected();
    return window.__cylinderTilingDebug.getCellTopColor(selectedCell.row, selectedCell.i);
  });
  assert.strictEqual(colorDuringHeatmap, "#d63a2f", `${name} an alpha~2.0 (max) cell should render at the heatmap's "expanded" red once enabled`);
  await page.evaluate(() => {
    const heatmapInput = document.getElementById("heatmapEnabled");
    heatmapInput.checked = false;
    heatmapInput.dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.waitForTimeout(100);
  const colorAfterHeatmap = await page.evaluate(() => {
    const selectedCell = window.__cylinderTilingDebug.getSelected();
    return window.__cylinderTilingDebug.getCellTopColor(selectedCell.row, selectedCell.i);
  });
  assert.strictEqual(colorAfterHeatmap, colorBeforeHeatmap, `${name} disabling the heatmap should restore the cell's original row-palette color`);

  // Clear All Actuators/Locks should reset the ring back to the uniform,
  // exactly-closed baseline state.
  await page.click("#clearRolesBtn");
  await page.waitForTimeout(150);
  const afterClear = await readMetrics(page);
  assert.strictEqual(afterClear.actuatorCount, "0", `${name} Clear All should remove every actuator`);
  assert.strictEqual(afterClear.lockedCount, "0", `${name} Clear All should remove every lock`);
  assert.ok(Math.abs(mmValue(afterClear.closure)) <= 0.01, `${name} clearing all roles should restore exact ring closure`);

  // Re-actuate once more so the JSON save/load round-trip below has
  // non-trivial per-cell state to carry through.
  await page.mouse.click(selection.x, selection.y);
  await page.waitForTimeout(100);
  await page.evaluate(() => {
    const roleSelect = document.getElementById("cellRoleSelect");
    roleSelect.value = "locked";
    roleSelect.dispatchEvent(new Event("change", { bubbles: true }));
  });
  await page.waitForTimeout(100);
  const beforeSaveRoles = await readMetrics(page);
  assert.strictEqual(beforeSaveRoles.lockedCount, "1", `${name} locking the selected cell should count it before saving`);

  // Export OBJ: the downloaded file should have one "o" group per visible
  // cell and internally consistent face indices (every face index must
  // reference a vertex that was actually written, and in range).
  const objDownloadPromise = page.waitForEvent("download");
  await page.click("#exportObjBtn");
  const objDownload = await objDownloadPromise;
  const objPath = path.join(root, `cylinder-tiling-${name}.obj`);
  await objDownload.saveAs(objPath);
  const objContent = fs.readFileSync(objPath, "utf8");
  const objLines = objContent.split("\n");
  const vertexCount = objLines.filter((line) => line.startsWith("v ")).length;
  const faceLines = objLines.filter((line) => line.startsWith("f "));
  const objectCount = objLines.filter((line) => line.startsWith("o ")).length;
  const expectedCellCount = Number(defaults.totalCells);
  assert.strictEqual(objectCount, expectedCellCount, `${name} OBJ export should have one object per visible cell (n*m=${expectedCellCount})`);
  assert.ok(vertexCount > 0, `${name} OBJ export should contain vertices`);
  assert.ok(faceLines.length > 0, `${name} OBJ export should contain faces`);
  let maxFaceIndex = 0;
  let minFaceIndex = Infinity;
  faceLines.forEach((line) => {
    line
      .slice(2)
      .trim()
      .split(/\s+/)
      .forEach((token) => {
        const index = parseInt(token, 10);
        maxFaceIndex = Math.max(maxFaceIndex, index);
        minFaceIndex = Math.min(minFaceIndex, index);
      });
  });
  assert.ok(minFaceIndex >= 1, `${name} OBJ face indices should be 1-based (got min ${minFaceIndex})`);
  assert.ok(maxFaceIndex <= vertexCount, `${name} OBJ face indices (max ${maxFaceIndex}) should not exceed the vertex count (${vertexCount})`);
  fs.unlinkSync(objPath);

  // Save/Load JSON: round-trip a changed alpha through a downloaded file.
  const downloadPromise = page.waitForEvent("download");
  await page.evaluate(() => {
    document.getElementById("alpha").value = "1.7";
    document.getElementById("alpha").dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.click("#saveJsonBtn");
  const download = await downloadPromise;
  const savedPath = path.join(root, `cylinder-tiling-${name}-saved.json`);
  await download.saveAs(savedPath);
  const saved = require(savedPath);
  assert.strictEqual(saved.format, "rad-cylinder-tiling.v1", `${name} saved JSON should carry the expected format tag`);
  assert.strictEqual(saved.alpha, 1.7, `${name} saved JSON should capture the current alpha`);
  assert.ok(Array.isArray(saved.cellRoles), `${name} saved JSON should include per-cell roles`);
  const savedLockedCount = saved.cellRoles.flat().filter((c) => c && c.role === "locked").length;
  assert.strictEqual(savedLockedCount, 1, `${name} saved JSON should capture the locked cell`);

  await page.evaluate(() => {
    document.getElementById("alpha").value = "0.4";
    document.getElementById("alpha").dispatchEvent(new Event("input", { bubbles: true }));
    document.getElementById("clearRolesBtn").click();
  });
  await page.waitForTimeout(100);
  const [fileChooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("#loadJsonBtn")]);
  await fileChooser.setFiles(savedPath);
  await page.waitForTimeout(100);
  const afterLoadRoles = await readMetrics(page);
  assert.strictEqual(afterLoadRoles.lockedCount, "1", `${name} loading the saved JSON should restore the locked cell`);
  const loadedAlpha = await page.evaluate(() => document.getElementById("alpha").value);
  assert.strictEqual(loadedAlpha, "1.7", `${name} loading the saved JSON should restore alpha=1.7`);

  // Target Diameter: fitting should set an alpha whose *actual rendered*
  // diameter (not just the fit tool's own estimate) matches the achieved
  // value it reports - this is a real regression check: an earlier version
  // computed the fit ignoring the backlash dead-zone the render pipeline
  // applies, so the alpha it picked rendered at a different diameter than
  // what the fit claimed.
  // Dragging the slider alone (no button click) must apply the fit live,
  // matching every other control in the app - it previously only updated
  // its own label until "Fit Alpha" was separately clicked, which read as
  // "target diameter doesn't do anything" since nothing else here requires
  // a second step.
  const alphaBeforeDrag = await page.evaluate(() => document.getElementById("alpha").value);
  await page.evaluate(() => {
    document.getElementById("targetDiameter").value = "125";
    document.getElementById("targetDiameter").dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.waitForTimeout(150);
  const afterDragOnly = await readMetrics(page);
  assert.notStrictEqual(afterDragOnly.alpha, alphaBeforeDrag, `${name} dragging Target Diameter alone (no button click) should already change alpha`);
  assert.strictEqual(mmValue(afterDragOnly.diameter).toFixed(1), "125.0", `${name} dragging Target Diameter alone should already re-render the ring at the new target`);

  await page.evaluate(() => {
    document.getElementById("targetDiameter").value = "140";
    document.getElementById("targetDiameter").dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.click("#fitDiameterBtn");
  await page.waitForTimeout(150);
  const fitResult = await readMetrics(page);
  const fitAchieved = await page.evaluate(() => document.getElementById("achievedDiameterMetric").textContent);
  assert.ok(fitAchieved.endsWith("mm"), `${name} fit should report an achieved diameter`);
  assert.strictEqual(mmValue(fitAchieved).toFixed(1), mmValue(fitResult.diameter).toFixed(1), `${name} the fit's reported achieved diameter should match the actual rendered ring diameter`);
  assert.strictEqual(fitResult.actuatorCount, "0", `${name} fitting should clear any actuators first (targets one shared alpha)`);
  assert.strictEqual(fitResult.lockedCount, "0", `${name} fitting should clear any locks first`);

  // An unreachable target (above the achievable range) should still report
  // an honest, non-hidden residual rather than silently clamping.
  await page.evaluate(() => {
    document.getElementById("targetDiameter").value = "259";
    document.getElementById("targetDiameter").dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.click("#fitDiameterBtn");
  await page.waitForTimeout(150);
  const outOfRangeResidual = await page.evaluate(() => document.getElementById("fitResidualMetric").textContent);
  assert.ok(mmValue(outOfRangeResidual) < -1, `${name} an unreachable target diameter should report a clearly nonzero residual instead of pretending it matched`);

  // Timeline: capture two keyframes at different alphas, jump between them,
  // play the sequence end to end, then save/load/delete.
  await page.evaluate(() => {
    document.getElementById("alpha").value = "1.3";
    document.getElementById("alpha").dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.click("#captureKeyframeBtn");
  await page.waitForTimeout(80);
  await page.evaluate(() => {
    document.getElementById("alpha").value = "1.9";
    document.getElementById("alpha").dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.click("#captureKeyframeBtn");
  await page.waitForTimeout(80);
  const afterTwoCaptures = await page.evaluate(() => ({
    count: document.getElementById("keyframeCountMetric").textContent,
    rowCount: document.querySelectorAll(".keyframe-row").length,
    playDisabled: document.getElementById("playSequenceBtn").disabled,
    saveDisabled: document.getElementById("saveSequenceBtn").disabled,
  }));
  assert.strictEqual(afterTwoCaptures.count, "2", `${name} capturing two keyframes should count two`);
  assert.strictEqual(afterTwoCaptures.rowCount, 2, `${name} should render two keyframe rows`);
  assert.ok(!afterTwoCaptures.playDisabled, `${name} Play Sequence should be enabled once keyframes exist`);
  assert.ok(!afterTwoCaptures.saveDisabled, `${name} Save Sequence JSON should be enabled once keyframes exist`);

  await page.click(".keyframe-row:first-child button:nth-of-type(1)");
  await page.waitForTimeout(100);
  const afterGoToFirst = await page.evaluate(() => document.getElementById("alpha").value);
  assert.strictEqual(afterGoToFirst, "1.3", `${name} jumping to the first keyframe should restore alpha=1.3`);

  await page.click("#playSequenceBtn");
  await page.waitForTimeout(100);
  const duringPlayback = await page.evaluate(() => document.getElementById("playSequenceBtn").textContent);
  assert.strictEqual(duringPlayback, "Stop", `${name} Play Sequence should relabel to Stop while playing`);
  await page.waitForTimeout(2900); // two keyframes at a 1200ms step interval, plus margin
  const afterPlayback = await page.evaluate(() => ({
    label: document.getElementById("playSequenceBtn").textContent,
    alpha: document.getElementById("alpha").value,
  }));
  assert.strictEqual(afterPlayback.label, "Play Sequence", `${name} playback should stop itself and relabel after the last keyframe`);
  assert.strictEqual(afterPlayback.alpha, "1.9", `${name} playback should end on the last keyframe's alpha`);

  const sequenceDownloadPromise = page.waitForEvent("download");
  await page.click("#saveSequenceBtn");
  const sequenceDownload = await sequenceDownloadPromise;
  const sequencePath = path.join(root, `cylinder-tiling-${name}-sequence.json`);
  await sequenceDownload.saveAs(sequencePath);
  const savedSequence = require(sequencePath);
  assert.strictEqual(savedSequence.format, "rad-cylinder-tiling-sequence.v1", `${name} saved sequence JSON should carry the expected format tag`);
  assert.strictEqual(savedSequence.keyframes.length, 2, `${name} saved sequence JSON should include both keyframes`);

  await page.click(".keyframe-row:first-child button:nth-of-type(2)"); // delete
  await page.waitForTimeout(80);
  const afterDelete = await page.evaluate(() => document.getElementById("keyframeCountMetric").textContent);
  assert.strictEqual(afterDelete, "1", `${name} deleting a keyframe should reduce the count`);

  const [sequenceFileChooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("#loadSequenceBtn")]);
  await sequenceFileChooser.setFiles(sequencePath);
  await page.waitForTimeout(100);
  const afterLoadSequence = await page.evaluate(() => document.getElementById("keyframeCountMetric").textContent);
  assert.strictEqual(afterLoadSequence, "2", `${name} loading a saved sequence should restore both keyframes`);

  // Measurement line: enabling it should populate real (non-degenerate)
  // dimension-line geometry - drawn offset below the part with extension
  // lines back to the true diameter endpoints, not straight through the
  // cell geometry (which turned out to render mostly occluded when tried
  // that way first).
  await page.evaluate(() => {
    const input = document.getElementById("showMeasurements");
    input.checked = true;
    input.dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.waitForTimeout(100);
  const measurementBox = await page.evaluate(() => {
    const canvas = document.querySelector("#threeMount canvas");
    return { w: canvas.width, h: canvas.height };
  });
  assert.ok(measurementBox.w > 0 && measurementBox.h > 0, `${name} canvas should still render with measurements enabled`);
  await page.evaluate(() => {
    const input = document.getElementById("showMeasurements");
    input.checked = false;
    input.dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.waitForTimeout(80);

  // Reset All: change a wide spread of settings (drive, ring/row counts,
  // an actuator, a selection), then confirm every one of them - and only
  // those, nothing left stale - lands back on its HTML-declared default.
  await page.evaluate(() => {
    document.getElementById("alpha").value = "1.9";
    document.getElementById("alpha").dispatchEvent(new Event("input", { bubbles: true }));
    document.getElementById("ringCount").value = "14";
    document.getElementById("ringCount").dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.waitForTimeout(100);
  await page.mouse.click(selection.x, selection.y);
  await page.waitForTimeout(100);
  await page.evaluate(() => {
    const roleSelect = document.getElementById("cellRoleSelect");
    if (roleSelect.disabled) return;
    roleSelect.value = "actuator";
    roleSelect.dispatchEvent(new Event("change", { bubbles: true }));
  });
  await page.waitForTimeout(100);
  await page.click("#resetAllBtn");
  await page.waitForTimeout(150);
  const afterResetAll = await readMetrics(page);
  assert.strictEqual(afterResetAll.alpha, "1.3", `${name} Reset All should restore the default alpha`);
  assert.strictEqual(afterResetAll.ringCount, "10", `${name} Reset All should restore the default ring count`);
  assert.strictEqual(afterResetAll.rowCount, "3", `${name} Reset All should restore the default row count`);
  assert.strictEqual(afterResetAll.actuatorCount, "0", `${name} Reset All should clear actuators`);
  assert.strictEqual(afterResetAll.lockedCount, "0", `${name} Reset All should clear locks`);
  assert.strictEqual(afterResetAll.selectedStatus, "no cell selected", `${name} Reset All should clear the current selection`);
  const afterResetAllSites = await page.evaluate(() => ({
    aSite: document.getElementById("aSite").value,
    bSite: document.getElementById("bSite").value,
  }));
  assert.strictEqual(afterResetAllSites.aSite, "east", `${name} Reset All should restore the default outgoing site`);
  assert.strictEqual(afterResetAllSites.bSite, "west", `${name} Reset All should restore the default incoming site`);

  assert.deepStrictEqual(errors, [], `${name} should not emit browser errors`);

  const screenshotPath = path.join(root, `cylinder-tiling-${name}-check.png`);
  const buffer = await page.locator("#threeMount").screenshot({ path: screenshotPath });
  const orangePixels = colorCount(buffer, (r, g, b) => r > 130 && g > 55 && g < 190 && b < 130);
  const slatePixels = colorCount(buffer, (r, g, b) => Math.abs(r - 68) < 30 && Math.abs(g - 81) < 30 && Math.abs(b - 95) < 30);
  assert.ok(orangePixels > 80, `${name} should render row 0's upper-cross color`);
  assert.ok(slatePixels > 40, `${name} should render row 0's lower-cross color`);

  await page.close();
}

// Regression test for a real bug: renderer.setSize(w, h, false) leaves the
// canvas's on-screen CSS size to fall back to its width/height attributes
// (the backing-buffer resolution, i.e. w*devicePixelRatio) whenever the
// display isn't at 100% OS scaling, unless the canvas's CSS size is pinned
// explicitly. Combined with resize() running every render frame, that
// mismatch compounded into a runaway feedback loop (canvas measuring
// itself too big, growing, being measured too big again, ...), observed
// growing to tens of millions of pixels within about a second on a
// 150%-scaled display. This checks the canvas's rendered CSS size stays
// locked to its container at several non-100% scale factors, across
// several seconds of real frames, instead of drifting or exploding.
async function checkDeviceScaleFactor(browser, scaleFactor) {
  const page = await browser.newPage({ viewport: { width: 1000, height: 700 }, deviceScaleFactor: scaleFactor });
  await page.goto(pageUrl);
  await page.waitForSelector("#threeMount canvas", { state: "visible" });
  await page.waitForTimeout(2500); // ~150 frames at 60fps - long enough for a feedback loop to blow up
  const sizes = await page.evaluate(() => {
    const canvas = document.querySelector("#threeMount canvas");
    const mount = document.getElementById("threeMount");
    const mountRect = mount.getBoundingClientRect();
    const canvasRect = canvas.getBoundingClientRect();
    return { mountRect: { w: mountRect.width, h: mountRect.height }, canvasRect: { w: canvasRect.width, h: canvasRect.height } };
  });
  const label = `deviceScaleFactor=${scaleFactor}`;
  assert.ok(sizes.canvasRect.w < 3000, `${label}: canvas CSS width should stay near its container (got ${sizes.canvasRect.w}px) - not run away`);
  assert.ok(sizes.canvasRect.h < 3000, `${label}: canvas CSS height should stay near its container (got ${sizes.canvasRect.h}px) - not run away`);
  assert.ok(
    Math.abs(sizes.canvasRect.w - sizes.mountRect.w) <= 2,
    `${label}: canvas CSS width (${sizes.canvasRect.w}) should match its container (${sizes.mountRect.w})`
  );
  assert.ok(
    Math.abs(sizes.canvasRect.h - sizes.mountRect.h) <= 2,
    `${label}: canvas CSS height (${sizes.canvasRect.h}) should match its container (${sizes.mountRect.h})`
  );
  await page.close();
}

(async () => {
  const browser = await launchBrowser();
  try {
    await checkViewport(browser, "desktop", { width: 1280, height: 820 });
    await checkViewport(browser, "mobile", { width: 390, height: 860 });
    for (const scaleFactor of [1.25, 1.5, 2]) {
      await checkDeviceScaleFactor(browser, scaleFactor);
    }
  } finally {
    await browser.close();
  }
  console.log("cylinder tiling page validation passed");
})();
