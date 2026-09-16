const assert = require("assert");
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

  await page.evaluate(() => {
    document.getElementById("alpha").value = "0.4";
    document.getElementById("alpha").dispatchEvent(new Event("input", { bubbles: true }));
  });
  const [fileChooser] = await Promise.all([page.waitForEvent("filechooser"), page.click("#loadJsonBtn")]);
  await fileChooser.setFiles(savedPath);
  await page.waitForTimeout(100);
  const loadedAlpha = await page.evaluate(() => document.getElementById("alpha").value);
  assert.strictEqual(loadedAlpha, "1.7", `${name} loading the saved JSON should restore alpha=1.7`);

  assert.deepStrictEqual(errors, [], `${name} should not emit browser errors`);

  const screenshotPath = path.join(root, `cylinder-tiling-${name}-check.png`);
  const buffer = await page.locator("#threeMount").screenshot({ path: screenshotPath });
  const orangePixels = colorCount(buffer, (r, g, b) => r > 130 && g > 55 && g < 190 && b < 130);
  const slatePixels = colorCount(buffer, (r, g, b) => Math.abs(r - 68) < 30 && Math.abs(g - 81) < 30 && Math.abs(b - 95) < 30);
  assert.ok(orangePixels > 80, `${name} should render row 0's upper-cross color`);
  assert.ok(slatePixels > 40, `${name} should render row 0's lower-cross color`);

  await page.close();
}

(async () => {
  const browser = await launchBrowser();
  try {
    await checkViewport(browser, "desktop", { width: 1280, height: 820 });
    await checkViewport(browser, "mobile", { width: 390, height: 860 });
  } finally {
    await browser.close();
  }
  console.log("cylinder tiling page validation passed");
})();
