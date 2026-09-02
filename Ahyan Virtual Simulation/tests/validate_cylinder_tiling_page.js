const assert = require("assert");
const path = require("path");
const { pathToFileURL } = require("url");
const { PNG } = require("pngjs");
const { chromium } = require("playwright");

const root = path.resolve(__dirname, "..");
const pageUrl = pathToFileURL(path.join(root, "web", "cylinder-tiling", "index.html")).href;

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
    aTop: document.getElementById("aTopAngle").value,
    ringCount: document.getElementById("ringCount").value,
    rowCount: document.getElementById("rowCount").value,
    pitch: document.getElementById("pitchMetric").textContent,
    turn: document.getElementById("turnMetric").textContent,
    diameter: document.getElementById("diameterMetric").textContent,
    radius: document.getElementById("radiusMetric").textContent,
    closure: document.getElementById("closureMetric").textContent,
    totalCells: document.getElementById("totalCellsMetric").textContent,
    height: document.getElementById("heightMetric").textContent,
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
  assert.ok(Math.abs(mmValue(defaults.closure)) <= 0.01, `${name} default ring closure residual should be ~0`);
  const defaultDiameter = mmValue(defaults.diameter);
  const defaultRadius = mmValue(defaults.radius);
  assert.ok(Math.abs(defaultDiameter - 2 * defaultRadius) <= 0.2, `${name} diameter should equal 2x radius`);
  assert.ok(defaultDiameter > 0, `${name} default diameter should be positive`);
  const turnDeg = Number(defaults.turn.replace(/ ?deg$/, ""));
  assert.ok(Math.abs(turnDeg - 36) <= 0.1, `${name} turn angle for n=10 should be 36 deg (360/10)`);
  const defaultHeight = mmValue(defaults.height);
  assert.ok(Math.abs(defaultHeight - 2 * 50) <= 0.5, `${name} default assembly height should be (m-1)*axialPitch`);

  // Sweep the drive angle across its full range; the ring should stay closed at every setting.
  for (const which of ["min", "max"]) {
    const swept = await page.evaluate((minOrMax) => {
      const input = document.getElementById("aTopAngle");
      input.value = input[minOrMax];
      input.dispatchEvent(new Event("input", { bubbles: true }));
      return document.getElementById("closureMetric").textContent;
    }, which);
    assert.ok(Math.abs(mmValue(swept)) <= 0.01, `${name} ring closure residual should stay ~0 at drive ${which}`);
  }
  await page.evaluate((value) => {
    const input = document.getElementById("aTopAngle");
    input.value = value;
    input.dispatchEvent(new Event("input", { bubbles: true }));
  }, defaults.aTop);

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
