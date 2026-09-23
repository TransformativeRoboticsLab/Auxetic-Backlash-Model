const assert = require("assert");
const path = require("path");
const { pathToFileURL } = require("url");
const { PNG } = require("pngjs");
const { chromium } = require("playwright");

const root = path.resolve(__dirname, "..");
const pageUrl = pathToFileURL(path.join(root, "web", "two-cell-attachment", "index.html")).href;

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

function parseCenter(value) {
  const match = value.match(/^\((-?\d+(?:\.\d+)?), (-?\d+(?:\.\d+)?)\)$/);
  assert.ok(match, `could not parse center metric: ${value}`);
  return [Number(match[1]), Number(match[2])];
}

async function launchBrowser() {
  try {
    return await chromium.launch({ channel: "chrome" });
  } catch (_error) {
    return chromium.launch();
  }
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

  const metrics = await page.evaluate(() => {
    const canvas = document.querySelector("#threeMount canvas");
    return {
      canvasWidth: canvas.width,
      canvasHeight: canvas.height,
      driveRangeIds: Array.from(document.querySelectorAll(".control-panel section:first-child input[type='range']")).map((input) => input.id),
      aTop: document.getElementById("aTopAngle").value,
      residual: document.getElementById("residualMetric").textContent,
      secondResidual: document.getElementById("secondResidualMetric").textContent,
      pitch: document.getElementById("pitchMetric").textContent,
      aCenter: document.getElementById("aCenterMetric").textContent,
      bCenter: document.getElementById("bCenterMetric").textContent,
      cCenter: document.getElementById("cCenterMetric").textContent,
      systemCenter: document.getElementById("systemCenterMetric").textContent,
      secondTarget: document.getElementById("secondTargetMetric").textContent,
      secondB: document.getElementById("secondBMetric").textContent,
      bcPrimaryResidual: document.getElementById("bcPrimaryResidualMetric").textContent,
      bcSecondResidual: document.getElementById("bcSecondResidualMetric").textContent,
      bcPitch: document.getElementById("bcPitchMetric").textContent,
      collisionState: document.getElementById("collisionState").textContent,
      allowedDrive: document.getElementById("allowedDriveMetric").textContent,
      driveClamp: document.getElementById("driveClampMetric").textContent,
      penetration: document.getElementById("penetrationMetric").textContent,
      clearance: document.getElementById("clearanceMetric").textContent,
      overflow: document.documentElement.scrollWidth - document.documentElement.clientWidth,
    };
  });
  assert.ok(metrics.canvasWidth > 300, `${name} canvas should have real width`);
  assert.ok(metrics.canvasHeight > 300, `${name} canvas should have real height`);
  assert.deepStrictEqual(metrics.driveRangeIds, ["aTopAngle"], `${name} should expose one uniform rotational drive slider`);
  assert.strictEqual(metrics.residual, "0.000 mm", `${name} shared-pin residual should be zero`);
  assert.strictEqual(metrics.secondResidual, "0.000 mm", `${name} second shared-pin residual should be zero`);
  assert.strictEqual(metrics.bcPrimaryResidual, "0.000 mm", `${name} B-C primary residual should be zero`);
  assert.strictEqual(metrics.bcSecondResidual, "0.000 mm", `${name} B-C second residual should be zero`);
  assert.notStrictEqual(metrics.collisionState, "blocked", `${name} default pose should be physically feasible`);
  assert.ok(metrics.allowedDrive.includes("deg"), `${name} should report the feasible drive angle range`);
  assert.notStrictEqual(metrics.driveClamp, "no feasible range", `${name} should have a feasible drive clamp state`);
  assert.strictEqual(metrics.penetration, "0.000 mm", `${name} rendered pose should have no material penetration`);
  assert.ok(metrics.clearance.endsWith("mm"), `${name} should report material clearance`);
  assert.ok(metrics.pitch.endsWith("mm"), `${name} should report center pitch`);
  assert.ok(metrics.aCenter.startsWith("("), `${name} should report Cell A moving center`);
  assert.ok(metrics.bCenter.startsWith("("), `${name} should report Cell B center`);
  assert.ok(metrics.cCenter.startsWith("("), `${name} should report Cell C center`);
  assert.strictEqual(metrics.systemCenter, "(0.0, 0.0)", `${name} should keep the three-cell system center fixed`);
  const aCenter = parseCenter(metrics.aCenter);
  const bCenter = parseCenter(metrics.bCenter);
  const cCenter = parseCenter(metrics.cCenter);
  const abVector = [bCenter[0] - aCenter[0], bCenter[1] - aCenter[1]];
  const bcVector = [cCenter[0] - bCenter[0], cCenter[1] - bCenter[1]];
  assert.ok(Math.abs(abVector[0] - bcVector[0]) <= 0.2, `${name} A-B and B-C center steps should match in x`);
  assert.ok(Math.abs(abVector[1] - bcVector[1]) <= 0.2, `${name} A-B and B-C center steps should match in y`);
  assert.ok(Math.abs(aCenter[0] + bCenter[0] + cCenter[0]) <= 0.2, `${name} three-cell center of mass should stay fixed in x`);
  assert.ok(Math.abs(aCenter[1] + bCenter[1] + cCenter[1]) <= 0.2, `${name} three-cell center of mass should stay fixed in y`);
  assert.strictEqual(metrics.pitch, metrics.bcPitch, `${name} A-B and B-C pitch readouts should match`);
  assert.ok(metrics.secondTarget.startsWith("("), `${name} should report the second pin target`);
  assert.ok(metrics.secondB.startsWith("("), `${name} should report the second pin on Cell B`);
  assert.ok(metrics.overflow <= 2, `${name} should not create horizontal overflow`);
  assert.deepStrictEqual(errors, [], `${name} should not emit browser errors`);

  const contractedMetrics = await page.evaluate(() => {
    const input = document.getElementById("aTopAngle");
    input.value = input.min;
    input.dispatchEvent(new Event("input", { bubbles: true }));
    return {
      aCenter: document.getElementById("aCenterMetric").textContent,
      bCenter: document.getElementById("bCenterMetric").textContent,
      cCenter: document.getElementById("cCenterMetric").textContent,
      systemCenter: document.getElementById("systemCenterMetric").textContent,
      residual: document.getElementById("residualMetric").textContent,
      secondResidual: document.getElementById("secondResidualMetric").textContent,
      bcPrimaryResidual: document.getElementById("bcPrimaryResidualMetric").textContent,
      bcSecondResidual: document.getElementById("bcSecondResidualMetric").textContent,
      penetration: document.getElementById("penetrationMetric").textContent,
    };
  });
  const contractedA = parseCenter(contractedMetrics.aCenter);
  const contractedB = parseCenter(contractedMetrics.bCenter);
  const contractedC = parseCenter(contractedMetrics.cCenter);
  const centerlineCross = cCenter[0] * contractedC[1] - cCenter[1] * contractedC[0];
  assert.ok(Math.abs(centerlineCross) <= 2.0, `${name} cell centers should stay on one straight centerline while contracting`);
  assert.ok(Math.abs(contractedB[0]) <= 0.15, `${name} contracted middle center should stay on the system center in x`);
  assert.ok(Math.abs(contractedB[1]) <= 0.15, `${name} contracted middle center should stay on the system center in y`);
  assert.ok(Math.abs(contractedA[0] + contractedB[0] + contractedC[0]) <= 0.2, `${name} contracted center of mass should remain fixed in x`);
  assert.ok(Math.abs(contractedA[1] + contractedB[1] + contractedC[1]) <= 0.2, `${name} contracted center of mass should remain fixed in y`);
  assert.strictEqual(contractedMetrics.systemCenter, "(0.0, 0.0)", `${name} contracted system center should remain fixed`);
  assert.strictEqual(contractedMetrics.residual, "0.000 mm", `${name} contracted primary pin residual should remain zero`);
  assert.strictEqual(contractedMetrics.secondResidual, "0.000 mm", `${name} contracted second pin residual should remain zero`);
  assert.strictEqual(contractedMetrics.bcPrimaryResidual, "0.000 mm", `${name} contracted B-C primary residual should remain zero`);
  assert.strictEqual(contractedMetrics.bcSecondResidual, "0.000 mm", `${name} contracted B-C second residual should remain zero`);
  assert.strictEqual(contractedMetrics.penetration, "0.000 mm", `${name} contracted pose should have no material penetration`);
  await page.evaluate((value) => {
    const input = document.getElementById("aTopAngle");
    input.value = value;
    input.dispatchEvent(new Event("input", { bubbles: true }));
  }, metrics.aTop || "64");

  const screenshotPath = path.join(root, `two-cell-attachment-${name}-check.png`);
  const buffer = await page.locator("#threeMount").screenshot({ path: screenshotPath });
  const orangePixels = colorCount(buffer, (r, g, b) => r > 130 && g > 55 && g < 190 && b < 130);
  const bluePixels = colorCount(buffer, (r, g, b) => b > 120 && r < 100 && g > 70 && g < 160);
  const purplePixels = colorCount(buffer, (r, g, b) => b > 120 && r > 80 && r < 170 && g > 50 && g < 140);
  const greenPixels = colorCount(buffer, (r, g, b) => g > 95 && r < 120 && b > 70 && b < 170);
  const pinPixels = colorCount(buffer, (r, g, b) => g > 100 && r < 85 && b < 130);
  assert.ok(orangePixels > 100, `${name} should render Cell A upper cross`);
  assert.ok(bluePixels > 100, `${name} should render Cell B upper cross`);
  assert.ok(purplePixels > 100, `${name} should render Cell B lower cross`);
  assert.ok(greenPixels > 100, `${name} should render Cell C upper cross`);
  assert.ok(pinPixels > 80, `${name} should render the shared pins`);
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
  console.log("two-cell attachment page validation passed");
})();
