const assert = require("assert");
const fs = require("fs");
const path = require("path");
const { pathToFileURL } = require("url");
const { PNG } = require("pngjs");
const { chromium } = require("playwright");

const root = path.resolve(__dirname, "..");
const pageUrl = pathToFileURL(path.join(root, "web", "one-cell-rotation", "index.html")).href;

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
    const mount = document.getElementById("threeMount");
    return {
      canvasWidth: canvas.width,
      canvasHeight: canvas.height,
      mountWidth: mount.clientWidth,
      mountHeight: mount.clientHeight,
      thetaText: document.getElementById("thetaMetric").textContent,
      bodyOverflow: document.documentElement.scrollWidth - document.documentElement.clientWidth,
    };
  });
  assert.ok(metrics.canvasWidth > 300, `${name} canvas should have real width`);
  assert.ok(metrics.canvasHeight > 300, `${name} canvas should have real height`);
  assert.ok(metrics.thetaText.includes("35"), `${name} should show the default driven angle`);
  assert.ok(metrics.bodyOverflow <= 2, `${name} should not create horizontal overflow`);
  assert.deepStrictEqual(errors, [], `${name} should not emit browser errors`);

  const screenshotPath = path.join(root, `one-cell-rotation-${name}-check.png`);
  const buffer = await page.locator("#threeMount").screenshot({ path: screenshotPath });
  const orangePixels = colorCount(buffer, (r, g, b) => r > 130 && g > 55 && g < 180 && b < 120);
  const darkPixels = colorCount(buffer, (r, g, b) => r < 105 && g < 115 && b < 130);
  assert.ok(orangePixels > 100, `${name} should render the driven upper cross`);
  assert.ok(darkPixels > 100, `${name} should render the fixed lower cross or hole geometry`);
  await page.close();
}

(async () => {
  const browser = await launchBrowser();
  try {
    await checkViewport(browser, "desktop", { width: 1280, height: 800 });
    await checkViewport(browser, "mobile", { width: 390, height: 840 });
  } finally {
    await browser.close();
  }
  console.log("one-cell rotation page validation passed");
})();
