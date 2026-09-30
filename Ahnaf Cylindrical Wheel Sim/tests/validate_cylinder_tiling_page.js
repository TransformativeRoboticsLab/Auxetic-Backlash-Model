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
  // Default 1.4 (effective 1.30, theta 31 deg) is inside the collision-free
  // window; the earlier 1.3 (effective 1.20) overlaps neighboring pads.
  assert.strictEqual(defaults.alpha, "1.4", `${name} default alpha should be 1.4`);
  assert.strictEqual(defaults.alphaCommand, "1.40", `${name} commanded alpha readout should mirror the slider`);
  // Dead-zone shifts the engaged effective alpha down by the gap width:
  // effective = 1.0 + reluDeadzone(alpha-1.0, backlash) = 1.0 + (0.4-0.1) = 1.3
  assert.strictEqual(defaults.alphaEffective, "1.30", `${name} at backlash=0.1 and |alpha-1|=0.4>0.1, effective alpha should be shifted by the gap width to 1.30`);
  assert.strictEqual(defaults.deadzone, "engaged", `${name} default alpha=1.4 should be outside the b=0.1 dead zone around 1.0`);
  const theta = Number(defaults.theta.replace(/ ?deg$/, ""));
  assert.ok(Math.abs(theta - (70 * 1.3 - 60)) <= 0.05, `${name} theta should follow theta = 70*effectiveAlpha - 60`);
  assert.ok(Math.abs(mmValue(defaults.closure)) <= 0.01, `${name} default ring closure residual should be ~0`);
  const defaultDiameter = mmValue(defaults.diameter);
  const defaultRadius = mmValue(defaults.radius);
  assert.ok(Math.abs(defaultDiameter - 2 * defaultRadius) <= 0.2, `${name} diameter should equal 2x radius`);
  assert.ok(defaultDiameter > 0, `${name} default diameter should be positive`);
  const turnDeg = Number(defaults.turn.replace(/ ?deg$/, ""));
  assert.ok(Math.abs(turnDeg - 36) <= 0.1, `${name} turn angle for n=10 should be 36 deg (360/10)`);
  // Rows are pinned by default, so the row spacing is set by the pins:
  // 2L*cos(theta/2) = 2*22.451*cos(15.5 deg) = 43.3 mm at theta = 31 deg
  // (L from Jacob's RADs unit cell STL).
  const defaultHeight = mmValue(defaults.height);
  const pinnedPitch = 2 * 22.451 * Math.cos(((70 * 1.3 - 60) * Math.PI) / 360);
  assert.ok(Math.abs(defaultHeight - 2 * pinnedPitch) <= 0.2, `${name} default assembly height should be (m-1) pinned row pitches (got ${defaultHeight})`);

  // Radial cell orientation: each cell's thickness-axis normal (read via the
  // debug hook, which reports the cell group's local-Z direction in world
  // space) should point away from the ring's own centroid - not merely
  // agree with whatever angle the orientation code itself was fed. A prior
  // version fed the ring-closure walk's chain heading (bottomRot) into the
  // orientation code directly; that heading leads the true polar position
  // angle by a construction-dependent phase (not a fixed 90deg), so cells
  // rendered rotated off-axis (arms pointing radially like spokes, normals
  // pointing tangentially) despite every self-consistency check on the
  // orientation code's own math passing. This test catches that class of
  // bug by comparing against the ring's actual center/centroid geometry,
  // an independent source of truth the orientation code does not feed from.
  const radialOrientation = await page.evaluate(() => {
    const debugHook = window.__cylinderTilingDebug;
    const ring = debugHook.getRingGeometry(0);
    const samples = [0, 2, 5];
    return samples.map((i) => {
      const normal = debugHook.getCellGroupNormal(0, i);
      const center = ring.centers[i];
      const expectedAngle = Math.atan2(center.y - ring.centroid.y, center.x - ring.centroid.x);
      const actualAngle = Math.atan2(normal.y, normal.x);
      let diffDeg = ((actualAngle - expectedAngle) * 180) / Math.PI;
      diffDeg = ((((diffDeg + 180) % 360) + 360) % 360) - 180;
      return { i, diffDeg };
    });
  });
  radialOrientation.forEach(({ i, diffDeg }) => {
    assert.ok(
      Math.abs(diffDeg) <= 1,
      `${name} cell ${i}'s orientation normal should point away from the ring's own centroid (off by ${diffDeg.toFixed(1)}deg)`
    );
  });

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
    document.getElementById("alpha").value = "1.4";
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

  // Pinned holes must actually line up. Read straight from the rendered
  // meshes (not from the pin-alignment math): at every joint, cell i's
  // upper east pad and cell i+1's lower west pad - and, with opposite sites,
  // cell i's lower east pad and cell i+1's upper west pad - must sit on one
  // pin axis (the bisector of the two cells' normals), spaced along it only
  // by the stacked plates: 4 mm * cos(bend/2). A regression here once left
  // the pinned pads 10.3 mm apart (9 mm of it along the axle) while an
  // earlier pin-gap readout, comparing the wrong pads, reported 1.9 mm.
  const jointGeometry = await page.evaluate(() => {
    const debugHook = window.__cylinderTilingDebug;
    const n = debugHook.getCellAlphas()[0].length;
    let worstLateral = 0;
    let worstAlong = 0;
    for (let i = 0; i < n; i += 1) {
      const next = (i + 1) % n;
      const nA = debugHook.getCellGroupNormal(0, i);
      const nB = debugHook.getCellGroupNormal(0, next);
      const axisLength = Math.hypot(nA.x + nB.x, nA.y + nB.y, nA.z + nB.z);
      const axis = { x: (nA.x + nB.x) / axisLength, y: (nA.y + nB.y) / axisLength, z: (nA.z + nB.z) / axisLength };
      const bendHalf = Math.acos(Math.min(1, nA.x * nB.x + nA.y * nB.y + nA.z * nB.z)) / 2;
      [
        ["upper", "lower"],
        ["lower", "upper"],
      ].forEach(([layerA, layerB]) => {
        const a = debugHook.getSiteWorld(0, i, layerA, "east");
        const b = debugHook.getSiteWorld(0, next, layerB, "west");
        const d = { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
        const along = d.x * axis.x + d.y * axis.y + d.z * axis.z;
        const lateral = Math.hypot(d.x - along * axis.x, d.y - along * axis.y, d.z - along * axis.z);
        worstLateral = Math.max(worstLateral, lateral);
        worstAlong = Math.max(worstAlong, Math.abs(Math.abs(along) - 4 * Math.cos(bendHalf)));
      });
    }
    return { worstLateral, worstAlong };
  });
  assert.ok(jointGeometry.worstLateral < 1e-6, `${name} pinned pads should share one pin axis (off by ${jointGeometry.worstLateral} mm)`);
  assert.ok(jointGeometry.worstAlong < 1e-6, `${name} pinned pads should be spaced only by the stacked plates along the pin (off by ${jointGeometry.worstAlong} mm)`);

  // The pin-alignment readout agrees (~0), and with rows pinned (default)
  // a cell's north pads meet the south pads of the cell above too.
  const pinAlignmentDefault = await page.evaluate(() => window.__cylinderTilingDebug.getPinAlignment());
  assert.ok(pinAlignmentDefault.maxCircumferential < 1e-6, `${name} pin-offset readout should be ~0 (got ${pinAlignmentDefault.maxCircumferential})`);
  assert.ok(pinAlignmentDefault.maxAxial < 1e-6, `${name} pinned rows' north/south pads should meet (got ${pinAlignmentDefault.maxAxial})`);

  // Unpinning the rows hands the spacing to the axial pitch slider: rows
  // become separate rings (no axial pins), and a large pitch opens a real
  // gap between them - the pads 70 - 42.6 = 27.4 mm apart along the axle.
  const setPinRows = (on) =>
    page.evaluate((value) => {
      const input = document.getElementById("pinRows");
      input.checked = value;
      input.dispatchEvent(new Event("change", { bubbles: true }));
    }, on);
  await setPinRows(false);
  await page.evaluate(() => {
    const input = document.getElementById("axialPitch");
    input.value = "70";
    input.dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.waitForTimeout(200);
  const unpinned = await page.evaluate(() => ({
    disabled: document.getElementById("axialPitch").disabled,
    height: document.getElementById("heightMetric").textContent,
    gap: window.__cylinderTilingDebug.getPinAlignment().maxAxial,
    pins: window.__cylinderTilingDebug.getPinInfo().visibleCount,
  }));
  assert.strictEqual(unpinned.disabled, false, `${name} the pitch slider should be usable once rows are unpinned`);
  assert.ok(Math.abs(mmValue(unpinned.height) - 140) < 0.2, `${name} unpinned height should follow the set pitch (got ${unpinned.height})`);
  assert.ok(Math.abs(unpinned.gap - (70 - pinnedPitch)) < 0.1, `${name} unpinned rows should show their real gap (got ${unpinned.gap})`);
  assert.strictEqual(unpinned.pins, 90, `${name} unpinned rows should have no axial pins (30 hub + 60 ring-joint pins)`);
  await page.evaluate(() => {
    const input = document.getElementById("axialPitch");
    input.value = input.defaultValue;
    input.dispatchEvent(new Event("input", { bubbles: true }));
  });
  await setPinRows(true);
  await page.waitForTimeout(200);
  const repinned = await page.evaluate(() => ({
    disabled: document.getElementById("axialPitch").disabled,
    gap: window.__cylinderTilingDebug.getPinAlignment().maxAxial,
  }));
  assert.ok(repinned.disabled && repinned.gap < 1e-6, `${name} re-pinning should rejoin the rows`);
  const circumferentialPinText = await page.evaluate(() => document.getElementById("circumferentialPinMetric").textContent);
  assert.ok(circumferentialPinText.endsWith("mm"), `${name} circumferential pin gap readout should be in mm`);

  // Backlash tilt check (exact pin-in-hole contact: theta_max =
  // acos(d/sqrt(D^2+t^2)) - atan(t/D), D=3.4, t=4.0, assumed d=3.0, doubled
  // for two plates on one floating bolt ~= 11.03deg). With radially facing
  // cells each circumferential joint bends 360/n = 36deg at n=10, which
  // backlash alone can't reach; the minimum ring is ceil(360/11.03) = 33.
  assert.ok(Math.abs((pinAlignmentDefault.tiltLimit * 180) / Math.PI - 11.03) < 0.05, `${name} joint tilt limit should be ~11.03deg`);
  assert.ok(Math.abs((pinAlignmentDefault.maxCircumferentialBend * 180) / Math.PI - 36) < 0.1, `${name} n=10 joints should each bend 36deg`);
  assert.strictEqual(pinAlignmentDefault.circumferentialOk, false, `${name} a 36deg joint bend should exceed the backlash tilt limit`);
  assert.strictEqual(pinAlignmentDefault.axialOk, true, `${name} uniform rows need no axial tilt`);
  assert.strictEqual(pinAlignmentDefault.minCellsForBacklashClosure, 33, `${name} minimum cells per ring for backlash-only closure`);
  const bendText = await page.evaluate(() => document.getElementById("circumferentialBendMetric").textContent);
  assert.ok(bendText.includes("exceeds backlash"), `${name} circumferential bend readout should flag exceeding backlash (got "${bendText}")`);

  // Constraints: no fusing through, no discontinuous motion. The default
  // pose sits inside the collision-free window; driving past either edge
  // holds the pose at contact instead of overlapping, and the pose can't
  // jump across a colliding band to the other, separate clear window.
  const readConstraint = () =>
    page.evaluate(() => {
      const debugHook = window.__cylinderTilingDebug;
      const realized = debugHook.getRealizedAlphas().flat();
      return {
        state: debugHook.getConstraintState(),
        envelope: debugHook.getEnvelope(),
        realizedMean: realized.reduce((a, b) => a + b, 0) / realized.length,
        collision: debugHook.getCollisionReport(),
      };
    });
  const driveAlpha = async (value) => {
    await page.evaluate((v) => {
      const input = document.getElementById("alpha");
      input.value = v;
      input.dispatchEvent(new Event("input", { bubbles: true }));
    }, value);
    await page.waitForTimeout(350);
  };
  const atDefault = await readConstraint();
  assert.strictEqual(atDefault.state, "free", `${name} default pose should be free (no overlap)`);
  assert.ok(atDefault.collision.clear, `${name} default pose should be collision-free`);
  const window1 = atDefault.envelope.find(([lo, hi]) => 1.3 >= lo && 1.3 <= hi);
  assert.ok(window1, `${name} collision-free range should contain the default effective alpha 1.30 (got ${JSON.stringify(atDefault.envelope)})`);
  // The window's lower edge is where neighbors' 4 mm pads touch at a joint:
  // 2L*sin(theta/2) = 2*4.0 -> theta = 20.5 deg -> effective alpha ~1.15.
  assert.ok(Math.abs(window1[0] - 1.15) <= 0.01, `${name} window should start where neighboring pads touch, ~1.15 (got ${window1[0]})`);

  // Driving the shared alpha past the window limits the drive to just
  // inside the window's edge (0.02 margin from contact), so the rest of the
  // structure isn't jammed against its stops.
  await driveAlpha("2.0");
  const drivenUp = await readConstraint();
  assert.strictEqual(drivenUp.state, "limited", `${name} driving past the window should be limited (got ${drivenUp.state})`);
  assert.ok(drivenUp.collision.clear, `${name} a limited pose must not overlap`);
  assert.ok(
    Math.abs(drivenUp.realizedMean - (window1[1] - 0.02)) < 0.005,
    `${name} should stop just inside the window's upper edge ${window1[1]} (got ${drivenUp.realizedMean})`
  );
  const driveText = await page.evaluate(() => document.getElementById("driveRealizedMetric").textContent);
  assert.ok(driveText.includes("out of reach"), `${name} the Drive panel should say the command is out of reach (got "${driveText}")`);

  await driveAlpha("0.4");
  const drivenDown = await readConstraint();
  assert.ok(drivenDown.collision.clear, `${name} driving down must not overlap either`);
  assert.ok(
    Math.abs(drivenDown.realizedMean - (window1[0] + 0.02)) < 0.005,
    `${name} should stop just inside the window's lower edge ${window1[0]}, not jump to the other clear range (got ${drivenDown.realizedMean})`
  );

  // With constraints off, the pose follows the command straight into overlap.
  await page.evaluate(() => {
    const input = document.getElementById("constraintsEnabled");
    input.checked = false;
    input.dispatchEvent(new Event("change", { bubbles: true }));
  });
  await driveAlpha("1.15");
  const unconstrained = await readConstraint();
  assert.strictEqual(unconstrained.state, "off", `${name} constraint state should read off`);
  assert.ok(Math.abs(unconstrained.realizedMean - 1.05) < 1e-6, `${name} unconstrained pose should follow the command exactly`);
  assert.ok(!unconstrained.collision.clear, `${name} effective alpha 1.05 should be reported as overlapping`);
  await page.evaluate(() => {
    const input = document.getElementById("constraintsEnabled");
    input.checked = true;
    input.dispatchEvent(new Event("change", { bubbles: true }));
  });
  await page.waitForTimeout(150);
  // Re-enabling limits the out-of-reach command back into the window.
  const reenabled = await readConstraint();
  assert.strictEqual(reenabled.state, "limited", `${name} re-enabling with the command out of reach should limit it (got ${reenabled.state})`);
  assert.ok(reenabled.collision.clear, `${name} re-enabling should leave an overlap-free pose`);
  await driveAlpha(defaults.alpha);
  const recovered = await readConstraint();
  assert.strictEqual(recovered.state, "free", `${name} moving back to the default should leave the overlap (got ${recovered.state})`);
  assert.ok(Math.abs(recovered.realizedMean - 1.3) < 1e-6, `${name} should reach the default pose again`);

  // Physical pins (pin-sized backlash, concept from Ahyan's two-cell V2):
  // one per hub, two per circumferential joint (opposite sites pin both
  // pad pairs) and two per axial joint (north/south, likewise) =
  // n*m + 2*n*m + 2*n*(m-1) = 30 + 60 + 40 at the n=10, m=3 default.
  const pinsDefault = await page.evaluate(() => window.__cylinderTilingDebug.getPinInfo());
  assert.strictEqual(pinsDefault.visibleCount, 130, `${name} should render one pin per hub and two per circumferential and axial joint`);
  assert.ok(Math.abs(pinsDefault.renderedRadius - 1.5) < 1e-6, `${name} default 3.0 mm pin should render at 1.5 mm radius`);
  // dphi = asin(radial clearance / arm length) = asin(0.2 / 22.451) = 0.51 deg.
  assert.ok(Math.abs(pinsDefault.deltaPhiDeg - 0.5104) < 0.001, `${name} pin in-plane dead zone should be asin(b/L) (got ${pinsDefault.deltaPhiDeg})`);

  // A thinner pin has more clearance: bigger dead zone, bigger tilt limit,
  // fewer cells needed to close a ring through backlash alone.
  const setPinDiameter = (value) =>
    page.evaluate((v) => {
      const input = document.getElementById("pinDiameter");
      input.value = v;
      input.dispatchEvent(new Event("input", { bubbles: true }));
    }, value);
  await setPinDiameter("2.44");
  await page.waitForTimeout(100);
  const thinPin = await page.evaluate(() => ({
    pins: window.__cylinderTilingDebug.getPinInfo(),
    alignment: window.__cylinderTilingDebug.getPinAlignment(),
  }));
  assert.ok(thinPin.pins.deltaPhiDeg > pinsDefault.deltaPhiDeg, `${name} thinner pin should widen the in-plane dead zone`);
  assert.ok(thinPin.alignment.tiltLimit > pinAlignmentDefault.tiltLimit, `${name} thinner pin should raise the joint tilt limit`);
  assert.ok(
    thinPin.alignment.minCellsForBacklashClosure < pinAlignmentDefault.minCellsForBacklashClosure,
    `${name} thinner pin should need fewer cells per ring for backlash-only closure`
  );

  // A pin that fills the hole leaves no clearance: no tilt allowance, and
  // the minimum-cells readout says so instead of dividing by zero.
  await setPinDiameter("3.4");
  await page.waitForTimeout(100);
  const fullPin = await page.evaluate(() => ({
    alignment: window.__cylinderTilingDebug.getPinAlignment(),
    text: document.getElementById("minCellsClosureMetric").textContent,
  }));
  assert.ok(fullPin.alignment.tiltLimit < 1e-9, `${name} a hole-filling pin should allow no tilt`);
  assert.ok(fullPin.text.includes("none"), `${name} min-cells readout should explain a hole-filling pin (got "${fullPin.text}")`);
  await setPinDiameter("3.0");

  await page.evaluate(() => {
    const input = document.getElementById("showPins");
    input.checked = false;
    input.dispatchEvent(new Event("change", { bubbles: true }));
  });
  await page.waitForTimeout(100);
  const pinsHidden = await page.evaluate(() => window.__cylinderTilingDebug.getPinInfo());
  assert.strictEqual(pinsHidden.visibleCount, 0, `${name} Show pins off should hide every pin`);
  await page.evaluate(() => {
    const input = document.getElementById("showPins");
    input.checked = true;
    input.dispatchEvent(new Event("change", { bubbles: true }));
  });
  await page.waitForTimeout(100);

  // Under differential dilation rows sit at different diameters, but they
  // must stay coaxial on the axle with each column stacked at one angle
  // (axial bolts force that). A regression here once left rows shifted
  // sideways and twisted relative to each other, inflating the axial gap.
  await page.evaluate(() => document.getElementById("presetBarrelBtn").click());
  await page.waitForTimeout(150);
  const barrelGeometry = await page.evaluate(() => {
    const debugHook = window.__cylinderTilingDebug;
    const rows = debugHook.getCellAlphas().length;
    const out = [];
    for (let row = 0; row < rows; row += 1) {
      const ring = debugHook.getRingGeometry(row);
      const n = ring.centers.length;
      const mean = ring.centers.reduce((acc, c) => ({ x: acc.x + c.x / n, y: acc.y + c.y / n }), { x: 0, y: 0 });
      out.push({ meanOffset: Math.hypot(mean.x, mean.y), cell0Angle: Math.atan2(ring.centers[0].y, ring.centers[0].x) });
    }
    return out;
  });
  barrelGeometry.forEach(({ meanOffset, cell0Angle }, row) => {
    assert.ok(meanOffset < 1e-6, `${name} Barrel row ${row} should be centered on the axle (offset ${meanOffset})`);
    assert.ok(Math.abs(cell0Angle) < 1e-6, `${name} Barrel row ${row} cell 0 should sit at polar angle 0 (got ${cell0Angle})`);
  });
  const pinAlignmentBarrel = await page.evaluate(() => window.__cylinderTilingDebug.getPinAlignment());
  assert.ok(pinAlignmentBarrel.maxAxialBend > 0, `${name} Barrel rows at different diameters should need some axial tilt`);
  // Pinned rows at different twists can't have their north/south pads meet
  // exactly: each pad sits L*sin(theta/2) off its hub along the ring, so
  // rows at alpha 1.24 and 1.73 leave L*|sin(t1/2) - sin(t2/2)|
  // for the axial joint to absorb - and nothing more than that (a regression
  // once left unaligned rows 43 mm apart).
  const barrelRows = await page.evaluate(() => window.__cylinderTilingDebug.getRealizedAlphas().map((row) => row[0]));
  const offsetFor = (a) => 22.451 * Math.sin(((70 * a - 60) * Math.PI) / 360);
  let expectedAxialOffset = 0;
  for (let row = 0; row + 1 < barrelRows.length; row += 1) {
    expectedAxialOffset = Math.max(expectedAxialOffset, Math.abs(offsetFor(barrelRows[row]) - offsetFor(barrelRows[row + 1])));
  }
  assert.ok(
    Math.abs(pinAlignmentBarrel.maxAxial - expectedAxialOffset) < 1,
    `${name} Barrel axial pad offset should come only from the rows' different twists (expected ~${expectedAxialOffset.toFixed(2)}, got ${pinAlignmentBarrel.maxAxial})`
  );
  await page.evaluate(() => document.getElementById("clearRolesBtn").click());
  await page.waitForTimeout(150);

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

  // Cell selection: click over the hub of the row-0 cell nearest the camera
  // (the ring is hollow and centered on the axle, so the canvas center looks
  // into empty space). Row 0 because the closure/diameter readouts below
  // report row 0, and rows solve independently.
  // Recomputed per click: Frame Cell below moves the camera, so a point
  // captured once would go stale.
  async function row0CellPoint() {
    const point = await page.evaluate(() => window.__cylinderTilingDebug.getNearestCellClientPoint(0));
    assert.ok(point, `${name} some row-0 cell should be on screen to click`);
    return point;
  }
  async function clickRow0Cell() {
    const point = await row0CellPoint();
    await page.mouse.click(point.x, point.y);
  }
  await clickRow0Cell();
  await page.waitForTimeout(100);
  async function shiftClickSelection() {
    // Locator.click's modifiers option reliably holds Shift for the whole
    // click across viewports; page.keyboard.down/up around page.mouse.click
    // was found to silently drop the modifier on the mobile viewport size.
    const point = await row0CellPoint();
    const box = await page.locator("#threeMount canvas").boundingBox();
    await page.locator("#threeMount canvas").click({
      position: { x: point.x - box.x, y: point.y - box.y },
      modifiers: ["Shift"],
    });
  }
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
  // should (a) show up in the actuator count, (b) actually open that cell
  // and resize its ring (the ring re-solves its radius, so closure stays
  // exact), and (c) pull a neighboring free cell's alpha away from the
  // shared baseline through the backlash-gated propagation - not just the
  // actuated cell itself. With constraints on, the realized pose may stop
  // short of the command at contact, so (b) checks the realized value.
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
  assert.ok(
    parseFloat(afterActuation.selectedCellAlpha) > parseFloat(preActuationAlpha) + 0.05,
    `${name} the actuated cell's realized alpha should rise toward its commanded max (${preActuationAlpha} -> ${afterActuation.selectedCellAlpha})`
  );
  assert.ok(Math.abs(mmValue(afterActuation.closure)) <= 0.01, `${name} the ring should still close exactly with non-uniform per-cell alpha`);
  assert.ok(
    Math.abs(mmValue(afterActuation.diameter) - mmValue(beforeActuation.diameter)) > 0.1,
    `${name} actuating one cell should resize its ring (${beforeActuation.diameter} -> ${afterActuation.diameter})`
  );
  // Both pad pairs of every joint are bolted, so the actuated cell drags
  // the whole ring toward 2.0 - past the reachable range, so every target
  // is limited to just inside it, with nothing overlapping or binding.
  await page.waitForFunction(() => window.__cylinderTilingDebug.getConstraintState() !== "moving", null, { timeout: 5000 });
  const actuatedContact = await page.evaluate(() => ({
    state: window.__cylinderTilingDebug.getConstraintState(),
    clear: window.__cylinderTilingDebug.getCollisionReport().clear,
    bind: window.__cylinderTilingDebug.getPinAlignment().maxCircumferential,
  }));
  assert.strictEqual(actuatedContact.state, "limited", `${name} an actuator commanded past reach should be limited (got ${actuatedContact.state})`);
  assert.ok(actuatedContact.clear, `${name} a limited actuator must not overlap its neighbors`);
  assert.ok(actuatedContact.bind <= 0.4 + 1e-3, `${name} bolted pins must stay within their 0.4 mm clearance (got ${actuatedContact.bind})`);

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
    Math.abs(propagation.neighborAlpha - 1.3) > 0.01,
    `${name} the actuated cell's free neighbor should shift away from the 1.30 no-actuator baseline (got ${propagation.neighborAlpha})`
  );

  // Alpha heatmap: read the actuated cell's actual material color through
  // the debug hook (precise; a screenshot pixel-color approach turned out
  // to be too fragile here, since lighting/shading shifts rendered pixels
  // away from the material's raw hex color, and several row-palette colors
  // incidentally fall inside any reasonably wide "reddish" band). The
  // heatmap colors the *realized* alpha; with constraints on, the actuator
  // (commanded 2.0) is held at contact somewhere above the 1.2 midpoint, so
  // its color should be on the red "expanded" side once the toggle is on,
  // and back to its row palette color once off.
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
  const realizedActuatorAlpha = await page.evaluate(() => {
    const selectedCell = window.__cylinderTilingDebug.getSelected();
    return window.__cylinderTilingDebug.getRealizedAlphas()[selectedCell.row][selectedCell.i];
  });
  assert.ok(realizedActuatorAlpha > 1.2, `${name} the actuator should be realized above the heatmap midpoint (got ${realizedActuatorAlpha})`);
  const [heatR, heatG, heatB] = [1, 3, 5].map((k) => parseInt(colorDuringHeatmap.slice(k, k + 2), 16));
  assert.notStrictEqual(colorDuringHeatmap, colorBeforeHeatmap, `${name} enabling the heatmap should recolor the cell`);
  assert.ok(heatR > heatG && heatR > heatB, `${name} an expanded cell should render on the heatmap's red side (got ${colorDuringHeatmap})`);
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

  // Batch selection (shift-click): toggle a cell into/out of a multi-select
  // set without disturbing the single-cell Selected Cell inspector, then
  // apply a role/alpha to every selected cell at once - the "select
  // individual cells to expand/contract independently" workflow, distinct
  // from actuating one cell at a time.
  await clickRow0Cell();
  await page.waitForTimeout(100);
  const batchTargetCell = await page.evaluate(() => window.__cylinderTilingDebug.getSelected());
  assert.ok(batchTargetCell, `${name} a cell should be selectable at the known screen point for batch testing`);

  await shiftClickSelection();
  await page.waitForTimeout(100);
  const batchCountAfterAdd = await page.evaluate(() => document.getElementById("batchSelectedCountMetric").textContent);
  assert.strictEqual(batchCountAfterAdd, "1", `${name} shift-clicking a cell should add it to the batch selection`);

  await shiftClickSelection();
  await page.waitForTimeout(100);
  const batchCountAfterToggleOff = await page.evaluate(() => document.getElementById("batchSelectedCountMetric").textContent);
  assert.strictEqual(batchCountAfterToggleOff, "0", `${name} shift-clicking the same cell again should remove it from the batch selection`);

  await shiftClickSelection();
  await page.waitForTimeout(100);
  await page.evaluate(() => {
    document.getElementById("batchRoleSelect").value = "actuator";
    const alphaInput = document.getElementById("batchAlpha");
    alphaInput.value = "1.8";
    alphaInput.dispatchEvent(new Event("input", { bubbles: true }));
    document.getElementById("applyBatchBtn").click();
  });
  await page.waitForTimeout(150);
  const batchApplyResult = await page.evaluate((cell) => {
    const roles = window.__cylinderTilingDebug.getCellRoles();
    const alphas = window.__cylinderTilingDebug.getCellAlphas();
    return { role: roles[cell.row][cell.i].role, alpha: alphas[cell.row][cell.i] };
  }, batchTargetCell);
  assert.strictEqual(batchApplyResult.role, "actuator", `${name} Apply to Selected should set the batch-selected cell's role`);
  assert.ok(
    Math.abs(batchApplyResult.alpha - 1.8) < 0.05,
    `${name} Apply to Selected should command the batch-selected cell's alpha (got ${batchApplyResult.alpha})`
  );

  // Clear Selection should empty the batch set but must not undo the role/
  // alpha it already applied.
  await page.click("#clearBatchSelectionBtn");
  await page.waitForTimeout(100);
  const batchCountAfterClear = await page.evaluate(() => document.getElementById("batchSelectedCountMetric").textContent);
  assert.strictEqual(batchCountAfterClear, "0", `${name} Clear Selection should empty the batch selection`);
  const roleAfterClearSelection = await page.evaluate(
    (cell) => window.__cylinderTilingDebug.getCellRoles()[cell.row][cell.i].role,
    batchTargetCell
  );
  assert.strictEqual(roleAfterClearSelection, "actuator", `${name} clearing the batch selection should not undo the role it already applied`);

  await page.click("#clearRolesBtn");
  await page.waitForTimeout(150);

  // Shape presets: each should command a different alpha per row and, with
  // constraints on, actually reach it - producing the named profile in the
  // realized ring diameters. Rows are bolted to each other (both pad pairs),
  // so adjacent rows can only differ by what the pins' clearance allows:
  // with the snug default 3.0 mm pin the profile is small (~1 mm), and a
  // thinner pin allows a bigger one.
  const presetShape = async (buttonId) => {
    await page.evaluate((id) => {
      document.getElementById("clearRolesBtn").click();
      document.getElementById(id).click();
    }, buttonId);
    await page.waitForFunction(() => window.__cylinderTilingDebug.getConstraintState() !== "moving", null, { timeout: 5000 });
    return page.evaluate(() => {
      const debugHook = window.__cylinderTilingDebug;
      const rows = debugHook.getCellAlphas().length;
      const diameters = [];
      for (let row = 0; row < rows; row += 1) {
        const ring = debugHook.getRingGeometry(row);
        diameters.push((2 * ring.centers.reduce((sum, c) => sum + Math.hypot(c.x, c.y), 0)) / ring.centers.length);
      }
      return {
        diameters,
        state: debugHook.getConstraintState(),
        clear: debugHook.getCollisionReport().clear,
        axialPinOffset: debugHook.getPinAlignment().maxAxial,
      };
    });
  };
  const barrel = await presetShape("presetBarrelBtn");
  assert.ok(barrel.diameters.length >= 3, `${name} preset test needs at least 3 rows (the default row count)`);
  const midRow = Math.floor((barrel.diameters.length - 1) / 2);
  const lastRow = barrel.diameters.length - 1;
  assert.strictEqual(barrel.state, "free", `${name} Barrel preset should be reachable under constraints (got ${barrel.state})`);
  assert.ok(barrel.clear, `${name} Barrel pose should be collision-free`);
  assert.ok(barrel.axialPinOffset <= 0.4 + 1e-3, `${name} Barrel rows must stay within the 0.4 mm pin clearance (got ${barrel.axialPinOffset})`);
  assert.ok(
    barrel.diameters[midRow] > barrel.diameters[0] + 0.3 && barrel.diameters[midRow] > barrel.diameters[lastRow] + 0.3,
    `${name} Barrel should bulge the middle row (diameters ${JSON.stringify(barrel.diameters)})`
  );
  const saddle = await presetShape("presetSaddleBtn");
  assert.strictEqual(saddle.state, "free", `${name} Saddle preset should be reachable under constraints (got ${saddle.state})`);
  assert.ok(
    saddle.diameters[midRow] < saddle.diameters[0] - 0.3 && saddle.diameters[midRow] < saddle.diameters[lastRow] - 0.3,
    `${name} Saddle should pinch the middle row (diameters ${JSON.stringify(saddle.diameters)})`
  );
  const cone = await presetShape("presetConeBtn");
  assert.strictEqual(cone.state, "free", `${name} Cone preset should be reachable under constraints (got ${cone.state})`);
  assert.ok(
    cone.diameters.every((d, row) => row === 0 || d > cone.diameters[row - 1]),
    `${name} Cone should widen row by row (diameters ${JSON.stringify(cone.diameters)})`
  );
  // A looser (thinner) pin lets bolted rows differ more: bigger barrel.
  await page.evaluate(() => {
    const input = document.getElementById("pinDiameter");
    input.value = "2.3";
    input.dispatchEvent(new Event("input", { bubbles: true }));
  });
  const looseBarrel = await presetShape("presetBarrelBtn");
  const bulge = (shape) => shape.diameters[midRow] - (shape.diameters[0] + shape.diameters[lastRow]) / 2;
  assert.ok(
    bulge(looseBarrel) > 2 * bulge(barrel),
    `${name} a 2.3 mm pin should allow a much bigger barrel than 3.0 mm (bulge ${bulge(looseBarrel)} vs ${bulge(barrel)})`
  );
  await page.evaluate(() => {
    const input = document.getElementById("pinDiameter");
    input.value = input.defaultValue;
    input.dispatchEvent(new Event("input", { bubbles: true }));
  });

  await page.click("#clearRolesBtn");
  await page.waitForTimeout(300);

  // Re-actuate once more so the JSON save/load round-trip below has
  // non-trivial per-cell state to carry through.
  await clickRow0Cell();
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
    document.getElementById("pinDiameter").value = "2.6";
    document.getElementById("pinDiameter").dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.click("#saveJsonBtn");
  const download = await downloadPromise;
  const savedPath = path.join(root, `cylinder-tiling-${name}-saved.json`);
  await download.saveAs(savedPath);
  const saved = require(savedPath);
  assert.strictEqual(saved.format, "rad-cylinder-tiling.v1", `${name} saved JSON should carry the expected format tag`);
  assert.strictEqual(saved.alpha, 1.7, `${name} saved JSON should capture the current alpha`);
  assert.strictEqual(saved.pinDiameter, 2.6, `${name} saved JSON should capture the pin diameter`);
  assert.ok(Array.isArray(saved.cellRoles), `${name} saved JSON should include per-cell roles`);
  const savedLockedCount = saved.cellRoles.flat().filter((c) => c && c.role === "locked").length;
  assert.strictEqual(savedLockedCount, 1, `${name} saved JSON should capture the locked cell`);

  await page.evaluate(() => {
    document.getElementById("alpha").value = "0.4";
    document.getElementById("alpha").dispatchEvent(new Event("input", { bubbles: true }));
    document.getElementById("pinDiameter").value = "3.2";
    document.getElementById("pinDiameter").dispatchEvent(new Event("input", { bubbles: true }));
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
  const loadedPin = await page.evaluate(() => document.getElementById("pinDiameter").value);
  assert.strictEqual(Number(loadedPin), 2.6, `${name} loading the saved JSON should restore the pin diameter`);

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
  await page.waitForFunction(() => window.__cylinderTilingDebug.getConstraintState() !== "moving", null, { timeout: 5000 });
  const afterDragOnly = await readMetrics(page);
  const dragAchieved = await page.evaluate(() => document.getElementById("achievedDiameterMetric").textContent);
  assert.notStrictEqual(afterDragOnly.alpha, alphaBeforeDrag, `${name} dragging Target Diameter alone (no button click) should already change alpha`);
  assert.ok(Math.abs(mmValue(dragAchieved) - 125) < 1, `${name} the fit should land within 1 mm of a reachable 125 mm target (got ${dragAchieved})`);
  assert.strictEqual(
    mmValue(afterDragOnly.diameter).toFixed(1),
    mmValue(dragAchieved).toFixed(1),
    `${name} dragging Target Diameter alone should already re-render the ring at the fitted diameter`
  );

  await page.evaluate(() => {
    document.getElementById("targetDiameter").value = "140";
    document.getElementById("targetDiameter").dispatchEvent(new Event("input", { bubbles: true }));
  });
  await page.click("#fitDiameterBtn");
  await page.waitForFunction(() => window.__cylinderTilingDebug.getConstraintState() !== "moving", null, { timeout: 5000 });
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
  await clickRow0Cell();
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
  assert.strictEqual(afterResetAll.alpha, "1.4", `${name} Reset All should restore the default alpha`);
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
