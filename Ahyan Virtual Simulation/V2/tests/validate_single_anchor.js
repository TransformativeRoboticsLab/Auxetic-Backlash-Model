const assert = require('assert');
const path = require('path');
const fs = require('fs');
const { pathToFileURL } = require('url');
const { chromium } = require('playwright');

const near = (a, b, tolerance, label) => assert.ok(Math.abs(a - b) <= tolerance, `${label}: ${a} vs ${b}`);
const grids = [[1,1],[1,2],[2,1],[1,10],[10,1],[1,20],[20,1],[1,30],[30,1],
  [2,2],[3,3],[3,7],[7,3],[7,7],[15,15],[30,30]];
const ratios = [.35,.5,.72,.9,.99,1];

function verify(free, anchored, options) {
  assert.ok(!anchored.driveDisabled && !anchored.animationDisabled, 'one anchor cannot jam the mechanism');
  near(anchored.uniformEffectiveDriveDeltaDeg, free.uniformEffectiveDriveDeltaDeg, 1e-8, 'uniform travel');
  const lockIndex = options.locks[0][0] * options.cols + options.locks[0][1];
  anchored.samples.forEach((sample, index) => {
    const reference = free.samples[index];
    for (const key of ['min','max','value']) near(sample.slider[key], reference.slider[key], 1e-8, `slider ${key}`);
    near(sample.effectiveDrive, reference.effectiveDrive, 1e-8, 'selected angle');
    near(sample.effectiveDrive, sample.slider.value, 1e-8, 'actual vs requested');
    const anchor = sample.cells[lockIndex], lock = anchor.lock;
    near(anchor.center[0], lock.targetX, 1e-9, 'anchor X');
    near(anchor.center[1], lock.targetY, 1e-9, 'anchor Y');
    if (lock.lockZ) near(anchor.center[2], lock.targetZ, 1e-9, 'anchor Z');
    const shift = anchor.center.map((value, k) => value - reference.cells[lockIndex].center[k]);
    sample.cells.forEach((cell, i) => {
      near(cell.delta, reference.cells[i].delta, 1e-8, 'passive angle unchanged');
      cell.center.forEach((value, k) => near(value - reference.cells[i].center[k], shift[k], 1e-6, 'rigid translation only'));
      cell.normal.forEach((value, k) => near(value, reference.cells[i].normal[k], 1e-8, 'orientation unchanged'));
    });
    assert.ok(sample.pinFit.clear, 'pin containment');
    assert.strictEqual(sample.pinFit.visiblePins, 2 * (options.rows * (options.cols - 1) + options.cols * (options.rows - 1)));
    sample.pins.forEach((pin, i) => {
      pin.origin.forEach((value, k) => near(value - reference.pins[i].origin[k], shift[k], 1e-5, 'pins translated with cells'));
      pin.holes.forEach(hole => {
        assert.ok(Math.hypot(pin.origin[0] - hole.center[0], pin.origin[1] - hole.center[1]) + pin.radius <= hole.radius + 1e-4,
          'independent pin/bore radial check');
        assert.ok(Math.abs(pin.origin[2] - hole.center[2]) + hole.thickness / 2 <= pin.halfLength + 1e-4,
          'pin spans plate thickness');
      });
    });
  });
}

(async () => {
  const browser = await chromium.launch({ channel: 'chrome' });
  let cases = 0;
  try {
    const page = await browser.newPage({ viewport: { width: 1280, height: 900 } });
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    await page.goto(pathToFileURL(path.resolve(__dirname, '../web/two-cell-attachment/index.html')).href);
    await page.waitForFunction(() => window.RADAttachmentDiagnostics);
    const run = options => page.evaluate(o => RADAttachmentDiagnostics.lockedStringBacklashScenario(o), options);
    for (const [g, [rows, cols]] of (process.env.RAD_ANCHOR_UI_ONLY ? [] : grids).entries()) {
      if (process.env.RAD_ANCHOR_SMOKE && !(rows === 7 && cols === 7)) continue;
      if (process.env.RAD_ANCHOR_LARGE_ONLY && rows * cols !== 900) continue;
      const middle = [Math.floor(rows / 2), Math.floor(cols / 2)];
      const end = [rows - 1, cols - 1];
      const arrangements = rows * cols === 900 ? [[middle, middle]]
        : [[middle, middle], [[0,0], end], [[0, Math.floor(cols / 2)], middle]];
      for (const [i, [anchor, source]] of arrangements.entries()) {
        const options = { rows, cols, sourceRow: source[0], sourceCol: source[1], ratio: rows * cols === 900 ? .72 : ratios[(g + i) % ratios.length],
          driveFractions: rows * cols === 900 ? [0,1] : [0,.25,.75,1,.5,0] };
        const free = await run({ ...options, locks: [] });
        if (rows * cols === 900) console.log('30x30 free reference complete');
        for (const lockMode of ['xy','xyz']) {
          const locked = { ...options, locks: [anchor], lockMode };
          verify(free, await run(locked), locked);
          cases++;
          if (rows * cols === 900) console.log(`30x30 ${lockMode} anchor passed`);
        }
      }
      console.log(`${rows}x${cols}: ${arrangements.length} anchor arrangements passed for XY and XYZ`);
    }
    if (process.env.RAD_ANCHOR_LARGE_ONLY) {
      assert.deepStrictEqual(errors, []);
      console.log(`${cases} large-lattice single-anchor scenarios passed`);
      return;
    }
    for (const lockMode of ['xy', 'xyz']) {
      const options = { rows: 3, cols: 4, sourceRow: 1, sourceCol: 2, ratio: .72,
        driveFractions: [0,.5,1,.5,0] };
      const locked = { ...options, locks: [[0,0]], lockMode, lockOffset: [73,-41,12] };
      const result = await run(locked);
      verify(await run({ ...options, locks: [] }), result, locked);
      result.samples.forEach(sample => sample.cells.forEach(cell =>
        near(cell.center[2], lockMode === 'xyz' ? 12 : 0, 1e-9, 'only XYZ fixes height')));
      cases++;
    }

    const input = async (id, value) => {
      await page.locator('#' + id).evaluate((element, v) => {
        element.value = String(v);
        element.dispatchEvent(new Event('input', { bubbles: true }));
      }, value);
      await page.waitForTimeout(100);
    };
    const state = () => page.evaluate(() => RADAttachmentDiagnostics.currentState());
    const unchanged = (before, after, label) => after.cells.forEach((cell, i) =>
      cell.center.forEach((v,k) => near(v, before.cells[i].center[k], 1e-5, label)));
    const fixed = snapshot => {
      snapshot.cells.filter(cell => cell.lock).forEach(cell => {
        near(cell.center[0], cell.lock.targetX, 1e-9, 'live anchor X');
        near(cell.center[1], cell.lock.targetY, 1e-9, 'live anchor Y');
        if (cell.lock.lockZ) near(cell.center[2], cell.lock.targetZ, 1e-9, 'live anchor Z');
      });
      assert.ok(snapshot.pinFit.clear, 'live pins remain in holes');
    };
    const select = async (row, col) => {
      await page.locator('#fitView').click();
      const before = await state();
      await page.mouse.click(...before.cells.find(c => c.row === row && c.col === col).screen);
      await page.waitForTimeout(100);
      const after = await state();
      assert.deepStrictEqual(after.selected, { row, col });
      unchanged(before, after, 'selection must not move cells');
    };
    for (const [rows, cols] of [[1,20],[20,1],[3,4]]) {
      for (const lockMode of ['xy', 'xyz']) {
        await page.locator('#clearXYLocks').click();
        await page.locator('#actuationMode').selectOption('uniform');
        await input('rowCount', rows); await input('colCount', cols);
        await page.locator('#applyGrid').click();
        await input('aTopAngle', 48);
        await page.locator('[data-view="top"]').click();
        await select(0,0);
        await page.locator('#positionLockMode').selectOption(lockMode);
        let before = await state();
        await page.locator('#lockSelectedXY').click();
        unchanged(before, await state(), 'locking must not move cells');
        for (const command of [64,40,48]) {
          await input('aTopAngle', command);
          fixed(await state());
          near(Number(await page.locator('#aTopAngle').inputValue()), command, 1e-8, 'uniform command');
        }
        before = await state();
        await page.locator('#actuationMode').selectOption('contact');
        unchanged(before, await state(), 'uniform/contact transition');
        await input('pinRadiusRatio', .72);
        for (const [r,c] of [[0,0],[rows-1,cols-1],[Math.floor(rows/2),Math.floor(cols/2)]]) {
          await select(r,c);
          const selected = await state();
          near(Number(await page.locator('#aTopAngle').inputValue()), selected.field[r][c], 1e-8, 'selected slider');
          for (const command of [64,40,55]) {
            await input('aTopAngle', command);
            const current = await state();
            fixed(current);
            near(current.field[r][c], command, 1e-8, 'sequential actual actuation');
            assert.ok(!await page.locator('#aTopAngle').isDisabled());
          }
        }
        before = await state();
        await page.locator('#clearXYLocks').click();
        unchanged(before, await state(), 'unlocking must not recenter');
        await page.locator('#lockSelectedXY').click();
        unchanged(before, await state(), 'reanchor after deformation');
        await page.locator('#animate').check();
        for (let n = 0; n < 3; n++) { await page.waitForTimeout(200); fixed(await state()); }
        await page.locator('#animate').uncheck();
        console.log(`Live ${rows}x${cols} ${lockMode}: uniform, sequential, unlock/reanchor and animation passed`);
      }
    }
    await page.locator('#clearXYLocks').click();
    await page.locator('#actuationMode').selectOption('uniform');
    await input('rowCount', 1); await input('colCount', 10);
    await page.locator('#applyGrid').click();
    await input('aTopAngle', 26);
    await page.locator('#positionLockMode').selectOption('xyz');
    await select(0,0); await page.locator('#lockSelectedXY').click();
    await select(0,9); await page.locator('#lockSelectedXY').click();
    assert.ok(await page.locator('#aTopAngle').isDisabled(), 'two zero-clearance anchors must still jam');
    fixed(await state());
    await page.locator('#unlockSelectedXY').click();
    assert.ok(!await page.locator('#aTopAngle').isDisabled(), 'removing second anchor restores travel');
    await input('aTopAngle', 64); fixed(await state());
    await page.locator('#actuationMode').selectOption('contact');
    await input('pinRadiusRatio', .72);
    await input('aTopAngle', 45); fixed(await state());
    await page.locator('#showPins').uncheck();
    await input('aTopAngle', 55); fixed(await state());
    await page.locator('#showPins').check();
    await input('rowCount', 3); await input('colCount', 4);
    await page.locator('#applyGrid').click();
    for (const command of [40,64,45]) { await input('aTopAngle', command); fixed(await state()); }
    assert.strictEqual((await state()).cells.filter(cell => cell.lock).length, 1, 'resize preserves valid anchor');
    const output = path.resolve(__dirname, '../research/constraint_audit/after-fixes');
    fs.mkdirSync(output, { recursive: true });
    await page.locator('[data-view="iso"]').click();
    await page.locator('#fitView').click();
    const canvasPixels = async () => {
      const png = await page.locator('canvas').screenshot();
      const colored = await page.evaluate(async data => {
        const image = new Image(); image.src = data; await image.decode();
        const canvas = document.createElement('canvas'); canvas.width = image.width; canvas.height = image.height;
        const context = canvas.getContext('2d'); context.drawImage(image, 0, 0);
        const pixels = context.getImageData(0,0,canvas.width,canvas.height).data;
        let count = 0;
        for (let i = 0; i < pixels.length; i += 4) if (pixels[i+2] - pixels[i] > 40 && pixels[i] < 180) count++;
        return count;
      }, 'data:image/png;base64,' + png.toString('base64'));
      assert.ok(colored > 100, `nonblank colored lattice pixels: ${colored}`);
    };
    await canvasPixels();
    await page.screenshot({ path: path.join(output, 'single-anchor-desktop.png') });
    await page.setViewportSize({ width: 480, height: 900 });
    await page.locator('#fitView').click();
    await canvasPixels();
    await page.screenshot({ path: path.join(output, 'single-anchor-mobile.png') });
    assert.deepStrictEqual(errors, []);
    console.log(`${cases} single-anchor scenarios passed (full range, free-motion equivalence, fixed axes, pin containment)`);
  } finally { await browser.close(); }
})().catch(error => { console.error(error); process.exitCode = 1; });
