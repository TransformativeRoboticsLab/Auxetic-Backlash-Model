(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  function cloneData(value) {
    return JSON.parse(JSON.stringify(value));
  }

  function ensureLockAlpha(state) {
    const { rows, cols, initialAlpha } = state.grid;
    if (!state.cells.lockAlpha) {
      state.cells.lockAlpha = RAD.matrix(rows, cols, Number(initialAlpha) || 1);
    }
    return state;
  }

  function readCell(state, cell) {
    const r = Array.isArray(cell) ? Number(cell[0]) : Number(cell?.r);
    const c = Array.isArray(cell) ? Number(cell[1]) : Number(cell?.c);
    if (!Number.isInteger(r) || !Number.isInteger(c)) throw new Error("event requires an integer cell");
    if (r < 0 || c < 0 || r >= state.grid.rows || c >= state.grid.cols) {
      throw new Error(`cell ${r},${c} is outside lattice ${state.grid.rows}x${state.grid.cols}`);
    }
    return { r, c };
  }

  function finiteDieOffRadius(dieOff) {
    let radius = 0;
    for (const row of dieOff || []) {
      for (const value of row || []) {
        if (Number.isFinite(value)) radius = Math.max(radius, value);
      }
    }
    return radius;
  }

  function localActuationEvent(cell, alpha = 0, z = 0) {
    return { kind: "actuate", cell, alpha: Number(alpha) || 0, z: Number(z) || 0 };
  }

  function lockEvent(cell) {
    return { kind: "lock", cell };
  }

  function releaseEvent(cell) {
    return { kind: "release", cell };
  }

  function clearActuationEvent(cell = null) {
    return { kind: "clear_actuation", cell };
  }

  function applyProgrammableEvent(state, event, options = {}) {
    const next = options.mutate ? state : cloneData(state);
    ensureLockAlpha(next);
    if (event.kind === "actuate") {
      const { r, c } = readCell(next, event.cell);
      const alpha = (Number(next.cells.commandAlpha[r][c]) || 0) + (Number(event.alpha) || 0);
      const z = (Number(next.cells.commandZ[r][c]) || 0) + (Number(event.z) || 0);
      next.cells.commandAlpha[r][c] = RAD.clampCommandAlpha ? RAD.clampCommandAlpha(next, alpha) : alpha;
      next.cells.commandZ[r][c] = RAD.clampCommandZ ? RAD.clampCommandZ(next, z) : z;
    } else if (event.kind === "lock") {
      const { r, c } = readCell(next, event.cell);
      const sim = RAD.simulate(next);
      next.cells.lockAlpha[r][c] = sim.alpha[r][c];
      next.cells.commandAlpha[r][c] = 0;
      next.cells.commandZ[r][c] = 0;
      next.cells.locked[r][c] = true;
    } else if (event.kind === "release") {
      const { r, c } = readCell(next, event.cell);
      next.cells.locked[r][c] = false;
    } else if (event.kind === "clear_actuation") {
      if (event.cell === null || event.cell === undefined) {
        for (let r = 0; r < next.grid.rows; r += 1) {
          for (let c = 0; c < next.grid.cols; c += 1) {
            next.cells.commandAlpha[r][c] = 0;
            next.cells.commandZ[r][c] = 0;
          }
        }
      } else {
        const { r, c } = readCell(next, event.cell);
        next.cells.commandAlpha[r][c] = 0;
        next.cells.commandZ[r][c] = 0;
      }
    } else {
      throw new Error(`unsupported programmable discontinuity event: ${event.kind}`);
    }
    return next;
  }

  function applyEventSequence(state, events, options = {}) {
    let current = options.mutate ? state : cloneData(state);
    for (const event of events || []) {
      current = applyProgrammableEvent(current, event, { mutate: true });
    }
    return current;
  }

  function maxMatrixDiff(a, b) {
    let out = 0;
    for (let r = 0; r < a.length; r += 1) {
      for (let c = 0; c < a[r].length; c += 1) {
        out = Math.max(out, Math.abs((Number(a[r][c]) || 0) - (Number(b?.[r]?.[c]) || 0)));
      }
    }
    return out;
  }

  function matrixEqual(a, b, tolerance = 1e-9) {
    return maxMatrixDiff(a, b) <= tolerance;
  }

  function boolMatrixEqual(a, b) {
    return JSON.stringify(a) === JSON.stringify(b);
  }

  function compareEventOrder(state, first, second, tolerance = 1e-9) {
    const firstThenSecond = applyEventSequence(state, [first, second]);
    const secondThenFirst = applyEventSequence(state, [second, first]);
    const firstSim = RAD.simulate(firstThenSecond);
    const secondSim = RAD.simulate(secondThenFirst);
    return {
      firstThenSecond,
      secondThenFirst,
      modeCommutes: boolMatrixEqual(firstThenSecond.cells.locked, secondThenFirst.cells.locked),
      commandCommutes:
        matrixEqual(firstThenSecond.cells.commandAlpha, secondThenFirst.cells.commandAlpha, tolerance) &&
        matrixEqual(firstThenSecond.cells.commandZ, secondThenFirst.cells.commandZ, tolerance),
      lockAlphaCommutes: matrixEqual(firstThenSecond.cells.lockAlpha, secondThenFirst.cells.lockAlpha, tolerance),
      finalAlphaError: maxMatrixDiff(firstSim.alpha, secondSim.alpha),
      finalHeightError: maxMatrixDiff(firstSim.height, secondSim.height),
    };
  }

  RAD.finiteDieOffRadius = finiteDieOffRadius;
  RAD.localActuationEvent = localActuationEvent;
  RAD.lockEvent = lockEvent;
  RAD.releaseEvent = releaseEvent;
  RAD.clearActuationEvent = clearActuationEvent;
  RAD.applyProgrammableEvent = applyProgrammableEvent;
  RAD.applyEventSequence = applyEventSequence;
  RAD.compareEventOrder = compareEventOrder;
})();
