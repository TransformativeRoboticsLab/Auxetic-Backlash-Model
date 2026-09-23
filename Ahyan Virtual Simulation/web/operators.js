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

  function boolMatrixDiffCount(a, b) {
    let count = 0;
    for (let r = 0; r < a.length; r += 1) {
      for (let c = 0; c < a[r].length; c += 1) {
        if (Boolean(a[r][c]) !== Boolean(b?.[r]?.[c])) count += 1;
      }
    }
    return count;
  }

  function stateDistance(a, b, simA = null, simB = null) {
    const firstSim = simA || RAD.simulate(a);
    const secondSim = simB || RAD.simulate(b);
    const commandAlphaError = maxMatrixDiff(a.cells.commandAlpha, b.cells.commandAlpha);
    const commandZError = maxMatrixDiff(a.cells.commandZ, b.cells.commandZ);
    const lockAlphaError = maxMatrixDiff(a.cells.lockAlpha, b.cells.lockAlpha);
    const finalAlphaError = maxMatrixDiff(firstSim.alpha, secondSim.alpha);
    const finalHeightError = maxMatrixDiff(firstSim.height, secondSim.height);
    return {
      modeChanges: boolMatrixDiffCount(a.cells.locked, b.cells.locked),
      commandAlphaError,
      commandZError,
      commandError: Math.max(commandAlphaError, commandZError),
      lockAlphaError,
      finalAlphaError,
      finalHeightError,
      finalError: Math.max(finalAlphaError, finalHeightError),
    };
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

  function compareSequenceOrder(state, events, tolerance = 1e-9) {
    const sequence = (events || []).filter(Boolean).map(cloneData);
    const baseFinal = applyEventSequence(state, sequence);
    const baseSim = RAD.simulate(baseFinal);
    const reverseFinal = applyEventSequence(state, sequence.slice().reverse());
    const reverse = stateDistance(baseFinal, reverseFinal, baseSim);
    let maxAdjacentAlphaError = 0;
    let maxAdjacentHeightError = 0;
    let maxAdjacentCommandError = 0;
    let maxAdjacentLockAlphaError = 0;
    let noncommutingAdjacentPairs = 0;
    const adjacent = [];
    for (let index = 0; index + 1 < sequence.length; index += 1) {
      const swapped = sequence.map(cloneData);
      [swapped[index], swapped[index + 1]] = [swapped[index + 1], swapped[index]];
      const swappedFinal = applyEventSequence(state, swapped);
      const distance = stateDistance(baseFinal, swappedFinal, baseSim);
      const sensitive =
        distance.modeChanges > 0 ||
        distance.commandError > tolerance ||
        distance.lockAlphaError > tolerance ||
        distance.finalAlphaError > tolerance ||
        distance.finalHeightError > tolerance;
      if (sensitive) noncommutingAdjacentPairs += 1;
      maxAdjacentAlphaError = Math.max(maxAdjacentAlphaError, distance.finalAlphaError);
      maxAdjacentHeightError = Math.max(maxAdjacentHeightError, distance.finalHeightError);
      maxAdjacentCommandError = Math.max(maxAdjacentCommandError, distance.commandError);
      maxAdjacentLockAlphaError = Math.max(maxAdjacentLockAlphaError, distance.lockAlphaError);
      adjacent.push({
        index,
        firstKind: sequence[index]?.kind || "",
        secondKind: sequence[index + 1]?.kind || "",
        sensitive,
        ...distance,
      });
    }
    const maxOrderError = Math.max(
      reverse.finalError,
      reverse.commandError,
      reverse.lockAlphaError,
      maxAdjacentAlphaError,
      maxAdjacentHeightError,
      maxAdjacentCommandError,
      maxAdjacentLockAlphaError
    );
    return {
      eventCount: sequence.length,
      adjacentPairCount: Math.max(0, sequence.length - 1),
      noncommutingAdjacentPairs,
      orderSensitive:
        reverse.modeChanges > 0 ||
        maxOrderError > tolerance ||
        noncommutingAdjacentPairs > 0,
      reverseAlphaError: reverse.finalAlphaError,
      reverseHeightError: reverse.finalHeightError,
      reverseModeChanges: reverse.modeChanges,
      reverseCommandError: reverse.commandError,
      reverseLockAlphaError: reverse.lockAlphaError,
      maxAdjacentAlphaError,
      maxAdjacentHeightError,
      maxAdjacentCommandError,
      maxAdjacentLockAlphaError,
      maxOrderError,
      adjacent,
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
  RAD.compareSequenceOrder = compareSequenceOrder;
})();
