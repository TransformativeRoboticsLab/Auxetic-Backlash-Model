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
    if (!state.cells.lockZ) state.cells.lockZ = RAD.matrix(rows, cols, 0);
    return state;
  }

  function ensureRemoved(state) {
    const { rows, cols } = state.grid;
    if (!state.cells.removed) state.cells.removed = RAD.matrix(rows, cols, false);
    return state;
  }

  function ensurePositionLocked(state) {
    const { rows, cols } = state.grid;
    if (!state.cells.positionLocked) state.cells.positionLocked = RAD.matrix(rows, cols, false);
    if (!state.cells.positionLockX) state.cells.positionLockX = RAD.matrix(rows, cols, null);
    if (!state.cells.positionLockY) state.cells.positionLockY = RAD.matrix(rows, cols, null);
    if (!state.cells.positionLockZ) state.cells.positionLockZ = RAD.matrix(rows, cols, null);
    return state;
  }

  function clearPositionLockTarget(state, r, c) {
    ensurePositionLocked(state);
    state.cells.positionLockX[r][c] = null;
    state.cells.positionLockY[r][c] = null;
    state.cells.positionLockZ[r][c] = null;
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

  function groupActuationEvent(cells, alpha = 0, z = 0) {
    return {
      kind: "group_actuate",
      cells: (cells || []).map((cell) => readCell({ grid: { rows: Infinity, cols: Infinity } }, cell)),
      alpha: Number(alpha) || 0,
      z: Number(z) || 0,
    };
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

  function removeCellEvent(cell) {
    return { kind: "remove_cell", cell };
  }

  function restoreCellEvent(cell) {
    return { kind: "restore_cell", cell };
  }

  function applyProgrammableEvent(state, event, options = {}) {
    const next = options.mutate ? state : cloneData(state);
    ensureLockAlpha(next);
    ensureRemoved(next);
    ensurePositionLocked(next);
    if (event.kind === "actuate") {
      const { r, c } = readCell(next, event.cell);
      if (next.cells.removed?.[r]?.[c]) return next;
      const alpha = (Number(next.cells.commandAlpha[r][c]) || 0) + (Number(event.alpha) || 0);
      const z = (Number(next.cells.commandZ[r][c]) || 0) + (Number(event.z) || 0);
      next.cells.commandAlpha[r][c] = RAD.clampCommandAlpha ? RAD.clampCommandAlpha(next, alpha) : alpha;
      next.cells.commandZ[r][c] = RAD.clampCommandZ ? RAD.clampCommandZ(next, z) : z;
    } else if (event.kind === "group_actuate") {
      for (const cell of event.cells || []) {
        const { r, c } = readCell(next, cell);
        if (next.cells.removed?.[r]?.[c]) continue;
        const alpha = (Number(next.cells.commandAlpha[r][c]) || 0) + (Number(event.alpha) || 0);
        const z = (Number(next.cells.commandZ[r][c]) || 0) + (Number(event.z) || 0);
        next.cells.commandAlpha[r][c] = RAD.clampCommandAlpha ? RAD.clampCommandAlpha(next, alpha) : alpha;
        next.cells.commandZ[r][c] = RAD.clampCommandZ ? RAD.clampCommandZ(next, z) : z;
      }
    } else if (event.kind === "lock") {
      const { r, c } = readCell(next, event.cell);
      if (next.cells.removed?.[r]?.[c]) return next;
      const sim = RAD.simulate(next);
      next.cells.lockAlpha[r][c] = sim.alpha[r][c];
      next.cells.lockZ[r][c] = sim.height[r][c];
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
    } else if (event.kind === "remove_cell") {
      const { r, c } = readCell(next, event.cell);
      next.cells.commandAlpha[r][c] = 0;
      next.cells.commandZ[r][c] = 0;
      next.cells.locked[r][c] = false;
      next.cells.lockAlpha[r][c] = Number(next.grid.initialAlpha) || 1;
      next.cells.lockZ[r][c] = 0;
      next.cells.positionLocked[r][c] = false;
      clearPositionLockTarget(next, r, c);
      next.cells.actuatorAllowed[r][c] = false;
      next.cells.removed[r][c] = true;
    } else if (event.kind === "restore_cell") {
      const { r, c } = readCell(next, event.cell);
      next.cells.removed[r][c] = false;
      if (next.cells.actuatorAllowed?.[r]?.[c] === false) next.cells.actuatorAllowed[r][c] = true;
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
    const rows = Math.max(a?.length || 0, b?.length || 0);
    for (let r = 0; r < rows; r += 1) {
      const rowA = a?.[r] || [];
      const rowB = b?.[r] || [];
      const cols = Math.max(rowA.length, rowB.length);
      for (let c = 0; c < cols; c += 1) {
        out = Math.max(out, Math.abs((Number(rowA[c]) || 0) - (Number(rowB[c]) || 0)));
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
    const lockZError = maxMatrixDiff(a.cells.lockZ, b.cells.lockZ);
    const finalAlphaError = maxMatrixDiff(firstSim.alpha, secondSim.alpha);
    const finalHeightError = maxMatrixDiff(firstSim.height, secondSim.height);
    return {
      modeChanges: boolMatrixDiffCount(a.cells.locked, b.cells.locked),
      topologyChanges: boolMatrixDiffCount(a.cells.removed || [], b.cells.removed || []),
      commandAlphaError,
      commandZError,
      commandError: Math.max(commandAlphaError, commandZError),
      lockAlphaError,
      lockZError,
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
      modeCommutes:
        boolMatrixEqual(firstThenSecond.cells.locked, secondThenFirst.cells.locked) &&
        boolMatrixEqual(firstThenSecond.cells.removed, secondThenFirst.cells.removed),
      commandCommutes:
        matrixEqual(firstThenSecond.cells.commandAlpha, secondThenFirst.cells.commandAlpha, tolerance) &&
        matrixEqual(firstThenSecond.cells.commandZ, secondThenFirst.cells.commandZ, tolerance),
      lockAlphaCommutes: matrixEqual(firstThenSecond.cells.lockAlpha, secondThenFirst.cells.lockAlpha, tolerance),
      lockZCommutes: matrixEqual(firstThenSecond.cells.lockZ, secondThenFirst.cells.lockZ, tolerance),
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
    let maxAdjacentLockZError = 0;
    let noncommutingAdjacentPairs = 0;
    const adjacent = [];
    for (let index = 0; index + 1 < sequence.length; index += 1) {
      const swapped = sequence.map(cloneData);
      [swapped[index], swapped[index + 1]] = [swapped[index + 1], swapped[index]];
      const swappedFinal = applyEventSequence(state, swapped);
      const distance = stateDistance(baseFinal, swappedFinal, baseSim);
      const sensitive =
        distance.modeChanges > 0 ||
        distance.topologyChanges > 0 ||
        distance.commandError > tolerance ||
        distance.lockAlphaError > tolerance ||
        distance.lockZError > tolerance ||
        distance.finalAlphaError > tolerance ||
        distance.finalHeightError > tolerance;
      if (sensitive) noncommutingAdjacentPairs += 1;
      maxAdjacentAlphaError = Math.max(maxAdjacentAlphaError, distance.finalAlphaError);
      maxAdjacentHeightError = Math.max(maxAdjacentHeightError, distance.finalHeightError);
      maxAdjacentCommandError = Math.max(maxAdjacentCommandError, distance.commandError);
      maxAdjacentLockAlphaError = Math.max(maxAdjacentLockAlphaError, distance.lockAlphaError);
      maxAdjacentLockZError = Math.max(maxAdjacentLockZError, distance.lockZError);
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
      reverse.lockZError,
      maxAdjacentAlphaError,
      maxAdjacentHeightError,
      maxAdjacentCommandError,
      maxAdjacentLockAlphaError,
      maxAdjacentLockZError
    );
    return {
      eventCount: sequence.length,
      adjacentPairCount: Math.max(0, sequence.length - 1),
      noncommutingAdjacentPairs,
      orderSensitive:
        reverse.modeChanges > 0 ||
        reverse.topologyChanges > 0 ||
        maxOrderError > tolerance ||
        noncommutingAdjacentPairs > 0,
      reverseAlphaError: reverse.finalAlphaError,
      reverseHeightError: reverse.finalHeightError,
      reverseModeChanges: reverse.modeChanges,
      reverseTopologyChanges: reverse.topologyChanges,
      reverseCommandError: reverse.commandError,
      reverseLockAlphaError: reverse.lockAlphaError,
      reverseLockZError: reverse.lockZError,
      maxAdjacentAlphaError,
      maxAdjacentHeightError,
      maxAdjacentCommandError,
      maxAdjacentLockAlphaError,
      maxAdjacentLockZError,
      maxOrderError,
      adjacent,
    };
  }

  function compareGroupActuationDecomposition(state, cells, alpha = 0, z = 0, tolerance = 1e-9) {
    const base = cloneData(state);
    ensureRemoved(base);
    const validated = (cells || []).map((cell) => readCell(base, cell));
    const appliedCells = validated.filter(({ r, c }) => !base.cells.removed?.[r]?.[c]);
    const skippedRemovedCells = validated.filter(({ r, c }) => base.cells.removed?.[r]?.[c]);
    const groupEvent = groupActuationEvent(validated, alpha, z);
    const localEvents = validated.map((cell) => localActuationEvent(cell, alpha, z));
    const simultaneousState = applyEventSequence(base, [groupEvent]);
    const sequencedState = applyEventSequence(base, localEvents);
    const distance = stateDistance(simultaneousState, sequencedState);
    const decomposes =
      distance.modeChanges === 0 &&
      distance.topologyChanges === 0 &&
      distance.commandError <= tolerance &&
      distance.lockAlphaError <= tolerance &&
      distance.finalAlphaError <= tolerance &&
      distance.finalHeightError <= tolerance;
    return {
      cells: validated,
      appliedCells,
      skippedRemovedCells,
      groupEvent,
      localEvents,
      simultaneousState,
      sequencedState,
      decomposes,
      ...distance,
    };
  }

  function cellKey(cell) {
    return `${cell.r},${cell.c}`;
  }

  function cellsAboveThreshold(values, tolerance = 1e-9, exclude = new Set()) {
    const out = [];
    for (let r = 0; r < values.length; r += 1) {
      for (let c = 0; c < values[r].length; c += 1) {
        if (exclude.has(`${r},${c}`)) continue;
        if (Math.abs(Number(values[r][c]) || 0) > tolerance) out.push({ r, c });
      }
    }
    return out;
  }

  function maxAbsMatrix(values) {
    let out = 0;
    for (const row of values || []) {
      for (const value of row || []) out = Math.max(out, Math.abs(Number(value) || 0));
    }
    return out;
  }

  function verticalContactPenalty(state, cells, contactStiffness = 1, tolerance = 1e-9) {
    const clearance = typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : 0;
    let energy = 0;
    const engagedCells = [];
    for (const cell of cells || []) {
      const { r, c } = readCell(state, cell);
      if (state.cells.removed?.[r]?.[c]) continue;
      const command = Math.abs(Number(state.cells.commandZ?.[r]?.[c]) || 0);
      const penetration = Math.max(0, command - clearance);
      if (penetration > tolerance) engagedCells.push({ r, c });
      energy += 0.5 * Math.max(0, Number(contactStiffness) || 0) * penetration ** 2;
    }
    return { energy, engagedCells };
  }

  function loadActiveCells(state, fixedCells = []) {
    const fixed = new Set((fixedCells || []).map(cellKey));
    const out = [];
    for (let r = 0; r < state.grid.rows; r += 1) {
      for (let c = 0; c < state.grid.cols; c += 1) {
        if (fixed.has(`${r},${c}`)) continue;
        if (state.cells.removed?.[r]?.[c]) continue;
        if (state.cells.locked?.[r]?.[c]) continue;
        out.push({ r, c });
      }
    }
    return out;
  }

  function heightLoadContactMetrics(state, height, fixedCells, externalZLoad, contactStiffness, tolerance) {
    const clearance = typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : 0;
    const activeCells = loadActiveCells(state, fixedCells);
    const force = Number(externalZLoad) || 0;
    const stiffness = Math.max(0, Number(contactStiffness) || 0);
    let loadWork = 0;
    let loadWorkMagnitude = 0;
    let contactPenaltyEnergy = 0;
    const contactEngagedCells = [];
    for (const { r, c } of activeCells) {
      const displacement = Number(height?.[r]?.[c]) || 0;
      loadWork += -force * displacement;
      loadWorkMagnitude += Math.abs(force) * Math.abs(displacement);
      const penetration = Math.max(0, Math.abs(displacement) - clearance);
      if (penetration > tolerance) contactEngagedCells.push({ r, c });
      contactPenaltyEnergy += 0.5 * stiffness * penetration ** 2;
    }
    return { activeCells, loadWork, loadWorkMagnitude, contactPenaltyEnergy, contactEngagedCells };
  }

  function compareGroupActuationUnderRemoval(state, cells, removedCells, alpha = 0, z = 0, tolerance = 1e-9) {
    const base = cloneData(state);
    ensureRemoved(base);
    const validated = (cells || []).map((cell) => readCell(base, cell));
    const removed = (removedCells || []).map((cell) => readCell(base, cell));
    const removalEvents = removed.map((cell) => removeCellEvent(cell));
    const intact = compareGroupActuationDecomposition(base, validated, alpha, z, tolerance);
    const removedBase = applyEventSequence(base, removalEvents);
    const afterRemoval = compareGroupActuationDecomposition(removedBase, validated, alpha, z, tolerance);
    const intactApplied = new Set(intact.appliedCells.map(cellKey));
    const afterApplied = new Set(afterRemoval.appliedCells.map(cellKey));
    const intactSkipped = new Set(intact.skippedRemovedCells.map(cellKey));
    const afterSkipped = new Set(afterRemoval.skippedRemovedCells.map(cellKey));
    const intactTopology = typeof RAD.topologyDiagnostics === "function" ? RAD.topologyDiagnostics(base) : null;
    const removedTopology = typeof RAD.topologyDiagnostics === "function" ? RAD.topologyDiagnostics(removedBase) : null;
    const distance = stateDistance(intact.simultaneousState, afterRemoval.simultaneousState);
    return {
      cells: validated,
      removedCells: removed,
      removalEvents,
      intact,
      removed: afterRemoval,
      lostAppliedCells: validated.filter((cell) => intactApplied.has(cellKey(cell)) && !afterApplied.has(cellKey(cell))),
      newlySkippedCells: validated.filter((cell) => afterSkipped.has(cellKey(cell)) && !intactSkipped.has(cellKey(cell))),
      intactTopology,
      removedTopology,
      componentCountDelta: (removedTopology?.componentCount || 0) - (intactTopology?.componentCount || 0),
      deletedEdgeDelta: (removedTopology?.deletedEdges || 0) - (intactTopology?.deletedEdges || 0),
      ...distance,
    };
  }

  function compareVerticalResidualUnderRemoval(
    state,
    cells,
    removedCells,
    z,
    alpha = 0,
    tolerance = 1e-9,
    contactStiffness = 1,
    options = {}
  ) {
    const base = cloneData(state);
    ensureRemoved(base);
    const validated = (cells || []).map((cell) => readCell(base, cell));
    const removed = (removedCells || []).map((cell) => readCell(base, cell));
    const fixedCells = (options.fixedCells || []).map((cell) => readCell(base, cell));
    const externalZLoad = Number(options.externalZLoad) || 0;
    const sourceKeys = new Set(validated.map(cellKey));
    const removalEvents = removed.map((cell) => removeCellEvent(cell));
    const groupEvent = groupActuationEvent(validated, alpha, z);
    const intactState = applyEventSequence(base, [groupEvent]);
    const removedBase = applyEventSequence(base, removalEvents);
    const removedState = applyEventSequence(removedBase, [groupEvent]);
    const intactSim = RAD.simulate(intactState);
    const removedSim = RAD.simulate(removedState);
    const residualDelta = RAD.matrix(base.grid.rows, base.grid.cols, (r, c) =>
      (Number(removedSim.zResidual?.[r]?.[c]) || 0) - (Number(intactSim.zResidual?.[r]?.[c]) || 0)
    );
    const intactReach = cellsAboveThreshold(intactSim.zResidual, tolerance);
    const removedReach = cellsAboveThreshold(removedSim.zResidual, tolerance);
    const intactNeighbors = cellsAboveThreshold(intactSim.zResidual, tolerance, sourceKeys);
    const removedNeighbors = cellsAboveThreshold(removedSim.zResidual, tolerance, sourceKeys);
    const lostResidualCells = intactReach.filter(({ r, c }) =>
      Math.abs(Number(removedSim.zResidual?.[r]?.[c]) || 0) <= tolerance
    );
    const topologyBlockedCells = lostResidualCells.filter(({ r, c }) => !removedState.cells.removed?.[r]?.[c]);
    const intactContact = verticalContactPenalty(intactState, validated, contactStiffness, tolerance);
    const removedContact = verticalContactPenalty(removedState, validated, contactStiffness, tolerance);
    const intactLoad = heightLoadContactMetrics(
      intactState,
      intactSim.height,
      fixedCells,
      externalZLoad,
      contactStiffness,
      tolerance
    );
    const removedLoad = heightLoadContactMetrics(
      removedState,
      removedSim.height,
      fixedCells,
      externalZLoad,
      contactStiffness,
      tolerance
    );
    const intactTopology = typeof RAD.topologyDiagnostics === "function" ? RAD.topologyDiagnostics(base) : null;
    const removedTopology = typeof RAD.topologyDiagnostics === "function" ? RAD.topologyDiagnostics(removedBase) : null;
    const maxAbsNeighborResidual = removedNeighbors.reduce(
      (max, { r, c }) => Math.max(max, Math.abs(Number(removedSim.zResidual?.[r]?.[c]) || 0)),
      0
    );
    return {
      cells: validated,
      removedCells: removed,
      fixedCells,
      removalEvents,
      intactState,
      removedState,
      intactTopology,
      removedTopology,
      intactZResidual: intactSim.zResidual,
      removedZResidual: removedSim.zResidual,
      residualDelta,
      intactZDieOff: intactSim.zDieOff,
      removedZDieOff: removedSim.zDieOff,
      affectedNeighborCells: removedNeighbors,
      lostResidualCells,
      topologyBlockedCells,
      removedSourceCells: validated.filter(({ r, c }) => removedState.cells.removed?.[r]?.[c]),
      intactReachCount: intactReach.length,
      removedReachCount: removedReach.length,
      intactNeighborReachCount: intactNeighbors.length,
      removedNeighborReachCount: removedNeighbors.length,
      intactDieOffRadius: finiteDieOffRadius(intactSim.zDieOff),
      removedDieOffRadius: finiteDieOffRadius(removedSim.zDieOff),
      maxAbsNeighborResidual,
      maxAbsResidualDelta: maxAbsMatrix(residualDelta),
      clearance: typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(base) : 0,
      contactStiffness: Math.max(0, Number(contactStiffness) || 0),
      intactContactEngagedCells: intactContact.engagedCells,
      removedContactEngagedCells: removedContact.engagedCells,
      intactContactPenaltyEnergy: intactContact.energy,
      removedContactPenaltyEnergy: removedContact.energy,
      contactPenaltyDelta: removedContact.energy - intactContact.energy,
      externalZLoad,
      intactLoadActiveCells: intactLoad.activeCells,
      removedLoadActiveCells: removedLoad.activeCells,
      intactLoadWork: intactLoad.loadWork,
      removedLoadWork: removedLoad.loadWork,
      loadWorkDelta: removedLoad.loadWork - intactLoad.loadWork,
      intactLoadWorkMagnitude: intactLoad.loadWorkMagnitude,
      removedLoadWorkMagnitude: removedLoad.loadWorkMagnitude,
      loadWorkMagnitudeDelta: removedLoad.loadWorkMagnitude - intactLoad.loadWorkMagnitude,
      intactHeightContactEngagedCells: intactLoad.contactEngagedCells,
      removedHeightContactEngagedCells: removedLoad.contactEngagedCells,
      intactHeightContactPenaltyEnergy: intactLoad.contactPenaltyEnergy,
      removedHeightContactPenaltyEnergy: removedLoad.contactPenaltyEnergy,
      heightContactPenaltyDelta: removedLoad.contactPenaltyEnergy - intactLoad.contactPenaltyEnergy,
      componentCountDelta: (removedTopology?.componentCount || 0) - (intactTopology?.componentCount || 0),
      deletedEdgeDelta: (removedTopology?.deletedEdges || 0) - (intactTopology?.deletedEdges || 0),
    };
  }

  RAD.finiteDieOffRadius = finiteDieOffRadius;
  RAD.localActuationEvent = localActuationEvent;
  RAD.groupActuationEvent = groupActuationEvent;
  RAD.lockEvent = lockEvent;
  RAD.releaseEvent = releaseEvent;
  RAD.clearActuationEvent = clearActuationEvent;
  RAD.removeCellEvent = removeCellEvent;
  RAD.restoreCellEvent = restoreCellEvent;
  RAD.applyProgrammableEvent = applyProgrammableEvent;
  RAD.applyEventSequence = applyEventSequence;
  RAD.compareEventOrder = compareEventOrder;
  RAD.compareSequenceOrder = compareSequenceOrder;
  RAD.compareGroupActuationDecomposition = compareGroupActuationDecomposition;
  RAD.compareGroupActuationUnderRemoval = compareGroupActuationUnderRemoval;
  RAD.compareVerticalResidualUnderRemoval = compareVerticalResidualUnderRemoval;
})();
