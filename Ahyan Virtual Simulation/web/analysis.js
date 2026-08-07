(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  function sequenceFrames(state) {
    const initial = state.experiment.initialSnapshot || RAD.snapshotState(state);
    const frames = [
      {
        index: 0,
        type: "initial",
        name: "initial",
        at: "",
        snapshot: initial,
      },
    ];
    for (const [eventIndex, event] of (state.experiment.eventList || []).entries()) {
      if (!event.snapshot) continue;
      frames.push({
        index: eventIndex + 1,
        type: event.type || "event",
        name: event.name || event.target || event.previewType || "",
        at: event.at || "",
        snapshot: event.snapshot,
      });
    }
    return frames;
  }

  function simulateSnapshot(snapshot) {
    const rows = snapshot?.grid?.rows || 1;
    const cols = snapshot?.grid?.cols || 1;
    const temp = RAD.createState(rows, cols);
    RAD.restoreSnapshot(temp, snapshot);
    return { state: temp, sim: RAD.simulate(temp) };
  }

  function frameCommandStats(snapshot) {
    const commandAlpha = snapshot?.cells?.commandAlpha || [];
    const commandZ = snapshot?.cells?.commandZ || [];
    const locked = snapshot?.cells?.locked || [];
    const actuatorAllowed = snapshot?.cells?.actuatorAllowed || [];
    let active = 0;
    let lockedCount = 0;
    let allowedCount = 0;
    let meanAbsAlphaCommand = 0;
    let meanAbsZCommand = 0;
    let maxAbsAlphaCommand = 0;
    let maxAbsZCommand = 0;
    let cells = 0;
    for (let r = 0; r < commandAlpha.length; r += 1) {
      for (let c = 0; c < commandAlpha[r].length; c += 1) {
        cells += 1;
        const alpha = Math.abs(Number(commandAlpha[r][c]) || 0);
        const z = Math.abs(Number(commandZ?.[r]?.[c]) || 0);
        if (alpha > 1e-9 || z > 1e-9) active += 1;
        if (locked?.[r]?.[c]) lockedCount += 1;
        if (actuatorAllowed?.[r]?.[c] !== false) allowedCount += 1;
        meanAbsAlphaCommand += alpha;
        meanAbsZCommand += z;
        maxAbsAlphaCommand = Math.max(maxAbsAlphaCommand, alpha);
        maxAbsZCommand = Math.max(maxAbsZCommand, z);
      }
    }
    return {
      active,
      lockedCount,
      allowedCount,
      meanAbsAlphaCommand: cells ? meanAbsAlphaCommand / cells : 0,
      meanAbsZCommand: cells ? meanAbsZCommand / cells : 0,
      maxAbsAlphaCommand,
      maxAbsZCommand,
    };
  }

  function frameDelta(prevSnapshot, snapshot) {
    if (!prevSnapshot) return { commandDelta: 0, heightDelta: 0 };
    const { sim: prevSim } = simulateSnapshot(prevSnapshot);
    const { sim } = simulateSnapshot(snapshot);
    const rows = Math.min(prevSnapshot.grid.rows, snapshot.grid.rows);
    const cols = Math.min(prevSnapshot.grid.cols, snapshot.grid.cols);
    let commandDelta = 0;
    let heightDelta = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        commandDelta += Math.abs((snapshot.cells.commandAlpha?.[r]?.[c] || 0) - (prevSnapshot.cells.commandAlpha?.[r]?.[c] || 0));
        commandDelta += Math.abs((snapshot.cells.commandZ?.[r]?.[c] || 0) - (prevSnapshot.cells.commandZ?.[r]?.[c] || 0));
        heightDelta += Math.abs((sim.height?.[r]?.[c] || 0) - (prevSim.height?.[r]?.[c] || 0));
      }
    }
    const cells = Math.max(1, rows * cols);
    return { commandDelta: commandDelta / cells, heightDelta: heightDelta / cells };
  }

  function worstTargetResidual(sim) {
    let row = 0;
    let col = 0;
    let value = 0;
    let abs = 0;
    for (let r = 0; r < sim.targetError.length; r += 1) {
      for (let c = 0; c < sim.targetError[r].length; c += 1) {
        const residual = Number(sim.targetError[r][c]) || 0;
        if (Math.abs(residual) > abs) {
          row = r;
          col = c;
          value = residual;
          abs = Math.abs(residual);
        }
      }
    }
    return { row, col, value, abs };
  }

  function clampIndex(value, max) {
    return Math.max(0, Math.min(max - 1, Number(value) || 0));
  }

  function cellKey(cell) {
    return `${cell.r},${cell.c}`;
  }

  function characterizationCells(state, r, c, scope) {
    const rows = state.grid.rows;
    const cols = state.grid.cols;
    const selected = { r: clampIndex(r, rows), c: clampIndex(c, cols) };
    if (scope === "lattice") {
      const active = [];
      for (let rr = 0; rr < rows; rr += 1) {
        for (let cc = 0; cc < cols; cc += 1) {
          if (Math.abs(state.cells.commandAlpha?.[rr]?.[cc] || 0) > 1e-9 || Math.abs(state.cells.commandZ?.[rr]?.[cc] || 0) > 1e-9) {
            active.push({ r: rr, c: cc });
          }
        }
      }
      return active.length ? active : [selected];
    }
    if (scope === "cluster") return RAD.brushCells(state, selected.r, selected.c, 1);
    if (scope === "pair") {
      const candidates = [
        { r: selected.r, c: selected.c + 1 },
        { r: selected.r, c: selected.c - 1 },
        { r: selected.r + 1, c: selected.c },
        { r: selected.r - 1, c: selected.c },
      ].filter((cell) => cell.r >= 0 && cell.r < rows && cell.c >= 0 && cell.c < cols);
      return [selected, candidates[0] || selected].filter((cell, index, cells) => cells.findIndex((other) => cellKey(other) === cellKey(cell)) === index);
    }
    return [selected];
  }

  function cloneForCharacterization(state) {
    const temp = RAD.createState(state.grid.rows, state.grid.cols);
    RAD.restoreSnapshot(temp, RAD.snapshotState(state));
    return temp;
  }

  function scopedState(state, cells) {
    const temp = cloneForCharacterization(state);
    const keep = new Set(cells.map(cellKey));
    for (let r = 0; r < temp.grid.rows; r += 1) {
      for (let c = 0; c < temp.grid.cols; c += 1) {
        if (!keep.has(`${r},${c}`)) {
          temp.cells.commandAlpha[r][c] = 0;
          temp.cells.commandZ[r][c] = 0;
        }
      }
    }
    return temp;
  }

  function finiteMax(field) {
    let value = 0;
    for (const row of field || []) {
      for (const entry of row || []) {
        if (Number.isFinite(entry)) value = Math.max(value, entry);
      }
    }
    return value;
  }

  function localResponseStats(state, sim, baselineSim, sourceCells) {
    const sourceSet = new Set(sourceCells.map(cellKey));
    let responseCells = 0;
    let alphaReachCells = 0;
    let zReachCells = 0;
    let maxAlphaDelta = 0;
    let maxHeightDelta = 0;
    let meanAbsAlphaDelta = 0;
    let meanAbsHeightDelta = 0;
    let activeSources = 0;
    const totalCells = Math.max(1, state.grid.rows * state.grid.cols);
    for (const cell of sourceCells) {
      if (Math.abs(state.cells.commandAlpha?.[cell.r]?.[cell.c] || 0) > 1e-9 || Math.abs(state.cells.commandZ?.[cell.r]?.[cell.c] || 0) > 1e-9) {
        activeSources += 1;
      }
    }
    for (let r = 0; r < state.grid.rows; r += 1) {
      for (let c = 0; c < state.grid.cols; c += 1) {
        const alphaDelta = (sim.alpha?.[r]?.[c] || 0) - (baselineSim.alpha?.[r]?.[c] || 0);
        const heightDelta = (sim.height?.[r]?.[c] || 0) - (baselineSim.height?.[r]?.[c] || 0);
        const alphaAbs = Math.abs(alphaDelta);
        const heightAbs = Math.abs(heightDelta);
        const isSource = sourceSet.has(`${r},${c}`);
        if (alphaAbs > 1e-8 || heightAbs > 1e-8) responseCells += 1;
        if (!isSource && Math.abs(sim.influence?.[r]?.[c] || 0) > 1e-8) alphaReachCells += 1;
        if (!isSource && Math.abs(sim.zResidual?.[r]?.[c] || 0) > 1e-8) zReachCells += 1;
        maxAlphaDelta = Math.max(maxAlphaDelta, alphaAbs);
        maxHeightDelta = Math.max(maxHeightDelta, heightAbs);
        meanAbsAlphaDelta += alphaAbs;
        meanAbsHeightDelta += heightAbs;
      }
    }
    return {
      activeSources,
      responseCells,
      alphaReachCells,
      zReachCells,
      maxAlphaDelta,
      maxHeightDelta,
      meanAbsAlphaDelta: meanAbsAlphaDelta / totalCells,
      meanAbsHeightDelta: meanAbsHeightDelta / totalCells,
      alphaDieOff: finiteMax(sim.dieOff),
      zDieOff: finiteMax(sim.zDieOff),
    };
  }

  function superpositionError(state, sourceCells, combinedSim, baselineSim) {
    if (sourceCells.length <= 1) return { rms: 0, max: 0, skipped: false, sourceCount: sourceCells.length };
    const activeSources = sourceCells.filter((cell) => Math.abs(state.cells.commandAlpha?.[cell.r]?.[cell.c] || 0) > 1e-9 || Math.abs(state.cells.commandZ?.[cell.r]?.[cell.c] || 0) > 1e-9);
    if (activeSources.length <= 1) return { rms: 0, max: 0, skipped: false, sourceCount: activeSources.length };
    if (activeSources.length > 16) return { rms: null, max: null, skipped: true, sourceCount: activeSources.length };
    const rows = state.grid.rows;
    const cols = state.grid.cols;
    const alphaSum = RAD.matrix(rows, cols, 0);
    const heightSum = RAD.matrix(rows, cols, 0);
    for (const source of activeSources) {
      const singleState = scopedState(state, [source]);
      const singleSim = RAD.simulate(singleState);
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          alphaSum[r][c] += (singleSim.alpha[r][c] || 0) - (baselineSim.alpha[r][c] || 0);
          heightSum[r][c] += (singleSim.height[r][c] || 0) - (baselineSim.height[r][c] || 0);
        }
      }
    }
    let squared = 0;
    let max = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        const alphaResidual = (combinedSim.alpha[r][c] || 0) - (baselineSim.alpha[r][c] || 0) - alphaSum[r][c];
        const heightResidual = (combinedSim.height[r][c] || 0) - (baselineSim.height[r][c] || 0) - heightSum[r][c];
        squared += alphaResidual * alphaResidual + heightResidual * heightResidual;
        max = Math.max(max, Math.abs(alphaResidual), Math.abs(heightResidual));
      }
    }
    return { rms: Math.sqrt(squared / Math.max(1, rows * cols * 2)), max, skipped: false, sourceCount: activeSources.length };
  }

  function responseVector(next, base, field) {
    const values = [];
    for (let r = 0; r < next[field].length; r += 1) {
      for (let c = 0; c < next[field][r].length; c += 1) values.push((next[field][r][c] || 0) - (base[field][r][c] || 0));
    }
    return values;
  }

  function matrixRankFromColumns(columns, tolerance = 1e-9) {
    if (!columns.length) return 0;
    const rows = columns[0].length;
    const matrix = Array.from({ length: rows }, (_, r) => columns.map((column) => Number(column[r]) || 0));
    let rank = 0;
    for (let col = 0; col < columns.length && rank < rows; col += 1) {
      let pivot = rank;
      for (let r = rank + 1; r < rows; r += 1) {
        if (Math.abs(matrix[r][col]) > Math.abs(matrix[pivot][col])) pivot = r;
      }
      if (Math.abs(matrix[pivot][col]) <= tolerance) continue;
      [matrix[rank], matrix[pivot]] = [matrix[pivot], matrix[rank]];
      const scale = matrix[rank][col];
      for (let c = col; c < columns.length; c += 1) matrix[rank][c] /= scale;
      for (let r = 0; r < rows; r += 1) {
        if (r === rank) continue;
        const factor = matrix[r][col];
        if (Math.abs(factor) <= tolerance) continue;
        for (let c = col; c < columns.length; c += 1) matrix[r][c] -= factor * matrix[rank][c];
      }
      rank += 1;
    }
    return rank;
  }

  function countReachableFromColumns(columns, tolerance = 1e-9) {
    if (!columns.length) return 0;
    let count = 0;
    for (let row = 0; row < columns[0].length; row += 1) {
      if (columns.some((column) => Math.abs(column[row] || 0) > tolerance)) count += 1;
    }
    return count;
  }

  function responseMatrixDiagnostic(state, sourceCells, baselineSim, tolerance = 1e-9) {
    const limits = RAD.commandLimits(state);
    const stepAlpha = Math.min(0.12, limits.alphaContract);
    const stepZ = Math.min(0.12, limits.z);
    const alphaColumns = [];
    const heightColumns = [];
    let columnCount = 0;
    for (const cell of sourceCells) {
      if (state.cells.locked?.[cell.r]?.[cell.c]) continue;
      if (stepAlpha > tolerance) {
        const alphaState = scopedState(state, []);
        alphaState.cells.commandAlpha[cell.r][cell.c] = -stepAlpha;
        const alphaSim = RAD.simulate(alphaState);
        alphaColumns.push(responseVector(alphaSim, baselineSim, "alpha"));
        heightColumns.push(responseVector(alphaSim, baselineSim, "height"));
        columnCount += 1;
      }
      if (stepZ > tolerance) {
        const zState = scopedState(state, []);
        zState.cells.commandZ[cell.r][cell.c] = stepZ;
        const zSim = RAD.simulate(zState);
        alphaColumns.push(responseVector(zSim, baselineSim, "alpha"));
        heightColumns.push(responseVector(zSim, baselineSim, "height"));
        columnCount += 1;
      }
    }
    const totalCells = state.grid.rows * state.grid.cols;
    const reachableAlphaCells = countReachableFromColumns(alphaColumns, tolerance);
    const reachableHeightCells = countReachableFromColumns(heightColumns, tolerance);
    return {
      diagnosticColumnCount: columnCount,
      responseRankAlpha: matrixRankFromColumns(alphaColumns, tolerance),
      responseRankHeight: matrixRankFromColumns(heightColumns, tolerance),
      reachableAlphaCells,
      reachableHeightCells,
      alphaUnderactuatedCells: Math.max(0, totalCells - reachableAlphaCells),
      heightUnderactuatedCells: Math.max(0, totalCells - reachableHeightCells),
    };
  }

  function nearestSourceDistance(r, c, sourceCells) {
    let best = Infinity;
    for (const source of sourceCells) best = Math.min(best, Math.abs(r - source.r) + Math.abs(c - source.c));
    return Number.isFinite(best) ? best : 0;
  }

  function fitShellDecay(shells, key, tolerance = 1e-8) {
    const entries = shells
      .filter((shell) => shell[key] > tolerance)
      .sort((a, b) => a.distance - b.distance);
    const first = entries[0]?.[key] || 0;
    const last = entries.at(-1)?.[key] || 0;
    const reach = entries.at(-1)?.distance || 0;
    if (entries.length < 2) {
      return { ratio: 0, length: 0, shells: entries.length, reach, first, last };
    }
    let sumX = 0;
    let sumY = 0;
    let sumXX = 0;
    let sumXY = 0;
    for (const entry of entries) {
      const x = entry.distance;
      const y = Math.log(Math.max(tolerance, entry[key]));
      sumX += x;
      sumY += y;
      sumXX += x * x;
      sumXY += x * y;
    }
    const n = entries.length;
    const denom = n * sumXX - sumX * sumX;
    const slope = Math.abs(denom) > tolerance ? (n * sumXY - sumX * sumY) / denom : 0;
    const ratio = Math.exp(slope);
    const length = slope < -tolerance ? -1 / slope : 0;
    return { ratio, length, shells: entries.length, reach, first, last };
  }

  function responseDecayProfile(state, sim, baselineSim, sourceCells, tolerance = 1e-8) {
    const shellMap = new Map();
    for (let r = 0; r < state.grid.rows; r += 1) {
      for (let c = 0; c < state.grid.cols; c += 1) {
        const distance = nearestSourceDistance(r, c, sourceCells);
        if (!shellMap.has(distance)) {
          shellMap.set(distance, { distance, count: 0, alphaMax: 0, zMax: 0 });
        }
        const shell = shellMap.get(distance);
        const alphaDelta = Math.abs((sim.alpha?.[r]?.[c] || 0) - (baselineSim.alpha?.[r]?.[c] || 0));
        const zDelta = Math.abs((sim.height?.[r]?.[c] || 0) - (baselineSim.height?.[r]?.[c] || 0));
        shell.alphaMax = Math.max(shell.alphaMax, alphaDelta);
        shell.zMax = Math.max(shell.zMax, zDelta);
        shell.count += 1;
      }
    }
    const shells = Array.from(shellMap.values()).sort((a, b) => a.distance - b.distance);
    const alpha = fitShellDecay(shells, "alphaMax", tolerance);
    const z = fitShellDecay(shells, "zMax", tolerance);
    return {
      decayModel: "log-linear shell max",
      alphaDecayRatio: alpha.ratio,
      zDecayRatio: z.ratio,
      alphaDecayLength: alpha.length,
      zDecayLength: z.length,
      alphaDecayShells: alpha.shells,
      zDecayShells: z.shells,
      alphaDecayReach: alpha.reach,
      zDecayReach: z.reach,
      alphaDecayFirst: alpha.first,
      zDecayFirst: z.first,
      alphaDecayLast: alpha.last,
      zDecayLast: z.last,
    };
  }

  function physicalPreviewComparison(state, baselineState, sim, baselineSim) {
    const unavailable = {
      physicalPreviewAvailable: false,
      physicalPreviewSuccess: false,
      physicalHeightRmsError: 0,
      physicalHeightMaxError: 0,
      physicalCenterRmsError: 0,
      physicalCenterMaxError: 0,
      physicalPreviewIterations: 0,
    };
    if (typeof RAD.simulatePhysicalRelaxation !== "function") return unavailable;
    try {
      const physicalBase = RAD.simulatePhysicalRelaxation(baselineState, { baseSim: baselineSim });
      const physical = RAD.simulatePhysicalRelaxation(state, { baseSim: sim });
      const rows = Math.min(state.grid.rows, baselineState.grid.rows);
      const cols = Math.min(state.grid.cols, baselineState.grid.cols);
      let heightSquared = 0;
      let heightMax = 0;
      let centerSquared = 0;
      let centerMax = 0;
      let heightCount = 0;
      let centerCount = 0;
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          const kinematicHeightDelta = (sim.height?.[r]?.[c] || 0) - (baselineSim.height?.[r]?.[c] || 0);
          const physicalHeightDelta = (physical.height?.[r]?.[c] || 0) - (physicalBase.height?.[r]?.[c] || 0);
          const heightError = physicalHeightDelta - kinematicHeightDelta;
          heightSquared += heightError * heightError;
          heightMax = Math.max(heightMax, Math.abs(heightError));
          heightCount += 1;

          const kinematicPoint = sim.centers?.[r]?.[c];
          const kinematicBasePoint = baselineSim.centers?.[r]?.[c];
          const physicalPoint = physical.centers?.[r]?.[c];
          const physicalBasePoint = physicalBase.centers?.[r]?.[c];
          if (kinematicPoint && kinematicBasePoint && physicalPoint && physicalBasePoint) {
            const dx = (physicalPoint.x - physicalBasePoint.x) - (kinematicPoint.x - kinematicBasePoint.x);
            const dy = (physicalPoint.y - physicalBasePoint.y) - (kinematicPoint.y - kinematicBasePoint.y);
            const dz = (physicalPoint.z - physicalBasePoint.z) - (kinematicPoint.z - kinematicBasePoint.z);
            const centerError = Math.hypot(dx, dy, dz);
            centerSquared += centerError * centerError;
            centerMax = Math.max(centerMax, centerError);
            centerCount += 1;
          }
        }
      }
      return {
        physicalPreviewAvailable: true,
        physicalPreviewSuccess: true,
        physicalHeightRmsError: Math.sqrt(heightSquared / Math.max(1, heightCount)),
        physicalHeightMaxError: heightMax,
        physicalCenterRmsError: Math.sqrt(centerSquared / Math.max(1, centerCount)),
        physicalCenterMaxError: centerMax,
        physicalPreviewIterations: Number(physical.metrics?.physicalIterations || physicalBase.metrics?.physicalIterations || 0),
      };
    } catch (error) {
      return {
        ...unavailable,
        physicalPreviewAvailable: true,
        physicalPreviewError: error?.message || String(error),
      };
    }
  }

  function characterizeLocalResponse(state, options = {}) {
    const scope = options.scope || state.experiment?.characterizationScope || "single";
    const selected = {
      r: clampIndex(options.r ?? state.selection?.r ?? 0, state.grid.rows),
      c: clampIndex(options.c ?? state.selection?.c ?? 0, state.grid.cols),
    };
    const sourceCells = characterizationCells(state, selected.r, selected.c, scope);
    const combinedState = scopedState(state, sourceCells);
    const baselineState = scopedState(state, []);
    const sim = RAD.simulate(combinedState);
    const baselineSim = RAD.simulate(baselineState);
    const stats = localResponseStats(combinedState, sim, baselineSim, sourceCells);
    const interaction = superpositionError(combinedState, sourceCells, sim, baselineSim);
    const matrixDiagnostic = responseMatrixDiagnostic(combinedState, sourceCells, baselineSim);
    const decayProfile = responseDecayProfile(combinedState, sim, baselineSim, sourceCells);
    const physicalComparison = physicalPreviewComparison(combinedState, baselineState, sim, baselineSim);
    const calibration = RAD.paperRadCalibration(combinedState);
    return {
      scope,
      selected,
      cells: sourceCells,
      regionCellCount: sourceCells.length,
      ...stats,
      ...matrixDiagnostic,
      ...decayProfile,
      ...physicalComparison,
      superpositionError: interaction.rms,
      maxSuperpositionError: interaction.max,
      superpositionSkipped: interaction.skipped,
      superpositionSources: interaction.sourceCount,
      nonadditive: !interaction.skipped && (interaction.rms || 0) > 1e-9,
      backlash: Number(combinedState.grid.backlash) || 0,
      zDeadZone: RAD.verticalDeadZone(combinedState),
      pinHoleClearance: RAD.pinHoleClearance(combinedState),
      backlashMm: calibration.configuredBacklashMm,
      pinHoleClearanceMm: calibration.pinHoleClearanceMm,
      model: combinedState.view.simulationMode || "kinematic",
    };
  }

  function analyzeExperimentSequence(state) {
    const frames = sequenceFrames(state);
    const frameMetrics = frames.map((frame, frameIndex) => {
      const { sim } = simulateSnapshot(frame.snapshot);
      const command = frameCommandStats(frame.snapshot);
      const delta = frameDelta(frames[frameIndex - 1]?.snapshot, frame.snapshot);
      const worstResidual = worstTargetResidual(sim);
      return {
        index: frame.index,
        type: frame.type,
        name: frame.name,
        at: frame.at,
        rows: frame.snapshot.grid.rows,
        cols: frame.snapshot.grid.cols,
        rmsTargetError: sim.metrics.rmsTargetError,
        maxTargetError: sim.metrics.maxTargetError,
        meanSignedTargetError: sim.metrics.meanSignedTargetError,
        meanAlpha: sim.metrics.meanAlpha,
        maxAbsHeight: sim.metrics.maxAbsHeight,
        minTheta: sim.metrics.minTheta,
        maxTheta: sim.metrics.maxTheta,
        meanSaturation: sim.metrics.meanSaturation,
        maxSaturation: sim.metrics.maxSaturation,
        saturatedActuators: sim.metrics.saturatedActuators,
        meanTravel: sim.metrics.meanTravel,
        activeCells: sim.metrics.activeCells,
        maxFiniteDieOff: sim.metrics.maxFiniteDieOff,
        meanAbsLinkStrain: sim.metrics.meanAbsLinkStrain,
        maxAbsLinkStrain: sim.metrics.maxAbsLinkStrain,
        meanSurfaceSlope: sim.metrics.meanSurfaceSlope,
        maxSurfaceSlope: sim.metrics.maxSurfaceSlope,
        commandDelta: delta.commandDelta,
        heightDelta: delta.heightDelta,
        worstResidual,
        ...command,
      };
    });
    const finalFrame = frameMetrics.at(-1);
    const bestFrame = frameMetrics.reduce((best, frame) => (frame.rmsTargetError < best.rmsTargetError ? frame : best), frameMetrics[0]);
    return {
      frameCount: Math.max(0, frames.length - 1),
      sampleCount: frames.length,
      frames: frameMetrics,
      finalFrame,
      bestFrame,
    };
  }

  function csvValue(value) {
    if (typeof value === "number") return Number.isFinite(value) ? String(value) : "";
    const text = String(value ?? "");
    return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
  }

  function exportSequenceMetricsCsv(state) {
    const frames = sequenceFrames(state);
    const header = [
      "frame",
      "type",
      "name",
      "at",
      "row",
      "col",
      "command_alpha",
      "command_z",
      "alpha",
      "theta_deg",
      "height_z",
      "target_z",
      "target_error",
      "influence",
      "die_off",
      "saturation",
      "locked",
      "actuator_allowed",
      "link_strain_max",
      "surface_slope",
      "normal_tilt_deg",
      "frame_worst_row",
      "frame_worst_col",
      "frame_worst_residual",
    ];
    const lines = [header.join(",")];
    for (const frame of frames) {
      const { state: temp, sim } = simulateSnapshot(frame.snapshot);
      const worstResidual = worstTargetResidual(sim);
      const rows = temp.grid.rows;
      const cols = temp.grid.cols;
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          const strain = sim.linkStrain.horizontal?.[r]?.[c] ?? sim.linkStrain.vertical?.[r]?.[c] ?? 0;
          lines.push(
            [
              frame.index,
              frame.type,
              frame.name,
              frame.at,
              r,
              c,
              temp.cells.commandAlpha[r][c],
              temp.cells.commandZ[r][c],
              sim.alpha[r][c],
              sim.theta[r][c],
              sim.height[r][c],
              sim.target[r][c],
              sim.targetError[r][c],
              sim.influence[r][c],
              sim.dieOff[r][c],
              sim.saturation[r][c],
              temp.cells.locked[r][c] ? 1 : 0,
              temp.cells.actuatorAllowed?.[r]?.[c] === false ? 0 : 1,
              strain,
              sim.slope.magnitude[r][c],
              sim.slope.tilt[r][c],
              worstResidual.row,
              worstResidual.col,
              worstResidual.value,
            ].map(csvValue).join(",")
          );
        }
      }
    }
    return lines.join("\n");
  }

  RAD.sequenceFrames = sequenceFrames;
  RAD.analyzeExperimentSequence = analyzeExperimentSequence;
  RAD.exportSequenceMetricsCsv = exportSequenceMetricsCsv;
  RAD.characterizeLocalResponse = characterizeLocalResponse;
  RAD.physicalPreviewComparison = physicalPreviewComparison;
  RAD.responseDecayProfile = responseDecayProfile;
})();
