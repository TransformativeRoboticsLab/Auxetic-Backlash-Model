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

  function pairwiseInteractionGraph(state, sourceCells, baselineSim, tolerance = 1e-9, maxPairs = 64) {
    const activeSources = sourceCells.filter((cell) => Math.abs(state.cells.commandAlpha?.[cell.r]?.[cell.c] || 0) > tolerance || Math.abs(state.cells.commandZ?.[cell.r]?.[cell.c] || 0) > tolerance);
    const totalPairCount = (activeSources.length * (activeSources.length - 1)) / 2;
    const alphaErrorMatrix = RAD.matrix(activeSources.length, activeSources.length, 0);
    const heightErrorMatrix = RAD.matrix(activeSources.length, activeSources.length, 0);
    const hotspotMap = RAD.matrix(state.grid.rows, state.grid.cols, 0);
    const degreeMap = RAD.matrix(state.grid.rows, state.grid.cols, 0);
    let hotspotMax = 0;
    let degreeMax = 0;
    const markHotspot = (source, value) => {
      const r = source?.r;
      const c = source?.c;
      if (!hotspotMap[r] || hotspotMap[r][c] === undefined) return;
      const strength = Math.abs(Number(value) || 0);
      hotspotMap[r][c] = Math.max(hotspotMap[r][c], strength);
      hotspotMax = Math.max(hotspotMax, hotspotMap[r][c]);
    };
    const markNonadditiveDegree = (...sources) => {
      const seen = new Set();
      for (const source of sources) {
        const r = source?.r;
        const c = source?.c;
        const key = `${r},${c}`;
        if (seen.has(key) || !degreeMap[r] || degreeMap[r][c] === undefined) continue;
        seen.add(key);
        degreeMap[r][c] += 1;
        degreeMax = Math.max(degreeMax, degreeMap[r][c]);
      }
    };
    if (activeSources.length < 2) {
      return {
        pairwiseInteractionModel: "pairwise superposition residual",
        pairwiseTotalPairs: totalPairCount,
        pairwiseEvaluatedPairs: 0,
        pairwiseNonadditivePairs: 0,
        pairwiseTruncated: false,
        pairwiseMaxAlphaError: 0,
        pairwiseMaxHeightError: 0,
        pairwiseMaxInteractionError: 0,
        pairwiseInteractions: [],
        pairwiseAlphaErrorMatrix: alphaErrorMatrix,
        pairwiseHeightErrorMatrix: heightErrorMatrix,
        pairwiseInteractionMap: hotspotMap,
        pairwiseInteractionMapMax: hotspotMax,
        pairwiseInteractionDegreeMap: degreeMap,
        pairwiseInteractionDegreeMax: degreeMax,
        pairwiseInteractionDensity: 0,
      };
    }

    const singles = activeSources.map((source) => {
      const singleState = scopedState(state, [source]);
      const singleSim = RAD.simulate(singleState);
      return {
        source,
        alphaDelta: responseVector(singleSim, baselineSim, "alpha"),
        heightDelta: responseVector(singleSim, baselineSim, "height"),
      };
    });
    const interactions = [];
    let evaluated = 0;
    let nonadditive = 0;
    let maxAlpha = 0;
    let maxHeight = 0;
    const limit = Math.max(0, Math.min(totalPairCount, Number(maxPairs) || 0));

    for (let i = 0; i < activeSources.length; i += 1) {
      for (let j = i + 1; j < activeSources.length; j += 1) {
        if (evaluated >= limit) {
          interactions.sort((a, b) => b.maxError - a.maxError);
          return {
            pairwiseInteractionModel: "pairwise superposition residual",
            pairwiseTotalPairs: totalPairCount,
            pairwiseEvaluatedPairs: evaluated,
            pairwiseNonadditivePairs: nonadditive,
            pairwiseTruncated: true,
            pairwiseMaxAlphaError: maxAlpha,
            pairwiseMaxHeightError: maxHeight,
            pairwiseMaxInteractionError: Math.max(maxAlpha, maxHeight),
            pairwiseInteractions: interactions.slice(0, 16),
            pairwiseAlphaErrorMatrix: alphaErrorMatrix,
            pairwiseHeightErrorMatrix: heightErrorMatrix,
            pairwiseInteractionMap: hotspotMap,
            pairwiseInteractionMapMax: hotspotMax,
            pairwiseInteractionDegreeMap: degreeMap,
            pairwiseInteractionDegreeMax: degreeMax,
            pairwiseInteractionDensity: evaluated ? nonadditive / evaluated : 0,
          };
        }
        const combinedState = scopedState(state, [activeSources[i], activeSources[j]]);
        const combinedSim = RAD.simulate(combinedState);
        const combinedAlpha = responseVector(combinedSim, baselineSim, "alpha");
        const combinedHeight = responseVector(combinedSim, baselineSim, "height");
        let alphaError = 0;
        let heightError = 0;
        for (let k = 0; k < combinedAlpha.length; k += 1) {
          alphaError = Math.max(alphaError, Math.abs(combinedAlpha[k] - singles[i].alphaDelta[k] - singles[j].alphaDelta[k]));
          heightError = Math.max(heightError, Math.abs(combinedHeight[k] - singles[i].heightDelta[k] - singles[j].heightDelta[k]));
        }
        const maxError = Math.max(alphaError, heightError);
        alphaErrorMatrix[i][j] = alphaError;
        alphaErrorMatrix[j][i] = alphaError;
        heightErrorMatrix[i][j] = heightError;
        heightErrorMatrix[j][i] = heightError;
        maxAlpha = Math.max(maxAlpha, alphaError);
        maxHeight = Math.max(maxHeight, heightError);
        markHotspot(activeSources[i], maxError);
        markHotspot(activeSources[j], maxError);
        if (maxError > tolerance) {
          nonadditive += 1;
          markNonadditiveDegree(activeSources[i], activeSources[j]);
        }
        interactions.push({
          firstIndex: i,
          secondIndex: j,
          first: activeSources[i],
          second: activeSources[j],
          manhattanDistance: Math.abs(activeSources[i].r - activeSources[j].r) + Math.abs(activeSources[i].c - activeSources[j].c),
          alphaSuperpositionError: alphaError,
          heightSuperpositionError: heightError,
          maxError,
          nonadditive: maxError > tolerance,
        });
        evaluated += 1;
      }
    }

    interactions.sort((a, b) => b.maxError - a.maxError);
    return {
      pairwiseInteractionModel: "pairwise superposition residual",
      pairwiseTotalPairs: totalPairCount,
      pairwiseEvaluatedPairs: evaluated,
      pairwiseNonadditivePairs: nonadditive,
      pairwiseTruncated: false,
      pairwiseMaxAlphaError: maxAlpha,
      pairwiseMaxHeightError: maxHeight,
      pairwiseMaxInteractionError: Math.max(maxAlpha, maxHeight),
      pairwiseInteractions: interactions.slice(0, 16),
      pairwiseAlphaErrorMatrix: alphaErrorMatrix,
      pairwiseHeightErrorMatrix: heightErrorMatrix,
      pairwiseInteractionMap: hotspotMap,
      pairwiseInteractionMapMax: hotspotMax,
      pairwiseInteractionDegreeMap: degreeMap,
      pairwiseInteractionDegreeMax: degreeMax,
      pairwiseInteractionDensity: evaluated ? nonadditive / evaluated : 0,
    };
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
    const pairwiseInteractions = pairwiseInteractionGraph(combinedState, sourceCells, baselineSim);
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
      ...pairwiseInteractions,
      backlash: Number(combinedState.grid.backlash) || 0,
      zDeadZone: RAD.verticalDeadZone(combinedState),
      pinHoleClearance: RAD.pinHoleClearance(combinedState),
      backlashMm: calibration.configuredBacklashMm,
      pinHoleClearanceMm: calibration.pinHoleClearanceMm,
      model: combinedState.view.simulationMode || "kinematic",
    };
  }

  const PROTOCOL_MEASUREMENT_FIELDS = Object.freeze([
    "alpha_delta_grid",
    "height_delta_grid",
    "center_displacement_grid",
    "actuator_command",
    "lock_state",
    "pin_hole_slip_mm",
    "actuator_force_n",
  ]);

  function protocolCell(cell) {
    return { row: cell.r, col: cell.c };
  }

  function protocolCommand(cell, alpha = 0, z = 0) {
    return { row: cell.r, col: cell.c, alpha, z };
  }

  function uniqueProtocolCells(cells) {
    const seen = new Set();
    const unique = [];
    for (const cell of cells) {
      const key = cellKey(cell);
      if (seen.has(key)) continue;
      seen.add(key);
      unique.push(cell);
    }
    return unique;
  }

  function calibrationExperimentProtocol(state, options = {}) {
    const limits = typeof RAD.commandLimits === "function" ? RAD.commandLimits(state) : { alphaContract: 0.55, z: 0.8 };
    const center = {
      r: clampIndex(options.r ?? state.selection?.r ?? 0, state.grid.rows),
      c: clampIndex(options.c ?? state.selection?.c ?? 0, state.grid.cols),
    };
    const pairCells = characterizationCells(state, center.r, center.c, "pair");
    const clusterCells = characterizationCells(state, center.r, center.c, "cluster");
    const primaryNeighbor = pairCells.find((cell) => cellKey(cell) !== cellKey(center)) || center;
    const secondaryNeighbor =
      clusterCells.find((cell) => cellKey(cell) !== cellKey(center) && cellKey(cell) !== cellKey(primaryNeighbor)) || primaryNeighbor;
    const alphaStep = Number(options.alphaStep ?? -Math.min(0.25, limits.alphaContract || 0.25));
    const alphaExpand = Math.abs(alphaStep) * 0.75;
    const zStep = Number(options.zStep ?? Math.min(0.3, limits.z || 0.3));
    const repeatCount = Math.max(1, Math.round(Number(options.repeatCount ?? 3)));
    const observationPair = uniqueProtocolCells([center, primaryNeighbor]);
    const observationCluster = uniqueProtocolCells([center, primaryNeighbor, secondaryNeighbor, ...clusterCells]);
    const step = (id, scope, commands, observationCells, lockedCells, purpose, expectedResponse) => ({
      id,
      scope,
      commands,
      observationCells: observationCells.map(protocolCell),
      lockedCells: lockedCells.map(protocolCell),
      measurementFields: [...PROTOCOL_MEASUREMENT_FIELDS],
      purpose,
      expectedResponse,
      repeatCount,
    });
    const steps = [
      step(
        "single_alpha_contract",
        "single",
        [protocolCommand(center, alphaStep, 0)],
        [center],
        [],
        "Measure the local rotating-square dilation response to contraction.",
        "Primary alpha change at the commanded cell with backlash-gated neighbor influence."
      ),
      step(
        "single_alpha_expand",
        "single",
        [protocolCommand(center, alphaExpand, 0)],
        [center],
        [],
        "Measure expansion-side travel and check for asymmetric backlash.",
        "Positive alpha response at the commanded cell with smaller expansion command."
      ),
      step(
        "single_z_lift",
        "single",
        [protocolCommand(center, 0, zStep)],
        observationPair,
        [],
        "Measure direct vertical actuation and residual neighbor lift.",
        "Commanded cell moves vertically; adjacent observation cell captures pin-hole residual coupling."
      ),
      step(
        "pair_z_residual",
        "pair",
        [protocolCommand(center, 0, zStep)],
        observationPair,
        [],
        "Quantify vertical die-off from one actuated cell into a neighboring cell.",
        "Neighbor height response should decay with clearance, backlash, and graph distance."
      ),
      step(
        "pair_superposition",
        "pair",
        [protocolCommand(center, alphaStep, 0.5 * zStep), protocolCommand(primaryNeighbor, alphaStep, 0.5 * zStep)],
        observationPair,
        [],
        "Measure whether adjacent cell commands add linearly or interact through backlash.",
        "Any deviation from summed single-cell responses identifies a programmable-discontinuity interaction."
      ),
      step(
        "cluster_mixed_actuation",
        "cluster",
        [
          protocolCommand(center, 0, zStep),
          protocolCommand(primaryNeighbor, alphaStep, 0),
          protocolCommand(secondaryNeighbor, 0.5 * alphaExpand, -0.5 * zStep),
        ],
        observationCluster,
        [],
        "Measure collective response of mixed horizontal and vertical actuation.",
        "Cluster field should reveal multi-operator coupling, residual height spread, and reachable directions."
      ),
      step(
        "locked_cell_control",
        "lock",
        [protocolCommand(center, alphaStep, zStep)],
        observationPair,
        [center],
        "Verify lock enforcement against commanded alpha and vertical motion.",
        "Locked cell should remain fixed while any neighbor residual exposes compliance leakage."
      ),
    ];
    const profile = typeof RAD.hardwareProfile === "function" ? RAD.hardwareProfile(state) : { name: "paper-reference" };
    return {
      schema: "rad-sim.calibration-experiment-protocol.v1",
      hardwareProfile: profile.name || "paper-reference",
      readiness: typeof RAD.calibrationReadiness === "function" ? RAD.calibrationReadiness(state) : null,
      measurementPlan: typeof RAD.calibrationMeasurementPlan === "function" ? RAD.calibrationMeasurementPlan(state) : [],
      grid: { rows: state.grid.rows, cols: state.grid.cols },
      centerCell: protocolCell(center),
      measurementFields: [...PROTOCOL_MEASUREMENT_FIELDS],
      notes: "Protocol defines repeatable simulator/bench measurements; it does not claim the current spring-hinge solver is calibrated.",
      steps,
    };
  }

  function exportCalibrationExperimentProtocol(state, options = {}) {
    return JSON.stringify(calibrationExperimentProtocol(state, options), null, 2);
  }

  function protocolStepState(state, step, includeCommands = true) {
    const temp = cloneForCharacterization(state);
    for (let r = 0; r < temp.grid.rows; r += 1) {
      for (let c = 0; c < temp.grid.cols; c += 1) {
        temp.cells.commandAlpha[r][c] = 0;
        temp.cells.commandZ[r][c] = 0;
        temp.cells.locked[r][c] = false;
      }
    }
    for (const locked of step.lockedCells || []) {
      if (temp.cells.locked?.[locked.row]?.[locked.col] !== undefined) temp.cells.locked[locked.row][locked.col] = true;
    }
    if (includeCommands) {
      for (const command of step.commands || []) {
        if (temp.cells.commandAlpha?.[command.row]?.[command.col] === undefined) continue;
        temp.cells.commandAlpha[command.row][command.col] += Number(command.alpha) || 0;
        temp.cells.commandZ[command.row][command.col] += Number(command.z) || 0;
      }
    }
    return temp;
  }

  function calibrationExperimentResultsTemplate(state, options = {}) {
    const protocol = options.protocol || calibrationExperimentProtocol(state, options);
    const steps = [];
    for (const step of protocol.steps || []) {
      const repeats = Math.max(1, Math.round(Number(step.repeatCount) || 1));
      for (let repeatIndex = 1; repeatIndex <= repeats; repeatIndex += 1) {
        steps.push({
          stepId: step.id,
          repeatIndex,
          cells: (step.observationCells || []).map((cell) => ({
            row: cell.row,
            col: cell.col,
            alphaDelta: null,
            heightDelta: null,
            centerDelta: null,
            pinHoleSlipMm: null,
            actuatorForceN: null,
          })),
          notes: "Replace null fields with measured bench data.",
        });
      }
    }
    return {
      schema: "rad-sim.calibration-experiment-results.v1",
      protocolSchema: protocol.schema,
      hardwareProfile: protocol.hardwareProfile,
      notes: "Fill optional measured fields with real bench measurements.",
      steps,
    };
  }

  function exportCalibrationExperimentResultsTemplate(state, options = {}) {
    return JSON.stringify(calibrationExperimentResultsTemplate(state, options), null, 2);
  }

  function numericOrNull(value) {
    if (value === null || value === undefined || value === "") return null;
    const numeric = Number(value);
    return Number.isFinite(numeric) ? numeric : null;
  }

  function rmse(values) {
    if (!values.length) return null;
    const mean = values.reduce((sum, value) => sum + value * value, 0) / values.length;
    return Math.sqrt(mean);
  }

  function meanOrNull(values) {
    if (!values.length) return null;
    return values.reduce((sum, value) => sum + value, 0) / values.length;
  }

  function calibrationLinearFit(pairs) {
    const sampleCount = pairs.length;
    if (!sampleCount) {
      return {
        sampleCount: 0,
        gain: null,
        bias: null,
        suggestedGain: null,
        suggestedBias: null,
        gainIdentifiable: false,
        rmsRawError: null,
        rmsResidual: null,
        meanPredicted: null,
        meanMeasured: null,
        predictedRange: null,
        measuredRange: null,
      };
    }
    const predicted = pairs.map((pair) => pair.predicted);
    const measured = pairs.map((pair) => pair.measured);
    const meanPredicted = meanOrNull(predicted);
    const meanMeasured = meanOrNull(measured);
    const variance = predicted.reduce((sum, value) => sum + (value - meanPredicted) ** 2, 0);
    const covariance = pairs.reduce(
      (sum, pair) => sum + (pair.predicted - meanPredicted) * (pair.measured - meanMeasured),
      0
    );
    const gainIdentifiable = variance > 1e-12;
    const gain = gainIdentifiable ? covariance / variance : null;
    const suggestedGain = gainIdentifiable ? gain : 1;
    const bias = gainIdentifiable
      ? meanMeasured - gain * meanPredicted
      : meanOrNull(pairs.map((pair) => pair.measured - pair.predicted));
    const residuals = pairs.map((pair) => pair.measured - (suggestedGain * pair.predicted + bias));
    const rawErrors = pairs.map((pair) => pair.measured - pair.predicted);
    return {
      sampleCount,
      gain,
      bias,
      suggestedGain,
      suggestedBias: bias,
      gainIdentifiable,
      rmsRawError: rmse(rawErrors),
      rmsResidual: rmse(residuals),
      meanPredicted,
      meanMeasured,
      predictedRange: [Math.min(...predicted), Math.max(...predicted)],
      measuredRange: [Math.min(...measured), Math.max(...measured)],
    };
  }

  function createCalibrationErrorField(state) {
    const rows = state.grid.rows;
    const cols = state.grid.cols;
    return {
      alphaError: RAD.matrix(rows, cols, 0),
      heightError: RAD.matrix(rows, cols, 0),
      combinedError: RAD.matrix(rows, cols, 0),
      sampleCount: RAD.matrix(rows, cols, 0),
      alphaSampleCount: RAD.matrix(rows, cols, 0),
      heightSampleCount: RAD.matrix(rows, cols, 0),
      maxAbsAlphaError: 0,
      maxAbsHeightError: 0,
      maxCombinedError: 0,
      worstCell: null,
    };
  }

  function finalizeCalibrationErrorField(field) {
    for (let r = 0; r < field.sampleCount.length; r += 1) {
      for (let c = 0; c < field.sampleCount[r].length; c += 1) {
        const alphaCount = field.alphaSampleCount[r][c];
        const heightCount = field.heightSampleCount[r][c];
        const alpha = alphaCount ? field.alphaError[r][c] / alphaCount : 0;
        const height = heightCount ? field.heightError[r][c] / heightCount : 0;
        const combined = Math.hypot(alpha, height);
        field.alphaError[r][c] = alpha;
        field.heightError[r][c] = height;
        field.combinedError[r][c] = combined;
        field.maxAbsAlphaError = Math.max(field.maxAbsAlphaError, Math.abs(alpha));
        field.maxAbsHeightError = Math.max(field.maxAbsHeightError, Math.abs(height));
        field.maxCombinedError = Math.max(field.maxCombinedError, combined);
        if (
          field.sampleCount[r][c] > 0 &&
          (!field.worstCell || combined > field.worstCell.combinedError)
        ) {
          field.worstCell = {
            row: r,
            col: c,
            alphaError: alpha,
            heightError: height,
            combinedError: combined,
            sampleCount: field.sampleCount[r][c],
            alphaSampleCount: alphaCount,
            heightSampleCount: heightCount,
          };
        }
      }
    }
    return field;
  }

  function compareCalibrationExperimentResults(state, results, options = {}) {
    const parsed = typeof results === "string" ? JSON.parse(results) : results;
    if (parsed?.schema !== "rad-sim.calibration-experiment-results.v1") {
      throw new Error("unsupported calibration experiment results schema");
    }
    const protocol = options.protocol || calibrationExperimentProtocol(state, options);
    const stepsById = new Map((protocol.steps || []).map((step) => [step.id, step]));
    const comparisons = [];
    const field = createCalibrationErrorField(state);
    const alphaPairs = [];
    const heightPairs = [];
    for (const measuredStep of parsed.steps || []) {
      const step = stepsById.get(measuredStep.stepId);
      if (!step) continue;
      const commandState = protocolStepState(state, step, true);
      const baselineState = protocolStepState(state, step, false);
      const sim = RAD.simulate(commandState);
      const baseline = RAD.simulate(baselineState);
      const alphaErrors = [];
      const heightErrors = [];
      const centerErrors = [];
      const slipValues = [];
      const forceValues = [];
      const observed = new Set();
      for (const cell of measuredStep.cells || []) {
        const row = Number(cell.row);
        const col = Number(cell.col);
        if (!Number.isInteger(row) || !Number.isInteger(col) || row < 0 || col < 0 || row >= state.grid.rows || col >= state.grid.cols) continue;
        observed.add(`${row},${col}`);
        const alpha = numericOrNull(cell.alphaDelta);
        const height = numericOrNull(cell.heightDelta);
        let measuredField = false;
        if (alpha !== null) {
          const predictedAlpha = (sim.alpha?.[row]?.[col] || 0) - (baseline.alpha?.[row]?.[col] || 0);
          const alphaError = alpha - predictedAlpha;
          alphaErrors.push(alphaError);
          alphaPairs.push({ predicted: predictedAlpha, measured: alpha });
          field.alphaError[row][col] += alphaError;
          field.alphaSampleCount[row][col] += 1;
          measuredField = true;
        }
        if (height !== null) {
          const predictedHeight = (sim.height?.[row]?.[col] || 0) - (baseline.height?.[row]?.[col] || 0);
          const heightError = height - predictedHeight;
          heightErrors.push(heightError);
          heightPairs.push({ predicted: predictedHeight, measured: height });
          field.heightError[row][col] += heightError;
          field.heightSampleCount[row][col] += 1;
          measuredField = true;
        }
        if (measuredField) field.sampleCount[row][col] += 1;
        if (Array.isArray(cell.centerDelta) && cell.centerDelta.length === 3) {
          const current = sim.centers?.[row]?.[col] || { x: 0, y: 0, z: 0 };
          const base = baseline.centers?.[row]?.[col] || { x: 0, y: 0, z: 0 };
          centerErrors.push(
            Math.hypot(
              Number(cell.centerDelta[0]) - (current.x - base.x),
              Number(cell.centerDelta[1]) - (current.y - base.y),
              Number(cell.centerDelta[2]) - (current.z - base.z)
            )
          );
        }
        const slip = numericOrNull(cell.pinHoleSlipMm);
        const force = numericOrNull(cell.actuatorForceN);
        if (slip !== null) slipValues.push(slip);
        if (force !== null) forceValues.push(force);
      }
      const missingObservationCount = (step.observationCells || []).filter((cell) => !observed.has(`${cell.row},${cell.col}`)).length;
      comparisons.push({
        stepId: measuredStep.stepId,
        repeatIndex: Number(measuredStep.repeatIndex) || 1,
        measuredCellCount: observed.size,
        missingObservationCount,
        alphaRmse: rmse(alphaErrors),
        heightRmse: rmse(heightErrors),
        centerRmse: rmse(centerErrors),
        meanSignedAlphaError: meanOrNull(alphaErrors),
        meanSignedHeightError: meanOrNull(heightErrors),
        meanAbsAlphaError: meanOrNull(alphaErrors.map(Math.abs)),
        meanAbsHeightError: meanOrNull(heightErrors.map(Math.abs)),
        maxAbsHeightError: heightErrors.length ? Math.max(...heightErrors.map(Math.abs)) : null,
        meanActuatorForceN: meanOrNull(forceValues),
        meanPinHoleSlipMm: meanOrNull(slipValues),
      });
    }
    return {
      schema: "rad-sim.calibration-experiment-comparison.v1",
      comparisons,
      field: finalizeCalibrationErrorField(field),
      fit: {
        alpha: calibrationLinearFit(alphaPairs),
        height: calibrationLinearFit(heightPairs),
      },
    };
  }

  function finiteAverage(values) {
    const finite = values.filter((value) => Number.isFinite(value));
    if (!finite.length) return null;
    return finite.reduce((sum, value) => sum + value, 0) / finite.length;
  }

  function summarizeCalibrationComparison(comparison) {
    const comparisons = comparison?.comparisons || [];
    const heightErrors = comparisons.map((item) => item.heightRmse).filter(Number.isFinite);
    const alphaErrors = comparisons.map((item) => item.alphaRmse).filter(Number.isFinite);
    const centerErrors = comparisons.map((item) => item.centerRmse).filter(Number.isFinite);
    const signedAlphaErrors = comparisons.map((item) => item.meanSignedAlphaError).filter(Number.isFinite);
    const signedHeightErrors = comparisons.map((item) => item.meanSignedHeightError).filter(Number.isFinite);
    const absAlphaErrors = comparisons.map((item) => item.meanAbsAlphaError).filter(Number.isFinite);
    const absHeightErrors = comparisons.map((item) => item.meanAbsHeightError).filter(Number.isFinite);
    const forceValues = comparisons.map((item) => item.meanActuatorForceN).filter(Number.isFinite);
    const slipValues = comparisons.map((item) => item.meanPinHoleSlipMm).filter(Number.isFinite);
    const field = comparison?.field || {};
    let worst = null;
    for (const item of comparisons) {
      const value = Number(item.maxAbsHeightError);
      if (!Number.isFinite(value)) continue;
      if (!worst || Math.abs(value) > Math.abs(worst.maxAbsHeightError)) worst = item;
    }
    return {
      schema: "rad-sim.calibration-experiment-comparison-summary.v1",
      stepCount: comparisons.length,
      measuredCellCount: comparisons.reduce((sum, item) => sum + (Number(item.measuredCellCount) || 0), 0),
      missingObservationCount: comparisons.reduce((sum, item) => sum + (Number(item.missingObservationCount) || 0), 0),
      alphaRmseMean: finiteAverage(alphaErrors),
      heightRmseMean: finiteAverage(heightErrors),
      centerRmseMean: finiteAverage(centerErrors),
      meanSignedAlphaError: finiteAverage(signedAlphaErrors),
      meanSignedHeightError: finiteAverage(signedHeightErrors),
      meanAbsAlphaError: finiteAverage(absAlphaErrors),
      meanAbsHeightError: finiteAverage(absHeightErrors),
      fit: comparison?.fit || null,
      maxAbsHeightError: worst ? Number(worst.maxAbsHeightError) : null,
      maxAbsAlphaError: Number.isFinite(field.maxAbsAlphaError) ? field.maxAbsAlphaError : null,
      maxCombinedError: Number.isFinite(field.maxCombinedError) ? field.maxCombinedError : null,
      worstCell: field.worstCell || null,
      worstStepId: worst?.stepId || null,
      meanActuatorForceN: finiteAverage(forceValues),
      meanPinHoleSlipMm: finiteAverage(slipValues),
    };
  }

  function calibrationComparisonReport(state) {
    const comparison = state.experiment?.calibrationComparison;
    if (!comparison) throw new Error("no imported calibration comparison is available");
    const summary = state.experiment?.calibrationComparisonSummary || summarizeCalibrationComparison(comparison);
    return {
      schema: "rad-sim.calibration-comparison-report.v1",
      savedAt: new Date().toISOString(),
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        cellSize: state.grid.cellSize,
        backlash: state.grid.backlash,
        couplingGain: state.grid.couplingGain,
        zCouplingGain: state.grid.zCouplingGain,
        pinRadius: state.grid.pinRadius,
        holeRadius: state.grid.holeRadius,
      },
      hardwareProfile: state.grid.hardwareProfile?.name || null,
      sourceResultsSchema: state.experiment?.calibrationResults?.schema || null,
      summary,
      comparison,
    };
  }

  function exportCalibrationComparisonReport(state) {
    return JSON.stringify(calibrationComparisonReport(state), null, 2);
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
  RAD.calibrationExperimentProtocol = calibrationExperimentProtocol;
  RAD.exportCalibrationExperimentProtocol = exportCalibrationExperimentProtocol;
  RAD.calibrationExperimentResultsTemplate = calibrationExperimentResultsTemplate;
  RAD.exportCalibrationExperimentResultsTemplate = exportCalibrationExperimentResultsTemplate;
  RAD.compareCalibrationExperimentResults = compareCalibrationExperimentResults;
  RAD.summarizeCalibrationComparison = summarizeCalibrationComparison;
  RAD.calibrationComparisonReport = calibrationComparisonReport;
  RAD.exportCalibrationComparisonReport = exportCalibrationComparisonReport;
  RAD.physicalPreviewComparison = physicalPreviewComparison;
  RAD.responseDecayProfile = responseDecayProfile;
})();
