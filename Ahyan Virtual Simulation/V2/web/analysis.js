(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});
  const SWEEP_SENSITIVITY_METRICS = Object.freeze([
    "meanAlphaReach",
    "meanZReach",
    "maxObservedNeighborZResidual",
    "maxSuperpositionError",
  ]);

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
    let positiveZReachCells = 0;
    let negativeZReachCells = 0;
    let maxAlphaDelta = 0;
    let maxHeightDelta = 0;
    let maxPositiveHeightDelta = 0;
    let maxNegativeHeightDelta = 0;
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
        if (heightDelta > 1e-8) positiveZReachCells += 1;
        if (heightDelta < -1e-8) negativeZReachCells += 1;
        maxAlphaDelta = Math.max(maxAlphaDelta, alphaAbs);
        maxHeightDelta = Math.max(maxHeightDelta, heightAbs);
        maxPositiveHeightDelta = Math.max(maxPositiveHeightDelta, heightDelta);
        maxNegativeHeightDelta = Math.min(maxNegativeHeightDelta, heightDelta);
        meanAbsAlphaDelta += alphaAbs;
        meanAbsHeightDelta += heightAbs;
      }
    }
    return {
      activeSources,
      responseCells,
      alphaReachCells,
      zReachCells,
      positiveZReachCells,
      negativeZReachCells,
      maxAlphaDelta,
      maxHeightDelta,
      maxPositiveHeightDelta,
      maxNegativeHeightDelta,
      meanAbsAlphaDelta: meanAbsAlphaDelta / totalCells,
      meanAbsHeightDelta: meanAbsHeightDelta / totalCells,
      alphaDieOff: finiteMax(sim.dieOff),
      zDieOff: finiteMax(sim.zDieOff),
    };
  }

  function superpositionError(state, sourceCells, combinedSim, baselineSim) {
    if (sourceCells.length <= 1) return { rms: 0, max: 0, alphaMax: 0, heightMax: 0, skipped: false, sourceCount: sourceCells.length };
    const activeSources = sourceCells.filter((cell) => Math.abs(state.cells.commandAlpha?.[cell.r]?.[cell.c] || 0) > 1e-9 || Math.abs(state.cells.commandZ?.[cell.r]?.[cell.c] || 0) > 1e-9);
    if (activeSources.length <= 1) return { rms: 0, max: 0, alphaMax: 0, heightMax: 0, skipped: false, sourceCount: activeSources.length };
    if (activeSources.length > 16) return { rms: null, max: null, alphaMax: null, heightMax: null, skipped: true, sourceCount: activeSources.length };
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
    let alphaMax = 0;
    let heightMax = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        const alphaResidual = (combinedSim.alpha[r][c] || 0) - (baselineSim.alpha[r][c] || 0) - alphaSum[r][c];
        const heightResidual = (combinedSim.height[r][c] || 0) - (baselineSim.height[r][c] || 0) - heightSum[r][c];
        squared += alphaResidual * alphaResidual + heightResidual * heightResidual;
        alphaMax = Math.max(alphaMax, Math.abs(alphaResidual));
        heightMax = Math.max(heightMax, Math.abs(heightResidual));
        max = Math.max(max, alphaMax, heightMax);
      }
    }
    return { rms: Math.sqrt(squared / Math.max(1, rows * cols * 2)), max, alphaMax, heightMax, skipped: false, sourceCount: activeSources.length };
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

  function matrixFromColumns(columns, rowCount) {
    return Array.from({ length: rowCount }, (_, row) => columns.map((column) => Number(column[row]) || 0));
  }

  function responseMatrixCellOrder(rows, cols) {
    const cells = [];
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) cells.push({ index: r * cols + c, row: r, col: c });
    }
    return cells;
  }

  function normalizeMatrixCell(state, cell) {
    return {
      r: clampIndex(cell?.r ?? cell?.row ?? 0, state.grid.rows),
      c: clampIndex(cell?.c ?? cell?.col ?? 0, state.grid.cols),
    };
  }

  function buildResponseMatrix(state, options = {}) {
    const tolerance = Number(options.tolerance ?? 1e-9);
    const includeAlpha = options.includeAlpha !== false;
    const includeZ = options.includeZ !== false;
    if (!includeAlpha && !includeZ) throw new Error("at least one command family must be included");
    const scope = options.scope || state.experiment?.characterizationScope || "single";
    const selected = {
      r: clampIndex(options.r ?? state.selection?.r ?? 0, state.grid.rows),
      c: clampIndex(options.c ?? state.selection?.c ?? 0, state.grid.cols),
    };
    const rawCells = options.actuatorCells || characterizationCells(state, selected.r, selected.c, scope);
    const sourceCells = uniqueProtocolCells(rawCells.map((cell) => normalizeMatrixCell(state, cell)));
    const limits = RAD.commandLimits(state);
    const stepAlpha = Math.min(Math.abs(Number(options.alphaStep ?? 0.12)), limits.alphaContract);
    const stepZ = Math.min(Math.abs(Number(options.zStep ?? 0.12)), limits.z);
    const baselineState = scopedState(state, []);
    const baselineSim = RAD.simulate(baselineState);
    const alphaColumns = [];
    const heightColumns = [];
    const commands = [];
    const pushColumn = (cell, family, commandAlpha, commandZ) => {
      const nextState = scopedState(state, []);
      nextState.cells.commandAlpha[cell.r][cell.c] = commandAlpha;
      nextState.cells.commandZ[cell.r][cell.c] = commandZ;
      const nextSim = RAD.simulate(nextState);
      alphaColumns.push(responseVector(nextSim, baselineSim, "alpha"));
      heightColumns.push(responseVector(nextSim, baselineSim, "height"));
      commands.push({
        index: commands.length,
        row: cell.r,
        col: cell.c,
        family,
        alpha: commandAlpha,
        z: commandZ,
        locked: Boolean(state.cells.locked?.[cell.r]?.[cell.c]),
      });
    };
    for (const cell of sourceCells) {
      if (includeAlpha && stepAlpha > tolerance) pushColumn(cell, "alpha", -stepAlpha, 0);
      if (includeZ && stepZ > tolerance) pushColumn(cell, "z", 0, stepZ);
    }
    const totalCells = state.grid.rows * state.grid.cols;
    const reachableAlphaCells = countReachableFromColumns(alphaColumns, tolerance);
    const reachableHeightCells = countReachableFromColumns(heightColumns, tolerance);
    return {
      schema: "rad-sim.response-matrix.v1",
      grid: { rows: state.grid.rows, cols: state.grid.cols },
      source: {
        scope,
        selected,
        actuatorCellCount: sourceCells.length,
        alphaStep: -stepAlpha,
        zStep: stepZ,
      },
      cellOrder: responseMatrixCellOrder(state.grid.rows, state.grid.cols),
      commands,
      alpha: matrixFromColumns(alphaColumns, totalCells),
      height: matrixFromColumns(heightColumns, totalCells),
      diagnostics: {
        columnCount: commands.length,
        alphaRank: matrixRankFromColumns(alphaColumns, tolerance),
        heightRank: matrixRankFromColumns(heightColumns, tolerance),
        reachableAlphaCells,
        reachableHeightCells,
        alphaUnderactuatedCells: Math.max(0, totalCells - reachableAlphaCells),
        heightUnderactuatedCells: Math.max(0, totalCells - reachableHeightCells),
        tolerance,
      },
    };
  }

  function exportResponseMatrix(state, options = {}) {
    return JSON.stringify(buildResponseMatrix(state, options), null, 2);
  }

  function topologyScenarioCell(cell, state) {
    return {
      r: clampIndex(cell?.r ?? cell?.row ?? 0, state.grid.rows),
      c: clampIndex(cell?.c ?? cell?.col ?? 0, state.grid.cols),
    };
  }

  function topologyScenarioState(state, scenario) {
    const next = cloneForCharacterization(state);
    RAD.clearCommands(next);
    for (let r = 0; r < next.grid.rows; r += 1) {
      for (let c = 0; c < next.grid.cols; c += 1) {
        next.cells.locked[r][c] = false;
        if (next.cells.removed) next.cells.removed[r][c] = false;
      }
    }
    for (const cell of scenario.lockedCells || []) {
      const { r, c } = topologyScenarioCell(cell, next);
      next.cells.locked[r][c] = true;
    }
    for (const cell of scenario.removedCells || []) {
      const { r, c } = topologyScenarioCell(cell, next);
      if (next.cells.removed) next.cells.removed[r][c] = true;
      next.cells.commandAlpha[r][c] = 0;
      next.cells.commandZ[r][c] = 0;
    }
    for (const command of scenario.commands || []) {
      const { r, c } = topologyScenarioCell(command, next);
      if (next.cells.removed?.[r]?.[c]) continue;
      next.cells.commandAlpha[r][c] += Number(command.alpha ?? command.commandAlpha) || 0;
      next.cells.commandZ[r][c] += Number(command.z ?? command.commandZ) || 0;
    }
    return next;
  }

  function topologyScenarioCandidates(state, scenario) {
    const cells = scenario.actuatorCells || scenario.commands || [];
    const rawCells = cells.length
      ? cells
      : Array.from({ length: state.grid.rows * state.grid.cols }, (_, index) => ({
          r: Math.floor(index / state.grid.cols),
          c: index % state.grid.cols,
        }));
    const seen = new Set();
    const candidates = [];
    for (const cell of rawCells) {
      const { r, c } = topologyScenarioCell(cell, state);
      const key = `${r},${c}`;
      if (seen.has(key)) continue;
      seen.add(key);
      if (state.cells.removed?.[r]?.[c] || state.cells.locked?.[r]?.[c]) continue;
      if (state.cells.actuatorAllowed?.[r]?.[c] === false) continue;
      candidates.push({ r, c });
    }
    return candidates;
  }

  function topologyComponentReach(topology, candidates) {
    const active = new Set();
    for (const cell of candidates) {
      const label = topology.componentLabels?.[cell.r]?.[cell.c];
      if (Number.isInteger(label) && label >= 0) active.add(label);
    }
    let reachable = 0;
    let blocked = 0;
    for (const row of topology.componentLabels || []) {
      for (const label of row || []) {
        if (!Number.isInteger(label) || label < 0) continue;
        if (active.has(label)) reachable += 1;
        else blocked += 1;
      }
    }
    return { reachable, blocked };
  }

  function componentResponseSummaries(topology, responseMatrix, tolerance) {
    const labels = topology.componentLabels || [];
    const rows = labels.length;
    const cols = labels[0]?.length || 0;
    const summaries = [];
    for (let label = 0; label < (topology.componentCount || 0); label += 1) {
      const rowIndices = [];
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          if (labels[r][c] === label) rowIndices.push(r * cols + c);
        }
      }
      const commandIndices = [];
      const actuatorCells = new Set();
      for (const command of responseMatrix.commands || []) {
        const commandLabel = labels?.[command.row]?.[command.col];
        if (commandLabel !== label) continue;
        commandIndices.push(command.index);
        actuatorCells.add(`${command.row},${command.col}`);
      }
      const alphaColumns = commandIndices.map((index) => rowIndices.map((row) => Number(responseMatrix.alpha?.[row]?.[index]) || 0));
      const heightColumns = commandIndices.map((index) => rowIndices.map((row) => Number(responseMatrix.height?.[row]?.[index]) || 0));
      const reachableAlpha = countReachableFromColumns(alphaColumns, tolerance);
      const reachableHeight = countReachableFromColumns(heightColumns, tolerance);
      summaries.push({
        label,
        cellCount: rowIndices.length,
        actuatorCellCount: actuatorCells.size,
        commandColumnCount: commandIndices.length,
        alphaRank: matrixRankFromColumns(alphaColumns, tolerance),
        heightRank: matrixRankFromColumns(heightColumns, tolerance),
        reachableAlphaCells: reachableAlpha,
        reachableHeightCells: reachableHeight,
        blockedAlphaCells: Math.max(0, rowIndices.length - reachableAlpha),
        blockedHeightCells: Math.max(0, rowIndices.length - reachableHeight),
      });
    }
    return summaries;
  }

  function topologyScenarioComparisons(rows) {
    if (!rows.length) return [];
    const baseline = rows[0];
    return rows.slice(1).map((row) => ({
      baseline: baseline.name,
      scenario: row.name,
      alphaRankDelta: (row.responseMatrix.alphaRank || 0) - (baseline.responseMatrix.alphaRank || 0),
      heightRankDelta: (row.responseMatrix.heightRank || 0) - (baseline.responseMatrix.heightRank || 0),
      reachableAlphaCellDelta: (row.responseMatrix.reachableAlphaCells || 0) - (baseline.responseMatrix.reachableAlphaCells || 0),
      reachableHeightCellDelta: (row.responseMatrix.reachableHeightCells || 0) - (baseline.responseMatrix.reachableHeightCells || 0),
      componentCountDelta: (row.topology.componentCount || 0) - (baseline.topology.componentCount || 0),
      componentBlockedCellDelta: (row.metrics.componentBlockedCells || 0) - (baseline.metrics.componentBlockedCells || 0),
      deletedEdgeDelta: (row.topology.deletedEdges || 0) - (baseline.topology.deletedEdges || 0),
    }));
  }

  function signedHeightReach(state, candidates, stepZ, tolerance) {
    const positive = RAD.matrix(state.grid.rows, state.grid.cols, false);
    const negative = RAD.matrix(state.grid.rows, state.grid.cols, false);
    const step = Math.abs(Number(stepZ) || 0);
    if (step <= tolerance) return { positive: 0, negative: 0 };
    const baseline = scopedState(state, []);
    const baselineSim = RAD.simulate(baseline);
    for (const cell of candidates) {
      for (const signedStep of [step, -step]) {
        const next = scopedState(state, []);
        next.cells.commandZ[cell.r][cell.c] = signedStep;
        const sim = RAD.simulate(next);
        for (let r = 0; r < state.grid.rows; r += 1) {
          for (let c = 0; c < state.grid.cols; c += 1) {
            const delta = (sim.height?.[r]?.[c] || 0) - (baselineSim.height?.[r]?.[c] || 0);
            if (delta > tolerance) positive[r][c] = true;
            if (delta < -tolerance) negative[r][c] = true;
          }
        }
      }
    }
    const count = (matrix) => matrix.reduce((sum, row) => sum + row.filter(Boolean).length, 0);
    return { positive: count(positive), negative: count(negative) };
  }

  function maxAbsMatrix(matrix) {
    let value = 0;
    for (const row of matrix || []) {
      for (const entry of row || []) value = Math.max(value, Math.abs(Number(entry) || 0));
    }
    return value;
  }

  function topologyExperimentReport(state, scenarios = null, options = {}) {
    const tolerance = Number(options.tolerance ?? 1e-9);
    const limits = RAD.commandLimits(state);
    const stepZ = Math.min(Math.abs(Number(options.zStep ?? 0.12)), limits.z);
    const center = {
      r: clampIndex(options.center?.r ?? Math.floor(state.grid.rows / 2), state.grid.rows),
      c: clampIndex(options.center?.c ?? Math.floor(state.grid.cols / 2), state.grid.cols),
    };
    const neighbors = RAD.brushCells(state, center.r, center.c, 1).filter((cell) => cell.r !== center.r || cell.c !== center.c);
    const group = [center, ...neighbors];
    const defaultScenarios = [
      {
        name: "intact-group-actuation",
        commands: group.map((cell) => ({ ...cell, alpha: -0.25, z: 0.2 })),
        actuatorCells: group,
      },
      {
        name: "locked-center-group-actuation",
        commands: group.map((cell) => ({ ...cell, alpha: -0.25, z: 0.2 })),
        lockedCells: [center],
        actuatorCells: group,
      },
      {
        name: "removed-center-group-actuation",
        commands: group.map((cell) => ({ ...cell, alpha: -0.25, z: 0.2 })),
        removedCells: [center],
        actuatorCells: group,
      },
    ];
    const rows = [];
    for (const scenario of scenarios || defaultScenarios) {
      const scenarioState = topologyScenarioState(state, scenario);
      const sim = RAD.simulate(scenarioState);
      const topology = RAD.topologyDiagnostics(scenarioState);
      const candidates = topologyScenarioCandidates(scenarioState, scenario);
      const responseMatrix = buildResponseMatrix(scenarioState, {
        actuatorCells: candidates,
        scope: "lattice",
        tolerance,
        zStep: options.zStep,
        alphaStep: options.alphaStep,
      });
      const componentReach = topologyComponentReach(topology, candidates);
      const componentResponseRank = componentResponseSummaries(topology, responseMatrix, tolerance);
      const signedReach = signedHeightReach(scenarioState, candidates, stepZ, tolerance);
      rows.push({
        name: scenario.name || "scenario",
        commands: (scenario.commands || []).map((command) => ({
          row: topologyScenarioCell(command, scenarioState).r,
          col: topologyScenarioCell(command, scenarioState).c,
          alpha: Number(command.alpha ?? command.commandAlpha) || 0,
          z: Number(command.z ?? command.commandZ) || 0,
        })),
        lockedCells: (scenario.lockedCells || []).map((cell) => {
          const normalized = topologyScenarioCell(cell, scenarioState);
          return { row: normalized.r, col: normalized.c };
        }),
        removedCells: (scenario.removedCells || []).map((cell) => {
          const normalized = topologyScenarioCell(cell, scenarioState);
          return { row: normalized.r, col: normalized.c };
        }),
        actuatorCells: candidates.map((cell) => ({ row: cell.r, col: cell.c })),
        topology,
        responseMatrix: responseMatrix.diagnostics,
        metrics: {
          componentReachableCells: componentReach.reachable,
          componentBlockedCells: componentReach.blocked,
          componentResponseRank,
          positiveHeightReachableCells: signedReach.positive,
          negativeHeightReachableCells: signedReach.negative,
          targetHeightRms: sim.metrics.rmsTargetError,
          maxAbsTargetHeightResidual: maxAbsMatrix(sim.targetError),
          tolerance,
        },
      });
    }
    return {
      schema: "rad-sim.browser-topology-experiment-report.v1",
      grid: { rows: state.grid.rows, cols: state.grid.cols },
      summary: {
        scenarioCount: rows.length,
        maxComponentCount: Math.max(0, ...rows.map((row) => row.topology.componentCount || 0)),
        maxDeletedEdges: Math.max(0, ...rows.map((row) => row.topology.deletedEdges || 0)),
        maxComponentBlockedCells: Math.max(0, ...rows.map((row) => row.metrics.componentBlockedCells || 0)),
        minHeightReachableCells: Math.min(...rows.map((row) => row.responseMatrix.reachableHeightCells || 0)),
        maxComponentHeightRank: Math.max(
          0,
          ...rows.flatMap((row) => row.metrics.componentResponseRank.map((component) => component.heightRank || 0))
        ),
        maxComponentBlockedHeightCells: Math.max(
          0,
          ...rows.flatMap((row) => row.metrics.componentResponseRank.map((component) => component.blockedHeightCells || 0))
        ),
      },
      scenarioComparisons: topologyScenarioComparisons(rows),
      scenarios: rows,
      notes: "Browser topology experiment report is a finite-response simulator diagnostic, not calibrated hardware validation.",
    };
  }

  function exportTopologyExperimentReport(state, scenarios = null, options = {}) {
    return JSON.stringify(topologyExperimentReport(state, scenarios, options), null, 2);
  }

  function reportCell(cell) {
    return { row: Number(cell.r ?? cell.row) || 0, col: Number(cell.c ?? cell.col) || 0 };
  }

  function activeReportCommands(state, sourceCells, tolerance = 1e-9) {
    return sourceCells.map((cell, index) => ({
      index,
      row: cell.r,
      col: cell.c,
      alpha: Number(state.cells.commandAlpha?.[cell.r]?.[cell.c]) || 0,
      z: Number(state.cells.commandZ?.[cell.r]?.[cell.c]) || 0,
      locked: Boolean(state.cells.locked?.[cell.r]?.[cell.c]),
      active:
        Math.abs(Number(state.cells.commandAlpha?.[cell.r]?.[cell.c]) || 0) > tolerance ||
        Math.abs(Number(state.cells.commandZ?.[cell.r]?.[cell.c]) || 0) > tolerance,
    }));
  }

  function lockedReportCells(state) {
    const cells = [];
    for (let r = 0; r < state.grid.rows; r += 1) {
      for (let c = 0; c < state.grid.cols; c += 1) {
        if (state.cells.locked?.[r]?.[c]) cells.push({ row: r, col: c });
      }
    }
    return cells;
  }

  function finiteMatrix(field) {
    return (field || []).map((row) =>
      (row || []).map((value) => (Number.isFinite(Number(value)) ? Number(value) : null))
    );
  }

  function deltaMatrix(next, base, field) {
    return (next?.[field] || []).map((row, r) =>
      (row || []).map((value, c) => (Number(value) || 0) - (Number(base?.[field]?.[r]?.[c]) || 0))
    );
  }

  function reportResponseFields(sim, baselineSim) {
    return {
      alphaDelta: deltaMatrix(sim, baselineSim, "alpha"),
      heightDelta: deltaMatrix(sim, baselineSim, "height"),
      actuatorInfluence: finiteMatrix(sim.influence),
      zResidual: finiteMatrix(sim.zResidual),
      alphaDieOff: finiteMatrix(sim.dieOff),
      zDieOff: finiteMatrix(sim.zDieOff),
    };
  }

  function reportPhysicalPreview(characterization) {
    return {
      schema: "rad-sim.browser-physical-preview.v1",
      model: "browser-spring-preview",
      physicalPreviewAvailable: Boolean(characterization.physicalPreviewAvailable),
      physicalSuccess: Boolean(characterization.physicalPreviewSuccess),
      heightRmsModelError: Number(characterization.physicalHeightRmsError) || 0,
      heightMaxModelError: Number(characterization.physicalHeightMaxError) || 0,
      centerRmsModelError: Number(characterization.physicalCenterRmsError) || 0,
      centerMaxModelError: Number(characterization.physicalCenterMaxError) || 0,
      iterations: Number(characterization.physicalPreviewIterations) || 0,
      error: characterization.physicalPreviewError || null,
      note: "Browser spring-preview validation is an interactive approximation; use Python spring_hinge_3d reports for the reference research artifact.",
    };
  }

  function reportSequenceOrder(state, commands, tolerance = 1e-9) {
    if (typeof RAD.compareSequenceOrder !== "function" || typeof RAD.localActuationEvent !== "function") {
      return { eventSequence: [], sequenceOrder: null };
    }
    const primary = commands.find((command) => command.active);
    if (!primary) {
      return { eventSequence: [], sequenceOrder: null };
    }
    const cell = { r: primary.row, c: primary.col };
    const alpha = primary.alpha;
    const z = primary.z;
    const eventSequence = [
      RAD.localActuationEvent(cell, alpha, z),
      RAD.lockEvent(cell),
      RAD.clearActuationEvent(cell),
    ];
    const sequenceOrder = RAD.compareSequenceOrder(scopedState(state, []), eventSequence, tolerance);
    return {
      eventSequence: eventSequence.map((event, index) => ({
        index,
        kind: event.kind,
        cell: event.cell ? reportCell(event.cell) : null,
        alpha: Number(event.alpha) || 0,
        z: Number(event.z) || 0,
      })),
      sequenceOrder,
    };
  }

  function reportLawCandidate(id, operatorClass, property, statement, supported, evidence) {
    return {
      id,
      operatorClass,
      property,
      statement,
      supportedByDiagnostic: Boolean(supported),
      status: "simulator-diagnostic",
      evidence,
    };
  }

  function frameworkOperatorLawCandidates(characterization, sequenceOrder, eventCount, totalCells, tolerance) {
    const maxSuperpositionError = Math.max(
      Number(characterization.alphaSuperpositionError) || 0,
      Number(characterization.heightSuperpositionError) || 0,
      Number(characterization.superpositionError) || 0
    );
    const localitySupported =
      (Number(characterization.alphaDieOff) || 0) < totalCells &&
      (Number(characterization.zDieOff) || 0) < totalCells;
    const rankLimited =
      (Number(characterization.alphaUnderactuatedCells) || 0) > 0 ||
      (Number(characterization.heightUnderactuatedCells) || 0) > 0;
    const laws = [
      reportLawCandidate(
        "bounded_locality",
        "actuation and clearance operators",
        "locality",
        "The measured response is bounded by finite alpha and height die-off radii at the report tolerance.",
        localitySupported,
        {
          alphaLocalityRadius: Number(characterization.alphaDieOff) || 0,
          zLocalityRadius: Number(characterization.zDieOff) || 0,
          alphaDecayRatio: Number(characterization.alphaDecayRatio) || 0,
          zDecayRatio: Number(characterization.zDecayRatio) || 0,
          tolerance,
        }
      ),
      reportLawCandidate(
        "rank_limited_reachability",
        "finite actuator set",
        "reachable set",
        "The selected actuator operators span a finite response subspace; unreached cells mark underactuated regions for the current command basis.",
        rankLimited,
        {
          alphaRank: Number(characterization.responseRankAlpha) || 0,
          heightRank: Number(characterization.responseRankHeight) || 0,
          reachableAlphaCells: Number(characterization.reachableAlphaCells) || 0,
          reachableHeightCells: Number(characterization.reachableHeightCells) || 0,
          alphaUnderactuatedCells: Number(characterization.alphaUnderactuatedCells) || 0,
          heightUnderactuatedCells: Number(characterization.heightUnderactuatedCells) || 0,
          totalCells,
        }
      ),
      reportLawCandidate(
        "composition_nonadditivity",
        "actuation composition",
        "nonadditivity",
        "Operator composition is non-additive when the combined response exceeds the sum of individual responses by more than tolerance.",
        maxSuperpositionError > tolerance,
        {
          alphaSuperpositionError: Number(characterization.alphaSuperpositionError) || 0,
          heightSuperpositionError: Number(characterization.heightSuperpositionError) || 0,
          maxSuperpositionError,
          tolerance,
        }
      ),
      reportLawCandidate(
        "event_order_noncommutativity",
        "lock and actuation sequence",
        "noncommutativity",
        "Lock, release, and actuation events are noncommutative when reversal or adjacent swaps change the final state by more than tolerance.",
        Boolean(sequenceOrder?.orderSensitive),
        {
          eventCount,
          noncommutingAdjacentPairs: Number(sequenceOrder?.noncommutingAdjacentPairs) || 0,
          maxOrderError: Number(sequenceOrder?.maxOrderError) || 0,
          tolerance,
        }
      ),
    ];
    return {
      schema: "rad-sim.framework-law-candidates.v1",
      method: "thresholded diagnostic predicates over locality, reachability, composition, and event-order metrics",
      laws,
    };
  }

  function formalizationTarget(id, statement, source, status, readyForLean, dependencies, evidence) {
    return {
      id,
      statement,
      source,
      status,
      readyForLean: Boolean(readyForLean),
      dependencies,
      evidence,
    };
  }

  function frameworkFormalizationTargets(characterization, sequenceOrder, lockedCount, commandCount, tolerance) {
    const orderSensitive = Boolean(sequenceOrder?.orderSensitive);
    const targets = [
      formalizationTarget(
        "dead_zone_zero_inside_backlash",
        "For b >= 0, f_b(x)=max(0,x-b)+min(x+b,0) equals 0 whenever -b <= x <= b.",
        "paper-supported",
        "pending-lean-tooling",
        false,
        ["real max/min lemmas", "nonnegative backlash premise"],
        { formula: "f(x)=max(0,x-b)+min(x+b,0)" }
      ),
      formalizationTarget(
        "dead_zone_piecewise_linear_outside_gap",
        "For b >= 0, f_b(x)=x-b when x >= b and f_b(x)=x+b when x <= -b.",
        "paper-supported",
        "pending-lean-tooling",
        false,
        ["real max/min lemmas", "case split on backlash thresholds"],
        { formula: "f(x)=max(0,x-b)+min(x+b,0)" }
      ),
      formalizationTarget(
        "lock_projection_idempotent",
        "Applying the same lock projection twice is equivalent to applying it once.",
        "simulator-operator",
        "pending-lean-tooling",
        false,
        ["finite grid state model", "lock projection definition"],
        { lockedCellCount: lockedCount }
      ),
      formalizationTarget(
        "finite_response_rank_bound",
        "The rank of a finite response matrix is bounded by its command-column count.",
        "linear-algebra-diagnostic",
        "pending-lean-tooling",
        false,
        ["finite matrix rank theorem", "response matrix column count"],
        {
          alphaRank: Number(characterization.responseRankAlpha) || 0,
          heightRank: Number(characterization.responseRankHeight) || 0,
          commandCount,
        }
      ),
      formalizationTarget(
        "group_operator_support_decomposition",
        "A list-supported group operator acts only on the declared support, empty group support is neutral, and disjoint group supports commute in the abstract mode model.",
        "simulator-operator",
        "lean-proved-discrete",
        true,
        ["LocalModeOperator.groupOperator", "list membership", "disjoint support premise"],
        {
          leanTheorems: [
            "LocalModeOperator.groupOperator_support_inside",
            "LocalModeOperator.groupOperator_support_outside",
            "LocalModeOperator.groupOperator_empty_apply",
            "LocalModeOperator.groupOperator_append_support_left",
            "LocalModeOperator.groupOperator_append_support_right",
            "LocalModeOperator.groupOperators_with_disjoint_lists_commute",
          ],
          activeOperatorCount: commandCount,
        }
      ),
      formalizationTarget(
        "removed_cell_clears_group_supported_constraints",
        "If a constraint touches a removed cell, cell-removal constraint deletion makes that constraint inactive even after a group operator is applied.",
        "simulator-operator",
        "lean-proved-discrete",
        true,
        ["CellGraph.removeCellConstraints", "CellConstraintMap.touches", "LocalModeOperator.groupOperator"],
        {
          leanTheorems: [
            "CellGraph.removed_cell_constraints_clear_group_operator",
            "CellGraph.removed_cell_constraints_not_active_after_group_operator",
          ],
          activeOperatorCount: commandCount,
        }
      ),
      formalizationTarget(
        "vertical_clearance_gates_residual_contact",
        "A vertical command whose magnitude is inside pin-hole clearance has zero discrete residual transmission, and the clearance-excess contact penalty is nonnegative.",
        "contact-mechanics-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.clearanceExcessNat", "Mechanics.verticalResidualStepNat", "Mechanics.contactPenaltyFromClearanceNat"],
        {
          leanTheorems: [
            "clearanceExcessNat_zero_inside",
            "verticalResidualStepNat_zero_inside",
            "contactPenaltyFromClearanceNat_nonnegative",
          ],
          formula: "0.5*k*max(|z|-clearance,0)^2 in simulator diagnostics",
          claimLimit: "discrete theorem target plus uncalibrated normalized simulator penalty",
        }
      ),
      formalizationTarget(
        "fixed_cell_load_work_proxy_zero",
        "The finite load-work magnitude proxy is nonnegative, and a fixed cell with zero displacement contributes zero load work.",
        "load-mechanics-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.loadWorkMagnitudeNat", "fixed-cell zero-displacement premise"],
        {
          leanTheorems: [
            "loadWorkMagnitudeNat_nonnegative",
            "loadWorkMagnitudeNat_zero_fixed",
          ],
          claimLimit: "finite magnitude theorem plus uncalibrated simulator load-work proxy",
        }
      ),
      formalizationTarget(
        "signed_vertical_load_work_residual_int",
        "The signed finite load-work scaffold records the simulator sign convention: zero load or zero displacement gives zero signed work, and signed measured-versus-simulated residuals vanish when the signed quantities are equal.",
        "signed-load-validation-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.signedLoadWorkInt", "Mechanics.signedWorkResidualInt", "vertical-load comparison report schema"],
        {
          leanTheorems: [
            "signedLoadWorkInt_zero_displacement",
            "signedLoadWorkInt_zero_load",
            "absoluteErrorInt_self",
            "signedWorkResidualInt_zero_when_equal",
            "signedEnergyResidualTripleInt_nonnegative",
            "signedEnergyResidualTripleInt_zero_when_equal",
          ],
          pythonFunction: "write_vertical_load_energy_comparison_artifacts",
          cliModule: "rad_sim.compare_vertical_load_measurements",
          reportSchema: "rad-sim.vertical-load-energy-comparison-report.v1",
          claimLimit: "signed discrete validation scaffold; external load sign and hardware work remain calibration assumptions",
        }
      ),
      formalizationTarget(
        "integer_scaled_mechanics_scaffold",
        "Integer-scaled mechanics quantities carry an explicit denominator while preserving proved numerator facts: scaled spring/contact energies are nonnegative, zero displacement or penetration gives zero numerator, and scaled signed residual triples vanish when measured and simulated signed quantities agree.",
        "integer-scaled-mechanics-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ScaledNatQuantity", "Mechanics.ScaledIntQuantity", "future rational/real mechanics targets"],
        {
          leanTheorems: [
            "scaledSpringEnergyNat_numerator_nonnegative",
            "scaledSpringEnergyNat_zero_displacement",
            "scaledSpringEnergyNat_preserves_denominator",
            "scaledContactPenaltyNat_numerator_nonnegative",
            "scaledContactPenaltyNat_zero_penetration",
            "scaledSignedLoadWorkInt_zero_displacement",
            "scaledSignedLoadWorkInt_zero_load",
            "scaledSignedLoadWorkInt_preserves_denominator",
            "scaledSignedWorkResidualInt_numerator_nonnegative",
            "scaledSignedWorkResidualInt_zero_when_equal",
            "scaledSignedEnergyResidualTripleInt_numerator_nonnegative",
            "scaledSignedEnergyResidualTripleInt_zero_when_equal",
          ],
          claimLimit: "integer-scaled numerator/denominator scaffold; not yet a field-valued rational or real mechanics proof",
          nextFormalStep: "replace denominator-carrying records with Rat/Real semantics once algebra tooling is available",
        }
      ),
      formalizationTarget(
        "measurement_unit_scale_invariants",
        "A finite measurement-unit scale maps normalized simulator values into denominator-carrying physical-unit records while preserving zero commands, zero equality residuals, and denominator metadata for both nonnegative and signed quantities.",
        "calibration-unit-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.MeasurementUnitScale", "calibrate_paper_rad_config", "bench artifact unit metadata"],
        {
          leanTheorems: [
            "measurementUnitScaleNat_zero",
            "measurementUnitScaleNat_preserves_denominator",
            "measurementUnitScaleInt_zero",
            "measurementUnitScaleInt_preserves_denominator",
            "measurementUnitScaleResidualNat_zero_when_equal",
            "measurementUnitScaleResidualNat_preserves_denominator",
            "measurementUnitScaleResidualInt_zero_when_equal",
            "measurementUnitScaleResidualInt_preserves_denominator",
            "hardwareProfileCoverageComplete_true_when_equal",
            "hardwareProfileCoverageMissing_zero_when_complete",
            "hardwareProfileCoverageMissing_preserves_total",
          ],
          unitScaleFunction: "physical_unit_scale_metadata",
          unitScaleSchema: "rad-sim.physical-unit-scale-metadata.v1",
          hardwareProfileSchema: "rad-sim.hardware-profile.v1",
          hardwareProfileFunctions: [
            "RADHardwareProfile.to_dict",
            "hardware_profile_from_json",
            "export_hardware_profile_json",
          ],
          pythonFunction: "calibrate_paper_rad_config",
          benchPacketFunction: "write_vertical_load_bench_packet_artifacts",
          comparisonFunction: "write_vertical_load_energy_comparison_artifacts",
          claimLimit: "finite unit-scale metadata invariant; not a calibrated physical-units proof",
          nextFormalStep: "connect MeasurementUnitScale to Rat/Real unit conversion and measured hardware profiles",
        }
      ),
      formalizationTarget(
        "mechanics_energy_certificate_nonnegative_proxy",
        "The mechanical energy certificate separates a nonnegative stored/proxy total from signed external-work terms in the spring-hinge physical preview.",
        "variational-mechanics-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.mechanicalStoredEnergyNat", "spring-hinge solver energy breakdown", "mechanics energy certificate schema"],
        {
          leanTheorems: [
            "mechanicalStoredEnergyNat_nonnegative",
            "mechanicalStoredEnergyNat_zero_components",
            "springEnergyNat_nonnegative",
            "hingeEnergyNat_nonnegative",
            "contactPenaltyNat_nonnegative",
            "loadWorkMagnitudeNat_nonnegative",
          ],
          pythonFunction: "mechanics_energy_certificate_to_dict",
          schema: "rad-sim.mechanics-energy-certificate.v1",
          claimLimit: "nonnegative proxy scaffold; signed work and contact calibration remain physical assumptions",
        }
      ),
      formalizationTarget(
        "calibration_parameter_estimate_residual_bookkeeping",
        "A calibration parameter-estimate report records finite sample counts, raw residuals, fitted residuals, and tolerances so zero fitted residuals with at least one sample satisfy the finite pass predicate.",
        "bench-calibration-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.CalibrationFitResidualNat", "calibration_experiment_comparison_report", "rad-sim.calibration-parameter-estimates.v1"],
        {
          leanTheorems: [
            "calibrationFitResidualPassNat_zero",
            "calibrationFitResidualNat_preserves_sample_count",
          ],
          pythonFunction: "calibration_experiment_comparison_report",
          browserFunction: "calibrationParameterEstimates",
          reportSchema: "rad-sim.calibration-comparison-report.v1",
          parameterEstimateSchema: "rad-sim.calibration-parameter-estimates.v1",
          claimLimit: "finite residual bookkeeping only; fitted parameters remain empirical until separately validated",
        }
      ),
      formalizationTarget(
        "calibration_model_profile_safe_update_bounds",
        "A calibration-derived model-profile update is marked safe only when it has at least one finite sample and its proposed bounded simulator parameter lies between the declared lower and upper bounds.",
        "bench-calibration-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.CalibrationModelProfileUpdateNat", "calibration_model_profile_from_report", "rad-sim.calibration-model-profile.v1"],
        {
          leanTheorems: [
            "calibrationModelProfileUpdateSafeNat_intro",
            "calibrationModelProfileUpdateSafeNat_bounds",
          ],
          pythonFunction: "calibration_model_profile_from_report",
          pythonApplyFunction: "apply_calibration_model_profile",
          browserFunction: "calibrationModelProfile",
          browserApplyFunction: "applyCalibrationModelProfile",
          profileSchema: "rad-sim.calibration-model-profile.v1",
          applicationSchema: "rad-sim.calibration-model-profile-application.v1",
          claimLimit: "finite safety gate for bounded simulator updates; not a proof that fitted parameters are physical laws",
        }
      ),
      formalizationTarget(
        "calibration_model_profile_selection_predicate",
        "A calibration model-profile candidate is selectable only when it has an applied update, non-increased missing observations, and a non-increased finite residual score.",
        "bench-calibration-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.CalibrationModelProfileCandidateNat", "calibration_model_profile_residual_comparison", "rad-sim.calibration-model-profile-selection.v1"],
        {
          leanTheorems: [
            "calibrationModelProfileCandidateImprovesNat_intro",
            "calibrationModelProfileCandidateImprovesNat_applied_updates",
            "calibrationModelProfileCandidateImprovesNat_residual_nonincrease",
            "calibrationModelProfileCandidateImprovesNat_missing_nonincrease",
            "calibrationModelProfileCandidateScoreNat_le_before",
          ],
          pythonFunction: "select_calibration_model_profile",
          browserFunction: "selectCalibrationModelProfile",
          comparisonSchema: "rad-sim.calibration-model-profile-residual-comparison.v1",
          selectionSchema: "rad-sim.calibration-model-profile-selection.v1",
          claimLimit: "finite residual-selection rule only; independent validation remains required before claiming physical calibration",
        }
      ),
      formalizationTarget(
        "calibration_model_profile_holdout_predicate",
        "A train/holdout calibration profile validation passes only when a profile has an applied update and both fit and holdout residual scores and missing-observation counts do not increase.",
        "bench-calibration-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.CalibrationModelProfileHoldoutNat", "calibration_model_profile_holdout_validation", "rad-sim.calibration-model-profile-holdout-validation.v1"],
        {
          leanTheorems: [
            "calibrationModelProfileHoldoutPassNat_intro",
            "calibrationModelProfileHoldoutPassNat_applied_updates",
            "calibrationModelProfileHoldoutPassNat_fit_residual_nonincrease",
            "calibrationModelProfileHoldoutPassNat_holdout_residual_nonincrease",
            "calibrationModelProfileHoldoutPassNat_missing_nonincrease",
            "calibrationModelProfileHoldoutScoreNat_le_before",
          ],
          pythonFunction: "calibration_model_profile_holdout_validation",
          csvExportFunction: "export_calibration_model_profile_holdout_validation_csv",
          browserFunction: "calibrationModelProfileHoldoutValidation",
          schema: "rad-sim.calibration-model-profile-holdout-validation.v1",
          claimLimit: "finite train/holdout residual gate; still not a proof of contact, friction, material, or actuator physics",
        }
      ),
      formalizationTarget(
        "calibration_train_holdout_split_metadata",
        "A proper finite calibration split has positive fit samples, positive holdout samples, zero overlap, and a frozen profile before holdout evaluation.",
        "bench-calibration-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.CalibrationTrainHoldoutSplitNat", "rad-sim.calibration-train-holdout-split.v1", "bench data provenance metadata"],
        {
          leanTheorems: [
            "calibrationTrainHoldoutSplitReadyNat_intro",
            "calibrationTrainHoldoutSplitReadyNat_has_fit_samples",
            "calibrationTrainHoldoutSplitReadyNat_has_holdout_samples",
            "calibrationTrainHoldoutSplitReadyNat_zero_overlap",
            "calibrationTrainHoldoutSplitReadyNat_profile_frozen",
          ],
          pythonFunction: "calibration_model_profile_holdout_validation",
          csvExportFunction: "export_calibration_model_profile_holdout_validation_csv",
          browserFunction: "calibrationModelProfileHoldoutValidation",
          schema: "rad-sim.calibration-train-holdout-split.v1",
          claimLimit: "finite split metadata predicate; actual file independence and frozen-profile timing remain bench-protocol obligations",
        }
      ),
      formalizationTarget(
        "calibration_train_holdout_file_provenance",
        "A documented finite calibration holdout provenance record has known fit and holdout dataset IDs, distinct dataset and source-file IDs, marked fit/holdout roles, a frozen profile, and a holdout profile ID matching the frozen profile.",
        "bench-calibration-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.CalibrationTrainHoldoutProvenanceNat", "rad-sim.calibration-dataset-provenance.v1", "rad-sim.calibration-train-holdout-split.v1"],
        {
          leanTheorems: [
            "calibrationTrainHoldoutProvenanceReadyNat_intro",
            "calibrationTrainHoldoutProvenanceReadyNat_distinct_files",
            "calibrationTrainHoldoutProvenanceReadyNat_profile_frozen",
            "calibrationTrainHoldoutProvenanceReadyNat_profile_matches",
          ],
          pythonFunction: "calibration_model_profile_holdout_validation",
          browserFunction: "calibrationModelProfileHoldoutValidation",
          datasetProvenanceSchema: "rad-sim.calibration-dataset-provenance.v1",
          splitSchema: "rad-sim.calibration-train-holdout-split.v1",
          claimLimit: "finite provenance predicate; true file provenance and collection timing remain external laboratory obligations",
        }
      ),
      formalizationTarget(
        "calibration_bench_protocol_coverage",
        "A calibration bench notebook is coverage-ready only when it contains at least one scenario, at least one instrument, at least one output artifact, at least one measurement column, and both fit and holdout dataset roles.",
        "bench-calibration-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.CalibrationBenchProtocolCoverageNat", "calibration_bench_notebook", "rad-sim.calibration-bench-notebook.v1"],
        {
          leanTheorems: [
            "calibrationBenchProtocolCoverageReadyNat_intro",
            "calibrationBenchProtocolCoverageReadyNat_has_scenarios",
            "calibrationBenchProtocolCoverageReadyNat_has_outputs",
            "calibrationBenchProtocolCoverageReadyNat_has_two_dataset_roles",
          ],
          pythonFunction: "calibration_bench_notebook",
          csvExportFunction: "export_calibration_bench_notebook_csv",
          browserFunction: "calibrationBenchNotebook",
          schema: "rad-sim.calibration-bench-notebook.v1",
          claimLimit: "finite protocol-coverage predicate; actual independent collection, calibration, and physical-law validity remain laboratory obligations",
        }
      ),
      formalizationTarget(
        "calibration_bench_packet_completeness",
        "A calibration bench packet is complete only when it contains notebook artifacts, protocol artifacts, fit and holdout templates, validation instructions, and a sufficiently populated artifact manifest.",
        "bench-calibration-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.CalibrationBenchPacketCompletenessNat", "calibration_bench_packet", "rad-sim.calibration-bench-packet.v1"],
        {
          leanTheorems: [
            "calibrationBenchPacketCompleteNat_intro",
            "calibrationBenchPacketCompleteNat_has_notebook",
            "calibrationBenchPacketCompleteNat_has_fit_template",
            "calibrationBenchPacketCompleteNat_has_holdout_template",
            "calibrationBenchPacketCompleteNat_has_manifest",
          ],
          pythonFunction: "calibration_bench_packet",
          writerFunction: "write_calibration_bench_packet_artifacts",
          browserFunction: "calibrationBenchPacket",
          schema: "rad-sim.calibration-bench-packet.v1",
          claimLimit: "finite packet-completeness predicate; file independence, bench execution, and physical calibration remain external evidence",
        }
      ),
      formalizationTarget(
        "calibration_bench_executed_validation_gate",
        "A filled calibration bench execution is validation-ready only when fit and holdout samples are present, at least one bounded model-profile update was applied, residual and independent validation flags are positive, and no required provenance evidence is missing.",
        "bench-validation-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.CalibrationBenchExecutedValidationNat", "calibration_bench_execution_validation", "rad-sim.calibration-bench-execution-validation.v1"],
        {
          leanTheorems: [
            "calibrationBenchExecutedValidationReadyNat_intro",
            "calibrationBenchExecutedValidationReadyNat_has_fit_samples",
            "calibrationBenchExecutedValidationReadyNat_has_holdout_samples",
            "calibrationBenchExecutedValidationReadyNat_has_applied_updates",
            "calibrationBenchExecutedValidationReadyNat_has_independent_validation",
            "calibrationBenchExecutedValidationReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "calibration_bench_execution_validation",
          writerFunction: "write_calibration_bench_execution_validation_artifacts",
          browserFunction: "calibrationBenchExecutionValidation",
          csvExportFunction: "export_calibration_bench_execution_validation_csv",
          schema: "rad-sim.calibration-bench-execution-validation.v1",
          claimLimit: "finite execution-evidence predicate; lab procedure truth, measurement accuracy, and physical-law validity remain external evidence",
        }
      ),
      formalizationTarget(
        "measured_vertical_load_work_validation_zero_residual",
        "A measured-versus-simulated vertical load-work validation has zero finite residual when the measured finite work quantity equals the simulated quantity.",
        "bench-validation-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.measuredWorkResidualNat", "vertical-load energy validation schema", "load-cell height measurements"],
        {
          leanTheorems: [
            "absoluteErrorNat_self",
            "measuredWorkResidualNat_zero_when_equal",
            "measuredEnergyResidualTripleNat_zero_when_equal",
            "loadWorkMagnitudeNat_nonnegative",
          ],
          pythonFunction: "validate_vertical_load_energy_measurements",
          templateFunction: "vertical_load_energy_measurement_template",
          protocolFunction: "vertical_load_energy_experiment_protocol",
          packetFunction: "vertical_load_bench_packet",
          compareFunction: "compare_vertical_load_energy_measurement_results",
          schema: "rad-sim.vertical-load-energy-validation.v1",
          resultsSchema: "rad-sim.vertical-load-energy-measurement-results.v1",
          protocolSchema: "rad-sim.vertical-load-energy-experiment-protocol.v1",
          packetSchema: "rad-sim.vertical-load-bench-packet.v1",
          reportSchema: "rad-sim.vertical-load-energy-comparison-report.v1",
          claimLimit: "bench-data comparison scaffold; real contact/load calibration remains experimental",
        }
      ),
      formalizationTarget(
        "vertical_load_comparison_pass_predicate",
        "A filled vertical-load scenario with no missing measurements and zero signed-work, magnitude, and contact-proxy errors passes the finite comparison predicate for any tolerance; a scenario with a missing measurement fails that predicate.",
        "bench-validation-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.verticalLoadScenarioPassNat", "Mechanics.measuredEnergyResidualTripleNat", "vertical-load comparison report schema"],
        {
          leanTheorems: [
            "measuredEnergyResidualTripleNat_nonnegative",
            "measuredEnergyResidualTripleNat_zero_when_equal",
            "verticalLoadScenarioPassNat_zero_errors",
            "verticalLoadScenarioPassNat_false_when_missing",
          ],
          pythonFunction: "write_vertical_load_energy_comparison_artifacts",
          cliModule: "rad_sim.compare_vertical_load_measurements",
          summaryCsv: "vertical_load_energy_comparison_summary.csv",
          reportSchema: "rad-sim.vertical-load-energy-comparison-report.v1",
          claimLimit: "finite pass/fail predicate; physical accuracy still requires bench calibration",
        }
      ),
      formalizationTarget(
        "physical_validation_readiness_gate",
        "A physical validation evidence packet is readiness-complete only when executed calibration is present, vertical-load comparison scenarios are present and passing, load/contact proxy terms are present, clearance is configured, and no required evidence is missing.",
        "bench-validation-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.PhysicalValidationReadinessNat", "physical_validation_readiness_report", "rad-sim.physical-validation-readiness.v1"],
        {
          leanTheorems: [
            "physicalValidationReadyNat_intro",
            "physicalValidationReadyNat_has_calibration",
            "physicalValidationReadyNat_has_vertical_load",
            "physicalValidationReadyNat_has_load_proxy",
            "physicalValidationReadyNat_has_contact_proxy",
            "physicalValidationReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "physical_validation_readiness_report",
          browserFunction: "physicalValidationReadinessReport",
          csvExportFunction: "export_physical_validation_readiness_csv",
          schema: "rad-sim.physical-validation-readiness.v1",
          claimLimit: "finite evidence-readiness predicate; contact, friction, stiffness, gravity, and material-law accuracy remain external physical validation",
        }
      ),
      formalizationTarget(
        "contact_state_abstraction_gate",
        "A contact-state abstraction is complete only when active bodies have matching pin, hole, and clearance records, contact-state and penalty records cover the active bodies, clearance is configured, and no abstraction evidence is missing.",
        "contact-state-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ContactStateAbstractionNat", "contact_state_abstraction_report", "rad-sim.contact-state-abstraction.v1"],
        {
          leanTheorems: [
            "contactStateAbstractionReadyNat_intro",
            "contactStateAbstractionReadyNat_has_bodies",
            "contactStateAbstractionReadyNat_pins_match_bodies",
            "contactStateAbstractionReadyNat_holes_match_bodies",
            "contactStateAbstractionReadyNat_has_contact_records",
            "contactStateAbstractionReadyNat_has_penalty_terms",
            "contactStateAbstractionReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "contact_state_abstraction_report",
          browserFunction: "contactStateAbstractionReport",
          csvExportFunction: "export_contact_state_abstraction_csv",
          schema: "rad-sim.contact-state-abstraction.v1",
          claimLimit: "finite contact-state bookkeeping predicate; rigid-body contact, friction, and stiffness remain uncalibrated",
        }
      ),
      formalizationTarget(
        "contact_graph_consistency_gate",
        "A contact graph is consistency-complete only when removed-cell incident edges are deleted, removed cells have no active contact, active bodies have contact records, group-support cells are represented, and no consistency evidence is missing.",
        "contact-state-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ContactGraphConsistencyNat", "contact_graph_consistency_report", "rad-sim.contact-graph-consistency.v1"],
        {
          leanTheorems: [
            "contactGraphConsistentNat_intro",
            "contactGraphConsistentNat_has_active_bodies",
            "contactGraphConsistentNat_contact_records_cover_bodies",
            "contactGraphConsistentNat_removed_edges_clear",
            "contactGraphConsistentNat_removed_contacts_clear",
            "contactGraphConsistentNat_support_records_cover_support",
            "contactGraphConsistentNat_zero_missing_evidence",
          ],
          pythonFunction: "contact_graph_consistency_report",
          browserFunction: "contactGraphConsistencyReport",
          csvExportFunction: "export_contact_graph_consistency_csv",
          schema: "rad-sim.contact-graph-consistency.v1",
          claimLimit: "finite graph/contact/support bookkeeping predicate; physical contact remains uncalibrated",
        }
      ),
      formalizationTarget(
        "physical_realization_map_gate",
        "A physical realization map is readiness-complete only when at least one abstract operator is present, every abstract operator has a realized simulator state effect, support record, claim label, contact-graph evidence, and no missing realization evidence.",
        "physical-realization-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.PhysicalRealizationMapNat", "physical_realization_map_report", "rad-sim.physical-realization-map.v1"],
        {
          leanTheorems: [
            "physicalRealizationMapReadyNat_intro",
            "physicalRealizationMapReadyNat_has_operators",
            "physicalRealizationMapReadyNat_realizes_all",
            "physicalRealizationMapReadyNat_has_support_records",
            "physicalRealizationMapReadyNat_has_state_effects",
            "physicalRealizationMapReadyNat_has_claim_labels",
            "physicalRealizationMapReadyNat_has_contact_graph",
            "physicalRealizationMapReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "physical_realization_map_report",
          browserFunction: "physicalRealizationMapReport",
          csvExportFunction: "export_physical_realization_map_csv",
          schema: "rad-sim.physical-realization-map.v1",
          claimLimit: "finite abstract-to-simulator realization predicate; hardware realization remains experimentally unvalidated",
        }
      ),
      formalizationTarget(
        "external_physics_engine_audit_gate",
        "An external physics engine audit is readiness-complete only when an independent rigid-body/contact tool candidate is available, required features and validation scenarios are recorded, contact-model evidence exists, and no audit evidence is missing.",
        "external-physics-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ExternalPhysicsEngineAuditNat", "external_physics_engine_audit_report", "rad-sim.external-physics-engine-audit.v1"],
        {
          leanTheorems: [
            "externalPhysicsEngineAuditReadyNat_intro",
            "externalPhysicsEngineAuditReadyNat_has_engines",
            "externalPhysicsEngineAuditReadyNat_has_available",
            "externalPhysicsEngineAuditReadyNat_has_features",
            "externalPhysicsEngineAuditReadyNat_has_scenarios",
            "externalPhysicsEngineAuditReadyNat_has_contact_model",
            "externalPhysicsEngineAuditReadyNat_has_independent_tool",
            "externalPhysicsEngineAuditReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "external_physics_engine_audit_report",
          browserFunction: "externalPhysicsEngineAuditReport",
          csvExportFunction: "export_external_physics_engine_audit_csv",
          schema: "rad-sim.external-physics-engine-audit.v1",
          claimLimit: "finite external-engine readiness predicate; external solver execution and hardware agreement remain unvalidated",
        }
      ),
      formalizationTarget(
        "mujoco_model_export_gate",
        "A MuJoCo model export is readiness-complete only when a coarse external-engine model, active body records, fixed-body records, gravity records, MJCF XML bytes, and zero missing export evidence are present.",
        "external-physics-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ExternalPhysicsModelExportNat", "mujoco_model_export_report", "rad-sim.mujoco-model-export.v1"],
        {
          leanTheorems: [
            "externalPhysicsModelExportReadyNat_intro",
            "externalPhysicsModelExportReadyNat_has_model",
            "externalPhysicsModelExportReadyNat_has_bodies",
            "externalPhysicsModelExportReadyNat_has_fixed_bodies",
            "externalPhysicsModelExportReadyNat_has_gravity",
            "externalPhysicsModelExportReadyNat_has_xml",
            "externalPhysicsModelExportReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "mujoco_model_export_report",
          browserFunction: "mujocoModelExportReport",
          xmlExportFunction: "export_mujoco_model_xml",
          schema: "rad-sim.mujoco-model-export.v1",
          claimLimit: "finite MJCF export-completeness predicate; exported geometry is a coarse proxy, not validated CAD/contact mechanics",
        }
      ),
      formalizationTarget(
        "mujoco_pin_hole_contact_geometry_gate",
        "A MuJoCo pin-hole contact-geometry inventory is readiness-complete only when pin, hole, clearance, contact pair, active contact-pair, MJCF geometry, and zero missing geometry evidence records are present.",
        "external-physics-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ExternalContactGeometryNat", "mujoco_pin_hole_contact_geometry_report", "rad-sim.mujoco-pin-hole-contact-geometry.v1"],
        {
          leanTheorems: [
            "externalContactGeometryReadyNat_intro",
            "externalContactGeometryReadyNat_has_pins",
            "externalContactGeometryReadyNat_holes_cover_pins",
            "externalContactGeometryReadyNat_clearance_covers_pins",
            "externalContactGeometryReadyNat_pairs_cover_pins",
            "externalContactGeometryReadyNat_active_pairs_cover_pins",
            "externalContactGeometryReadyNat_has_xml",
            "externalContactGeometryReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "mujoco_pin_hole_contact_geometry_report",
          browserFunction: "mujocoPinHoleContactGeometryReport",
          csvExportFunction: "export_mujoco_pin_hole_contact_geometry_csv",
          schema: "rad-sim.mujoco-pin-hole-contact-geometry.v1",
          claimLimit: "finite pin-hole contact-geometry inventory; proxy geometry, friction, compliance, and hardware contact remain unvalidated",
        }
      ),
      formalizationTarget(
        "mujoco_contact_parameter_profile_gate",
        "A MuJoCo contact-parameter profile is readiness-complete only when contact pairs have parameter, friction, solver, stiffness, damping, XML-attribute, and zero missing evidence records.",
        "external-physics-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ExternalContactParameterProfileNat", "mujoco_contact_parameter_report", "rad-sim.mujoco-contact-parameter-profile.v1"],
        {
          leanTheorems: [
            "externalContactParameterProfileReadyNat_intro",
            "externalContactParameterProfileReadyNat_has_pairs",
            "externalContactParameterProfileReadyNat_parameters_cover_pairs",
            "externalContactParameterProfileReadyNat_friction_covers_pairs",
            "externalContactParameterProfileReadyNat_solver_covers_pairs",
            "externalContactParameterProfileReadyNat_stiffness_covers_pairs",
            "externalContactParameterProfileReadyNat_damping_covers_pairs",
            "externalContactParameterProfileReadyNat_has_xml_attributes",
            "externalContactParameterProfileReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "mujoco_contact_parameter_report",
          browserFunction: "mujocoContactParameterReport",
          csvExportFunction: "export_mujoco_contact_parameter_csv",
          schema: "rad-sim.mujoco-contact-parameter-profile.v1",
          claimLimit: "finite contact-parameter completeness predicate; values remain uncalibrated unless backed by bench measurements",
        }
      ),
      formalizationTarget(
        "contact_parameter_calibration_packet_completeness",
        "A contact-parameter calibration packet is readiness-complete only when contact-pair records, parameter records, measurement columns, fit and holdout template rows, artifact-manifest entries, attached profile evidence, and zero missing evidence are present.",
        "bench-calibration-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ContactParameterCalibrationPacketNat", "contact_parameter_calibration_packet", "rad-sim.contact-parameter-calibration-packet.v1"],
        {
          leanTheorems: [
            "contactParameterCalibrationPacketCompleteNat_intro",
            "contactParameterCalibrationPacketCompleteNat_has_pairs",
            "contactParameterCalibrationPacketCompleteNat_parameters_cover_pairs",
            "contactParameterCalibrationPacketCompleteNat_has_measurement_columns",
            "contactParameterCalibrationPacketCompleteNat_fit_rows_cover_pairs",
            "contactParameterCalibrationPacketCompleteNat_holdout_rows_cover_pairs",
            "contactParameterCalibrationPacketCompleteNat_has_manifest",
            "contactParameterCalibrationPacketCompleteNat_has_profile_evidence",
            "contactParameterCalibrationPacketCompleteNat_zero_missing_evidence",
          ],
          pythonFunction: "contact_parameter_calibration_packet",
          browserFunction: "contactParameterCalibrationPacket",
          csvExportFunction: "export_contact_parameter_calibration_packet_csv",
          schema: "rad-sim.contact-parameter-calibration-packet.v1",
          claimLimit: "finite contact-calibration packet completeness predicate; fitted contact physics remains external bench evidence",
        }
      ),
      formalizationTarget(
        "contact_parameter_bench_validation_gate",
        "A contact-parameter bench validation is readiness-complete only when a ready packet, fit rows, holdout rows, completed measurements, residual checks, fit pass, holdout pass, independent holdout pass, and zero missing evidence are present.",
        "bench-calibration-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ContactParameterBenchValidationNat", "compare_contact_parameter_calibration_results", "rad-sim.contact-parameter-bench-validation.v1"],
        {
          leanTheorems: [
            "contactParameterBenchValidationReadyNat_intro",
            "contactParameterBenchValidationReadyNat_has_packet",
            "contactParameterBenchValidationReadyNat_has_fit_rows",
            "contactParameterBenchValidationReadyNat_has_holdout_rows",
            "contactParameterBenchValidationReadyNat_measurements_cover_rows",
            "contactParameterBenchValidationReadyNat_has_residuals",
            "contactParameterBenchValidationReadyNat_has_fit_pass",
            "contactParameterBenchValidationReadyNat_has_holdout_pass",
            "contactParameterBenchValidationReadyNat_has_independent_holdout",
            "contactParameterBenchValidationReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "compare_contact_parameter_calibration_results",
          browserFunction: "compareContactParameterCalibrationResults",
          csvExportFunction: "export_contact_parameter_bench_validation_csv",
          schema: "rad-sim.contact-parameter-bench-validation.v1",
          resultsSchema: "rad-sim.contact-parameter-calibration-results.v1",
          claimLimit: "finite contact-parameter bench comparison predicate; physical contact laws remain empirical and require independent measurements",
        }
      ),
      formalizationTarget(
        "contact_parameter_interval_calibration_gate",
        "A contact-parameter interval calibration is readiness-complete only when a passing bench validation, parameter intervals, accepted interval records, uncertainty records, holdout agreement records, simulator parameters inside bounds, and zero missing evidence are present.",
        "bench-calibration-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ContactParameterIntervalCalibrationNat", "contact_parameter_interval_calibration_report", "rad-sim.contact-parameter-interval-calibration.v1"],
        {
          leanTheorems: [
            "contactParameterIntervalCalibrationReadyNat_intro",
            "contactParameterIntervalCalibrationReadyNat_has_validation",
            "contactParameterIntervalCalibrationReadyNat_has_intervals",
            "contactParameterIntervalCalibrationReadyNat_accepts_all_intervals",
            "contactParameterIntervalCalibrationReadyNat_has_uncertainty",
            "contactParameterIntervalCalibrationReadyNat_has_holdout_agreement",
            "contactParameterIntervalCalibrationReadyNat_parameters_inside_bounds",
            "contactParameterIntervalCalibrationReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "contact_parameter_interval_calibration_report",
          browserFunction: "contactParameterIntervalCalibrationReport",
          csvExportFunction: "export_contact_parameter_interval_calibration_csv",
          schema: "rad-sim.contact-parameter-interval-calibration.v1",
          upstreamSchema: "rad-sim.contact-parameter-bench-validation.v1",
          claimLimit: "finite empirical interval-calibration predicate; interval acceptance is not a constitutive contact proof",
        }
      ),
      formalizationTarget(
        "mujoco_external_run_gate",
        "A MuJoCo external run is readiness-complete only when the engine is available, a ready model export is attached, body records are present, result records cover body records, requested steps are positive and completed, and no run evidence is missing.",
        "external-physics-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ExternalPhysicsRunNat", "mujoco_external_run_report", "rad-sim.mujoco-external-run.v1"],
        {
          leanTheorems: [
            "externalPhysicsRunReadyNat_intro",
            "externalPhysicsRunReadyNat_has_engine",
            "externalPhysicsRunReadyNat_has_export",
            "externalPhysicsRunReadyNat_has_bodies",
            "externalPhysicsRunReadyNat_results_cover_bodies",
            "externalPhysicsRunReadyNat_has_steps",
            "externalPhysicsRunReadyNat_steps_completed",
            "externalPhysicsRunReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "mujoco_external_run_report",
          browserFunction: "mujocoExternalRunReport",
          schema: "rad-sim.mujoco-external-run.v1",
          claimLimit: "finite external-run completeness predicate; external solver output still needs simulator and bench comparison",
        }
      ),
      formalizationTarget(
        "mujoco_external_comparison_gate",
        "A MuJoCo external comparison is readiness-complete only when a complete external run is attached, simulator records are present, external and matched records cover simulator records, a tolerance record exists, the comparison is within tolerance, and no comparison evidence is missing.",
        "external-physics-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ExternalPhysicsComparisonNat", "mujoco_external_comparison_report", "rad-sim.mujoco-external-comparison.v1"],
        {
          leanTheorems: [
            "externalPhysicsComparisonReadyNat_intro",
            "externalPhysicsComparisonReadyNat_has_run",
            "externalPhysicsComparisonReadyNat_has_simulator_records",
            "externalPhysicsComparisonReadyNat_external_covers_simulator",
            "externalPhysicsComparisonReadyNat_matched_covers_simulator",
            "externalPhysicsComparisonReadyNat_has_tolerance",
            "externalPhysicsComparisonReadyNat_within_tolerance",
            "externalPhysicsComparisonReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "mujoco_external_comparison_report",
          browserFunction: "mujocoExternalComparisonReport",
          csvExportFunction: "export_mujoco_external_comparison_csv",
          schema: "rad-sim.mujoco-external-comparison.v1",
          claimLimit: "finite simulator-vs-external comparison predicate; bench validation remains external evidence",
        }
      ),
      formalizationTarget(
        "equilibrium_relation_gate",
        "An equilibrium relation is readiness-complete only when physical realization evidence is ready, solver success and contact-state evidence are present, energy terms are covered by nonnegative terms, the residual is within tolerance, and no equilibrium evidence is missing.",
        "variational-mechanics-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.EquilibriumRelationNat", "equilibrium_relation_report", "rad-sim.equilibrium-relation.v1"],
        {
          leanTheorems: [
            "equilibriumRelationReadyNat_intro",
            "equilibriumRelationReadyNat_has_realization",
            "equilibriumRelationReadyNat_has_solver",
            "equilibriumRelationReadyNat_has_contact",
            "equilibriumRelationReadyNat_has_energy_terms",
            "equilibriumRelationReadyNat_energy_terms_nonnegative",
            "equilibriumRelationReadyNat_residual_within_tolerance",
            "equilibriumRelationReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "equilibrium_relation_report",
          browserFunction: "equilibriumRelationReport",
          csvExportFunction: "export_equilibrium_relation_csv",
          schema: "rad-sim.equilibrium-relation.v1",
          claimLimit: "finite equilibrium evidence predicate; continuous minimization and physical material/contact accuracy remain unproved",
        }
      ),
      formalizationTarget(
        "reachable_equilibrium_controllability_gate",
        "A reachable-equilibrium controllability report is readiness-complete only when equilibrium evidence is ready, an actuator basis and response columns exist, target cells and reachable responses are recorded, topology evidence is present, optional full-target reachability policy is satisfied, and no evidence is missing.",
        "reachability-controllability-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ReachableEquilibriumControllabilityNat", "reachable_equilibrium_controllability_report", "rad-sim.reachable-equilibrium-controllability.v1"],
        {
          leanTheorems: [
            "reachableEquilibriumControllabilityReadyNat_intro",
            "reachableEquilibriumControllabilityReadyNat_has_equilibrium",
            "reachableEquilibriumControllabilityReadyNat_has_actuator_basis",
            "reachableEquilibriumControllabilityReadyNat_has_response_columns",
            "reachableEquilibriumControllabilityReadyNat_has_target_cells",
            "reachableEquilibriumControllabilityReadyNat_has_reachable_responses",
            "reachableEquilibriumControllabilityReadyNat_has_topology",
            "reachableEquilibriumControllabilityReadyNat_target_policy",
            "reachableEquilibriumControllabilityReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "reachable_equilibrium_controllability_report",
          browserFunction: "reachableEquilibriumControllabilityReport",
          csvExportFunction: "export_reachable_equilibrium_controllability_csv",
          schema: "rad-sim.reachable-equilibrium-controllability.v1",
          claimLimit: "finite response-matrix reachability predicate; nonlinear controllability and hardware feasibility remain unproved",
        }
      ),
      formalizationTarget(
        "reachable_equilibrium_bench_protocol_gate",
        "A reachable-equilibrium bench protocol is readiness-complete only when reachable-equilibrium evidence is attached, protocol steps and actuator-column trials exist, target observation cells and measurement columns are recorded, topology-blocked target policy is represented, pass/fail criteria are present, and no protocol evidence is missing.",
        "physical-experiment-protocol-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ReachableEquilibriumBenchProtocolNat", "reachable_equilibrium_bench_protocol", "rad-sim.reachable-equilibrium-bench-protocol.v1"],
        {
          leanTheorems: [
            "reachableEquilibriumBenchProtocolReadyNat_intro",
            "reachableEquilibriumBenchProtocolReadyNat_has_reachability",
            "reachableEquilibriumBenchProtocolReadyNat_has_steps",
            "reachableEquilibriumBenchProtocolReadyNat_has_actuator_trials",
            "reachableEquilibriumBenchProtocolReadyNat_has_targets",
            "reachableEquilibriumBenchProtocolReadyNat_has_measurement_columns",
            "reachableEquilibriumBenchProtocolReadyNat_has_topology_policy",
            "reachableEquilibriumBenchProtocolReadyNat_has_pass_fail",
            "reachableEquilibriumBenchProtocolReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "reachable_equilibrium_bench_protocol",
          browserFunction: "reachableEquilibriumBenchProtocol",
          csvExportFunction: "export_reachable_equilibrium_bench_protocol_csv",
          schema: "rad-sim.reachable-equilibrium-bench-protocol.v1",
          claimLimit: "bench protocol coverage predicate; physical measurements and nonlinear controllability remain unproved",
        }
      ),
      formalizationTarget(
        "reachable_equilibrium_bench_validation_gate",
        "A reachable-equilibrium bench validation is readiness-complete only when the protocol is ready, result rows and target measurements exist, completed measurements cover target measurements, reachability, topology-leakage, and group-sequence checks are present, the comparison passes, and no validation evidence is missing.",
        "physical-experiment-validation-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ReachableEquilibriumBenchValidationNat", "compare_reachable_equilibrium_bench_results", "rad-sim.reachable-equilibrium-bench-comparison.v1"],
        {
          leanTheorems: [
            "reachableEquilibriumBenchValidationReadyNat_intro",
            "reachableEquilibriumBenchValidationReadyNat_has_protocol",
            "reachableEquilibriumBenchValidationReadyNat_has_rows",
            "reachableEquilibriumBenchValidationReadyNat_has_targets",
            "reachableEquilibriumBenchValidationReadyNat_measurements_complete",
            "reachableEquilibriumBenchValidationReadyNat_has_reachability_checks",
            "reachableEquilibriumBenchValidationReadyNat_has_topology_checks",
            "reachableEquilibriumBenchValidationReadyNat_has_group_checks",
            "reachableEquilibriumBenchValidationReadyNat_has_pass",
            "reachableEquilibriumBenchValidationReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "compare_reachable_equilibrium_bench_results",
          templateFunction: "reachable_equilibrium_bench_results_template",
          browserFunction: "compareReachableEquilibriumBenchResults",
          csvExportFunction: "export_reachable_equilibrium_bench_comparison_csv",
          schema: "rad-sim.reachable-equilibrium-bench-comparison.v1",
          resultsSchema: "rad-sim.reachable-equilibrium-bench-results.v1",
          claimLimit: "finite filled-result validation predicate; physical laws and nonlinear controllability remain unproved",
        }
      ),
      formalizationTarget(
        "reachable_equilibrium_amplitude_calibration_gate",
        "A reachable-equilibrium amplitude calibration is readiness-complete only when bench validation is ready, repeated-trial groups, amplitude estimates, residual-field cells, uncertainty bands, topology-leakage bands, group-sequence residuals, and zero missing evidence are present.",
        "physical-amplitude-calibration-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ReachableEquilibriumAmplitudeCalibrationNat", "reachable_equilibrium_amplitude_calibration_report", "rad-sim.reachable-equilibrium-amplitude-calibration.v1"],
        {
          leanTheorems: [
            "reachableEquilibriumAmplitudeCalibrationReadyNat_intro",
            "reachableEquilibriumAmplitudeCalibrationReadyNat_has_validation",
            "reachableEquilibriumAmplitudeCalibrationReadyNat_has_repeated_trials",
            "reachableEquilibriumAmplitudeCalibrationReadyNat_has_amplitudes",
            "reachableEquilibriumAmplitudeCalibrationReadyNat_has_residual_field",
            "reachableEquilibriumAmplitudeCalibrationReadyNat_has_uncertainty",
            "reachableEquilibriumAmplitudeCalibrationReadyNat_has_topology_bands",
            "reachableEquilibriumAmplitudeCalibrationReadyNat_has_group_residuals",
            "reachableEquilibriumAmplitudeCalibrationReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "reachable_equilibrium_amplitude_calibration_report",
          browserFunction: "reachableEquilibriumAmplitudeCalibrationReport",
          csvExportFunction: "export_reachable_equilibrium_amplitude_calibration_csv",
          schema: "rad-sim.reachable-equilibrium-amplitude-calibration.v1",
          claimLimit: "finite repeated-trial amplitude artifact; statistical model, physical laws, and nonlinear controllability remain unproved",
        }
      ),
      formalizationTarget(
        "reachable_equilibrium_empirical_profile_gate",
        "A reachable-equilibrium empirical profile is readiness-complete only when amplitude calibration evidence is ready, bounded and safe update proposals exist, tolerance and uncertainty records are present, a holdout hook is represented, and no profile evidence is missing.",
        "physical-empirical-profile-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ReachableEquilibriumEmpiricalProfileNat", "reachable_equilibrium_empirical_profile_from_amplitude", "rad-sim.reachable-equilibrium-empirical-profile.v1"],
        {
          leanTheorems: [
            "reachableEquilibriumEmpiricalProfileReadyNat_intro",
            "reachableEquilibriumEmpiricalProfileReadyNat_has_amplitude",
            "reachableEquilibriumEmpiricalProfileReadyNat_has_bounded_proposals",
            "reachableEquilibriumEmpiricalProfileReadyNat_has_safe_proposals",
            "reachableEquilibriumEmpiricalProfileReadyNat_has_tolerances",
            "reachableEquilibriumEmpiricalProfileReadyNat_has_uncertainty",
            "reachableEquilibriumEmpiricalProfileReadyNat_has_holdout_hooks",
            "reachableEquilibriumEmpiricalProfileReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "reachable_equilibrium_empirical_profile_from_amplitude",
          browserFunction: "reachableEquilibriumEmpiricalProfileFromAmplitude",
          csvExportFunction: "export_reachable_equilibrium_empirical_profile_csv",
          schema: "rad-sim.reachable-equilibrium-empirical-profile.v1",
          claimLimit: "bounded empirical profile proposal artifact; no v1 LatticeConfig mutation or physical-law proof",
        }
      ),
      formalizationTarget(
        "reachable_equilibrium_profile_inverse_gate",
        "A profile-aware inverse diagnostic is readiness-complete only when the empirical profile is ready, safe proposals and target cells exist, the inverse solve succeeded, residual and score records are present, the profile is used read-only, and no evidence is missing.",
        "profile-aware-inverse-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ReachableEquilibriumProfileInverseNat", "reachable_equilibrium_profile_inverse_report", "rad-sim.reachable-equilibrium-profile-inverse.v1"],
        {
          leanTheorems: [
            "reachableEquilibriumProfileInverseReadyNat_intro",
            "reachableEquilibriumProfileInverseReadyNat_has_profile",
            "reachableEquilibriumProfileInverseReadyNat_has_safe_proposals",
            "reachableEquilibriumProfileInverseReadyNat_has_targets",
            "reachableEquilibriumProfileInverseReadyNat_has_solve",
            "reachableEquilibriumProfileInverseReadyNat_has_residuals",
            "reachableEquilibriumProfileInverseReadyNat_has_scores",
            "reachableEquilibriumProfileInverseReadyNat_read_only_profile_use",
            "reachableEquilibriumProfileInverseReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "reachable_equilibrium_profile_inverse_report",
          browserFunction: "reachableEquilibriumProfileInverseReport",
          csvExportFunction: "export_reachable_equilibrium_profile_inverse_csv",
          schema: "rad-sim.reachable-equilibrium-profile-inverse.v1",
          claimLimit: "profile-aware inverse residual certificate; no optimizer mutation, physical-law proof, or fabrication calibration",
        }
      ),
      formalizationTarget(
        "reachable_equilibrium_profile_inverse_acceptance_gate",
        "A profile-aware inverse acceptance decision is ready only when the profile-inverse report is ready, score records exist, the scaled residual score, band-failure count, and actuator count fit declared limits, the plan is accepted for preview, and no evidence is missing.",
        "profile-aware-inverse-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ReachableEquilibriumProfileInverseAcceptanceNat", "reachable_equilibrium_profile_inverse_acceptance_report", "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1"],
        {
          leanTheorems: [
            "reachableEquilibriumProfileInverseAcceptanceReadyNat_intro",
            "reachableEquilibriumProfileInverseAcceptanceReadyNat_has_profile_inverse",
            "reachableEquilibriumProfileInverseAcceptanceReadyNat_has_score_records",
            "reachableEquilibriumProfileInverseAcceptanceReadyNat_score_within_limit",
            "reachableEquilibriumProfileInverseAcceptanceReadyNat_band_failures_within_limit",
            "reachableEquilibriumProfileInverseAcceptanceReadyNat_actuators_within_limit",
            "reachableEquilibriumProfileInverseAcceptanceReadyNat_accepted",
            "reachableEquilibriumProfileInverseAcceptanceReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "reachable_equilibrium_profile_inverse_acceptance_report",
          browserFunction: "reachableEquilibriumProfileInverseAcceptanceReport",
          csvExportFunction: "export_reachable_equilibrium_profile_inverse_acceptance_csv",
          schema: "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1",
          claimLimit: "thresholded preview acceptance gate; not a hardware execution authorization or nonlinear controllability proof",
        }
      ),
      formalizationTarget(
        "reachable_equilibrium_profile_inverse_preview_packet_gate",
        "An accepted profile-aware inverse preview packet is ready only when acceptance evidence is ready, command records exist, preview event records cover commands, target and residual records exist, the packet is marked read-only, and no evidence is missing.",
        "profile-aware-inverse-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ReachableEquilibriumProfileInversePreviewPacketNat", "reachable_equilibrium_profile_inverse_preview_packet", "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1"],
        {
          leanTheorems: [
            "reachableEquilibriumProfileInversePreviewPacketReadyNat_intro",
            "reachableEquilibriumProfileInversePreviewPacketReadyNat_has_acceptance",
            "reachableEquilibriumProfileInversePreviewPacketReadyNat_has_commands",
            "reachableEquilibriumProfileInversePreviewPacketReadyNat_events_cover_commands",
            "reachableEquilibriumProfileInversePreviewPacketReadyNat_has_targets",
            "reachableEquilibriumProfileInversePreviewPacketReadyNat_has_residuals",
            "reachableEquilibriumProfileInversePreviewPacketReadyNat_read_only",
            "reachableEquilibriumProfileInversePreviewPacketReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "reachable_equilibrium_profile_inverse_preview_packet",
          browserFunction: "reachableEquilibriumProfileInversePreviewPacket",
          csvExportFunction: "export_reachable_equilibrium_profile_inverse_preview_packet_csv",
          schema: "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1",
          claimLimit: "read-only preview/lab handoff artifact; not hardware execution authorization",
        }
      ),
      formalizationTarget(
        "reachable_equilibrium_profile_inverse_preview_replay_gate",
        "A profile-aware inverse preview replay is ready only when the packet is ready, command records replay exactly, preview events cover the commands, simulation and target-residual records exist, residuals agree with packet metadata, replay is read-only, and no evidence is missing.",
        "profile-aware-inverse-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ReachableEquilibriumProfileInversePreviewReplayNat", "reachable_equilibrium_profile_inverse_preview_replay_report", "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1"],
        {
          leanTheorems: [
            "reachableEquilibriumProfileInversePreviewReplayReadyNat_intro",
            "reachableEquilibriumProfileInversePreviewReplayReadyNat_has_packet",
            "reachableEquilibriumProfileInversePreviewReplayReadyNat_has_commands",
            "reachableEquilibriumProfileInversePreviewReplayReadyNat_replays_all_commands",
            "reachableEquilibriumProfileInversePreviewReplayReadyNat_has_event_coverage",
            "reachableEquilibriumProfileInversePreviewReplayReadyNat_has_simulation",
            "reachableEquilibriumProfileInversePreviewReplayReadyNat_has_target_residuals",
            "reachableEquilibriumProfileInversePreviewReplayReadyNat_has_residual_agreement",
            "reachableEquilibriumProfileInversePreviewReplayReadyNat_read_only",
            "reachableEquilibriumProfileInversePreviewReplayReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "reachable_equilibrium_profile_inverse_preview_replay_report",
          browserFunction: "reachableEquilibriumProfileInversePreviewReplayReport",
          csvExportFunction: "export_reachable_equilibrium_profile_inverse_preview_replay_csv",
          schema: "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1",
          claimLimit: "deterministic kinematic replay certificate; not physical contact validation or hardware execution",
        }
      ),
      formalizationTarget(
        "reachable_equilibrium_profile_inverse_preview_physical_gate",
        "A profile-aware inverse preview physical check is ready only when the replay certificate is ready, the physical-preview solver succeeds, command and target-residual records exist, model-comparison and energy records exist, the check is read-only, and no evidence is missing.",
        "profile-aware-inverse-scaffold",
        "lean-proved-discrete",
        true,
        ["Mechanics.ReachableEquilibriumProfileInversePreviewPhysicalNat", "reachable_equilibrium_profile_inverse_preview_physical_report", "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1"],
        {
          leanTheorems: [
            "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_intro",
            "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_replay",
            "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_solver",
            "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_commands",
            "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_targets",
            "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_model_comparison",
            "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_has_energy",
            "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_read_only",
            "reachableEquilibriumProfileInversePreviewPhysicalReadyNat_zero_missing_evidence",
          ],
          pythonFunction: "reachable_equilibrium_profile_inverse_preview_physical_report",
          browserFunction: "reachableEquilibriumProfileInversePreviewPhysicalReport",
          csvExportFunction: "export_reachable_equilibrium_profile_inverse_preview_physical_csv",
          schema: "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1",
          claimLimit: "normalized spring-hinge/browser spring-preview check; not a calibrated contact, gravity, or hardware-execution proof",
        }
      ),
      formalizationTarget(
        "spring_hinge_removed_topology_load_comparison",
        "A removed-cell vertical residual/load case can be compared against the 3D spring-hinge physical-preview solver using the same deleted spring/hinge topology and fixed-cell load support.",
        "spring-hinge-diagnostic",
        "simulator-diagnostic",
        false,
        ["compare_vertical_residual_spring_hinge_3d", "solve_spring_hinge_3d", "LoadCase.fixed_cells"],
        {
          pythonFunction: "compare_vertical_residual_spring_hinge_3d",
          reportFunction: "build_vertical_load_physical_preview_report",
          csvExportFunction: "export_vertical_load_physical_preview_report_csv",
          comparisonSchema: "rad-sim.vertical-removal-physical-comparison.v1",
          reportSchema: "rad-sim.vertical-load-physical-preview-report.v1",
          solver: "solve_spring_hinge_3d",
          claimLimit: "numerical physical preview; not a calibrated rigid-body contact proof",
        }
      ),
      formalizationTarget(
        "browser_spring_preview_removed_edge_deletion",
        "The browser spring-preview relaxation treats removed cells as deleted graph nodes: removed cells do not average with active neighbors, and incident link-strain entries are marked as absent spring edges.",
        "spring-preview-diagnostic",
        "simulator-diagnostic",
        false,
        ["CellGraph.removeCell", "web/physics.js simulatePhysicalRelaxation", "link-strain absent-edge counts"],
        {
          leanTheorems: [
            "CellGraph.removed_cell_not_present",
            "CellGraph.removal_deletes_one_step_from_removed",
            "CellGraph.removal_deletes_one_step_to_removed",
          ],
          browserFunction: "RAD.simulatePhysicalRelaxation",
          browserMetrics: [
            "physicalActiveSpringEdges",
            "physicalSkippedSpringEdges",
          ],
          validation: "tests/validate_web_modules.js",
          claimLimit: "browser relaxation topology consistency only; not calibrated rigid-body contact",
        }
      ),
      formalizationTarget(
        "noncommutativity_witness_from_order_error",
        "If sequence-order distance is greater than tolerance, the corresponding event compositions are not equal.",
        "simulator-diagnostic",
        orderSensitive ? "pending-lean-tooling" : "pending-numeric-witness",
        false,
        ["state distance definition", "event composition semantics"],
        {
          orderSensitive,
          maxOrderError: Number(sequenceOrder?.maxOrderError) || 0,
          tolerance,
        }
      ),
      formalizationTarget(
        "bounded_locality_witness",
        "If all response magnitudes outside a reported die-off radius are below tolerance, the diagnostic has a finite locality witness.",
        "simulator-diagnostic",
        "requires-calibrated-premise",
        false,
        ["normed response field", "thresholded locality definition"],
        {
          alphaLocalityRadius: Number(characterization.alphaDieOff) || 0,
          zLocalityRadius: Number(characterization.zDieOff) || 0,
          tolerance,
        }
      ),
    ];
    return {
      schema: "rad-sim.formalization-targets.v1",
      method: "candidate theorem manifest; no Lean proof is emitted until Lean/Lake are available and the premises are first-principles enough to formalize",
      tooling: {
        engine: "Lean",
        leanPath: null,
        lakePath: null,
        available: false,
        status: "browser-cannot-inspect-path",
      },
      targets,
    };
  }

  function programmableDiscontinuityReport(state, options = {}) {
    const tolerance = Number(options.tolerance ?? 1e-9);
    const includeFields = options.includeFields !== false;
    const includeResponseMatrix = options.includeResponseMatrix !== false;
    const scope = options.scope || state.experiment?.characterizationScope || "single";
    const selected = {
      r: clampIndex(options.r ?? state.selection?.r ?? 0, state.grid.rows),
      c: clampIndex(options.c ?? state.selection?.c ?? 0, state.grid.cols),
    };
    const sourceCells = uniqueProtocolCells(
      characterizationCells(state, selected.r, selected.c, scope).map((cell) => normalizeMatrixCell(state, cell))
    );
    const combinedState = scopedState(state, sourceCells);
    const baselineState = scopedState(state, []);
    const sim = RAD.simulate(combinedState);
    const baselineSim = RAD.simulate(baselineState);
    const characterization = characterizeLocalResponse(state, { scope, r: selected.r, c: selected.c });
    const matrix = buildResponseMatrix(state, { scope, r: selected.r, c: selected.c, tolerance });
    const commands = activeReportCommands(combinedState, sourceCells, tolerance);
    const activeOperatorCount = commands.filter((command) => command.active).length;
    const order = reportSequenceOrder(combinedState, commands, tolerance);
    const calibration = RAD.paperRadCalibration(combinedState);
    const responseMatrix = includeResponseMatrix
      ? matrix
      : {
          schema: matrix.schema,
          grid: matrix.grid,
          source: matrix.source,
          commands: matrix.commands,
          diagnostics: matrix.diagnostics,
        };
    const report = {
      schema: "rad-sim.programmable-discontinuity-report.v1",
      savedAt: new Date().toISOString(),
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        totalCells: state.grid.rows * state.grid.cols,
      },
      config: {
        cellSize: state.grid.cellSize,
        initialAlpha: state.grid.initialAlpha,
        backlash: state.grid.backlash,
        couplingGain: state.grid.couplingGain,
        zCouplingGain: state.grid.zCouplingGain,
        pinRadius: state.grid.pinRadius,
        holeRadius: state.grid.holeRadius,
        pinHoleClearance: RAD.pinHoleClearance(state),
        backlashMm: calibration.configuredBacklashMm,
        pinHoleClearanceMm: calibration.pinHoleClearanceMm,
        alphaContractLimit: state.grid.alphaContractLimit,
        alphaExpandLimit: state.grid.alphaExpandLimit,
        zTravelLimit: state.grid.zTravelLimit,
      },
      operators: {
        scope,
        selected: reportCell(selected),
        commands,
        lockedCells: lockedReportCells(combinedState),
        eventSequence: order.eventSequence,
        activeOperatorCount,
      },
      paperSupportedAssumptions: [
        {
          name: "backlash dead-zone activation",
          formula: "f(x)=max(0,x-b)+min(x+b,0)",
          implementation: "RAD.backlashActivation",
        },
        {
          name: "normalized backlash",
          formula: "b_norm=b/L",
          implementation: "state.grid.backlash is dimensionless in v1",
        },
        {
          name: "rotating-square angle/dilation relation",
          formula: "theta_degrees=70*alpha-60",
          implementation: "RAD.alphaToTheta",
        },
      ],
      simulatorDiagnostics: [
        {
          name: "response matrix reachability",
          interpretation: "finite command columns approximate local reachable alpha/height directions",
        },
        {
          name: "superposition residual",
          interpretation: "nonzero residual marks non-additive operator composition caused by thresholds, locks, saturation, or coupling",
        },
        {
          name: "shellwise locality fit",
          interpretation: "log-linear shell maxima estimate die-off but are not a constitutive law",
        },
        {
          name: "event-order sensitivity",
          interpretation: "reversal and adjacent-swap differences test noncommutativity of lock and actuation operators",
        },
      ],
      operatorLawCandidates: frameworkOperatorLawCandidates(
        characterization,
        order.sequenceOrder,
        order.eventSequence.length,
        state.grid.rows * state.grid.cols,
        tolerance
      ),
      formalizationTargets: frameworkFormalizationTargets(
        characterization,
        order.sequenceOrder,
        lockedReportCells(combinedState).length,
        matrix.commands.length,
        tolerance
      ),
      locality: {
        alphaLocalityRadius: characterization.alphaDieOff,
        zLocalityRadius: characterization.zDieOff,
        alphaDecayRatio: characterization.alphaDecayRatio,
        zDecayRatio: characterization.zDecayRatio,
        alphaDecayLength: characterization.alphaDecayLength,
        zDecayLength: characterization.zDecayLength,
        decayProfile: {
          model: characterization.decayModel,
          alphaShells: characterization.alphaDecayShells,
          zShells: characterization.zDecayShells,
          alphaReach: characterization.alphaDecayReach,
          zReach: characterization.zDecayReach,
          alphaFirst: characterization.alphaDecayFirst,
          zFirst: characterization.zDecayFirst,
          alphaLast: characterization.alphaDecayLast,
          zLast: characterization.zDecayLast,
        },
      },
      reachability: {
        reachableAlphaCells: characterization.reachableAlphaCells,
        reachableHeightCells: characterization.reachableHeightCells,
        alphaUnderactuatedCells: characterization.alphaUnderactuatedCells,
        heightUnderactuatedCells: characterization.heightUnderactuatedCells,
        alphaRank: characterization.responseRankAlpha,
        heightRank: characterization.responseRankHeight,
      },
      composition: {
        nonadditive: characterization.nonadditive,
        alphaSuperpositionError: characterization.alphaSuperpositionError,
        heightSuperpositionError: characterization.heightSuperpositionError,
        superpositionRmsError: characterization.superpositionError,
        orderSensitive: Boolean(order.sequenceOrder?.orderSensitive),
        noncommutingAdjacentPairs: Number(order.sequenceOrder?.noncommutingAdjacentPairs) || 0,
        maxOrderError: Number(order.sequenceOrder?.maxOrderError) || 0,
        nonadditivePairCount: characterization.pairwiseNonadditivePairs,
        maxPairwiseInteractionError: characterization.pairwiseMaxInteractionError,
        maxPairwiseHotspotError: characterization.pairwiseInteractionMapMax,
        maxPairwiseInteractionDegree: characterization.pairwiseInteractionDegreeMax,
        pairwiseInteractionDensity: characterization.pairwiseInteractionDensity,
        pairwiseInteractionsTruncated: characterization.pairwiseTruncated,
      },
      combinedResponse: {
        scope,
        selected: reportCell(selected),
        cells: sourceCells.map(reportCell),
        activeSources: characterization.activeSources,
        responseCells: characterization.responseCells,
        alphaReach: characterization.alphaReachCells,
        zReach: characterization.zReachCells,
        effectiveAlphaDieOff: characterization.alphaDieOff,
        effectiveZDieOff: characterization.zDieOff,
        maxAbsAlphaDelta: characterization.maxAlphaDelta,
        maxAbsHeightDelta: characterization.maxHeightDelta,
        meanAbsAlphaDelta: characterization.meanAbsAlphaDelta,
        meanAbsHeightDelta: characterization.meanAbsHeightDelta,
        positiveZReachCells: characterization.positiveZReachCells,
        negativeZReachCells: characterization.negativeZReachCells,
        maxPositiveHeightDelta: characterization.maxPositiveHeightDelta,
        maxNegativeHeightDelta: characterization.maxNegativeHeightDelta,
      },
      responseMatrix,
      pairwiseInteractions: {
        model: characterization.pairwiseInteractionModel,
        commandCount: activeOperatorCount,
        totalPairCount: characterization.pairwiseTotalPairs,
        evaluatedPairCount: characterization.pairwiseEvaluatedPairs,
        nonadditivePairCount: characterization.pairwiseNonadditivePairs,
        maxAlphaError: characterization.pairwiseMaxAlphaError,
        maxHeightError: characterization.pairwiseMaxHeightError,
        maxInteractionError: characterization.pairwiseMaxInteractionError,
        maxHotspotError: characterization.pairwiseInteractionMapMax,
        maxInteractionDegree: characterization.pairwiseInteractionDegreeMax,
        interactionDensity: characterization.pairwiseInteractionDensity,
        truncated: characterization.pairwiseTruncated,
        interactions: characterization.pairwiseInteractions,
      },
      physicalValidation: reportPhysicalPreview(characterization),
      sequenceOrder: order.sequenceOrder,
      tolerance,
    };
    if (includeFields) {
      report.combinedResponse.fields = reportResponseFields(sim, baselineSim);
      report.pairwiseInteractions.fields = {
        alphaErrorMatrix: characterization.pairwiseAlphaErrorMatrix,
        heightErrorMatrix: characterization.pairwiseHeightErrorMatrix,
        interactionHotspotMap: characterization.pairwiseInteractionMap,
        interactionDegreeMap: characterization.pairwiseInteractionDegreeMap,
      };
    }
    return report;
  }

  function exportProgrammableDiscontinuityReport(state, options = {}) {
    return JSON.stringify(programmableDiscontinuityReport(state, options), null, 2);
  }

  function formalizationTargetManifest(state, options = {}) {
    return programmableDiscontinuityReport(state, {
      ...options,
      includeResponseMatrix: false,
      includeFields: false,
    }).formalizationTargets;
  }

  function exportFormalizationTargetManifest(state, options = {}) {
    return JSON.stringify(formalizationTargetManifest(state, options), null, 2);
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

  function springHingePhysicalPreviewReportDescriptor() {
    return {
      schema: "rad-sim.vertical-load-physical-preview-report.v1",
      comparisonSchema: "rad-sim.vertical-removal-physical-comparison.v1",
      mechanicsCertificateSchema: "rad-sim.mechanics-energy-certificate.v1",
      energyValidationSchema: "rad-sim.vertical-load-energy-validation.v1",
      energyMeasurementSchema: "rad-sim.vertical-load-energy-measurement-results.v1",
      energyProtocolSchema: "rad-sim.vertical-load-energy-experiment-protocol.v1",
      benchPacketSchema: "rad-sim.vertical-load-bench-packet.v1",
      energyComparisonSchema: "rad-sim.vertical-load-energy-comparison-report.v1",
      browserRunsSolver: false,
      pythonOnly: true,
      pythonFunctions: [
        "compare_vertical_residual_spring_hinge_3d",
        "mechanics_energy_certificate_to_dict",
        "vertical_load_energy_measurement_template",
        "vertical_load_energy_experiment_protocol",
        "export_vertical_load_energy_experiment_protocol_json",
        "vertical_load_bench_packet",
        "export_vertical_load_bench_packet_json",
        "vertical_load_energy_measurement_results_from_json",
        "compare_vertical_load_energy_measurement_results",
        "export_vertical_load_energy_measurement_template_json",
        "export_vertical_load_energy_comparison_report_json",
        "validate_vertical_load_energy_measurements",
        "export_vertical_load_energy_validation_json",
        "build_vertical_load_physical_preview_report",
        "export_vertical_load_physical_preview_report_json",
        "export_vertical_load_physical_preview_report_csv",
      ],
      requiredSolver: "solve_spring_hinge_3d",
      reportIncludes: [
        "solver convergence",
        "unanchored load components",
        "spring and hinge graph deletion counts",
        "kinematic-vs-physical height errors",
        "nonnegative mechanics energy certificate",
        "measured-vs-simulated load-work validation",
        "fillable vertical-load bench measurement template",
        "step-by-step vertical-load bench protocol",
        "combined vertical-load bench packet",
        "bench-table CSV rows",
        "claim labels",
      ],
      claimLimit: "Browser records the schema/descriptor only; SciPy physical-preview reports are generated by Python.",
    };
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
      alphaSuperpositionError: interaction.alphaMax,
      heightSuperpositionError: interaction.heightMax,
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
  const REACHABLE_EQUILIBRIUM_BENCH_FIELDS = Object.freeze([
    ...PROTOCOL_MEASUREMENT_FIELDS,
    "baseline_equilibrium_residual",
    "response_column_id",
    "target_reachability_flag",
    "topology_component_label",
    "topology_blocked_flag",
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

  function finiteOrNull(value) {
    const numeric = Number(value);
    return Number.isFinite(numeric) ? numeric : null;
  }

  function atlasObservationCells(sim, baseline, step) {
    return (step.observationCells || []).map((cell) => {
      const row = Number(cell.row) || 0;
      const col = Number(cell.col) || 0;
      return {
        row,
        col,
        alphaDelta: (Number(sim.alpha?.[row]?.[col]) || 0) - (Number(baseline.alpha?.[row]?.[col]) || 0),
        heightDelta: (Number(sim.height?.[row]?.[col]) || 0) - (Number(baseline.height?.[row]?.[col]) || 0),
        zResidual: Number(sim.zResidual?.[row]?.[col]) || 0,
        alphaDieOff: finiteOrNull(sim.dieOff?.[row]?.[col]),
        zDieOff: finiteOrNull(sim.zDieOff?.[row]?.[col]),
      };
    });
  }

  function atlasScopeCounts(entries) {
    return entries.reduce((counts, entry) => {
      counts[entry.scope] = (counts[entry.scope] || 0) + 1;
      return counts;
    }, {});
  }

  function meanOrZero(values) {
    const numeric = values.map(Number).filter(Number.isFinite);
    if (!numeric.length) return 0;
    return numeric.reduce((sum, value) => sum + value, 0) / numeric.length;
  }

  function uniqueParameterValues(values, fallback) {
    const raw = Array.isArray(values) && values.length ? values : fallback;
    const seen = new Set();
    const out = [];
    for (const value of raw) {
      const numeric = Math.max(0, Number(value) || 0);
      const key = numeric.toPrecision(12);
      if (seen.has(key)) continue;
      seen.add(key);
      out.push(numeric);
    }
    return out.length ? out : [0];
  }

  function defaultSweepValues(value, fallbackStep = 0.05) {
    const current = Math.max(0, Number(value) || 0);
    if (current <= 1e-9) return [0, fallbackStep, 2 * fallbackStep];
    return [0.5 * current, current, 2 * current];
  }

  function stateWithSweepSettings(state, backlash, clearance) {
    const temp = cloneForCharacterization(state);
    const pinRadius = Math.max(0, Number(temp.grid.pinRadius ?? 0.18));
    const pinHoleClearance = Math.max(0, Number(clearance) || 0);
    temp.grid.backlash = Math.max(0, Number(backlash) || 0);
    temp.grid.pinRadius = pinRadius;
    temp.grid.holeRadius = Math.max(pinRadius, pinRadius + pinHoleClearance);
    return temp;
  }

  function atlasSampleSummary(atlas) {
    const entries = atlas.entries || [];
    const observedZResiduals = [];
    const neighborZResiduals = [];
    for (const entry of entries) {
      const commandCells = new Set((entry.commands || []).map((command) => `${command.row},${command.col}`));
      for (const cell of entry.observationCells || []) {
        const residual = Math.abs(Number(cell.zResidual) || 0);
        observedZResiduals.push(residual);
        if (!commandCells.has(`${cell.row},${cell.col}`)) neighborZResiduals.push(residual);
      }
    }
    const physicalEntries = entries.filter((entry) => entry.physicalValidation?.physicalPreviewAvailable);
    const physicalHeightErrors = physicalEntries.map((entry) => entry.physicalValidation?.heightRmsModelError);
    return {
      meanAlphaReach: meanOrZero(entries.map((entry) => entry.alphaReach)),
      meanZReach: meanOrZero(entries.map((entry) => entry.zReach)),
      maxAlphaReach: entries.reduce((max, entry) => Math.max(max, Number(entry.alphaReach) || 0), 0),
      maxZReach: entries.reduce((max, entry) => Math.max(max, Number(entry.zReach) || 0), 0),
      maxSuperpositionError: entries.reduce(
        (max, entry) => Math.max(max, Number(entry.alphaSuperpositionError) || 0, Number(entry.heightSuperpositionError) || 0),
        0
      ),
      maxObservedZResidual: observedZResiduals.reduce((max, value) => Math.max(max, value), 0),
      maxObservedNeighborZResidual: neighborZResiduals.reduce((max, value) => Math.max(max, value), 0),
      meanAbsAlphaDelta: meanOrZero(entries.map((entry) => entry.meanAbsAlphaDelta)),
      meanAbsHeightDelta: meanOrZero(entries.map((entry) => entry.meanAbsHeightDelta)),
      physicalSuccessRate: physicalEntries.length
        ? physicalEntries.filter((entry) => entry.physicalValidation?.physicalSuccess).length / physicalEntries.length
        : null,
      meanPhysicalHeightRmsError: physicalEntries.length ? meanOrZero(physicalHeightErrors) : null,
    };
  }

  function sweepTrend(samples, parameter) {
    const buckets = new Map();
    for (const sample of samples) {
      const value = Number(sample.settings?.[parameter]) || 0;
      const key = value.toPrecision(12);
      if (!buckets.has(key)) buckets.set(key, { value, samples: [] });
      buckets.get(key).samples.push(sample);
    }
    return Array.from(buckets.values())
      .sort((a, b) => a.value - b.value)
      .map((bucket) => ({
        value: bucket.value,
        sampleCount: bucket.samples.length,
        meanAlphaReach: meanOrZero(bucket.samples.map((sample) => sample.summary.meanAlphaReach)),
        meanZReach: meanOrZero(bucket.samples.map((sample) => sample.summary.meanZReach)),
        maxObservedZResidual: bucket.samples.reduce((max, sample) => Math.max(max, sample.summary.maxObservedZResidual || 0), 0),
        maxObservedNeighborZResidual: bucket.samples.reduce((max, sample) => Math.max(max, sample.summary.maxObservedNeighborZResidual || 0), 0),
        maxSuperpositionError: bucket.samples.reduce((max, sample) => Math.max(max, sample.summary.maxSuperpositionError || 0), 0),
      }));
  }

  function endpointSlope(trend, metric) {
    if (!Array.isArray(trend) || trend.length < 2) return null;
    const first = trend[0];
    const last = trend[trend.length - 1];
    const dx = (Number(last.value) || 0) - (Number(first.value) || 0);
    if (Math.abs(dx) <= 1e-12) return null;
    return ((Number(last[metric]) || 0) - (Number(first[metric]) || 0)) / dx;
  }

  function sweepSensitivityFromTrends(trends) {
    const trendMap = {
      backlash: trends.byBacklash || [],
      pinHoleClearance: trends.byPinHoleClearance || [],
    };
    const metrics = {};
    let dominant = null;
    for (const [parameter, trend] of Object.entries(trendMap)) {
      metrics[parameter] = {};
      for (const metric of SWEEP_SENSITIVITY_METRICS) {
        const slope = endpointSlope(trend, metric);
        metrics[parameter][metric] = slope;
        if (slope === null) continue;
        const candidate = { parameter, metric, slope, absSlope: Math.abs(slope) };
        if (!dominant || candidate.absSlope > dominant.absSlope) dominant = candidate;
      }
    }
    return {
      method: "endpoint finite difference over each parameter trend",
      metrics,
      dominant,
    };
  }

  function trendMonotonicity(trend, metric, tolerance = 1e-12) {
    if (!Array.isArray(trend) || trend.length < 2) return "insufficient";
    let increases = 0;
    let decreases = 0;
    let flats = 0;
    for (let i = 1; i < trend.length; i += 1) {
      const delta = (Number(trend[i][metric]) || 0) - (Number(trend[i - 1][metric]) || 0);
      if (Math.abs(delta) <= tolerance) flats += 1;
      else if (delta > 0) increases += 1;
      else decreases += 1;
    }
    if (increases && !decreases) return "increasing";
    if (decreases && !increases) return "decreasing";
    if (flats && !increases && !decreases) return "flat";
    return "mixed";
  }

  function operatorLawStatement(parameter, metric, monotonicity) {
    const parameterLabels = {
      backlash: "backlash dead-zone width",
      pinHoleClearance: "pin-hole clearance",
    };
    const metricLabels = {
      meanAlphaReach: "mean dilation reach",
      meanZReach: "mean vertical reach",
      maxObservedNeighborZResidual: "neighbor vertical residual motion",
      maxSuperpositionError: "non-additive superposition residual",
    };
    const parameterLabel = parameterLabels[parameter] || parameter;
    const metricLabel = metricLabels[metric] || metric;
    if (monotonicity === "increasing") {
      return `Increasing ${parameterLabel} increases ${metricLabel} over the sampled simulator sweep.`;
    }
    if (monotonicity === "decreasing") {
      return `Increasing ${parameterLabel} decreases ${metricLabel} over the sampled simulator sweep.`;
    }
    if (monotonicity === "flat") {
      return `Changing ${parameterLabel} leaves ${metricLabel} approximately flat over the sampled simulator sweep.`;
    }
    if (monotonicity === "mixed") {
      return `${parameterLabel} has a mixed sampled relationship with ${metricLabel}; no monotone candidate is supported.`;
    }
    return `${parameterLabel} has insufficient sampled data to propose a ${metricLabel} law candidate.`;
  }

  function sweepOperatorLawCandidates(trends, sensitivity) {
    const trendMap = {
      backlash: trends.byBacklash || [],
      pinHoleClearance: trends.byPinHoleClearance || [],
    };
    const laws = [];
    for (const [parameter, trend] of Object.entries(trendMap)) {
      for (const metric of SWEEP_SENSITIVITY_METRICS) {
        const monotonicity = trendMonotonicity(trend, metric);
        laws.push({
          parameter,
          metric,
          monotonicity,
          slope: sensitivity?.metrics?.[parameter]?.[metric] ?? null,
          supportedBySweep: ["increasing", "decreasing", "flat"].includes(monotonicity),
          status: "simulator-diagnostic",
          statement: operatorLawStatement(parameter, metric, monotonicity),
        });
      }
    }
    return {
      schema: "rad-sim.operator-law-candidates.v1",
      method: "adjacent monotonicity over sampled parameter trend plus endpoint sensitivity",
      laws,
    };
  }

  function atlasStepEntry(state, step, tolerance = 1e-9) {
    const commandState = protocolStepState(state, step, true);
    const baselineState = protocolStepState(state, step, false);
    const sim = RAD.simulate(commandState);
    const baseline = RAD.simulate(baselineState);
    const sourceCells = uniqueProtocolCells(
      (step.commands || []).map((command) => normalizeMatrixCell(state, { row: command.row, col: command.col }))
    );
    const stats = localResponseStats(commandState, sim, baseline, sourceCells);
    const interaction = superpositionError(commandState, sourceCells, sim, baseline);
    const physical = physicalPreviewComparison(commandState, baselineState, sim, baseline);
    return {
      stepId: step.id,
      scope: step.scope,
      purpose: step.purpose,
      expectedResponse: step.expectedResponse,
      commands: step.commands || [],
      lockedCells: step.lockedCells || [],
      observationCells: atlasObservationCells(sim, baseline, step),
      alphaReach: stats.alphaReachCells,
      zReach: stats.zReachCells,
      effectiveAlphaDieOff: stats.alphaDieOff,
      effectiveZDieOff: stats.zDieOff,
      maxAbsAlphaDelta: stats.maxAlphaDelta,
      maxAbsHeightDelta: stats.maxHeightDelta,
      meanAbsAlphaDelta: stats.meanAbsAlphaDelta,
      meanAbsHeightDelta: stats.meanAbsHeightDelta,
      alphaSuperpositionError: interaction.alphaMax,
      heightSuperpositionError: interaction.heightMax,
      superpositionRmsError: interaction.rms,
      physicalValidation: {
        schema: "rad-sim.browser-physical-preview.v1",
        model: "browser-spring-preview",
        physicalPreviewAvailable: Boolean(physical.physicalPreviewAvailable),
        physicalSuccess: Boolean(physical.physicalPreviewSuccess),
        heightRmsModelError: Number(physical.physicalHeightRmsError) || 0,
        centerRmsModelError: Number(physical.physicalCenterRmsError) || 0,
        heightMaxModelError: Number(physical.physicalHeightMaxError) || 0,
        centerMaxModelError: Number(physical.physicalCenterMaxError) || 0,
        iterations: Number(physical.physicalPreviewIterations) || 0,
        error: physical.physicalPreviewError || null,
      },
      tolerance,
    };
  }

  function responseAtlas(state, options = {}) {
    const tolerance = Number(options.tolerance ?? 1e-9);
    const protocol = options.protocol || calibrationExperimentProtocol(state, options);
    const entries = (protocol.steps || []).map((step) => atlasStepEntry(state, step, tolerance));
    const maxSuperpositionError = entries.reduce(
      (max, entry) => Math.max(max, entry.alphaSuperpositionError || 0, entry.heightSuperpositionError || 0),
      0
    );
    const physicalEntries = entries.filter((entry) => entry.physicalValidation?.physicalPreviewAvailable);
    return {
      schema: "rad-sim.response-atlas.v1",
      savedAt: new Date().toISOString(),
      hardwareProfile: protocol.hardwareProfile,
      grid: protocol.grid,
      centerCell: protocol.centerCell,
      notes: "Browser-generated single, pair, cluster, and lock response atlas for comparing simulator settings before bench data exists.",
      summary: {
        entryCount: entries.length,
        scopeCounts: atlasScopeCounts(entries),
        maxAlphaReach: entries.reduce((max, entry) => Math.max(max, entry.alphaReach || 0), 0),
        maxZReach: entries.reduce((max, entry) => Math.max(max, entry.zReach || 0), 0),
        maxSuperpositionError,
        physicalEntryCount: physicalEntries.length,
        physicalSuccessCount: physicalEntries.filter((entry) => entry.physicalValidation?.physicalSuccess).length,
      },
      assumptions: {
        source: "Calibration protocol commands run through the browser simulator.",
        physical: "Browser spring-preview validation is interactive model-disagreement evidence, not calibrated hardware truth.",
      },
      protocol,
      entries,
      tolerance,
    };
  }

  function exportResponseAtlas(state, options = {}) {
    return JSON.stringify(responseAtlas(state, options), null, 2);
  }

  function responseAtlasSweep(state, options = {}) {
    const tolerance = Number(options.tolerance ?? 1e-9);
    const currentClearance = typeof RAD.pinHoleClearance === "function"
      ? RAD.pinHoleClearance(state)
      : Math.max(0, Number(state.grid.holeRadius ?? 0.225) - Number(state.grid.pinRadius ?? 0.18));
    const backlashValues = uniqueParameterValues(
      options.backlashValues,
      defaultSweepValues(state.grid.backlash, 0.05)
    );
    const clearanceValues = uniqueParameterValues(
      options.clearanceValues,
      defaultSweepValues(currentClearance, 0.04)
    );
    const samples = [];
    for (const backlash of backlashValues) {
      for (const clearance of clearanceValues) {
        const sampleState = stateWithSweepSettings(state, backlash, clearance);
        const atlas = responseAtlas(sampleState, { ...options, tolerance });
        samples.push({
          settings: {
            backlash: sampleState.grid.backlash,
            pinHoleClearance: Math.max(0, Number(sampleState.grid.holeRadius) - Number(sampleState.grid.pinRadius)),
            pinRadius: sampleState.grid.pinRadius,
            holeRadius: sampleState.grid.holeRadius,
          },
          summary: atlasSampleSummary(atlas),
          atlas,
        });
      }
    }
    const trends = {
      byBacklash: sweepTrend(samples, "backlash"),
      byPinHoleClearance: sweepTrend(samples, "pinHoleClearance"),
    };
    const sensitivity = sweepSensitivityFromTrends(trends);
    return {
      schema: "rad-sim.response-atlas-sweep.v1",
      savedAt: new Date().toISOString(),
      grid: { rows: state.grid.rows, cols: state.grid.cols },
      centerCell: {
        row: clampIndex(options.r ?? state.selection?.r ?? 0, state.grid.rows),
        col: clampIndex(options.c ?? state.selection?.c ?? 0, state.grid.cols),
      },
      physical: samples.some((sample) => sample.summary.physicalSuccessRate !== null),
      parameters: {
        backlashValues,
        pinHoleClearanceValues: clearanceValues,
      },
      notes: "Browser-generated sweep for comparing how backlash and pin-hole clearance change locality, residual vertical motion, and operator interaction before calibrated bench data exists.",
      summary: {
        sampleCount: samples.length,
        maxAlphaReach: samples.reduce((max, sample) => Math.max(max, sample.summary.maxAlphaReach || 0), 0),
        maxZReach: samples.reduce((max, sample) => Math.max(max, sample.summary.maxZReach || 0), 0),
        maxObservedZResidual: samples.reduce((max, sample) => Math.max(max, sample.summary.maxObservedZResidual || 0), 0),
        maxObservedNeighborZResidual: samples.reduce((max, sample) => Math.max(max, sample.summary.maxObservedNeighborZResidual || 0), 0),
        maxSuperpositionError: samples.reduce((max, sample) => Math.max(max, sample.summary.maxSuperpositionError || 0), 0),
      },
      trends,
      sensitivity,
      operatorLawCandidates: sweepOperatorLawCandidates(trends, sensitivity),
      assumptions: {
        source: "Each browser sample rebuilds the response atlas from the current calibration protocol commands.",
        interpretation: "Backlash and clearance are treated as programmable-discontinuity dead-zone parameters; trend summaries are simulator diagnostics.",
        physical: "Browser spring-preview metrics are model-disagreement diagnostics, not calibrated hardware validation.",
      },
      samples,
      tolerance,
    };
  }

  function exportResponseAtlasSweep(state, options = {}) {
    return JSON.stringify(responseAtlasSweep(state, options), null, 2);
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
      provenance: {
        schema: "rad-sim.calibration-dataset-provenance.v1",
        datasetId: options.datasetId || null,
        datasetRole: options.datasetRole || "unassigned",
        sourceFileId: options.sourceFileId || null,
        collectedAt: options.collectedAt || null,
        operator: options.operator || null,
        profileId: options.profileId || null,
        profileFrozenAt: options.profileFrozenAt || null,
        notes: options.provenanceNotes || "Set datasetRole to fit or holdout before using this file for independent validation.",
      },
      notes: "Fill optional measured fields with real bench measurements.",
      steps,
    };
  }

  function exportCalibrationExperimentResultsTemplate(state, options = {}) {
    return JSON.stringify(calibrationExperimentResultsTemplate(state, options), null, 2);
  }

  function calibrationBenchNotebook(state, options = {}) {
    const protocol = options.protocol || calibrationExperimentProtocol(state, options);
    const fitDatasetId = options.fitDatasetId || "fit-run-001";
    const holdoutDatasetId = options.holdoutDatasetId || "holdout-run-001";
    const profileId = options.profileId || "calibration-profile-v1";
    const profileFrozenAt = options.profileFrozenAt || "set-after-fit-before-holdout";
    const fitFilename = `${fitDatasetId}.json`;
    const holdoutFilename = `${holdoutDatasetId}.json`;
    const datasetPlan = {
      fit: {
        datasetId: fitDatasetId,
        role: "fit",
        suggestedFilename: fitFilename,
        sourceFileId: fitFilename,
        profileId: null,
        profileFrozenAt: null,
        provenanceSchema: "rad-sim.calibration-dataset-provenance.v1",
        purpose: "Estimate diagnostic response parameters and freeze a bounded model profile.",
      },
      holdout: {
        datasetId: holdoutDatasetId,
        role: "holdout",
        suggestedFilename: holdoutFilename,
        sourceFileId: holdoutFilename,
        profileId,
        profileFrozenAt,
        provenanceSchema: "rad-sim.calibration-dataset-provenance.v1",
        purpose: "Replay the frozen profile against data collected after the fit file is closed.",
      },
    };
    const scenarios = (protocol.steps || []).map((step) => ({
      stepId: step.id,
      scope: step.scope,
      repeatCount: step.repeatCount,
      commandCells: (step.commands || []).map((command) => ({ row: command.row, col: command.col })),
      commands: (step.commands || []).map((command) => ({
        row: command.row,
        col: command.col,
        alpha: Number(command.alpha) || 0,
        z: Number(command.z) || 0,
      })),
      lockedCells: (step.lockedCells || []).map((cell) => ({ row: cell.row, col: cell.col })),
      observationCells: (step.observationCells || []).map((cell) => ({ row: cell.row, col: cell.col })),
      measurementFields: [...(step.measurementFields || PROTOCOL_MEASUREMENT_FIELDS)],
      datasetRoles: ["fit", "holdout"],
      purpose: step.purpose,
      expectedResponse: step.expectedResponse,
      claimLabel: "bench protocol step; no physical law claimed until measured",
    }));
    return {
      schema: "rad-sim.calibration-bench-notebook.v1",
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        cellSize: state.grid.cellSize,
        backlash: state.grid.backlash,
        couplingGain: state.grid.couplingGain,
        zCouplingGain: state.grid.zCouplingGain,
        pinRadius: state.grid.pinRadius,
        holeRadius: state.grid.holeRadius,
        pinHoleClearance: typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : Math.max(0, (state.grid.holeRadius || 0) - (state.grid.pinRadius || 0)),
      },
      protocol: {
        schema: protocol.schema,
        hardwareProfile: protocol.hardwareProfile,
        centerCell: protocol.centerCell,
        stepCount: (protocol.steps || []).length,
        repeatCountTotal: (protocol.steps || []).reduce((sum, step) => sum + (Number(step.repeatCount) || 0), 0),
        notes: protocol.notes,
      },
      datasetPlan,
      instruments: [
        { id: "calipers-or-micrometer", measures: ["pin radius", "hole radius", "plate thickness", "joint stack height"], claimLabel: "required geometry measurement" },
        { id: "motion-tracking-camera", measures: ["cell center displacement", "height delta", "neighbor residual z"], claimLabel: "required response measurement" },
        { id: "force-gauge-or-load-cell", measures: ["actuator force", "contact/load proxy"], claimLabel: "required force measurement for contact calibration" },
        { id: "fixture-and-clamps", measures: ["locked cell enforcement", "fixed boundary repeatability"], claimLabel: "required boundary-condition control" },
        { id: "actuator-controller", measures: ["commanded alpha", "commanded z", "command repeat timing"], claimLabel: "required command provenance" },
        { id: "scale-reference-marker", measures: ["pixel-to-mm scale", "coordinate registration"], claimLabel: "required unit-scale evidence" },
      ],
      measurementColumns: [
        "datasetId",
        "datasetRole",
        "sourceFileId",
        "collectedAt",
        "operator",
        "profileId",
        "profileFrozenAt",
        ...PROTOCOL_MEASUREMENT_FIELDS,
      ],
      phases: [
        { id: "setup", description: "Measure hardware profile fields and register the cell coordinate frame.", requiredEvidence: ["hardware profile", "scale marker", "fixture notes"] },
        { id: "fit-collection", description: "Collect the fit dataset without using any holdout measurements.", requiredEvidence: [fitDatasetId, "datasetRole=fit"] },
        { id: "profile-freeze", description: "Generate and freeze a model profile before holdout collection starts.", requiredEvidence: [profileId, profileFrozenAt] },
        { id: "holdout-collection", description: "Collect a separate holdout file after the model profile is frozen.", requiredEvidence: [holdoutDatasetId, "datasetRole=holdout"] },
        {
          id: "validation",
          description: "Run residual validation and independent split checks before making physical claims.",
          requiredEvidence: ["rad-sim.calibration-model-profile-holdout-validation.v1", "rad-sim.calibration-train-holdout-split.v1"],
        },
      ],
      scenarios,
      outputs: [
        { id: "protocol-json", schema: protocol.schema, browserFunction: "exportCalibrationExperimentProtocol", browserAction: "Save Protocol" },
        { id: "fit-results-template-json", schema: "rad-sim.calibration-experiment-results.v1", datasetRole: "fit", browserFunction: "exportCalibrationExperimentResultsTemplate", browserAction: "Save Results Template" },
        { id: "holdout-results-template-json", schema: "rad-sim.calibration-experiment-results.v1", datasetRole: "holdout", browserFunction: "exportCalibrationExperimentResultsTemplate", browserAction: "Save Results Template" },
        { id: "model-profile-json", schema: "rad-sim.calibration-model-profile.v1", browserFunction: "exportCalibrationModelProfile", browserAction: "Save Model Profile" },
        { id: "holdout-validation-json", schema: "rad-sim.calibration-model-profile-holdout-validation.v1", browserFunction: "exportCalibrationModelProfileHoldoutValidation", browserAction: "Save Holdout Check" },
        { id: "holdout-validation-csv", schema: "rad-sim.calibration-model-profile-holdout-validation.v1", browserFunction: "exportCalibrationModelProfileHoldoutValidationCsv", browserAction: "Save Holdout CSV" },
        { id: "bench-notebook-csv", schema: "rad-sim.calibration-bench-notebook.v1", browserFunction: "exportCalibrationBenchNotebookCsv", browserAction: "Save Bench CSV" },
      ],
      passFailCriteria: {
        residualValidationPass: "fit and holdout residual scores do not increase after an applied profile update",
        independentValidationPass: "residual pass plus complete train/holdout provenance metadata",
        missingEvidence: "must be empty before claiming independent hardware validation",
        physicalClaimLimit: "passing this protocol supports calibration bookkeeping only; contact, friction, stiffness, and material laws still require model-specific validation",
      },
      formalization: {
        targetId: "calibration_bench_protocol_coverage",
        leanStructure: "Mechanics.CalibrationBenchProtocolCoverageNat",
        leanPredicate: "calibrationBenchProtocolCoverageReadyNat",
        schema: "rad-sim.calibration-bench-notebook.v1",
      },
      claimLabels: {
        notebook: "bench protocol artifact",
        scenarioRows: "experimentally unvalidated physical procedure until performed",
        formalCoverage: "Lean-proven finite coverage predicate, not physical accuracy",
      },
      limitations: [
        "The notebook enforces artifact coverage and provenance fields but cannot prove the lab actually collected independent files.",
        "The protocol is normalized and must be mapped to measured hardware units before physical claims.",
        "Rigid-body contact, friction, stiffness, gravity sag, and actuator-force laws remain uncalibrated.",
      ],
    };
  }

  function exportCalibrationBenchNotebook(state, options = {}) {
    return JSON.stringify(calibrationBenchNotebook(state, options), null, 2);
  }

  function exportCalibrationBenchNotebookCsv(notebookOrState, options = {}) {
    const notebook =
      notebookOrState?.schema === "rad-sim.calibration-bench-notebook.v1"
        ? notebookOrState
        : calibrationBenchNotebook(notebookOrState, options);
    const header = [
      "schema",
      "step_id",
      "scope",
      "repeat_count",
      "command_cells",
      "locked_cells",
      "observation_cells",
      "measurement_columns",
      "fit_dataset_id",
      "holdout_dataset_id",
      "purpose",
      "expected_response",
      "claim_label",
    ];
    const fit = notebook.datasetPlan?.fit || {};
    const holdout = notebook.datasetPlan?.holdout || {};
    const rows = [header];
    for (const scenario of notebook.scenarios || []) {
      rows.push([
        notebook.schema,
        scenario.stepId,
        scenario.scope,
        scenario.repeatCount,
        JSON.stringify(scenario.commandCells || []),
        JSON.stringify(scenario.lockedCells || []),
        JSON.stringify(scenario.observationCells || []),
        (scenario.measurementFields || []).join(";"),
        fit.datasetId || "",
        holdout.datasetId || "",
        scenario.purpose || "",
        scenario.expectedResponse || "",
        scenario.claimLabel || "",
      ]);
    }
    return rows.map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function calibrationBenchPacket(state, options = {}) {
    const protocol = options.protocol || calibrationExperimentProtocol(state, options);
    const fitDatasetId = options.fitDatasetId || "fit-run-001";
    const holdoutDatasetId = options.holdoutDatasetId || "holdout-run-001";
    const profileId = options.profileId || "calibration-profile-v1";
    const profileFrozenAt = options.profileFrozenAt || "set-after-fit-before-holdout";
    const notebook = calibrationBenchNotebook(state, {
      ...options,
      protocol,
      fitDatasetId,
      holdoutDatasetId,
      profileId,
      profileFrozenAt,
    });
    const fitTemplate = calibrationExperimentResultsTemplate(state, {
      protocol,
      datasetId: fitDatasetId,
      datasetRole: "fit",
      sourceFileId: `${fitDatasetId}.json`,
      profileId: null,
      profileFrozenAt: null,
      provenanceNotes: "Fit file: collect before generating and freezing a profile.",
    });
    const holdoutTemplate = calibrationExperimentResultsTemplate(state, {
      protocol,
      datasetId: holdoutDatasetId,
      datasetRole: "holdout",
      sourceFileId: `${holdoutDatasetId}.json`,
      profileId,
      profileFrozenAt,
      provenanceNotes: "Holdout file: collect only after the fit-derived profile is frozen.",
    });
    const filenames = {
      packet: "calibration_bench_packet.json",
      notebook: "calibration_bench_notebook.json",
      notebook_csv: "calibration_bench_notebook.csv",
      protocol: "calibration_experiment_protocol.json",
      fit_template: "calibration_fit_results_template.json",
      holdout_template: "calibration_holdout_results_template.json",
      readme: "README.md",
    };
    return {
      schema: "rad-sim.calibration-bench-packet.v1",
      method: "bench handoff packet bundling calibration protocol, fit and holdout measurement templates, notebook scenario table, and validation instructions",
      schemas: {
        notebook: "rad-sim.calibration-bench-notebook.v1",
        protocol: protocol.schema,
        resultsTemplate: "rad-sim.calibration-experiment-results.v1",
        datasetProvenance: "rad-sim.calibration-dataset-provenance.v1",
        modelProfile: "rad-sim.calibration-model-profile.v1",
        holdoutValidation: "rad-sim.calibration-model-profile-holdout-validation.v1",
        splitMetadata: "rad-sim.calibration-train-holdout-split.v1",
      },
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        cellSize: state.grid.cellSize,
        backlash: state.grid.backlash,
        couplingGain: state.grid.couplingGain,
        zCouplingGain: state.grid.zCouplingGain,
        pinRadius: state.grid.pinRadius,
        holeRadius: state.grid.holeRadius,
        pinHoleClearance: typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : Math.max(0, (state.grid.holeRadius || 0) - (state.grid.pinRadius || 0)),
      },
      datasetPlan: notebook.datasetPlan,
      filenames,
      benchNotebook: notebook,
      benchNotebookCsv: exportCalibrationBenchNotebookCsv(notebook),
      experimentProtocol: protocol,
      fitResultsTemplate: fitTemplate,
      holdoutResultsTemplate: holdoutTemplate,
      artifactManifest: [
        { id: "packet", filename: filenames.packet, schema: "rad-sim.calibration-bench-packet.v1", purpose: "single JSON bundle for review and archival" },
        { id: "notebook", filename: filenames.notebook, schema: "rad-sim.calibration-bench-notebook.v1", purpose: "human-readable protocol plan with phases and outputs" },
        { id: "notebook-csv", filename: filenames.notebook_csv, schema: "rad-sim.calibration-bench-notebook.v1", purpose: "one-row-per-scenario table for lab notebook or spreadsheet" },
        { id: "protocol", filename: filenames.protocol, schema: protocol.schema, purpose: "machine-readable protocol used by simulator comparisons" },
        { id: "fit-template", filename: filenames.fit_template, schema: "rad-sim.calibration-experiment-results.v1", purpose: "blank fit measurement rows" },
        { id: "holdout-template", filename: filenames.holdout_template, schema: "rad-sim.calibration-experiment-results.v1", purpose: "blank holdout measurement rows tied to frozen profile ID" },
      ],
      validationInstructions: {
        fitStep: "Fill the fit template first and use it to generate a model profile.",
        freezeStep: "Assign profileId and profileFrozenAt before collecting holdout measurements.",
        holdoutStep: "Fill the holdout template in a separate raw file after freezing the profile.",
        compareFunction: "calibrationModelProfileHoldoutValidation",
        csvFunction: "exportCalibrationModelProfileHoldoutValidationCsv",
        acceptanceRule: "Treat independentValidationPass as meaningful only when residual validation passes and splitMetadata.missingEvidence is empty.",
      },
      formalization: {
        targetId: "calibration_bench_packet_completeness",
        leanStructure: "Mechanics.CalibrationBenchPacketCompletenessNat",
        leanPredicate: "calibrationBenchPacketCompleteNat",
        schema: "rad-sim.calibration-bench-packet.v1",
      },
      claimLabels: {
        packetAssembly: "bench protocol artifact",
        fitTemplate: "experimentally unvalidated physical procedure until filled",
        holdoutTemplate: "experimentally unvalidated physical procedure until independently filled",
        formalCompleteness: "Lean-proven finite packet-completeness predicate, not physical accuracy",
      },
      limitations: [
        "The packet separates fit and holdout files but cannot prove laboratory independence by itself.",
        "The templates are blank until real bench measurements replace null values.",
        "Passing simulator residual checks is not a proof of contact, friction, stiffness, gravity, or material physics.",
      ],
    };
  }

  function exportCalibrationBenchPacket(state, options = {}) {
    return JSON.stringify(calibrationBenchPacket(state, options), null, 2);
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
      topCells: [],
    };
  }

  function finalizeCalibrationErrorField(field) {
    const topCells = [];
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
        if (field.sampleCount[r][c] > 0) {
          topCells.push({
            row: r,
            col: c,
            alphaError: alpha,
            heightError: height,
            combinedError: combined,
            sampleCount: field.sampleCount[r][c],
            alphaSampleCount: alphaCount,
            heightSampleCount: heightCount,
          });
        }
      }
    }
    field.topCells = topCells
      .sort((a, b) => b.combinedError - a.combinedError || a.row - b.row || a.col - b.col)
      .slice(0, 12);
    return field;
  }

  function calibrationFitResidualField(state, alphaPairs, heightPairs, fit) {
    const field = createCalibrationErrorField(state);
    const touched = new Map();
    const alphaGain = Number.isFinite(fit?.alpha?.suggestedGain) ? fit.alpha.suggestedGain : 1;
    const alphaBias = Number.isFinite(fit?.alpha?.suggestedBias) ? fit.alpha.suggestedBias : 0;
    const heightGain = Number.isFinite(fit?.height?.suggestedGain) ? fit.height.suggestedGain : 1;
    const heightBias = Number.isFinite(fit?.height?.suggestedBias) ? fit.height.suggestedBias : 0;
    for (const pair of alphaPairs) {
      const error = pair.measured - (alphaGain * pair.predicted + alphaBias);
      field.alphaError[pair.row][pair.col] += error;
      field.alphaSampleCount[pair.row][pair.col] += 1;
      touched.set(`${pair.row},${pair.col}`, [pair.row, pair.col]);
    }
    for (const pair of heightPairs) {
      const error = pair.measured - (heightGain * pair.predicted + heightBias);
      field.heightError[pair.row][pair.col] += error;
      field.heightSampleCount[pair.row][pair.col] += 1;
      touched.set(`${pair.row},${pair.col}`, [pair.row, pair.col]);
    }
    for (const [row, col] of touched.values()) field.sampleCount[row][col] += 1;
    return finalizeCalibrationErrorField(field);
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
          alphaPairs.push({ row, col, predicted: predictedAlpha, measured: alpha });
          field.alphaError[row][col] += alphaError;
          field.alphaSampleCount[row][col] += 1;
          measuredField = true;
        }
        if (height !== null) {
          const predictedHeight = (sim.height?.[row]?.[col] || 0) - (baseline.height?.[row]?.[col] || 0);
          const heightError = height - predictedHeight;
          heightErrors.push(heightError);
          heightPairs.push({ row, col, predicted: predictedHeight, measured: height });
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
    const fit = {
      alpha: calibrationLinearFit(alphaPairs),
      height: calibrationLinearFit(heightPairs),
    };
    return {
      schema: "rad-sim.calibration-experiment-comparison.v1",
      comparisons,
      field: finalizeCalibrationErrorField(field),
      fit,
      fitResidualField: calibrationFitResidualField(state, alphaPairs, heightPairs, fit),
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
    const fitResidualField = comparison?.fitResidualField || {};
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
      fitResidualMaxCombinedError: Number.isFinite(fitResidualField.maxCombinedError) ? fitResidualField.maxCombinedError : null,
      fitResidualWorstCell: fitResidualField.worstCell || null,
      topCells: Array.isArray(field.topCells) ? field.topCells : [],
      fitResidualTopCells: Array.isArray(fitResidualField.topCells) ? fitResidualField.topCells : [],
      maxAbsHeightError: worst ? Number(worst.maxAbsHeightError) : null,
      maxAbsAlphaError: Number.isFinite(field.maxAbsAlphaError) ? field.maxAbsAlphaError : null,
      maxCombinedError: Number.isFinite(field.maxCombinedError) ? field.maxCombinedError : null,
      worstCell: field.worstCell || null,
      worstStepId: worst?.stepId || null,
      meanActuatorForceN: finiteAverage(forceValues),
      meanPinHoleSlipMm: finiteAverage(slipValues),
    };
  }

  function calibrationParameterEstimates(state, results, comparison, summary = null, options = {}) {
    const parsed = typeof results === "string" ? JSON.parse(results) : results || {};
    const protocol = options.protocol || calibrationExperimentProtocol(state, options);
    const stepsById = new Map((protocol.steps || []).map((step) => [step.id, step]));
    const resolvedSummary = summary || summarizeCalibrationComparison(comparison || {});
    const fit = comparison?.fit || {};
    const zRatios = [];
    const slipValues = [];
    const forceValues = [];
    const forceHeightPairs = [];
    for (const measuredStep of parsed.steps || []) {
      const step = stepsById.get(measuredStep.stepId);
      if (!step) continue;
      const commandCells = new Set((step.commands || []).map((command) => `${command.row},${command.col}`));
      const directHeights = [];
      const neighborHeights = [];
      for (const cell of measuredStep.cells || []) {
        const height = numericOrNull(cell.heightDelta);
        const force = numericOrNull(cell.actuatorForceN);
        const slip = numericOrNull(cell.pinHoleSlipMm);
        if (slip !== null) slipValues.push(slip);
        if (force !== null) forceValues.push(Math.abs(force));
        if (height !== null) {
          const heightMagnitude = Math.abs(height);
          if (force !== null) forceHeightPairs.push({ height: heightMagnitude, force: Math.abs(force) });
          if (commandCells.has(`${cell.row},${cell.col}`)) directHeights.push(heightMagnitude);
          else neighborHeights.push(heightMagnitude);
        }
      }
      const direct = finiteAverage(directHeights);
      if (direct !== null && direct > 1e-12) {
        for (const neighbor of neighborHeights) zRatios.push(neighbor / direct);
      }
    }
    const forceHeightDenominator = forceHeightPairs.reduce((sum, pair) => sum + pair.height * pair.height, 0);
    const forceHeightSlope =
      forceHeightDenominator > 1e-12
        ? forceHeightPairs.reduce((sum, pair) => sum + pair.height * pair.force, 0) / forceHeightDenominator
        : null;
    const alphaFit = fit.alpha || {};
    const heightFit = fit.height || {};
    const zEstimate = finiteAverage(zRatios);
    return {
      schema: "rad-sim.calibration-parameter-estimates.v1",
      method:
        "conservative parameter bookkeeping from completed calibration experiment measurements; estimates are not treated as physical laws without separate validation",
      claimLabels: {
        alphaResponseGain: "simulator-derived empirical fit",
        heightResponseGain: "simulator-derived empirical fit",
        zCouplingGain: "experimentally unvalidated empirical estimate",
        verticalFreePlay: "experimentally unvalidated physical assumption",
        backlash: "not identifiable from current measurement set",
        contactProxy: "experimentally unvalidated physical assumption",
      },
      estimates: {
        alphaResponse: {
          suggestedGain: alphaFit.suggestedGain ?? null,
          suggestedBias: alphaFit.suggestedBias ?? null,
          sampleCount: Number(alphaFit.sampleCount) || 0,
          gainIdentifiable: Boolean(alphaFit.gainIdentifiable),
          rmsRawError: alphaFit.rmsRawError ?? null,
          rmsResidual: alphaFit.rmsResidual ?? null,
          status: Number(alphaFit.sampleCount) > 0 ? "estimated-from-results" : "missing-measurements",
        },
        heightResponse: {
          suggestedGain: heightFit.suggestedGain ?? null,
          suggestedBias: heightFit.suggestedBias ?? null,
          sampleCount: Number(heightFit.sampleCount) || 0,
          gainIdentifiable: Boolean(heightFit.gainIdentifiable),
          rmsRawError: heightFit.rmsRawError ?? null,
          rmsResidual: heightFit.rmsResidual ?? null,
          status: Number(heightFit.sampleCount) > 0 ? "estimated-from-results" : "missing-measurements",
        },
        zCouplingGain: {
          estimate: zEstimate,
          configured: state.grid.zCouplingGain,
          sampleCount: zRatios.length,
          method: "mean neighbor/direct measured height ratio on z-actuated protocol steps",
          status: zRatios.length ? "estimated-from-results" : "missing-pair-z-measurements",
        },
        verticalFreePlay: {
          configuredModelUnits: typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : Math.max(0, (state.grid.holeRadius || 0) - (state.grid.pinRadius || 0)),
          meanMeasuredSlipMm: finiteAverage(slipValues),
          sampleCount: slipValues.length,
          status: slipValues.length ? "measured-slip-summary" : "configured-clearance-only",
        },
        backlash: {
          configured: state.grid.backlash,
          estimate: null,
          sampleCount: 0,
          status: "not-identifiable-from-current-results",
          requiredExperiment: "one-cell alpha command sweep across the dead-zone with reversal/hysteresis measurements",
        },
        contactProxy: {
          meanActuatorForceN: finiteAverage(forceValues),
          forceHeightSlopeNPerModelUnit: forceHeightSlope,
          sampleCount: forceHeightPairs.length,
          status: forceHeightPairs.length ? "force-displacement-proxy" : "missing-force-measurements",
        },
      },
      quality: {
        comparisonCount: comparison?.comparisons?.length || 0,
        measuredCellCount: resolvedSummary.measuredCellCount ?? null,
        missingObservationCount: resolvedSummary.missingObservationCount ?? null,
        fitResidualMaxCombinedError: resolvedSummary.fitResidualMaxCombinedError ?? null,
        maxCombinedError: resolvedSummary.maxCombinedError ?? null,
      },
      suggestedModelUpdates: {
        zCouplingGain: zEstimate,
        alphaResponseGain: alphaFit.suggestedGain ?? null,
        heightResponseGain: heightFit.suggestedGain ?? null,
      },
      limitations: [
        "The z-coupling estimate is a measured response ratio, not a derivation of the simulator z-coupling law.",
        "Backlash is not identifiable without a command sweep across the dead-zone and reversal path.",
        "The contact proxy uses force/displacement summaries only and does not calibrate friction or rigid-body contact.",
      ],
    };
  }

  function profileFinite(value) {
    const numeric = Number(value);
    return Number.isFinite(numeric) ? numeric : null;
  }

  function profileUpdate({ name, configField, current, proposed, sampleCount, sourceStatus, claimLabel, safeToApply, reason }) {
    return {
      name,
      configField,
      current,
      proposed,
      delta: proposed !== null && current !== null ? proposed - current : null,
      sampleCount: Number(sampleCount) || 0,
      sourceStatus: sourceStatus || "unknown",
      claimLabel,
      safeToApply: Boolean(safeToApply),
      reason,
    };
  }

  function calibrationModelProfile(state, reportOrEstimates = null, options = {}) {
    const report = reportOrEstimates || (state.experiment?.calibrationComparison ? calibrationComparisonReport(state) : null);
    const parameterEstimates = report?.parameterEstimates || report;
    if (!parameterEstimates || parameterEstimates.schema !== "rad-sim.calibration-parameter-estimates.v1") {
      throw new Error("calibration model profile requires parameter estimates");
    }
    const minSamples = Math.max(1, Math.floor(Number(options.minSamples ?? 1)));
    const claims = parameterEstimates.claimLabels || {};
    const estimates = parameterEstimates.estimates || {};
    const quality = parameterEstimates.quality || {};
    const z = estimates.zCouplingGain || {};
    const zEstimate = profileFinite(z.estimate);
    const zSamples = Number(z.sampleCount) || 0;
    const zSafe = zEstimate !== null && zEstimate >= 0 && zEstimate <= 1 && zSamples >= minSamples && z.status === "estimated-from-results";
    const updates = [
      profileUpdate({
        name: "zCouplingGain",
        configField: "zCouplingGain",
        current: profileFinite(state.grid.zCouplingGain),
        proposed: zEstimate,
        sampleCount: zSamples,
        sourceStatus: z.status,
        claimLabel: claims.zCouplingGain || "experimentally unvalidated empirical estimate",
        safeToApply: zSafe,
        reason: zSafe
          ? "bounded measured neighbor/direct z-response ratio"
          : "requires a finite z-coupling estimate in [0, 1] with enough pair measurements",
      }),
    ];
    const diagnostics = [
      ["alphaResponseGain", "alphaResponse", "suggestedGain", claims.alphaResponseGain || "simulator-derived empirical fit"],
      ["heightResponseGain", "heightResponse", "suggestedGain", claims.heightResponseGain || "simulator-derived empirical fit"],
      ["verticalFreePlay", "verticalFreePlay", "meanMeasuredSlipMm", claims.verticalFreePlay || "experimentally unvalidated physical assumption"],
      ["backlash", "backlash", "estimate", claims.backlash || "not identifiable from current measurement set"],
      ["contactProxy", "contactProxy", "forceHeightSlopeNPerModelUnit", claims.contactProxy || "experimentally unvalidated physical assumption"],
    ];
    for (const [name, sectionName, valueField, claimLabel] of diagnostics) {
      const section = estimates[sectionName] || {};
      updates.push(
        profileUpdate({
          name,
          configField: null,
          current: null,
          proposed: profileFinite(section[valueField]),
          sampleCount: Number(section.sampleCount) || 0,
          sourceStatus: section.status,
          claimLabel,
          safeToApply: false,
          reason: "diagnostic in v1; no corresponding browser grid field is mutated",
        })
      );
    }
    const safeUpdates = updates.filter((update) => update.safeToApply);
    return {
      schema: "rad-sim.calibration-model-profile.v1",
      sourceReportSchema: options.sourceReportSchema || report?.schema || parameterEstimates.schema,
      method:
        "claim-labeled fitted-profile artifact; safe updates may adjust bounded simulator parameters, while non-configurable fitted response gains remain diagnostic evidence",
      configBefore: {
        backlash: state.grid.backlash,
        couplingGain: state.grid.couplingGain,
        zCouplingGain: state.grid.zCouplingGain,
        pinRadius: state.grid.pinRadius,
        holeRadius: state.grid.holeRadius,
        pinHoleClearance: typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : Math.max(0, (state.grid.holeRadius || 0) - (state.grid.pinRadius || 0)),
      },
      recommendedUpdates: updates,
      safeUpdateCount: safeUpdates.length,
      diagnosticOnlyCount: updates.length - safeUpdates.length,
      blockers: updates
        .filter((update) => !update.safeToApply)
        .map((update) => ({ name: update.name, reason: update.reason, sourceStatus: update.sourceStatus })),
      quality: {
        measuredCellCount: quality.measuredCellCount ?? null,
        missingObservationCount: quality.missingObservationCount ?? null,
        fitResidualMaxCombinedError: quality.fitResidualMaxCombinedError ?? null,
        maxCombinedError: quality.maxCombinedError ?? null,
        minSamples,
      },
      claimLabels: {
        profile: "simulator-derived empirical law",
        appliedConfigUpdates: "experimentally unvalidated empirical estimate",
        diagnosticOnlyUpdates: "not a simulator mutation in v1",
      },
      limitations: [
        "Applying this profile updates only bounded simulator parameters with direct grid fields.",
        "Alpha and height response gains are recorded for model review but do not change the kinematic law in v1.",
        "Backlash, pin-hole clearance, friction, and rigid-body contact remain separate calibration targets.",
      ],
    };
  }

  function applyCalibrationModelProfile(state, profile) {
    const resolvedProfile = profile || state.experiment?.calibrationModelProfile || calibrationModelProfile(state);
    if (resolvedProfile?.schema !== "rad-sim.calibration-model-profile.v1") {
      throw new Error("unsupported calibration model profile schema");
    }
    const applied = [];
    const skipped = [];
    for (const update of resolvedProfile.recommendedUpdates || []) {
      if (!update?.safeToApply || update.proposed === null || update.proposed === undefined) {
        skipped.push({ name: update?.name, configField: update?.configField, reason: update?.reason || "not marked safe to apply" });
        continue;
      }
      if (update.configField === "zCouplingGain") {
        const proposed = profileFinite(update.proposed);
        if (proposed === null || proposed < 0 || proposed > 1) {
          skipped.push({ name: update.name, configField: update.configField, reason: "zCouplingGain must stay in [0, 1]" });
          continue;
        }
        state.grid.zCouplingGain = proposed;
        applied.push({ name: update.name, configField: update.configField, value: proposed });
      } else {
        skipped.push({ name: update.name, configField: update.configField, reason: "no safe v1 grid mutation is defined for this field" });
      }
    }
    const audit = {
      schema: "rad-sim.calibration-model-profile-application.v1",
      sourceProfileSchema: resolvedProfile.schema,
      appliedUpdates: applied,
      skippedUpdates: skipped,
      configAfter: {
        zCouplingGain: state.grid.zCouplingGain,
      },
      claimLabel: "simulator configuration update from empirical calibration profile",
    };
    state.experiment.calibrationModelProfile = resolvedProfile;
    state.experiment.calibrationModelProfileApplication = audit;
    if (!Array.isArray(state.experiment.calibrationModelProfileHistory)) state.experiment.calibrationModelProfileHistory = [];
    state.experiment.calibrationModelProfileHistory.push({
      at: new Date().toISOString(),
      kind: "application",
      profileSchema: resolvedProfile.schema,
      appliedUpdates: applied,
      skippedUpdates: skipped,
      configAfter: audit.configAfter,
      claimLabel: audit.claimLabel,
    });
    if (state.experiment.calibrationModelProfileHistory.length > 50) {
      state.experiment.calibrationModelProfileHistory.splice(0, state.experiment.calibrationModelProfileHistory.length - 50);
    }
    return audit;
  }

  function exportCalibrationModelProfile(state) {
    return JSON.stringify(calibrationModelProfile(state), null, 2);
  }

  function cloneCalibrationState(state) {
    return RAD.deserialize(RAD.serialize(state));
  }

  function comparisonSummaryMetrics(reportOrSummary) {
    const summary = reportOrSummary?.summary || reportOrSummary || {};
    const metrics = {
      maxCombinedError: profileFinite(summary.maxCombinedError),
      fitResidualMaxCombinedError: profileFinite(summary.fitResidualMaxCombinedError),
      heightRmseMean: profileFinite(summary.heightRmseMean),
      alphaRmseMean: profileFinite(summary.alphaRmseMean),
      missingObservationCount: Number(summary.missingObservationCount) || 0,
    };
    const finiteResiduals = [
      metrics.maxCombinedError,
      metrics.fitResidualMaxCombinedError,
      metrics.heightRmseMean,
      metrics.alphaRmseMean,
    ].filter((value) => value !== null);
    metrics.residualScore = finiteResiduals.length ? finiteResiduals.reduce((sum, value) => sum + value, 0) : null;
    return metrics;
  }

  function residualDelta(before, after) {
    const delta = {};
    for (const key of ["maxCombinedError", "fitResidualMaxCombinedError", "heightRmseMean", "alphaRmseMean", "residualScore"]) {
      delta[key] = before[key] !== null && after[key] !== null ? after[key] - before[key] : null;
    }
    delta.missingObservationCount = (Number(after.missingObservationCount) || 0) - (Number(before.missingObservationCount) || 0);
    return delta;
  }

  function calibrationDatasetProvenance(results) {
    const provenance = results?.provenance || {};
    return {
      schema: "rad-sim.calibration-dataset-provenance.v1",
      datasetId: provenance.datasetId || results?.datasetId || null,
      datasetRole: provenance.datasetRole || results?.datasetRole || "unassigned",
      sourceFileId: provenance.sourceFileId || results?.sourceFileId || null,
      collectedAt: provenance.collectedAt || results?.collectedAt || null,
      operator: provenance.operator || results?.operator || null,
      profileId: provenance.profileId || results?.profileId || null,
      profileFrozenAt: provenance.profileFrozenAt || results?.profileFrozenAt || null,
      notes: provenance.notes || "",
    };
  }

  function calibrationProfileIdentifier(profile) {
    if (profile?.profileId) return String(profile.profileId);
    if (profile?.sourceReportSchema || profile?.safeUpdateCount !== undefined) {
      return `${profile.sourceReportSchema || "profile"}:safe=${profile.safeUpdateCount}`;
    }
    return null;
  }

  function calibrationSplitMetadata(fitResults, holdoutResults, profile, fitComparison, holdoutComparison) {
    const fitDataset = calibrationDatasetProvenance(fitResults);
    const holdoutDataset = calibrationDatasetProvenance(holdoutResults);
    const datasetOverlap = Boolean(fitDataset.datasetId && fitDataset.datasetId === holdoutDataset.datasetId);
    const sourceOverlap = Boolean(fitDataset.sourceFileId && fitDataset.sourceFileId === holdoutDataset.sourceFileId);
    const knownOverlapCount = datasetOverlap || sourceOverlap ? 1 : 0;
    const profileId = calibrationProfileIdentifier(profile);
    const requiredEvidence = {
      fitDatasetId: Boolean(fitDataset.datasetId),
      holdoutDatasetId: Boolean(holdoutDataset.datasetId),
      distinctDatasetIds: Boolean(fitDataset.datasetId && holdoutDataset.datasetId && fitDataset.datasetId !== holdoutDataset.datasetId),
      distinctSourceFileIds: Boolean(fitDataset.sourceFileId && holdoutDataset.sourceFileId && fitDataset.sourceFileId !== holdoutDataset.sourceFileId),
      fitRoleMarked: fitDataset.datasetRole === "fit",
      holdoutRoleMarked: holdoutDataset.datasetRole === "holdout",
      profileFrozenBeforeHoldout: Boolean(holdoutDataset.profileFrozenAt),
      holdoutProfileMatchesFrozenProfile: Boolean(profileId && holdoutDataset.profileId && holdoutDataset.profileId === profileId),
    };
    const missingEvidence = Object.entries(requiredEvidence)
      .filter(([, value]) => !value)
      .map(([key]) => key);
    const independenceReady = missingEvidence.length === 0 && knownOverlapCount === 0;
    return {
      schema: "rad-sim.calibration-train-holdout-split.v1",
      fitStepCount: Array.isArray(fitResults?.steps) ? fitResults.steps.length : 0,
      holdoutStepCount: Array.isArray(holdoutResults?.steps) ? holdoutResults.steps.length : 0,
      fitMeasuredCellCount: fitComparison.after?.metrics?.measuredCellCount ?? null,
      holdoutMeasuredCellCount: holdoutComparison.after?.metrics?.measuredCellCount ?? null,
      fitDataset,
      holdoutDataset,
      profileId,
      knownOverlapCount,
      profileFrozenBeforeHoldout: Boolean(holdoutDataset.profileFrozenAt),
      holdoutProfileMatchesFrozenProfile: requiredEvidence.holdoutProfileMatchesFrozenProfile,
      independenceReady,
      independenceStatus: independenceReady ? "documented-independent-holdout" : "requires-external-bench-protocol",
      missingEvidence,
      requiredEvidence,
      claimLabel: "experimentally unvalidated split metadata",
    };
  }

  function calibrationModelProfileResidualComparisonForResults(state, results, profile = null, options = {}) {
    if (!results) throw new Error("calibration profile residual comparison requires calibration results");
    const protocol = options.protocol || calibrationExperimentProtocol(state, options);
    const beforeState = cloneCalibrationState(state);
    const beforeComparison = compareCalibrationExperimentResults(beforeState, results, { ...options, protocol });
    const beforeSummary = summarizeCalibrationComparison(beforeComparison);
    beforeState.experiment.calibrationResults = results;
    beforeState.experiment.calibrationComparison = beforeComparison;
    beforeState.experiment.calibrationComparisonSummary = beforeSummary;
    const beforeReport = calibrationComparisonReport(beforeState);
    const resolvedProfile = profile || beforeReport.modelProfile;
    const afterState = cloneCalibrationState(state);
    afterState.experiment.calibrationResults = results;
    const application = applyCalibrationModelProfile(afterState, resolvedProfile);
    const afterComparison = compareCalibrationExperimentResults(afterState, results, { ...options, protocol });
    const afterSummary = summarizeCalibrationComparison(afterComparison);
    afterState.experiment.calibrationComparison = afterComparison;
    afterState.experiment.calibrationComparisonSummary = afterSummary;
    const afterReport = calibrationComparisonReport(afterState);
    const beforeMetrics = comparisonSummaryMetrics(beforeReport);
    const afterMetrics = comparisonSummaryMetrics(afterReport);
    const delta = residualDelta(beforeMetrics, afterMetrics);
    const minImprovement = Number(options.minImprovement) || 0;
    const improved =
      delta.residualScore !== null &&
      delta.residualScore <= -minImprovement &&
      delta.missingObservationCount <= 0;
    return {
      schema: "rad-sim.calibration-model-profile-residual-comparison.v1",
      profileSchema: resolvedProfile?.schema || null,
      applicationSchema: application.schema,
      application,
      before: {
        grid: beforeReport.grid,
        summary: beforeReport.summary,
        metrics: beforeMetrics,
      },
      after: {
        grid: afterReport.grid,
        summary: afterReport.summary,
        metrics: afterMetrics,
      },
      delta,
      improved: Boolean(improved),
      selectionScore: afterMetrics.residualScore,
      claimLabels: {
        comparison: "simulator-derived empirical law",
        profileApplication: application.claimLabel,
      },
      limitations: [
        "Residual improvement is evaluated against the same measurement file used to fit the profile.",
        "A lower simulator residual is not independent physical validation.",
        "Profile comparison does not calibrate backlash, friction, or rigid-body contact.",
      ],
    };
  }

  function calibrationModelProfileResidualComparison(state, profile = null, options = {}) {
    const results = state.experiment?.calibrationResults;
    if (!results) throw new Error("calibration profile residual comparison requires imported calibration results");
    return calibrationModelProfileResidualComparisonForResults(state, results, profile, options);
  }

  function calibrationModelProfileHoldoutValidation(state, holdoutResults = null, profile = null, options = {}) {
    const fitResults = state.experiment?.calibrationResults;
    const resolvedHoldoutResults = holdoutResults || state.experiment?.calibrationHoldoutResults;
    if (!fitResults) throw new Error("holdout validation requires imported fit calibration results");
    if (!resolvedHoldoutResults) throw new Error("holdout validation requires imported holdout calibration results");
    const protocol = options.protocol || calibrationExperimentProtocol(state, options);
    const fitState = cloneCalibrationState(state);
    fitState.experiment.calibrationResults = fitResults;
    fitState.experiment.calibrationComparison = compareCalibrationExperimentResults(fitState, fitResults, { ...options, protocol });
    fitState.experiment.calibrationComparisonSummary = summarizeCalibrationComparison(fitState.experiment.calibrationComparison);
    const fitReport = calibrationComparisonReport(fitState);
    const resolvedProfile = profile || fitReport.modelProfile;
    const fitComparison = calibrationModelProfileResidualComparisonForResults(
      state,
      fitResults,
      resolvedProfile,
      { ...options, protocol }
    );
    const holdoutComparison = calibrationModelProfileResidualComparisonForResults(
      state,
      resolvedHoldoutResults,
      resolvedProfile,
      { ...options, protocol }
    );
    const passes = (comparison) => {
      const appliedUpdates = Array.isArray(comparison?.application?.appliedUpdates)
        ? comparison.application.appliedUpdates.length
        : 0;
      const scoreDelta = profileFinite(comparison?.delta?.residualScore);
      const missingDelta = Number(comparison?.delta?.missingObservationCount) || 0;
      const minImprovement = Number(options.minImprovement) || 0;
      return appliedUpdates > 0 && scoreDelta !== null && scoreDelta <= -minImprovement && missingDelta <= 0;
    };
    const fitPass = passes(fitComparison);
    const holdoutPass = passes(holdoutComparison);
    const splitMetadata = calibrationSplitMetadata(
      fitResults,
      resolvedHoldoutResults,
      resolvedProfile,
      fitComparison,
      holdoutComparison
    );
    const residualValidationPass = Boolean(fitPass && holdoutPass);
    return {
      schema: "rad-sim.calibration-model-profile-holdout-validation.v1",
      profileSchema: resolvedProfile?.schema || null,
      modelProfile: resolvedProfile,
      fitComparisonSchema: fitComparison.schema,
      holdoutComparisonSchema: holdoutComparison.schema,
      fit: fitComparison,
      holdout: holdoutComparison,
      holdoutPass,
      fitPass,
      residualValidationPass,
      independentValidationPass: Boolean(residualValidationPass && splitMetadata.independenceReady),
      splitMetadata,
      provenanceWarnings: splitMetadata.missingEvidence.map((field) => `missing or invalid split evidence: ${field}`),
      rule: {
        name: "fit-profile-holdout-residual-nonincrease",
        minImprovement: Number(options.minImprovement) || 0,
        requiresFitImprovement: true,
        requiresHoldoutImprovement: true,
        requiresAppliedUpdate: true,
        requiresMissingObservationNonincrease: true,
      },
      claimLabels: {
        fit: "simulator-derived empirical law",
        holdout: "holdout simulator validation",
        independentValidationPass: "experimentally unvalidated until holdout measurements are independently collected",
      },
      limitations: [
        "Holdout validation is meaningful only if the holdout file was not used to fit the profile.",
        "Passing holdout residual checks supports simulator calibration but does not prove contact, friction, or material laws.",
        "The current v1 profile mutates only bounded simulator fields such as z-coupling.",
      ],
    };
  }

  function selectCalibrationModelProfile(comparisons, options = {}) {
    const minImprovement = Number(options.minImprovement) || 0;
    const candidates = (comparisons || []).map((candidate) => JSON.parse(JSON.stringify(candidate)));
    const scored = [];
    candidates.forEach((candidate, index) => {
      const afterMetrics = candidate?.after?.metrics || {};
      const delta = candidate?.delta || {};
      const score = profileFinite(afterMetrics.residualScore);
      const scoreDelta = profileFinite(delta.residualScore);
      const missingDelta = Number(delta.missingObservationCount) || 0;
      const appliedUpdates = Array.isArray(candidate?.application?.appliedUpdates)
        ? candidate.application.appliedUpdates.length
        : 0;
      const eligible =
        score !== null &&
        scoreDelta !== null &&
        scoreDelta <= -minImprovement &&
        missingDelta <= 0 &&
        appliedUpdates > 0;
      candidate.selectionEligible = Boolean(eligible);
      candidate.selectionRejectionReason = eligible
        ? null
        : "requires applied updates, non-increased missing count, and residual-score improvement";
      if (eligible) scored.push({ score, missingDelta, index, candidate });
    });
    scored.sort((a, b) => a.score - b.score || a.missingDelta - b.missingDelta || a.index - b.index);
    const selected = scored.length ? scored[0] : null;
    return {
      schema: "rad-sim.calibration-model-profile-selection.v1",
      candidateCount: candidates.length,
      eligibleCount: scored.length,
      selectedIndex: selected ? selected.index : null,
      selected: selected ? selected.candidate : null,
      candidates,
      rule: {
        name: "lowest-after-residual-score-with-improvement",
        minImprovement,
        requiresAppliedUpdate: true,
        requiresMissingObservationNonincrease: true,
      },
      claimLabel: "simulator-derived empirical profile selection rule",
      limitations: [
        "Selection is based on simulator residuals against available measurements.",
        "Independent validation data is still required before treating a profile as physically calibrated.",
      ],
    };
  }

  function exportCalibrationModelProfileResidualComparison(state, profile = null, options = {}) {
    return JSON.stringify(calibrationModelProfileResidualComparison(state, profile, options), null, 2);
  }

  function exportCalibrationModelProfileHoldoutValidation(state, holdoutResults = null, profile = null, options = {}) {
    return JSON.stringify(calibrationModelProfileHoldoutValidation(state, holdoutResults, profile, options), null, 2);
  }

  function holdoutValidationMetric(comparison, stage, key) {
    const metrics = comparison?.[stage]?.metrics || {};
    return metrics[key] ?? "";
  }

  function exportCalibrationModelProfileHoldoutValidationCsv(validationOrState) {
    const validation =
      validationOrState?.schema === "rad-sim.calibration-model-profile-holdout-validation.v1"
        ? validationOrState
        : validationOrState?.experiment?.calibrationModelProfileHoldoutValidation;
    if (!validation) throw new Error("no calibration holdout validation is available");
    const split = validation.splitMetadata || {};
    const header = [
      "schema",
      "dataset",
      "pass",
      "independent_validation_pass",
      "residual_validation_pass",
      "profile_schema",
      "fit_dataset_id",
      "holdout_dataset_id",
      "applied_updates",
      "z_coupling_before",
      "z_coupling_after",
      "before_residual_score",
      "after_residual_score",
      "residual_score_delta",
      "before_missing_observations",
      "after_missing_observations",
      "missing_observation_delta",
      "before_measured_cells",
      "after_measured_cells",
      "profile_frozen_before_holdout",
      "known_overlap_count",
      "independence_status",
      "independence_ready",
      "missing_evidence",
      "claim_label",
    ];
    const lines = [header.join(",")];
    for (const [dataset, key, passKey] of [
      ["fit", "fit", "fitPass"],
      ["holdout", "holdout", "holdoutPass"],
    ]) {
      const comparison = validation[key] || {};
      const application = comparison.application || {};
      const appliedUpdates = Array.isArray(application.appliedUpdates) ? application.appliedUpdates.length : 0;
      const configBefore = application.configBefore || {};
      const configAfter = application.configAfter || {};
      const delta = comparison.delta || {};
      lines.push(
        [
          validation.schema,
          dataset,
          validation[passKey],
          validation.independentValidationPass,
          validation.residualValidationPass,
          validation.profileSchema || "",
          split.fitDataset?.datasetId ?? "",
          split.holdoutDataset?.datasetId ?? "",
          appliedUpdates,
          configBefore.zCouplingGain ?? "",
          configAfter.zCouplingGain ?? "",
          holdoutValidationMetric(comparison, "before", "residualScore"),
          holdoutValidationMetric(comparison, "after", "residualScore"),
          delta.residualScore ?? "",
          holdoutValidationMetric(comparison, "before", "missingObservationCount"),
          holdoutValidationMetric(comparison, "after", "missingObservationCount"),
          delta.missingObservationCount ?? "",
          holdoutValidationMetric(comparison, "before", "measuredCellCount"),
          holdoutValidationMetric(comparison, "after", "measuredCellCount"),
          split.profileFrozenBeforeHoldout ?? "",
          split.knownOverlapCount ?? "",
          split.independenceStatus ?? "",
          split.independenceReady ?? "",
          Array.isArray(split.missingEvidence) ? split.missingEvidence.join(";") : "",
          split.claimLabel ?? "",
        ].map(csvValue).join(",")
      );
    }
    return lines.join("\n");
  }

  function calibrationBenchExecutionValidation(state, holdoutResults = null, profile = null, options = {}) {
    const fitResults = state.experiment?.calibrationResults;
    const resolvedHoldoutResults = holdoutResults || state.experiment?.calibrationHoldoutResults;
    if (!fitResults) throw new Error("bench execution validation requires imported fit calibration results");
    if (!resolvedHoldoutResults) throw new Error("bench execution validation requires imported holdout calibration results");
    const protocol = options.protocol || calibrationExperimentProtocol(state, options);
    const fitState = cloneCalibrationState(state);
    fitState.experiment.calibrationResults = fitResults;
    fitState.experiment.calibrationComparison = compareCalibrationExperimentResults(fitState, fitResults, { ...options, protocol });
    fitState.experiment.calibrationComparisonSummary = summarizeCalibrationComparison(fitState.experiment.calibrationComparison);
    const fitReport = calibrationComparisonReport(fitState);
    const resolvedProfile = profile || fitReport.modelProfile;
    const holdoutValidation = calibrationModelProfileHoldoutValidation(
      state,
      resolvedHoldoutResults,
      resolvedProfile,
      { ...options, protocol }
    );
    const split = holdoutValidation.splitMetadata || {};
    const missingEvidence = Array.isArray(split.missingEvidence) ? split.missingEvidence : [];
    const appliedUpdates = Array.isArray(holdoutValidation.fit?.application?.appliedUpdates)
      ? holdoutValidation.fit.application.appliedUpdates.length
      : 0;
    const executionValidationPass = Boolean(holdoutValidation.independentValidationPass && missingEvidence.length === 0);
    const status = executionValidationPass
      ? "documented-independent-validation"
      : holdoutValidation.residualValidationPass
        ? "residual-pass-needs-provenance"
        : "residual-validation-failed";
    return {
      schema: "rad-sim.calibration-bench-execution-validation.v1",
      method: "filled fit/holdout calibration packet ingestion using bounded model-profile residual replay and split-provenance checks",
      inputSchemas: {
        protocol: protocol.schema,
        fitResults: "rad-sim.calibration-experiment-results.v1",
        holdoutResults: "rad-sim.calibration-experiment-results.v1",
        datasetProvenance: "rad-sim.calibration-dataset-provenance.v1",
        holdoutValidation: "rad-sim.calibration-model-profile-holdout-validation.v1",
      },
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        cellSize: state.grid.cellSize,
        backlash: state.grid.backlash,
        couplingGain: state.grid.couplingGain,
        zCouplingGain: state.grid.zCouplingGain,
        pinHoleClearance: typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : Math.max(0, (state.grid.holeRadius || 0) - (state.grid.pinRadius || 0)),
      },
      protocol: {
        schema: protocol.schema,
        hardwareProfile: protocol.hardwareProfile,
        stepCount: (protocol.steps || []).length,
        repeatCountTotal: (protocol.steps || []).reduce((sum, step) => sum + (Number(step.repeatCount) || 0), 0),
        centerCell: protocol.centerCell,
      },
      fitDataset: fitResults.provenance || null,
      holdoutDataset: resolvedHoldoutResults.provenance || null,
      fitComparisonReport: fitReport,
      modelProfile: resolvedProfile,
      holdoutValidation,
      holdoutValidationCsv: exportCalibrationModelProfileHoldoutValidationCsv(holdoutValidation),
      summary: {
        status,
        fitStepCount: (fitResults.steps || []).length,
        holdoutStepCount: (resolvedHoldoutResults.steps || []).length,
        appliedUpdateCount: appliedUpdates,
        residualValidationPass: Boolean(holdoutValidation.residualValidationPass),
        independentValidationPass: Boolean(holdoutValidation.independentValidationPass),
        executionValidationPass,
        missingEvidenceCount: missingEvidence.length,
        independenceStatus: split.independenceStatus || null,
        knownOverlapCount: split.knownOverlapCount ?? null,
      },
      formalization: {
        targetId: "calibration_bench_executed_validation_gate",
        leanStructure: "Mechanics.CalibrationBenchExecutedValidationNat",
        leanPredicate: "calibrationBenchExecutedValidationReadyNat",
        schema: "rad-sim.calibration-bench-execution-validation.v1",
      },
      claimLabels: {
        executionValidation: "bench-data residual validation bookkeeping",
        modelProfile: "simulator-derived empirical law",
        independentValidationPass: "experimentally unvalidated unless provenance is externally true",
        physicalAccuracy: "experimentally unvalidated physical assumption",
      },
      limitations: [
        "The report can detect missing provenance fields but cannot prove the lab followed the protocol.",
        "Residual agreement supports simulator-profile bookkeeping only; it does not validate contact, friction, stiffness, gravity, or material laws.",
        "The current bounded profile application mutates only implemented simulator fields.",
      ],
    };
  }

  function exportCalibrationBenchExecutionValidation(state, holdoutResults = null, profile = null, options = {}) {
    return JSON.stringify(calibrationBenchExecutionValidation(state, holdoutResults, profile, options), null, 2);
  }

  function exportCalibrationBenchExecutionValidationCsv(validationOrState) {
    const validation =
      validationOrState?.schema === "rad-sim.calibration-bench-execution-validation.v1"
        ? validationOrState
        : validationOrState?.experiment?.calibrationBenchExecutionValidation;
    if (!validation) throw new Error("no calibration bench execution validation is available");
    const summary = validation.summary || {};
    const fitDataset = validation.fitDataset || {};
    const holdoutDataset = validation.holdoutDataset || {};
    const split = validation.holdoutValidation?.splitMetadata || {};
    const header = [
      "schema",
      "status",
      "execution_validation_pass",
      "residual_validation_pass",
      "independent_validation_pass",
      "fit_dataset_id",
      "holdout_dataset_id",
      "applied_updates",
      "fit_step_count",
      "holdout_step_count",
      "missing_evidence_count",
      "missing_evidence",
      "independence_status",
      "known_overlap_count",
      "claim_label",
    ];
    return [
      header,
      [
        validation.schema,
        summary.status || "",
        summary.executionValidationPass,
        summary.residualValidationPass,
        summary.independentValidationPass,
        fitDataset.datasetId || "",
        holdoutDataset.datasetId || "",
        summary.appliedUpdateCount ?? "",
        summary.fitStepCount ?? "",
        summary.holdoutStepCount ?? "",
        summary.missingEvidenceCount ?? "",
        Array.isArray(split.missingEvidence) ? split.missingEvidence.join(";") : "",
        summary.independenceStatus || "",
        summary.knownOverlapCount ?? "",
        validation.claimLabels?.executionValidation || "",
      ],
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function verticalLoadProxyTermCounts(verticalLoadComparisonReport) {
    if (!verticalLoadComparisonReport || typeof verticalLoadComparisonReport !== "object") return { load: 0, contact: 0 };
    let load = 0;
    let contact = 0;
    for (const scenario of verticalLoadComparisonReport.scenarios || []) {
      const validation = scenario?.validation || {};
      for (const side of ["intact", "removed"]) {
        const sidePayload = validation?.[side] || {};
        const errors = sidePayload.errors || {};
        const measured = sidePayload.measured || {};
        const simulated = sidePayload.simulatedObserved || {};
        if (
          Object.prototype.hasOwnProperty.call(errors, "signedLoadWork") &&
          Object.prototype.hasOwnProperty.call(measured, "signedLoadWork") &&
          Object.prototype.hasOwnProperty.call(simulated, "signedLoadWork") &&
          Object.prototype.hasOwnProperty.call(errors, "loadWorkMagnitude") &&
          Object.prototype.hasOwnProperty.call(measured, "loadWorkMagnitude") &&
          Object.prototype.hasOwnProperty.call(simulated, "loadWorkMagnitude")
        ) {
          load += 1;
        }
        if (
          Object.prototype.hasOwnProperty.call(errors, "heightContactPenalty") &&
          Object.prototype.hasOwnProperty.call(measured, "heightContactPenalty") &&
          Object.prototype.hasOwnProperty.call(simulated, "heightContactPenalty")
        ) {
          contact += 1;
        }
      }
    }
    return { load, contact };
  }

  function physicalValidationReadinessReport(state, options = {}) {
    const calibrationExecution =
      options.calibrationExecutionValidation || state.experiment?.calibrationBenchExecutionValidation;
    const verticalLoadComparison =
      options.verticalLoadComparisonReport || state.experiment?.verticalLoadComparisonReport;
    const calibrationSummary = calibrationExecution?.summary || {};
    const verticalSummary = verticalLoadComparison?.summary || {};
    const proxyCounts = verticalLoadProxyTermCounts(verticalLoadComparison);
    const pinRadius = Number(state.grid.pinRadius) || 0;
    const holeRadius = Number(state.grid.holeRadius) || 0;
    const clearance = typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : Math.max(0, holeRadius - pinRadius);
    const calibrationPass = Boolean(calibrationSummary.executionValidationPass);
    const verticalScenarios = Number(verticalSummary.scenarioCount || 0);
    const verticalMissing = Number(verticalSummary.missingMeasurementCount || 0);
    const verticalPass = Boolean(verticalSummary.allScenariosPassTolerance);
    const clearanceConfigured = pinRadius > 0 && holeRadius >= pinRadius && clearance >= 0;
    const missingEvidence = [];
    if (!calibrationExecution || typeof calibrationExecution !== "object") missingEvidence.push("calibrationBenchExecutionValidation");
    else if (!calibrationPass) missingEvidence.push("calibrationExecutionValidationPass");
    if (!verticalLoadComparison || typeof verticalLoadComparison !== "object") {
      missingEvidence.push("verticalLoadComparisonReport");
    } else {
      if (verticalScenarios <= 0) missingEvidence.push("verticalLoadScenarioCount");
      if (verticalMissing > 0) missingEvidence.push("verticalLoadMissingMeasurements");
      if (!verticalPass) missingEvidence.push("verticalLoadTolerancePass");
    }
    if (proxyCounts.load <= 0) missingEvidence.push("signedLoadWorkAndMagnitudeTerms");
    if (proxyCounts.contact <= 0) missingEvidence.push("heightContactPenaltyTerms");
    if (!clearanceConfigured) missingEvidence.push("pinHoleClearanceConfiguration");
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.physical-validation-readiness.v1",
      method: "conservative evidence gate combining executed calibration validation with filled vertical-load work/contact comparison reports",
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        pinRadius,
        holeRadius,
        pinHoleClearance: clearance,
        zCouplingGain: state.grid.zCouplingGain,
      },
      inputSchemas: {
        calibrationExecution: "rad-sim.calibration-bench-execution-validation.v1",
        verticalLoadComparison: "rad-sim.vertical-load-energy-comparison-report.v1",
      },
      evidence: {
        calibrationExecutionValidationPass: calibrationPass,
        fitStepCount: calibrationSummary.fitStepCount || 0,
        holdoutStepCount: calibrationSummary.holdoutStepCount || 0,
        verticalLoadScenarioCount: verticalScenarios,
        verticalLoadMissingMeasurementCount: verticalMissing,
        verticalLoadAllScenariosPassTolerance: verticalPass,
        loadProxyTermCount: proxyCounts.load,
        contactProxyTermCount: proxyCounts.contact,
        pinHoleClearanceConfigured: clearanceConfigured,
      },
      summary: {
        status: ready ? "ready-for-physical-claim-review" : "needs-bench-evidence",
        physicalValidationReady: ready,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "physical_validation_readiness_gate",
        leanStructure: "Mechanics.PhysicalValidationReadinessNat",
        leanPredicate: "physicalValidationReadyNat",
        schema: "rad-sim.physical-validation-readiness.v1",
      },
      claimLabels: {
        readinessGate: "bench-data evidence bookkeeping",
        loadProxy: "Lean-proven finite scaffold plus simulator-derived empirical law",
        contactProxy: "experimentally unvalidated physical assumption",
        physicalAccuracy: "experimentally unvalidated physical assumption",
      },
      limitations: [
        "Readiness means the required evidence artifacts are present and internally passing.",
        "It does not prove calibrated rigid-body contact, friction, stiffness, gravity, or material law.",
        "Physical claims still require lab audit of fixtures, sensors, raw files, and uncertainty.",
      ],
    };
  }

  function exportPhysicalValidationReadiness(state, options = {}) {
    return JSON.stringify(physicalValidationReadinessReport(state, options), null, 2);
  }

  function exportPhysicalValidationReadinessCsv(reportOrState, options = {}) {
    const report =
      reportOrState?.schema === "rad-sim.physical-validation-readiness.v1"
        ? reportOrState
        : reportOrState?.experiment?.physicalValidationReadiness || physicalValidationReadinessReport(reportOrState, options);
    const summary = report.summary || {};
    const evidence = report.evidence || {};
    const header = [
      "schema",
      "status",
      "physical_validation_ready",
      "missing_evidence_count",
      "missing_evidence",
      "calibration_execution_pass",
      "fit_step_count",
      "holdout_step_count",
      "vertical_load_scenarios",
      "vertical_load_missing_measurements",
      "vertical_load_pass",
      "load_proxy_terms",
      "contact_proxy_terms",
      "pin_hole_clearance_configured",
      "claim_label",
    ];
    return [
      header,
      [
        report.schema,
        summary.status || "",
        summary.physicalValidationReady,
        summary.missingEvidenceCount,
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
        evidence.calibrationExecutionValidationPass,
        evidence.fitStepCount,
        evidence.holdoutStepCount,
        evidence.verticalLoadScenarioCount,
        evidence.verticalLoadMissingMeasurementCount,
        evidence.verticalLoadAllScenariosPassTolerance,
        evidence.loadProxyTermCount,
        evidence.contactProxyTermCount,
        evidence.pinHoleClearanceConfigured,
        report.claimLabels?.readinessGate || "",
      ],
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function contactStateAbstractionReport(state, options = {}) {
    const contactStiffness = Math.max(0, Number(options.contactStiffness ?? 1));
    const tolerance = Math.max(0, Number(options.tolerance ?? 1e-9));
    const verticalLoadComparison = options.verticalLoadComparisonReport || state.experiment?.verticalLoadComparisonReport;
    const proxyCounts = verticalLoadProxyTermCounts(verticalLoadComparison);
    const pinRadius = Number(state.grid.pinRadius) || 0;
    const holeRadius = Number(state.grid.holeRadius) || 0;
    const clearance = typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : Math.max(0, holeRadius - pinRadius);
    const cells = [];
    let activeBodyCount = 0;
    let removedBodyCount = 0;
    let engagedContactCount = 0;
    let totalContactPenalty = 0;
    for (let r = 0; r < state.grid.rows; r += 1) {
      for (let c = 0; c < state.grid.cols; c += 1) {
        const removed = Boolean(state.cells.removed?.[r]?.[c]);
        const zCommand = Number(state.cells.commandZ?.[r]?.[c]) || 0;
        const penetration = removed ? 0 : Math.max(0, Math.abs(zCommand) - clearance);
        const contactPenalty = 0.5 * contactStiffness * penetration * penetration;
        const mode = removed ? "removed" : penetration > tolerance ? "engaged" : "free-clearance";
        if (removed) removedBodyCount += 1;
        else activeBodyCount += 1;
        if (mode === "engaged") engagedContactCount += 1;
        totalContactPenalty += contactPenalty;
        cells.push({
          row: r,
          col: c,
          removed,
          zCommand,
          pinRadius,
          holeRadius,
          clearance,
          penetration,
          contactPenalty,
          mode,
          contactState: {
            pinPresent: !removed,
            holePresent: !removed,
            unilateralContactActive: mode === "engaged",
            clearanceExceeded: penetration > tolerance,
          },
        });
      }
    }
    const clearanceConfigured = pinRadius > 0 && holeRadius >= pinRadius && clearance >= 0;
    const missingEvidence = [];
    if (activeBodyCount <= 0) missingEvidence.push("activeBodies");
    if (pinRadius <= 0) missingEvidence.push("pinRadius");
    if (holeRadius < pinRadius) missingEvidence.push("holeRadiusAtLeastPinRadius");
    if (!clearanceConfigured) missingEvidence.push("clearanceConfiguration");
    if (contactStiffness <= 0) missingEvidence.push("positiveContactStiffness");
    if (!cells.length) missingEvidence.push("contactStateRecords");
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.contact-state-abstraction.v1",
      method: "normalized per-cell pin-hole clearance and unilateral contact-state bookkeeping derived from z commands, removed topology, and optional vertical-load comparison evidence",
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        pinRadius,
        holeRadius,
        pinHoleClearance: clearance,
        contactStiffness,
      },
      topology: {
        activeBodyCount,
        removedBodyCount,
        pinCount: activeBodyCount,
        holeCount: activeBodyCount,
        clearancePairCount: activeBodyCount,
        contactStateRecordCount: cells.length,
        penaltyTermCount: activeBodyCount,
      },
      summary: {
        status: ready ? "contact-state-abstraction-ready" : "needs-contact-abstraction-inputs",
        contactStateAbstractionReady: ready,
        engagedContactCount,
        totalContactPenalty,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      benchEvidence: {
        verticalLoadComparisonAttached: Boolean(verticalLoadComparison && typeof verticalLoadComparison === "object"),
        loadProxyTermCount: proxyCounts.load,
        contactProxyTermCount: proxyCounts.contact,
      },
      cells,
      formalization: {
        targetId: "contact_state_abstraction_gate",
        leanStructure: "Mechanics.ContactStateAbstractionNat",
        leanPredicate: "contactStateAbstractionReadyNat",
        schema: "rad-sim.contact-state-abstraction.v1",
      },
      claimLabels: {
        contactState: "graph/contact bookkeeping abstraction",
        contactPenalty: "simulator-derived empirical law",
        physicalContact: "experimentally unvalidated physical assumption",
      },
      limitations: [
        "The contact mode is derived from normalized z command and clearance, not collision detection.",
        "The penalty is a unilateral contact proxy until contact stiffness and friction are measured.",
        "Removed cells are deleted topology records, not simulated detached rigid bodies.",
      ],
    };
  }

  function exportContactStateAbstraction(state, options = {}) {
    return JSON.stringify(contactStateAbstractionReport(state, options), null, 2);
  }

  function exportContactStateAbstractionCsv(reportOrState, options = {}) {
    const report =
      reportOrState?.schema === "rad-sim.contact-state-abstraction.v1"
        ? reportOrState
        : reportOrState?.experiment?.contactStateAbstraction || contactStateAbstractionReport(reportOrState, options);
    const rows = [[
      "row",
      "col",
      "removed",
      "mode",
      "z_command",
      "pin_radius",
      "hole_radius",
      "clearance",
      "penetration",
      "contact_penalty",
    ]];
    for (const cell of report.cells || []) {
      rows.push([
        cell.row,
        cell.col,
        cell.removed,
        cell.mode,
        cell.zCommand,
        cell.pinRadius,
        cell.holeRadius,
        cell.clearance,
        cell.penetration,
        cell.contactPenalty,
      ]);
    }
    return rows.map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function contactGraphKey(row, col) {
    return `${row},${col}`;
  }

  function normalizeContactGraphCell(cell) {
    const row = Number(cell?.r ?? cell?.row);
    const col = Number(cell?.c ?? cell?.col);
    if (!Number.isInteger(row) || !Number.isInteger(col)) return null;
    return { row, col };
  }

  function graphEdgesForState(state) {
    const edges = [];
    let activeEdgeCount = 0;
    let deletedEdgeCount = 0;
    let removedIncidentActiveEdgeCount = 0;
    for (let row = 0; row < state.grid.rows; row += 1) {
      for (let col = 0; col < state.grid.cols; col += 1) {
        for (const [nr, nc] of [[row, col + 1], [row + 1, col]]) {
          if (nr >= state.grid.rows || nc >= state.grid.cols) continue;
          const touchesRemoved = Boolean(state.cells.removed?.[row]?.[col] || state.cells.removed?.[nr]?.[nc]);
          const active = !touchesRemoved;
          if (active) activeEdgeCount += 1;
          else deletedEdgeCount += 1;
          if (active && touchesRemoved) removedIncidentActiveEdgeCount += 1;
          edges.push({
            from: { row, col },
            to: { row: nr, col: nc },
            active,
            touchesRemoved,
            status: active ? "active" : "deleted-by-removal",
          });
        }
      }
    }
    return { edges, activeEdgeCount, deletedEdgeCount, removedIncidentActiveEdgeCount };
  }

  function supportFromState(state, tolerance = 1e-9) {
    const support = [];
    for (let row = 0; row < state.grid.rows; row += 1) {
      for (let col = 0; col < state.grid.cols; col += 1) {
        const alpha = Math.abs(Number(state.cells.commandAlpha?.[row]?.[col]) || 0);
        const z = Math.abs(Number(state.cells.commandZ?.[row]?.[col]) || 0);
        const locked = Boolean(state.cells.locked?.[row]?.[col]);
        if (alpha > tolerance || z > tolerance || locked) support.push({ row, col });
      }
    }
    return support;
  }

  function contactGraphConsistencyReport(state, options = {}) {
    const contactStiffness = Math.max(0, Number(options.contactStiffness ?? 1));
    const tolerance = Math.max(0, Number(options.tolerance ?? 1e-9));
    const contactReport =
      options.contactReport ||
      state.experiment?.contactStateAbstraction ||
      contactStateAbstractionReport(state, { contactStiffness, tolerance });
    const contactCells = Array.isArray(contactReport?.cells) ? contactReport.cells : [];
    const contactByCell = new Map();
    for (const cell of contactCells) {
      const normalized = normalizeContactGraphCell(cell);
      if (normalized) contactByCell.set(contactGraphKey(normalized.row, normalized.col), cell);
    }
    const graph = graphEdgesForState(state);
    let activeBodyCount = 0;
    let removedBodyCount = 0;
    let activeContactRecordCount = 0;
    let missingActiveContactRecordCount = 0;
    let removedActiveContactCount = 0;
    for (let row = 0; row < state.grid.rows; row += 1) {
      for (let col = 0; col < state.grid.cols; col += 1) {
        const removed = Boolean(state.cells.removed?.[row]?.[col]);
        if (removed) removedBodyCount += 1;
        else activeBodyCount += 1;
        const record = contactByCell.get(contactGraphKey(row, col));
        if (!record) {
          if (!removed) missingActiveContactRecordCount += 1;
          continue;
        }
        const activeContact = Boolean(record.contactState?.unilateralContactActive);
        if (!removed) activeContactRecordCount += 1;
        if (removed && (activeContact || record.mode === "engaged")) removedActiveContactCount += 1;
      }
    }
    const rawSupport = Array.isArray(options.groupSupport) ? options.groupSupport : supportFromState(state, tolerance);
    const validSupport = [];
    const skippedSupport = [];
    let supportContactRecordCount = 0;
    let activeSupportCellCount = 0;
    let removedSupportCellCount = 0;
    for (const cell of rawSupport) {
      const normalized = normalizeContactGraphCell(cell);
      if (!normalized || normalized.row < 0 || normalized.row >= state.grid.rows || normalized.col < 0 || normalized.col >= state.grid.cols) {
        if (normalized) skippedSupport.push(normalized);
        continue;
      }
      validSupport.push(normalized);
      if (contactByCell.has(contactGraphKey(normalized.row, normalized.col))) supportContactRecordCount += 1;
      if (state.cells.removed?.[normalized.row]?.[normalized.col]) removedSupportCellCount += 1;
      else activeSupportCellCount += 1;
    }
    const missingEvidence = [];
    if (!contactReport || typeof contactReport !== "object") missingEvidence.push("contactStateAbstraction");
    else if (!contactReport.summary?.contactStateAbstractionReady) missingEvidence.push("contactStateAbstractionReady");
    if (activeBodyCount <= 0) missingEvidence.push("activeBodies");
    if (missingActiveContactRecordCount > 0) missingEvidence.push("activeContactRecords");
    if (graph.removedIncidentActiveEdgeCount > 0) missingEvidence.push("removedIncidentActiveEdges");
    if (removedActiveContactCount > 0) missingEvidence.push("removedActiveContacts");
    if (supportContactRecordCount < validSupport.length) missingEvidence.push("supportContactRecords");
    const contactGraphConsistent = missingEvidence.length === 0;
    return {
      schema: "rad-sim.contact-graph-consistency.v1",
      method: "checks that graph edges touching removed cells are deleted, contact records are inactive on removed cells, and support cells are represented in the contact-state abstraction",
      inputSchemas: {
        contactState: "rad-sim.contact-state-abstraction.v1",
        state: "rad-sim lattice state with removed/contact support masks",
      },
      graph: {
        activeBodyCount,
        removedBodyCount,
        totalEdgeCount: graph.edges.length,
        activeEdgeCount: graph.activeEdgeCount,
        deletedEdgeCount: graph.deletedEdgeCount,
        removedIncidentActiveEdgeCount: graph.removedIncidentActiveEdgeCount,
        edges: graph.edges,
      },
      contact: {
        contactRecordCount: contactByCell.size,
        activeContactRecordCount,
        missingActiveContactRecordCount,
        removedActiveContactCount,
        engagedContactCount: Number(contactReport?.summary?.engagedContactCount) || 0,
      },
      support: {
        supportCells: validSupport,
        skippedSupportCells: skippedSupport,
        supportCellCount: validSupport.length,
        activeSupportCellCount,
        removedSupportCellCount,
        supportContactRecordCount,
      },
      summary: {
        status: contactGraphConsistent ? "contact-graph-consistent" : "needs-contact-graph-evidence",
        contactGraphConsistent,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "contact_graph_consistency_gate",
        leanStructure: "Mechanics.ContactGraphConsistencyNat",
        leanPredicate: "contactGraphConsistentNat",
        schema: "rad-sim.contact-graph-consistency.v1",
      },
      claimLabels: {
        graphDeletion: "Lean-proven graph deletion theorem scaffold",
        contactRecords: "graph/contact bookkeeping abstraction",
        support: "simulator-derived event support diagnostic",
        physicalContact: "experimentally unvalidated physical assumption",
      },
      limitations: [
        "Consistency only checks graph/contact bookkeeping against removed topology.",
        "It does not prove that inferred contact modes match hardware contact.",
        "Removed support is reported as lost support, not treated as an error by itself.",
      ],
    };
  }

  function exportContactGraphConsistency(state, options = {}) {
    return JSON.stringify(contactGraphConsistencyReport(state, options), null, 2);
  }

  function exportContactGraphConsistencyCsv(reportOrState, options = {}) {
    const report =
      reportOrState?.schema === "rad-sim.contact-graph-consistency.v1"
        ? reportOrState
        : reportOrState?.experiment?.contactGraphConsistency || contactGraphConsistencyReport(reportOrState, options);
    const graph = report.graph || {};
    const contact = report.contact || {};
    const support = report.support || {};
    const summary = report.summary || {};
    return [
      [
        "schema",
        "status",
        "contact_graph_consistent",
        "active_bodies",
        "removed_bodies",
        "active_edges",
        "deleted_edges",
        "removed_incident_active_edges",
        "contact_records",
        "removed_active_contacts",
        "support_cells",
        "support_contact_records",
        "removed_support_cells",
        "missing_evidence",
      ],
      [
        report.schema,
        summary.status || "",
        summary.contactGraphConsistent,
        graph.activeBodyCount ?? "",
        graph.removedBodyCount ?? "",
        graph.activeEdgeCount ?? "",
        graph.deletedEdgeCount ?? "",
        graph.removedIncidentActiveEdgeCount ?? "",
        contact.contactRecordCount ?? "",
        contact.removedActiveContactCount ?? "",
        support.supportCellCount ?? "",
        support.supportContactRecordCount ?? "",
        support.removedSupportCellCount ?? "",
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
      ],
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function realizationCommandSupport(state, tolerance = 1e-9) {
    const support = [];
    for (let row = 0; row < state.grid.rows; row += 1) {
      for (let col = 0; col < state.grid.cols; col += 1) {
        const alpha = Math.abs(Number(state.cells.commandAlpha?.[row]?.[col]) || 0);
        const z = Math.abs(Number(state.cells.commandZ?.[row]?.[col]) || 0);
        if (alpha > tolerance || z > tolerance) support.push({ row, col });
      }
    }
    return support;
  }

  function inferRealizationEvents(state, tolerance = 1e-9) {
    const events = [];
    for (let row = 0; row < state.grid.rows; row += 1) {
      for (let col = 0; col < state.grid.cols; col += 1) {
        const alpha = Number(state.cells.commandAlpha?.[row]?.[col]) || 0;
        const z = Number(state.cells.commandZ?.[row]?.[col]) || 0;
        if (Math.abs(alpha) > tolerance || Math.abs(z) > tolerance) {
          events.push({ kind: "actuate", cell: { r: row, c: col }, alpha, z });
        }
        if (state.cells.locked?.[row]?.[col]) events.push({ kind: "lock", cell: { r: row, c: col } });
        if (state.cells.removed?.[row]?.[col]) events.push({ kind: "remove_cell", cell: { r: row, c: col } });
      }
    }
    return events;
  }

  function eventOperatorClass(kind) {
    return {
      actuate: "local actuation operator",
      group_actuate: "group actuation operator",
      lock: "constraint activation operator",
      release: "constraint release operator",
      clear_actuation: "command reset operator",
      remove_cell: "graph deletion operator",
      restore_cell: "graph restoration operator",
    }[kind] || "unknown programmable operator";
  }

  function eventHardwareChannel(kind) {
    return {
      actuate: "linear alpha/z actuator channel",
      group_actuate: "simultaneous multi-cell actuator channel",
      lock: "cell lock or constraint latch",
      release: "lock release channel",
      clear_actuation: "actuator command reset",
      remove_cell: "removable cell / deleted topology",
      restore_cell: "restored cell / topology insertion",
    }[kind] || "unmapped mechanism channel";
  }

  function realizationEventSupport(state, event, tolerance = 1e-9) {
    const candidates = [];
    const primary = normalizeContactGraphCell(event?.cell);
    if (primary) candidates.push(primary);
    for (const cell of event?.cells || []) {
      const normalized = normalizeContactGraphCell(cell);
      if (normalized) candidates.push(normalized);
    }
    if (event?.kind === "clear_actuation" && !primary && candidates.length === 0) {
      candidates.push(...realizationCommandSupport(state, tolerance));
    }
    const valid = [];
    const invalid = [];
    for (const cell of candidates) {
      if (cell.row >= 0 && cell.row < state.grid.rows && cell.col >= 0 && cell.col < state.grid.cols) valid.push(cell);
      else invalid.push(cell);
    }
    return { valid, invalid };
  }

  function realizationMaxMatrixDiff(a, b) {
    let out = 0;
    for (let r = 0; r < a.length; r += 1) {
      for (let c = 0; c < a[r].length; c += 1) {
        out = Math.max(out, Math.abs((Number(a[r][c]) || 0) - (Number(b?.[r]?.[c]) || 0)));
      }
    }
    return out;
  }

  function realizationBoolMatrixDiff(a, b) {
    let count = 0;
    for (let r = 0; r < a.length; r += 1) {
      for (let c = 0; c < a[r].length; c += 1) {
        if (Boolean(a[r][c]) !== Boolean(b?.[r]?.[c])) count += 1;
      }
    }
    return count;
  }

  function realizationStateEffect(before, after, tolerance = 1e-9) {
    const commandAlphaDelta = realizationMaxMatrixDiff(after.cells.commandAlpha, before.cells.commandAlpha);
    const commandZDelta = realizationMaxMatrixDiff(after.cells.commandZ, before.cells.commandZ);
    const lockDeltaCount = realizationBoolMatrixDiff(after.cells.locked, before.cells.locked);
    const removedDeltaCount = realizationBoolMatrixDiff(after.cells.removed || [], before.cells.removed || []);
    const active =
      commandAlphaDelta > tolerance ||
      commandZDelta > tolerance ||
      lockDeltaCount > 0 ||
      removedDeltaCount > 0;
    return { active, commandAlphaDelta, commandZDelta, lockDeltaCount, removedDeltaCount };
  }

  function realizationEventRecord(event, index, before, after, support, invalid, error, tolerance) {
    const effect = realizationStateEffect(before, after, tolerance);
    const removedSupportCells = support.filter((cell) => Boolean(before.cells.removed?.[cell.row]?.[cell.col]));
    const evidenceGaps = [];
    if (!support.length) evidenceGaps.push("support");
    if (invalid.length) evidenceGaps.push("validSupport");
    if (error) evidenceGaps.push("eventApplication");
    if (!effect.active) evidenceGaps.push("stateEffect");
    if (!["actuate", "group_actuate", "lock", "release", "clear_actuation", "remove_cell", "restore_cell"].includes(event?.kind)) {
      evidenceGaps.push("operatorKind");
    }
    return {
      index,
      kind: event?.kind || "",
      operatorClass: eventOperatorClass(event?.kind),
      hardwareChannel: eventHardwareChannel(event?.kind),
      cell: normalizeContactGraphCell(event?.cell),
      supportCells: support,
      invalidSupportCells: invalid,
      removedSupportCells,
      alpha: Number(event?.alpha) || 0,
      z: Number(event?.z) || 0,
      stateEffect: effect,
      realized: Boolean(effect.active && !error && invalid.length === 0),
      evidenceGaps,
      claimLabel: "simulator-derived realization map; hardware channel remains experimentally unvalidated",
      error,
    };
  }

  function physicalRealizationMapReport(state, options = {}) {
    const tolerance = Math.max(0, Number(options.tolerance ?? 1e-9));
    const contactStiffness = Math.max(0, Number(options.contactStiffness ?? 1));
    const sequence = Array.isArray(options.eventSequence) ? options.eventSequence : inferRealizationEvents(state, tolerance);
    let current = cloneForCharacterization(state);
    const records = [];
    const supportByKey = new Map();
    for (const [index, event] of sequence.entries()) {
      const { valid, invalid } = realizationEventSupport(current, event, tolerance);
      for (const cell of valid) supportByKey.set(contactGraphKey(cell.row, cell.col), cell);
      let error = null;
      let after = current;
      try {
        if (typeof RAD.applyProgrammableEvent !== "function") throw new Error("RAD.applyProgrammableEvent unavailable");
        after = RAD.applyProgrammableEvent(current, event);
      } catch (err) {
        error = err?.message || String(err);
      }
      records.push(realizationEventRecord(event, index, current, after, valid, invalid, error, tolerance));
      current = after;
    }
    const supportCells = Array.from(supportByKey.values());
    const contactGraph =
      options.contactGraphReport ||
      contactGraphConsistencyReport(current, {
        contactReport: options.contactReport,
        groupSupport: supportCells,
        contactStiffness,
        tolerance,
      });
    const realizedOperatorCount = records.filter((record) => record.realized).length;
    const supportRecordCount = records.filter((record) => record.supportCells.length > 0).length;
    const stateEffectRecordCount = records.filter((record) => record.stateEffect.active).length;
    const claimLabelCount = records.filter((record) => record.claimLabel).length;
    const missingEvidence = [];
    if (!records.length) missingEvidence.push("abstractOperators");
    if (realizedOperatorCount < records.length) missingEvidence.push("realizationRecords");
    if (supportRecordCount < records.length) missingEvidence.push("supportRecords");
    if (stateEffectRecordCount < records.length) missingEvidence.push("stateEffectRecords");
    if (claimLabelCount < records.length) missingEvidence.push("claimLabels");
    if (!contactGraph?.summary?.contactGraphConsistent) missingEvidence.push("contactGraphConsistency");
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.physical-realization-map.v1",
      method: "finite map from abstract programmable discontinuity events to simulator state effects, support records, mechanism channels, and contact-graph consistency evidence",
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        backlash: state.grid.backlash,
        pinHoleClearance: typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : 0,
      },
      operators: records,
      support: { supportCells, supportCellCount: supportCells.length },
      contactGraph: {
        attached: Boolean(contactGraph && typeof contactGraph === "object"),
        schema: contactGraph?.schema || "",
        status: contactGraph?.summary?.status || "",
        contactGraphConsistent: Boolean(contactGraph?.summary?.contactGraphConsistent),
      },
      summary: {
        status: ready ? "physical-realization-map-ready" : "needs-realization-map-evidence",
        physicalRealizationMapReady: ready,
        abstractOperatorCount: records.length,
        realizedOperatorCount,
        supportRecordCount,
        stateEffectRecordCount,
        claimLabelCount,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "physical_realization_map_gate",
        leanStructure: "Mechanics.PhysicalRealizationMapNat",
        leanPredicate: "physicalRealizationMapReadyNat",
        schema: "rad-sim.physical-realization-map.v1",
      },
      claimLabels: {
        map: "simulator-derived physical realization map",
        operatorSemantics: "Lean-proven finite operator/event scaffold",
        contactGraph: "Lean-proven graph/contact bookkeeping scaffold",
        hardwareRealization: "experimentally unvalidated physical assumption",
      },
      limitations: [
        "The map proves finite bookkeeping, not that a real actuator produces the same state change.",
        "Mechanism channels are simulator labels until hardware experiments identify them.",
        "Contact-graph consistency is necessary evidence, not physical contact validation.",
      ],
    };
  }

  function exportPhysicalRealizationMap(state, options = {}) {
    return JSON.stringify(physicalRealizationMapReport(state, options), null, 2);
  }

  function exportPhysicalRealizationMapCsv(reportOrState, options = {}) {
    const report =
      reportOrState?.schema === "rad-sim.physical-realization-map.v1"
        ? reportOrState
        : reportOrState?.experiment?.physicalRealizationMap || physicalRealizationMapReport(reportOrState, options);
    const rows = [[
      "index",
      "kind",
      "operator_class",
      "hardware_channel",
      "realized",
      "support_cells",
      "removed_support_cells",
      "command_alpha_delta",
      "command_z_delta",
      "lock_delta_count",
      "removed_delta_count",
      "evidence_gaps",
    ]];
    for (const record of report.operators || []) {
      rows.push([
        record.index,
        record.kind,
        record.operatorClass,
        record.hardwareChannel,
        record.realized,
        (record.supportCells || []).map((cell) => `${cell.row}:${cell.col}`).join(";"),
        (record.removedSupportCells || []).map((cell) => `${cell.row}:${cell.col}`).join(";"),
        record.stateEffect?.commandAlphaDelta ?? "",
        record.stateEffect?.commandZDelta ?? "",
        record.stateEffect?.lockDeltaCount ?? "",
        record.stateEffect?.removedDeltaCount ?? "",
        Array.isArray(record.evidenceGaps) ? record.evidenceGaps.join(";") : "",
      ]);
    }
    return rows.map((row) => row.map(csvValue).join(",")).join("\n");
  }

  const EXTERNAL_ENGINE_CANDIDATES = [
    {
      engine: "MuJoCo",
      pythonPackage: "mujoco",
      role: "primary rigid-body contact and gravity validation candidate",
      supportsRigidBody: true,
      supportsContact: true,
      supportsGravity: true,
      supportsJointConstraints: true,
      supportsHeadless: true,
      independentFromRadSolver: true,
      notes: "Good first target for pin-hole clearance and fixed-cell gravity tests.",
    },
    {
      engine: "PyBullet",
      pythonPackage: "pybullet",
      role: "secondary rigid-body dynamics and fast prototyping candidate",
      supportsRigidBody: true,
      supportsContact: true,
      supportsGravity: true,
      supportsJointConstraints: true,
      supportsHeadless: true,
      independentFromRadSolver: true,
      notes: "Useful for quick independent checks, but contact tuning must be documented.",
    },
    {
      engine: "Project Chrono",
      pythonPackage: "pychrono",
      role: "higher-fidelity multibody/contact candidate",
      supportsRigidBody: true,
      supportsContact: true,
      supportsGravity: true,
      supportsJointConstraints: true,
      supportsHeadless: true,
      independentFromRadSolver: true,
      notes: "Potential later target for detailed rigid-body and contact studies.",
    },
  ];

  function externalEngineAvailable(candidate, availability) {
    return Boolean(availability?.[candidate.pythonPackage] || availability?.[candidate.engine]);
  }

  function externalPhysicsEngineAuditReport(state, options = {}) {
    const requiredFeatures = Array.from(new Set(options.requiredFeatures || [
      "supportsRigidBody",
      "supportsContact",
      "supportsGravity",
      "supportsJointConstraints",
      "supportsHeadless",
    ]));
    const availability = options.engineAvailability || state.experiment?.externalPhysicsEngineAvailability || {};
    const contact =
      options.contactReport ||
      state.experiment?.contactStateAbstraction ||
      contactStateAbstractionReport(state, options);
    const realization =
      options.physicalRealizationReport ||
      state.experiment?.physicalRealizationMap ||
      physicalRealizationMapReport(state, { ...options, contactReport: contact });

    let availableEngineCount = 0;
    let feasibleEngineCount = 0;
    const engines = EXTERNAL_ENGINE_CANDIDATES.map((candidate) => {
      const available = externalEngineAvailable(candidate, availability);
      const missingFeatures = requiredFeatures.filter((feature) => !candidate[feature]);
      const feasible = available && missingFeatures.length === 0 && Boolean(candidate.independentFromRadSolver);
      if (available) availableEngineCount += 1;
      if (feasible) feasibleEngineCount += 1;
      return {
        ...candidate,
        available,
        requiredFeatures,
        missingFeatures,
        feasible,
      };
    });

    const contactSummary = contact?.summary || {};
    const contactTopology = contact?.topology || {};
    const realizationSummary = realization?.summary || {};
    const activeBodyCount = Math.max(0, Math.round(Number(contactTopology.activeBodyCount || 0)));
    const contactStateRecordCount = Math.max(0, Math.round(Number(contactTopology.contactStateRecordCount || 0)));
    const penaltyTermCount = Math.max(0, Math.round(Number(contactTopology.penaltyTermCount || 0)));
    const abstractOperatorCount = Math.max(0, Math.round(Number(realizationSummary.abstractOperatorCount || 0)));
    const scenarioRecordCount = (activeBodyCount > 0 ? 1 : 0) + (abstractOperatorCount > 0 ? 1 : 0);

    const missingEvidence = [];
    if (!engines.length) missingEvidence.push("engineCandidates");
    if (availableEngineCount === 0) missingEvidence.push("availableExternalEngine");
    if (requiredFeatures.length === 0) missingEvidence.push("requiredFeatures");
    if (scenarioRecordCount === 0) missingEvidence.push("validationScenario");
    if (!contactSummary.contactStateAbstractionReady) missingEvidence.push("contactStateAbstraction");
    if (contactStateRecordCount === 0) missingEvidence.push("contactStateRecords");
    if (penaltyTermCount === 0) missingEvidence.push("contactPenaltyTerms");
    if (feasibleEngineCount === 0) missingEvidence.push("independentToolRecords");

    const ready = missingEvidence.length === 0;
    const report = {
      schema: "rad-sim.external-physics-engine-audit.v1",
      method: "finite readiness audit for independent rigid-body/contact engines before external validation is claimed",
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        backlash: state.grid.backlash,
        pinHoleClearance: typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : 0,
      },
      requiredFeatures,
      engines,
      scenario: {
        activeBodyCount,
        contactStateRecordCount,
        penaltyTermCount,
        abstractOperatorCount,
        scenarioRecordCount,
        contactSchema: contact?.schema || "",
        realizationSchema: realization?.schema || "",
      },
      summary: {
        status: ready ? "external-physics-engine-audit-ready" : "needs-external-physics-evidence",
        externalPhysicsEngineAuditReady: ready,
        engineCandidateCount: engines.length,
        availableEngineCount,
        feasibleEngineCount,
        requiredFeatureRecordCount: requiredFeatures.length,
        scenarioRecordCount,
        contactModelRecordCount: contactStateRecordCount + penaltyTermCount,
        independentToolRecordCount: feasibleEngineCount,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "external_physics_engine_audit_gate",
        leanStructure: "Mechanics.ExternalPhysicsEngineAuditNat",
        leanPredicate: "externalPhysicsEngineAuditReadyNat",
        schema: "rad-sim.external-physics-engine-audit.v1",
      },
      claimLabels: {
        audit: "Lean-proven finite external-engine readiness predicate",
        contactModel: "simulator-derived contact abstraction requiring external validation",
        externalTool: "independent package availability and feature audit, not a completed physics run",
        physicalAccuracy: "experimentally unvalidated physical assumption until bench/external comparisons pass",
      },
      limitations: [
        "The browser cannot inspect installed Python packages; pass engineAvailability explicitly or import a saved audit.",
        "This audit does not run MuJoCo, PyBullet, or Chrono.",
        "A ready audit only proves required evidence fields exist before external validation work begins.",
      ],
    };
    if (options.store === true) {
      state.experiment.externalPhysicsEngineAudit = JSON.parse(JSON.stringify(report));
    }
    return report;
  }

  function exportExternalPhysicsEngineAudit(state, options = {}) {
    return JSON.stringify(externalPhysicsEngineAuditReport(state, options), null, 2);
  }

  function exportExternalPhysicsEngineAuditCsv(reportOrState, options = {}) {
    const report =
      reportOrState?.schema === "rad-sim.external-physics-engine-audit.v1"
        ? reportOrState
        : reportOrState?.experiment?.externalPhysicsEngineAudit ||
          externalPhysicsEngineAuditReport(reportOrState, options);
    const rows = [[
      "engine",
      "python_package",
      "available",
      "feasible",
      "independent_from_rad_solver",
      "missing_features",
      "role",
      "notes",
    ]];
    for (const engine of report.engines || []) {
      rows.push([
        engine.engine,
        engine.pythonPackage,
        engine.available,
        engine.feasible,
        engine.independentFromRadSolver,
        Array.isArray(engine.missingFeatures) ? engine.missingFeatures.join(";") : "",
        engine.role,
        engine.notes,
      ]);
    }
    return rows.map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function mujocoCellName(row, col) {
    return `rad_cell_${Number(row)}_${Number(col)}`;
  }

  function vec3(values) {
    const source = Array.isArray(values) ? values : [];
    return [
      Number(source[0] || 0),
      Number(source[1] || 0),
      Number(source[2] || 0),
    ];
  }

  function xmlVec(values) {
    return vec3(values).map((value) => Number(value).toPrecision(9).replace(/\.?0+$/, "")).join(" ");
  }

  function normalizeFixedCells(state, fixedCells) {
    const source = Array.isArray(fixedCells) && fixedCells.length ? fixedCells : [{ row: 0, col: 0 }];
    const out = new Set();
    for (const raw of source) {
      const cell = normalizeContactGraphCell(raw);
      if (!cell) continue;
      if (cell.row < 0 || cell.row >= state.grid.rows || cell.col < 0 || cell.col >= state.grid.cols) continue;
      out.add(contactGraphKey(cell.row, cell.col));
    }
    return out;
  }

  function normalizeExternalForces(state, externalForces) {
    const out = new Map();
    if (!externalForces || typeof externalForces !== "object") return out;
    const entries = Array.isArray(externalForces)
      ? externalForces.map((entry) => [entry.cell || entry, entry.force || entry.externalForce || [0, 0, 0]])
      : Object.entries(externalForces).map(([key, value]) => {
          const parts = key.split(/[,:]/).map((item) => Number(item));
          return [{ row: parts[0], col: parts[1] }, value];
        });
    for (const [rawCell, rawForce] of entries) {
      const cell = normalizeContactGraphCell(rawCell);
      if (!cell) continue;
      if (cell.row < 0 || cell.row >= state.grid.rows || cell.col < 0 || cell.col >= state.grid.cols) continue;
      out.set(contactGraphKey(cell.row, cell.col), vec3(rawForce));
    }
    return out;
  }

  function mujocoExportState(state, eventSequence) {
    const base = cloneForCharacterization(state);
    if (!Array.isArray(eventSequence) || !eventSequence.length || typeof RAD.applyEventSequence !== "function") {
      return base;
    }
    return RAD.applyEventSequence(base, eventSequence);
  }

  function simulationCenterVector(sim, state, row, col) {
    const center = sim.centers?.[row]?.[col] || {};
    return [
      Number(center.x ?? center[0] ?? col * state.grid.cellSize),
      Number(center.y ?? center[1] ?? row * state.grid.cellSize),
      Number(center.z ?? center[2] ?? sim.height?.[row]?.[col] ?? state.cells.z?.[row]?.[col] ?? 0),
    ];
  }

  function contactDefaults(options = {}) {
    const asArray = (value, fallback) => {
      const source = Array.isArray(value) ? value : fallback;
      return fallback.map((defaultValue, index) => Number(source[index] ?? defaultValue));
    };
    return {
      contactStiffness: Math.max(0, Number(options.contactStiffness ?? 1000)),
      contactDamping: Math.max(0, Number(options.contactDamping ?? 2)),
      friction: asArray(options.friction, [0.4, 0.02, 0.001]),
      solref: asArray(options.solref, [0.02, 1]),
      solimp: asArray(options.solimp, [0.9, 0.95, 0.001, 0.5, 2]),
      margin: Math.max(0, Number(options.contactMargin ?? options.margin ?? 0)),
      gap: Math.max(0, Number(options.contactGap ?? options.gap ?? 0)),
      condim: Math.max(1, Math.round(Number(options.condim ?? 3))),
      calibrated: Boolean(options.calibratedContact || options.calibrated),
    };
  }

  function contactXmlAttributes(defaults) {
    return ` condim="${defaults.condim}" friction="${xmlVec(defaults.friction)}" solref="${xmlVec(defaults.solref)}" solimp="${xmlVec(defaults.solimp)}" margin="${defaults.margin}" gap="${defaults.gap}"`;
  }

  function contactParameterProfileFromPairs(contactPairs, defaults) {
    const xmlAttributes = contactXmlAttributes(defaults);
    const parameters = (contactPairs || []).map((pair) => ({
      row: pair.row,
      col: pair.col,
      pivot: pair.pivot,
      pin: pair.pin,
      hole: pair.hole,
      contactStiffness: defaults.contactStiffness,
      contactDamping: defaults.contactDamping,
      friction: defaults.friction,
      solref: defaults.solref,
      solimp: defaults.solimp,
      margin: defaults.margin,
      gap: defaults.gap,
      condim: defaults.condim,
      calibrated: defaults.calibrated,
      xmlAttributes,
    }));
    const contactPairRecordCount = (contactPairs || []).length;
    const frictionRecordCount = parameters.filter((item) => item.friction.length).length;
    const solverParameterRecordCount = parameters.filter((item) => item.solref.length && item.solimp.length).length;
    const stiffnessRecordCount = parameters.filter((item) => item.contactStiffness > 0).length;
    const dampingRecordCount = parameters.filter((item) => item.contactDamping > 0).length;
    const calibratedRecordCount = parameters.filter((item) => item.calibrated).length;
    const missingEvidence = [];
    if (contactPairRecordCount <= 0) missingEvidence.push("contactPairRecords");
    if (parameters.length < contactPairRecordCount) missingEvidence.push("parameterRecords");
    if (frictionRecordCount < contactPairRecordCount) missingEvidence.push("frictionRecords");
    if (solverParameterRecordCount < contactPairRecordCount) missingEvidence.push("solverParameterRecords");
    if (stiffnessRecordCount < contactPairRecordCount) missingEvidence.push("stiffnessRecords");
    if (dampingRecordCount < contactPairRecordCount) missingEvidence.push("dampingRecords");
    if (!xmlAttributes) missingEvidence.push("xmlContactAttributes");
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.mujoco-contact-parameter-profile.v1",
      method: "finite contact-parameter profile for MuJoCo proxy contact records; parameters are exported but remain uncalibrated unless calibrated=true",
      defaults,
      parameters,
      xmlAttributes,
      summary: {
        status: ready ? "mujoco-contact-parameter-profile-ready" : "needs-mujoco-contact-parameters",
        mujocoContactParameterProfileReady: ready,
        contactPairRecordCount,
        parameterRecordCount: parameters.length,
        frictionRecordCount,
        solverParameterRecordCount,
        stiffnessRecordCount,
        dampingRecordCount,
        calibratedRecordCount,
        xmlAttributeByteCount: xmlAttributes.length,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "mujoco_contact_parameter_profile_gate",
        leanStructure: "Mechanics.ExternalContactParameterProfileNat",
        leanPredicate: "externalContactParameterProfileReadyNat",
        schema: "rad-sim.mujoco-contact-parameter-profile.v1",
      },
      claimLabels: {
        parameters: "simulator-derived external contact parameter profile",
        calibration: defaults.calibrated ? "calibrated contact parameter" : "experimentally unvalidated contact parameter assumption",
        externalEngine: "MuJoCo XML contact attributes for proxy geometry",
      },
      limitations: [
        "Contact stiffness and damping are proxy metadata, not derived from MuJoCo constitutive laws.",
        "Friction, solref, and solimp values require bench calibration before physical claims.",
        "A complete parameter profile does not prove external-run or hardware agreement.",
      ],
    };
  }

  function mujocoContactGeometryFromBodies(state, bodies, options = {}) {
    const defaults = contactDefaults(options);
    const contactAttributes = contactXmlAttributes(defaults);
    const pinRadius = Math.max(1e-9, Number(state.grid.pinRadius || 0));
    const holeRadius = Math.max(pinRadius, Number(state.grid.holeRadius || pinRadius));
    const clearance = holeRadius - pinRadius;
    const plateThickness = Math.max(1e-9, Number(options.plateThickness ?? 0.04));
    const halfHeight = plateThickness * 0.55;
    const pins = [];
    const holes = [];
    const clearanceRecords = [];
    const contactPairs = [];
    const bodyGeomXml = {};
    for (const body of bodies || []) {
      const name = body.name || "";
      if (!name) continue;
      const center = vec3(body.position || [0, 0, 0]);
      const alpha = Math.max(Number(body.alpha || state.grid.initialAlpha || 1), 0.1);
      const pivotSpan = Math.max(state.grid.cellSize * 0.18, 0.32 * state.grid.cellSize * alpha);
      const pivots = [
        ["nw", [-pivotSpan, -pivotSpan, 0]],
        ["ne", [pivotSpan, -pivotSpan, 0]],
        ["se", [pivotSpan, pivotSpan, 0]],
        ["sw", [-pivotSpan, pivotSpan, 0]],
      ];
      bodyGeomXml[name] = [];
      for (const [pivot, local] of pivots) {
        const pin = `${name}_${pivot}_pin`;
        const hole = `${name}_${pivot}_hole`;
        const position = [center[0] + local[0], center[1] + local[1], center[2] + local[2]];
        pins.push({ row: body.row, col: body.col, pivot, name: pin, radius: pinRadius, halfHeight, position });
        holes.push({ row: body.row, col: body.col, pivot, name: hole, radius: holeRadius, halfHeight, position });
        clearanceRecords.push({
          row: body.row,
          col: body.col,
          pivot,
          pinRadius,
          holeRadius,
          clearance,
          clearanceRatio: holeRadius > 0 ? clearance / holeRadius : 0,
        });
        contactPairs.push({ row: body.row, col: body.col, pivot, pin, hole, clearance, active: clearance >= 0 });
        bodyGeomXml[name].push(`      <geom name="${pin}" type="cylinder" pos="${xmlVec(local)}" size="${pinRadius} ${halfHeight}"${contactAttributes} rgba="0.05 0.05 0.05 1"/>`);
        bodyGeomXml[name].push(`      <geom name="${hole}_clearance" type="cylinder" pos="${xmlVec(local)}" size="${holeRadius} ${halfHeight * 1.05}" contype="0" conaffinity="0" rgba="0.1 0.7 0.9 0.18"/>`);
      }
    }
    const xmlFragment = Object.values(bodyGeomXml).flat().join("\n");
    const contactParameterProfile = contactParameterProfileFromPairs(contactPairs, defaults);
    const missingEvidence = [];
    if (!pins.length) missingEvidence.push("pinRecords");
    if (!holes.length) missingEvidence.push("holeRecords");
    if (!clearanceRecords.length) missingEvidence.push("clearanceRecords");
    if (!contactPairs.length) missingEvidence.push("contactPairRecords");
    if (clearance < 0) missingEvidence.push("nonnegativeClearance");
    if (!xmlFragment) missingEvidence.push("mjcfContactGeometry");
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.mujoco-pin-hole-contact-geometry.v1",
      method: "finite pin-hole contact-geometry inventory for the MuJoCo RAD handoff; records geometry proxies but does not calibrate contact",
      pins,
      holes,
      clearanceRecords,
      contactPairs,
      contactParameterProfile,
      bodyGeomXml,
      xmlFragment,
      summary: {
        status: ready ? "mujoco-pin-hole-contact-geometry-ready" : "needs-mujoco-pin-hole-contact-geometry",
        mujocoPinHoleContactGeometryReady: ready,
        pinRecordCount: pins.length,
        holeRecordCount: holes.length,
        clearanceRecordCount: clearanceRecords.length,
        contactPairRecordCount: contactPairs.length,
        activeContactPairCount: contactPairs.filter((pair) => pair.active).length,
        contactParameterRecordCount: contactParameterProfile.summary.parameterRecordCount,
        frictionRecordCount: contactParameterProfile.summary.frictionRecordCount,
        solverParameterRecordCount: contactParameterProfile.summary.solverParameterRecordCount,
        minClearance: clearanceRecords.length ? Math.min(...clearanceRecords.map((record) => record.clearance)) : 0,
        maxClearance: clearanceRecords.length ? Math.max(...clearanceRecords.map((record) => record.clearance)) : 0,
        xmlFragmentByteCount: xmlFragment.length,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "mujoco_pin_hole_contact_geometry_gate",
        leanStructure: "Mechanics.ExternalContactGeometryNat",
        leanPredicate: "externalContactGeometryReadyNat",
        schema: "rad-sim.mujoco-pin-hole-contact-geometry.v1",
      },
      claimLabels: {
        geometry: "simulator-derived pin-hole proxy geometry",
        clearance: "configured geometric clearance record",
        contact: "external-engine contact candidate; stiffness/friction remain uncalibrated",
      },
      limitations: [
        "Pin and clearance cylinders are proxy geometry, not measured CAD surfaces.",
        "The hole clearance shell is non-colliding in the exported MJCF until contact parameters are calibrated.",
        "Friction, compliance, assembly offsets, and wear are not represented.",
      ],
    };
  }

  function mujocoModelExportReport(state, options = {}) {
    const eventSequence = Array.isArray(options.eventSequence) ? options.eventSequence : [];
    const exportState = mujocoExportState(state, eventSequence);
    const sim = RAD.simulate(exportState);
    const gravity = vec3(options.gravity || [0, 0, -9.81]);
    const timestep = Math.max(1e-9, Number(options.timestep ?? 0.002));
    const plateThickness = Math.max(1e-9, Number(options.plateThickness ?? 0.04));
    const bodyDensity = Math.max(1e-9, Number(options.bodyDensity ?? 900));
    const fixed = normalizeFixedCells(exportState, options.fixedCells || options.loadCase?.fixedCells);
    const forces = normalizeExternalForces(exportState, options.externalForces || options.loadCase?.externalForces);
    const bodies = [];
    const removedCells = [];
    const halfThickness = plateThickness / 2;
    const minHalf = Math.max(exportState.grid.cellSize * 0.05, 1e-6);
    for (let row = 0; row < exportState.grid.rows; row += 1) {
      for (let col = 0; col < exportState.grid.cols; col += 1) {
        if (exportState.cells.removed?.[row]?.[col]) {
          removedCells.push({ row, col });
          continue;
        }
        const center = simulationCenterVector(sim, exportState, row, col);
        const alpha = Number(sim.alpha?.[row]?.[col] ?? exportState.cells.alpha?.[row]?.[col] ?? exportState.grid.initialAlpha);
        const halfSide = Math.max(minHalf, 0.24 * exportState.grid.cellSize * Math.max(alpha, 0.1));
        const key = contactGraphKey(row, col);
        bodies.push({
          row,
          col,
          name: mujocoCellName(row, col),
          position: center,
          halfExtents: [halfSide, halfSide, halfThickness],
          alpha,
          thetaDegrees: Number(sim.theta?.[row]?.[col] ?? exportState.cells.theta?.[row]?.[col] ?? 0),
          commandAlpha: Number(exportState.cells.commandAlpha?.[row]?.[col] || 0),
          commandZ: Number(exportState.cells.commandZ?.[row]?.[col] || 0),
          locked: Boolean(exportState.cells.locked?.[row]?.[col]),
          fixed: fixed.has(key),
          hasExternalForce: forces.has(key),
        });
      }
    }
    const contactGeometry = mujocoContactGeometryFromBodies(exportState, bodies, { ...options, plateThickness });
    const size = Math.max(exportState.grid.rows, exportState.grid.cols, 1) * exportState.grid.cellSize * 2;
    const xml = [
      '<mujoco model="rad_lattice_external_validation">',
      '  <compiler angle="radian"/>',
      `  <option timestep="${timestep}" gravity="${xmlVec(gravity)}"/>`,
      "  <worldbody>",
      `    <geom name="floor" type="plane" size="${size} ${size} 0.05" rgba="0.65 0.65 0.65 0.25"/>`,
      ...bodies.flatMap((body) => {
        const lines = [`    <body name="${body.name}" pos="${xmlVec(body.position)}">`];
        if (!body.fixed) lines.push(`      <joint name="${body.name}_free" type="free" damping="0.05"/>`);
        lines.push(`      <geom name="${body.name}_plate" type="box" size="${xmlVec(body.halfExtents)}" density="${bodyDensity}" rgba="${body.fixed ? "0.20 0.45 0.85 1" : "0.85 0.55 0.18 1"}"/>`);
        for (const line of contactGeometry.bodyGeomXml?.[body.name] || []) lines.push(line);
        lines.push("    </body>");
        return lines;
      }),
      "  </worldbody>",
      "</mujoco>",
    ].join("\n");
    const loads = [];
    for (const [key, force] of forces.entries()) {
      const [row, col] = key.split(":").map((value) => Number(value));
      if (exportState.cells.removed?.[row]?.[col]) continue;
      loads.push({ row, col, force });
    }
    const fixedBodyCount = bodies.filter((body) => body.fixed).length;
    const gravityEnabled = gravity.some((value) => Math.abs(value) > 0);
    const missingEvidence = [];
    if (!bodies.length) missingEvidence.push("bodyRecords");
    if (fixedBodyCount <= 0) missingEvidence.push("fixedBodyRecords");
    if (!gravityEnabled) missingEvidence.push("gravityRecord");
    if (!xml.trim()) missingEvidence.push("mjcfXml");
    const ready = missingEvidence.length === 0;
    const report = {
      schema: "rad-sim.mujoco-model-export.v1",
      method: "coarse MJCF export for independent rigid-body/contact validation; not fabrication-accurate CAD",
      engine: "MuJoCo",
      xml,
      grid: {
        rows: exportState.grid.rows,
        cols: exportState.grid.cols,
        cellSize: exportState.grid.cellSize,
        backlash: exportState.grid.backlash,
        pinHoleClearance: typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(exportState) : 0,
      },
      settings: { gravity, timestep, plateThickness, bodyDensity, contactParameters: contactGeometry.contactParameterProfile.defaults },
      bodies,
      removedCells,
      contactGeometry,
      loads,
      summary: {
        status: ready ? "mujoco-model-export-ready" : "needs-mujoco-export-evidence",
        mujocoModelExportReady: ready,
        bodyRecordCount: bodies.length,
        fixedBodyCount,
        removedBodyCount: removedCells.length,
        loadRecordCount: loads.length,
        pinRecordCount: contactGeometry.summary.pinRecordCount,
        holeRecordCount: contactGeometry.summary.holeRecordCount,
        clearanceRecordCount: contactGeometry.summary.clearanceRecordCount,
        contactPairRecordCount: contactGeometry.summary.contactPairRecordCount,
        contactParameterRecordCount: contactGeometry.summary.contactParameterRecordCount,
        frictionRecordCount: contactGeometry.summary.frictionRecordCount,
        solverParameterRecordCount: contactGeometry.summary.solverParameterRecordCount,
        gravityRecordCount: gravityEnabled ? 1 : 0,
        xmlByteCount: xml.length,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "mujoco_model_export_gate",
        leanStructure: "Mechanics.ExternalPhysicsModelExportNat",
        leanPredicate: "externalPhysicsModelExportReadyNat",
        schema: "rad-sim.mujoco-model-export.v1",
      },
      claimLabels: {
        export: "simulator-derived coarse external-engine model export",
        geometry: "normalized proxy geometry, not fabrication-accurate CAD",
        physics: "experimentally unvalidated until MuJoCo run and bench comparison pass",
      },
      limitations: [
        "Each cell is represented as one coarse box body rather than the full RAD pin/plate assembly.",
        "The browser can export and compare MuJoCo artifacts but cannot execute MuJoCo.",
        "Pin-hole clearance is recorded in metadata but not yet represented as exact contact geometry in MJCF.",
      ],
    };
    if (options.store === true) state.experiment.mujocoModelExport = JSON.parse(JSON.stringify(report));
    return report;
  }

  function exportMujocoModelXml(state, options = {}) {
    return mujocoModelExportReport(state, options).xml;
  }

  function exportMujocoModelReport(state, options = {}) {
    return JSON.stringify(mujocoModelExportReport(state, options), null, 2);
  }

  function mujocoPinHoleContactGeometryReport(state, options = {}) {
    const exportReport = options.exportReport || state.experiment?.mujocoModelExport;
    if (exportReport?.contactGeometry) return exportReport.contactGeometry;
    if (Array.isArray(exportReport?.bodies)) {
      return mujocoContactGeometryFromBodies(state, exportReport.bodies, options);
    }
    return mujocoModelExportReport(state, options).contactGeometry;
  }

  function exportMujocoPinHoleContactGeometry(state, options = {}) {
    return JSON.stringify(mujocoPinHoleContactGeometryReport(state, options), null, 2);
  }

  function exportMujocoPinHoleContactGeometryCsv(report) {
    const rows = [["row", "col", "pivot", "pin", "hole", "pin_radius", "hole_radius", "clearance", "active"]];
    const pairByKey = new Map();
    for (const pair of report.contactPairs || []) {
      pairByKey.set(`${pair.row}:${pair.col}:${pair.pivot}`, pair);
    }
    for (const record of report.clearanceRecords || []) {
      const pair = pairByKey.get(`${record.row}:${record.col}:${record.pivot}`) || {};
      rows.push([
        record.row,
        record.col,
        record.pivot,
        pair.pin || "",
        pair.hole || "",
        record.pinRadius,
        record.holeRadius,
        record.clearance,
        pair.active ?? "",
      ]);
    }
    rows.push([
      "summary",
      "",
      "",
      report.summary?.pinRecordCount ?? "",
      report.summary?.holeRecordCount ?? "",
      "",
      "",
      report.summary?.minClearance ?? "",
      report.summary?.mujocoPinHoleContactGeometryReady ?? "",
    ]);
    return rows.map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function mujocoContactParameterReport(state, options = {}) {
    const geometry =
      options.contactGeometryReport ||
      mujocoPinHoleContactGeometryReport(state, options);
    if (
      geometry?.contactParameterProfile &&
      options.contactStiffness === undefined &&
      options.contactDamping === undefined &&
      options.friction === undefined &&
      options.solref === undefined &&
      options.solimp === undefined &&
      options.contactMargin === undefined &&
      options.contactGap === undefined &&
      options.condim === undefined &&
      options.calibratedContact === undefined
    ) {
      return geometry.contactParameterProfile;
    }
    return contactParameterProfileFromPairs(geometry?.contactPairs || [], contactDefaults(options));
  }

  function exportMujocoContactParameter(state, options = {}) {
    return JSON.stringify(mujocoContactParameterReport(state, options), null, 2);
  }

  function exportMujocoContactParameterCsv(report) {
    const rows = [["row", "col", "pivot", "pin", "hole", "contact_stiffness", "contact_damping", "friction", "solref", "solimp", "condim", "calibrated"]];
    for (const parameter of report.parameters || []) {
      rows.push([
        parameter.row,
        parameter.col,
        parameter.pivot,
        parameter.pin,
        parameter.hole,
        parameter.contactStiffness,
        parameter.contactDamping,
        Array.isArray(parameter.friction) ? parameter.friction.join(";") : "",
        Array.isArray(parameter.solref) ? parameter.solref.join(";") : "",
        Array.isArray(parameter.solimp) ? parameter.solimp.join(";") : "",
        parameter.condim,
        parameter.calibrated,
      ]);
    }
    rows.push([
      "summary",
      "",
      "",
      "",
      "",
      report.summary?.stiffnessRecordCount ?? "",
      report.summary?.dampingRecordCount ?? "",
      report.summary?.frictionRecordCount ?? "",
      report.summary?.solverParameterRecordCount ?? "",
      "",
      "",
      report.summary?.calibratedRecordCount ?? "",
    ]);
    return rows.map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function contactParameterMeasurementColumns() {
    return [
      "dataset_id",
      "dataset_role",
      "scenario_id",
      "repeat_index",
      "row",
      "col",
      "pivot",
      "pin_radius_mm",
      "hole_radius_mm",
      "clearance_mm",
      "normal_load_n",
      "tangential_load_n",
      "imposed_z_mm",
      "measured_pin_hole_slip_mm",
      "measured_normal_force_n",
      "measured_tangent_force_n",
      "measured_rebound_ratio",
      "measured_contact_duration_s",
      "measured_static_friction_coeff",
      "measured_dynamic_friction_coeff",
      "fitted_contact_stiffness",
      "fitted_contact_damping",
      "fitted_solref_timeconst",
      "fitted_solref_damping_ratio",
      "fitted_solimp_width",
      "fixture_notes",
    ];
  }

  function contactParameterTemplateRows(parameterReport, { datasetId, datasetRole, repeatCount }) {
    const rows = [];
    for (const parameter of parameterReport.parameters || []) {
      const row = Number(parameter.row || 0);
      const col = Number(parameter.col || 0);
      const pivot = String(parameter.pivot || "");
      for (let repeatIndex = 1; repeatIndex <= repeatCount; repeatIndex += 1) {
        rows.push({
          datasetId,
          datasetRole,
          scenarioId: `${datasetRole}-r${row}c${col}-${pivot}-contact-repeat-${repeatIndex}`,
          repeatIndex,
          row,
          col,
          pivot,
          pin: parameter.pin || "",
          hole: parameter.hole || "",
          assumedContactStiffness: parameter.contactStiffness ?? null,
          assumedContactDamping: parameter.contactDamping ?? null,
          assumedFriction: Array.isArray(parameter.friction) ? [...parameter.friction] : [],
          assumedSolref: Array.isArray(parameter.solref) ? [...parameter.solref] : [],
          assumedSolimp: Array.isArray(parameter.solimp) ? [...parameter.solimp] : [],
          assumedCondim: parameter.condim ?? null,
          requiredMeasurements: {
            pinRadiusMm: null,
            holeRadiusMm: null,
            clearanceMm: null,
            normalLoadN: null,
            tangentialLoadN: null,
            imposedZMm: null,
            measuredPinHoleSlipMm: null,
            measuredNormalForceN: null,
            measuredTangentForceN: null,
            measuredReboundRatio: null,
            measuredContactDurationS: null,
            measuredStaticFrictionCoeff: null,
            measuredDynamicFrictionCoeff: null,
            fittedContactStiffness: null,
            fittedContactDamping: null,
            fittedSolrefTimeconst: null,
            fittedSolrefDampingRatio: null,
            fittedSolimpWidth: null,
            fixtureNotes: "",
          },
        });
      }
    }
    return rows;
  }

  function contactParameterCalibrationPacket(state, options = {}) {
    const repeatCount = Math.max(1, Math.floor(Number(options.repeatCount ?? 3)));
    const fitDatasetId = options.fitDatasetId || "contact-fit-run-001";
    const holdoutDatasetId = options.holdoutDatasetId || "contact-holdout-run-001";
    const profileId = options.profileId || "contact-parameter-profile-v1";
    const profileFrozenAt = options.profileFrozenAt || "set-after-contact-fit-before-holdout";
    const geometry =
      options.contactGeometryReport ||
      state.experiment?.mujocoPinHoleContactGeometry ||
      mujocoPinHoleContactGeometryReport(state, options);
    const parameters =
      options.contactParameterReport ||
      state.experiment?.mujocoContactParameterProfile ||
      mujocoContactParameterReport(state, { ...options, contactGeometryReport: geometry });
    const summary = parameters.summary || {};
    const contactPairs = Number(summary.contactPairRecordCount || 0);
    const parameterRecords = Number(summary.parameterRecordCount || 0);
    const profileReady = Boolean(summary.mujocoContactParameterProfileReady);
    const measurementColumns = contactParameterMeasurementColumns();
    const fitTemplateRows = contactParameterTemplateRows(parameters, {
      datasetId: fitDatasetId,
      datasetRole: "fit",
      repeatCount,
    });
    const holdoutTemplateRows = contactParameterTemplateRows(parameters, {
      datasetId: holdoutDatasetId,
      datasetRole: "holdout",
      repeatCount,
    });
    const filenames = {
      packet: "contact_parameter_calibration_packet.json",
      csv: "contact_parameter_calibration_template.csv",
      contact_geometry: "mujoco_pin_hole_contact_geometry.json",
      contact_parameters: "mujoco_contact_parameter_profile.json",
      readme: "README.md",
    };
    const artifactManifest = [
      { id: "packet", filename: filenames.packet, schema: "rad-sim.contact-parameter-calibration-packet.v1", purpose: "single JSON bundle for contact-parameter calibration review" },
      { id: "template-csv", filename: filenames.csv, schema: "rad-sim.contact-parameter-calibration-template.v1", purpose: "blank fit and holdout measurement rows for contact pairs" },
      { id: "contact-geometry", filename: filenames.contact_geometry, schema: "rad-sim.mujoco-pin-hole-contact-geometry.v1", purpose: "proxy pin-hole geometry that defines the contact pair inventory" },
      { id: "contact-parameters", filename: filenames.contact_parameters, schema: "rad-sim.mujoco-contact-parameter-profile.v1", purpose: "explicit proxy stiffness, damping, friction, and solver settings" },
      { id: "readme", filename: filenames.readme, schema: "text/markdown", purpose: "human instructions and claim limitations" },
    ];
    const missingEvidence = [];
    if (contactPairs <= 0) missingEvidence.push("contactPairRecords");
    if (parameterRecords < contactPairs) missingEvidence.push("parameterRecords");
    if (!measurementColumns.length) missingEvidence.push("measurementColumns");
    if (fitTemplateRows.length < contactPairs) missingEvidence.push("fitTemplateRows");
    if (holdoutTemplateRows.length < contactPairs) missingEvidence.push("holdoutTemplateRows");
    if (artifactManifest.length < 5) missingEvidence.push("artifactManifest");
    if (!profileReady) missingEvidence.push("contactParameterProfile");
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.contact-parameter-calibration-packet.v1",
      method: "bench handoff packet linking proxy pin-hole geometry, MuJoCo contact parameter assumptions, and blank fit/holdout measurement rows",
      schemas: {
        template: "rad-sim.contact-parameter-calibration-template.v1",
        contactGeometry: "rad-sim.mujoco-pin-hole-contact-geometry.v1",
        contactParameters: "rad-sim.mujoco-contact-parameter-profile.v1",
      },
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        cellSize: state.grid.cellSize,
        backlash: state.grid.backlash,
        pinRadius: state.grid.pinRadius,
        holeRadius: state.grid.holeRadius,
        pinHoleClearance: typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : Math.max(0, (state.grid.holeRadius || 0) - (state.grid.pinRadius || 0)),
      },
      profilePlan: {
        profileId,
        profileFrozenAt,
        fitDatasetId,
        holdoutDatasetId,
        repeatCount,
        freezeRule: "Fit contact parameters from fit rows, freeze the profile ID, then collect holdout rows without modifying the fitted parameters.",
      },
      measurementColumns,
      fitTemplateRows,
      holdoutTemplateRows,
      contactGeometry: geometry,
      contactParameterProfile: parameters,
      filenames,
      artifactManifest,
      validationInstructions: {
        fitStep: "Measure pin radius, hole radius, slip, force, friction, and rebound fields in fit rows.",
        freezeStep: "Estimate stiffness, damping, friction, solref, and solimp once, then freeze the profile.",
        holdoutStep: "Repeat the same contact-pair measurements in a separate holdout dataset.",
        acceptanceRule: "Treat calibratedContact=true as meaningful only after fit and holdout rows are filled, residuals are compared, and missing measurement fields are zero.",
      },
      summary: {
        status: ready ? "contact-parameter-calibration-packet-ready" : "needs-contact-parameter-calibration-evidence",
        contactParameterCalibrationPacketReady: ready,
        contactPairRecordCount: contactPairs,
        parameterRecordCount: parameterRecords,
        measurementColumnCount: measurementColumns.length,
        fitTemplateRowCount: fitTemplateRows.length,
        holdoutTemplateRowCount: holdoutTemplateRows.length,
        artifactManifestEntryCount: artifactManifest.length,
        profileEvidenceReady: profileReady,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "contact_parameter_calibration_packet_completeness",
        leanStructure: "Mechanics.ContactParameterCalibrationPacketNat",
        leanPredicate: "contactParameterCalibrationPacketCompleteNat",
        schema: "rad-sim.contact-parameter-calibration-packet.v1",
      },
      claimLabels: {
        packetAssembly: "bench protocol artifact",
        parameterValues: "experimentally unvalidated physical assumption until fit/holdout measurements pass",
        formalCompleteness: "Lean-proven finite packet-completeness predicate, not calibrated contact mechanics",
      },
      limitations: [
        "The packet creates measurement rows but cannot prove that the measurements were collected.",
        "Fitted stiffness, damping, friction, and solver parameters remain empirical until independently validated.",
        "MuJoCo agreement is not hardware validation without bench measurements and provenance.",
      ],
    };
  }

  function exportContactParameterCalibrationPacket(state, options = {}) {
    return JSON.stringify(contactParameterCalibrationPacket(state, options), null, 2);
  }

  function exportContactParameterCalibrationPacketCsv(packetOrState, options = {}) {
    const packet =
      packetOrState?.schema === "rad-sim.contact-parameter-calibration-packet.v1"
        ? packetOrState
        : contactParameterCalibrationPacket(packetOrState, options);
    const rows = [[
      "dataset_id",
      "dataset_role",
      "scenario_id",
      "repeat_index",
      "row",
      "col",
      "pivot",
      "pin",
      "hole",
      "assumed_contact_stiffness",
      "assumed_contact_damping",
      "assumed_friction",
      "assumed_solref",
      "assumed_solimp",
      "assumed_condim",
      "pin_radius_mm",
      "hole_radius_mm",
      "clearance_mm",
      "normal_load_n",
      "tangential_load_n",
      "imposed_z_mm",
      "measured_pin_hole_slip_mm",
      "measured_normal_force_n",
      "measured_tangent_force_n",
      "measured_static_friction_coeff",
      "measured_dynamic_friction_coeff",
      "fitted_contact_stiffness",
      "fitted_contact_damping",
      "fitted_solref_timeconst",
      "fitted_solref_damping_ratio",
      "fitted_solimp_width",
      "fixture_notes",
    ]];
    for (const row of [...(packet.fitTemplateRows || []), ...(packet.holdoutTemplateRows || [])]) {
      const measured = row.requiredMeasurements || {};
      rows.push([
        row.datasetId || "",
        row.datasetRole || "",
        row.scenarioId || "",
        row.repeatIndex ?? "",
        row.row ?? "",
        row.col ?? "",
        row.pivot || "",
        row.pin || "",
        row.hole || "",
        row.assumedContactStiffness ?? "",
        row.assumedContactDamping ?? "",
        Array.isArray(row.assumedFriction) ? row.assumedFriction.join(";") : "",
        Array.isArray(row.assumedSolref) ? row.assumedSolref.join(";") : "",
        Array.isArray(row.assumedSolimp) ? row.assumedSolimp.join(";") : "",
        row.assumedCondim ?? "",
        measured.pinRadiusMm ?? "",
        measured.holeRadiusMm ?? "",
        measured.clearanceMm ?? "",
        measured.normalLoadN ?? "",
        measured.tangentialLoadN ?? "",
        measured.imposedZMm ?? "",
        measured.measuredPinHoleSlipMm ?? "",
        measured.measuredNormalForceN ?? "",
        measured.measuredTangentForceN ?? "",
        measured.measuredStaticFrictionCoeff ?? "",
        measured.measuredDynamicFrictionCoeff ?? "",
        measured.fittedContactStiffness ?? "",
        measured.fittedContactDamping ?? "",
        measured.fittedSolrefTimeconst ?? "",
        measured.fittedSolrefDampingRatio ?? "",
        measured.fittedSolimpWidth ?? "",
        measured.fixtureNotes ?? "",
      ]);
    }
    rows.push([
      "summary",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      packet.summary?.parameterRecordCount ?? "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      packet.summary?.contactParameterCalibrationPacketReady ?? "",
    ]);
    return rows.map((row) => row.map(csvValue).join(",")).join("\n");
  }

  const CONTACT_PARAMETER_RESULT_FIELDS = [
    "pinRadiusMm",
    "holeRadiusMm",
    "clearanceMm",
    "normalLoadN",
    "tangentialLoadN",
    "imposedZMm",
    "measuredPinHoleSlipMm",
    "measuredNormalForceN",
    "measuredTangentForceN",
    "measuredReboundRatio",
    "measuredContactDurationS",
    "measuredStaticFrictionCoeff",
    "measuredDynamicFrictionCoeff",
    "fittedContactStiffness",
    "fittedContactDamping",
    "fittedSolrefTimeconst",
    "fittedSolrefDampingRatio",
    "fittedSolimpWidth",
  ];

  function cloneJson(value) {
    return JSON.parse(JSON.stringify(value ?? null));
  }

  function contactParameterCalibrationResultsTemplate(packetOrState, options = {}) {
    const packet =
      packetOrState?.schema === "rad-sim.contact-parameter-calibration-packet.v1"
        ? packetOrState
        : contactParameterCalibrationPacket(packetOrState, options);
    const measurements = [
      ...(Array.isArray(packet.fitTemplateRows) ? packet.fitTemplateRows : []),
      ...(Array.isArray(packet.holdoutTemplateRows) ? packet.holdoutTemplateRows : []),
    ].map((row) => cloneJson(row));
    const datasetId = options.datasetId || "contact-parameter-results-001";
    return {
      schema: "rad-sim.contact-parameter-calibration-results.v1",
      sourcePacketSchema: packet.schema || "",
      sourcePacketTargetId: packet.formalization?.targetId || "",
      datasetId,
      datasetRole: options.datasetRole || "fit-holdout-results",
      sourceFileId: options.sourceFileId || `${datasetId}.json`,
      operator: options.operator || "",
      collectedAt: options.collectedAt || "",
      profilePlan: cloneJson(packet.profilePlan || {}),
      measurementColumns: Array.isArray(packet.measurementColumns) ? [...packet.measurementColumns] : [],
      provenanceNotes: "Fill requiredMeasurements from bench or external contact fitting; do not edit assumedContact*, assumedFriction, assumedSolref, or assumedSolimp.",
      measurements,
    };
  }

  function exportContactParameterCalibrationResultsTemplate(packetOrState, options = {}) {
    return JSON.stringify(contactParameterCalibrationResultsTemplate(packetOrState, options), null, 2);
  }

  function contactParameterCalibrationResultsFromJson(text) {
    const raw = JSON.parse(text);
    if (!raw || typeof raw !== "object" || Array.isArray(raw)) {
      throw new Error("contact-parameter calibration results must be a JSON object");
    }
    if (raw.schema !== "rad-sim.contact-parameter-calibration-results.v1") {
      throw new Error("unsupported contact-parameter calibration results schema");
    }
    if (!Array.isArray(raw.measurements)) {
      throw new Error("contact-parameter calibration results require measurements");
    }
    return raw;
  }

  function contactParameterResultNumber(row, field) {
    return numericOrNull(row.requiredMeasurements?.[field]);
  }

  function contactParameterExpectedValue(packet, row, field) {
    const grid = packet.grid || {};
    if (field === "pinRadiusMm") return numericOrNull(grid.pinRadius);
    if (field === "holeRadiusMm") return numericOrNull(grid.holeRadius);
    if (field === "clearanceMm") return numericOrNull(grid.pinHoleClearance);
    if (field === "fittedContactStiffness") return numericOrNull(row.assumedContactStiffness);
    if (field === "fittedContactDamping") return numericOrNull(row.assumedContactDamping);
    if (field === "measuredStaticFrictionCoeff") return numericOrNull(row.assumedFriction?.[0]);
    if (field === "measuredDynamicFrictionCoeff") return numericOrNull(row.assumedFriction?.[1]);
    if (field === "fittedSolrefTimeconst") return numericOrNull(row.assumedSolref?.[0]);
    if (field === "fittedSolrefDampingRatio") return numericOrNull(row.assumedSolref?.[1]);
    if (field === "fittedSolimpWidth") return numericOrNull(row.assumedSolimp?.[0]);
    return null;
  }

  function contactParameterKey(row) {
    return [
      Number(row.row || 0),
      Number(row.col || 0),
      String(row.pivot || ""),
      String(row.pin || ""),
      String(row.hole || ""),
    ].join("|");
  }

  function meanNumber(values) {
    return values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : null;
  }

  function compareContactParameterCalibrationResults(packetOrState, results = null, options = {}) {
    const packet =
      packetOrState?.schema === "rad-sim.contact-parameter-calibration-packet.v1"
        ? packetOrState
        : contactParameterCalibrationPacket(packetOrState, options);
    if (packet.schema !== "rad-sim.contact-parameter-calibration-packet.v1") {
      throw new Error("unsupported contact-parameter calibration packet schema");
    }
    const filledResults = results || contactParameterCalibrationResultsTemplate(packet, options);
    if (filledResults.schema !== "rad-sim.contact-parameter-calibration-results.v1") {
      throw new Error("unsupported contact-parameter calibration results schema");
    }
    const tolerance = Math.max(0, Number(options.tolerance ?? 1e-9));
    const holdoutTolerance = Math.max(0, Number(options.holdoutTolerance ?? 1e-6));
    const expectedFitRows = Array.isArray(packet.fitTemplateRows) ? packet.fitTemplateRows.length : 0;
    const expectedHoldoutRows = Array.isArray(packet.holdoutTemplateRows) ? packet.holdoutTemplateRows.length : 0;
    const packetReady = Boolean(packet.summary?.contactParameterCalibrationPacketReady);
    const resultRows = Array.isArray(filledResults.measurements) ? filledResults.measurements.filter((row) => row && typeof row === "object") : [];
    const fitRows = resultRows.filter((row) => row.datasetRole === "fit");
    const holdoutRows = resultRows.filter((row) => row.datasetRole === "holdout");
    const compareFields = [
      "pinRadiusMm",
      "holeRadiusMm",
      "clearanceMm",
      "measuredStaticFrictionCoeff",
      "measuredDynamicFrictionCoeff",
      "fittedContactStiffness",
      "fittedContactDamping",
      "fittedSolrefTimeconst",
      "fittedSolrefDampingRatio",
      "fittedSolimpWidth",
    ];
    const roleValues = new Map();
    const fitResiduals = [];
    const holdoutResiduals = [];
    const rowReports = [];
    let missingMeasurements = 0;
    let completedMeasurementRows = 0;
    let parameterResidualCount = 0;
    for (const row of resultRows) {
      const missingFields = CONTACT_PARAMETER_RESULT_FIELDS.filter((field) => contactParameterResultNumber(row, field) === null);
      missingMeasurements += missingFields.length;
      if (missingFields.length === 0) completedMeasurementRows += 1;
      const residuals = [];
      let maxRowResidual = 0;
      const role = String(row.datasetRole || "");
      const key = contactParameterKey(row);
      for (const field of compareFields) {
        const observed = contactParameterResultNumber(row, field);
        const expected = contactParameterExpectedValue(packet, row, field);
        if (observed === null || expected === null) continue;
        const residual = Math.abs(observed - expected);
        parameterResidualCount += 1;
        maxRowResidual = Math.max(maxRowResidual, residual);
        const pass = residual <= (role === "holdout" ? holdoutTolerance : tolerance);
        residuals.push({ field, observed, expected, absResidual: residual, pass });
        const mapKey = `${role}|${key}|${field}`;
        if (!roleValues.has(mapKey)) roleValues.set(mapKey, []);
        roleValues.get(mapKey).push(observed);
        if (role === "holdout") holdoutResiduals.push(residual);
        else fitResiduals.push(residual);
      }
      rowReports.push({
        datasetId: row.datasetId || "",
        datasetRole: role,
        scenarioId: row.scenarioId || "",
        repeatIndex: row.repeatIndex ?? "",
        row: row.row ?? "",
        col: row.col ?? "",
        pivot: row.pivot || "",
        pin: row.pin || "",
        hole: row.hole || "",
        missingFields,
        missingMeasurementCount: missingFields.length,
        residuals,
        maxParameterResidual: maxRowResidual,
        pass: missingFields.length === 0 && residuals.every((item) => item.pass),
      });
    }
    const holdoutPairDeltas = [];
    let comparedPairFields = 0;
    for (const [mapKey, values] of roleValues.entries()) {
      if (!mapKey.startsWith("fit|")) continue;
      const holdoutKey = `holdout|${mapKey.slice("fit|".length)}`;
      const fitValue = meanNumber(values);
      const holdoutValue = meanNumber(roleValues.get(holdoutKey) || []);
      if (fitValue === null || holdoutValue === null) continue;
      comparedPairFields += 1;
      holdoutPairDeltas.push(Math.abs(fitValue - holdoutValue));
    }
    const maxFitResidual = fitResiduals.length ? Math.max(...fitResiduals) : 0;
    const maxHoldoutResidual = holdoutResiduals.length ? Math.max(...holdoutResiduals) : 0;
    const maxFitHoldoutDelta = holdoutPairDeltas.length ? Math.max(...holdoutPairDeltas) : 0;
    const fitParameterPass =
      fitRows.length >= expectedFitRows &&
      fitRows.length > 0 &&
      rowReports.filter((row) => row.datasetRole === "fit").every((row) => row.missingMeasurementCount === 0) &&
      fitResiduals.length > 0 &&
      maxFitResidual <= tolerance;
    const holdoutParameterPass =
      holdoutRows.length >= expectedHoldoutRows &&
      holdoutRows.length > 0 &&
      rowReports.filter((row) => row.datasetRole === "holdout").every((row) => row.missingMeasurementCount === 0) &&
      holdoutResiduals.length > 0 &&
      maxHoldoutResidual <= holdoutTolerance;
    const independentHoldoutPass = comparedPairFields > 0 && maxFitHoldoutDelta <= holdoutTolerance;
    const missingEvidence = [];
    if (!packetReady) missingEvidence.push("contactParameterCalibrationPacket");
    if (expectedFitRows <= 0) missingEvidence.push("fitTemplateRows");
    if (expectedHoldoutRows <= 0) missingEvidence.push("holdoutTemplateRows");
    if (!resultRows.length) missingEvidence.push("resultRows");
    if (fitRows.length < expectedFitRows) missingEvidence.push("completedFitRows");
    if (holdoutRows.length < expectedHoldoutRows) missingEvidence.push("completedHoldoutRows");
    if (missingMeasurements > 0) missingEvidence.push("completedMeasurements");
    if (!fitResiduals.length) missingEvidence.push("fitParameterResiduals");
    if (!holdoutResiduals.length) missingEvidence.push("holdoutParameterResiduals");
    if (comparedPairFields <= 0) missingEvidence.push("independentHoldoutPairs");
    const passFlag = missingEvidence.length === 0 && fitParameterPass && holdoutParameterPass && independentHoldoutPass;
    return {
      schema: "rad-sim.contact-parameter-bench-validation.v1",
      method: "comparison of filled fit and holdout contact-parameter measurements against the frozen proxy MuJoCo contact profile",
      sourcePacket: {
        schema: packet.schema || "",
        ready: packetReady,
        profilePlan: cloneJson(packet.profilePlan || {}),
      },
      sourceResults: {
        schema: filledResults.schema || "",
        datasetId: filledResults.datasetId || "",
        datasetRole: filledResults.datasetRole || "",
        sourceFileId: filledResults.sourceFileId || "",
      },
      metrics: {
        resultRowCount: resultRows.length,
        expectedFitRowCount: expectedFitRows,
        expectedHoldoutRowCount: expectedHoldoutRows,
        fitRowCount: fitRows.length,
        holdoutRowCount: holdoutRows.length,
        completedMeasurementRowCount: completedMeasurementRows,
        missingMeasurementCount: missingMeasurements,
        comparedPairFieldCount: comparedPairFields,
        parameterResidualCount,
        fitParameterResidualCount: fitResiduals.length,
        holdoutParameterResidualCount: holdoutResiduals.length,
        maxFitParameterResidual: maxFitResidual,
        maxHoldoutParameterResidual: maxHoldoutResidual,
        maxFitHoldoutParameterDelta: maxFitHoldoutDelta,
        tolerance,
        holdoutTolerance,
      },
      summary: {
        status: passFlag ? "contact-parameter-bench-validation-pass" : "needs-contact-parameter-bench-review",
        contactParameterBenchValidationPass: passFlag,
        fitParameterPass,
        holdoutParameterPass,
        independentHoldoutPass,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      rows: rowReports,
      formalization: {
        targetId: "contact_parameter_bench_validation_gate",
        leanStructure: "Mechanics.ContactParameterBenchValidationNat",
        leanPredicate: "contactParameterBenchValidationReadyNat",
        schema: "rad-sim.contact-parameter-bench-validation.v1",
      },
      claimLabels: {
        comparison: "bench comparison artifact",
        parameterFit: "empirical comparison only; not a first-principles contact law",
        holdout: "independent holdout check if collected after profile freeze",
        formalCompleteness: "Lean-proven finite pass/fail predicate over result counts and flags",
      },
      limitations: [
        "Passing this gate requires filled rows but does not prove friction, damping, or contact constitutive laws.",
        "Fit and holdout rows must come from independent collection after the profile is frozen.",
        "Radius and clearance fields are normalized proxy comparisons until the hardware unit scale is calibrated.",
      ],
    };
  }

  function exportContactParameterBenchValidation(packetOrState, results = null, options = {}) {
    return JSON.stringify(compareContactParameterCalibrationResults(packetOrState, results, options), null, 2);
  }

  function exportContactParameterBenchValidationCsv(report) {
    const metrics = report.metrics || {};
    const summary = report.summary || {};
    const rows = [[
      "dataset_id",
      "dataset_role",
      "scenario_id",
      "repeat_index",
      "row",
      "col",
      "pivot",
      "pin",
      "hole",
      "missing_measurement_count",
      "max_parameter_residual",
      "pass",
    ]];
    for (const row of report.rows || []) {
      rows.push([
        row.datasetId || "",
        row.datasetRole || "",
        row.scenarioId || "",
        row.repeatIndex ?? "",
        row.row ?? "",
        row.col ?? "",
        row.pivot || "",
        row.pin || "",
        row.hole || "",
        row.missingMeasurementCount ?? "",
        row.maxParameterResidual ?? "",
        row.pass ?? "",
      ]);
    }
    rows.push([
      "summary",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      metrics.missingMeasurementCount ?? "",
      metrics.maxFitHoldoutParameterDelta ?? "",
      summary.contactParameterBenchValidationPass ?? "",
    ]);
    return rows.map((row) => row.map(csvValue).join(",")).join("\n");
  }

  const CONTACT_PARAMETER_INTERVAL_FIELDS = [
    ["pinRadius", "pinRadiusMm", "normalized radius"],
    ["holeRadius", "holeRadiusMm", "normalized radius"],
    ["pinHoleClearance", "clearanceMm", "normalized radius"],
    ["staticFrictionCoeff", "measuredStaticFrictionCoeff", "coefficient"],
    ["dynamicFrictionCoeff", "measuredDynamicFrictionCoeff", "coefficient"],
    ["contactStiffness", "fittedContactStiffness", "proxy stiffness"],
    ["contactDamping", "fittedContactDamping", "proxy damping"],
    ["solrefTimeconst", "fittedSolrefTimeconst", "MuJoCo solref"],
    ["solrefDampingRatio", "fittedSolrefDampingRatio", "MuJoCo solref"],
    ["solimpWidth", "fittedSolimpWidth", "MuJoCo solimp"],
  ];

  function contactParameterIntervalCalibrationReport(packetOrState, results = null, options = {}) {
    const packet =
      packetOrState?.schema === "rad-sim.contact-parameter-calibration-packet.v1"
        ? packetOrState
        : contactParameterCalibrationPacket(packetOrState, options);
    if (packet.schema !== "rad-sim.contact-parameter-calibration-packet.v1") {
      throw new Error("unsupported contact-parameter calibration packet schema");
    }
    const filledResults = results || contactParameterCalibrationResultsTemplate(packet, options);
    if (filledResults.schema !== "rad-sim.contact-parameter-calibration-results.v1") {
      throw new Error("unsupported contact-parameter calibration results schema");
    }
    const validation =
      options.benchValidation ||
      compareContactParameterCalibrationResults(packet, filledResults, {
        tolerance: options.tolerance ?? 1e-9,
        holdoutTolerance: options.holdoutTolerance ?? 1e-6,
      });
    const validationPass = Boolean(validation.summary?.contactParameterBenchValidationPass);
    const resultRows = Array.isArray(filledResults.measurements)
      ? filledResults.measurements.filter((row) => row && typeof row === "object")
      : [];
    const margin = Math.max(0, Number(options.confidenceMargin ?? 0));
    const parameters = [];
    let acceptedCount = 0;
    let uncertaintyCount = 0;
    let holdoutAgreementCount = 0;
    for (const [parameterName, measurementField, unit] of CONTACT_PARAMETER_INTERVAL_FIELDS) {
      const observedValues = [];
      const fitValues = [];
      const holdoutValues = [];
      const residuals = [];
      const assumedValues = [];
      for (const row of resultRows) {
        const observed = contactParameterResultNumber(row, measurementField);
        const expected = contactParameterExpectedValue(packet, row, measurementField);
        if (expected !== null) assumedValues.push(expected);
        if (observed === null) continue;
        observedValues.push(observed);
        if (row.datasetRole === "holdout") holdoutValues.push(observed);
        else fitValues.push(observed);
        if (expected !== null) residuals.push(Math.abs(observed - expected));
      }
      const assumedValue = meanNumber(assumedValues);
      const meanObserved = meanNumber(observedValues);
      const fitMean = meanNumber(fitValues);
      const holdoutMean = meanNumber(holdoutValues);
      const sampleMin = observedValues.length ? Math.min(...observedValues) : null;
      const sampleMax = observedValues.length ? Math.max(...observedValues) : null;
      const maxResidual = residuals.length ? Math.max(...residuals) : null;
      const fitHoldoutDelta = fitMean !== null && holdoutMean !== null ? Math.abs(fitMean - holdoutMean) : null;
      let halfWidth = null;
      let lowerBound = null;
      let upperBound = null;
      if (meanObserved !== null && sampleMin !== null && sampleMax !== null) {
        const spreadHalfWidth = Math.max(Math.abs(meanObserved - sampleMin), Math.abs(sampleMax - meanObserved));
        halfWidth = Math.max(spreadHalfWidth, maxResidual ?? 0, fitHoldoutDelta ?? 0, margin);
        lowerBound = meanObserved - halfWidth;
        upperBound = meanObserved + halfWidth;
      }
      const accepted = assumedValue !== null && lowerBound !== null && upperBound !== null && lowerBound <= assumedValue && assumedValue <= upperBound;
      if (accepted) acceptedCount += 1;
      if (halfWidth !== null) uncertaintyCount += 1;
      if (fitHoldoutDelta !== null) holdoutAgreementCount += 1;
      parameters.push({
        parameter: parameterName,
        measurementField,
        unit,
        assumedValue,
        meanObserved,
        fitMean,
        holdoutMean,
        sampleMin,
        sampleMax,
        halfWidth,
        lowerBound,
        upperBound,
        maxResidualToAssumption: maxResidual,
        fitHoldoutDelta,
        sampleCount: observedValues.length,
        fitSampleCount: fitValues.length,
        holdoutSampleCount: holdoutValues.length,
        acceptedAssumptionInsideBounds: accepted,
        claimLabel: "empirical interval estimate from filled fit/holdout rows",
      });
    }
    const missingEvidence = [];
    if (!validationPass) missingEvidence.push("contactParameterBenchValidationPass");
    if (!resultRows.length) missingEvidence.push("resultRows");
    if (parameters.length < CONTACT_PARAMETER_INTERVAL_FIELDS.length) missingEvidence.push("parameterIntervals");
    if (acceptedCount < parameters.length) missingEvidence.push("simulatorParametersInsideBounds");
    if (uncertaintyCount < parameters.length) missingEvidence.push("uncertaintyBounds");
    if (holdoutAgreementCount < parameters.length) missingEvidence.push("holdoutAgreement");
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.contact-parameter-interval-calibration.v1",
      method: "conservative interval calibration over filled fit and holdout pin-hole/contact-parameter measurements",
      sourcePacket: {
        schema: packet.schema || "",
        profilePlan: cloneJson(packet.profilePlan || {}),
      },
      sourceResults: {
        schema: filledResults.schema || "",
        datasetId: filledResults.datasetId || "",
        datasetRole: filledResults.datasetRole || "",
        sourceFileId: filledResults.sourceFileId || "",
      },
      sourceBenchValidation: {
        schema: validation.schema || "",
        passed: validationPass,
        targetId: validation.formalization?.targetId || "",
      },
      parameters,
      summary: {
        status: ready ? "contact-parameter-interval-calibration-ready" : "needs-contact-parameter-interval-review",
        contactParameterIntervalCalibrationReady: ready,
        benchValidationPass: validationPass,
        parameterIntervalCount: parameters.length,
        acceptedParameterIntervalCount: acceptedCount,
        uncertaintyRecordCount: uncertaintyCount,
        holdoutAgreementRecordCount: holdoutAgreementCount,
        simulatorParametersInsideBounds: acceptedCount === parameters.length,
        confidenceMargin: margin,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "contact_parameter_interval_calibration_gate",
        leanStructure: "Mechanics.ContactParameterIntervalCalibrationNat",
        leanPredicate: "contactParameterIntervalCalibrationReadyNat",
        schema: "rad-sim.contact-parameter-interval-calibration.v1",
      },
      claimLabels: {
        intervals: "empirical bounded-parameter artifact",
        acceptance: "frozen simulator parameter lies inside measured interval",
        formalCompleteness: "Lean-proven finite readiness predicate over interval evidence",
        physics: "not a first-principles contact mechanics proof",
      },
      limitations: [
        "Intervals are only as valid as the bench/result rows and measurement provenance.",
        "A parameter inside an empirical interval can still fail outside the tested load, velocity, wear, or assembly regime.",
        "This artifact bounds proxy MuJoCo/simulator parameters; it does not derive friction or damping from continuum mechanics.",
      ],
    };
  }

  function exportContactParameterIntervalCalibration(packetOrState, results = null, options = {}) {
    return JSON.stringify(contactParameterIntervalCalibrationReport(packetOrState, results, options), null, 2);
  }

  function exportContactParameterIntervalCalibrationCsv(report) {
    const summary = report.summary || {};
    const rows = [[
      "parameter",
      "measurement_field",
      "unit",
      "assumed_value",
      "mean_observed",
      "fit_mean",
      "holdout_mean",
      "lower_bound",
      "upper_bound",
      "half_width",
      "max_residual_to_assumption",
      "fit_holdout_delta",
      "sample_count",
      "accepted_assumption_inside_bounds",
    ]];
    for (const parameter of report.parameters || []) {
      rows.push([
        parameter.parameter || "",
        parameter.measurementField || "",
        parameter.unit || "",
        parameter.assumedValue ?? "",
        parameter.meanObserved ?? "",
        parameter.fitMean ?? "",
        parameter.holdoutMean ?? "",
        parameter.lowerBound ?? "",
        parameter.upperBound ?? "",
        parameter.halfWidth ?? "",
        parameter.maxResidualToAssumption ?? "",
        parameter.fitHoldoutDelta ?? "",
        parameter.sampleCount ?? "",
        parameter.acceptedAssumptionInsideBounds ?? "",
      ]);
    }
    rows.push([
      "summary",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      summary.confidenceMargin ?? "",
      "",
      "",
      summary.parameterIntervalCount ?? "",
      summary.contactParameterIntervalCalibrationReady ?? "",
    ]);
    return rows.map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function mujocoExternalRunReport(state, options = {}) {
    const exportReport = options.exportReport || state.experiment?.mujocoModelExport || mujocoModelExportReport(state, options);
    const exportSummary = exportReport?.summary || {};
    const missingEvidence = [];
    if (!exportSummary.mujocoModelExportReady) missingEvidence.push("mujocoModelExport");
    missingEvidence.push("browserCannotExecuteMujoco");
    return {
      schema: "rad-sim.mujoco-external-run.v1",
      method: "browser placeholder for an external MuJoCo run; use Python to execute MuJoCo",
      engine: "MuJoCo",
      modelExport: exportReport,
      solver: {
        engineAvailable: false,
        ran: false,
        stepsRequested: Math.max(1, Math.round(Number(options.steps ?? 120))),
        stepsCompleted: 0,
        error: "MuJoCo cannot be executed from this static browser UI.",
      },
      results: { bodies: [] },
      summary: {
        status: "needs-mujoco-external-run-evidence",
        mujocoExternalRunComplete: false,
        engineAvailable: false,
        bodyResultCount: 0,
        expectedBodyResultCount: Math.max(0, Math.round(Number(exportSummary.bodyRecordCount || 0))),
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "mujoco_external_run_gate",
        leanStructure: "Mechanics.ExternalPhysicsRunNat",
        leanPredicate: "externalPhysicsRunReadyNat",
        schema: "rad-sim.mujoco-external-run.v1",
      },
      claimLabels: {
        run: "external-engine run not completed in browser",
        physicalAccuracy: "experimentally unvalidated",
      },
    };
  }

  function exportMujocoExternalRun(state, options = {}) {
    return JSON.stringify(mujocoExternalRunReport(state, options), null, 2);
  }

  function mujocoExternalComparisonReport(state, runReport, options = {}) {
    const eventSequence = Array.isArray(options.eventSequence) ? options.eventSequence : [];
    const comparisonState = mujocoExportState(state, eventSequence);
    const sim = RAD.simulate(comparisonState);
    const runSummary = runReport?.summary || {};
    const resultMap = new Map();
    for (const result of runReport?.results?.bodies || []) {
      resultMap.set(contactGraphKey(Number(result.row), Number(result.col)), result);
    }
    const tolerance = Math.max(0, Number(options.tolerance ?? 1e-6));
    const records = [];
    let missingBodyRecords = 0;
    let maxPositionError = 0;
    let squaredError = 0;
    let matchedCount = 0;
    for (let row = 0; row < comparisonState.grid.rows; row += 1) {
      for (let col = 0; col < comparisonState.grid.cols; col += 1) {
        if (comparisonState.cells.removed?.[row]?.[col]) continue;
        const expected = simulationCenterVector(sim, comparisonState, row, col);
        const result = resultMap.get(contactGraphKey(row, col));
        if (!result) {
          missingBodyRecords += 1;
          records.push({ row, col, matched: false, simulatorPosition: expected, externalPosition: null, positionError: null });
          continue;
        }
        const external = vec3(result.finalPosition || [0, 0, 0]);
        const dx = external[0] - expected[0];
        const dy = external[1] - expected[1];
        const dz = external[2] - expected[2];
        const error = Math.sqrt(dx * dx + dy * dy + dz * dz);
        maxPositionError = Math.max(maxPositionError, error);
        squaredError += error * error;
        matchedCount += 1;
        records.push({ row, col, matched: true, simulatorPosition: expected, externalPosition: external, positionError: error, withinTolerance: error <= tolerance });
      }
    }
    const missingEvidence = [];
    if (!runSummary.mujocoExternalRunComplete) missingEvidence.push("mujocoExternalRun");
    if (!records.length) missingEvidence.push("comparisonRecords");
    if (missingBodyRecords) missingEvidence.push("matchedBodyRecords");
    if (maxPositionError > tolerance) missingEvidence.push("positionTolerance");
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.mujoco-external-comparison.v1",
      method: "finite simulator-vs-MuJoCo body-center comparison for external validation review",
      engine: "MuJoCo",
      summary: {
        status: ready ? "mujoco-external-comparison-ready" : "needs-mujoco-external-comparison-evidence",
        mujocoExternalComparisonReady: ready,
        comparisonRecordCount: records.length,
        matchedBodyRecordCount: matchedCount,
        missingBodyRecordCount: missingBodyRecords,
        rmsPositionError: matchedCount ? Math.sqrt(squaredError / matchedCount) : 0,
        maxPositionError,
        tolerance,
        withinTolerance: maxPositionError <= tolerance,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      records,
      runSummary,
      formalization: {
        targetId: "mujoco_external_comparison_gate",
        leanStructure: "Mechanics.ExternalPhysicsComparisonNat",
        leanPredicate: "externalPhysicsComparisonReadyNat",
        schema: "rad-sim.mujoco-external-comparison.v1",
      },
      claimLabels: {
        comparison: "external-engine-derived comparison if the run report is complete",
        simulatorAgreement: "numerical agreement diagnostic, not physical bench validation",
        physicalAccuracy: "experimentally unvalidated until bench measurements agree",
      },
      limitations: [
        "This compares body-center positions only.",
        "A passing external comparison does not validate contact, friction, or actuator mechanics without bench data.",
        "The comparison inherits the coarse MJCF export assumptions.",
      ],
    };
  }

  function exportMujocoExternalComparison(state, runReport, options = {}) {
    return JSON.stringify(mujocoExternalComparisonReport(state, runReport, options), null, 2);
  }

  function exportMujocoExternalComparisonCsv(report) {
    const rows = [["row", "col", "matched", "simulator_position", "external_position", "position_error", "within_tolerance"]];
    for (const record of report.records || []) {
      rows.push([
        record.row,
        record.col,
        record.matched,
        Array.isArray(record.simulatorPosition) ? record.simulatorPosition.join(";") : "",
        Array.isArray(record.externalPosition) ? record.externalPosition.join(";") : "",
        record.positionError ?? "",
        record.withinTolerance ?? "",
      ]);
    }
    rows.push(["summary", "", report.summary?.mujocoExternalComparisonReady ?? "", "", "", report.summary?.maxPositionError ?? "", report.summary?.withinTolerance ?? ""]);
    return rows.map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function finalStateForEquilibrium(state, sequence) {
    if (!Array.isArray(sequence) || !sequence.length || typeof RAD.applyEventSequence !== "function") {
      return cloneForCharacterization(state);
    }
    return RAD.applyEventSequence(cloneForCharacterization(state), sequence);
  }

  function browserSpringEnergyProxy(physical) {
    let springEnergy = 0;
    let activeSpringEdges = 0;
    for (const row of physical?.linkStrain?.horizontal || []) {
      for (const value of row || []) {
        if (value === null || value === undefined) continue;
        activeSpringEdges += 1;
        springEnergy += 0.5 * Number(value) * Number(value);
      }
    }
    for (const row of physical?.linkStrain?.vertical || []) {
      for (const value of row || []) {
        if (value === null || value === undefined) continue;
        activeSpringEdges += 1;
        springEnergy += 0.5 * Number(value) * Number(value);
      }
    }
    return { springEnergy, activeSpringEdges };
  }

  function browserLockEnergyProxy(state, physical) {
    let energy = 0;
    for (let row = 0; row < state.grid.rows; row += 1) {
      for (let col = 0; col < state.grid.cols; col += 1) {
        if (!state.cells.locked?.[row]?.[col]) continue;
        const z = Number(physical?.height?.[row]?.[col]) || 0;
        energy += 0.5 * z * z;
      }
    }
    return energy;
  }

  function equilibriumRelationReport(state, options = {}) {
    const tolerance = Math.max(0, Number(options.tolerance ?? 1e-6));
    const residualTolerance = Math.max(0, Number(options.residualTolerance ?? tolerance));
    const contactStiffness = Math.max(0, Number(options.contactStiffness ?? 1));
    const sequence = Array.isArray(options.eventSequence) ? options.eventSequence : null;
    const finalState = finalStateForEquilibrium(state, sequence);
    const realization =
      options.realizationReport ||
      physicalRealizationMapReport(state, {
        eventSequence: sequence,
        contactReport: options.contactReport,
        contactGraphReport: options.contactGraphReport,
        contactStiffness,
        tolerance,
      });
    const contact =
      options.contactReport ||
      contactStateAbstractionReport(finalState, { contactStiffness, tolerance });
    const physical =
      typeof RAD.simulatePhysicalRelaxation === "function"
        ? RAD.simulatePhysicalRelaxation(finalState, { iterations: options.iterations ?? 18 })
        : RAD.simulate(finalState);
    const spring = browserSpringEnergyProxy(physical);
    const hingeEnergy = Math.max(0, Number(physical?.slope?.mean || 0) ** 2);
    const lockPenaltyEnergy = browserLockEnergyProxy(finalState, physical);
    const contactPenaltyEnergy = Math.max(0, Number(contact?.summary?.totalContactPenalty) || 0);
    const loadWorkMagnitudeProxy = Math.abs(Number(options.externalLoadWorkProxy) || 0);
    const storedEnergy = spring.springEnergy + hingeEnergy + lockPenaltyEnergy + contactPenaltyEnergy;
    const objectiveEnergy = storedEnergy + loadWorkMagnitudeProxy;
    const objectiveBalanceError = Math.abs(objectiveEnergy - (storedEnergy + loadWorkMagnitudeProxy));
    const residual = Math.max(
      0,
      Number(physical?.metrics?.physicalRmsHeightDelta ?? physical?.metrics?.modelErrorHeight ?? 0) || 0
    );
    const nonnegativeTerms = [
      spring.springEnergy,
      hingeEnergy,
      lockPenaltyEnergy,
      contactPenaltyEnergy,
      loadWorkMagnitudeProxy,
    ];
    const nonnegativeTermCount = nonnegativeTerms.filter((value) => value >= -tolerance).length;
    const missingEvidence = [];
    if (!realization?.summary?.physicalRealizationMapReady) missingEvidence.push("physicalRealizationMap");
    if (!physical) missingEvidence.push("solverSuccess");
    if (residual > residualTolerance) missingEvidence.push("equilibriumResidual");
    if (nonnegativeTermCount < nonnegativeTerms.length) missingEvidence.push("nonnegativeEnergyTerms");
    if (objectiveBalanceError > tolerance) missingEvidence.push("energyBalance");
    if (!contact?.summary?.contactStateAbstractionReady) missingEvidence.push("contactStateAbstraction");
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.equilibrium-relation.v1",
      method: "finite equilibrium evidence gate tying a realized operator state to browser spring-preview status, nonnegative stored/proxy energy terms, contact-state evidence, and residual tolerance",
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        backlash: state.grid.backlash,
        pinHoleClearance: typeof RAD.pinHoleClearance === "function" ? RAD.pinHoleClearance(state) : 0,
      },
      solver: {
        model: physical?.metrics?.model || "browser-spring-preview",
        success: Boolean(physical),
        message: "browser spring-preview proxy",
        iterations: Number(physical?.metrics?.physicalIterations) || 0,
        springEdges: spring.activeSpringEdges,
        hingeTriples: 0,
        removedCells: (finalState.cells.removed || []).flat().filter(Boolean).length,
      },
      residual: {
        targetRmsError: residual,
        tolerance: residualTolerance,
        passesTolerance: residual <= residualTolerance,
        objectiveBalanceError,
      },
      energy: {
        objectiveEnergy,
        storedEnergy,
        springEnergy: spring.springEnergy,
        hingeEnergy,
        lockPenaltyEnergy,
        externalPotentialEnergy: -loadWorkMagnitudeProxy,
        contactPenaltyEnergy,
        loadWorkMagnitudeProxy,
        nonnegativeTermCount,
        requiredNonnegativeTermCount: nonnegativeTerms.length,
      },
      realization: {
        attached: Boolean(realization),
        schema: realization?.schema || "",
        physicalRealizationMapReady: Boolean(realization?.summary?.physicalRealizationMapReady),
        abstractOperatorCount: Number(realization?.summary?.abstractOperatorCount) || 0,
      },
      contact: {
        attached: Boolean(contact),
        schema: contact?.schema || "",
        contactStateAbstractionReady: Boolean(contact?.summary?.contactStateAbstractionReady),
        contactPenaltyEnergy,
      },
      summary: {
        status: ready ? "equilibrium-relation-ready" : "needs-equilibrium-evidence",
        equilibriumRelationReady: ready,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "equilibrium_relation_gate",
        leanStructure: "Mechanics.EquilibriumRelationNat",
        leanPredicate: "equilibriumRelationReadyNat",
        schema: "rad-sim.equilibrium-relation.v1",
      },
      claimLabels: {
        equilibriumRelation: "Lean-proven finite evidence gate",
        energyTerms: "browser spring-preview proxy",
        solver: "interactive numerical approximation",
        physicalAccuracy: "experimentally unvalidated physical assumption",
      },
      limitations: [
        "The browser equilibrium relation is a finite evidence gate over a spring-preview proxy.",
        "Use Python solve_spring_hinge_3d reports for stronger research artifacts.",
        "Residual tolerance is a user-chosen numerical threshold, not a hardware-calibrated bound.",
      ],
    };
  }

  function exportEquilibriumRelation(state, options = {}) {
    return JSON.stringify(equilibriumRelationReport(state, options), null, 2);
  }

  function exportEquilibriumRelationCsv(reportOrState, options = {}) {
    const report =
      reportOrState?.schema === "rad-sim.equilibrium-relation.v1"
        ? reportOrState
        : reportOrState?.experiment?.equilibriumRelation || equilibriumRelationReport(reportOrState, options);
    const solver = report.solver || {};
    const residual = report.residual || {};
    const energy = report.energy || {};
    const summary = report.summary || {};
    return [
      [
        "schema",
        "status",
        "equilibrium_relation_ready",
        "solver_success",
        "target_rms_error",
        "residual_tolerance",
        "objective_balance_error",
        "stored_energy",
        "spring_energy",
        "hinge_energy",
        "lock_penalty_energy",
        "contact_penalty_energy",
        "load_work_magnitude_proxy",
        "missing_evidence",
      ],
      [
        report.schema,
        summary.status || "",
        summary.equilibriumRelationReady,
        solver.success,
        residual.targetRmsError ?? "",
        residual.tolerance ?? "",
        residual.objectiveBalanceError ?? "",
        energy.storedEnergy ?? "",
        energy.springEnergy ?? "",
        energy.hingeEnergy ?? "",
        energy.lockPenaltyEnergy ?? "",
        energy.contactPenaltyEnergy ?? "",
        energy.loadWorkMagnitudeProxy ?? "",
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
      ],
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function reachabilityEventActuatorCells(event) {
    if (event?.kind === "actuate") {
      const cell = normalizeContactGraphCell(event.cell);
      return cell ? [cell] : [];
    }
    if (event?.kind === "group_actuate") {
      return (event.cells || []).map(normalizeContactGraphCell).filter(Boolean);
    }
    return [];
  }

  function reachabilityUniqueCells(cells) {
    const out = [];
    const seen = new Set();
    for (const cell of cells || []) {
      const normalized = normalizeContactGraphCell(cell);
      if (!normalized) continue;
      const key = contactGraphKey(normalized.row, normalized.col);
      if (seen.has(key)) continue;
      seen.add(key);
      out.push(normalized);
    }
    return out;
  }

  function reachabilityActuatorBasis(state, finalState, sequence, explicit, tolerance) {
    const cells = [];
    if (Array.isArray(explicit)) {
      cells.push(...explicit);
    } else {
      cells.push(...realizationCommandSupport(state, tolerance));
      cells.push(...realizationCommandSupport(finalState, tolerance));
      for (const event of sequence || []) cells.push(...reachabilityEventActuatorCells(event));
    }
    return reachabilityUniqueCells(cells).filter(
      (cell) =>
        cell.row >= 0 &&
        cell.row < finalState.grid.rows &&
        cell.col >= 0 &&
        cell.col < finalState.grid.cols &&
        !finalState.cells.removed?.[cell.row]?.[cell.col]
    );
  }

  function reachabilityTargetCells(state, explicit) {
    if (Array.isArray(explicit)) return reachabilityUniqueCells(explicit);
    const cells = [];
    for (let row = 0; row < state.grid.rows; row += 1) {
      for (let col = 0; col < state.grid.cols; col += 1) {
        if (!state.cells.removed?.[row]?.[col]) cells.push({ row, col });
      }
    }
    return cells;
  }

  function reachabilityMapFromMatrix(matrix, state, tolerance) {
    const out = RAD.matrix(state.grid.rows, state.grid.cols, false);
    for (let row = 0; row < state.grid.rows; row += 1) {
      for (let col = 0; col < state.grid.cols; col += 1) {
        const index = row * state.grid.cols + col;
        out[row][col] = (matrix?.[index] || []).some((value) => Math.abs(Number(value) || 0) > tolerance);
      }
    }
    return out;
  }

  function topologyReachabilityMap(topology, basis, state) {
    const activeComponents = new Set();
    for (const cell of basis) {
      const label = topology.componentLabels?.[cell.row]?.[cell.col];
      if (Number.isInteger(label) && label >= 0) activeComponents.add(label);
    }
    return RAD.matrix(state.grid.rows, state.grid.cols, (row, col) => {
      const label = topology.componentLabels?.[row]?.[col];
      return Number.isInteger(label) && label >= 0 && activeComponents.has(label);
    });
  }

  function countTargetMap(targets, map) {
    return targets.filter((cell) => Boolean(map?.[cell.row]?.[cell.col])).length;
  }

  function reachableEquilibriumControllabilityReport(state, options = {}) {
    const tolerance = Math.max(0, Number(options.tolerance ?? 1e-9));
    const sequence = Array.isArray(options.eventSequence) ? options.eventSequence : [];
    const finalState = finalStateForEquilibrium(state, sequence);
    const basis = reachabilityActuatorBasis(state, finalState, sequence, options.actuatorCells, tolerance);
    const targets = reachabilityTargetCells(finalState, options.targetCells);
    const matrix = buildResponseMatrix(finalState, {
      actuatorCells: basis.map((cell) => ({ r: cell.row, c: cell.col })),
      alphaStep: options.alphaStep ?? 0.12,
      zStep: options.zStep ?? 0.12,
      includeAlpha: options.includeAlpha,
      includeZ: options.includeZ,
      tolerance,
    });
    const topology = typeof RAD.topologyDiagnostics === "function" ? RAD.topologyDiagnostics(finalState) : null;
    const alphaMap = reachabilityMapFromMatrix(matrix.alpha, finalState, tolerance);
    const heightMap = reachabilityMapFromMatrix(matrix.height, finalState, tolerance);
    const topologyMap = topology ? topologyReachabilityMap(topology, basis, finalState) : RAD.matrix(finalState.grid.rows, finalState.grid.cols, false);
    const targetAlphaReachableCells = countTargetMap(targets, alphaMap);
    const targetHeightReachableCells = countTargetMap(targets, heightMap);
    const topologyReachableTargetCells = countTargetMap(targets, topologyMap);
    const targetAlphaUnderactuatedCells = Math.max(0, targets.length - targetAlphaReachableCells);
    const targetHeightUnderactuatedCells = Math.max(0, targets.length - targetHeightReachableCells);
    const topologyBlockedTargetCells = Math.max(0, targets.length - topologyReachableTargetCells);
    const equilibrium =
      options.equilibriumReport ||
      equilibriumRelationReport(state, {
        eventSequence: sequence,
        tolerance,
        residualTolerance: options.residualTolerance,
        contactStiffness: options.contactStiffness,
      });
    const requireFull = Boolean(options.requireFullTargetReachability);
    const missingEvidence = [];
    if (!equilibrium?.summary?.equilibriumRelationReady) missingEvidence.push("equilibriumRelation");
    if (!basis.length) missingEvidence.push("actuatorBasis");
    if (!matrix.commands?.length) missingEvidence.push("responseColumns");
    if (!targets.length) missingEvidence.push("targetCells");
    if ((matrix.diagnostics?.reachableAlphaCells || 0) + (matrix.diagnostics?.reachableHeightCells || 0) <= 0) {
      missingEvidence.push("reachableResponse");
    }
    if (requireFull && (targetAlphaUnderactuatedCells > 0 || targetHeightUnderactuatedCells > 0 || topologyBlockedTargetCells > 0)) {
      missingEvidence.push("targetReachability");
    }
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.reachable-equilibrium-controllability.v1",
      method: "finite actuator-basis reachability report over an equilibrium state, response matrix, and removed-cell topology graph",
      grid: { rows: state.grid.rows, cols: state.grid.cols, totalCells: state.grid.rows * state.grid.cols },
      actuatorBasis: {
        cells: basis,
        cellCount: basis.length,
        commandColumnCount: matrix.commands?.length || 0,
        alphaStep: Number(options.alphaStep ?? 0.12),
        zStep: Number(options.zStep ?? 0.12),
      },
      targets: {
        cells: targets,
        targetCellCount: targets.length,
        requireFullTargetReachability: requireFull,
      },
      response: {
        alphaRank: matrix.diagnostics?.alphaRank || 0,
        heightRank: matrix.diagnostics?.heightRank || 0,
        reachableAlphaCells: matrix.diagnostics?.reachableAlphaCells || 0,
        reachableHeightCells: matrix.diagnostics?.reachableHeightCells || 0,
        alphaUnderactuatedCells: matrix.diagnostics?.alphaUnderactuatedCells || 0,
        heightUnderactuatedCells: matrix.diagnostics?.heightUnderactuatedCells || 0,
        targetAlphaReachableCells,
        targetHeightReachableCells,
        targetAlphaUnderactuatedCells,
        targetHeightUnderactuatedCells,
        tolerance,
        alphaReachableMap: alphaMap,
        heightReachableMap: heightMap,
      },
      topology: {
        componentCount: topology?.componentCount || 0,
        deletedEdgeCount: topology?.deletedEdges || 0,
        removedCellCount: topology?.removedCells || 0,
        topologyReachableTargetCells,
        topologyBlockedTargetCells,
        componentLabels: topology?.componentLabels || [],
        topologyReachableMap: topologyMap,
      },
      equilibrium: {
        attached: Boolean(equilibrium),
        schema: equilibrium?.schema || "",
        equilibriumRelationReady: Boolean(equilibrium?.summary?.equilibriumRelationReady),
        status: equilibrium?.summary?.status || "",
      },
      summary: {
        status: ready ? "reachable-equilibrium-controllability-ready" : "needs-reachability-controllability-evidence",
        reachableEquilibriumControllabilityReady: ready,
        targetFullyReachable:
          targets.length > 0 &&
          targetAlphaUnderactuatedCells === 0 &&
          targetHeightUnderactuatedCells === 0 &&
          topologyBlockedTargetCells === 0,
        targetUnderactuatedCells: Math.max(targetAlphaUnderactuatedCells, targetHeightUnderactuatedCells, topologyBlockedTargetCells),
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "reachable_equilibrium_controllability_gate",
        leanStructure: "Mechanics.ReachableEquilibriumControllabilityNat",
        leanPredicate: "reachableEquilibriumControllabilityReadyNat",
        schema: "rad-sim.reachable-equilibrium-controllability.v1",
      },
      claimLabels: {
        reachability: "finite response-matrix simulator diagnostic",
        equilibrium: "Lean-proven finite equilibrium evidence gate",
        topology: "graph deletion and component bookkeeping",
        controllability: "linearized finite actuator-basis proxy",
        physicalAccuracy: "experimentally unvalidated physical assumption",
      },
      limitations: [
        "Reachability is measured through finite response columns, not full nonlinear controllability.",
        "Topology reachability is graph-component bookkeeping, not a proof of mechanical motion feasibility.",
        "Full target reachability is optional because underactuation is itself useful design information.",
      ],
    };
  }

  function exportReachableEquilibriumControllability(state, options = {}) {
    return JSON.stringify(reachableEquilibriumControllabilityReport(state, options), null, 2);
  }

  function exportReachableEquilibriumControllabilityCsv(reportOrState, options = {}) {
    const report =
      reportOrState?.schema === "rad-sim.reachable-equilibrium-controllability.v1"
        ? reportOrState
        : reportOrState?.experiment?.reachableEquilibriumControllability || reachableEquilibriumControllabilityReport(reportOrState, options);
    const response = report.response || {};
    const topology = report.topology || {};
    const basis = report.actuatorBasis || {};
    const targets = report.targets || {};
    const summary = report.summary || {};
    return [
      [
        "schema",
        "status",
        "reachable_equilibrium_controllability_ready",
        "target_fully_reachable",
        "actuator_basis_cells",
        "command_columns",
        "target_cells",
        "alpha_rank",
        "height_rank",
        "reachable_alpha_cells",
        "reachable_height_cells",
        "target_alpha_underactuated_cells",
        "target_height_underactuated_cells",
        "topology_blocked_target_cells",
        "component_count",
        "deleted_edges",
        "missing_evidence",
      ],
      [
        report.schema,
        summary.status || "",
        summary.reachableEquilibriumControllabilityReady,
        summary.targetFullyReachable,
        basis.cellCount ?? "",
        basis.commandColumnCount ?? "",
        targets.targetCellCount ?? "",
        response.alphaRank ?? "",
        response.heightRank ?? "",
        response.reachableAlphaCells ?? "",
        response.reachableHeightCells ?? "",
        response.targetAlphaUnderactuatedCells ?? "",
        response.targetHeightUnderactuatedCells ?? "",
        topology.topologyBlockedTargetCells ?? "",
        topology.componentCount ?? "",
        topology.deletedEdgeCount ?? "",
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
      ],
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function reachableEquilibriumCellList(raw) {
    return reachabilityUniqueCells(Array.isArray(raw) ? raw : []);
  }

  function reachableBenchCommand(cell, alpha = 0, z = 0) {
    return { row: cell.row, col: cell.col, alpha, z };
  }

  function reachableBenchBoolMapGet(map, cell) {
    return Boolean(map?.[cell.row]?.[cell.col]);
  }

  function reachableBenchIntMapGet(map, cell, fallback = -1) {
    const value = map?.[cell.row]?.[cell.col];
    return Number.isInteger(value) ? value : fallback;
  }

  function reachableEquilibriumTargetSummaries(report) {
    const targets = reachableEquilibriumCellList(report?.targets?.cells || []);
    const alphaMap = report?.response?.alphaReachableMap || [];
    const heightMap = report?.response?.heightReachableMap || [];
    const topologyMap = report?.topology?.topologyReachableMap || [];
    const componentLabels = report?.topology?.componentLabels || [];
    return targets.map((cell) => {
      const topologyReachable = reachableBenchBoolMapGet(topologyMap, cell);
      return {
        row: cell.row,
        col: cell.col,
        alphaReachable: reachableBenchBoolMapGet(alphaMap, cell),
        heightReachable: reachableBenchBoolMapGet(heightMap, cell),
        topologyReachable,
        componentLabel: reachableBenchIntMapGet(componentLabels, cell),
        blockedByTopology: !topologyReachable,
      };
    });
  }

  function reachableEquilibriumBenchStep({
    id,
    scope,
    commands = [],
    observationCells = [],
    targetCells = [],
    purpose,
    expectedResponse,
    passFailCriterion,
    repeatCount,
    claimLabel,
  }) {
    return {
      id,
      scope,
      commands,
      observationCells: reachabilityUniqueCells(observationCells),
      targetCells: reachabilityUniqueCells(targetCells),
      requiredMeasurements: {
        fields: [...REACHABLE_EQUILIBRIUM_BENCH_FIELDS],
        coordinateFrame: "registered cell-center grid with pixel/mm or probe/mm scale",
      },
      purpose,
      expectedResponse,
      passFailCriterion,
      repeatCount,
      claimLabel,
    };
  }

  function reachableEquilibriumBenchProtocol(state, options = {}) {
    const alphaStep = Number(options.alphaStep ?? 0.12);
    const zStep = Number(options.zStep ?? 0.12);
    const repeatCount = Math.max(1, Math.round(Number(options.repeatCount ?? 3)));
    const report =
      options.controllabilityReport ||
      reachableEquilibriumControllabilityReport(state, {
        eventSequence: options.eventSequence,
        actuatorCells: options.actuatorCells,
        targetCells: options.targetCells,
        alphaStep,
        zStep,
        residualTolerance: options.residualTolerance,
        tolerance: options.tolerance,
      });
    const basis = reachableEquilibriumCellList(report?.actuatorBasis?.cells || []);
    const targets = reachableEquilibriumCellList(report?.targets?.cells || []);
    const targetSummaries = reachableEquilibriumTargetSummaries(report);
    const blockedTargets = targetSummaries
      .filter((target) => target.blockedByTopology)
      .map((target) => ({ row: target.row, col: target.col }));
    const observationCells = reachabilityUniqueCells([...basis, ...targets]);
    const steps = [
      reachableEquilibriumBenchStep({
        id: "baseline_equilibrium_capture",
        scope: "baseline",
        commands: [],
        observationCells,
        targetCells: targets,
        purpose: "Record the undeformed or post-event equilibrium state before any actuator-column perturbation.",
        expectedResponse: "Measured alpha, height, lock, contact, and component labels establish the zero command reference.",
        passFailCriterion:
          "baseline equilibrium residual is at or below the configured measurement tolerance before response columns are collected",
        repeatCount,
        claimLabel: "required boundary-condition and equilibrium evidence",
      }),
    ];
    basis.forEach((cell, index) => {
      steps.push(
        reachableEquilibriumBenchStep({
          id: `basis_${index}_alpha_response_column`,
          scope: "actuator-column",
          commands: [reachableBenchCommand(cell, alphaStep, 0)],
          observationCells,
          targetCells: targets,
          purpose: "Measure the physical alpha response column for one candidate actuator-basis cell.",
          expectedResponse:
            "Target alpha deltas should match the simulator-predicted reachable map within calibrated tolerance after backlash.",
          passFailCriterion: "alpha target deltas, topology labels, and lock states are recorded for every target cell",
          repeatCount,
          claimLabel: "simulator-derived empirical law until bench data exists",
        }),
        reachableEquilibriumBenchStep({
          id: `basis_${index}_z_response_column`,
          scope: "actuator-column",
          commands: [reachableBenchCommand(cell, 0, zStep)],
          observationCells,
          targetCells: targets,
          purpose: "Measure the physical vertical response column and residual neighbor movement for one actuator-basis cell.",
          expectedResponse:
            "Height deltas should show direct actuation plus die-off through pin-hole clearance and graph connectivity.",
          passFailCriterion: "height target deltas and pin-hole slip measurements are recorded for every target cell",
          repeatCount,
          claimLabel: "experimentally unvalidated physical assumption",
        })
      );
    });
    if (basis.length) {
      steps.push(
        reachableEquilibriumBenchStep({
          id: "group_target_reachability_probe",
          scope: "group",
          commands: basis.map((cell) => reachableBenchCommand(cell, alphaStep, zStep)),
          observationCells,
          targetCells: targets,
          purpose: "Probe whether simultaneous basis actuation reaches the target set predicted by finite response-column evidence.",
          expectedResponse:
            "The measured group field should be compared against the linear response-matrix span and a sequenced application.",
          passFailCriterion: "group response is saved with simultaneous/sequenced order metadata and target residuals",
          repeatCount,
          claimLabel: "finite response-matrix controllability proxy",
        })
      );
    }
    if (blockedTargets.length) {
      steps.push(
        reachableEquilibriumBenchStep({
          id: "topology_blocked_target_control",
          scope: "topology-control",
          commands: basis.map((cell) => reachableBenchCommand(cell, 0, zStep)),
          observationCells: blockedTargets,
          targetCells: blockedTargets,
          purpose: "Check whether removed-cell graph deletion really blocks the target component under physical actuation.",
          expectedResponse: "Topology-blocked targets should remain below calibrated noise unless the real sheet has bypass coupling.",
          passFailCriterion:
            "blocked targets are explicitly measured and any nonzero response is labeled as model-mismatch evidence",
          repeatCount,
          claimLabel: "topology/connectivity physical validation trial",
        })
      );
    }
    const topologyPolicySatisfied = blockedTargets.length === 0 || steps.some((step) => step.id === "topology_blocked_target_control");
    const missingEvidence = [];
    if (!report?.summary?.reachableEquilibriumControllabilityReady) missingEvidence.push("reachableEquilibriumControllability");
    if (!basis.length) missingEvidence.push("actuatorBasis");
    if (!targets.length) missingEvidence.push("targetCells");
    if (!steps.length) missingEvidence.push("protocolSteps");
    if (!REACHABLE_EQUILIBRIUM_BENCH_FIELDS.length) missingEvidence.push("measurementFields");
    if (!topologyPolicySatisfied) missingEvidence.push("topologyBlockedControl");
    const passFailCriteria = {
      baseline: "equilibrium residual and fixture registration are recorded before perturbation",
      responseColumns: "each actuator-basis alpha and z response column has target measurements",
      topologyControls: "topology-blocked targets are measured whenever the report predicts a block",
      claimLimit: "passing the protocol validates a bench observation table, not nonlinear controllability",
    };
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.reachable-equilibrium-bench-protocol.v1",
      method: "bench protocol derived from finite reachable-equilibrium controllability evidence",
      grid: { rows: state.grid.rows, cols: state.grid.cols },
      sourceReport: {
        schema: report?.schema || "",
        status: report?.summary?.status || "",
        targetFullyReachable: Boolean(report?.summary?.targetFullyReachable),
        targetUnderactuatedCells: Number(report?.summary?.targetUnderactuatedCells || 0),
      },
      actuatorBasis: basis,
      targetSummaries,
      topology: {
        componentCount: Number(report?.topology?.componentCount || 0),
        topologyBlockedTargetCells: blockedTargets.length,
        topologyPolicySatisfied,
      },
      measurementFields: [...REACHABLE_EQUILIBRIUM_BENCH_FIELDS],
      repeatCount,
      steps,
      passFailCriteria,
      summary: {
        status: ready ? "reachable-equilibrium-bench-protocol-ready" : "needs-reachable-equilibrium-bench-evidence",
        benchProtocolReady: ready,
        stepCount: steps.length,
        actuatorColumnTrialCount: Math.max(0, 2 * basis.length),
        targetObservationCellCount: targets.length,
        topologyBlockedControlCount: blockedTargets.length ? 1 : 0,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "reachable_equilibrium_bench_protocol_gate",
        leanStructure: "Mechanics.ReachableEquilibriumBenchProtocolNat",
        leanPredicate: "reachableEquilibriumBenchProtocolReadyNat",
        schema: "rad-sim.reachable-equilibrium-bench-protocol.v1",
      },
      claimLabels: {
        protocol: "bench protocol artifact",
        reachability: "finite simulator evidence to be tested physically",
        blockedTargets: "graph-topology hypothesis requiring measurement",
        physicalAccuracy: "experimentally unvalidated physical assumption",
      },
      limitations: [
        "The protocol defines measurements but does not certify they were collected.",
        "No physical tolerance is calibrated until hardware scale, backlash, force, and clearance are measured.",
        "A passing bench protocol still does not prove continuous nonlinear controllability.",
      ],
    };
  }

  function exportReachableEquilibriumBenchProtocol(state, options = {}) {
    return JSON.stringify(reachableEquilibriumBenchProtocol(state, options), null, 2);
  }

  function exportReachableEquilibriumBenchProtocolCsv(protocolOrState, options = {}) {
    const protocol =
      protocolOrState?.schema === "rad-sim.reachable-equilibrium-bench-protocol.v1"
        ? protocolOrState
        : protocolOrState?.experiment?.reachableEquilibriumBenchProtocol || reachableEquilibriumBenchProtocol(protocolOrState, options);
    const summary = protocol.summary || {};
    const topology = protocol.topology || {};
    return [
      [
        "schema",
        "protocol_ready",
        "step_id",
        "scope",
        "command_count",
        "observation_cells",
        "target_cells",
        "topology_blocked_target_cells",
        "repeat_count",
        "expected_response",
        "claim_label",
        "missing_evidence",
      ],
      ...(protocol.steps || []).map((step) => [
        protocol.schema,
        summary.benchProtocolReady,
        step.id || "",
        step.scope || "",
        Array.isArray(step.commands) ? step.commands.length : "",
        JSON.stringify(step.observationCells || []),
        JSON.stringify(step.targetCells || []),
        topology.topologyBlockedTargetCells ?? "",
        step.repeatCount ?? "",
        step.expectedResponse || "",
        step.claimLabel || "",
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
      ]),
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function reachableBenchMeasurementModes(step) {
    const scope = String(step?.scope || "");
    const id = String(step?.id || "");
    if (scope === "baseline") return ["baseline"];
    if (scope === "topology-control") return ["topology"];
    if (scope === "group") return ["alpha", "height", "topology"];
    if (id.endsWith("_alpha_response_column")) return ["alpha"];
    if (id.endsWith("_z_response_column")) return ["height"];
    return ["alpha", "height"];
  }

  function reachableBenchTargetSummaryMap(protocol) {
    const out = new Map();
    for (const raw of protocol?.targetSummaries || []) {
      const cell = normalizeContactGraphCell(raw);
      if (!cell) continue;
      out.set(contactGraphKey(cell.row, cell.col), raw);
    }
    return out;
  }

  function reachableEquilibriumBenchResultsTemplate(protocol, options = {}) {
    const targetSummary = reachableBenchTargetSummaryMap(protocol);
    const datasetId = options.datasetId || "reachable-equilibrium-run-001";
    const rows = [];
    for (const step of protocol?.steps || []) {
      const repeatCount = Math.max(1, Math.round(Number(step.repeatCount ?? protocol.repeatCount ?? 1)));
      const applicationModes = step.scope === "group" ? ["simultaneous", "sequenced"] : [String(step.scope || "single")];
      const targets = reachabilityUniqueCells(step.targetCells || []);
      for (let repeatIndex = 0; repeatIndex < repeatCount; repeatIndex += 1) {
        for (const applicationMode of applicationModes) {
          rows.push({
            stepId: step.id || "",
            scope: step.scope || "",
            repeatIndex,
            applicationMode,
            commands: step.commands || [],
            expectedResponse: step.expectedResponse || "",
            measurementMode: reachableBenchMeasurementModes(step),
            targetMeasurements: targets.map((cell) => {
              const summary = targetSummary.get(contactGraphKey(cell.row, cell.col)) || {};
              return {
                row: cell.row,
                col: cell.col,
                predictedAlphaReachable: Boolean(summary.alphaReachable),
                predictedHeightReachable: Boolean(summary.heightReachable),
                predictedTopologyBlocked: Boolean(summary.blockedByTopology),
                expectedTopologyComponentLabel: Number.isInteger(summary.componentLabel) ? summary.componentLabel : -1,
                measuredAlphaDelta: null,
                measuredHeightDelta: null,
                measuredTopologyComponentLabel: null,
                measuredPinHoleSlipMm: null,
                measuredActuatorForceN: null,
                notes: "",
              };
            }),
          });
        }
      }
    }
    return {
      schema: "rad-sim.reachable-equilibrium-bench-results.v1",
      protocolSchema: protocol?.schema || "",
      datasetId,
      datasetRole: options.datasetRole || "bench",
      sourceFileId: options.sourceFileId || `${datasetId}.json`,
      operator: options.operator || "",
      collectedAt: options.collectedAt || "",
      measurementFields: [...REACHABLE_EQUILIBRIUM_BENCH_FIELDS],
      provenanceNotes:
        "Fill measuredAlphaDelta and measuredHeightDelta from the registered cell-center frame; leave predictions unchanged.",
      measurements: rows,
    };
  }

  function exportReachableEquilibriumBenchResultsTemplate(protocolOrState, options = {}) {
    const protocol =
      protocolOrState?.schema === "rad-sim.reachable-equilibrium-bench-protocol.v1"
        ? protocolOrState
        : protocolOrState?.experiment?.reachableEquilibriumBenchProtocol || reachableEquilibriumBenchProtocol(protocolOrState, options);
    return JSON.stringify(reachableEquilibriumBenchResultsTemplate(protocol, options), null, 2);
  }

  function reachableEquilibriumBenchResultsFromJson(text) {
    const parsed = JSON.parse(text);
    if (parsed?.schema !== "rad-sim.reachable-equilibrium-bench-results.v1") {
      throw new Error("unsupported reachable equilibrium bench results schema");
    }
    if (!Array.isArray(parsed.measurements)) throw new Error("reachable equilibrium bench results require measurements");
    return parsed;
  }

  function benchMeasurementNumber(value) {
    const numeric = Number(value);
    return Number.isFinite(numeric) ? numeric : null;
  }

  function benchTargetMeasurementMap(row) {
    const out = new Map();
    for (const target of row?.targetMeasurements || []) {
      const cell = normalizeContactGraphCell(target);
      if (!cell) continue;
      out.set(contactGraphKey(cell.row, cell.col), target);
    }
    return out;
  }

  function compareReachableEquilibriumBenchResults(protocol, results, options = {}) {
    if (results?.schema !== "rad-sim.reachable-equilibrium-bench-results.v1") {
      throw new Error("unsupported reachable equilibrium bench results schema");
    }
    const responseTolerance = Math.max(0, Number(options.responseTolerance ?? 1e-6));
    const blockedTolerance = Math.max(0, Number(options.blockedTolerance ?? 1e-6));
    const groupSequenceTolerance = Math.max(0, Number(options.groupSequenceTolerance ?? 1e-6));
    const stepById = new Map((protocol?.steps || []).map((step) => [String(step.id || ""), step]));
    let targetMeasurementCount = 0;
    let missingMeasurementCount = 0;
    let alphaCheckCount = 0;
    let heightCheckCount = 0;
    let alphaReachabilityMismatchCount = 0;
    let heightReachabilityMismatchCount = 0;
    let topologyLeakageCount = 0;
    let componentMismatchCount = 0;
    let maxAbsMeasuredAlphaDelta = 0;
    let maxAbsMeasuredHeightDelta = 0;
    const rows = [];
    for (const row of results.measurements || []) {
      const step = stepById.get(String(row.stepId || "")) || {};
      const modes = Array.isArray(row.measurementMode) ? row.measurementMode.map(String) : reachableBenchMeasurementModes(step);
      let rowMissing = 0;
      let rowAlphaMismatch = 0;
      let rowHeightMismatch = 0;
      let rowTopologyLeakage = 0;
      let rowComponentMismatch = 0;
      for (const target of row.targetMeasurements || []) {
        targetMeasurementCount += 1;
        const alpha = benchMeasurementNumber(target.measuredAlphaDelta);
        const height = benchMeasurementNumber(target.measuredHeightDelta);
        if (alpha === null) rowMissing += 1;
        else maxAbsMeasuredAlphaDelta = Math.max(maxAbsMeasuredAlphaDelta, Math.abs(alpha));
        if (height === null) rowMissing += 1;
        else maxAbsMeasuredHeightDelta = Math.max(maxAbsMeasuredHeightDelta, Math.abs(height));
        if (modes.includes("baseline")) {
          if (alpha !== null) {
            alphaCheckCount += 1;
            if (Math.abs(alpha) > responseTolerance) rowAlphaMismatch += 1;
          }
          if (height !== null) {
            heightCheckCount += 1;
            if (Math.abs(height) > responseTolerance) rowHeightMismatch += 1;
          }
        }
        if (modes.includes("alpha") && alpha !== null) {
          alphaCheckCount += 1;
          if (Boolean(target.predictedAlphaReachable) !== Math.abs(alpha) > responseTolerance) rowAlphaMismatch += 1;
        }
        if (modes.includes("height") && height !== null) {
          heightCheckCount += 1;
          if (Boolean(target.predictedHeightReachable) !== Math.abs(height) > responseTolerance) rowHeightMismatch += 1;
        }
        if (modes.includes("topology") && target.predictedTopologyBlocked) {
          const leaked =
            (alpha !== null && Math.abs(alpha) > blockedTolerance) ||
            (height !== null && Math.abs(height) > blockedTolerance);
          if (leaked) rowTopologyLeakage += 1;
        }
        if (target.measuredTopologyComponentLabel !== null && target.measuredTopologyComponentLabel !== undefined) {
          const measured = Number(target.measuredTopologyComponentLabel);
          const expected = Number(target.expectedTopologyComponentLabel);
          if (!Number.isInteger(measured) || measured !== expected) rowComponentMismatch += 1;
        }
      }
      missingMeasurementCount += rowMissing;
      alphaReachabilityMismatchCount += rowAlphaMismatch;
      heightReachabilityMismatchCount += rowHeightMismatch;
      topologyLeakageCount += rowTopologyLeakage;
      componentMismatchCount += rowComponentMismatch;
      rows.push({
        stepId: row.stepId || "",
        repeatIndex: row.repeatIndex || 0,
        applicationMode: row.applicationMode || "",
        missingMeasurementCount: rowMissing,
        alphaReachabilityMismatchCount: rowAlphaMismatch,
        heightReachabilityMismatchCount: rowHeightMismatch,
        topologyLeakageCount: rowTopologyLeakage,
        componentMismatchCount: rowComponentMismatch,
      });
    }
    const rowsByMode = new Map();
    for (const row of results.measurements || []) {
      rowsByMode.set(`${row.stepId || ""}|${row.repeatIndex || 0}|${row.applicationMode || ""}`, row);
    }
    let groupSequencePairCount = 0;
    let groupSequenceMaxError = 0;
    for (const step of protocol?.steps || []) {
      if (step.scope !== "group") continue;
      const repeatCount = Math.max(1, Math.round(Number(step.repeatCount ?? protocol.repeatCount ?? 1)));
      for (let repeatIndex = 0; repeatIndex < repeatCount; repeatIndex += 1) {
        const simultaneous = rowsByMode.get(`${step.id || ""}|${repeatIndex}|simultaneous`);
        const sequenced = rowsByMode.get(`${step.id || ""}|${repeatIndex}|sequenced`);
        if (!simultaneous || !sequenced) continue;
        const simTargets = benchTargetMeasurementMap(simultaneous);
        const seqTargets = benchTargetMeasurementMap(sequenced);
        for (const [key, simTarget] of simTargets.entries()) {
          const seqTarget = seqTargets.get(key);
          if (!seqTarget) continue;
          groupSequencePairCount += 1;
          for (const field of ["measuredAlphaDelta", "measuredHeightDelta"]) {
            const simValue = benchMeasurementNumber(simTarget[field]);
            const seqValue = benchMeasurementNumber(seqTarget[field]);
            if (simValue === null || seqValue === null) continue;
            groupSequenceMaxError = Math.max(groupSequenceMaxError, Math.abs(simValue - seqValue));
          }
        }
      }
    }
    const groupSequencePass = groupSequenceMaxError <= groupSequenceTolerance;
    const protocolReady = Boolean(protocol?.summary?.benchProtocolReady);
    const missingEvidence = [];
    if (!protocolReady) missingEvidence.push("reachableEquilibriumBenchProtocol");
    if (!(results.measurements || []).length) missingEvidence.push("resultRows");
    if (targetMeasurementCount <= 0) missingEvidence.push("targetMeasurements");
    if (missingMeasurementCount > 0) missingEvidence.push("completedMeasurements");
    if (alphaCheckCount + heightCheckCount <= 0) missingEvidence.push("reachabilityChecks");
    if (groupSequencePairCount <= 0) missingEvidence.push("groupSequenceChecks");
    const pass =
      missingEvidence.length === 0 &&
      alphaReachabilityMismatchCount === 0 &&
      heightReachabilityMismatchCount === 0 &&
      topologyLeakageCount === 0 &&
      componentMismatchCount === 0 &&
      groupSequencePass;
    return {
      schema: "rad-sim.reachable-equilibrium-bench-comparison.v1",
      method: "comparison of filled reachable-equilibrium bench measurements against simulator reachability flags and topology-blocked controls",
      sourceProtocol: { schema: protocol?.schema || "", ready: protocolReady },
      sourceResults: {
        schema: results.schema,
        datasetId: results.datasetId || "",
        datasetRole: results.datasetRole || "",
        sourceFileId: results.sourceFileId || "",
      },
      metrics: {
        resultRowCount: (results.measurements || []).length,
        targetMeasurementCount,
        missingMeasurementCount,
        alphaCheckCount,
        heightCheckCount,
        alphaReachabilityMismatchCount,
        heightReachabilityMismatchCount,
        topologyLeakageCount,
        componentMismatchCount,
        groupSequencePairCount,
        groupSequenceMaxError,
        maxAbsMeasuredAlphaDelta,
        maxAbsMeasuredHeightDelta,
        responseTolerance,
        blockedTolerance,
        groupSequenceTolerance,
      },
      summary: {
        status: pass ? "reachable-equilibrium-bench-comparison-pass" : "needs-reachable-equilibrium-bench-review",
        benchComparisonPass: pass,
        topologyBlockedLeakagePass: topologyLeakageCount === 0,
        groupSequencePass: groupSequencePass && groupSequencePairCount > 0,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      rows,
      formalization: {
        targetId: "reachable_equilibrium_bench_validation_gate",
        leanStructure: "Mechanics.ReachableEquilibriumBenchValidationNat",
        leanPredicate: "reachableEquilibriumBenchValidationReadyNat",
        schema: "rad-sim.reachable-equilibrium-bench-comparison.v1",
      },
      claimLabels: {
        comparison: "bench result comparison artifact",
        reachability: "simulator-derived empirical law compared to measurement",
        topologyLeakage: "experimentally measured model-mismatch indicator",
        physicalAccuracy: "requires calibrated hardware and repeated trials",
      },
      limitations: [
        "Reachable flags are qualitative nonzero-response checks, not fitted amplitude laws.",
        "Topology leakage can indicate unmodeled bypass coupling, fixture compliance, or measurement noise.",
        "Passing this comparison does not prove nonlinear controllability or contact mechanics.",
      ],
    };
  }

  function exportReachableEquilibriumBenchComparison(protocol, results, options = {}) {
    return JSON.stringify(compareReachableEquilibriumBenchResults(protocol, results, options), null, 2);
  }

  function exportReachableEquilibriumBenchComparisonCsv(report) {
    const summary = report.summary || {};
    const metrics = report.metrics || {};
    return [
      [
        "schema",
        "comparison_pass",
        "step_id",
        "repeat_index",
        "application_mode",
        "missing_measurements",
        "alpha_mismatches",
        "height_mismatches",
        "topology_leakage",
        "component_mismatches",
        "group_sequence_max_error",
        "missing_evidence",
      ],
      ...(report.rows || []).map((row) => [
        report.schema,
        summary.benchComparisonPass,
        row.stepId || "",
        row.repeatIndex ?? "",
        row.applicationMode || "",
        row.missingMeasurementCount ?? "",
        row.alphaReachabilityMismatchCount ?? "",
        row.heightReachabilityMismatchCount ?? "",
        row.topologyLeakageCount ?? "",
        row.componentMismatchCount ?? "",
        metrics.groupSequenceMaxError ?? "",
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
      ]),
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function benchNumericStats(values, confidenceSigma) {
    const finite = (values || []).map(Number).filter(Number.isFinite);
    if (!finite.length) {
      return { sampleCount: 0, mean: null, sampleStd: null, meanAbs: null, maxAbs: null, uncertaintyHalfWidth: null };
    }
    const mean = finite.reduce((sum, value) => sum + value, 0) / finite.length;
    const variance =
      finite.length > 1
        ? finite.reduce((sum, value) => sum + (value - mean) ** 2, 0) / (finite.length - 1)
        : 0;
    const sampleStd = Math.sqrt(variance);
    return {
      sampleCount: finite.length,
      mean,
      sampleStd,
      meanAbs: finite.reduce((sum, value) => sum + Math.abs(value), 0) / finite.length,
      maxAbs: Math.max(...finite.map(Math.abs)),
      uncertaintyHalfWidth: Number(confidenceSigma) * sampleStd / Math.sqrt(finite.length),
    };
  }

  function benchQualitativeResidual(meanValue, predictedReachable, tolerance) {
    const value = benchMeasurementNumber(meanValue);
    if (value === null) return tolerance;
    const observed = Math.abs(value) > tolerance;
    if (observed === Boolean(predictedReachable)) return 0;
    return observed ? Math.abs(value) : tolerance;
  }

  function reachableEquilibriumAmplitudeCalibrationReport(protocol, results, options = {}) {
    const responseTolerance = Math.max(0, Number(options.responseTolerance ?? 1e-6));
    const blockedTolerance = Math.max(0, Number(options.blockedTolerance ?? 1e-6));
    const groupSequenceTolerance = Math.max(0, Number(options.groupSequenceTolerance ?? 1e-6));
    const confidenceSigma = Math.max(0, Number(options.confidenceSigma ?? 2));
    const comparison =
      options.comparisonReport ||
      compareReachableEquilibriumBenchResults(protocol, results, {
        responseTolerance,
        blockedTolerance,
        groupSequenceTolerance,
      });
    const grouped = new Map();
    const stepById = new Map((protocol?.steps || []).map((step) => [String(step.id || ""), step]));
    for (const row of results?.measurements || []) {
      const step = stepById.get(String(row.stepId || "")) || {};
      const modes = Array.isArray(row.measurementMode) ? row.measurementMode.map(String) : reachableBenchMeasurementModes(step);
      for (const target of row.targetMeasurements || []) {
        const cell = normalizeContactGraphCell(target);
        if (!cell) continue;
        const key = `${row.stepId || ""}|${row.applicationMode || ""}|${cell.row}|${cell.col}`;
        if (!grouped.has(key)) {
          grouped.set(key, {
            stepId: String(row.stepId || ""),
            scope: String(row.scope || ""),
            applicationMode: String(row.applicationMode || ""),
            row: cell.row,
            col: cell.col,
            measurementMode: modes,
            predictedAlphaReachable: Boolean(target.predictedAlphaReachable),
            predictedHeightReachable: Boolean(target.predictedHeightReachable),
            predictedTopologyBlocked: Boolean(target.predictedTopologyBlocked),
            expectedTopologyComponentLabel: Number.isInteger(target.expectedTopologyComponentLabel) ? target.expectedTopologyComponentLabel : -1,
            alphaValues: [],
            heightValues: [],
            slipValues: [],
            forceValues: [],
          });
        }
        const entry = grouped.get(key);
        const alpha = benchMeasurementNumber(target.measuredAlphaDelta);
        const height = benchMeasurementNumber(target.measuredHeightDelta);
        const slip = benchMeasurementNumber(target.measuredPinHoleSlipMm);
        const force = benchMeasurementNumber(target.measuredActuatorForceN);
        if (alpha !== null) entry.alphaValues.push(alpha);
        if (height !== null) entry.heightValues.push(height);
        if (slip !== null) entry.slipValues.push(slip);
        if (force !== null) entry.forceValues.push(force);
      }
    }
    const amplitudeEstimates = [];
    const residualField = [];
    const topologyResponseBands = [];
    for (const raw of grouped.values()) {
      const alpha = benchNumericStats(raw.alphaValues, confidenceSigma);
      const height = benchNumericStats(raw.heightValues, confidenceSigma);
      const pinHoleSlipMm = benchNumericStats(raw.slipValues, confidenceSigma);
      const actuatorForceN = benchNumericStats(raw.forceValues, confidenceSigma);
      amplitudeEstimates.push({
        stepId: raw.stepId,
        scope: raw.scope,
        applicationMode: raw.applicationMode,
        row: raw.row,
        col: raw.col,
        measurementMode: raw.measurementMode,
        predictedAlphaReachable: raw.predictedAlphaReachable,
        predictedHeightReachable: raw.predictedHeightReachable,
        predictedTopologyBlocked: raw.predictedTopologyBlocked,
        expectedTopologyComponentLabel: raw.expectedTopologyComponentLabel,
        alpha,
        height,
        pinHoleSlipMm,
        actuatorForceN,
      });
      const alphaResidual = benchQualitativeResidual(alpha.mean, raw.predictedAlphaReachable, responseTolerance);
      const heightResidual = benchQualitativeResidual(height.mean, raw.predictedHeightReachable, responseTolerance);
      residualField.push({
        stepId: raw.stepId,
        applicationMode: raw.applicationMode,
        row: raw.row,
        col: raw.col,
        alphaQualitativeResidual: alphaResidual,
        heightQualitativeResidual: heightResidual,
        combinedResidual: Math.max(alphaResidual, heightResidual),
      });
      const alphaBand = alpha.meanAbs === null ? null : alpha.meanAbs + (alpha.uncertaintyHalfWidth || 0);
      const heightBand = height.meanAbs === null ? null : height.meanAbs + (height.uncertaintyHalfWidth || 0);
      topologyResponseBands.push({
        stepId: raw.stepId,
        applicationMode: raw.applicationMode,
        row: raw.row,
        col: raw.col,
        predictedTopologyBlocked: raw.predictedTopologyBlocked,
        alphaResponseBand: alphaBand,
        heightResponseBand: heightBand,
        blockedTolerance,
        blockedLeakageBandPass:
          !raw.predictedTopologyBlocked ||
          ((alphaBand || 0) <= blockedTolerance && (heightBand || 0) <= blockedTolerance),
      });
    }
    const rowsByMode = new Map();
    for (const row of results?.measurements || []) {
      rowsByMode.set(`${row.stepId || ""}|${row.repeatIndex || 0}|${row.applicationMode || ""}`, row);
    }
    const groupSequenceResiduals = [];
    for (const step of protocol?.steps || []) {
      if (step.scope !== "group") continue;
      const repeatCount = Math.max(1, Math.round(Number(step.repeatCount ?? protocol.repeatCount ?? 1)));
      for (let repeatIndex = 0; repeatIndex < repeatCount; repeatIndex += 1) {
        const simultaneous = rowsByMode.get(`${step.id || ""}|${repeatIndex}|simultaneous`);
        const sequenced = rowsByMode.get(`${step.id || ""}|${repeatIndex}|sequenced`);
        if (!simultaneous || !sequenced) continue;
        const simTargets = benchTargetMeasurementMap(simultaneous);
        const seqTargets = benchTargetMeasurementMap(sequenced);
        for (const [key, simTarget] of simTargets.entries()) {
          const seqTarget = seqTargets.get(key);
          if (!seqTarget) continue;
          const cell = normalizeContactGraphCell(simTarget);
          if (!cell) continue;
          const simAlpha = benchMeasurementNumber(simTarget.measuredAlphaDelta);
          const seqAlpha = benchMeasurementNumber(seqTarget.measuredAlphaDelta);
          const simHeight = benchMeasurementNumber(simTarget.measuredHeightDelta);
          const seqHeight = benchMeasurementNumber(seqTarget.measuredHeightDelta);
          const alphaError = simAlpha === null || seqAlpha === null ? 0 : Math.abs(simAlpha - seqAlpha);
          const heightError = simHeight === null || seqHeight === null ? 0 : Math.abs(simHeight - seqHeight);
          groupSequenceResiduals.push({
            stepId: step.id || "",
            repeatIndex,
            row: cell.row,
            col: cell.col,
            alphaError,
            heightError,
            combinedError: Math.max(alphaError, heightError),
          });
        }
      }
    }
    const repeatedTrialGroupCount = amplitudeEstimates.filter(
      (estimate) => Math.max(estimate.alpha.sampleCount, estimate.height.sampleCount) >= 2
    ).length;
    const failedTopologyBandCount = topologyResponseBands.filter((band) => !band.blockedLeakageBandPass).length;
    const maxQualitativeResidual = residualField.reduce((max, item) => Math.max(max, item.combinedResidual), 0);
    const maxGroupSequenceResidual = groupSequenceResiduals.reduce((max, item) => Math.max(max, item.combinedError), 0);
    const comparisonPass = Boolean(comparison?.summary?.benchComparisonPass);
    const missingEvidence = [];
    if (!comparisonPass) missingEvidence.push("benchComparison");
    if (!amplitudeEstimates.length) missingEvidence.push("amplitudeEstimates");
    if (repeatedTrialGroupCount <= 0) missingEvidence.push("repeatedTrials");
    if (!residualField.length) missingEvidence.push("residualField");
    if (!topologyResponseBands.length) missingEvidence.push("topologyBands");
    if (failedTopologyBandCount > 0) missingEvidence.push("topologyLeakageBand");
    if (!groupSequenceResiduals.length) missingEvidence.push("groupSequenceResiduals");
    if (maxGroupSequenceResidual > groupSequenceTolerance) missingEvidence.push("groupSequenceTolerance");
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.reachable-equilibrium-amplitude-calibration.v1",
      method: "repeated-trial amplitude calibration report over filled reachable-equilibrium bench measurements",
      sourceProtocol: { schema: protocol?.schema || "", stepCount: (protocol?.steps || []).length },
      sourceResults: {
        schema: results?.schema || "",
        datasetId: results?.datasetId || "",
        datasetRole: results?.datasetRole || "",
        sourceFileId: results?.sourceFileId || "",
      },
      sourceComparison: { schema: comparison?.schema || "", benchComparisonPass: comparisonPass },
      parameters: { responseTolerance, blockedTolerance, groupSequenceTolerance, confidenceSigma },
      amplitudeEstimates,
      residualField,
      topologyResponseBands,
      groupSequenceResiduals,
      metrics: {
        amplitudeEstimateCount: amplitudeEstimates.length,
        repeatedTrialGroupCount,
        residualFieldCellCount: residualField.length,
        topologyBandCount: topologyResponseBands.length,
        failedTopologyBandCount,
        groupSequenceResidualCount: groupSequenceResiduals.length,
        maxQualitativeResidual,
        maxGroupSequenceResidual,
      },
      summary: {
        status: ready ? "reachable-equilibrium-amplitude-calibration-ready" : "needs-reachable-equilibrium-amplitude-review",
        amplitudeCalibrationReady: ready,
        topologyLeakageBandsPass: failedTopologyBandCount === 0,
        groupSequenceResidualPass: maxGroupSequenceResidual <= groupSequenceTolerance && groupSequenceResiduals.length > 0,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "reachable_equilibrium_amplitude_calibration_gate",
        leanStructure: "Mechanics.ReachableEquilibriumAmplitudeCalibrationNat",
        leanPredicate: "reachableEquilibriumAmplitudeCalibrationReadyNat",
        schema: "rad-sim.reachable-equilibrium-amplitude-calibration.v1",
      },
      claimLabels: {
        amplitudes: "bench-measured empirical law",
        uncertainty: "finite repeated-trial statistic",
        topologyBands: "experimentally measured leakage bound",
        physicalAccuracy: "requires external hardware calibration and uncertainty analysis",
      },
      limitations: [
        "The report estimates measured amplitudes but does not derive a constitutive mechanics law.",
        "Uncertainty bands use a simple sigma multiplier and are not a full statistical model.",
        "Group residuals compare two command orderings but do not prove operator commutativity.",
      ],
    };
  }

  function exportReachableEquilibriumAmplitudeCalibration(protocol, results, options = {}) {
    return JSON.stringify(reachableEquilibriumAmplitudeCalibrationReport(protocol, results, options), null, 2);
  }

  function exportReachableEquilibriumAmplitudeCalibrationCsv(report) {
    const summary = report.summary || {};
    const metrics = report.metrics || {};
    return [
      [
        "schema",
        "calibration_ready",
        "step_id",
        "application_mode",
        "row",
        "col",
        "alpha_samples",
        "alpha_mean",
        "alpha_uncertainty",
        "height_samples",
        "height_mean",
        "height_uncertainty",
        "predicted_topology_blocked",
        "max_group_sequence_residual",
        "missing_evidence",
      ],
      ...(report.amplitudeEstimates || []).map((estimate) => [
        report.schema,
        summary.amplitudeCalibrationReady,
        estimate.stepId || "",
        estimate.applicationMode || "",
        estimate.row ?? "",
        estimate.col ?? "",
        estimate.alpha?.sampleCount ?? "",
        estimate.alpha?.mean ?? "",
        estimate.alpha?.uncertaintyHalfWidth ?? "",
        estimate.height?.sampleCount ?? "",
        estimate.height?.mean ?? "",
        estimate.height?.uncertaintyHalfWidth ?? "",
        estimate.predictedTopologyBlocked,
        metrics.maxGroupSequenceResidual ?? "",
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
      ]),
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function amplitudeSectionValues(report, channel) {
    const reachableKey = channel === "alpha" ? "predictedAlphaReachable" : "predictedHeightReachable";
    return (report?.amplitudeEstimates || [])
      .filter((estimate) => Boolean(estimate?.[reachableKey]))
      .map((estimate) => Number(estimate?.[channel]?.meanAbs))
      .filter(Number.isFinite);
  }

  function amplitudeUncertaintyValues(report) {
    const values = [];
    for (const estimate of report?.amplitudeEstimates || []) {
      for (const channel of ["alpha", "height"]) {
        const value = Number(estimate?.[channel]?.uncertaintyHalfWidth);
        if (Number.isFinite(value)) values.push(value);
      }
    }
    return values;
  }

  function amplitudeTopologyBandValues(report) {
    const values = [];
    for (const band of report?.topologyResponseBands || []) {
      if (!band?.predictedTopologyBlocked) continue;
      const alpha = Number.isFinite(Number(band.alphaResponseBand)) ? Number(band.alphaResponseBand) : 0;
      const height = Number.isFinite(Number(band.heightResponseBand)) ? Number(band.heightResponseBand) : 0;
      values.push(Math.max(alpha, height));
    }
    return values;
  }

  function boundedProfileProposal(name, proposed, lower, upper, sampleCount, { minSamples, claimLabel, reason }) {
    const finite = Number.isFinite(Number(proposed));
    const bounded = finite && Number(proposed) >= lower && Number(proposed) <= upper;
    const safeToApply = bounded && sampleCount >= minSamples;
    return {
      name,
      configField: null,
      current: null,
      proposed: finite ? Number(proposed) : null,
      lowerBound: lower,
      upperBound: upper,
      sampleCount,
      bounded,
      safeToApply,
      claimLabel,
      reason: safeToApply ? reason : `${reason}; requires finite bounded value and at least ${minSamples} samples`,
    };
  }

  function reachableEquilibriumEmpiricalProfileFromAmplitude(amplitudeReport, options = {}) {
    const safetyFactor = Math.max(1, Number(options.safetyFactor ?? 1.5));
    const minSamples = Math.max(0, Math.round(Number(options.minSamples ?? 1)));
    const requireHoldout = Boolean(options.requireHoldout);
    const alphaValues = amplitudeSectionValues(amplitudeReport, "alpha");
    const heightValues = amplitudeSectionValues(amplitudeReport, "height");
    const uncertaintyValues = amplitudeUncertaintyValues(amplitudeReport);
    const topologyValues = amplitudeTopologyBandValues(amplitudeReport);
    const mean = (values) => (values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : null);
    const alphaScale = mean(alphaValues);
    const heightScale = mean(heightValues);
    const topologyLeakageTolerance = topologyValues.length ? safetyFactor * Math.max(...topologyValues) : 0;
    const groupSequenceResidual = Number(amplitudeReport?.metrics?.maxGroupSequenceResidual);
    const groupSequenceTolerance = Number.isFinite(groupSequenceResidual) ? safetyFactor * groupSequenceResidual : null;
    const uncertaintyBudget = uncertaintyValues.length ? safetyFactor * Math.max(...uncertaintyValues) : null;
    const groupSequenceResidualCount = Number(amplitudeReport?.metrics?.groupSequenceResidualCount || 0);
    const proposals = [
      boundedProfileProposal("alphaResponseScale", alphaScale, 0, 10, alphaValues.length, {
        minSamples,
        claimLabel: "bench-measured empirical alpha response scale",
        reason: "mean absolute measured alpha response over reachable targets",
      }),
      boundedProfileProposal("heightResponseScale", heightScale, 0, 10, heightValues.length, {
        minSamples,
        claimLabel: "bench-measured empirical height response scale",
        reason: "mean absolute measured height response over reachable targets",
      }),
      boundedProfileProposal("topologyLeakageTolerance", topologyLeakageTolerance, 0, 10, topologyValues.length, {
        minSamples: 0,
        claimLabel: "bench-measured empirical topology leakage bound",
        reason: "safety-factor-scaled maximum blocked-target response band",
      }),
      boundedProfileProposal("groupSequenceTolerance", groupSequenceTolerance, 0, 10, groupSequenceResidualCount, {
        minSamples,
        claimLabel: "bench-measured empirical group ordering tolerance",
        reason: "safety-factor-scaled maximum simultaneous/sequenced residual",
      }),
      boundedProfileProposal("uncertaintyBudget", uncertaintyBudget, 0, 10, uncertaintyValues.length, {
        minSamples,
        claimLabel: "finite repeated-trial uncertainty metadata",
        reason: "safety-factor-scaled maximum alpha/height uncertainty half-width",
      }),
    ];
    let holdoutValidation = {
      attached: Boolean(options.holdoutAmplitudeReport),
      schema: options.holdoutAmplitudeReport?.schema || "",
      holdoutReady: false,
      holdoutPass: !requireHoldout && !options.holdoutAmplitudeReport,
      reason: "holdout not required and not attached",
    };
    if (options.holdoutAmplitudeReport) {
      const holdout = options.holdoutAmplitudeReport;
      const holdoutGroup = Number(holdout?.metrics?.maxGroupSequenceResidual);
      const groupLimit = groupSequenceTolerance ?? 0;
      const holdoutReady = Boolean(holdout?.summary?.amplitudeCalibrationReady);
      const holdoutTopology = Boolean(holdout?.summary?.topologyLeakageBandsPass);
      const holdoutPass = holdoutReady && holdoutTopology && Number.isFinite(holdoutGroup) && holdoutGroup <= groupLimit;
      holdoutValidation = {
        attached: true,
        schema: holdout.schema || "",
        holdoutReady,
        holdoutPass,
        maxGroupSequenceResidual: Number.isFinite(holdoutGroup) ? holdoutGroup : null,
        groupSequenceTolerance: groupLimit,
        topologyLeakageBandsPass: holdoutTopology,
        reason: holdoutPass
          ? "holdout amplitude report fits within empirical profile tolerances"
          : "holdout amplitude report is missing readiness or exceeds empirical tolerances",
      };
    }
    const safeProposalCount = proposals.filter((proposal) => proposal.safeToApply).length;
    const boundedProposalCount = proposals.filter((proposal) => proposal.bounded).length;
    const missingEvidence = [];
    if (!amplitudeReport?.summary?.amplitudeCalibrationReady) missingEvidence.push("amplitudeCalibration");
    if (!alphaValues.length) missingEvidence.push("alphaResponseScale");
    if (!heightValues.length) missingEvidence.push("heightResponseScale");
    if (!Number.isFinite(groupSequenceResidual)) missingEvidence.push("groupSequenceResidual");
    if (!uncertaintyValues.length) missingEvidence.push("uncertaintyBudget");
    if (safeProposalCount <= 0) missingEvidence.push("safeBoundedProposals");
    if (requireHoldout && !holdoutValidation.holdoutPass) missingEvidence.push("holdoutValidation");
    const ready = missingEvidence.length === 0;
    return {
      schema: "rad-sim.reachable-equilibrium-empirical-profile.v1",
      sourceReportSchema: amplitudeReport?.schema || "",
      method: "bounded empirical profile proposals derived from repeated-trial reachable-equilibrium amplitude calibration",
      parameters: { safetyFactor, minSamples, requireHoldout },
      recommendedUpdates: proposals,
      holdoutValidation,
      metrics: {
        proposalCount: proposals.length,
        boundedProposalCount,
        safeProposalCount,
        alphaSampleCount: alphaValues.length,
        heightSampleCount: heightValues.length,
        topologyLeakageBandCount: topologyValues.length,
        uncertaintySampleCount: uncertaintyValues.length,
        groupSequenceResidualCount,
      },
      summary: {
        status: ready ? "reachable-equilibrium-empirical-profile-ready" : "needs-reachable-equilibrium-profile-review",
        empiricalProfileReady: ready,
        safeBoundedProposalCount: safeProposalCount,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "reachable_equilibrium_empirical_profile_gate",
        leanStructure: "Mechanics.ReachableEquilibriumEmpiricalProfileNat",
        leanPredicate: "reachableEquilibriumEmpiricalProfileReadyNat",
        schema: "rad-sim.reachable-equilibrium-empirical-profile.v1",
      },
      claimLabels: {
        profile: "bench-measured empirical law",
        safeProposals: "bounded finite-sample update proposals",
        holdout: "optional empirical validation hook",
        physicalAccuracy: "not a constitutive mechanics proof",
      },
      limitations: [
        "Profile proposals are calibration metadata; no LatticeConfig field is mutated in v1.",
        "Bounds prevent nonsensical profile values but do not prove physical correctness.",
        "Holdout validation checks tolerance consistency, not statistical independence unless the lab supplied independent data.",
      ],
    };
  }

  function exportReachableEquilibriumEmpiricalProfile(amplitudeReport, options = {}) {
    return JSON.stringify(reachableEquilibriumEmpiricalProfileFromAmplitude(amplitudeReport, options), null, 2);
  }

  function exportReachableEquilibriumEmpiricalProfileCsv(profile) {
    const summary = profile.summary || {};
    return [
      [
        "schema",
        "profile_ready",
        "name",
        "proposed",
        "lower_bound",
        "upper_bound",
        "sample_count",
        "bounded",
        "safe_to_apply",
        "claim_label",
        "missing_evidence",
      ],
      ...(profile.recommendedUpdates || []).map((update) => [
        profile.schema,
        summary.empiricalProfileReady,
        update.name || "",
        update.proposed ?? "",
        update.lowerBound ?? "",
        update.upperBound ?? "",
        update.sampleCount ?? "",
        update.bounded,
        update.safeToApply,
        update.claimLabel || "",
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
      ]),
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function exportCalibrationModelProfileSelection(comparisons, options = {}) {
    const resolved = Array.isArray(comparisons)
      ? comparisons
      : (comparisons?.experiment?.calibrationModelProfileHistory || [])
          .map((entry) => entry.residualComparison)
          .filter(Boolean);
    return JSON.stringify(selectCalibrationModelProfile(resolved, options), null, 2);
  }

  function calibrationComparisonReport(state) {
    const comparison = state.experiment?.calibrationComparison;
    if (!comparison) throw new Error("no imported calibration comparison is available");
    const summary = state.experiment?.calibrationComparisonSummary || summarizeCalibrationComparison(comparison);
    const parameterEstimates = calibrationParameterEstimates(
      state,
      state.experiment?.calibrationResults || {},
      comparison,
      summary
    );
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
      parameterEstimates,
      modelProfile: calibrationModelProfile(state, parameterEstimates, {
        sourceReportSchema: "rad-sim.calibration-comparison-report.v1",
      }),
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
  RAD.buildResponseMatrix = buildResponseMatrix;
  RAD.exportResponseMatrix = exportResponseMatrix;
  RAD.topologyExperimentReport = topologyExperimentReport;
  RAD.exportTopologyExperimentReport = exportTopologyExperimentReport;
  RAD.programmableDiscontinuityReport = programmableDiscontinuityReport;
  RAD.exportProgrammableDiscontinuityReport = exportProgrammableDiscontinuityReport;
  RAD.formalizationTargetManifest = formalizationTargetManifest;
  RAD.exportFormalizationTargetManifest = exportFormalizationTargetManifest;
  RAD.calibrationExperimentProtocol = calibrationExperimentProtocol;
  RAD.exportCalibrationExperimentProtocol = exportCalibrationExperimentProtocol;
  RAD.responseAtlas = responseAtlas;
  RAD.exportResponseAtlas = exportResponseAtlas;
  RAD.responseAtlasSweep = responseAtlasSweep;
  RAD.exportResponseAtlasSweep = exportResponseAtlasSweep;
  RAD.calibrationExperimentResultsTemplate = calibrationExperimentResultsTemplate;
  RAD.exportCalibrationExperimentResultsTemplate = exportCalibrationExperimentResultsTemplate;
  RAD.calibrationBenchNotebook = calibrationBenchNotebook;
  RAD.exportCalibrationBenchNotebook = exportCalibrationBenchNotebook;
  RAD.exportCalibrationBenchNotebookCsv = exportCalibrationBenchNotebookCsv;
  RAD.calibrationBenchPacket = calibrationBenchPacket;
  RAD.exportCalibrationBenchPacket = exportCalibrationBenchPacket;
  RAD.calibrationBenchExecutionValidation = calibrationBenchExecutionValidation;
  RAD.exportCalibrationBenchExecutionValidation = exportCalibrationBenchExecutionValidation;
  RAD.exportCalibrationBenchExecutionValidationCsv = exportCalibrationBenchExecutionValidationCsv;
  RAD.physicalValidationReadinessReport = physicalValidationReadinessReport;
  RAD.exportPhysicalValidationReadiness = exportPhysicalValidationReadiness;
  RAD.exportPhysicalValidationReadinessCsv = exportPhysicalValidationReadinessCsv;
  RAD.contactStateAbstractionReport = contactStateAbstractionReport;
  RAD.exportContactStateAbstraction = exportContactStateAbstraction;
  RAD.exportContactStateAbstractionCsv = exportContactStateAbstractionCsv;
  RAD.contactGraphConsistencyReport = contactGraphConsistencyReport;
  RAD.exportContactGraphConsistency = exportContactGraphConsistency;
  RAD.exportContactGraphConsistencyCsv = exportContactGraphConsistencyCsv;
  RAD.physicalRealizationMapReport = physicalRealizationMapReport;
  RAD.exportPhysicalRealizationMap = exportPhysicalRealizationMap;
  RAD.exportPhysicalRealizationMapCsv = exportPhysicalRealizationMapCsv;
  RAD.externalPhysicsEngineAuditReport = externalPhysicsEngineAuditReport;
  RAD.exportExternalPhysicsEngineAudit = exportExternalPhysicsEngineAudit;
  RAD.exportExternalPhysicsEngineAuditCsv = exportExternalPhysicsEngineAuditCsv;
  RAD.mujocoModelExportReport = mujocoModelExportReport;
  RAD.exportMujocoModelXml = exportMujocoModelXml;
  RAD.exportMujocoModelReport = exportMujocoModelReport;
  RAD.mujocoPinHoleContactGeometryReport = mujocoPinHoleContactGeometryReport;
  RAD.exportMujocoPinHoleContactGeometry = exportMujocoPinHoleContactGeometry;
  RAD.exportMujocoPinHoleContactGeometryCsv = exportMujocoPinHoleContactGeometryCsv;
  RAD.mujocoContactParameterReport = mujocoContactParameterReport;
  RAD.exportMujocoContactParameter = exportMujocoContactParameter;
  RAD.exportMujocoContactParameterCsv = exportMujocoContactParameterCsv;
  RAD.contactParameterCalibrationPacket = contactParameterCalibrationPacket;
  RAD.exportContactParameterCalibrationPacket = exportContactParameterCalibrationPacket;
  RAD.exportContactParameterCalibrationPacketCsv = exportContactParameterCalibrationPacketCsv;
  RAD.contactParameterCalibrationResultsTemplate = contactParameterCalibrationResultsTemplate;
  RAD.exportContactParameterCalibrationResultsTemplate = exportContactParameterCalibrationResultsTemplate;
  RAD.contactParameterCalibrationResultsFromJson = contactParameterCalibrationResultsFromJson;
  RAD.compareContactParameterCalibrationResults = compareContactParameterCalibrationResults;
  RAD.exportContactParameterBenchValidation = exportContactParameterBenchValidation;
  RAD.exportContactParameterBenchValidationCsv = exportContactParameterBenchValidationCsv;
  RAD.contactParameterIntervalCalibrationReport = contactParameterIntervalCalibrationReport;
  RAD.exportContactParameterIntervalCalibration = exportContactParameterIntervalCalibration;
  RAD.exportContactParameterIntervalCalibrationCsv = exportContactParameterIntervalCalibrationCsv;
  RAD.mujocoExternalRunReport = mujocoExternalRunReport;
  RAD.exportMujocoExternalRun = exportMujocoExternalRun;
  RAD.mujocoExternalComparisonReport = mujocoExternalComparisonReport;
  RAD.exportMujocoExternalComparison = exportMujocoExternalComparison;
  RAD.exportMujocoExternalComparisonCsv = exportMujocoExternalComparisonCsv;
  RAD.equilibriumRelationReport = equilibriumRelationReport;
  RAD.exportEquilibriumRelation = exportEquilibriumRelation;
  RAD.exportEquilibriumRelationCsv = exportEquilibriumRelationCsv;
  RAD.reachableEquilibriumControllabilityReport = reachableEquilibriumControllabilityReport;
  RAD.exportReachableEquilibriumControllability = exportReachableEquilibriumControllability;
  RAD.exportReachableEquilibriumControllabilityCsv = exportReachableEquilibriumControllabilityCsv;
  RAD.reachableEquilibriumBenchProtocol = reachableEquilibriumBenchProtocol;
  RAD.exportReachableEquilibriumBenchProtocol = exportReachableEquilibriumBenchProtocol;
  RAD.exportReachableEquilibriumBenchProtocolCsv = exportReachableEquilibriumBenchProtocolCsv;
  RAD.reachableEquilibriumBenchResultsTemplate = reachableEquilibriumBenchResultsTemplate;
  RAD.exportReachableEquilibriumBenchResultsTemplate = exportReachableEquilibriumBenchResultsTemplate;
  RAD.reachableEquilibriumBenchResultsFromJson = reachableEquilibriumBenchResultsFromJson;
  RAD.compareReachableEquilibriumBenchResults = compareReachableEquilibriumBenchResults;
  RAD.exportReachableEquilibriumBenchComparison = exportReachableEquilibriumBenchComparison;
  RAD.exportReachableEquilibriumBenchComparisonCsv = exportReachableEquilibriumBenchComparisonCsv;
  RAD.reachableEquilibriumAmplitudeCalibrationReport = reachableEquilibriumAmplitudeCalibrationReport;
  RAD.exportReachableEquilibriumAmplitudeCalibration = exportReachableEquilibriumAmplitudeCalibration;
  RAD.exportReachableEquilibriumAmplitudeCalibrationCsv = exportReachableEquilibriumAmplitudeCalibrationCsv;
  RAD.reachableEquilibriumEmpiricalProfileFromAmplitude = reachableEquilibriumEmpiricalProfileFromAmplitude;
  RAD.exportReachableEquilibriumEmpiricalProfile = exportReachableEquilibriumEmpiricalProfile;
  RAD.exportReachableEquilibriumEmpiricalProfileCsv = exportReachableEquilibriumEmpiricalProfileCsv;
  RAD.compareCalibrationExperimentResults = compareCalibrationExperimentResults;
  RAD.summarizeCalibrationComparison = summarizeCalibrationComparison;
  RAD.calibrationParameterEstimates = calibrationParameterEstimates;
  RAD.calibrationModelProfile = calibrationModelProfile;
  RAD.applyCalibrationModelProfile = applyCalibrationModelProfile;
  RAD.exportCalibrationModelProfile = exportCalibrationModelProfile;
  RAD.calibrationModelProfileResidualComparison = calibrationModelProfileResidualComparison;
  RAD.calibrationModelProfileHoldoutValidation = calibrationModelProfileHoldoutValidation;
  RAD.selectCalibrationModelProfile = selectCalibrationModelProfile;
  RAD.exportCalibrationModelProfileResidualComparison = exportCalibrationModelProfileResidualComparison;
  RAD.exportCalibrationModelProfileHoldoutValidation = exportCalibrationModelProfileHoldoutValidation;
  RAD.exportCalibrationModelProfileHoldoutValidationCsv = exportCalibrationModelProfileHoldoutValidationCsv;
  RAD.exportCalibrationModelProfileSelection = exportCalibrationModelProfileSelection;
  RAD.calibrationComparisonReport = calibrationComparisonReport;
  RAD.exportCalibrationComparisonReport = exportCalibrationComparisonReport;
  RAD.physicalPreviewComparison = physicalPreviewComparison;
  RAD.springHingePhysicalPreviewReportDescriptor = springHingePhysicalPreviewReportDescriptor;
  RAD.responseDecayProfile = responseDecayProfile;
})();
