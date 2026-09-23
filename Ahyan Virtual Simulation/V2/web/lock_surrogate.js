(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  function dataset() {
    return window.RAD_LOCK_DATASET || null;
  }

  function knownFailedLockCells(cells) {
    const key = [...cells].sort((a, b) => a - b).join(",");
    const known = {
      "1,2,3": [2],
      "4,5,6": [5],
      "7,8,9": [8],
      "10,11,12": [10, 11],
    };
    if (known[key]) return [...known[key]];
    const ordered = [...cells].sort((a, b) => a - b);
    const failed = [];
    let start = 0;
    while (start < ordered.length) {
      let end = start + 1;
      while (end < ordered.length && ordered[end] === ordered[end - 1] + 1) end += 1;
      const run = ordered.slice(start, end);
      if (run.length >= 3) failed.push(...run.slice(1, -1));
      start = end;
    }
    return failed;
  }

  function realizedLockCells(cells) {
    const failed = new Set(knownFailedLockCells(cells));
    return cells.filter((cell) => !failed.has(cell));
  }

  function featureVector(cells, angles, stateIndex, strandCells = 12) {
    const nominal = Array(strandCells).fill(0);
    const realized = Array(strandCells).fill(0);
    const angleVector = Array(strandCells).fill(0);
    const realizedCells = realizedLockCells(cells);
    for (const cell of cells) {
      const index = Number(cell) - 1;
      if (index >= 0 && index < strandCells) nominal[index] = 1;
    }
    for (const cell of realizedCells) {
      const index = Number(cell) - 1;
      if (index >= 0 && index < strandCells) realized[index] = 1;
    }
    for (let i = 0; i < cells.length; i += 1) {
      const index = Number(cells[i]) - 1;
      if (index >= 0 && index < strandCells) angleVector[index] = Number(angles[i] || 30) / 40;
    }
    return [
      ...nominal,
      ...realized,
      ...nominal.map((value, index) => value - realized[index]),
      ...angleVector,
      Number(stateIndex) || 1,
    ];
  }

  function rowLockCells(state, row) {
    const cells = [];
    for (let c = 0; c < state.grid.cols; c += 1) {
      if (state.cells.locked?.[row]?.[c]) cells.push(c + 1);
    }
    return cells;
  }

  function samePattern(a, b) {
    if (a.length !== b.length) return false;
    return a.every((value, index) => value === b[index]);
  }

  function repeatedRowPattern(state) {
    if (state.grid.cols !== 12) return null;
    if (state.cells.removed?.some((row) => row.some(Boolean))) return null;
    const first = rowLockCells(state, 0);
    for (let r = 1; r < state.grid.rows; r += 1) {
      if (!samePattern(first, rowLockCells(state, r))) return null;
    }
    return first;
  }

  function distance(a, b) {
    let sum = 0;
    for (let i = 0; i < a.length; i += 1) sum += (Number(a[i]) - Number(b[i])) ** 2;
    return Math.sqrt(sum);
  }

  function predictEmpiricalLockCoordinates(state, options = {}) {
    const data = dataset();
    if (!data?.records?.length) return null;
    const cells = options.lockCells || repeatedRowPattern(state);
    if (!cells) return null;
    const stateIndex = Number(options.stateIndex ?? state.grid.lockStateIndex ?? 1) || 1;
    const angles = options.lockAngles || cells.map(() => 30);
    const query = featureVector(cells, angles, stateIndex, data.strandCells || 12);
    const ranked = data.records
      .map((record) => ({
        record,
        distance: distance(record.featureVector || [], query),
      }))
      .sort((a, b) => a.distance - b.distance);
    if (!ranked.length) return null;
    const nearest = ranked[0];
    return {
      schema: "rad-sim.browser-lock-coordinate-prediction.v1",
      method: "nearest",
      lockCells: [...cells],
      lockAngles: [...angles],
      stateIndex,
      failedLockCells: knownFailedLockCells(cells),
      realizedLockCells: realizedLockCells(cells),
      cellCoordinates: nearest.record.cellCoordinates,
      nearestRecord: {
        filename: nearest.record.filename,
        distance: nearest.distance,
        lockCells: nearest.record.lockCells,
        realizedLockCells: nearest.record.realizedLockCells,
        stateIndex: nearest.record.stateIndex,
      },
    };
  }

  function cloneCenter(point) {
    return { x: Number(point.x) || 0, y: Number(point.y) || 0, z: Number(point.z) || 0 };
  }

  function rowMean(row, axis) {
    return row.reduce((sum, point) => sum + point[axis], 0) / Math.max(1, row.length);
  }

  function extent(values) {
    return Math.max(...values) - Math.min(...values);
  }

  function meanCoordinate(coordinates, axisIndex) {
    return coordinates.reduce((sum, point) => sum + (Number(point[axisIndex]) || 0), 0) / Math.max(1, coordinates.length);
  }

  function vectorLength(vector) {
    return Math.hypot(vector[0], vector[1], vector[2]);
  }

  function cross(a, b) {
    return [
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0],
    ];
  }

  function dot(a, b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  function rotateBetweenVectors(point, source, target) {
    const sourceLength = vectorLength(source);
    const targetLength = vectorLength(target);
    if (sourceLength < 1e-12 || targetLength < 1e-12) return point.slice();
    const a = source.map((value) => value / sourceLength);
    const b = target.map((value) => value / targetLength);
    const cosine = Math.max(-1, Math.min(1, dot(a, b)));
    if (Math.abs(cosine - 1) < 1e-12) return point.slice();
    if (Math.abs(cosine + 1) < 1e-12) return point.map((value) => -value);
    const v = cross(a, b);
    const s2 = Math.max(1e-12, dot(v, v));
    const vxPoint = cross(v, point);
    const vxVxPoint = cross(v, vxPoint);
    return point.map((value, index) => value + vxPoint[index] + vxVxPoint[index] * ((1 - cosine) / s2));
  }

  function positionAnchorIndices(positionLocks, count) {
    if (!Array.isArray(positionLocks)) return [];
    const anchors = [];
    for (let c = 0; c < Math.min(count, positionLocks.length); c += 1) {
      if (positionLocks[c] === true) anchors.push(c);
    }
    return anchors.length >= 2 ? [anchors[0], anchors[anchors.length - 1]] : [];
  }

  function cloneBasePoint(point) {
    return { x: Number(point.x) || 0, y: Number(point.y) || 0, z: Number(point.z) || 0 };
  }

  function scaledRowTemplate(prediction, baseRow, positionLocks) {
    const coordinates = prediction.cellCoordinates || [];
    if (coordinates.length !== baseRow.length) return null;
    const anchors = positionAnchorIndices(positionLocks, coordinates.length);
    if (anchors.length === 2) {
      const [first, last] = anchors;
      const measuredFirst = coordinates[first].map((value) => Number(value) || 0);
      const measuredLast = coordinates[last].map((value) => Number(value) || 0);
      const baseFirst = cloneBasePoint(baseRow[first]);
      const baseLast = cloneBasePoint(baseRow[last]);
      const sourceVector = measuredLast.map((value, index) => value - measuredFirst[index]);
      const targetVector = [baseLast.x - baseFirst.x, baseLast.y - baseFirst.y, baseLast.z - baseFirst.z];
      const sourceLength = vectorLength(sourceVector);
      const targetLength = vectorLength(targetVector);
      if (sourceLength > 1e-12 && targetLength > 1e-12) {
        const scale = targetLength / sourceLength;
        const centers = coordinates.map((point) => {
          const relative = point.map((value, index) => (Number(value) || 0) - measuredFirst[index]);
          const rotated = rotateBetweenVectors(relative, sourceVector, targetVector);
          return {
            x: baseFirst.x + rotated[0] * scale,
            y: baseFirst.y + rotated[1] * scale,
            z: baseFirst.z + rotated[2] * scale,
          };
        });
        centers[first] = baseFirst;
        centers[last] = baseLast;
        return {
          centers,
          endpointAnchored: true,
          positionLockedAnchors: anchors,
          scale,
        };
      }
    }
    const measuredX = coordinates.map((point) => Number(point[0]) || 0);
    const baseX = baseRow.map((point) => Number(point.x) || 0);
    const scale = extent(baseX) / Math.max(1e-9, extent(measuredX));
    const meanX = rowMean(baseRow, "x");
    const meanY = rowMean(baseRow, "y");
    const meanZ = rowMean(baseRow, "z");
    const measuredMean = [meanCoordinate(coordinates, 0), meanCoordinate(coordinates, 1), meanCoordinate(coordinates, 2)];
    return {
      centers: coordinates.map((point) => ({
        x: meanX + ((Number(point[0]) || 0) - measuredMean[0]) * scale,
        y: meanY + ((Number(point[1]) || 0) - measuredMean[1]) * scale,
        z: meanZ + ((Number(point[2]) || 0) - measuredMean[2]) * scale,
      })),
      endpointAnchored: false,
      positionLockedAnchors: [],
      scale,
    };
  }

  function recomputeSlope(state, height) {
    const { rows, cols, cellSize } = state.grid;
    const magnitude = RAD.matrix(rows, cols, 0);
    const normal = RAD.matrix(rows, cols, null);
    const tilt = RAD.matrix(rows, cols, 0);
    let total = 0;
    let max = 0;
    let totalTilt = 0;
    let maxTilt = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        const center = height[r][c];
        const left = height[r][Math.max(0, c - 1)];
        const right = height[r][Math.min(cols - 1, c + 1)];
        const down = height[Math.max(0, r - 1)][c];
        const up = height[Math.min(rows - 1, r + 1)][c];
        const dzdx = (right - left) / Math.max(1e-9, cellSize * (c > 0 && c + 1 < cols ? 2 : 1));
        const dzdy = (up - down) / Math.max(1e-9, cellSize * (r > 0 && r + 1 < rows ? 2 : 1));
        const value = Math.hypot(dzdx, dzdy);
        const normalLength = Math.hypot(dzdx, dzdy, 1);
        magnitude[r][c] = value;
        normal[r][c] = { x: -dzdx / normalLength, y: -dzdy / normalLength, z: 1 / normalLength };
        tilt[r][c] = (Math.atan(value) * 180) / Math.PI;
        total += value;
        totalTilt += tilt[r][c];
        max = Math.max(max, value);
        maxTilt = Math.max(maxTilt, tilt[r][c]);
      }
    }
    return {
      magnitude,
      normal,
      tilt,
      mean: total / Math.max(1, rows * cols),
      max,
      meanTilt: totalTilt / Math.max(1, rows * cols),
      maxTilt,
    };
  }

  function computeModelError(state, original, centers, height) {
    const { rows, cols } = state.grid;
    const modelErrorHeight = RAD.matrix(rows, cols, 0);
    const modelErrorCenter = RAD.matrix(rows, cols, 0);
    let heightSq = 0;
    let centerSq = 0;
    let heightMax = 0;
    let centerMax = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        const dh = height[r][c] - original.height[r][c];
        const before = original.centers[r][c];
        const after = centers[r][c];
        const dc = Math.hypot(after.x - before.x, after.y - before.y, after.z - before.z);
        modelErrorHeight[r][c] = dh;
        modelErrorCenter[r][c] = dc;
        heightSq += dh * dh;
        centerSq += dc * dc;
        heightMax = Math.max(heightMax, Math.abs(dh));
        centerMax = Math.max(centerMax, dc);
      }
    }
    const count = Math.max(1, rows * cols);
    return {
      modelErrorHeight,
      modelErrorCenter,
      physicalRmsHeightDelta: Math.sqrt(heightSq / count),
      physicalMaxHeightDelta: heightMax,
      physicalRmsCenterDelta: Math.sqrt(centerSq / count),
      physicalMaxCenterDelta: centerMax,
    };
  }

  function applyEmpiricalLockSurrogate(state, sim) {
    if (state.grid.empiricalLockModel === false) return sim;
    const pattern = repeatedRowPattern(state);
    if (!pattern) return sim;
    const prediction = predictEmpiricalLockCoordinates(state, { lockCells: pattern });
    if (!prediction?.cellCoordinates) return sim;
    const centers = sim.centers.map((row) => row.map(cloneCenter));
    let endpointAnchoredRows = 0;
    let positionLockedAnchorCells = 0;
    let scaleTotal = 0;
    for (let r = 0; r < state.grid.rows; r += 1) {
      const rowTemplate = scaledRowTemplate(prediction, sim.centers[r], state.cells.positionLocked?.[r]);
      if (!rowTemplate) return sim;
      if (rowTemplate.endpointAnchored) endpointAnchoredRows += 1;
      positionLockedAnchorCells += rowTemplate.positionLockedAnchors.length;
      scaleTotal += rowTemplate.scale;
      for (let c = 0; c < state.grid.cols; c += 1) {
        centers[r][c] = {
          x: rowTemplate.centers[c].x,
          y: rowTemplate.centers[c].y,
          z: rowTemplate.centers[c].z,
        };
      }
    }
    const height = RAD.matrix(state.grid.rows, state.grid.cols, (r, c) => centers[r][c].z);
    const slope = recomputeSlope(state, height);
    const modelError = computeModelError(state, sim, centers, height);
    const maxAbsHeight = Math.max(...height.flat().map((value) => Math.abs(value)));
    return {
      ...sim,
      centers,
      height,
      slope,
      modelErrorHeight: modelError.modelErrorHeight,
      modelErrorCenter: modelError.modelErrorCenter,
      empiricalLockPrediction: prediction,
      metrics: {
        ...sim.metrics,
        model: `${sim.metrics?.model || "kinematic"}+lock-data`,
        empiricalLockDatasetApplied: true,
        empiricalLockDatasetSchema: dataset()?.schema || "none",
        empiricalLockNearestFile: prediction.nearestRecord.filename,
        empiricalLockNearestDistance: prediction.nearestRecord.distance,
        empiricalLockRealizedLocks: prediction.realizedLockCells.length,
        empiricalLockEndpointAnchoredRows: endpointAnchoredRows,
        empiricalLockPositionAnchorCells: positionLockedAnchorCells,
        empiricalLockScaleAppliedToMeasuredData: scaleTotal / Math.max(1, state.grid.rows),
        physicalRmsHeightDelta: modelError.physicalRmsHeightDelta,
        physicalMaxHeightDelta: modelError.physicalMaxHeightDelta,
        physicalRmsCenterDelta: modelError.physicalRmsCenterDelta,
        physicalMaxCenterDelta: modelError.physicalMaxCenterDelta,
        maxAbsHeight,
        meanSurfaceSlope: slope.mean,
        maxSurfaceSlope: slope.max,
        meanNormalTilt: slope.meanTilt,
        maxNormalTilt: slope.maxTilt,
      },
    };
  }

  function lockDatasetSummary() {
    const data = dataset();
    if (!data) return { available: false, configurationCount: 0 };
    return {
      available: true,
      schema: data.schema,
      configurationCount: data.configurationCount,
      strandCells: data.strandCells,
      inputDimension: data.inputDimension,
      outputDimension: data.outputDimension,
    };
  }

  RAD.failedLockCellsForNominal = knownFailedLockCells;
  RAD.predictEmpiricalLockCoordinates = predictEmpiricalLockCoordinates;
  RAD.applyEmpiricalLockSurrogate = applyEmpiricalLockSurrogate;
  RAD.lockDatasetSummary = lockDatasetSummary;
})();
