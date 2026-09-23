(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  function clonePoint(point) {
    return { x: Number(point.x) || 0, y: Number(point.y) || 0, z: Number(point.z) || 0 };
  }

  function cellRemoved(state, r, c) {
    return state.cells?.removed?.[r]?.[c] === true;
  }

  function positionLocked(state, r, c) {
    return state.cells?.positionLocked?.[r]?.[c] === true && !cellRemoved(state, r, c);
  }

  function positionLockTarget(state, r, c) {
    if (typeof RAD.positionLockTarget === "function") return RAD.positionLockTarget(state, r, c);
    const reference =
      typeof RAD.referenceCenter === "function"
        ? RAD.referenceCenter(state, r, c)
        : { x: 0, y: 0, z: 0 };
    return clonePoint(reference);
  }

  function neighborCells(state, r, c) {
    const { rows, cols } = state.grid;
    const out = [];
    if (r > 0 && !cellRemoved(state, r - 1, c)) out.push([r - 1, c]);
    if (r + 1 < rows && !cellRemoved(state, r + 1, c)) out.push([r + 1, c]);
    if (c > 0 && !cellRemoved(state, r, c - 1)) out.push([r, c - 1]);
    if (c + 1 < cols && !cellRemoved(state, r, c + 1)) out.push([r, c + 1]);
    return out;
  }

  function recomputeLinkStrain(state, centers) {
    const { rows, cols, cellSize } = state.grid;
    const horizontal = RAD.matrix(rows, Math.max(0, cols - 1), 0);
    const vertical = RAD.matrix(Math.max(0, rows - 1), cols, 0);
    let totalAbs = 0;
    let maxAbs = 0;
    let count = 0;
    let skipped = 0;
    function add(a, b) {
      const strain = Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z) / cellSize - 1;
      totalAbs += Math.abs(strain);
      maxAbs = Math.max(maxAbs, Math.abs(strain));
      count += 1;
      return strain;
    }
    function activeLink(r0, c0, r1, c1) {
      if (cellRemoved(state, r0, c0) || cellRemoved(state, r1, c1)) {
        skipped += 1;
        return false;
      }
      return true;
    }
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c + 1 < cols; c += 1) {
        horizontal[r][c] = activeLink(r, c, r, c + 1) ? add(centers[r][c], centers[r][c + 1]) : null;
      }
    }
    for (let r = 0; r + 1 < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        vertical[r][c] = activeLink(r, c, r + 1, c) ? add(centers[r][c], centers[r + 1][c]) : null;
      }
    }
    return { horizontal, vertical, meanAbs: totalAbs / Math.max(1, count), maxAbs, activeEdges: count, skippedEdges: skipped };
  }

  function recomputeSlope(state, height) {
    const { rows, cols, cellSize } = state.grid;
    const gradient = RAD.matrix(rows, cols, 0);
    const normalTilt = RAD.matrix(rows, cols, 0);
    const normal = RAD.matrix(rows, cols, null);
    let total = 0;
    let max = 0;
    let totalTilt = 0;
    let maxTilt = 0;
    let activeCount = 0;
    function sampleHeight(rr, cc, fallback) {
      return cellRemoved(state, rr, cc) ? fallback : height[rr][cc];
    }
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) {
          gradient[r][c] = null;
          normalTilt[r][c] = null;
          continue;
        }
        const center = height[r][c];
        const left = sampleHeight(r, Math.max(0, c - 1), center);
        const right = sampleHeight(r, Math.min(cols - 1, c + 1), center);
        const down = sampleHeight(Math.max(0, r - 1), c, center);
        const up = sampleHeight(Math.min(rows - 1, r + 1), c, center);
        const dzdx = (right - left) / (cellSize * (c > 0 && c + 1 < cols ? 2 : 1));
        const dzdy = (up - down) / (cellSize * (r > 0 && r + 1 < rows ? 2 : 1));
        const slope = Math.hypot(dzdx, dzdy);
        const tilt = (Math.atan(slope) * 180) / Math.PI;
        const normalLength = Math.hypot(dzdx, dzdy, 1);
        gradient[r][c] = slope;
        normalTilt[r][c] = tilt;
        normal[r][c] = { x: -dzdx / normalLength, y: -dzdy / normalLength, z: 1 / normalLength };
        total += slope;
        max = Math.max(max, slope);
        totalTilt += tilt;
        maxTilt = Math.max(maxTilt, tilt);
        activeCount += 1;
      }
    }
    return {
      gradient,
      magnitude: gradient,
      normal,
      normalTilt,
      tilt: normalTilt,
      mean: total / Math.max(1, activeCount),
      max,
      meanTilt: totalTilt / Math.max(1, activeCount),
      maxTilt,
      activeCells: activeCount,
    };
  }

  function computeModelError(state, base, centers, height) {
    const rows = base.height.length;
    const cols = base.height[0]?.length || 0;
    const modelErrorHeight = RAD.matrix(rows, cols, 0);
    const modelErrorCenter = RAD.matrix(rows, cols, 0);
    let heightSquared = 0;
    let centerSquared = 0;
    let heightMax = 0;
    let centerMax = 0;
    let count = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) continue;
        const heightDelta = height[r][c] - base.height[r][c];
        const baseCenter = base.centers[r][c];
        const center = centers[r][c];
        const centerDelta = Math.hypot(center.x - baseCenter.x, center.y - baseCenter.y, center.z - baseCenter.z);
        modelErrorHeight[r][c] = heightDelta;
        modelErrorCenter[r][c] = centerDelta;
        heightSquared += heightDelta * heightDelta;
        centerSquared += centerDelta * centerDelta;
        heightMax = Math.max(heightMax, Math.abs(heightDelta));
        centerMax = Math.max(centerMax, centerDelta);
        count += 1;
      }
    }
    return {
      modelErrorHeight,
      modelErrorCenter,
      physicalRmsHeightDelta: Math.sqrt(heightSquared / Math.max(1, count)),
      physicalMaxHeightDelta: heightMax,
      physicalRmsCenterDelta: Math.sqrt(centerSquared / Math.max(1, count)),
      physicalMaxCenterDelta: centerMax,
    };
  }

  function simulatePhysicalRelaxation(state, options = {}) {
    const base = options.baseSim || RAD.simulate(state);
    const { rows, cols } = state.grid;
    const iterations = Math.max(1, Math.min(80, Math.floor(options.iterations ?? 18)));
    const springGain = Math.max(0, Math.min(0.45, Number(options.springGain ?? 0.16)));
    const anchorGain = Math.max(0, Math.min(0.95, Number(options.anchorGain ?? 0.32)));
    let centers = base.centers.map((row) => row.map(clonePoint));
    const target = base.centers.map((row) => row.map(clonePoint));
    const fixedTargets = RAD.matrix(rows, cols, null);
    let positionLockedCells = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (!positionLocked(state, r, c)) continue;
        const fixture = positionLockTarget(state, r, c);
        fixedTargets[r][c] = fixture;
        centers[r][c] = clonePoint(fixture);
        target[r][c] = clonePoint(fixture);
        positionLockedCells += 1;
      }
    }

    for (let step = 0; step < iterations; step += 1) {
      const next = centers.map((row) => row.map(clonePoint));
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          if (cellRemoved(state, r, c)) continue;
          if (fixedTargets[r][c]) {
            next[r][c] = clonePoint(fixedTargets[r][c]);
            continue;
          }
          const neighbors = neighborCells(state, r, c);
          if (!neighbors.length) continue;
          const current = centers[r][c];
          const average = neighbors.reduce(
            (sum, [nr, nc]) => {
              sum.x += centers[nr][nc].x;
              sum.y += centers[nr][nc].y;
              sum.z += centers[nr][nc].z;
              return sum;
            },
            { x: 0, y: 0, z: 0 }
          );
          average.x /= neighbors.length;
          average.y /= neighbors.length;
          average.z /= neighbors.length;
          const commanded = Math.abs(state.cells.commandAlpha[r][c]) + Math.abs(state.cells.commandZ[r][c]);
          const locked = state.cells.locked[r][c];
          const localAnchor = Math.min(0.98, anchorGain + (locked ? 0.55 : 0) + Math.min(0.35, commanded * 0.35));
          next[r][c] = {
            x: current.x + springGain * (average.x - current.x) + localAnchor * (target[r][c].x - current.x),
            y: current.y + springGain * (average.y - current.y) + localAnchor * (target[r][c].y - current.y),
            z: current.z + springGain * (average.z - current.z) + localAnchor * (target[r][c].z - current.z),
          };
        }
      }
      centers = next;
    }

    const height = RAD.matrix(rows, cols, (r, c) => centers[r][c].z);
    const slope = recomputeSlope(state, height);
    const linkStrain = recomputeLinkStrain(state, centers);
    const modelError = computeModelError(state, base, centers, height);
    let maxAbsHeight = 0;
    let activeCells = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) continue;
        activeCells += 1;
        maxAbsHeight = Math.max(maxAbsHeight, Math.abs(height[r][c]));
      }
    }
    return {
      ...base,
      height,
      centers,
      modelErrorHeight: modelError.modelErrorHeight,
      modelErrorCenter: modelError.modelErrorCenter,
      slope,
      linkStrain,
      metrics: {
        ...base.metrics,
        model: "spring-preview",
        physicalPreview: true,
        physicalIterations: iterations,
        physicalRmsHeightDelta: modelError.physicalRmsHeightDelta,
        physicalMaxHeightDelta: modelError.physicalMaxHeightDelta,
        physicalRmsCenterDelta: modelError.physicalRmsCenterDelta,
        physicalMaxCenterDelta: modelError.physicalMaxCenterDelta,
        maxAbsHeight,
        meanAbsLinkStrain: linkStrain.meanAbs,
        maxAbsLinkStrain: linkStrain.maxAbs,
        physicalActiveCells: activeCells,
        physicalPositionLockedCells: positionLockedCells,
        physicalActiveSpringEdges: linkStrain.activeEdges,
        physicalSkippedSpringEdges: linkStrain.skippedEdges,
        meanSurfaceSlope: slope.mean,
        maxSurfaceSlope: slope.max,
        meanNormalTilt: slope.meanTilt,
        maxNormalTilt: slope.maxTilt,
      },
    };
  }

  function simulateActive(state) {
    const base = RAD.simulate(state);
    let active;
    if (state.view?.simulationMode === "springPreview") {
      active = simulatePhysicalRelaxation(state, { baseSim: base });
    } else if (state.view?.simulationMode === "constraintSolved") {
      active = { ...base, metrics: { ...base.metrics, model: "constraint-solved", physicalPreview: false, constraintSolved: true } };
    } else {
      active = { ...base, metrics: { ...base.metrics, model: "kinematic", physicalPreview: false } };
    }
    return typeof RAD.applyEmpiricalLockSurrogate === "function"
      ? RAD.applyEmpiricalLockSurrogate(state, active)
      : active;
  }

  RAD.simulatePhysicalRelaxation = simulatePhysicalRelaxation;
  RAD.simulateActive = simulateActive;
})();
