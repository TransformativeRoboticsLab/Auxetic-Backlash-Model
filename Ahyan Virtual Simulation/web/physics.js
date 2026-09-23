(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  function clonePoint(point) {
    return { x: Number(point.x) || 0, y: Number(point.y) || 0, z: Number(point.z) || 0 };
  }

  function neighborCells(rows, cols, r, c) {
    const out = [];
    if (r > 0) out.push([r - 1, c]);
    if (r + 1 < rows) out.push([r + 1, c]);
    if (c > 0) out.push([r, c - 1]);
    if (c + 1 < cols) out.push([r, c + 1]);
    return out;
  }

  function recomputeLinkStrain(state, centers) {
    const { rows, cols, cellSize } = state.grid;
    const horizontal = RAD.matrix(rows, Math.max(0, cols - 1), 0);
    const vertical = RAD.matrix(Math.max(0, rows - 1), cols, 0);
    let totalAbs = 0;
    let maxAbs = 0;
    let count = 0;
    function add(a, b) {
      const strain = Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z) / cellSize - 1;
      totalAbs += Math.abs(strain);
      maxAbs = Math.max(maxAbs, Math.abs(strain));
      count += 1;
      return strain;
    }
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c + 1 < cols; c += 1) horizontal[r][c] = add(centers[r][c], centers[r][c + 1]);
    }
    for (let r = 0; r + 1 < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) vertical[r][c] = add(centers[r][c], centers[r + 1][c]);
    }
    return { horizontal, vertical, meanAbs: totalAbs / Math.max(1, count), maxAbs };
  }

  function recomputeSlope(state, height) {
    const { rows, cols, cellSize } = state.grid;
    const gradient = RAD.matrix(rows, cols, 0);
    const normalTilt = RAD.matrix(rows, cols, 0);
    let total = 0;
    let max = 0;
    let totalTilt = 0;
    let maxTilt = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        const left = height[r][Math.max(0, c - 1)];
        const right = height[r][Math.min(cols - 1, c + 1)];
        const down = height[Math.max(0, r - 1)][c];
        const up = height[Math.min(rows - 1, r + 1)][c];
        const dzdx = (right - left) / (cellSize * (c > 0 && c + 1 < cols ? 2 : 1));
        const dzdy = (up - down) / (cellSize * (r > 0 && r + 1 < rows ? 2 : 1));
        const slope = Math.hypot(dzdx, dzdy);
        const tilt = (Math.atan(slope) * 180) / Math.PI;
        gradient[r][c] = slope;
        normalTilt[r][c] = tilt;
        total += slope;
        max = Math.max(max, slope);
        totalTilt += tilt;
        maxTilt = Math.max(maxTilt, tilt);
      }
    }
    return {
      gradient,
      normalTilt,
      mean: total / Math.max(1, rows * cols),
      max,
      meanTilt: totalTilt / Math.max(1, rows * cols),
      maxTilt,
    };
  }

  function computeModelError(base, centers, height) {
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

    for (let step = 0; step < iterations; step += 1) {
      const next = centers.map((row) => row.map(clonePoint));
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          const neighbors = neighborCells(rows, cols, r, c);
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
          if (locked) next[r][c].z = 0;
        }
      }
      centers = next;
    }

    const height = RAD.matrix(rows, cols, (r, c) => centers[r][c].z);
    const slope = recomputeSlope(state, height);
    const linkStrain = recomputeLinkStrain(state, centers);
    const modelError = computeModelError(base, centers, height);
    let maxAbsHeight = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
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
        meanSurfaceSlope: slope.mean,
        maxSurfaceSlope: slope.max,
        meanNormalTilt: slope.meanTilt,
        maxNormalTilt: slope.maxTilt,
      },
    };
  }

  function simulateActive(state) {
    const base = RAD.simulate(state);
    if (state.view?.simulationMode === "springPreview") {
      return simulatePhysicalRelaxation(state, { baseSim: base });
    }
    return { ...base, metrics: { ...base.metrics, model: "kinematic", physicalPreview: false } };
  }

  RAD.simulatePhysicalRelaxation = simulatePhysicalRelaxation;
  RAD.simulateActive = simulateActive;
})();
