(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  function matrix(rows, cols, fill) {
    return typeof RAD.matrix === "function"
      ? RAD.matrix(rows, cols, fill)
      : Array.from({ length: rows }, (_, r) => Array.from({ length: cols }, (_, c) => (typeof fill === "function" ? fill(r, c) : fill)));
  }

  function cloneCenter(center) {
    return {
      x: Number(center?.x) || 0,
      y: Number(center?.y) || 0,
      z: Number(center?.z) || 0,
    };
  }

  function cellRemoved(state, r, c) {
    return typeof RAD.cellRemoved === "function" ? RAD.cellRemoved(state, r, c) : state.cells?.removed?.[r]?.[c] === true;
  }

  function activeNeighbors(state, r, c) {
    const { rows, cols } = state.grid;
    const out = [];
    if (c + 1 < cols && !cellRemoved(state, r, c + 1)) out.push([r, c + 1]);
    if (r + 1 < rows && !cellRemoved(state, r + 1, c)) out.push([r + 1, c]);
    return out;
  }

  function boundaryState(state) {
    if (typeof RAD.ensureBoundarySchema === "function") return RAD.ensureBoundarySchema(state);
    return state.boundary || {};
  }

  function numericOrNull(value) {
    const number = Number(value);
    return Number.isFinite(number) ? number : null;
  }

  function boundaryConstraintsActive(state) {
    const boundary = boundaryState(state);
    if (boundary.wallType !== "inactive" && boundary.mode === "walls") {
      if (["xMin", "xMax", "yMin", "yMax"].some((key) => numericOrNull(boundary[key]) !== null)) return true;
    }
    if (boundary.wallType !== "inactive" && boundary.mode === "channel" && Number(boundary.channelWidth) > 0) return true;
    const { rows, cols } = state.grid;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (state.cells?.positionLocked?.[r]?.[c] === true && !cellRemoved(state, r, c)) return true;
      }
    }
    return false;
  }

  function referenceBounds(state, centers = null) {
    const { rows, cols } = state.grid;
    let xMin = Infinity;
    let xMax = -Infinity;
    let yMin = Infinity;
    let yMax = -Infinity;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) continue;
        const center = centers?.[r]?.[c] || (typeof RAD.referenceCenter === "function" ? RAD.referenceCenter(state, r, c) : { x: c, y: r, z: 0 });
        xMin = Math.min(xMin, center.x);
        xMax = Math.max(xMax, center.x);
        yMin = Math.min(yMin, center.y);
        yMax = Math.max(yMax, center.y);
      }
    }
    if (!Number.isFinite(xMin)) return { xMin: -1, xMax: 1, yMin: -1, yMax: 1 };
    return { xMin, xMax, yMin, yMax };
  }

  function fitBoundaryToReference(state, padding = 0) {
    const pad = Math.max(0, Number(padding) || 0);
    const bounds = referenceBounds(state);
    const boundary = boundaryState(state);
    boundary.mode = "walls";
    boundary.wallType = boundary.wallType === "inactive" ? "rigid" : boundary.wallType;
    boundary.xMin = bounds.xMin - pad;
    boundary.xMax = bounds.xMax + pad;
    boundary.yMin = bounds.yMin - pad;
    boundary.yMax = bounds.yMax + pad;
    return boundary;
  }

  function clearBoundaryConstraints(state) {
    state.boundary = typeof RAD.defaultBoundaryState === "function" ? RAD.defaultBoundaryState() : { mode: "free", wallType: "rigid" };
    return state.boundary;
  }

  function applyWallConstraint(center, boundary) {
    let x = center.x;
    let y = center.y;
    let contact = false;
    const xMin = numericOrNull(boundary.xMin);
    const xMax = numericOrNull(boundary.xMax);
    const yMin = numericOrNull(boundary.yMin);
    const yMax = numericOrNull(boundary.yMax);
    if (xMin !== null && x < xMin) {
      x = xMin;
      contact = true;
    }
    if (xMax !== null && x > xMax) {
      x = xMax;
      contact = true;
    }
    if (yMin !== null && y < yMin) {
      y = yMin;
      contact = true;
    }
    if (yMax !== null && y > yMax) {
      y = yMax;
      contact = true;
    }
    return { x, y, contact };
  }

  function applyChannelConstraint(center, boundary) {
    const width = Math.max(0, Number(boundary.channelWidth) || 0);
    if (width <= 0) return { x: center.x, y: center.y, contact: false };
    const half = width / 2;
    const axis = boundary.channelAxis === "y" ? "y" : "x";
    const freeAxis = axis === "x" ? "x" : "y";
    const blockedAxis = axis === "x" ? "y" : "x";
    const out = { x: center.x, y: center.y, contact: false };
    out[freeAxis] = center[freeAxis];
    if (center[blockedAxis] < -half) {
      out[blockedAxis] = -half;
      out.contact = true;
    } else if (center[blockedAxis] > half) {
      out[blockedAxis] = half;
      out.contact = true;
    }
    return out;
  }

  function realizedCenter(state, boundary, preferred, r, c) {
    let next = cloneCenter(preferred);
    let wallContact = false;
    if (state.cells?.positionLocked?.[r]?.[c] === true) {
      const target = typeof RAD.positionLockTarget === "function" ? RAD.positionLockTarget(state, r, c) : preferred;
      next.x = target.x;
      next.y = target.y;
      return { center: next, wallContact: false, positionLocked: true };
    }
    if (boundary.wallType !== "inactive" && boundary.mode === "walls") {
      const constrained = applyWallConstraint(next, boundary);
      wallContact = constrained.contact;
      next.x = constrained.x;
      next.y = constrained.y;
    } else if (boundary.wallType !== "inactive" && boundary.mode === "channel") {
      const constrained = applyChannelConstraint(next, boundary);
      wallContact = constrained.contact;
      next.x = constrained.x;
      next.y = constrained.y;
    }
    if (boundary.wallType === "soft" && wallContact) {
      const stiffness = 0.55;
      next.x = preferred.x + (next.x - preferred.x) * stiffness;
      next.y = preferred.y + (next.y - preferred.y) * stiffness;
    }
    return { center: next, wallContact, positionLocked: false };
  }

  function solveBoundaryConstraints(state, alpha, preferredCenters = null) {
    const { rows, cols } = state.grid;
    const boundary = boundaryState(state);
    const preferred = preferredCenters || (typeof RAD.computeCenters === "function" ? RAD.computeCenters(state, alpha) : matrix(rows, cols, (r, c) => ({ x: c, y: r, z: 0 })));
    const centers = matrix(rows, cols, (r, c) => cloneCenter(preferred[r][c]));
    const compressionResidual = matrix(rows, cols, 0);
    const inducedHeight = matrix(rows, cols, 0);
    const constraintDisplacement = matrix(rows, cols, 0);
    const wallContact = matrix(rows, cols, false);
    const positionLockContact = matrix(rows, cols, false);
    const edgeCounts = matrix(rows, cols, 0);
    let constrainedCellCount = 0;

    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) continue;
        const realized = realizedCenter(state, boundary, preferred[r][c], r, c);
        centers[r][c].x = realized.center.x;
        centers[r][c].y = realized.center.y;
        wallContact[r][c] = realized.wallContact;
        positionLockContact[r][c] = realized.positionLocked;
        constraintDisplacement[r][c] = Math.hypot(centers[r][c].x - preferred[r][c].x, centers[r][c].y - preferred[r][c].y);
        if (constraintDisplacement[r][c] > 1e-9) constrainedCellCount += 1;
      }
    }

    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) continue;
        for (const [nr, nc] of activeNeighbors(state, r, c)) {
          const preferredLength = Math.hypot(preferred[nr][nc].x - preferred[r][c].x, preferred[nr][nc].y - preferred[r][c].y);
          const realizedLength = Math.hypot(centers[nr][nc].x - centers[r][c].x, centers[nr][nc].y - centers[r][c].y);
          const compression = Math.max(0, preferredLength - realizedLength);
          compressionResidual[r][c] += compression;
          compressionResidual[nr][nc] += compression;
          edgeCounts[r][c] += 1;
          edgeCounts[nr][nc] += 1;
        }
      }
    }

    let maxCompressionResidual = 0;
    let meanCompressionResidual = 0;
    let maxInducedHeight = 0;
    let activeCells = 0;
    const zGain = Math.max(0, Number(boundary.zGain) || 0);
    const threshold = Math.max(0, Number(boundary.zThreshold) || 0);
    const power = Math.max(0.5, Number(boundary.zPower) || 1);
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (cellRemoved(state, r, c)) continue;
        activeCells += 1;
        const edgeResidual = compressionResidual[r][c] / Math.max(1, edgeCounts[r][c]);
        const residual = Math.max(edgeResidual, constraintDisplacement[r][c]);
        compressionResidual[r][c] = residual;
        const induced = zGain * Math.pow(Math.max(0, residual - threshold), power);
        inducedHeight[r][c] = state.cells?.positionLocked?.[r]?.[c] === true ? 0 : induced;
        maxCompressionResidual = Math.max(maxCompressionResidual, residual);
        meanCompressionResidual += residual;
        maxInducedHeight = Math.max(maxInducedHeight, Math.abs(inducedHeight[r][c]));
      }
    }

    return {
      schema: "rad-sim.constraint-realization.v1",
      boundary: { ...boundary },
      preferredCenters: preferred,
      centers,
      compressionResidual,
      inducedHeight,
      constraintDisplacement,
      wallContact,
      positionLockContact,
      active: boundaryConstraintsActive(state),
      metrics: {
        constrainedCellCount,
        meanCompressionResidual: meanCompressionResidual / Math.max(1, activeCells),
        maxCompressionResidual,
        maxInducedHeight,
      },
    };
  }

  RAD.boundaryConstraintsActive = boundaryConstraintsActive;
  RAD.referenceBounds = referenceBounds;
  RAD.fitBoundaryToReference = fitBoundaryToReference;
  RAD.clearBoundaryConstraints = clearBoundaryConstraints;
  RAD.solveBoundaryConstraints = solveBoundaryConstraints;
})();
