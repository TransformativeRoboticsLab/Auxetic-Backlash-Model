(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
  }

  function cloneState(state) {
    return JSON.parse(JSON.stringify(state));
  }

  function cloneData(value) {
    return value === undefined ? null : JSON.parse(JSON.stringify(value));
  }

  function csvValue(value) {
    const text = String(value ?? "");
    return /[",\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
  }

  function planningStateFrom(state) {
    const next = cloneState(state);
    const locked = JSON.parse(JSON.stringify(next.cells.locked));
    const actuatorAllowed = JSON.parse(JSON.stringify(next.cells.actuatorAllowed || RAD.matrix(next.grid.rows, next.grid.cols, true)));
    RAD.clearCommands(next);
    next.cells.locked = locked;
    next.cells.actuatorAllowed = actuatorAllowed;
    return next;
  }

  function objective(sim, state) {
    const actuatorPenalty = state.target.actuatorPenalty ?? 0.012;
    const travelPenalty = state.target.travelPenalty ?? 0.018;
    const n = state.grid.rows * state.grid.cols;
    return sim.metrics.rmsTargetError + actuatorPenalty * (sim.metrics.recommendedActuators / n) + travelPenalty * sim.metrics.meanTravel;
  }

  function rmsDelta(a, b) {
    let sum = 0;
    let count = 0;
    let peak = 0;
    for (let r = 0; r < a.length; r += 1) {
      for (let c = 0; c < a[r].length; c += 1) {
        const delta = a[r][c] - b[r][c];
        sum += delta * delta;
        peak = Math.max(peak, Math.abs(delta));
        count += 1;
      }
    }
    return { rms: Math.sqrt(sum / Math.max(1, count)), peak };
  }

  function rmsTargetResidual(target, height) {
    let sum = 0;
    let count = 0;
    let maxAbs = 0;
    let bias = 0;
    for (let r = 0; r < target.length; r += 1) {
      for (let c = 0; c < target[r].length; c += 1) {
        const residual = target[r][c] - height[r][c];
        sum += residual * residual;
        bias += residual;
        maxAbs = Math.max(maxAbs, Math.abs(residual));
        count += 1;
      }
    }
    return {
      rms: Math.sqrt(sum / Math.max(1, count)),
      maxAbs,
      bias: bias / Math.max(1, count),
    };
  }

  function countResponsiveCells(a, b, threshold) {
    let count = 0;
    for (let r = 0; r < a.length; r += 1) {
      for (let c = 0; c < a[r].length; c += 1) {
        if (Math.abs(a[r][c] - b[r][c]) >= threshold) count += 1;
      }
    }
    return count;
  }

  function flattenDelta(next, base) {
    const values = [];
    let normSquared = 0;
    let peak = 0;
    for (let r = 0; r < next.length; r += 1) {
      for (let c = 0; c < next[r].length; c += 1) {
        const value = next[r][c] - base[r][c];
        values.push(value);
        normSquared += value * value;
        peak = Math.max(peak, Math.abs(value));
      }
    }
    return { values, norm: Math.sqrt(normSquared), peak };
  }

  function dot(a, b) {
    let value = 0;
    for (let i = 0; i < Math.min(a.length, b.length); i += 1) value += a[i] * b[i];
    return value;
  }

  function targetResidualVector(sim) {
    const values = [];
    for (let r = 0; r < sim.target.length; r += 1) {
      for (let c = 0; c < sim.target[r].length; c += 1) values.push(sim.target[r][c] - sim.height[r][c]);
    }
    return values;
  }

  function profileUpdateMap(profile) {
    const updates = Array.isArray(profile?.recommendedUpdates) ? profile.recommendedUpdates : [];
    const mapped = new Map();
    for (const update of updates) {
      if (update && typeof update.name === "string") mapped.set(update.name, update);
    }
    return mapped;
  }

  function profileUpdateValue(updates, name) {
    const value = Number(updates.get(name)?.proposed);
    return Number.isFinite(value) ? value : null;
  }

  function profileScale(value, uncertainty, tolerance) {
    const values = [Math.abs(Number(tolerance) || 0), 1e-12];
    if (Number.isFinite(Number(value))) values.push(Math.abs(Number(value)));
    if (Number.isFinite(Number(uncertainty))) values.push(Math.abs(Number(uncertainty)));
    return Math.max(...values);
  }

  function residualStats(target, actual) {
    if (!Array.isArray(target) || !Array.isArray(actual)) return { rms: 0, maxAbs: 0, field: null };
    let sum = 0;
    let count = 0;
    let maxAbs = 0;
    const field = RAD.matrix(target.length, target[0]?.length || 0, 0);
    for (let r = 0; r < target.length; r += 1) {
      for (let c = 0; c < target[r].length; c += 1) {
        const residual = (Number(target[r][c]) || 0) - (Number(actual?.[r]?.[c]) || 0);
        field[r][c] = residual;
        sum += residual ** 2;
        maxAbs = Math.max(maxAbs, Math.abs(residual));
        count += 1;
      }
    }
    return { rms: Math.sqrt(sum / Math.max(1, count)), maxAbs, field };
  }

  function targetCellCount(target, baseline, tolerance) {
    if (!Array.isArray(target) || !Array.isArray(baseline)) return 0;
    let count = 0;
    for (let r = 0; r < target.length; r += 1) {
      for (let c = 0; c < target[r].length; c += 1) {
        if (Math.abs((Number(target[r][c]) || 0) - (Number(baseline?.[r]?.[c]) || 0)) > tolerance) count += 1;
      }
    }
    return count;
  }

  function bandFailureCount(field, band, tolerance) {
    if (!Array.isArray(field)) return 0;
    let count = 0;
    for (let r = 0; r < field.length; r += 1) {
      for (let c = 0; c < field[r].length; c += 1) {
        if (Math.abs(Number(field[r][c]) || 0) > band + tolerance) count += 1;
      }
    }
    return count;
  }

  function targetReachabilityReport(state, baseSim, columns, options = {}) {
    const { rows, cols } = state.grid;
    const responseThreshold = options.responseThreshold ?? 0.012;
    const targetThreshold = options.targetThreshold ?? responseThreshold;
    const heightReachableMap = RAD.matrix(rows, cols, false);
    const positiveHeightReachableMap = RAD.matrix(rows, cols, false);
    const negativeHeightReachableMap = RAD.matrix(rows, cols, false);
    const topology = baseSim.topology || (typeof RAD.topologyDiagnostics === "function" ? RAD.topologyDiagnostics(state) : null);
    const topologyReachableMap = RAD.matrix(rows, cols, false);
    const topologyBlockedHeightMap = RAD.matrix(rows, cols, 0);
    const activeComponents = new Set();
    if (topology?.componentLabels) {
      for (const column of columns || []) {
        const label = topology.componentLabels?.[column.r]?.[column.c];
        if (Number.isInteger(label) && label >= 0) activeComponents.add(label);
      }
      for (let rr = 0; rr < rows; rr += 1) {
        for (let cc = 0; cc < cols; cc += 1) {
          topologyReachableMap[rr][cc] = activeComponents.has(topology.componentLabels[rr][cc]);
        }
      }
    }
    const targetHeightMap = RAD.matrix(rows, cols, 0);
    const underactuatedHeightMap = RAD.matrix(rows, cols, 0);
    let heightReachableCells = 0;
    let positiveHeightReachableCells = 0;
    let negativeHeightReachableCells = 0;
    let targetHeightCells = 0;
    let upwardTargetHeightCells = 0;
    let downwardTargetHeightCells = 0;
    let underactuatedHeightCells = 0;
    let positiveUnderactuatedHeightCells = 0;
    let negativeUnderactuatedHeightCells = 0;
    let topologyReachableCells = 0;
    let topologyBlockedHeightCells = 0;
    let unreachableSquared = 0;
    let topologyBlockedSquared = 0;
    let maxUnreachableHeightResidual = 0;
    let maxTopologyBlockedHeightResidual = 0;
    let worstUnderactuatedCell = null;
    let worstTopologyBlockedCell = null;

    for (let index = 0; index < rows * cols; index += 1) {
      const r = Math.floor(index / cols);
      const c = index % cols;
      if (topologyReachableMap[r][c]) topologyReachableCells += 1;
      const positiveReachable = columns.some((column) => (column.heightDelta?.[index] || 0) >= responseThreshold);
      const negativeReachable = columns.some((column) => (column.heightDelta?.[index] || 0) <= -responseThreshold);
      const reachable = positiveReachable || negativeReachable;
      heightReachableMap[r][c] = reachable;
      positiveHeightReachableMap[r][c] = positiveReachable;
      negativeHeightReachableMap[r][c] = negativeReachable;
      if (reachable) heightReachableCells += 1;
      if (positiveReachable) positiveHeightReachableCells += 1;
      if (negativeReachable) negativeHeightReachableCells += 1;
      const residual = (baseSim.target?.[r]?.[c] || 0) - (baseSim.height?.[r]?.[c] || 0);
      const upwardRequested = residual >= targetThreshold;
      const downwardRequested = residual <= -targetThreshold;
      const requested = upwardRequested || downwardRequested;
      if (requested) {
        targetHeightCells += 1;
        targetHeightMap[r][c] = residual;
      }
      if (upwardRequested) upwardTargetHeightCells += 1;
      if (downwardRequested) downwardTargetHeightCells += 1;
      const signedUnderactuated = (upwardRequested && !positiveReachable) || (downwardRequested && !negativeReachable);
      if (signedUnderactuated) {
        const magnitude = Math.abs(residual);
        underactuatedHeightCells += 1;
        if (upwardRequested) positiveUnderactuatedHeightCells += 1;
        if (downwardRequested) negativeUnderactuatedHeightCells += 1;
        underactuatedHeightMap[r][c] = magnitude;
        unreachableSquared += residual * residual;
        if (magnitude > maxUnreachableHeightResidual) {
          maxUnreachableHeightResidual = magnitude;
          worstUnderactuatedCell = { row: r, col: c, residual, absResidual: magnitude };
        }
      }
      if (requested && !topologyReachableMap[r][c]) {
        const magnitude = Math.abs(residual);
        topologyBlockedHeightCells += 1;
        topologyBlockedHeightMap[r][c] = magnitude;
        topologyBlockedSquared += residual * residual;
        if (magnitude > maxTopologyBlockedHeightResidual) {
          maxTopologyBlockedHeightResidual = magnitude;
          worstTopologyBlockedCell = { row: r, col: c, residual, absResidual: magnitude };
        }
      }
    }

    return {
      model: "finite-response-height-reachability",
      directionalModel: "sign-compatible finite-response-height-reachability",
      responseThreshold,
      targetThreshold,
      heightReachableMap,
      positiveHeightReachableMap,
      negativeHeightReachableMap,
      topologyReachableMap,
      targetHeightMap,
      underactuatedHeightMap,
      topologyBlockedHeightMap,
      heightReachableCells,
      positiveHeightReachableCells,
      negativeHeightReachableCells,
      topologyReachableCells,
      topologyComponentCount: topology?.componentCount || 0,
      topologyActuatedComponentCount: activeComponents.size,
      targetHeightCells,
      upwardTargetHeightCells,
      downwardTargetHeightCells,
      underactuatedHeightCells,
      positiveUnderactuatedHeightCells,
      negativeUnderactuatedHeightCells,
      topologyBlockedHeightCells,
      unreachableHeightRms: Math.sqrt(unreachableSquared / Math.max(1, underactuatedHeightCells)),
      topologyBlockedHeightRms: Math.sqrt(topologyBlockedSquared / Math.max(1, topologyBlockedHeightCells)),
      maxUnreachableHeightResidual,
      maxTopologyBlockedHeightResidual,
      worstUnderactuatedCell,
      worstTopologyBlockedCell,
    };
  }

  function commandFromResidual(state, residual) {
    return {
      commandZ: RAD.clampCommandZ(state, residual * 0.82),
      commandAlpha: RAD.clampCommandAlpha(state, -Math.sign(residual || 1) * Math.min(0.42, Math.abs(residual) * 0.24)),
    };
  }

  function buildInverseDesignPlan(state, options = {}) {
    const projectedState = planningStateFrom(state);
    const baseSim = RAD.simulate(projectedState);
    const baseScore = objective(baseSim, projectedState);
    const { rows, cols } = state.grid;
    const maxActuators = options.maxActuators ?? state.inverse?.maxActuators ?? 24;
    const maxEvaluations = Math.max(1, Math.min(maxActuators, rows * cols));
    const commands = [];
    const history = [];
    const selected = new Set();
    let currentScore = baseScore;
    let currentSim = baseSim;

    for (let step = 0; step < maxEvaluations; step += 1) {
      let best = null;
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          const key = `${r},${c}`;
          if (projectedState.cells.removed?.[r]?.[c] || projectedState.cells.locked[r][c] || projectedState.cells.actuatorAllowed?.[r]?.[c] === false || selected.has(key)) continue;
          const residual = currentSim.target[r][c] - currentSim.height[r][c];
          const command = commandFromResidual(projectedState, residual);
          if (Math.abs(command.commandZ) + Math.abs(command.commandAlpha) < 1e-9) continue;
          const trialState = cloneState(projectedState);
          trialState.cells.commandZ[r][c] = command.commandZ;
          trialState.cells.commandAlpha[r][c] = command.commandAlpha;
          const trialSim = RAD.simulate(trialState);
          const trialScore = objective(trialSim, trialState);
          const improvement = currentScore - trialScore;
          if (!best || improvement > best.improvement) {
            best = { r, c, residual, command, trialScore, improvement, trialSim };
          }
        }
      }
      if (!best || best.improvement <= 1e-6) break;
      const key = `${best.r},${best.c}`;
      selected.add(key);
      projectedState.cells.commandZ[best.r][best.c] = best.command.commandZ;
      projectedState.cells.commandAlpha[best.r][best.c] = best.command.commandAlpha;
      currentScore = best.trialScore;
      currentSim = best.trialSim;
      const accepted = {
        step: step + 1,
        r: best.r,
        c: best.c,
        commandZ: best.command.commandZ,
        commandAlpha: best.command.commandAlpha,
        improvement: best.improvement,
        residual: best.residual,
        score: best.trialScore,
      };
      commands.push(accepted);
      history.push({
        step: accepted.step,
        score: accepted.score,
        improvement: accepted.improvement,
        rmsError: best.trialSim.metrics.rmsTargetError,
        actuators: best.trialSim.metrics.recommendedActuators,
      });
    }

    const candidates = [];

    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        const key = `${r},${c}`;
        if (projectedState.cells.removed?.[r]?.[c] || projectedState.cells.locked[r][c] || projectedState.cells.actuatorAllowed?.[r]?.[c] === false || selected.has(key)) continue;
        const residual = currentSim.target[r][c] - currentSim.height[r][c];
        const command = commandFromResidual(projectedState, residual);
        const trialState = cloneState(projectedState);
        trialState.cells.commandZ[r][c] = command.commandZ;
        trialState.cells.commandAlpha[r][c] = command.commandAlpha;
        const trialSim = RAD.simulate(trialState);
        const trialScore = objective(trialSim, trialState);
        const improvement = currentScore - trialScore;
        const reach = Number.isFinite(currentSim.dieOff[r][c]) ? currentSim.dieOff[r][c] : rows + cols;
        candidates.push({
          r,
          c,
          residual,
          improvement,
          score: trialScore,
          reach,
          commandZ: trialState.cells.commandZ[r][c],
          commandAlpha: trialState.cells.commandAlpha[r][c],
        });
      }
    }

    candidates.sort((a, b) => b.improvement - a.improvement || Math.abs(b.residual) - Math.abs(a.residual));
    const projectedSim = RAD.simulate(projectedState);
    const projectedScore = objective(projectedSim, projectedState);
    const plan = {
      baseScore,
      projectedScore,
      projectedError: projectedSim.metrics.rmsTargetError,
      projectedActuators: projectedSim.metrics.recommendedActuators,
      totalImprovement: baseScore - projectedScore,
      maxActuators,
      candidateCount: candidates.length,
      strategy: "iterative-greedy",
      steps: history.length,
      history,
      candidates: candidates.slice(0, Math.min(48, candidates.length)),
      commands,
    };
    state.inverse = { ...(state.inverse || {}), maxActuators, lastScore: baseScore, plan, preview: null, physicalValidation: null };
    return plan;
  }

  function analyzeActuatorSensitivity(state, options = {}) {
    const projectedState = planningStateFrom(state);
    const baseSim = RAD.simulate(projectedState);
    const baseScore = objective(baseSim, projectedState);
    const { rows, cols } = projectedState.grid;
    const limits = RAD.commandLimits(projectedState);
    const stepZ = Math.min(options.stepZ ?? 0.12, limits.z);
    const stepAlpha = Math.min(options.stepAlpha ?? 0.12, limits.alphaContract);
    const responseThreshold = options.responseThreshold ?? 0.012;
    const map = RAD.matrix(rows, cols, 0);
    const candidates = [];
    let totalGain = 0;
    let maxGain = 0;

    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (projectedState.cells.removed?.[r]?.[c] || projectedState.cells.locked[r][c] || projectedState.cells.actuatorAllowed?.[r]?.[c] === false) continue;

        const zState = cloneState(projectedState);
        zState.cells.commandZ[r][c] = stepZ;
        const zSim = RAD.simulate(zState);
        const zResponse = rmsDelta(zSim.height, baseSim.height);
        const zScore = objective(zSim, zState);

        const alphaState = cloneState(projectedState);
        alphaState.cells.commandAlpha[r][c] = -stepAlpha;
        const alphaSim = RAD.simulate(alphaState);
        const alphaHeightResponse = rmsDelta(alphaSim.height, baseSim.height);
        const alphaDilationResponse = rmsDelta(alphaSim.alpha, baseSim.alpha);
        const alphaScore = objective(alphaSim, alphaState);

        const zGain = zResponse.rms / Math.max(1e-9, Math.abs(stepZ));
        const alphaGain = alphaHeightResponse.rms / Math.max(1e-9, Math.abs(stepAlpha));
        const dilationGain = alphaDilationResponse.rms / Math.max(1e-9, Math.abs(stepAlpha));
        const combinedGain = Math.hypot(zGain, alphaGain, dilationGain * 0.25);
        const reach = Math.max(
          countResponsiveCells(zSim.height, baseSim.height, responseThreshold),
          countResponsiveCells(alphaSim.height, baseSim.height, responseThreshold),
          countResponsiveCells(alphaSim.alpha, baseSim.alpha, responseThreshold * 0.5)
        );
        const residual = baseSim.target[r][c] - baseSim.height[r][c];
        const bestObjectiveGain = Math.max(0, baseScore - Math.min(zScore, alphaScore));
        map[r][c] = combinedGain;
        totalGain += combinedGain;
        maxGain = Math.max(maxGain, combinedGain);
        candidates.push({
          r,
          c,
          combinedGain,
          zGain,
          alphaGain,
          dilationGain,
          zPeak: zResponse.peak,
          alphaPeak: alphaHeightResponse.peak,
          reach,
          residual,
          objectiveGain: bestObjectiveGain,
          commandZ: RAD.clampCommandZ(projectedState, residual * 0.82 || stepZ),
          commandAlpha: RAD.clampCommandAlpha(projectedState, -Math.sign(residual || 1) * Math.min(0.42, Math.abs(residual || stepAlpha) * 0.24)),
        });
      }
    }

    candidates.sort((a, b) => b.combinedGain - a.combinedGain || b.objectiveGain - a.objectiveGain);
    const sensitivity = {
      strategy: "finite-difference-command-response",
      stepZ,
      stepAlpha,
      responseThreshold,
      map,
      candidates: candidates.slice(0, Math.min(64, candidates.length)),
      controllableCells: candidates.length,
      meanGain: totalGain / Math.max(1, candidates.length),
      maxGain,
      baseScore,
    };
    state.inverse = { ...(state.inverse || {}), sensitivity, preview: null, physicalValidation: null };
    return sensitivity;
  }

  function buildResponseJacobian(state, options = {}) {
    const projectedState = planningStateFrom(state);
    const baseSim = RAD.simulate(projectedState);
    const { rows, cols } = projectedState.grid;
    const limits = RAD.commandLimits(projectedState);
    const stepZ = Math.min(options.stepZ ?? 0.12, limits.z);
    const stepAlpha = Math.min(options.stepAlpha ?? 0.12, limits.alphaContract);
    const residual = targetResidualVector(baseSim);
    const coverageMap = RAD.matrix(rows, cols, 0);
    const columns = [];
    let maxCoverage = 0;
    let meanCoverage = 0;
    let actuatorCount = 0;

    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        if (projectedState.cells.removed?.[r]?.[c] || projectedState.cells.locked[r][c] || projectedState.cells.actuatorAllowed?.[r]?.[c] === false) continue;
        actuatorCount += 1;
        const localColumns = [];

        if (stepZ > 1e-9) {
          const zState = cloneState(projectedState);
          zState.cells.commandZ[r][c] = stepZ;
          const zSim = RAD.simulate(zState);
          const heightDelta = flattenDelta(zSim.height, baseSim.height);
          const alphaDelta = flattenDelta(zSim.alpha, baseSim.alpha);
          localColumns.push({
            type: "z+",
            r,
            c,
            step: stepZ,
            commandZ: stepZ,
            commandAlpha: 0,
            heightDelta: heightDelta.values,
            alphaDelta: alphaDelta.values,
            heightNorm: heightDelta.norm,
            alphaNorm: alphaDelta.norm,
            peakHeight: heightDelta.peak,
            targetAlignment: dot(heightDelta.values, residual) / Math.max(1e-9, heightDelta.norm),
            responsiveCells: countResponsiveCells(zSim.height, baseSim.height, options.responseThreshold ?? 0.012),
          });

          const zNegativeState = cloneState(projectedState);
          zNegativeState.cells.commandZ[r][c] = -stepZ;
          const zNegativeSim = RAD.simulate(zNegativeState);
          const zNegativeHeightDelta = flattenDelta(zNegativeSim.height, baseSim.height);
          const zNegativeAlphaDelta = flattenDelta(zNegativeSim.alpha, baseSim.alpha);
          localColumns.push({
            type: "z-",
            r,
            c,
            step: -stepZ,
            commandZ: -stepZ,
            commandAlpha: 0,
            heightDelta: zNegativeHeightDelta.values,
            alphaDelta: zNegativeAlphaDelta.values,
            heightNorm: zNegativeHeightDelta.norm,
            alphaNorm: zNegativeAlphaDelta.norm,
            peakHeight: zNegativeHeightDelta.peak,
            targetAlignment: dot(zNegativeHeightDelta.values, residual) / Math.max(1e-9, zNegativeHeightDelta.norm),
            responsiveCells: countResponsiveCells(zNegativeSim.height, baseSim.height, options.responseThreshold ?? 0.012),
          });
        }

        if (stepAlpha > 1e-9) {
          const alphaState = cloneState(projectedState);
          alphaState.cells.commandAlpha[r][c] = -stepAlpha;
          const alphaSim = RAD.simulate(alphaState);
          const heightDelta = flattenDelta(alphaSim.height, baseSim.height);
          const alphaDelta = flattenDelta(alphaSim.alpha, baseSim.alpha);
          localColumns.push({
            type: "alpha-",
            r,
            c,
            step: -stepAlpha,
            commandZ: 0,
            commandAlpha: -stepAlpha,
            heightDelta: heightDelta.values,
            alphaDelta: alphaDelta.values,
            heightNorm: heightDelta.norm,
            alphaNorm: alphaDelta.norm,
            peakHeight: heightDelta.peak,
            targetAlignment: dot(heightDelta.values, residual) / Math.max(1e-9, heightDelta.norm),
            responsiveCells: Math.max(
              countResponsiveCells(alphaSim.height, baseSim.height, options.responseThreshold ?? 0.012),
              countResponsiveCells(alphaSim.alpha, baseSim.alpha, (options.responseThreshold ?? 0.012) * 0.5)
            ),
          });

          const alphaExpandState = cloneState(projectedState);
          alphaExpandState.cells.commandAlpha[r][c] = stepAlpha;
          const alphaExpandSim = RAD.simulate(alphaExpandState);
          const alphaExpandHeightDelta = flattenDelta(alphaExpandSim.height, baseSim.height);
          const alphaExpandAlphaDelta = flattenDelta(alphaExpandSim.alpha, baseSim.alpha);
          localColumns.push({
            type: "alpha+",
            r,
            c,
            step: stepAlpha,
            commandZ: 0,
            commandAlpha: stepAlpha,
            heightDelta: alphaExpandHeightDelta.values,
            alphaDelta: alphaExpandAlphaDelta.values,
            heightNorm: alphaExpandHeightDelta.norm,
            alphaNorm: alphaExpandAlphaDelta.norm,
            peakHeight: alphaExpandHeightDelta.peak,
            targetAlignment: dot(alphaExpandHeightDelta.values, residual) / Math.max(1e-9, alphaExpandHeightDelta.norm),
            responsiveCells: Math.max(
              countResponsiveCells(alphaExpandSim.height, baseSim.height, options.responseThreshold ?? 0.012),
              countResponsiveCells(alphaExpandSim.alpha, baseSim.alpha, (options.responseThreshold ?? 0.012) * 0.5)
            ),
          });
        }

        const localCoverage = localColumns.reduce((sum, column) => sum + column.heightNorm + 0.25 * column.alphaNorm, 0);
        coverageMap[r][c] = localCoverage;
        meanCoverage += localCoverage;
        maxCoverage = Math.max(maxCoverage, localCoverage);
        columns.push(...localColumns);
      }
    }

    const targetReachability = targetReachabilityReport(projectedState, baseSim, columns, {
      responseThreshold: options.responseThreshold ?? 0.012,
      targetThreshold: options.targetThreshold,
    });
    const norms = columns.map((column) => column.heightNorm).filter((value) => value > 1e-9).sort((a, b) => a - b);
    const conditionEstimate = norms.length > 1 ? norms[norms.length - 1] / norms[0] : 0;
    columns.sort((a, b) => b.targetAlignment - a.targetAlignment || b.heightNorm - a.heightNorm);
    const jacobian = {
      strategy: "finite-difference-response-jacobian",
      stepZ,
      stepAlpha,
      rows: rows * cols,
      columnCount: columns.length,
      actuatorCount,
      coverageMap,
      meanCoverage: meanCoverage / Math.max(1, actuatorCount),
      maxCoverage,
      conditionEstimate,
      targetReachability,
      baseError: baseSim.metrics.rmsTargetError,
      columns: columns.slice(0, Math.min(96, columns.length)),
    };
    state.inverse = { ...(state.inverse || {}), jacobian, preview: null, physicalValidation: null };
    return jacobian;
  }

  function residualNorm(values) {
    return Math.sqrt(values.reduce((sum, value) => sum + value * value, 0) / Math.max(1, values.length));
  }

  function solveLinearizedTargetFit(state, options = {}) {
    const projectedState = planningStateFrom(state);
    const baseSim = RAD.simulate(projectedState);
    const maxActuators = options.maxActuators ?? state.inverse?.maxActuators ?? 24;
    const maxColumns = Math.max(1, Math.min(options.maxColumns ?? maxActuators * 2, projectedState.grid.rows * projectedState.grid.cols * 4));
    const damping = options.damping ?? 0.018;
    const minGain = options.minGain ?? 1e-5;
    const jacobian = state.inverse?.jacobian?.columns?.length ? state.inverse.jacobian : buildResponseJacobian(state, options);
    jacobian.targetReachability = targetReachabilityReport(projectedState, baseSim, jacobian.columns || [], {
      responseThreshold: options.responseThreshold ?? jacobian.targetReachability?.responseThreshold ?? 0.012,
      targetThreshold: options.targetThreshold ?? jacobian.targetReachability?.targetThreshold,
    });
    const limits = RAD.commandLimits(projectedState);
    const residual = targetResidualVector(baseSim);
    const accepted = [];
    const usedCells = new Set();
    const commandsByCell = new Map();
    const columns = (jacobian.columns || []).filter((column) => column.heightNorm > 1e-9);

    for (let step = 0; step < maxColumns; step += 1) {
      let best = null;
      for (const column of columns) {
        const key = `${column.r},${column.c}`;
        if (!usedCells.has(key) && usedCells.size >= maxActuators) continue;
        const current = commandsByCell.get(key) || { commandZ: 0, commandAlpha: 0 };
        const rawCoefficient = dot(column.heightDelta, residual) / (column.heightNorm ** 2 + damping);
        const coefficient = Math.max(0, rawCoefficient);
        if (coefficient <= 1e-9) continue;
        const trialZ = RAD.clampCommandZ(projectedState, current.commandZ + (column.commandZ || 0) * coefficient);
        const trialAlpha = RAD.clampCommandAlpha(projectedState, current.commandAlpha + (column.commandAlpha || 0) * coefficient);
        const appliedZ = trialZ - current.commandZ;
        const appliedAlpha = trialAlpha - current.commandAlpha;
        const appliedScale =
          Math.abs(column.commandZ || 0) > 1e-9
            ? Math.abs(appliedZ / column.commandZ)
            : Math.abs(column.commandAlpha || 0) > 1e-9
              ? Math.abs(appliedAlpha / column.commandAlpha)
              : 0;
        if (appliedScale <= 1e-9) continue;
        const predictedGain = dot(column.heightDelta, residual) * appliedScale - 0.5 * column.heightNorm ** 2 * appliedScale ** 2;
        if (!best || predictedGain > best.predictedGain) {
          best = { column, key, current, coefficient: appliedScale, trialZ, trialAlpha, appliedZ, appliedAlpha, predictedGain };
        }
      }
      if (!best || best.predictedGain <= minGain) break;
      commandsByCell.set(best.key, { commandZ: best.trialZ, commandAlpha: best.trialAlpha });
      usedCells.add(best.key);
      for (let i = 0; i < residual.length; i += 1) residual[i] -= best.column.heightDelta[i] * best.coefficient;
      accepted.push({
        step: accepted.length + 1,
        r: best.column.r,
        c: best.column.c,
        type: best.column.type,
        coefficient: best.coefficient,
        commandZ: best.trialZ,
        commandAlpha: best.trialAlpha,
        predictedGain: best.predictedGain,
      });
    }

    for (const [key, command] of commandsByCell.entries()) {
      const [r, c] = key.split(",").map(Number);
      projectedState.cells.commandZ[r][c] = command.commandZ;
      projectedState.cells.commandAlpha[r][c] = command.commandAlpha;
    }
    const projectedSim = RAD.simulate(projectedState);
    const solution = {
      strategy: "linearized-jacobian-greedy-fit",
      damping,
      maxActuators,
      baseError: baseSim.metrics.rmsTargetError,
      predictedError: residualNorm(residual),
      projectedError: projectedSim.metrics.rmsTargetError,
      projectedActuators: projectedSim.metrics.recommendedActuators,
      targetReachability: jacobian.targetReachability,
      underactuatedHeightCells: jacobian.targetReachability?.underactuatedHeightCells || 0,
      unreachableHeightRms: jacobian.targetReachability?.unreachableHeightRms || 0,
      maxUnreachableHeightResidual: jacobian.targetReachability?.maxUnreachableHeightResidual || 0,
      steps: accepted.length,
      commands: Array.from(commandsByCell.entries()).map(([key, command], index) => {
        const [r, c] = key.split(",").map(Number);
        return { step: index + 1, r, c, commandZ: command.commandZ, commandAlpha: command.commandAlpha };
      }),
      history: accepted,
    };
    state.inverse = { ...(state.inverse || {}), jacobian, linearSolution: solution, preview: null, physicalValidation: null };
    return solution;
  }

  function applyLinearizedTargetFit(state) {
    const solution = state.inverse?.linearSolution?.commands?.length ? state.inverse.linearSolution : solveLinearizedTargetFit(state);
    if (!state.experiment.initialSnapshot) state.experiment.initialSnapshot = RAD.snapshotState(state);
    const locked = JSON.parse(JSON.stringify(state.cells.locked));
    RAD.clearCommands(state);
    state.cells.locked = locked;
    for (const command of solution.commands) {
      state.cells.commandZ[command.r][command.c] = command.commandZ;
      state.cells.commandAlpha[command.r][command.c] = command.commandAlpha;
    }
    RAD.recordEvent(state, {
      type: "linear-fit-applied",
      actuators: solution.commands.length,
      baseError: solution.baseError,
      projectedError: solution.projectedError,
    });
    state.inverse = { ...(state.inverse || {}), physicalValidation: null };
    return solution;
  }

  function inverseCommandSet(state, source = "auto") {
    const linear = state.inverse?.linearSolution?.commands || [];
    const plan = state.inverse?.plan?.commands || [];
    if (source === "linear") return { source: "linear", commands: linear };
    if (source === "plan") return { source: "plan", commands: plan };
    if (linear.length) return { source: "linear", commands: linear };
    return { source: "plan", commands: plan };
  }

  function validateInversePlanPhysical(state, options = {}) {
    const selected = inverseCommandSet(state, options.source || "auto");
    const commands = selected.commands || [];
    const planningState = planningStateFrom(state);
    const baseKinematic = RAD.simulate(planningState);
    const commandState = cloneState(planningState);
    for (const command of commands) {
      commandState.cells.commandZ[command.r][command.c] = command.commandZ;
      commandState.cells.commandAlpha[command.r][command.c] = command.commandAlpha;
    }
    const projectedKinematic = RAD.simulate(commandState);
    const physicalAvailable = typeof RAD.simulatePhysicalRelaxation === "function";
    const basePhysical = physicalAvailable
      ? RAD.simulatePhysicalRelaxation(planningState, {
          baseSim: baseKinematic,
          iterations: options.iterations ?? 24,
        })
      : baseKinematic;
    const projectedPhysical = physicalAvailable
      ? RAD.simulatePhysicalRelaxation(commandState, {
          baseSim: projectedKinematic,
          iterations: options.iterations ?? 24,
        })
      : projectedKinematic;
    const baseResidual = rmsTargetResidual(basePhysical.target, basePhysical.height);
    const projectedResidual = rmsTargetResidual(projectedPhysical.target, projectedPhysical.height);
    const validation = {
      strategy: "spring-preview-inverse-validation",
      source: selected.source,
      physicalAvailable,
      physicalSuccess: physicalAvailable,
      iterations: projectedPhysical.metrics?.physicalIterations || 0,
      commandCount: commands.length,
      baseKinematicError: baseKinematic.metrics.rmsTargetError,
      projectedKinematicError: projectedKinematic.metrics.rmsTargetError,
      physicalBaseError: baseResidual.rms,
      physicalProjectedError: projectedResidual.rms,
      physicalErrorDelta: baseResidual.rms - projectedResidual.rms,
      physicalMaxAbsResidual: projectedResidual.maxAbs,
      physicalMeanResidual: projectedResidual.bias,
      heightModelRms: projectedPhysical.metrics?.physicalRmsHeightDelta || 0,
      heightModelMax: projectedPhysical.metrics?.physicalMaxHeightDelta || 0,
      centerModelRms: projectedPhysical.metrics?.physicalRmsCenterDelta || 0,
      centerModelMax: projectedPhysical.metrics?.physicalMaxCenterDelta || 0,
      modelAgreementScore: 1 / (1 + (projectedPhysical.metrics?.physicalRmsCenterDelta || 0)),
      commands: commands.map((command) => ({
        r: command.r,
        c: command.c,
        commandZ: command.commandZ,
        commandAlpha: command.commandAlpha,
      })),
    };
    state.inverse = { ...(state.inverse || {}), physicalValidation: validation };
    return validation;
  }

  function inverseDesignReport(state) {
    const sim = RAD.simulateActive(state);
    const inverse = state.inverse || {};
    const jacobian = inverse.jacobian?.columns?.length ? inverse.jacobian : null;
    const linearSolution = inverse.linearSolution?.commands?.length || inverse.linearSolution?.steps ? inverse.linearSolution : null;
    const plan = inverse.plan?.commands?.length || inverse.plan?.candidates?.length ? inverse.plan : null;
    const targetReachability =
      jacobian?.targetReachability ||
      linearSolution?.targetReachability ||
      null;
    let lockedCount = 0;
    let allowedCount = 0;
    let activeCommandCount = 0;
    for (let r = 0; r < state.grid.rows; r += 1) {
      for (let c = 0; c < state.grid.cols; c += 1) {
        if (state.cells.removed?.[r]?.[c]) continue;
        if (state.cells.locked?.[r]?.[c]) lockedCount += 1;
        if (state.cells.actuatorAllowed?.[r]?.[c] !== false) allowedCount += 1;
        if (Math.abs(state.cells.commandAlpha?.[r]?.[c] || 0) > 1e-9 || Math.abs(state.cells.commandZ?.[r]?.[c] || 0) > 1e-9) activeCommandCount += 1;
      }
    }
    return {
      schema: "rad-sim.inverse-design-report.v1",
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        backlash: state.grid.backlash,
        couplingGain: state.grid.couplingGain,
        zCouplingGain: state.grid.zCouplingGain,
        zTravelLimit: state.grid.zTravelLimit,
        alphaContractLimit: state.grid.alphaContractLimit,
        alphaExpandLimit: state.grid.alphaExpandLimit,
      },
      target: cloneData(state.target),
      selection: cloneData(state.selection),
      actuatorState: {
        lockedCount,
        allowedCount,
        activeCommandCount,
        commandAlpha: cloneData(state.cells.commandAlpha),
        commandZ: cloneData(state.cells.commandZ),
        locked: cloneData(state.cells.locked),
        actuatorAllowed: cloneData(state.cells.actuatorAllowed),
        removed: cloneData(state.cells.removed || RAD.matrix(state.grid.rows, state.grid.cols, false)),
      },
      currentMetrics: {
        rmsTargetError: sim.metrics.rmsTargetError,
        meanSignedTargetError: sim.metrics.meanSignedTargetError,
        maxNegativeTargetError: sim.metrics.maxNegativeTargetError,
        maxPositiveTargetError: sim.metrics.maxPositiveTargetError,
        activeCells: sim.metrics.activeCells,
        recommendedActuators: sim.metrics.recommendedActuators,
        maxSaturation: sim.metrics.maxSaturation,
      },
      inverse: {
        plan: cloneData(plan),
        jacobian: cloneData(jacobian),
        linearSolution: cloneData(linearSolution),
        physicalValidation: cloneData(inverse.physicalValidation || null),
        targetReachability: cloneData(targetReachability),
      },
      assumptions: {
        inverseModel: "browser finite-difference columns and greedy linearized target fit",
        physicalValidation: "spring-preview validation is a check of proposed commands, not a calibrated optimizer",
      },
    };
  }

  function exportInverseDesignReport(state) {
    return JSON.stringify(inverseDesignReport(state), null, 2);
  }

  function reachableEquilibriumProfileInverseReport(state, empiricalProfile = null, options = {}) {
    const tolerance = Math.max(0, Number(options.tolerance ?? 1e-9));
    const profile = empiricalProfile || state.experiment?.reachableEquilibriumEmpiricalProfile || {};
    const updates = profileUpdateMap(profile);
    const alphaScale = profileUpdateValue(updates, "alphaResponseScale");
    const heightScale = profileUpdateValue(updates, "heightResponseScale");
    const topologyLeakageTolerance = profileUpdateValue(updates, "topologyLeakageTolerance");
    const groupSequenceTolerance = profileUpdateValue(updates, "groupSequenceTolerance");
    const uncertaintyBudget = profileUpdateValue(updates, "uncertaintyBudget");
    const alphaBand = profileScale(alphaScale, uncertaintyBudget, tolerance);
    let heightBand = profileScale(heightScale, uncertaintyBudget, tolerance);
    if (Number.isFinite(Number(topologyLeakageTolerance))) heightBand = Math.max(heightBand, Math.abs(Number(topologyLeakageTolerance)));

    const baselineState = planningStateFrom(state);
    const baselineSim = RAD.simulate(baselineState);
    const sim = RAD.simulateActive(state);
    const targetHeight = options.targetHeight || sim.target;
    const targetAlpha = options.targetAlpha || null;
    const alphaStats = residualStats(targetAlpha, sim.alpha);
    const heightStats = residualStats(targetHeight, sim.height);
    const alphaTargets = targetCellCount(targetAlpha, baselineSim.alpha, tolerance);
    const heightTargets = targetCellCount(targetHeight, baselineSim.height, tolerance);
    const alphaFailures = bandFailureCount(alphaStats.field, alphaBand, tolerance);
    const heightFailures = bandFailureCount(heightStats.field, heightBand, tolerance);
    const safeProposalCount = Array.from(updates.values()).filter((update) => update.safeToApply === true).length;
    const empiricalReady = Boolean(profile?.summary?.empiricalProfileReady);
    let activeActuatorCount = 0;
    for (let r = 0; r < state.grid.rows; r += 1) {
      for (let c = 0; c < state.grid.cols; c += 1) {
        if (state.cells.removed?.[r]?.[c]) continue;
        if (
          Math.abs(Number(state.cells.commandAlpha?.[r]?.[c]) || 0) > tolerance ||
          Math.abs(Number(state.cells.commandZ?.[r]?.[c]) || 0) > tolerance
        ) {
          activeActuatorCount += 1;
        }
      }
    }
    const inverseSuccess = Boolean(state.inverse?.linearSolution?.commands?.length || state.inverse?.plan?.commands?.length || activeActuatorCount > 0);
    const missingEvidence = [];
    if (!empiricalReady) missingEvidence.push("empiricalProfileReady");
    if (safeProposalCount <= 0) missingEvidence.push("safeBoundedProposals");
    for (const [name, value] of [
      ["alphaResponseScale", alphaScale],
      ["heightResponseScale", heightScale],
      ["topologyLeakageTolerance", topologyLeakageTolerance],
      ["groupSequenceTolerance", groupSequenceTolerance],
      ["uncertaintyBudget", uncertaintyBudget],
    ]) {
      if (!Number.isFinite(Number(value))) missingEvidence.push(name);
    }
    if (!targetAlpha && !targetHeight) missingEvidence.push("inverseTarget");
    if (alphaTargets + heightTargets <= 0) missingEvidence.push("targetCells");
    if (!inverseSuccess) missingEvidence.push("inverseSolveSuccess");
    const scoreTerms = [
      targetAlpha ? alphaStats.rms / alphaBand : null,
      targetHeight ? heightStats.rms / heightBand : null,
      alphaFailures + heightFailures,
    ].filter((value) => Number.isFinite(Number(value)));
    const profileWeightedResidualScore = scoreTerms.length
      ? Math.sqrt(scoreTerms.reduce((sum, value) => sum + Number(value) ** 2, 0) / scoreTerms.length)
      : 0;
    const report = {
      schema: "rad-sim.reachable-equilibrium-profile-inverse.v1",
      method: "profile-aware diagnostic scoring of the current browser inverse state; the empirical profile is read-only calibration metadata",
      profile: {
        schema: profile?.schema || "",
        empiricalProfileReady: empiricalReady,
        safeBoundedProposalCount: safeProposalCount,
        proposalNames: Array.from(updates.keys()).sort(),
      },
      parameters: {
        alphaResponseScale: alphaScale,
        heightResponseScale: heightScale,
        topologyLeakageTolerance,
        groupSequenceTolerance,
        uncertaintyBudget,
        alphaResidualBand: alphaBand,
        heightResidualBand: heightBand,
        tolerance,
      },
      target: {
        targetCellCount: alphaTargets + heightTargets,
        alphaTargetCellCount: alphaTargets,
        heightTargetCellCount: heightTargets,
      },
      inverse: {
        success: inverseSuccess,
        activeActuatorCount,
        projectedError: state.inverse?.linearSolution?.projectedError ?? sim.metrics.rmsTargetError,
        targetReachability: cloneData(state.inverse?.linearSolution?.targetReachability || state.inverse?.jacobian?.targetReachability || null),
      },
      residuals: {
        alphaRmsResidual: alphaStats.rms,
        heightRmsResidual: heightStats.rms,
        maxAbsAlphaResidual: alphaStats.maxAbs,
        maxAbsHeightResidual: heightStats.maxAbs,
        normalizedAlphaRmsResidual: alphaStats.rms / alphaBand,
        normalizedHeightRmsResidual: heightStats.rms / heightBand,
        normalizedMaxAbsAlphaResidual: alphaStats.maxAbs / alphaBand,
        normalizedMaxAbsHeightResidual: heightStats.maxAbs / heightBand,
        alphaProfileBandFailureCount: alphaFailures,
        heightProfileBandFailureCount: heightFailures,
      },
      summary: {
        status: missingEvidence.length ? "needs-profile-aware-inverse-review" : "reachable-equilibrium-profile-inverse-ready",
        profileInverseReady: missingEvidence.length === 0,
        profileWeightedResidualScore,
        profileBandFailureCount: alphaFailures + heightFailures,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "reachable_equilibrium_profile_inverse_gate",
        leanStructure: "Mechanics.ReachableEquilibriumProfileInverseNat",
        leanPredicate: "reachableEquilibriumProfileInverseReadyNat",
        schema: "rad-sim.reachable-equilibrium-profile-inverse.v1",
      },
      claimLabels: {
        profile: "bench-measured empirical law",
        inverse: "browser finite-difference inverse plan",
        diagnostic: "read-only calibration-aware residual certificate",
        physicalAccuracy: "not a physical-law proof and not an optimizer mutation",
      },
      limitations: [
        "The report scores residuals against empirical bands but does not prove reachability.",
        "Browser inverse targets are height-surface targets unless an explicit alpha target grid is provided.",
        "Spring-preview or rigid-body contact validation remains a separate check.",
      ],
    };
    if (options.store === true) state.experiment.reachableEquilibriumProfileInverse = cloneData(report);
    return report;
  }

  function exportReachableEquilibriumProfileInverse(state, empiricalProfile = null, options = {}) {
    return JSON.stringify(reachableEquilibriumProfileInverseReport(state, empiricalProfile, options), null, 2);
  }

  function exportReachableEquilibriumProfileInverseCsv(report) {
    const summary = report.summary || {};
    const profile = report.profile || {};
    const target = report.target || {};
    const inverse = report.inverse || {};
    const residuals = report.residuals || {};
    return [
      [
        "schema",
        "profile_inverse_ready",
        "empirical_profile_ready",
        "target_cells",
        "active_actuators",
        "profile_weighted_residual_score",
        "alpha_rms_residual",
        "height_rms_residual",
        "alpha_band_failures",
        "height_band_failures",
        "missing_evidence",
      ],
      [
        report.schema,
        summary.profileInverseReady,
        profile.empiricalProfileReady,
        target.targetCellCount,
        inverse.activeActuatorCount,
        summary.profileWeightedResidualScore,
        residuals.alphaRmsResidual,
        residuals.heightRmsResidual,
        residuals.alphaProfileBandFailureCount,
        residuals.heightProfileBandFailureCount,
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
      ],
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function numberOrNull(value) {
    const numeric = Number(value);
    return Number.isFinite(numeric) ? numeric : null;
  }

  function reachableEquilibriumProfileInverseAcceptanceReport(profileInverseReport, options = {}) {
    const maxWeightedResidualScore = Math.max(0, Number(options.maxWeightedResidualScore ?? 1.0));
    const maxBandFailures = Math.max(0, Math.round(Number(options.maxBandFailures ?? 0)));
    const maxActiveActuators = Math.max(0, Math.round(Number(options.maxActiveActuators ?? 24)));
    const allowUnderactuated = Boolean(options.allowUnderactuated);
    const requirePhysicalValidation = Boolean(options.requirePhysicalValidation);
    const scoreScale = Math.max(1, Math.round(Number(options.scoreScale ?? 1000)));
    const summary = profileInverseReport?.summary || {};
    const inverse = profileInverseReport?.inverse || {};
    const residuals = profileInverseReport?.residuals || {};
    const physical = profileInverseReport?.physicalValidation || null;
    const score = numberOrNull(summary.profileWeightedResidualScore);
    const bandFailures = Math.max(0, Math.round(Number(summary.profileBandFailureCount || 0)));
    const activeActuators = Math.max(0, Math.round(Number(inverse.activeActuatorCount || 0)));
    const underactuatedCells =
      Math.max(0, Math.round(Number(inverse.alphaUnderactuatedCells || 0))) +
      Math.max(0, Math.round(Number(inverse.heightUnderactuatedCells || 0))) +
      Math.max(0, Math.round(Number(inverse.topologyBlockedAlphaCells || 0))) +
      Math.max(0, Math.round(Number(inverse.topologyBlockedHeightCells || 0)));
    const profileInverseReady = Boolean(summary.profileInverseReady);
    const inverseSuccess = Boolean(inverse.success);
    const physicalPass = physical ? Boolean(physical.physicalSuccess) : !requirePhysicalValidation;
    const missingEvidence = [];
    if (profileInverseReport?.schema !== "rad-sim.reachable-equilibrium-profile-inverse.v1") missingEvidence.push("profileInverseReportSchema");
    if (!profileInverseReady) missingEvidence.push("profileInverseReady");
    if (!inverseSuccess) missingEvidence.push("inverseSolveSuccess");
    if (score === null) missingEvidence.push("profileWeightedResidualScore");
    if (requirePhysicalValidation && !physical) missingEvidence.push("physicalValidation");
    const failedCriteria = [];
    if (score !== null && score > maxWeightedResidualScore) failedCriteria.push("profileWeightedResidualScore");
    if (bandFailures > maxBandFailures) failedCriteria.push("profileBandFailureCount");
    if (activeActuators > maxActiveActuators) failedCriteria.push("activeActuatorBudget");
    if (underactuatedCells > 0 && !allowUnderactuated) failedCriteria.push("underactuatedOrTopologyBlockedCells");
    if (requirePhysicalValidation && !physicalPass) failedCriteria.push("physicalValidation");
    const decision = missingEvidence.length
      ? "reject-missing-evidence"
      : failedCriteria.length
        ? "review-required"
        : "accept-for-preview";
    const accepted = decision === "accept-for-preview";
    return {
      schema: "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1",
      sourceReportSchema: profileInverseReport?.schema || "",
      method: "thresholded acceptance gate for a read-only profile-aware inverse residual certificate",
      criteria: {
        maxWeightedResidualScore,
        maxBandFailures,
        maxActiveActuators,
        allowUnderactuated,
        requirePhysicalValidation,
        scoreScale,
      },
      metrics: {
        profileInverseReady,
        inverseSolveSuccess: inverseSuccess,
        profileWeightedResidualScore: score,
        profileWeightedResidualScoreScaled: Math.ceil((score || 0) * scoreScale),
        maxWeightedResidualScoreScaled: Math.floor(maxWeightedResidualScore * scoreScale),
        profileBandFailureCount: bandFailures,
        activeActuatorCount: activeActuators,
        underactuatedOrTopologyBlockedCellCount: underactuatedCells,
        physicalValidationPass: physicalPass,
        alphaProfileBandFailureCount: Math.max(0, Math.round(Number(residuals.alphaProfileBandFailureCount || 0))),
        heightProfileBandFailureCount: Math.max(0, Math.round(Number(residuals.heightProfileBandFailureCount || 0))),
      },
      decision: {
        decision,
        acceptedForPreview: accepted,
        reviewRequired: !accepted,
        failedCriteria,
      },
      summary: {
        status: accepted ? "profile-aware-inverse-accepted" : "profile-aware-inverse-not-accepted",
        profileInverseAcceptanceReady: accepted,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "reachable_equilibrium_profile_inverse_acceptance_gate",
        leanStructure: "Mechanics.ReachableEquilibriumProfileInverseAcceptanceNat",
        leanPredicate: "reachableEquilibriumProfileInverseAcceptanceReadyNat",
        schema: "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1",
      },
      claimLabels: {
        gate: "finite thresholded acceptance predicate",
        profile: "read-only empirical metadata",
        inverse: "accepted only for preview, not automatic hardware execution",
        physicalAccuracy: "not a physical-law or convergence proof",
      },
      limitations: [
        "Acceptance means the residual certificate fits declared thresholds; it does not prove nonlinear reachability.",
        "Underactuated and topology-blocked targets default to review unless explicitly allowed.",
        "Physical validation remains optional unless required by the caller.",
      ],
    };
  }

  function exportReachableEquilibriumProfileInverseAcceptance(profileInverseReport, options = {}) {
    return JSON.stringify(reachableEquilibriumProfileInverseAcceptanceReport(profileInverseReport, options), null, 2);
  }

  function exportReachableEquilibriumProfileInverseAcceptanceCsv(report) {
    const summary = report.summary || {};
    const decision = report.decision || {};
    const metrics = report.metrics || {};
    return [
      [
        "schema",
        "acceptance_ready",
        "decision",
        "accepted_for_preview",
        "profile_inverse_ready",
        "residual_score",
        "band_failures",
        "active_actuators",
        "underactuated_or_topology_blocked_cells",
        "failed_criteria",
        "missing_evidence",
      ],
      [
        report.schema,
        summary.profileInverseAcceptanceReady,
        decision.decision,
        decision.acceptedForPreview,
        metrics.profileInverseReady,
        metrics.profileWeightedResidualScore,
        metrics.profileBandFailureCount,
        metrics.activeActuatorCount,
        metrics.underactuatedOrTopologyBlockedCellCount,
        Array.isArray(decision.failedCriteria) ? decision.failedCriteria.join(";") : "",
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
      ],
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function activeCommandRecords(state, tolerance = 1e-9) {
    const records = [];
    for (let r = 0; r < state.grid.rows; r += 1) {
      for (let c = 0; c < state.grid.cols; c += 1) {
        if (state.cells.removed?.[r]?.[c]) continue;
        const alpha = Number(state.cells.commandAlpha?.[r]?.[c]) || 0;
        const z = Number(state.cells.commandZ?.[r]?.[c]) || 0;
        if (Math.abs(alpha) <= tolerance && Math.abs(z) <= tolerance) continue;
        records.push({
          index: records.length,
          row: r,
          col: c,
          alpha,
          z,
          commandType: "profile-inverse-actuator-command",
        });
      }
    }
    return records;
  }

  function reachableEquilibriumProfileInversePreviewPacket(state, profileInverseReport = null, acceptanceReport = null, options = {}) {
    const tolerance = Math.max(0, Number(options.tolerance ?? 1e-9));
    const profileInverse = profileInverseReport || state.experiment?.reachableEquilibriumProfileInverse || null;
    const acceptance =
      acceptanceReport ||
      state.experiment?.reachableEquilibriumProfileInverseAcceptance ||
      (profileInverse ? reachableEquilibriumProfileInverseAcceptanceReport(profileInverse, options.acceptanceOptions || {}) : null);
    const commands = activeCommandRecords(state, tolerance);
    const previewEvents = commands.map((command) => ({
      index: command.index,
      type: "set-actuator-command",
      row: command.row,
      col: command.col,
      alpha: command.alpha,
      z: command.z,
      source: "accepted-profile-aware-inverse",
    }));
    const profileSummary = profileInverse?.summary || {};
    const target = profileInverse?.target || {};
    const residuals = profileInverse?.residuals || {};
    const acceptanceSummary = acceptance?.summary || {};
    const acceptanceDecision = acceptance?.decision || {};
    const missingEvidence = [];
    if (!profileInverse || profileInverse.schema !== "rad-sim.reachable-equilibrium-profile-inverse.v1") missingEvidence.push("profileInverseReportSchema");
    if (!acceptance || acceptance.schema !== "rad-sim.reachable-equilibrium-profile-inverse-acceptance.v1") missingEvidence.push("acceptanceReportSchema");
    if (!profileSummary.profileInverseReady) missingEvidence.push("profileInverseReady");
    if (!acceptanceSummary.profileInverseAcceptanceReady) missingEvidence.push("profileInverseAcceptanceReady");
    if (acceptanceDecision.acceptedForPreview !== true) missingEvidence.push("acceptedForPreview");
    if (!commands.length) missingEvidence.push("commandRecords");
    if (Math.max(0, Math.round(Number(target.targetCellCount || 0))) <= 0) missingEvidence.push("targetRecords");
    if (!Number.isFinite(Number(residuals.heightRmsResidual)) && !Number.isFinite(Number(residuals.alphaRmsResidual))) missingEvidence.push("residualRecords");
    const ready = missingEvidence.length === 0;
    const packet = {
      schema: "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1",
      packetId: options.packetId || "profile-inverse-preview",
      sourceReportSchemas: {
        profileInverse: profileInverse?.schema || "",
        acceptance: acceptance?.schema || "",
        inverse: "rad-sim.browser-inverse-state.v1",
      },
      method: "read-only packet of accepted profile-aware inverse commands for visual preview, bench review, or supplemental artifacts",
      acceptance: {
        decision: acceptanceDecision.decision || "",
        acceptedForPreview: Boolean(acceptanceDecision.acceptedForPreview),
        failedCriteria: Array.isArray(acceptanceDecision.failedCriteria) ? acceptanceDecision.failedCriteria : [],
      },
      target: {
        targetCellCount: target.targetCellCount || 0,
        alphaTargetCellCount: target.alphaTargetCellCount || 0,
        heightTargetCellCount: target.heightTargetCellCount || 0,
      },
      residuals: {
        profileWeightedResidualScore: profileSummary.profileWeightedResidualScore,
        profileBandFailureCount: profileSummary.profileBandFailureCount,
        alphaRmsResidual: residuals.alphaRmsResidual,
        heightRmsResidual: residuals.heightRmsResidual,
      },
      commands,
      previewEvents,
      reviewProtocol: {
        mode: "preview-only",
        operator: "accepted-profile-aware-inverse",
        requiresHumanReviewBeforeHardware: true,
        notes: options.notes || "",
        steps: [
          "Load the packet in the browser or notebook.",
          "Inspect target residual, unreachable masks, and actuator budget.",
          "Run spring-preview or rigid-body contact validation before physical actuation.",
          "Record any manual overrides as a separate event sequence.",
        ],
      },
      summary: {
        status: ready ? "profile-aware-inverse-preview-packet-ready" : "profile-aware-inverse-preview-packet-not-ready",
        profileInversePreviewPacketReady: ready,
        commandCount: commands.length,
        eventCount: previewEvents.length,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "reachable_equilibrium_profile_inverse_preview_packet_gate",
        leanStructure: "Mechanics.ReachableEquilibriumProfileInversePreviewPacketNat",
        leanPredicate: "reachableEquilibriumProfileInversePreviewPacketReadyNat",
        schema: "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1",
      },
      claimLabels: {
        packet: "finite preview handoff artifact",
        commands: "accepted inverse command proposal",
        hardware: "requires separate human and physical validation",
        physicalAccuracy: "not a hardware execution proof",
      },
      limitations: [
        "The packet is read-only and does not apply commands.",
        "Accepted preview packets still require physical validation before hardware use.",
        "The commands come from the existing inverse solve and are not re-optimized here.",
      ],
    };
    if (options.store === true) state.experiment.reachableEquilibriumProfileInversePreviewPacket = cloneData(packet);
    return packet;
  }

  function exportReachableEquilibriumProfileInversePreviewPacket(state, profileInverseReport = null, acceptanceReport = null, options = {}) {
    return JSON.stringify(reachableEquilibriumProfileInversePreviewPacket(state, profileInverseReport, acceptanceReport, options), null, 2);
  }

  function exportReachableEquilibriumProfileInversePreviewPacketCsv(packet) {
    const summary = packet.summary || {};
    const target = packet.target || {};
    const acceptance = packet.acceptance || {};
    const residuals = packet.residuals || {};
    const rows = [
      [
        "schema",
        "packet_ready",
        "decision",
        "target_cells",
        "command_count",
        "event_count",
        "profile_weighted_residual_score",
        "profile_band_failures",
        "missing_evidence",
      ],
      [
        packet.schema,
        summary.profileInversePreviewPacketReady,
        acceptance.decision,
        target.targetCellCount,
        summary.commandCount,
        summary.eventCount,
        residuals.profileWeightedResidualScore,
        residuals.profileBandFailureCount,
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
      ],
      ...((packet.commands || []).map((command) => [
        "command",
        command.index,
        command.row,
        command.col,
        command.alpha,
        command.z,
        command.commandType,
        "",
        "",
      ])),
    ];
    return rows.map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function stateFromPreviewPacket(state, packet, tolerance = 1e-9) {
    const replayState = planningStateFrom(state);
    const invalidCommands = [];
    const commands = [];
    const records = Array.isArray(packet?.commands) ? packet.commands : [];
    for (let index = 0; index < records.length; index += 1) {
      const raw = records[index];
      if (!raw || typeof raw !== "object") {
        invalidCommands.push({ index, reason: "command record is not an object" });
        continue;
      }
      const row = Math.trunc(Number(raw.row));
      const col = Math.trunc(Number(raw.col));
      const alpha = Number(raw.alpha);
      const z = Number(raw.z);
      if (!Number.isInteger(row) || !Number.isInteger(col)) {
        invalidCommands.push({ index, reason: "row/col must be integers" });
        continue;
      }
      if (row < 0 || row >= replayState.grid.rows || col < 0 || col >= replayState.grid.cols) {
        invalidCommands.push({ index, row, col, reason: "cell outside grid" });
        continue;
      }
      if (replayState.cells.removed?.[row]?.[col]) {
        invalidCommands.push({ index, row, col, reason: "cell is removed" });
        continue;
      }
      if (replayState.cells.locked?.[row]?.[col]) {
        invalidCommands.push({ index, row, col, reason: "cell is locked" });
        continue;
      }
      if (!Number.isFinite(alpha) && !Number.isFinite(z)) {
        invalidCommands.push({ index, row, col, reason: "missing finite alpha/z command" });
        continue;
      }
      const commandAlpha = Number.isFinite(alpha) ? alpha : 0;
      const commandZ = Number.isFinite(z) ? z : 0;
      if (Math.abs(commandAlpha) <= tolerance && Math.abs(commandZ) <= tolerance) continue;
      replayState.cells.commandAlpha[row][col] += commandAlpha;
      replayState.cells.commandZ[row][col] += commandZ;
      commands.push({ index: commands.length, row, col, alpha: commandAlpha, z: commandZ });
    }
    return { replayState, commands, invalidCommands };
  }

  function finiteCount(matrix) {
    let count = 0;
    for (const row of matrix || []) {
      for (const value of row || []) {
        if (Number.isFinite(Number(value))) count += 1;
      }
    }
    return count;
  }

  function reachableEquilibriumProfileInversePreviewReplayReport(state, packet = null, options = {}) {
    const tolerance = Math.max(0, Number(options.tolerance ?? 1e-9));
    const residualAgreementTolerance = Math.max(0, Number(options.residualAgreementTolerance ?? 1e-9));
    const resolvedPacket = packet || state.experiment?.reachableEquilibriumProfileInversePreviewPacket || null;
    const packetSummary = resolvedPacket?.summary || {};
    const packetResiduals = resolvedPacket?.residuals || {};
    const packetTarget = resolvedPacket?.target || {};
    const { replayState, commands, invalidCommands } = stateFromPreviewPacket(state, resolvedPacket, tolerance);
    const baselineSim = RAD.simulate(planningStateFrom(state));
    const sim = RAD.simulate(replayState);
    const targetHeight = options.targetHeight || sim.target;
    const targetAlpha = options.targetAlpha || null;
    const alphaStats = residualStats(targetAlpha, sim.alpha);
    const heightStats = residualStats(targetHeight, sim.height);
    const packetAlphaRms = numberOrNull(packetResiduals.alphaRmsResidual);
    const packetHeightRms = numberOrNull(packetResiduals.heightRmsResidual);
    const residualDiffs = [];
    if (packetAlphaRms !== null && targetAlpha) residualDiffs.push(Math.abs(alphaStats.rms - packetAlphaRms));
    if (packetHeightRms !== null && targetHeight) residualDiffs.push(Math.abs(heightStats.rms - packetHeightRms));
    const maxResidualDisagreement = residualDiffs.length ? Math.max(...residualDiffs) : 0;
    const comparableResiduals = residualDiffs.length > 0;
    const residualAgreementPass =
      comparableResiduals && maxResidualDisagreement <= Math.max(tolerance, residualAgreementTolerance);
    const eventCount = Array.isArray(resolvedPacket?.previewEvents) ? resolvedPacket.previewEvents.length : 0;
    const expectedCommandCount = Math.max(0, Math.round(Number(packetSummary.commandCount || 0)));
    const commandEventCoveragePass = eventCount >= commands.length;
    let targetCells =
      targetCellCount(targetAlpha, baselineSim.alpha, tolerance) +
      targetCellCount(targetHeight, baselineSim.height, tolerance);
    if (targetCells <= 0) targetCells = Math.max(0, Math.round(Number(packetTarget.targetCellCount || 0)));
    const missingEvidence = [];
    if (!resolvedPacket || resolvedPacket.schema !== "rad-sim.reachable-equilibrium-profile-inverse-preview-packet.v1") missingEvidence.push("packetSchema");
    if (!packetSummary.profileInversePreviewPacketReady) missingEvidence.push("packetReady");
    if (!commands.length) missingEvidence.push("replayedCommands");
    if (invalidCommands.length) missingEvidence.push("validCommandRecords");
    if (expectedCommandCount && expectedCommandCount !== commands.length) missingEvidence.push("commandCountAgreement");
    if (!commandEventCoveragePass) missingEvidence.push("eventCoverage");
    if (!targetAlpha && !targetHeight) missingEvidence.push("targetArrays");
    if (targetCells <= 0) missingEvidence.push("targetRecords");
    if (!comparableResiduals) missingEvidence.push("residualComparison");
    if (comparableResiduals && !residualAgreementPass) missingEvidence.push("residualAgreement");
    const ready = missingEvidence.length === 0;
    const report = {
      schema: "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1",
      sourcePacketSchema: resolvedPacket?.schema || "",
      method: "deterministic kinematic replay of a read-only accepted inverse preview packet",
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        backlash: state.grid.backlash,
        couplingGain: state.grid.couplingGain,
        zCouplingGain: state.grid.zCouplingGain,
      },
      commands: {
        expectedCommandCount,
        replayedCommandCount: commands.length,
        invalidCommandCount: invalidCommands.length,
        invalidCommands,
        eventCount,
        commandEventCoveragePass,
      },
      simulation: {
        model: sim.metrics?.model || "kinematic",
        meanAlpha: sim.metrics?.meanAlpha,
        maxAbsHeight: sim.metrics?.maxAbsHeight,
        finiteDieOffCells: finiteCount(sim.dieOff),
      },
      residuals: {
        targetCellCount: targetCells,
        alphaRmsResidual: targetAlpha ? alphaStats.rms : null,
        heightRmsResidual: targetHeight ? heightStats.rms : null,
        packetAlphaRmsResidual: packetAlphaRms,
        packetHeightRmsResidual: packetHeightRms,
        maxResidualDisagreement,
        residualAgreementPass,
      },
      summary: {
        status: ready ? "profile-aware-inverse-preview-replay-ready" : "profile-aware-inverse-preview-replay-not-ready",
        profileInversePreviewReplayReady: ready,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "reachable_equilibrium_profile_inverse_preview_replay_gate",
        leanStructure: "Mechanics.ReachableEquilibriumProfileInversePreviewReplayNat",
        leanPredicate: "reachableEquilibriumProfileInversePreviewReplayReadyNat",
        schema: "rad-sim.reachable-equilibrium-profile-inverse-preview-replay.v1",
      },
      claimLabels: {
        replay: "deterministic simulator replay certificate",
        packet: "read-only accepted inverse command packet",
        physicalAccuracy: "kinematic replay only; not hardware execution or contact validation",
      },
      limitations: [
        "Replay uses the current browser kinematic simulator and state target.",
        "Residual agreement is only checked when target fields are available.",
        "Passing replay does not imply spring-preview or rigid-body contact validity.",
      ],
    };
    if (options.store === true) state.experiment.reachableEquilibriumProfileInversePreviewReplay = cloneData(report);
    return report;
  }

  function exportReachableEquilibriumProfileInversePreviewReplay(state, packet = null, options = {}) {
    return JSON.stringify(reachableEquilibriumProfileInversePreviewReplayReport(state, packet, options), null, 2);
  }

  function exportReachableEquilibriumProfileInversePreviewReplayCsv(report) {
    const summary = report.summary || {};
    const commands = report.commands || {};
    const residuals = report.residuals || {};
    return [
      [
        "schema",
        "replay_ready",
        "replayed_commands",
        "invalid_commands",
        "event_count",
        "target_cells",
        "height_rms_residual",
        "max_residual_disagreement",
        "residual_agreement_pass",
        "missing_evidence",
      ],
      [
        report.schema,
        summary.profileInversePreviewReplayReady,
        commands.replayedCommandCount,
        commands.invalidCommandCount,
        commands.eventCount,
        residuals.targetCellCount,
        residuals.heightRmsResidual,
        residuals.maxResidualDisagreement,
        residuals.residualAgreementPass,
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
      ],
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function heightModelComparisonStats(physicalHeight, kinematicHeight) {
    if (!Array.isArray(physicalHeight) || !Array.isArray(kinematicHeight)) {
      return { rms: 0, maxAbs: 0, finite: false };
    }
    let sum = 0;
    let count = 0;
    let maxAbs = 0;
    for (let r = 0; r < physicalHeight.length; r += 1) {
      for (let c = 0; c < physicalHeight[r].length; c += 1) {
        const delta = Number(physicalHeight[r][c]) - Number(kinematicHeight?.[r]?.[c]);
        if (!Number.isFinite(delta)) return { rms: 0, maxAbs: 0, finite: false };
        sum += delta ** 2;
        maxAbs = Math.max(maxAbs, Math.abs(delta));
        count += 1;
      }
    }
    return { rms: Math.sqrt(sum / Math.max(1, count)), maxAbs, finite: true };
  }

  function centerModelComparisonStats(physicalCenters, kinematicCenters) {
    if (!Array.isArray(physicalCenters) || !Array.isArray(kinematicCenters)) {
      return { rms: 0, maxAbs: 0, finite: false };
    }
    let sum = 0;
    let count = 0;
    let maxAbs = 0;
    for (let r = 0; r < physicalCenters.length; r += 1) {
      for (let c = 0; c < physicalCenters[r].length; c += 1) {
        const physical = physicalCenters[r][c] || {};
        const kinematic = kinematicCenters?.[r]?.[c] || {};
        const delta = Math.hypot(
          Number(physical.x) - Number(kinematic.x),
          Number(physical.y) - Number(kinematic.y),
          Number(physical.z) - Number(kinematic.z)
        );
        if (!Number.isFinite(delta)) return { rms: 0, maxAbs: 0, finite: false };
        sum += delta ** 2;
        maxAbs = Math.max(maxAbs, delta);
        count += 1;
      }
    }
    return { rms: Math.sqrt(sum / Math.max(1, count)), maxAbs, finite: true };
  }

  function reachableEquilibriumProfileInversePreviewPhysicalReport(state, packet = null, options = {}) {
    const tolerance = Math.max(0, Number(options.tolerance ?? 1e-9));
    const replay = reachableEquilibriumProfileInversePreviewReplayReport(state, packet, {
      targetAlpha: options.targetAlpha,
      targetHeight: options.targetHeight,
      tolerance,
      residualAgreementTolerance: options.residualAgreementTolerance ?? 1e-9,
    });
    const resolvedPacket = packet || state.experiment?.reachableEquilibriumProfileInversePreviewPacket || null;
    const { replayState, commands, invalidCommands } = stateFromPreviewPacket(state, resolvedPacket, tolerance);
    const kinematic = RAD.simulate(replayState);
    const physicalAvailable = typeof RAD.simulatePhysicalRelaxation === "function";
    const physical = physicalAvailable
      ? RAD.simulatePhysicalRelaxation(replayState, {
          baseSim: kinematic,
          iterations: options.iterations,
          springGain: options.springGain,
          anchorGain: options.anchorGain,
        })
      : null;
    const targetHeight = options.targetHeight || kinematic.target;
    const targetAlpha = options.targetAlpha || null;
    const heightResidual = physical ? residualStats(targetHeight, physical.height) : { rms: null, maxAbs: null };
    const alphaResidual = physical ? residualStats(targetAlpha, physical.alpha) : { rms: null, maxAbs: null };
    const heightStats = physical ? heightModelComparisonStats(physical.height, kinematic.height) : { rms: 0, maxAbs: 0, finite: false };
    const centerStats = physical ? centerModelComparisonStats(physical.centers, kinematic.centers) : { rms: 0, maxAbs: 0, finite: false };
    let targetCells =
      targetCellCount(targetAlpha, kinematic.alpha, tolerance) +
      targetCellCount(targetHeight, kinematic.height, tolerance);
    if (targetCells <= 0) targetCells = Math.max(0, Math.round(Number(resolvedPacket?.target?.targetCellCount || 0)));
    const maxHeightModelError = Number.isFinite(Number(options.maxHeightModelError)) ? Number(options.maxHeightModelError) : null;
    const maxCenterModelError = Number.isFinite(Number(options.maxCenterModelError)) ? Number(options.maxCenterModelError) : null;
    let modelErrorPass = true;
    if (maxHeightModelError !== null) modelErrorPass = modelErrorPass && heightStats.maxAbs <= maxHeightModelError + tolerance;
    if (maxCenterModelError !== null) modelErrorPass = modelErrorPass && centerStats.maxAbs <= maxCenterModelError + tolerance;
    const physicalSuccess = Boolean(physical?.metrics?.physicalPreview);
    const finiteModelComparison = Boolean(heightStats.finite && centerStats.finite);
    const missingEvidence = [];
    if (!replay.summary?.profileInversePreviewReplayReady) missingEvidence.push("previewReplayReady");
    if (!commands.length) missingEvidence.push("commandRecords");
    if (invalidCommands.length) missingEvidence.push("validCommandRecords");
    if (!physicalSuccess) missingEvidence.push("physicalSolverSuccess");
    if (!targetAlpha && !targetHeight) missingEvidence.push("targetArrays");
    if (targetCells <= 0) missingEvidence.push("targetRecords");
    if (!finiteModelComparison) missingEvidence.push("modelComparisonRecords");
    if (!Number.isFinite(Number(physical?.metrics?.meanAbsLinkStrain))) missingEvidence.push("energyRecords");
    if (!modelErrorPass) missingEvidence.push("modelErrorThreshold");
    const ready = missingEvidence.length === 0;
    const report = {
      schema: "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1",
      sourcePacketSchema: resolvedPacket?.schema || "",
      sourceReplaySchema: replay.schema || "",
      method: "read-only spring-preview physical check of an accepted inverse command packet",
      grid: {
        rows: state.grid.rows,
        cols: state.grid.cols,
        backlash: state.grid.backlash,
        couplingGain: state.grid.couplingGain,
        zCouplingGain: state.grid.zCouplingGain,
      },
      commands: {
        replayedCommandCount: commands.length,
        invalidCommandCount: invalidCommands.length,
      },
      physical: {
        model: physical?.metrics?.model || "spring-preview",
        physicalAvailable,
        physicalSuccess,
        iterations: physical?.metrics?.physicalIterations || 0,
        physicalEnergyProxy: physical?.metrics?.meanAbsLinkStrain ?? null,
        storedEnergy: null,
        springEdges: physical?.metrics?.physicalActiveSpringEdges || 0,
        skippedSpringEdges: physical?.metrics?.physicalSkippedSpringEdges || 0,
        removedCells: physical?.metrics?.removedCells || 0,
      },
      comparison: {
        heightRmsModelError: heightStats.rms,
        centerRmsModelError: centerStats.rms,
        maxAbsHeightModelError: heightStats.maxAbs,
        maxAbsCenterModelError: centerStats.maxAbs,
        modelAgreementScore: 1 / (1 + centerStats.rms),
        finiteModelComparison,
        modelErrorThresholdPass: modelErrorPass,
        maxHeightModelErrorLimit: maxHeightModelError,
        maxCenterModelErrorLimit: maxCenterModelError,
      },
      residuals: {
        targetCellCount: targetCells,
        physicalAlphaRmsResidual: targetAlpha ? alphaResidual.rms : null,
        physicalHeightRmsResidual: targetHeight ? heightResidual.rms : null,
        physicalMaxAbsHeightResidual: targetHeight ? heightResidual.maxAbs : null,
      },
      summary: {
        status: ready ? "profile-aware-inverse-preview-physical-ready" : "profile-aware-inverse-preview-physical-not-ready",
        profileInversePreviewPhysicalReady: ready,
        missingEvidenceCount: missingEvidence.length,
        missingEvidence,
      },
      formalization: {
        targetId: "reachable_equilibrium_profile_inverse_preview_physical_gate",
        leanStructure: "Mechanics.ReachableEquilibriumProfileInversePreviewPhysicalNat",
        leanPredicate: "reachableEquilibriumProfileInversePreviewPhysicalReadyNat",
        schema: "rad-sim.reachable-equilibrium-profile-inverse-preview-physical.v1",
      },
      claimLabels: {
        gate: "finite physical-preview evidence predicate",
        physics: "simulator-derived spring preview",
        physicalAccuracy: "experimentally unvalidated physical assumption",
        hardware: "not hardware execution authorization",
      },
      limitations: [
        "The browser spring preview is normalized and uncalibrated.",
        "Passing this gate does not prove rigid-body contact, friction, gravity, or material accuracy.",
        "The report replays a command packet and does not re-optimize inverse design commands.",
      ],
    };
    if (options.store === true) state.experiment.reachableEquilibriumProfileInversePreviewPhysical = cloneData(report);
    return report;
  }

  function exportReachableEquilibriumProfileInversePreviewPhysical(state, packet = null, options = {}) {
    return JSON.stringify(reachableEquilibriumProfileInversePreviewPhysicalReport(state, packet, options), null, 2);
  }

  function exportReachableEquilibriumProfileInversePreviewPhysicalCsv(report) {
    const summary = report.summary || {};
    const commands = report.commands || {};
    const physical = report.physical || {};
    const comparison = report.comparison || {};
    const residuals = report.residuals || {};
    return [
      [
        "schema",
        "physical_ready",
        "physical_success",
        "replayed_commands",
        "target_cells",
        "physical_height_rms_residual",
        "height_rms_model_error",
        "center_rms_model_error",
        "physical_energy_proxy",
        "missing_evidence",
      ],
      [
        report.schema,
        summary.profileInversePreviewPhysicalReady,
        physical.physicalSuccess,
        commands.replayedCommandCount,
        residuals.targetCellCount,
        residuals.physicalHeightRmsResidual,
        comparison.heightRmsModelError,
        comparison.centerRmsModelError,
        physical.physicalEnergyProxy,
        Array.isArray(summary.missingEvidence) ? summary.missingEvidence.join(";") : "",
      ],
    ].map((row) => row.map(csvValue).join(",")).join("\n");
  }

  function previewFromCommands(state, commands, meta = {}) {
    const previewState = planningStateFrom(state);
    const baseSim = RAD.simulate(previewState);
    for (const command of commands) {
      previewState.cells.commandZ[command.r][command.c] = command.commandZ;
      previewState.cells.commandAlpha[command.r][command.c] = command.commandAlpha;
    }
    const previewSim = RAD.simulate(previewState);
    const { rows, cols } = state.grid;
    const contribution = RAD.matrix(rows, cols, 0);
    const errorReduction = RAD.matrix(rows, cols, 0);
    let maxContribution = 0;
    let totalErrorReduction = 0;
    for (let r = 0; r < rows; r += 1) {
      for (let c = 0; c < cols; c += 1) {
        contribution[r][c] = previewSim.height[r][c] - baseSim.height[r][c];
        errorReduction[r][c] = Math.abs(baseSim.targetError[r][c]) - Math.abs(previewSim.targetError[r][c]);
        maxContribution = Math.max(maxContribution, Math.abs(contribution[r][c]));
        totalErrorReduction += errorReduction[r][c];
      }
    }
    return {
      ...meta,
      commands: commands.map((command) => ({
        r: command.r,
        c: command.c,
        commandZ: command.commandZ,
        commandAlpha: command.commandAlpha,
      })),
      baseError: baseSim.metrics.rmsTargetError,
      previewError: previewSim.metrics.rmsTargetError,
      errorDelta: baseSim.metrics.rmsTargetError - previewSim.metrics.rmsTargetError,
      totalErrorReduction,
      maxContribution,
      contribution,
      errorReduction,
    };
  }

  function setInversePreview(state, candidate) {
    if (!candidate) {
      state.inverse = { ...(state.inverse || {}), preview: null };
      return null;
    }
    const preview = previewFromCommands(state, [candidate], {
      type: "candidate",
      r: candidate.r,
      c: candidate.c,
      commandZ: candidate.commandZ,
      commandAlpha: candidate.commandAlpha,
      improvement: candidate.improvement || 0,
      residual: candidate.residual || 0,
    });
    state.inverse = { ...(state.inverse || {}), preview };
    return preview;
  }

  function setInversePlanStepPreview(state, stepCount) {
    const plan = state.inverse?.plan?.commands?.length ? state.inverse.plan : buildInverseDesignPlan(state);
    const count = Math.max(0, Math.min(plan.commands.length, Number(stepCount) || 0));
    if (count === 0) return setInversePreview(state, null);
    const commands = plan.commands.slice(0, count);
    const last = commands[commands.length - 1];
    const preview = previewFromCommands(state, commands, {
      type: "plan-step",
      step: count,
      r: last.r,
      c: last.c,
      commandZ: last.commandZ,
      commandAlpha: last.commandAlpha,
      improvement: commands.reduce((sum, command) => sum + Math.max(0, command.improvement || 0), 0),
      residual: last.residual || 0,
    });
    state.inverse = { ...(state.inverse || {}), preview };
    return preview;
  }

  function applyInverseDesignPlan(state) {
    const plan = state.inverse?.plan?.commands?.length ? state.inverse.plan : buildInverseDesignPlan(state);
    if (!state.experiment.initialSnapshot) state.experiment.initialSnapshot = RAD.snapshotState(state);
    const locked = JSON.parse(JSON.stringify(state.cells.locked));
    RAD.clearCommands(state);
    state.cells.locked = locked;
    for (const command of plan.commands) {
      state.cells.commandZ[command.r][command.c] = command.commandZ;
      state.cells.commandAlpha[command.r][command.c] = command.commandAlpha;
    }
    RAD.recordEvent(state, {
      type: "inverse-plan-applied",
      target: state.target.type,
      actuators: plan.commands.length,
      score: plan.baseScore,
    });
    state.inverse = { ...(state.inverse || {}), physicalValidation: null };
    return plan;
  }

  RAD.buildInverseDesignPlan = buildInverseDesignPlan;
  RAD.analyzeActuatorSensitivity = analyzeActuatorSensitivity;
  RAD.buildResponseJacobian = buildResponseJacobian;
  RAD.solveLinearizedTargetFit = solveLinearizedTargetFit;
  RAD.applyLinearizedTargetFit = applyLinearizedTargetFit;
  RAD.validateInversePlanPhysical = validateInversePlanPhysical;
  RAD.inverseDesignReport = inverseDesignReport;
  RAD.exportInverseDesignReport = exportInverseDesignReport;
  RAD.reachableEquilibriumProfileInverseReport = reachableEquilibriumProfileInverseReport;
  RAD.exportReachableEquilibriumProfileInverse = exportReachableEquilibriumProfileInverse;
  RAD.exportReachableEquilibriumProfileInverseCsv = exportReachableEquilibriumProfileInverseCsv;
  RAD.reachableEquilibriumProfileInverseAcceptanceReport = reachableEquilibriumProfileInverseAcceptanceReport;
  RAD.exportReachableEquilibriumProfileInverseAcceptance = exportReachableEquilibriumProfileInverseAcceptance;
  RAD.exportReachableEquilibriumProfileInverseAcceptanceCsv = exportReachableEquilibriumProfileInverseAcceptanceCsv;
  RAD.reachableEquilibriumProfileInversePreviewPacket = reachableEquilibriumProfileInversePreviewPacket;
  RAD.exportReachableEquilibriumProfileInversePreviewPacket = exportReachableEquilibriumProfileInversePreviewPacket;
  RAD.exportReachableEquilibriumProfileInversePreviewPacketCsv = exportReachableEquilibriumProfileInversePreviewPacketCsv;
  RAD.reachableEquilibriumProfileInversePreviewReplayReport = reachableEquilibriumProfileInversePreviewReplayReport;
  RAD.exportReachableEquilibriumProfileInversePreviewReplay = exportReachableEquilibriumProfileInversePreviewReplay;
  RAD.exportReachableEquilibriumProfileInversePreviewReplayCsv = exportReachableEquilibriumProfileInversePreviewReplayCsv;
  RAD.reachableEquilibriumProfileInversePreviewPhysicalReport = reachableEquilibriumProfileInversePreviewPhysicalReport;
  RAD.exportReachableEquilibriumProfileInversePreviewPhysical = exportReachableEquilibriumProfileInversePreviewPhysical;
  RAD.exportReachableEquilibriumProfileInversePreviewPhysicalCsv = exportReachableEquilibriumProfileInversePreviewPhysicalCsv;
  RAD.setInversePreview = setInversePreview;
  RAD.setInversePlanStepPreview = setInversePlanStepPreview;
  RAD.applyInverseDesignPlan = applyInverseDesignPlan;
})();
