(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
  }

  function cloneState(state) {
    return JSON.parse(JSON.stringify(state));
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

  function targetReachabilityReport(state, baseSim, columns, options = {}) {
    const { rows, cols } = state.grid;
    const responseThreshold = options.responseThreshold ?? 0.012;
    const targetThreshold = options.targetThreshold ?? responseThreshold;
    const heightReachableMap = RAD.matrix(rows, cols, false);
    const targetHeightMap = RAD.matrix(rows, cols, 0);
    const underactuatedHeightMap = RAD.matrix(rows, cols, 0);
    let heightReachableCells = 0;
    let targetHeightCells = 0;
    let underactuatedHeightCells = 0;
    let unreachableSquared = 0;
    let maxUnreachableHeightResidual = 0;
    let worstUnderactuatedCell = null;

    for (let index = 0; index < rows * cols; index += 1) {
      const r = Math.floor(index / cols);
      const c = index % cols;
      const reachable = columns.some((column) => Math.abs(column.heightDelta?.[index] || 0) >= responseThreshold);
      heightReachableMap[r][c] = reachable;
      if (reachable) heightReachableCells += 1;
      const residual = (baseSim.target?.[r]?.[c] || 0) - (baseSim.height?.[r]?.[c] || 0);
      const requested = Math.abs(residual) >= targetThreshold;
      if (requested) {
        targetHeightCells += 1;
        targetHeightMap[r][c] = residual;
      }
      if (requested && !reachable) {
        const magnitude = Math.abs(residual);
        underactuatedHeightCells += 1;
        underactuatedHeightMap[r][c] = magnitude;
        unreachableSquared += residual * residual;
        if (magnitude > maxUnreachableHeightResidual) {
          maxUnreachableHeightResidual = magnitude;
          worstUnderactuatedCell = { row: r, col: c, residual, absResidual: magnitude };
        }
      }
    }

    return {
      model: "finite-response-height-reachability",
      responseThreshold,
      targetThreshold,
      heightReachableMap,
      targetHeightMap,
      underactuatedHeightMap,
      heightReachableCells,
      targetHeightCells,
      underactuatedHeightCells,
      unreachableHeightRms: Math.sqrt(unreachableSquared / Math.max(1, underactuatedHeightCells)),
      maxUnreachableHeightResidual,
      worstUnderactuatedCell,
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
          if (projectedState.cells.locked[r][c] || projectedState.cells.actuatorAllowed?.[r]?.[c] === false || selected.has(key)) continue;
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
        if (projectedState.cells.locked[r][c] || projectedState.cells.actuatorAllowed?.[r]?.[c] === false || selected.has(key)) continue;
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
        if (projectedState.cells.locked[r][c] || projectedState.cells.actuatorAllowed?.[r]?.[c] === false) continue;

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
        if (projectedState.cells.locked[r][c] || projectedState.cells.actuatorAllowed?.[r]?.[c] === false) continue;
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
  RAD.setInversePreview = setInversePreview;
  RAD.setInversePlanStepPreview = setInversePlanStepPreview;
  RAD.applyInverseDesignPlan = applyInverseDesignPlan;
})();
