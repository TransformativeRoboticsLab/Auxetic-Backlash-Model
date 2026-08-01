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
})();
