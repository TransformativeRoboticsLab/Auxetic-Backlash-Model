(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  const CORE_DOM_IDS = [
    "threeMount",
    "statusLine",
    "runSystemCheck",
    "systemCheckPanel",
    "systemCheckSummary",
    "systemCheckList",
    "alphaCommand",
    "zCommand",
    "locked",
    "positionLocked",
    "applyCommand",
    "operatorCommitLock",
    "cellGrid",
    "overlayMode",
    "simulationMode",
  ];

  const CORE_APIS = [
    "createState",
    "simulate",
    "simulateActive",
    "serialize",
    "deserialize",
    "applyPreset",
    "applyEventSequence",
    "localActuationEvent",
    "lockEvent",
    "clearActuationEvent",
    "referenceCenter",
    "matrix",
  ];

  function byId(id, root = document) {
    return root.getElementById ? root.getElementById(id) : root.querySelector(`#${id}`);
  }

  function requireIds(ids, root = document) {
    const missing = [];
    const elements = {};
    for (const id of ids) {
      const element = byId(id, root);
      if (element) elements[id] = element;
      else missing.push(id);
    }
    return { ok: missing.length === 0, missing, elements };
  }

  function setPressed(element, active) {
    if (!element) return;
    element.classList.toggle("is-active", Boolean(active));
    element.setAttribute("aria-pressed", String(Boolean(active)));
  }

  function setText(element, value) {
    if (element) element.textContent = String(value);
  }

  function clampNumber(value, min, max, fallback = 0) {
    const numeric = Number(value);
    if (!Number.isFinite(numeric)) return fallback;
    return Math.max(min, Math.min(max, numeric));
  }

  function cloneData(value) {
    return JSON.parse(JSON.stringify(value));
  }

  function installCollapsibleSections(options = {}) {
    const root = options.root || document;
    const initiallyCollapsed = new Set(options.initiallyCollapsed || []);
    root.querySelectorAll(".control-panel section").forEach((section) => {
      const heading = section.querySelector("h2");
      if (!heading || heading.querySelector("button")) return;
      const title = heading.textContent.trim();
      section.classList.add("control-section");
      heading.classList.add("section-heading");
      const button = document.createElement("button");
      button.type = "button";
      button.className = "section-toggle";
      button.setAttribute("aria-label", `Toggle ${title} controls`);
      heading.appendChild(button);
      const setCollapsed = (collapsed) => {
        section.classList.toggle("is-collapsed", collapsed);
        button.textContent = collapsed ? "+" : "-";
        button.setAttribute("aria-expanded", String(!collapsed));
      };
      setCollapsed(initiallyCollapsed.has(title));
      heading.addEventListener("click", (event) => {
        if (event.target.closest("button") || event.target === heading) {
          setCollapsed(!section.classList.contains("is-collapsed"));
        }
      });
      button.addEventListener("click", (event) => {
        event.stopPropagation();
        setCollapsed(!section.classList.contains("is-collapsed"));
      });
    });
  }

  function runCheck(name, fn) {
    const started = performance.now();
    try {
      const detail = fn();
      return {
        name,
        ok: detail === undefined ? true : detail !== false,
        detail: typeof detail === "string" ? detail : "",
        durationMs: performance.now() - started,
      };
    } catch (error) {
      return {
        name,
        ok: false,
        detail: error.message || String(error),
        durationMs: performance.now() - started,
      };
    }
  }

  function assert(condition, message) {
    if (!condition) throw new Error(message);
  }

  function checkCoreApis() {
    const missing = CORE_APIS.filter((name) => typeof RAD[name] !== "function");
    assert(!missing.length, `missing APIs: ${missing.join(", ")}`);
    return `${CORE_APIS.length} APIs`;
  }

  function checkSimulationSmoke() {
    const state = RAD.createState(4, 4);
    state.cells.commandAlpha[1][1] = -0.25;
    state.cells.commandZ[1][1] = 0.2;
    const sim = RAD.simulateActive(state);
    assert(sim.alpha.length === 4 && sim.alpha[0].length === 4, "wrong alpha shape");
    assert(Number.isFinite(sim.height[1][1]), "height is not finite");
    assert(Number.isFinite(sim.metrics.meanAlpha), "mean alpha metric is not finite");
    return `mean alpha ${sim.metrics.meanAlpha.toFixed(3)}`;
  }

  function checkSerializationRoundtrip() {
    const state = RAD.createState(3, 3);
    state.cells.commandAlpha[1][1] = -0.31;
    state.cells.commandZ[1][1] = 0.27;
    state.cells.positionLocked[0][0] = true;
    const restored = RAD.deserialize(RAD.serialize(state));
    assert(restored.grid.rows === 3 && restored.grid.cols === 3, "grid did not roundtrip");
    assert(Math.abs(restored.cells.commandAlpha[1][1] + 0.31) < 1e-9, "alpha command did not roundtrip");
    assert(restored.cells.positionLocked[0][0] === true, "position lock did not roundtrip");
    return "JSON schema ok";
  }

  function checkEventLockHeight() {
    const state = RAD.createState(3, 3);
    RAD.clearCommands(state);
    const preLock = RAD.applyEventSequence(state, [
      RAD.localActuationEvent({ r: 1, c: 1 }, -0.3, 0.2),
    ]);
    const expectedHeight = RAD.simulate(preLock).height[1][1];
    const locked = RAD.applyEventSequence(state, [
      RAD.localActuationEvent({ r: 1, c: 1 }, -0.3, 0.2),
      RAD.lockEvent({ r: 1, c: 1 }),
      RAD.clearActuationEvent(),
    ]);
    const sim = RAD.simulate(locked);
    assert(locked.cells.locked[1][1] === true, "cell is not locked");
    assert(Math.abs(locked.cells.lockZ[1][1] - expectedHeight) < 1e-9, "lockZ did not commit current height");
    assert(Math.abs(sim.height[1][1] - expectedHeight) < 1e-9, "locked height changed after commands cleared");
    assert(Math.abs(sim.height[1][1]) > 1e-6, "locked height snapped to ground");
    return `z ${sim.height[1][1].toFixed(3)}`;
  }

  function checkPositionFixture() {
    const state = RAD.createState(3, 3);
    state.view.simulationMode = "springPreview";
    state.cells.positionLocked[1][1] = true;
    state.cells.commandAlpha[1][1] = -0.4;
    state.cells.commandZ[1][1] = 0.5;
    const sim = RAD.simulateActive(state);
    const ref = RAD.referenceCenter(state, 1, 1);
    const center = sim.centers[1][1];
    assert(Math.abs(center.x - ref.x) < 1e-9, "fixture x moved");
    assert(Math.abs(center.y - ref.y) < 1e-9, "fixture y moved");
    assert(Math.abs(center.z - ref.z) < 1e-9, "fixture z moved");
    return "fixture held";
  }

  function checkPresetDiversity() {
    const signatures = ["center", "dome", "saddle", "wave"].map((name) => {
      const state = RAD.createState(5, 5);
      RAD.applyPreset(state, name);
      const sim = RAD.simulate(state);
      return `${name}:${sim.height.flat().map((value) => value.toFixed(3)).join(",")}`;
    });
    assert(new Set(signatures).size > 2, "presets are not producing distinct states");
    return `${signatures.length} presets`;
  }

  function checkRendererContext(context) {
    const renderer = context?.renderer;
    assert(renderer && typeof renderer.getDiagnostics === "function", "renderer diagnostics unavailable");
    const diagnostics = renderer.getDiagnostics();
    assert(Number.isFinite(Number(diagnostics.drawCalls || 0)), "draw call metric invalid");
    return `${diagnostics.drawCalls || 0} draw calls`;
  }

  function runSelfCheck(context = {}) {
    const dom = requireIds(CORE_DOM_IDS);
    const checks = [
      runCheck("DOM anchors", () => {
        assert(dom.ok, `missing DOM ids: ${dom.missing.join(", ")}`);
        return `${CORE_DOM_IDS.length} anchors`;
      }),
      runCheck("Core APIs", checkCoreApis),
      runCheck("Simulation smoke", checkSimulationSmoke),
      runCheck("JSON roundtrip", checkSerializationRoundtrip),
      runCheck("Event lock preserves height", checkEventLockHeight),
      runCheck("Position lock fixture", checkPositionFixture),
      runCheck("Preset diversity", checkPresetDiversity),
      runCheck("Renderer diagnostics", () => checkRendererContext(context)),
    ];
    const passed = checks.filter((check) => check.ok).length;
    return {
      schema: "rad-sim.browser-self-check.v1",
      at: new Date().toISOString(),
      passed,
      total: checks.length,
      ok: passed === checks.length,
      checks,
      summary: `${passed}/${checks.length} checks passed`,
    };
  }

  function renderSelfCheck(panel, summary, list, result) {
    if (!panel || !summary || !list) return;
    panel.hidden = false;
    panel.classList.toggle("is-ok", result.ok);
    panel.classList.toggle("is-failing", !result.ok);
    summary.textContent = result.ok ? `System check passed: ${result.summary}` : `System check failed: ${result.summary}`;
    list.replaceChildren(
      ...result.checks.map((check) => {
        const item = document.createElement("li");
        item.className = check.ok ? "is-ok" : "is-failing";
        const label = document.createElement("span");
        label.textContent = check.name;
        const detail = document.createElement("small");
        detail.textContent = check.detail || `${check.durationMs.toFixed(1)} ms`;
        item.append(label, detail);
        return item;
      })
    );
  }

  function installSelfCheckPanel(options = {}) {
    const button = options.button || byId("runSystemCheck");
    const panel = options.panel || byId("systemCheckPanel");
    const summary = options.summary || byId("systemCheckSummary");
    const list = options.list || byId("systemCheckList");
    const getContext = options.getContext || (() => ({}));
    let lastResult = null;
    const run = () => {
      lastResult = runSelfCheck(getContext());
      renderSelfCheck(panel, summary, list, lastResult);
      if (button) {
        setPressed(button, lastResult.ok);
        button.textContent = lastResult.ok ? "Checks OK" : "Check Failed";
        button.title = lastResult.summary;
      }
      document.body.classList.toggle("system-check-ok", lastResult.ok);
      document.body.classList.toggle("system-check-failing", !lastResult.ok);
      return lastResult;
    };
    if (button) button.addEventListener("click", run);
    return {
      run,
      get lastResult() {
        return lastResult;
      },
    };
  }

  RAD.Primitives = {
    CORE_DOM_IDS,
    CORE_APIS,
    byId,
    requireIds,
    setPressed,
    setText,
    clampNumber,
    cloneData,
    installCollapsibleSections,
    runSelfCheck,
    installSelfCheckPanel,
  };
})();
