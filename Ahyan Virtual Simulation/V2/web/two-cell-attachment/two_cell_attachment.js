(function () {
  "use strict";

  const THREE = window.THREE;
  const mount = document.getElementById("threeMount");
  const aTopInput = document.getElementById("aTopAngle");
  const actuationModeInput = document.getElementById("actuationMode");
  const animateInput = document.getElementById("animate");
  const animationHalfPeriodSecondsInput = document.getElementById("animationHalfPeriodSeconds");
  const rowCountInput = document.getElementById("rowCount");
  const colCountInput = document.getElementById("colCount");
  const applyGridButton = document.getElementById("applyGrid");
  const showPinsInput = document.getElementById("showPins");
  const showCenterLinesInput = document.getElementById("showCenterLines");
  const pinRadiusRatioInput = document.getElementById("pinRadiusRatio");
  const plateThicknessInput = document.getElementById("plateThickness");
  const lockSelectedXYButton = document.getElementById("lockSelectedXY");
  const positionLockModeInput = document.getElementById("positionLockMode");
  const unlockSelectedXYButton = document.getElementById("unlockSelectedXY");
  const clearXYLocksButton = document.getElementById("clearXYLocks");
  const multiSelectModeInput = document.getElementById("multiSelectMode");
  const aTopOut = document.getElementById("aTopOut");
  const residualMetric = document.getElementById("residualMetric");
  const pitchMetric = document.getElementById("pitchMetric");
  const aCenterMetric = document.getElementById("aCenterMetric");
  const bCenterMetric = document.getElementById("bCenterMetric");
  const cCenterMetric = document.getElementById("cCenterMetric");
  const systemCenterMetric = document.getElementById("systemCenterMetric");
  const secondResidualMetric = document.getElementById("secondResidualMetric");
  const secondTargetMetric = document.getElementById("secondTargetMetric");
  const secondBMetric = document.getElementById("secondBMetric");
  const bcPrimaryResidualMetric = document.getElementById("bcPrimaryResidualMetric");
  const bcSecondResidualMetric = document.getElementById("bcSecondResidualMetric");
  const bcPitchMetric = document.getElementById("bcPitchMetric");
  const collisionState = document.getElementById("collisionState");
  const allowedDriveMetric = document.getElementById("allowedDriveMetric");
  const driveClampMetric = document.getElementById("driveClampMetric");
  const penetrationMetric = document.getElementById("penetrationMetric");
  const clearanceMetric = document.getElementById("clearanceMetric");
  const effectiveBLowerMetric = document.getElementById("effectiveBLowerMetric");
  const effectiveBTopMetric = document.getElementById("effectiveBTopMetric");
  const effectiveCLowerMetric = document.getElementById("effectiveCLowerMetric");
  const effectiveCTopMetric = document.getElementById("effectiveCTopMetric");
  const latticeSizeMetric = document.getElementById("latticeSizeMetric");
  const cellCountMetric = document.getElementById("cellCountMetric");
  const selectedCellMetric = document.getElementById("selectedCellMetric");
  const dieOffRadiusMetric = document.getElementById("dieOffRadiusMetric");
  const maxPropagatedMetric = document.getElementById("maxPropagatedMetric");
  const selectedCellsMetric = document.getElementById("selectedCellsMetric");
  const selectedLockMetric = document.getElementById("selectedLockMetric");
  const activeLocksMetric = document.getElementById("activeLocksMetric");
  const maxLiftMetric = document.getElementById("maxLiftMetric");
  const maxTiltMetric = document.getElementById("maxTiltMetric");
  const maxLockDriftMetric = document.getElementById("maxLockDriftMetric");
  const pinRadiusRatioOut = document.getElementById("pinRadiusRatioOut");
  const holeRadiusMetric = document.getElementById("holeRadiusMetric");
  const pinRadiusMetric = document.getElementById("pinRadiusMetric");
  const diametralBacklashMetric = document.getElementById("diametralBacklashMetric");
  const radialClearanceMetric = document.getElementById("radialClearanceMetric");
  const freeAngleMetric = document.getElementById("freeAngleMetric");
  const relativeBacklashMetric = document.getElementById("relativeBacklashMetric");
  const exactTiltMetric = document.getElementById("exactTiltMetric");
  const computedContactCommandMetric = document.getElementById("computedContactCommandMetric");
  const contactStateMetric = document.getElementById("contactStateMetric");
  const transmittedCommandMetric = document.getElementById("transmittedCommandMetric");

  if (!THREE || !mount) {
    if (mount) mount.textContent = "Three.js did not load.";
    return;
  }

  const CAD = Object.freeze({
    cellWidthMm: 55.604331,
    nominalHoleDiameterMm: 3.4,
    bodyThicknessMm: 4.0,
    padRadiusMm: 4.9,
    hubRadiusMm: 4.6,
    armWidthMm: 5.4,
    siteRadiusMm: 22.1,
  });
  const BODY_EPSILON_MM = 1e-4;
  const FEASIBLE_PIN_RESIDUAL_MM = 0.01;
  const DRIVE_SAMPLE_STEP_DEG = 1;
  const A_TOP_FULL_MIN_DEG = Number(aTopInput.min);
  const A_TOP_FULL_MAX_DEG = Number(aTopInput.max);
  const B_LOWER_REFERENCE_DEG = 0;
  const B_TOP_REFERENCE_DEG = -25;
  const B_LOWER_MIN_DEG = -90;
  const B_LOWER_MAX_DEG = 90;
  const B_TOP_MIN_DEG = -80;
  const B_TOP_MAX_DEG = 80;
  const LONG_SINGLE_STRAND_MIN_CELLS = 15;
  const LONG_SINGLE_STRAND_MIN_DRIVE_DEG = 33;
  const MIN_CONTINUOUS_CENTER_PITCH_RATIO = 0.75;
  const MAX_DRIVE_STEP_CENTER_JUMP_MM = 4;
  const MAX_DRIVE_STEP_ROTATION_JUMP_DEG = 5;
  const ANIMATION_EASE_BLEND = 0.18;

  const SITE_VECTORS = Object.freeze({
    east: new THREE.Vector2(CAD.siteRadiusMm, 0),
    north: new THREE.Vector2(0, CAD.siteRadiusMm),
    west: new THREE.Vector2(-CAD.siteRadiusMm, 0),
    south: new THREE.Vector2(0, -CAD.siteRadiusMm),
  });
  const FIXED_AB_SITES = Object.freeze({
    aUpper: "east",
    bLower: "south",
    bUpper: "south",
    aLower: "north",
  });

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(42, 1, 1, 6000);
  camera.up.set(0, 0, 1);

  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.shadowMap.enabled = true;
  mount.appendChild(renderer.domElement);

  const cameraState = {
    radius: 360,
    azimuth: -0.82,
    elevation: 0.58,
    target: new THREE.Vector3(0, 0, 0),
  };
  let driveEnvelopeCache = { key: "", value: null };
  let latestLatticeRadius = 180;
  let currentView = "iso";
  let lastMechanismUpdateKey = "";
  let lastControlStateKey = "";
  let renderedPinRadius = null;
  let latestContactField = null;
  let latestRestSnapshot = null;
  let latestBoundaryReport = null;
  let latestPinFitReport = null;
  let contactMemory = null;
  let freeTranslation = [0, 0, 0];
  let selectedDriveCell = { row: 0, col: 0 };
  let lastActuationMode = actuationModeInput.value;
  let selectedDriveCommandDirty = false;
  const selectedCellKeys = new Set(["0,0"]);
  const cellDriveCommands = new Map();
  const xyLocks = new Map();
  const pointer = new THREE.Vector2();
  const raycaster = new THREE.Raycaster();
  const selectableMeshes = [];

  const upperMaterial = new THREE.MeshStandardMaterial({ color: 0x8f969f, roughness: 0.52, metalness: 0.18 });
  const lowerMaterial = new THREE.MeshStandardMaterial({ color: 0x5d6570, roughness: 0.58, metalness: 0.16 });

  const materials = {
    upper: upperMaterial,
    lower: lowerMaterial,
    edge: new THREE.LineBasicMaterial({ color: 0x121820, transparent: true, opacity: 0.46 }),
    pin: new THREE.MeshStandardMaterial({ color: 0x1f2933, roughness: 0.36, metalness: 0.35 }),
    centerLine: new THREE.LineBasicMaterial({ color: 0x6b7480, transparent: true, opacity: 0.62 }),
    source: new THREE.LineBasicMaterial({ color: 0x22a36f, transparent: true, opacity: 0.95 }),
    selectionSet: new THREE.LineBasicMaterial({ color: 0x4bb6ff, transparent: true, opacity: 0.9 }),
    lock: new THREE.MeshStandardMaterial({ color: 0xd9273f, roughness: 0.44, metalness: 0.1 }),
  };

  const latticeRoot = new THREE.Group();
  const pinRoot = new THREE.Group();
  const guideRoot = new THREE.Group();
  const selectionRoot = new THREE.Group();
  scene.add(latticeRoot);
  scene.add(pinRoot);
  scene.add(guideRoot);
  scene.add(selectionRoot);

  let cellViews = [];
  let connections = [];
  let connectionPinMeshes = [];
  const latticeLine = new THREE.LineSegments(new THREE.BufferGeometry(), materials.centerLine);
  guideRoot.add(latticeLine);
  const selectionRing = new THREE.LineLoop(new THREE.BufferGeometry(), materials.source);
  selectionRoot.add(selectionRing);
  const selectionMarkerRoot = new THREE.Group();
  selectionRoot.add(selectionMarkerRoot);
  const lockMarkerRoot = new THREE.Group();
  selectionRoot.add(lockMarkerRoot);

  function rotate2(vector, angle) {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return new THREE.Vector2(c * vector.x - s * vector.y, s * vector.x + c * vector.y);
  }

  function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
  }

  function clampInt(value, min, max) {
    const parsed = Number.parseInt(value, 10);
    if (!Number.isFinite(parsed)) return min;
    return Math.max(min, Math.min(max, parsed));
  }

  function cleanZero(value) {
    return Math.abs(value) < 0.05 ? 0 : value;
  }

  function formatPair(point) {
    return `(${cleanZero(point.x).toFixed(1)}, ${cleanZero(point.y).toFixed(1)})`;
  }

  function selectedCellLabel() {
    const active = `r${selectedDriveCell.row + 1} c${selectedDriveCell.col + 1}`;
    return selectedCellKeys.size > 1 ? `${active} (+${selectedCellKeys.size - 1})` : active;
  }

  function selectedCellsLabel() {
    const labels = [...selectedCellKeys]
      .sort()
      .map((key) => {
        const [row, col] = key.split(",").map(Number);
        return `r${row + 1} c${col + 1}`;
      });
    if (!labels.length) return selectedCellLabel();
    if (labels.length <= 4) return labels.join(", ");
    return `${labels.slice(0, 4).join(", ")} +${labels.length - 4}`;
  }

  function holeRadiusMm() {
    return CAD.nominalHoleDiameterMm * 0.5;
  }

  function pinRadiusRatio() {
    if (actuationModeInput.value === "uniform") return 1;
    const value = Number(pinRadiusRatioInput.value);
    return clamp(Number.isFinite(value) ? value : 0.72, Number(pinRadiusRatioInput.min), Number(pinRadiusRatioInput.max));
  }

  function pinRadiusMm() {
    return holeRadiusMm() * pinRadiusRatio();
  }

  function pinDiameterMm() {
    return pinRadiusMm() * 2;
  }

  function holeDiameterMm() {
    return holeRadiusMm() * 2;
  }

  function diametralBacklashMm() {
    return Math.max(0, holeDiameterMm() - pinDiameterMm());
  }

  function pinHoleClearanceMm() {
    return diametralBacklashMm() * 0.5;
  }

  function plateThicknessMm() {
    const value = Number(plateThicknessInput.value);
    return clamp(Number.isFinite(value) ? value : CAD.bodyThicknessMm, Number(plateThicknessInput.min), Number(plateThicknessInput.max));
  }

  function normalizedBacklash() {
    return clamp(diametralBacklashMm() / CAD.siteRadiusMm, 0, 1);
  }

  function researchAngularBacklashRad() {
    return Math.asin(normalizedBacklash());
  }

  function researchAngularBacklashDeg() {
    return radToDeg(researchAngularBacklashRad());
  }

  function contactTransmit(commandDeg) {
    const backlash = researchAngularBacklashDeg();
    return Math.sign(commandDeg) * Math.max(0, Math.abs(commandDeg) - backlash);
  }

  function boundaryActuationScale() {
    if (xyLocks.size < 2) return 1;
    const isString = cellViews.length === 1 || cellViews[0]?.length === 1;
    if (isString && [...xyLocks.values()].filter((lock) => lock.lockZ !== false).length < 2) return 1;
    return exactPinTiltLimitRad() <= 1e-9 ? 0 : 1;
  }

  function jammedLockDriveDeg(envelope) {
    if (boundaryActuationScale() > 0) return null;
    const drives = [...xyLocks.values()]
      .map((lock) => Number(lock.targetDriveDeg))
      .filter(Number.isFinite)
      .sort((a, b) => a - b);
    const fallback = Number.isFinite(Number(aTopInput.min)) ? Number(aTopInput.min) : A_TOP_FULL_MIN_DEG;
    return clampDriveCommand(drives.length ? drives[Math.floor(drives.length * 0.5)] : fallback, envelope);
  }

  function exactPinTiltLimitRad() {
    const d = pinDiameterMm();
    const D = holeDiameterMm();
    const t = plateThicknessMm();
    if (d <= 0 || D <= 0 || t <= 0 || d > D) return 0;
    const ratio = clamp(d / Math.sqrt(D * D + t * t), -1, 1);
    return Math.max(0, Math.acos(ratio) - Math.atan(t / D));
  }

  function exactPinTiltLimitDeg() {
    return radToDeg(exactPinTiltLimitRad());
  }

  function verticalRiseLimitForSpan(spanMm, phiRad = exactPinTiltLimitRad()) {
    return Math.max(0, Math.tan(Math.max(0, phiRad)) * Math.max(0, spanMm));
  }

  function cellKey(row, col) {
    return `${row},${col}`;
  }

  function lockForCell(row, col) {
    return xyLocks.get(cellKey(row, col)) || null;
  }

  function clampDriveCommand(value, envelope) {
    if (envelope && !envelopeContainsDrive(value, envelope)) return nearestFeasibleDrive(value, envelope);
    return clamp(value, A_TOP_FULL_MIN_DEG, A_TOP_FULL_MAX_DEG);
  }

  function driveCommandForCell(row, col, envelope) {
    const key = cellKey(row, col);
    const stored = cellDriveCommands.has(key)
      ? cellDriveCommands.get(key)
      : latestContactField?.field?.[row]?.[col] ?? contactRestDriveDeg(envelope);
    return clampDriveCommand(stored, envelope);
  }

  function setDriveCommandForCell(row, col, value, envelope) {
    cellDriveCommands.set(cellKey(row, col), clampDriveCommand(value, envelope));
  }

  function selectedDriveCommand(envelope) {
    return driveCommandForCell(selectedDriveCell.row, selectedDriveCell.col, envelope);
  }

  function seedContactCommandsFromUniformDrive(angleDeg, envelope) {
    contactMemory = { restDrive: angleDeg, field: uniformDriveField(angleDeg), pins: [] };
    cellDriveCommands.clear();
    cellViews.forEach((rowViews) => {
      rowViews.forEach((cell) => {
        setDriveCommandForCell(cell.row, cell.col, angleDeg, envelope);
      });
    });
  }

  function selectedCells() {
    const cells = [...selectedCellKeys].map((key) => {
      const [row, col] = key.split(",").map(Number);
      return cellAt(row, col);
    });
    const filtered = cells.filter(Boolean);
    if (filtered.length) return filtered;
    const active = cellAt(selectedDriveCell.row, selectedDriveCell.col);
    return active ? [active] : [];
  }

  function contactRestDriveDeg(envelope) {
    if (contactMemory) return contactMemory.restDrive;
    if (cellDriveCommands.size) {
      const commands = [...cellDriveCommands.values()].map((value) => clampDriveCommand(value, envelope)).sort((a, b) => a - b);
      return commands[Math.floor(commands.length * 0.5)];
    }
    return envelope ? envelope.max : A_TOP_FULL_MAX_DEG;
  }

  function neighborCells(row, col, rows, cols) {
    const cells = [];
    if (row > 0) cells.push([row - 1, col]);
    if (row + 1 < rows) cells.push([row + 1, col]);
    if (col > 0) cells.push([row, col - 1]);
    if (col + 1 < cols) cells.push([row, col + 1]);
    return cells;
  }

  function contactDriveField(envelope) {
    const rows = cellViews.length;
    const cols = rows ? cellViews[0].length : 0;
    const restDrive = contactRestDriveDeg(envelope);
    let field = Array.from({ length: rows }, () => Array.from({ length: cols }, () => restDrive));
    const deltaField = Array.from({ length: rows }, () => Array.from({ length: cols }, () => 0));
    const distance = Array.from({ length: rows }, () => Array.from({ length: cols }, () => Infinity));
    const sourceRow = clampInt(selectedDriveCell.row, 0, Math.max(0, rows - 1));
    const sourceCol = clampInt(selectedDriveCell.col, 0, Math.max(0, cols - 1));
    const lockScale = boundaryActuationScale();
    if (lockScale <= 1e-9 && xyLocks.size >= 2) {
      return { field, deltaField, distance, reach: 0, dieOffRadius: 0, maxAbs: 0, sourceRow, sourceCol, restDrive, constraintScale: 0 };
    }
    const previous = contactMemory?.field?.length === rows && contactMemory.field[0]?.length === cols
      ? contactMemory.field : field;
    const command = selectedDriveCommand(envelope);
    field = Math.abs(command - previous[sourceRow][sourceCol]) < 1e-10
      ? previous.map((row) => row.slice())
      : window.RADBacklashState.advance(previous, sourceRow, sourceCol, command, researchAngularBacklashDeg());
    field.forEach((row, r) => row.forEach((value, c) => {
      deltaField[r][c] = value - restDrive;
      distance[r][c] = Math.abs(r - sourceRow) + Math.abs(c - sourceCol);
    }));

    let reach = 0;
    let dieOffRadius = 0;
    let maxAbs = 0;
    field.forEach((row, rowIndex) => {
      row.forEach((value, colIndex) => {
        if (Math.abs(value - restDrive) > 1e-6) reach += 1;
        if (Math.abs(value - restDrive) > 1e-6 && Number.isFinite(distance[rowIndex][colIndex])) {
          dieOffRadius = Math.max(dieOffRadius, distance[rowIndex][colIndex]);
        }
        maxAbs = Math.max(maxAbs, Math.abs(value - restDrive));
      });
    });
    return { field, deltaField, distance, reach, dieOffRadius, maxAbs, sourceRow, sourceCol, restDrive };
  }

  function scaledContactField(contactField, scale) {
    if (!contactField) return null;
    const clampedScale = clamp(scale, 0, 1);
    const field = contactField.field.map((row) => row.map((value) => contactField.restDrive + (value - contactField.restDrive) * clampedScale));
    const deltaField = contactField.deltaField.map((row) => row.map((value) => value * clampedScale));
    let reach = 0;
    let dieOffRadius = 0;
    let maxAbs = 0;
    field.forEach((row, rowIndex) => {
      row.forEach((value, colIndex) => {
        const delta = value - contactField.restDrive;
        if (Math.abs(delta) > 1e-6) {
          reach += 1;
          if (Number.isFinite(contactField.distance[rowIndex][colIndex])) {
            dieOffRadius = Math.max(dieOffRadius, contactField.distance[rowIndex][colIndex]);
          }
        }
        maxAbs = Math.max(maxAbs, Math.abs(delta));
      });
    });
    return {
      ...contactField,
      field,
      deltaField,
      reach,
      dieOffRadius,
      maxAbs,
      constraintScale: clampedScale,
    };
  }

  function degToRad(value) {
    return (value * Math.PI) / 180;
  }

  function radToDeg(value) {
    return (value * 180) / Math.PI;
  }

  function displayFrameForPose(pose) {
    const solvedAngle = Math.atan2(pose.bCenter.y, pose.bCenter.x);
    const safeSolvedAngle = Number.isFinite(solvedAngle) ? solvedAngle : 0;
    const displayRotation = -safeSolvedAngle;
    const displayedBVector = rotate2(pose.bCenter, displayRotation);
    const shift = displayedBVector.clone().multiplyScalar(-1);
    return {
      rotation: displayRotation,
      centerlineAngle: 0,
      mapPoint(point) {
        return rotate2(point, displayRotation).add(shift);
      },
    };
  }

  function pointSegmentDistance(point, start, end) {
    const segment = end.clone().sub(start);
    const lengthSquared = segment.lengthSq();
    if (lengthSquared <= 1e-12) return point.distanceTo(start);
    const t = clamp(point.clone().sub(start).dot(segment) / lengthSquared, 0, 1);
    return point.distanceTo(start.clone().add(segment.multiplyScalar(t)));
  }

  function orientation(a, b, c) {
    return Math.sign((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
  }

  function onSegment(a, b, c) {
    return (
      Math.min(a.x, b.x) - 1e-9 <= c.x &&
      c.x <= Math.max(a.x, b.x) + 1e-9 &&
      Math.min(a.y, b.y) - 1e-9 <= c.y &&
      c.y <= Math.max(a.y, b.y) + 1e-9
    );
  }

  function segmentsIntersect(a, b, c, d) {
    const o1 = orientation(a, b, c);
    const o2 = orientation(a, b, d);
    const o3 = orientation(c, d, a);
    const o4 = orientation(c, d, b);
    if (o1 !== o2 && o3 !== o4) return true;
    if (o1 === 0 && onSegment(a, b, c)) return true;
    if (o2 === 0 && onSegment(a, b, d)) return true;
    if (o3 === 0 && onSegment(c, d, a)) return true;
    return o4 === 0 && onSegment(c, d, b);
  }

  function segmentSegmentDistance(a, b, c, d) {
    if (segmentsIntersect(a, b, c, d)) return 0;
    return Math.min(
      pointSegmentDistance(a, c, d),
      pointSegmentDistance(b, c, d),
      pointSegmentDistance(c, a, b),
      pointSegmentDistance(d, a, b)
    );
  }

  function primitiveDistance(first, second) {
    if (first.type === "disk" && second.type === "disk") {
      return first.center.distanceTo(second.center) - first.radius - second.radius;
    }
    if (first.type === "capsule" && second.type === "capsule") {
      return segmentSegmentDistance(first.start, first.end, second.start, second.end) - first.radius - second.radius;
    }
    const capsule = first.type === "capsule" ? first : second;
    const disk = first.type === "disk" ? first : second;
    return pointSegmentDistance(disk.center, capsule.start, capsule.end) - disk.radius - capsule.radius;
  }

  function primitiveBounds(primitive) {
    if (primitive.type === "disk") {
      return {
        minX: primitive.center.x - primitive.radius,
        maxX: primitive.center.x + primitive.radius,
        minY: primitive.center.y - primitive.radius,
        maxY: primitive.center.y + primitive.radius,
      };
    }
    return {
      minX: Math.min(primitive.start.x, primitive.end.x) - primitive.radius,
      maxX: Math.max(primitive.start.x, primitive.end.x) + primitive.radius,
      minY: Math.min(primitive.start.y, primitive.end.y) - primitive.radius,
      maxY: Math.max(primitive.start.y, primitive.end.y) + primitive.radius,
    };
  }

  function combineBounds(bounds) {
    return bounds.reduce(
      (combined, next) => ({
        minX: Math.min(combined.minX, next.minX),
        maxX: Math.max(combined.maxX, next.maxX),
        minY: Math.min(combined.minY, next.minY),
        maxY: Math.max(combined.maxY, next.maxY),
      }),
      { minX: Infinity, maxX: -Infinity, minY: Infinity, maxY: -Infinity }
    );
  }

  function boundsOverlap(first, second) {
    return first.minX <= second.maxX && first.maxX >= second.minX && first.minY <= second.maxY && first.maxY >= second.minY;
  }

  function zBandsOverlap(first, second) {
    return Math.max(first.zMin, second.zMin) < Math.min(first.zMax, second.zMax) - BODY_EPSILON_MM;
  }

  function localBodyPrimitives() {
    const armRadius = CAD.armWidthMm * 0.5;
    const r = CAD.siteRadiusMm;
    return [
      { type: "capsule", label: "horizontal arm", start: new THREE.Vector2(-r, 0), end: new THREE.Vector2(r, 0), radius: armRadius },
      { type: "capsule", label: "vertical arm", start: new THREE.Vector2(0, -r), end: new THREE.Vector2(0, r), radius: armRadius },
      { type: "disk", label: "center hub", center: new THREE.Vector2(0, 0), radius: CAD.hubRadiusMm },
      { type: "disk", label: "east pad", center: SITE_VECTORS.east, radius: CAD.padRadiusMm },
      { type: "disk", label: "north pad", center: SITE_VECTORS.north, radius: CAD.padRadiusMm },
      { type: "disk", label: "west pad", center: SITE_VECTORS.west, radius: CAD.padRadiusMm },
      { type: "disk", label: "south pad", center: SITE_VECTORS.south, radius: CAD.padRadiusMm },
    ];
  }

  const LOCAL_BODY_PRIMITIVES = localBodyPrimitives();

  function transformPoint(point, origin, angle) {
    return rotate2(point, angle).add(origin);
  }

  function bodyPrimitives(body) {
    return LOCAL_BODY_PRIMITIVES.map((primitive) => {
      const common = {
        part: body.part,
        cell: body.cell,
        zMin: body.zMin,
        zMax: body.zMax,
        label: primitive.label,
      };
      if (primitive.type === "disk") {
        return {
          ...common,
          type: "disk",
          center: transformPoint(primitive.center, body.origin, body.rotation),
          radius: primitive.radius,
        };
      }
      return {
        ...common,
        type: "capsule",
        start: transformPoint(primitive.start, body.origin, body.rotation),
        end: transformPoint(primitive.end, body.origin, body.rotation),
        radius: primitive.radius,
      };
    });
  }

  function collisionReportForBodies(rawBodies) {
    const bodies = rawBodies.map((body) => {
      const primitives = bodyPrimitives(body);
      return { ...body, primitives, bounds: combineBounds(primitives.map(primitiveBounds)) };
    });
    let maxPenetration = 0;
    let minClearance = Infinity;
    const collisions = [];
    for (let i = 0; i < bodies.length; i += 1) {
      for (let j = i + 1; j < bodies.length; j += 1) {
        const firstBody = bodies[i];
        const secondBody = bodies[j];
        if (firstBody.cell === secondBody.cell || !zBandsOverlap(firstBody, secondBody) || !boundsOverlap(firstBody.bounds, secondBody.bounds)) continue;
        for (const first of firstBody.primitives) {
          for (const second of secondBody.primitives) {
            const clearance = primitiveDistance(first, second);
            minClearance = Math.min(minClearance, clearance);
            if (clearance < 0) {
              const penetration = -clearance;
              maxPenetration = Math.max(maxPenetration, penetration);
              collisions.push({
                first: `${first.part} ${first.label}`,
                second: `${second.part} ${second.label}`,
                penetration,
              });
            }
          }
        }
      }
    }
    return {
      maxPenetration,
      minClearance: Number.isFinite(minClearance) ? minClearance : 0,
      collisionCount: collisions.length,
      collisions,
      clear: maxPenetration <= BODY_EPSILON_MM,
    };
  }

  function normalizeNear(angle, reference) {
    let normalized = angle;
    while (normalized - reference > Math.PI) normalized -= Math.PI * 2;
    while (normalized - reference < -Math.PI) normalized += Math.PI * 2;
    return normalized;
  }

  function angleDelta(first, second) {
    return Math.atan2(Math.sin(first - second), Math.cos(first - second));
  }

  function oppositeVector(vector) {
    return vector.clone().multiplyScalar(-1);
  }

  function buildPose(aTopRad, bLowerRad, bTopRad) {
    const aAttachLocal = SITE_VECTORS[FIXED_AB_SITES.aUpper];
    const bAttachLocal = SITE_VECTORS[FIXED_AB_SITES.bLower];
    const aLowerAttachLocal = SITE_VECTORS[FIXED_AB_SITES.aLower];
    const bTopAttachLocal = SITE_VECTORS[FIXED_AB_SITES.bUpper];
    const aAttachWorld = rotate2(aAttachLocal, aTopRad);
    const bAttachOffset = rotate2(bAttachLocal, bLowerRad);
    const bCenter = aAttachWorld.clone().sub(bAttachOffset);
    const bTopWorldRad = bLowerRad + bTopRad;
    const aLowerAttachWorld = aLowerAttachLocal.clone();
    const bTopAttachOffset = rotate2(bTopAttachLocal, bTopWorldRad);
    const bTopAttachWorld = bCenter.clone().add(bTopAttachOffset);
    const cCenter = bCenter.clone().multiplyScalar(2);
    const cLowerRad = 0;
    const cTopWorldRad = aTopRad;
    const bcPrimaryBLocal = oppositeVector(bAttachLocal);
    const bcPrimaryCLocal = oppositeVector(aAttachLocal);
    const bcSecondBLocal = oppositeVector(bTopAttachLocal);
    const bcSecondCLocal = oppositeVector(aLowerAttachLocal);
    const bcPrimaryBWorld = bCenter.clone().add(rotate2(bcPrimaryBLocal, bLowerRad));
    const bcPrimaryCWorld = cCenter.clone().add(rotate2(bcPrimaryCLocal, cTopWorldRad));
    const bcSecondBWorld = bCenter.clone().add(rotate2(bcSecondBLocal, bTopWorldRad));
    const bcSecondCWorld = cCenter.clone().add(rotate2(bcSecondCLocal, cLowerRad));
    return {
      aTopRad,
      bLowerRad,
      bTopRad,
      bTopWorldRad,
      cLowerRad,
      cTopWorldRad,
      aAttachWorld,
      bAttachOffset,
      bCenter,
      cCenter,
      aLowerAttachWorld,
      bTopAttachWorld,
      bcPrimaryBWorld,
      bcPrimaryCWorld,
      bcSecondBWorld,
      bcSecondCWorld,
      bodies: [
        { cell: "A", part: "A lower", origin: new THREE.Vector2(0, 0), rotation: 0, zMin: -CAD.bodyThicknessMm, zMax: 0 },
        { cell: "A", part: "A upper", origin: new THREE.Vector2(0, 0), rotation: aTopRad, zMin: 0, zMax: CAD.bodyThicknessMm },
        { cell: "B", part: "B lower", origin: bCenter, rotation: bLowerRad, zMin: -CAD.bodyThicknessMm, zMax: 0 },
        { cell: "B", part: "B upper", origin: bCenter, rotation: bTopWorldRad, zMin: 0, zMax: CAD.bodyThicknessMm },
        { cell: "C", part: "C lower", origin: cCenter, rotation: cLowerRad, zMin: -CAD.bodyThicknessMm, zMax: 0 },
        { cell: "C", part: "C upper", origin: cCenter, rotation: cTopWorldRad, zMin: 0, zMax: CAD.bodyThicknessMm },
      ],
    };
  }

  function collisionReport(pose) {
    return collisionReportForBodies(pose.bodies);
  }

  function pinResiduals(pose) {
    const abPrimary = pose.aAttachWorld.distanceTo(pose.bCenter.clone().add(pose.bAttachOffset));
    const abSecond = pose.aLowerAttachWorld.distanceTo(pose.bTopAttachWorld);
    const bcPrimary = pose.bcPrimaryBWorld.distanceTo(pose.bcPrimaryCWorld);
    const bcSecond = pose.bcSecondBWorld.distanceTo(pose.bcSecondCWorld);
    return {
      abPrimary,
      abSecond,
      bcPrimary,
      bcSecond,
      max: Math.max(abPrimary, abSecond, bcPrimary, bcSecond),
    };
  }

  function exactSecondPinCandidates(aTopRad, requestedBLowerRad, requestedBTopRad) {
    const aAttachLocal = SITE_VECTORS[FIXED_AB_SITES.aUpper];
    const bAttachLocal = SITE_VECTORS[FIXED_AB_SITES.bLower];
    const aLowerAttachLocal = SITE_VECTORS[FIXED_AB_SITES.aLower];
    const bTopAttachLocal = SITE_VECTORS[FIXED_AB_SITES.bUpper];
    const primaryWorld = rotate2(aAttachLocal, aTopRad);
    const targetWorld = aLowerAttachLocal.clone();
    const targetFromPrimary = targetWorld.clone().sub(primaryWorld);
    const targetDistance = targetFromPrimary.length();
    const lowerArmRadius = bAttachLocal.length();
    const topArmRadius = bTopAttachLocal.length();
    if (targetDistance < 1e-9 || lowerArmRadius < 1e-9 || topArmRadius < 1e-9) return [];
    if (Math.abs(lowerArmRadius - topArmRadius) > 1e-6 || targetDistance > lowerArmRadius + topArmRadius + 1e-6) return [];

    const lowerMin = degToRad(B_LOWER_MIN_DEG);
    const lowerMax = degToRad(B_LOWER_MAX_DEG);
    const topMin = degToRad(B_TOP_MIN_DEG);
    const topMax = degToRad(B_TOP_MAX_DEG);
    const targetAngle = Math.atan2(targetFromPrimary.y, targetFromPrimary.x);
    const lowerLocalAngle = Math.atan2(bAttachLocal.y, bAttachLocal.x);
    const topLocalAngle = Math.atan2(bTopAttachLocal.y, bTopAttachLocal.x);
    const triangleAngle = Math.acos(clamp(-targetDistance / (2 * lowerArmRadius), -1, 1));
    const candidates = [];

    [-1, 1].forEach((sign) => {
      const bLowerRad = normalizeNear(targetAngle + sign * triangleAngle - lowerLocalAngle, requestedBLowerRad);
      const bAttachOffset = rotate2(bAttachLocal, bLowerRad);
      const bCenter = primaryWorld.clone().sub(bAttachOffset);
      const targetFromB = targetWorld.clone().sub(bCenter);
      const bTopWorldRad = Math.atan2(targetFromB.y, targetFromB.x) - topLocalAngle;
      const bTopRad = normalizeNear(bTopWorldRad - bLowerRad, requestedBTopRad);
      if (bLowerRad < lowerMin - 1e-8 || bLowerRad > lowerMax + 1e-8) return;
      if (bTopRad < topMin - 1e-8 || bTopRad > topMax + 1e-8) return;
      candidates.push(buildPose(aTopRad, bLowerRad, bTopRad));
    });

    return candidates;
  }

  function evaluateExactLoopClosure(aTopRad, requestedBLowerRad, requestedBTopRad) {
    let best = null;
    exactSecondPinCandidates(aTopRad, requestedBLowerRad, requestedBTopRad).forEach((pose) => {
      const report = collisionReport(pose);
      const residual = pinResiduals(pose).max;
      const angleScore = Math.abs(angleDelta(pose.bLowerRad, requestedBLowerRad)) + 0.65 * Math.abs(angleDelta(pose.bTopRad, requestedBTopRad));
      const score = angleScore + report.maxPenetration * 10;
      const continuous = isContinuousAssemblyPose(pose);
      const candidate = { pose, report, residual, score, angleScore, continuous };
      const candidateFeasible = continuous && report.clear && residual <= FEASIBLE_PIN_RESIDUAL_MM;
      const bestFeasible = best && best.continuous && best.report.clear && best.residual <= FEASIBLE_PIN_RESIDUAL_MM;
      if (
        !best ||
        (candidateFeasible && !bestFeasible) ||
        (candidateFeasible === bestFeasible &&
          (residual < best.residual - 1e-9 ||
            (Math.abs(residual - best.residual) <= 1e-9 && report.maxPenetration < best.report.maxPenetration - 1e-9) ||
            (Math.abs(residual - best.residual) <= 1e-9 &&
              Math.abs(report.maxPenetration - best.report.maxPenetration) <= 1e-9 &&
              score < best.score)))
      ) {
        best = candidate;
      }
    });

    if (!best) {
      const fallbackPose = buildPose(aTopRad, requestedBLowerRad, requestedBTopRad);
      const fallbackReport = collisionReport(fallbackPose);
      return {
        pose: fallbackPose,
        report: fallbackReport,
        residual: pinResiduals(fallbackPose).max,
        score: Infinity,
        angleScore: Infinity,
        status: "blocked",
      };
    }

    const attached = best.residual <= FEASIBLE_PIN_RESIDUAL_MM;
    const clear = best.report.clear;
    return {
      ...best,
      status: best.continuous && clear && attached ? (best.angleScore <= 1e-9 ? "clear" : "adjusted") : "blocked",
    };
  }

  function isContinuousAssemblyPose(pose) {
    return pose.bCenter.length() >= CAD.siteRadiusMm * MIN_CONTINUOUS_CENTER_PITCH_RATIO;
  }

  function solveNonPenetratingPose(aTopRad, requestedBLowerRad, requestedBTopRad) {
    return evaluateExactLoopClosure(aTopRad, requestedBLowerRad, requestedBTopRad);
  }

  function isDriveAngleFeasible(angleDeg, requestedBLowerRad, requestedBTopRad) {
    const solution = evaluateExactLoopClosure(degToRad(angleDeg), requestedBLowerRad, requestedBTopRad);
    if (solution.status === "blocked") return false;
    return latticeCollisionReportForDrive(angleDeg, requestedBLowerRad, requestedBTopRad).clear;
  }

  function latticeMinimumDriveDeg() {
    const rows = cellViews.length;
    const cols = rows ? cellViews[0].length : 0;
    const isSingleWidthStrand = Math.min(rows, cols) === 1;
    const strandLength = Math.max(rows, cols);
    if (isSingleWidthStrand && strandLength >= LONG_SINGLE_STRAND_MIN_CELLS) {
      return LONG_SINGLE_STRAND_MIN_DRIVE_DEG;
    }
    return A_TOP_FULL_MIN_DEG;
  }

  function snapshotEntries(snapshot) {
    return snapshot.flatMap((row) =>
      row.map((entry) => ({
        center: entry.center,
        topRotation: entry.topRotation,
        bottomRotation: entry.bottomRotation,
      }))
    );
  }

  function maxSnapshotEntryJump(previous, current) {
    let centerJump = 0;
    let rotationJump = 0;
    for (let index = 0; index < Math.min(previous.length, current.length); index += 1) {
      centerJump = Math.max(centerJump, previous[index].center.distanceTo(current[index].center));
      rotationJump = Math.max(
        rotationJump,
        Math.abs(angleDelta(previous[index].topRotation, current[index].topRotation)),
        Math.abs(angleDelta(previous[index].bottomRotation, current[index].bottomRotation))
      );
    }
    return { centerJump, rotationJumpDeg: radToDeg(rotationJump) };
  }

  function continuousDriveSamples(feasible, requestedBLowerRad, requestedBTopRad) {
    const continuous = [];
    let previous = null;
    feasible.forEach((angleDeg) => {
      const snapshot = snapshotEntries(buildLatticeSnapshot(uniformDriveField(angleDeg), requestedBLowerRad, requestedBTopRad));
      if (previous && angleDeg - previous.angleDeg <= DRIVE_SAMPLE_STEP_DEG + 1e-9) {
        const jump = maxSnapshotEntryJump(previous.snapshot, snapshot);
        if (jump.centerJump > MAX_DRIVE_STEP_CENTER_JUMP_MM || jump.rotationJumpDeg > MAX_DRIVE_STEP_ROTATION_JUMP_DEG) {
          previous = null;
          return;
        }
      }
      continuous.push(angleDeg);
      previous = { angleDeg, snapshot };
    });
    return continuous;
  }

  function computeDriveEnvelope(requestedBLowerRad, requestedBTopRad) {
    const rawFeasible = [];
    const lowerBound = latticeMinimumDriveDeg();
    for (let angleDeg = A_TOP_FULL_MIN_DEG; angleDeg <= A_TOP_FULL_MAX_DEG; angleDeg += DRIVE_SAMPLE_STEP_DEG) {
      if (angleDeg < lowerBound) continue;
      if (isDriveAngleFeasible(angleDeg, requestedBLowerRad, requestedBTopRad)) rawFeasible.push(angleDeg);
    }
    const feasible = continuousDriveSamples(rawFeasible, requestedBLowerRad, requestedBTopRad);
    if (!feasible.length) return null;

    const ranges = [];
    let start = feasible[0];
    let previous = feasible[0];
    for (let index = 1; index < feasible.length; index += 1) {
      const value = feasible[index];
      if (value - previous > DRIVE_SAMPLE_STEP_DEG + 1e-9) {
        ranges.push([start, previous]);
        start = value;
      }
      previous = value;
    }
    ranges.push([start, previous]);

    return {
      feasible,
      ranges,
      min: feasible[0],
      max: feasible[feasible.length - 1],
    };
  }

  function cachedDriveEnvelope(requestedBLowerRad, requestedBTopRad) {
    const rows = cellViews.length;
    const cols = rows ? cellViews[0].length : 0;
    const key = [
      FIXED_AB_SITES.aUpper,
      FIXED_AB_SITES.bLower,
      FIXED_AB_SITES.bUpper,
      FIXED_AB_SITES.aLower,
      rows,
      cols,
      requestedBLowerRad.toFixed(6),
      requestedBTopRad.toFixed(6),
    ].join("|");
    if (driveEnvelopeCache.key !== key) {
      driveEnvelopeCache = { key, value: computeDriveEnvelope(requestedBLowerRad, requestedBTopRad) };
    }
    return driveEnvelopeCache.value;
  }

  function nearestFeasibleDrive(requestedDeg, envelope) {
    let best = envelope.feasible[0];
    let bestDistance = Math.abs(requestedDeg - best);
    envelope.feasible.forEach((value) => {
      const distance = Math.abs(requestedDeg - value);
      if (distance < bestDistance) {
        best = value;
        bestDistance = distance;
      }
    });
    return best;
  }

  function envelopeContainsDrive(angleDeg, envelope) {
    if (!envelope) return false;
    return envelope.ranges.some((range) => angleDeg >= range[0] - 1e-9 && angleDeg <= range[1] + 1e-9);
  }

  function widestRange(envelope) {
    return envelope.ranges.reduce((best, range) => (range[1] - range[0] > best[1] - best[0] ? range : best), envelope.ranges[0]);
  }

  function mildlyEasedPingPong(timeMs, minDeg, maxDeg, halfPeriodSeconds, easeBlend) {
    const span = Math.max(0, maxDeg - minDeg);
    const duration = Math.max(0.25, halfPeriodSeconds);
    if (span <= 1e-9) return minDeg;
    const phase = (timeMs * 0.001) / duration;
    const cycle = phase % 2;
    const linearT = cycle <= 1 ? cycle : 2 - cycle;
    const sineEaseT = 0.5 - 0.5 * Math.cos(Math.PI * linearT);
    const t = linearT + clamp(easeBlend, 0, 0.35) * (sineEaseT - linearT);
    return minDeg + span * t;
  }

  function formatDriveRanges(envelope) {
    if (!envelope) return "none";
    return envelope.ranges
      .map((range) => (range[0] === range[1] ? `${range[0].toFixed(0)} deg` : `${range[0].toFixed(0)} to ${range[1].toFixed(0)} deg`))
      .join("; ");
  }

  function cylinderZ(radius, depth, material, segments = 48) {
    const mesh = new THREE.Mesh(new THREE.CylinderGeometry(radius, radius, depth, segments), material);
    mesh.rotation.x = Math.PI / 2;
    mesh.castShadow = true;
    mesh.receiveShadow = true;
    return mesh;
  }

  function replaceCylinderGeometry(mesh, radius, depth, segments = 36) {
    mesh.geometry.dispose();
    mesh.geometry = new THREE.CylinderGeometry(radius, radius, depth, segments);
  }

  function updatePinGeometry() {
    const radius = pinRadiusMm();
    if (renderedPinRadius !== null && Math.abs(renderedPinRadius - radius) < 1e-6) return;
    renderedPinRadius = radius;
    cellViews.forEach((row) => {
      row.forEach((cell) => {
        replaceCylinderGeometry(cell.centerPin, radius, CAD.bodyThicknessMm * 2.55, 36);
      });
    });
    connectionPinMeshes.forEach((pin) => {
      replaceCylinderGeometry(pin, radius, CAD.bodyThicknessMm * 2.95, 36);
    });
  }

  function annularDisk(outerRadius, innerRadius, depth, material, segments = 48) {
    const shape = new THREE.Shape();
    shape.absarc(0, 0, outerRadius, 0, Math.PI * 2, false);
    addCircularHole(shape, 0, 0, innerRadius);
    const geometry = new THREE.ExtrudeGeometry(shape, {
      depth,
      steps: 1,
      bevelEnabled: false,
      curveSegments: segments,
    });
    geometry.translate(0, 0, -depth * 0.5);
    const mesh = new THREE.Mesh(geometry, material);
    mesh.castShadow = true;
    mesh.receiveShadow = true;
    return mesh;
  }

  function addCircularHole(shape, x, y, radius) {
    const hole = new THREE.Path();
    hole.absarc(x, y, radius, 0, Math.PI * 2, true);
    shape.holes.push(hole);
  }

  function extrudedShapeMesh(shape, depth, material, curveSegments = 32) {
    const geometry = new THREE.ExtrudeGeometry(shape, {
      depth,
      steps: 1,
      bevelEnabled: false,
      curveSegments,
    });
    geometry.translate(0, 0, -depth * 0.5);
    const mesh = new THREE.Mesh(geometry, material);
    mesh.castShadow = true;
    mesh.receiveShadow = true;
    return mesh;
  }

  function addEdges(parent, mesh) {
    const edges = new THREE.LineSegments(new THREE.EdgesGeometry(mesh.geometry, 24), materials.edge);
    edges.position.copy(mesh.position);
    edges.rotation.copy(mesh.rotation);
    edges.scale.copy(mesh.scale);
    parent.add(edges);
  }

  function addArm(parent, rotationZ, material) {
    const halfWidth = CAD.armWidthMm * 0.5;
    const holeRadius = holeRadiusMm();
    const halfLength = CAD.siteRadiusMm + holeRadius + 0.15;
    const shape = new THREE.Shape();
    shape.moveTo(-halfLength, -halfWidth);
    shape.lineTo(halfLength, -halfWidth);
    shape.lineTo(halfLength, halfWidth);
    shape.lineTo(-halfLength, halfWidth);
    shape.lineTo(-halfLength, -halfWidth);
    addCircularHole(shape, -CAD.siteRadiusMm, 0, holeRadius);
    addCircularHole(shape, 0, 0, holeRadius);
    addCircularHole(shape, CAD.siteRadiusMm, 0, holeRadius);
    const arm = extrudedShapeMesh(shape, CAD.bodyThicknessMm, material, 32);
    arm.rotation.z = rotationZ;
    parent.add(arm);
    addEdges(parent, arm);
  }

  function addPad(parent, x, y, material) {
    const pad = annularDisk(CAD.padRadiusMm, holeRadiusMm(), CAD.bodyThicknessMm, material, 44);
    pad.position.set(x, y, 0);
    parent.add(pad);
    addEdges(parent, pad);
  }

  function createCrossPart(name, material) {
    const group = new THREE.Group();
    group.name = name;
    addArm(group, 0, material);
    addArm(group, Math.PI / 2, material);
    const hub = annularDisk(CAD.hubRadiusMm, holeRadiusMm(), CAD.bodyThicknessMm, material, 48);
    group.add(hub);
    addEdges(group, hub);
    addPad(group, CAD.siteRadiusMm, 0, material);
    addPad(group, -CAD.siteRadiusMm, 0, material);
    addPad(group, 0, CAD.siteRadiusMm, material);
    addPad(group, 0, -CAD.siteRadiusMm, material);
    return group;
  }

  function materialPairForCell() {
    return { top: materials.upper.clone(), bottom: materials.lower.clone() };
  }

  function driveColor(angleDeg, envelope) {
    const min = envelope ? envelope.min : A_TOP_FULL_MIN_DEG;
    const max = envelope ? envelope.max : A_TOP_FULL_MAX_DEG;
    const t = clamp((angleDeg - min) / Math.max(1e-6, max - min), 0, 1);
    return new THREE.Color(0xe60018).lerp(new THREE.Color(0x1478ff), t);
  }

  function applyCellColor(cell, angleDeg, envelope) {
    const color = driveColor(angleDeg, envelope);
    cell.topMaterial.color.copy(color);
    cell.bottomMaterial.color.copy(color).multiplyScalar(0.72);
  }

  function cellName(row, col) {
    if (row === 0 && col === 0) return "Cell A";
    if (row === 0 && col === 1) return "Cell B";
    if (row === 0 && col === 2) return "Cell C";
    return `Cell R${row + 1}C${col + 1}`;
  }

  function createCell(row, col) {
    const group = new THREE.Group();
    const name = cellName(row, col);
    const pair = materialPairForCell(row, col);
    group.name = name;
    const bottom = createCrossPart(`${name} lower cross`, pair.bottom);
    const top = createCrossPart(`${name} upper cross`, pair.top);
    const centerPin = cylinderZ(pinRadiusMm(), CAD.bodyThicknessMm * 2.55, materials.pin, 36);
    bottom.position.z = -CAD.bodyThicknessMm * 0.5;
    top.position.z = CAD.bodyThicknessMm * 0.5;
    group.add(bottom);
    group.add(top);
    group.add(centerPin);
    group.traverse((object) => {
      if (object.isMesh) {
        object.userData.cellRow = row;
        object.userData.cellCol = col;
        object.userData.cellSelectable = true;
        selectableMeshes.push(object);
      }
    });
    return {
      row,
      col,
      parity: (row + col) % 2,
      group,
      top,
      bottom,
      centerPin,
      topMaterial: pair.top,
      bottomMaterial: pair.bottom,
      lockMarker: null,
      center: new THREE.Vector2(),
      z: 0,
      tiltX: 0,
      tiltY: 0,
      topRotation: 0,
      bottomRotation: 0,
    };
  }

  function createLockMarker() {
    const geometry = new THREE.TorusGeometry(CAD.siteRadiusMm * 1.08, 0.85, 8, 64);
    const marker = new THREE.Mesh(geometry, materials.lock);
    marker.rotation.x = Math.PI / 2;
    marker.visible = false;
    lockMarkerRoot.add(marker);
    return marker;
  }

  function removeChildren(group) {
    while (group.children.length) {
      group.remove(group.children[group.children.length - 1]);
    }
  }

  function createConnectionPin() {
    const pin = cylinderZ(pinRadiusMm(), CAD.bodyThicknessMm * 2.95, materials.pin, 36);
    pinRoot.add(pin);
    return pin;
  }

  function ensureConnectionPinMeshCount(count) {
    while (connectionPinMeshes.length < count) {
      connectionPinMeshes.push(createConnectionPin());
    }
    for (let index = count; index < connectionPinMeshes.length; index += 1) {
      connectionPinMeshes[index].visible = false;
    }
  }

  function resizeFloorGridForDimensions(rows, cols) {
    const widestCellCount = Math.max(rows, cols);
    const desiredHalfSpan = widestCellCount * CAD.cellWidthMm * 0.58 + CAD.cellWidthMm;
    floorGrid.scale.setScalar(Math.max(1, desiredHalfSpan / 130));
  }

  function rebuildLattice() {
    contactMemory = null;
    freeTranslation = [0, 0, 0];
    const rows = clampInt(rowCountInput.value, Number(rowCountInput.min), Number(rowCountInput.max));
    const cols = clampInt(colCountInput.value, Number(colCountInput.min), Number(colCountInput.max));
    rowCountInput.value = String(rows);
    colCountInput.value = String(cols);
    resizeFloorGridForDimensions(rows, cols);
    lastMechanismUpdateKey = "";
    removeChildren(latticeRoot);
    removeChildren(pinRoot);
    removeChildren(lockMarkerRoot);
    cellViews = [];
    connections = [];
    connectionPinMeshes = [];
    selectableMeshes.length = 0;
    selectedDriveCell.row = Math.min(selectedDriveCell.row, rows - 1);
    selectedDriveCell.col = Math.min(selectedDriveCell.col, cols - 1);
    [...selectedCellKeys].forEach((key) => {
      const [row, col] = key.split(",").map(Number);
      if (row >= rows || col >= cols) selectedCellKeys.delete(key);
    });
    selectedCellKeys.add(cellKey(selectedDriveCell.row, selectedDriveCell.col));
    [...xyLocks.keys()].forEach((key) => {
      const [row, col] = key.split(",").map(Number);
      if (row >= rows || col >= cols) xyLocks.delete(key);
    });
    [...cellDriveCommands.keys()].forEach((key) => {
      const [row, col] = key.split(",").map(Number);
      if (row >= rows || col >= cols) cellDriveCommands.delete(key);
    });

    for (let row = 0; row < rows; row += 1) {
      const rowViews = [];
      for (let col = 0; col < cols; col += 1) {
        const cell = createCell(row, col);
        cell.lockMarker = createLockMarker();
        rowViews.push(cell);
        latticeRoot.add(cell.group);
      }
      cellViews.push(rowViews);
    }

    for (let row = 0; row < rows; row += 1) {
      for (let col = 0; col < cols; col += 1) {
        if (col + 1 < cols) {
          connections.push({ axis: "x", from: cellViews[row][col], to: cellViews[row][col + 1] });
        }
        if (row + 1 < rows) {
          connections.push({ axis: "y", from: cellViews[row][col], to: cellViews[row + 1][col] });
        }
      }
    }

    latticeSizeMetric.textContent = `${rows} x ${cols}`;
    cellCountMetric.textContent = String(rows * cols);
  }

  function cellAt(row, col) {
    return cellViews[row] && cellViews[row][col] ? cellViews[row][col] : null;
  }

  function siteWorld(cell, layer, site) {
    const rotation = layer === "top" ? cell.topRotation : cell.bottomRotation;
    return cell.center.clone().add(rotate2(SITE_VECTORS[site], rotation));
  }

  function siteWorldZ(cell, layer, site) {
    const rotation = layer === "top" ? cell.topRotation : cell.bottomRotation;
    const offset = rotate2(SITE_VECTORS[site], rotation);
    return (cell.z || 0) + (cell.tiltX || 0) * offset.y - (cell.tiltY || 0) * offset.x;
  }

  function cellNormal(cell) {
    return new THREE.Vector3(0, 0, 1).applyEuler(new THREE.Euler(cell.tiltX || 0, cell.tiltY || 0, 0, "XYZ")).normalize();
  }

  function cellPlaneZAt(cell, point) {
    return (cell.z || 0) + (cell.tiltX || 0) * (point.y - cell.center.y) - (cell.tiltY || 0) * (point.x - cell.center.x);
  }

  function holeRefPoint(ref) {
    return siteWorld(ref.cell, ref.layer, ref.site);
  }

  function makeHoleRef(cell, layer, site) {
    return { cell, layer, site };
  }

  function connectionPinPairs(connection) {
    const from = connection.from;
    const to = connection.to;
    if (connection.axis === "x") {
      if (from.parity === 0) {
        return [
          [makeHoleRef(from, "top", "east"), makeHoleRef(to, "bottom", "south")],
          [makeHoleRef(from, "bottom", "north"), makeHoleRef(to, "top", "south")],
        ];
      }
      return [
        [makeHoleRef(from, "bottom", "north"), makeHoleRef(to, "top", "west")],
        [makeHoleRef(from, "top", "north"), makeHoleRef(to, "bottom", "south")],
      ];
    }

    if (from.parity === 0) {
      return [
        [makeHoleRef(from, "top", "north"), makeHoleRef(to, "bottom", "east")],
        [makeHoleRef(from, "bottom", "west"), makeHoleRef(to, "top", "east")],
      ];
    }
    return [
      [makeHoleRef(from, "bottom", "west"), makeHoleRef(to, "top", "south")],
      [makeHoleRef(from, "top", "west"), makeHoleRef(to, "bottom", "east")],
    ];
  }

  function connectionPinPoints(connection) {
    return connectionPinPairs(connection).map(([first, second]) => holeRefPoint(first).add(holeRefPoint(second)).multiplyScalar(0.5));
  }

  function explicitConnectionPinPoints() {
    const points = [];
    connections.forEach((connection) => {
      connectionPinPoints(connection).forEach((point) => points.push(point));
    });
    return points;
  }

  function explicitConnectionPinPairs() {
    const pairs = [];
    connections.forEach((connection) => {
      connectionPinPairs(connection).forEach((pair) => pairs.push(pair));
    });
    return pairs;
  }

  function allArmHoleSites() {
    const sites = [];
    const siteNames = ["east", "north", "west", "south"];
    cellViews.forEach((row) => {
      row.forEach((cell) => {
        siteNames.forEach((site) => {
          sites.push({
            cell,
            layer: "top",
            site,
            point: siteWorld(cell, "top", site),
          });
          sites.push({
            cell,
            layer: "bottom",
            site,
            point: siteWorld(cell, "bottom", site),
          });
        });
      });
    });
    return sites;
  }

  function sameCell(first, second) {
    return first.cell.row === second.cell.row && first.cell.col === second.cell.col;
  }

  function overlappingHolePinPoints() {
    const tolerance = CAD.nominalHoleDiameterMm * 0.38;
    const binSize = Math.max(tolerance, 0.5);
    const bins = new Map();
    const groups = [];

    allArmHoleSites().forEach((site) => {
      const bx = Math.round(site.point.x / binSize);
      const by = Math.round(site.point.y / binSize);
      let matched = null;

      for (let dx = -1; dx <= 1 && !matched; dx += 1) {
        for (let dy = -1; dy <= 1 && !matched; dy += 1) {
          const nearby = bins.get(`${bx + dx}|${by + dy}`) || [];
          matched = nearby.find((group) => group.point.distanceTo(site.point) <= tolerance);
        }
      }

      if (!matched) {
        matched = { point: site.point.clone(), sites: [] };
        groups.push(matched);
        const key = `${bx}|${by}`;
        if (!bins.has(key)) bins.set(key, []);
        bins.get(key).push(matched);
      }

      matched.sites.push(site);
      matched.point.multiplyScalar(matched.sites.length - 1).add(site.point).multiplyScalar(1 / matched.sites.length);
    });

    return groups
      .filter((group) => group.sites.some((first, firstIndex) => group.sites.some((second, secondIndex) => firstIndex !== secondIndex && !sameCell(first, second))))
      .map((group) => group.point);
  }

  function dedupePinPoints(points) {
    const dedupeTolerance = CAD.nominalHoleDiameterMm * 0.28;
    const deduped = [];
    points.forEach((point) => {
      if (!deduped.some((existing) => existing.distanceTo(point) <= dedupeTolerance)) {
        deduped.push(point);
      }
    });
    return deduped;
  }

  function cellDriveMagnitude(ref, contactField) {
    return Math.abs(contactField?.deltaField?.[ref.cell.row]?.[ref.cell.col] || 0);
  }

  function snapshotHoleRefPoint(ref, snapshot) {
    const entry = snapshot?.[ref.cell.row]?.[ref.cell.col];
    if (!entry) return holeRefPoint(ref);
    const rotation = ref.layer === "top" ? entry.topRotation : entry.bottomRotation;
    return entry.center.clone().add(rotate2(SITE_VECTORS[ref.site], rotation));
  }

  function contactPinPointForPair(first, second, contactField) {
    const firstPoint = holeRefPoint(first);
    const secondPoint = holeRefPoint(second);
    const firstDrive = cellDriveMagnitude(first, contactField);
    const secondDrive = cellDriveMagnitude(second, contactField);
    const driverPoint = firstDrive >= secondDrive ? firstPoint : secondPoint;
    const restFirstPoint = snapshotHoleRefPoint(first, latestRestSnapshot);
    const restSecondPoint = snapshotHoleRefPoint(second, latestRestSnapshot);
    const restPinPoint = restFirstPoint.add(restSecondPoint).multiplyScalar(0.5);
    if (xyLocks.size <= 1) restPinPoint.add(new THREE.Vector2(freeTranslation[0], freeTranslation[1]));
    const displacement = driverPoint.clone().sub(restPinPoint);
    const separation = displacement.length();
    const clearance = pinHoleClearanceMm();
    if (separation <= clearance || separation <= 1e-9) {
      return restPinPoint;
    }
    return restPinPoint.add(displacement.multiplyScalar((separation - clearance) / separation));
  }

  function contactPinPoints(contactField) {
    return explicitConnectionPinPairs().map(([first, second]) => contactPinPointForPair(first, second, contactField));
  }

  function clampPointToDisk(point, center, radius) {
    if (radius <= 1e-9) return center.clone();
    const offset = point.clone().sub(center);
    const distance = offset.length();
    if (distance <= radius || distance <= 1e-9) return point;
    return center.clone().add(offset.multiplyScalar(radius / distance));
  }

  function feasiblePinCenter(firstPoint, secondPoint, targetPoint) {
    const clearance = Math.max(0, holeRadiusMm() - pinRadiusMm());
    const separation = firstPoint.distanceTo(secondPoint);
    let point;
    if (separation >= 2 * clearance || clearance < 1e-12) {
      point = firstPoint.clone().add(secondPoint).multiplyScalar(0.5);
    } else {
      const candidates = [targetPoint.clone(), clampPointToDisk(targetPoint.clone(), firstPoint, clearance),
        clampPointToDisk(targetPoint.clone(), secondPoint, clearance)]
        .filter((p) => p.distanceTo(firstPoint) <= clearance + 1e-10 && p.distanceTo(secondPoint) <= clearance + 1e-10);
      if (separation > 1e-12) {
        const midpoint = firstPoint.clone().add(secondPoint).multiplyScalar(0.5);
        const delta = secondPoint.clone().sub(firstPoint).multiplyScalar(1 / separation);
        const normal = new THREE.Vector2(-delta.y, delta.x);
        const height = Math.sqrt(Math.max(0, clearance * clearance - separation * separation / 4));
        candidates.push(midpoint.clone().addScaledVector(normal, height), midpoint.clone().addScaledVector(normal, -height));
      }
      point = candidates.sort((a, b) => a.distanceToSquared(targetPoint) - b.distanceToSquared(targetPoint))[0];
    }
    const firstResidual = Math.max(0, point.distanceTo(firstPoint) - clearance);
    const secondResidual = Math.max(0, point.distanceTo(secondPoint) - clearance);
    return {
      point,
      radialResidual: Math.max(firstResidual, secondResidual),
      centerSeparation: firstPoint.distanceTo(secondPoint),
      clearance,
    };
  }

  function connectionPinItems(contactField = null) {
    return explicitConnectionPinPairs().map(([first, second], pinIndex) => {
      const firstPoint = holeRefPoint(first);
      const secondPoint = holeRefPoint(second);
      const previousPin = contactField && contactMemory?.pins?.[pinIndex];
      const savedPin = previousPin && xyLocks.size <= 1
        ? previousPin.map((value, axis) => value + freeTranslation[axis] - (contactMemory.freeTranslation?.[axis] || 0))
        : previousPin;
      const targetPoint = savedPin ? new THREE.Vector2(savedPin[0], savedPin[1])
        : contactField ? contactPinPointForPair(first, second, contactField) : firstPoint.clone().add(secondPoint).multiplyScalar(0.5);
      const fit = feasiblePinCenter(firstPoint, secondPoint, targetPoint);
      const bores = [first, second].map((ref) => {
        const offset = rotate2(SITE_VECTORS[ref.site], ref.layer === "top" ? ref.cell.topRotation : ref.cell.bottomRotation);
        const center = new THREE.Vector3(offset.x, offset.y, (ref.layer === "top" ? 1 : -1) * plateThicknessMm() * 0.5)
          .applyQuaternion(ref.cell.group.quaternion).add(ref.cell.group.position);
        return { center, normal: cellNormal(ref.cell) };
      });
      const axis = bores[0].normal.clone().add(bores[1].normal);
      if (axis.lengthSq() <= 1e-9) axis.set(0, 0, 1);
      axis.normalize();
      const midpoint = bores[0].center.clone().add(bores[1].center).multiplyScalar(0.5);
      const planar = bores.every((bore) => bore.normal.z > 1 - 1e-12);
      const origin = savedPin && !planar
        ? new THREE.Vector3(...savedPin) : new THREE.Vector3(fit.point.x, fit.point.y, midpoint.z);
      const thickness = plateThicknessMm();
      // Maximize the actual elliptical section at both bore faces. The norm
      // is convex through the thickness, so the two faces bound the full bore.
      function boreFit(bore) {
        const cosine = clamp(axis.dot(bore.normal), 1e-9, 1);
        const tangent = Math.sqrt(Math.max(0, 1 - cosine * cosine)) / cosine;
        const along = bore.center.clone().sub(origin).dot(bore.normal) / cosine;
        const offset = origin.clone().addScaledVector(axis, along).sub(bore.center);
        const u = axis.clone().addScaledVector(bore.normal, -cosine);
        if (u.lengthSq() < 1e-20) {
          u.crossVectors(bore.normal, Math.abs(bore.normal.x) < 0.9 ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 1, 0));
        }
        u.normalize();
        const v = new THREE.Vector3().crossVectors(bore.normal, u).normalize();
        let excess = -Infinity, correction = new THREE.Vector3();
        for (const face of [-1, 1]) {
          const shift = axis.clone().multiplyScalar(face * thickness * 0.5 / cosine)
            .addScaledVector(bore.normal, -face * thickness * 0.5).add(offset);
          const extent = window.RADPinGeometry.ellipseExtent(shift.dot(u), shift.dot(v), pinRadiusMm() / cosine, pinRadiusMm());
          if (extent.radius - holeRadiusMm() > excess) {
            excess = extent.radius - holeRadiusMm();
            correction = u.clone().multiplyScalar(extent.x).addScaledVector(v, extent.y).normalize().multiplyScalar(Math.max(0, excess));
          }
        }
        return { offset, excess, correction, along, cosine, tangent };
      }
      for (let iteration = 0; iteration < 48; iteration += 1) {
        const fits = bores.map(boreFit);
        if (fits.every((fit) => fit.excess <= 1e-8)) break;
        const correction = fits[0].correction.clone().add(fits[1].correction).multiplyScalar(0.5);
        if (correction.lengthSq() < 1e-24) break;
        origin.sub(correction);
      }
      // Keep the fixed-length pin centered axially without changing its line.
      origin.addScaledVector(axis, midpoint.clone().sub(origin).dot(axis));
      const fits = bores.map(boreFit);
      const radialResidual = Math.max(0, ...fits.map((item) => item.offset.length() - pinHoleClearanceMm()));
      const tiltResidual = Math.max(0, ...fits.map((item) => item.excess));
      const axialResidual = Math.max(0, ...fits.map((item) => Math.abs(item.along) + thickness * 0.5 / item.cosine + pinRadiusMm() * item.tangent - CAD.bodyThicknessMm * 2.95 * 0.5));
      return {
        point: new THREE.Vector2(origin.x, origin.y),
        z: origin.z,
        axis,
        radialResidual,
        tiltResidual,
        axialResidual,
        centerSeparation: fit.centerSeparation,
        clearance: fit.clearance,
      };
    });
  }

  function constrainedLayoutSolve(baseCenters, pitch) {
    const rows = baseCenters.length;
    const cols = rows ? baseCenters[0].length : 0;
    const centers = baseCenters.map((row) => row.map((point) => point.clone()));
    const z = Array.from({ length: rows }, () => Array(cols).fill(0));
    const tiltX = z.map((row) => row.slice());
    const tiltY = z.map((row) => row.slice());
    const locks = [...xyLocks.values()].filter((lock) => lock.row < rows && lock.col < cols);
    const report = { activeLocks: locks.length, maxLift: 0, maxTiltDeg: 0, maxLockDrift: 0,
      maxEdgeCompression: 0, maxAdmissibleEdgeRise: 0, backlashLimitedLift: false,
      maxEdgeLengthError: 0, iterations: 0, constrained: locks.length > 0 };
    if (locks.length < 2) return { centers, z, tiltX, tiltY, report };

    const index = (row, col) => row * cols + col;
    const reference = baseCenters.flatMap((row, r) => row.map((point, c) =>
      (locks.length > 1 && latestRestSnapshot?.[r]?.[c]?.center) || point));
    const fixed = new Map(locks.map((lock) => [index(lock.row, lock.col), lock]));
    const positions = reference.map((point) => {
      let dx = 0, dy = 0, dz = 0, weightSum = 0;
      locks.forEach((lock) => {
        const anchor = reference[index(lock.row, lock.col)];
        const weight = 1 / Math.max(1, point.distanceToSquared(anchor));
        dx += weight * (lock.targetX - anchor.x);
        dy += weight * (lock.targetY - anchor.y);
        dz += weight * lock.targetZ;
        weightSum += weight;
      });
      return new THREE.Vector3(point.x + dx / weightSum, point.y + dy / weightSum, dz / weightSum);
    });
    fixed.forEach((lock, i) => positions[i].set(lock.targetX, lock.targetY, lock.targetZ));
    const edges = connections.map(({ from, to }) => ({
      a: index(from.row, from.col), b: index(to.row, to.col),
      length: 0.5 * (pitch[from.row][from.col] + pitch[to.row][to.col]),
    }));
    const influence = positions.map((_, i) =>
      Math.abs(latestContactField?.deltaField?.[Math.floor(i / cols)]?.[i % cols] || 0));
    const maxInfluence = Math.max(1e-12, ...influence);
    const support = new Set(influence.map((value, i) => value > 1e-9 ? i : -1).filter((i) => i >= 0));
    // A transition ring allows a local arch to meet the unactuated straight tails.
    const active = new Set(support);
    const excitation = influence.slice();
    edges.forEach(({ a, b }) => {
      if (active.has(a)) support.add(b);
      if (active.has(b)) support.add(a);
      excitation[a] = Math.max(excitation[a], influence[b]);
      excitation[b] = Math.max(excitation[b], influence[a]);
    });
    const isString = rows === 1 || cols === 1;
    const hasFreeZLock = locks.some((lock) => lock.lockZ === false);
    const rigidStringTilt = isString && locks.length >= 2 && exactPinTiltLimitRad() <= 1e-9
      && locks.filter((lock) => lock.lockZ !== false).length < 2;
    const canBend = locks.length >= 2 && (exactPinTiltLimitRad() > 1e-9 || rigidStringTilt);
    const flatAnchors = locks.every((lock) => Math.abs(lock.targetZ) < 1e-9);
    const seed = positions.map(() => 0);
    edges.forEach(({ a, b, length }) => {
      const distance = positions[a].distanceTo(positions[b]);
      report.maxEdgeCompression = Math.max(report.maxEdgeCompression, length - distance);
      const rise = Math.sqrt(Math.max(0, length * length - distance * distance));
      seed[a] = Math.max(seed[a], rise * influence[a] / maxInfluence);
      seed[b] = Math.max(seed[b], rise * influence[b] / maxInfluence);
    });
    const mobility = positions.map((point, i) => {
      if (canBend && support.has(i) && (!fixed.has(i) || fixed.get(i).lockZ === false)) point.z += seed[i];
      if (fixed.has(i)) return [0, 0, canBend && fixed.get(i).lockZ === false ? 4 : 0];
      return [cols === 1 ? 0 : 0.25, rows === 1 ? 0 : 0.25,
        canBend ? 4 : 0];
    });
    const globalString = canBend && isString;
    if (rigidStringTilt) {
      const anchors = [...fixed.keys()].sort((a, b) => a - b);
      const pitches = pitch.flat();
      const distances = [0];
      for (let i = 1; i < positions.length; i += 1)
        distances.push(distances[i - 1] + 0.5 * (pitches[i - 1] + pitches[i]));
      const first = anchors[0], last = anchors.at(-1);
      const start = positions[first].clone(), end = positions[last];
      const span = distances[last] - distances[first];
      const horizontal = new THREE.Vector3(end.x - start.x, end.y - start.y, 0);
      const slope = Math.sqrt(Math.max(0, 1 - horizontal.lengthSq() / (span * span)));
      const heightAnchor = anchors.find((i) => fixed.get(i).lockZ !== false) ?? first;
      const height = fixed.get(heightAnchor).targetZ;
      positions.forEach((point, i) => {
        point.copy(start).addScaledVector(horizontal, (distances[i] - distances[first]) / span);
        point.z = height + slope * (distances[i] - distances[heightAnchor]);
      });
    } else if (globalString) {
      let feasibleBranch = true;
      const allAnchors = [...fixed.keys()].sort((a, b) => a - b);
      const anchors = allAnchors.filter((i) => fixed.get(i).lockZ !== false
        || i === allAnchors[0] || i === allAnchors.at(-1));
      const pitches = pitch.flat();
      for (let span = 1; span < anchors.length; span += 1) {
        const first = anchors[span - 1], last = anchors[span];
        const start = positions[first], end = positions[last];
        const axis = end.clone().sub(start);
        const chord = axis.length();
        if (chord < 1e-9) { feasibleBranch = false; continue; }
        axis.divideScalar(chord);
        const up = new THREE.Vector3(0, 0, 1).addScaledVector(axis, -axis.z).normalize();
        const lengths = Array.from({ length: last - first }, (_, i) =>
          0.5 * (pitches[first + i] + pitches[first + i + 1]));
        const total = lengths.reduce((sum, length) => sum + length, 0);
        // A circular polygon preserves each edge length, including unequal
        // pitches. Curvature is global; dilation die-off does not fix cell Z.
        const turn = (curvature) => lengths.reduce((sum, length) =>
          sum + 2 * Math.asin(clamp(0.5 * curvature * length, 0, 1)), 0);
        if (total <= chord + 1e-9) {
          let distance = 0;
          for (let i = first + 1; i < last; i += 1) {
            distance += lengths[i - first - 1];
            positions[i].copy(start).addScaledVector(axis, distance);
          }
          continue;
        }
        let low = 0, high = 2 / Math.max(...lengths);
        for (let i = 0; i < 60; i += 1) {
          const mid = 0.5 * (low + high);
          if (turn(mid) < Math.PI) low = mid; else high = mid;
        }
        const spanAt = (curvature) => 2 * Math.sin(0.5 * turn(curvature)) / curvature;
        // Reject geometry outside this continuous minor-arc branch rather than
        // switching branches or inventing an extension of the mechanism.
        if (spanAt(high) > chord || up.lengthSq() < 1e-12) {
          if (hasFreeZLock && last === first + 1 && (fixed.get(first).lockZ === false || fixed.get(last).lockZ === false)) {
            const rise = Math.sqrt(Math.max(0, total * total - chord * chord));
            if (fixed.get(last).lockZ === false) end.z += rise; else start.z += rise;
          } else feasibleBranch = false;
          continue;
        }
        low = 0;
        for (let i = 0; i < 60; i += 1) {
          const mid = 0.5 * (low + high);
          if (spanAt(mid) > chord) low = mid; else high = mid;
        }
        const curvature = 0.5 * (low + high);
        const halfTurn = 0.5 * turn(curvature);
        let angle = -halfTurn;
        for (let i = first + 1; i < last; i += 1) {
          angle += 2 * Math.asin(clamp(0.5 * curvature * lengths[i - first - 1], 0, 1));
          positions[i].copy(start)
            .addScaledVector(axis, (Math.sin(angle) + Math.sin(halfTurn)) / curvature)
            .addScaledVector(up, (Math.cos(angle) - Math.cos(halfTurn)) / curvature);
        }
      }
      // Free tails follow the tangent of the anchored span; they are not
      // additional positional or orientation constraints at an internal lock.
      const first = allAnchors[0], last = allAnchors.at(-1);
      const before = positions[first + 1].clone().sub(positions[first]).normalize();
      const after = positions[last].clone().sub(positions[last - 1]).normalize();
      for (let i = first - 1; i >= 0; i -= 1)
        positions[i].copy(positions[i + 1]).addScaledVector(before, -0.5 * (pitches[i] + pitches[i + 1]));
      for (let i = last + 1; i < positions.length; i += 1)
        positions[i].copy(positions[i - 1]).addScaledVector(after, 0.5 * (pitches[i - 1] + pitches[i]));
      fixed.forEach((lock, i) => {
        positions[i].x = lock.targetX;
        positions[i].y = lock.targetY;
        if (lock.lockZ !== false) positions[i].z = lock.targetZ;
      });
      report.maxEdgeLengthError = edges.reduce((error, edge) =>
        Math.max(error, Math.abs(positions[edge.a].distanceTo(positions[edge.b]) - edge.length)), 0);
      if (!feasibleBranch) report.maxEdgeLengthError = Infinity;
    }
    // Project the shared 3D edge constraints together, preserving XYZ anchors
    // exactly. Sheets retain their existing branch; neither is a force equilibrium.
    const corrections = positions.map(() => new THREE.Vector3());
    const degree = positions.map(() => 0);
    edges.forEach(({ a, b }) => { degree[a] += 1; degree[b] += 1; });
    const relaxation = 1 / Math.max(1, ...degree);
    for (let iteration = 0; !rigidStringTilt && (!globalString || hasFreeZLock) && iteration < 3200; iteration += 1) {
      corrections.forEach((correction) => correction.set(0, 0, 0));
      for (let offset = 0; offset < edges.length; offset += 1) {
        const edge = edges[iteration % 2 ? edges.length - 1 - offset : offset];
        const a = positions[edge.a], b = positions[edge.b];
        const distance = a.distanceTo(b);
        if (distance < 1e-12) continue;
        const nx = (b.x - a.x) / distance, ny = (b.y - a.y) / distance, nz = (b.z - a.z) / distance;
        const wa = mobility[edge.a], wb = mobility[edge.b];
        const denominator = nx * nx * (wa[0] + wb[0]) + ny * ny * (wa[1] + wb[1]) + nz * nz * (wa[2] + wb[2]);
        if (denominator < 1e-16) continue;
        const correction = (distance - edge.length) / denominator;
        corrections[edge.a].x += correction * wa[0] * nx;
        corrections[edge.a].y += correction * wa[1] * ny;
        corrections[edge.a].z += correction * wa[2] * nz;
        corrections[edge.b].x -= correction * wb[0] * nx;
        corrections[edge.b].y -= correction * wb[1] * ny;
        corrections[edge.b].z -= correction * wb[2] * nz;
      }
      // Jacobi updates are independent of traversal order and grid orientation.
      positions.forEach((point, i) => point.addScaledVector(corrections[i], relaxation));
      report.maxEdgeLengthError = edges.reduce((error, edge) =>
        Math.max(error, Math.abs(positions[edge.a].distanceTo(positions[edge.b]) - edge.length)), 0);
      report.iterations = iteration + 1;
      if (report.maxEdgeLengthError < 1e-6) break;
    }

    for (let row = 0; row < rows; row += 1) {
      for (let col = 0; col < cols; col += 1) {
        const i = index(row, col), point = positions[i];
        centers[row][col].set(point.x, point.y);
        z[row][col] = point.z;
        const tangentX = cols > 1
          ? positions[index(row, Math.min(cols - 1, col + 1))].clone().sub(positions[index(row, Math.max(0, col - 1))])
          : new THREE.Vector3(1, 0, 0);
        const tangentY = rows > 1
          ? positions[index(Math.min(rows - 1, row + 1), col)].clone().sub(positions[index(Math.max(0, row - 1), col)])
          : new THREE.Vector3(0, 1, 0);
        const normal = new THREE.Vector3().crossVectors(tangentX, tangentY).normalize();
        if (normal.lengthSq() < 1e-12) normal.set(0, 0, 1);
        if (normal.z < 0) normal.negate();
        tiltX[row][col] = Math.atan2(-normal.y, normal.z);
        tiltY[row][col] = Math.asin(clamp(normal.x, -1, 1));
        report.maxLift = Math.max(report.maxLift, point.z);
        report.maxTiltDeg = Math.max(report.maxTiltDeg, radToDeg(Math.acos(clamp(normal.z, -1, 1))));
      }
    }
    fixed.forEach((lock, i) => {
      report.maxLockDrift = Math.max(report.maxLockDrift,
        Math.hypot(positions[i].x - lock.targetX, positions[i].y - lock.targetY,
          lock.lockZ === false ? 0 : positions[i].z - lock.targetZ));
    });
    report.maxEdgeLengthError = Math.max(report.maxEdgeLengthError, report.maxLockDrift,
      ...edges.map((edge) => Math.abs(positions[edge.a].distanceTo(positions[edge.b]) - edge.length)));
    return { centers, z, tiltX, tiltY, report };
  }

  function buildCenterLineGeometry() {
    const points = [];
    for (let row = 0; row < cellViews.length; row += 1) {
      for (let col = 0; col < cellViews[row].length; col += 1) {
        const cell = cellViews[row][col];
        if (col + 1 < cellViews[row].length) {
          const next = cellViews[row][col + 1];
          points.push(new THREE.Vector3(cell.center.x, cell.center.y, cell.z + 0.3), new THREE.Vector3(next.center.x, next.center.y, next.z + 0.3));
        }
        if (row + 1 < cellViews.length) {
          const next = cellViews[row + 1][col];
          points.push(new THREE.Vector3(cell.center.x, cell.center.y, cell.z + 0.3), new THREE.Vector3(next.center.x, next.center.y, next.z + 0.3));
        }
      }
    }
    latticeLine.geometry.dispose();
    latticeLine.geometry = new THREE.BufferGeometry().setFromPoints(points);
  }

  function updateConnectionPins() {
    const visible = showPinsInput.checked;
    const contactMode = actuationModeInput.value === "contact" && latestContactField;
    const pinItems = connectionPinItems(contactMode ? latestContactField : null);
    ensureConnectionPinMeshCount(pinItems.length);
    updatePinGeometry();
    let visibleCount = 0;
    let maxRadialResidual = 0;
    let maxTiltResidual = 0;
    let maxAxialResidual = 0;
    const pinAxisLocal = new THREE.Vector3(0, 1, 0);
    maxAxialResidual = Math.max(0, plateThicknessMm() - CAD.bodyThicknessMm * 2.55 * 0.5);
    connectionPinMeshes.forEach((pin, index) => {
      const item = pinItems[index];
      pin.visible = Boolean(visible && item);
      if (item) {
        visibleCount += 1;
        maxRadialResidual = Math.max(maxRadialResidual, item.radialResidual || 0);
        maxTiltResidual = Math.max(maxTiltResidual, item.tiltResidual || 0);
        maxAxialResidual = Math.max(maxAxialResidual, item.axialResidual || 0);
        pin.position.set(item.point.x, item.point.y, item.z);
        pin.quaternion.setFromUnitVectors(pinAxisLocal, item.axis || new THREE.Vector3(0, 0, 1));
      }
    });
    latestPinFitReport = {
      expectedPins: explicitConnectionPinPairs().length,
      visiblePins: visible ? visibleCount : 0,
      maxRadialResidual,
      maxTiltResidual,
      maxAxialResidual,
      clear: maxRadialResidual <= BODY_EPSILON_MM && maxTiltResidual <= BODY_EPSILON_MM && maxAxialResidual <= BODY_EPSILON_MM && (latestBoundaryReport?.maxEdgeLengthError || 0) <= BODY_EPSILON_MM,
    };
  }

  function updateContactReadout(contactField = latestContactField) {
    const ratio = pinRadiusRatio();
    const holeRadius = holeRadiusMm();
    const pinRadius = pinRadiusMm();
    const diametral = diametralBacklashMm();
    const clearance = pinHoleClearanceMm();
    const normalized = normalizedBacklash();
    const backlash = researchAngularBacklashDeg();
    const command =
      contactField && actuationModeInput.value === "contact"
        ? contactField.deltaField?.[contactField.sourceRow]?.[contactField.sourceCol] || 0
        : 0;
    const transmitted = Math.sign(command) * Math.max(0, Math.abs(command) - backlash);
    const floating = Math.abs(command) <= backlash + 1e-9;

    pinRadiusRatioInput.value = ratio.toFixed(2);
    pinRadiusRatioOut.textContent = ratio.toFixed(2);
    plateThicknessInput.value = plateThicknessMm().toFixed(2);
    holeRadiusMetric.textContent = `${holeRadius.toFixed(3)} mm`;
    pinRadiusMetric.textContent = `${pinRadius.toFixed(3)} mm`;
    diametralBacklashMetric.textContent = `${diametral.toFixed(3)} mm`;
    radialClearanceMetric.textContent = `${clearance.toFixed(3)} mm`;
    freeAngleMetric.textContent = `${backlash.toFixed(2)} deg`;
    relativeBacklashMetric.textContent = normalized.toFixed(4);
    exactTiltMetric.textContent = `${exactPinTiltLimitDeg().toFixed(2)} deg`;
    computedContactCommandMetric.textContent = `${command.toFixed(1)} deg`;
    contactStateMetric.textContent = contactField && actuationModeInput.value === "contact" ? (floating ? "floating in clearance" : "contact engaged") : "uniform mode";
    transmittedCommandMetric.textContent = `${transmitted.toFixed(1)} deg`;
  }

  function updateSelectionMarker() {
    const cell = cellAt(selectedDriveCell.row, selectedDriveCell.col);
    removeChildren(selectionMarkerRoot);
    if (!cell) {
      selectionRing.visible = false;
      return;
    }
    const points = [];
    const radius = CAD.siteRadiusMm * 1.38;
    for (let index = 0; index <= 64; index += 1) {
      const angle = (index / 64) * Math.PI * 2;
      points.push(new THREE.Vector3(Math.cos(angle) * radius, Math.sin(angle) * radius, 0));
    }
    selectionRing.geometry.dispose();
    selectionRing.geometry = new THREE.BufferGeometry().setFromPoints(points);
    selectionRing.position.set(cell.center.x, cell.center.y, (cell.z || 0) + CAD.bodyThicknessMm * 1.9);
    selectionRing.visible = true;
    const activeKey = cellKey(selectedDriveCell.row, selectedDriveCell.col);
    selectedCellKeys.forEach((key) => {
      if (key === activeKey) return;
      const [row, col] = key.split(",").map(Number);
      const selectedCell = cellAt(row, col);
      if (!selectedCell) return;
      const marker = new THREE.LineLoop(new THREE.BufferGeometry().setFromPoints(points), materials.selectionSet);
      marker.position.set(selectedCell.center.x, selectedCell.center.y, (selectedCell.z || 0) + CAD.bodyThicknessMm * 1.8);
      selectionMarkerRoot.add(marker);
    });
  }

  function updateLockMarkers() {
    cellViews.forEach((row) => {
      row.forEach((cell) => {
        const lock = lockForCell(cell.row, cell.col);
        if (!cell.lockMarker) return;
        cell.lockMarker.visible = Boolean(lock);
        if (lock) cell.lockMarker.position.set(lock.targetX, lock.targetY,
          (lock.lockZ === false ? cell.z : lock.targetZ || 0) + CAD.bodyThicknessMm * 2.25);
      });
    });
  }

  function updateBoundaryReadout(report = latestBoundaryReport) {
    const selectedLock = lockForCell(selectedDriveCell.row, selectedDriveCell.col);
    const lockedSelectedCount = selectedCells().filter((cell) => lockForCell(cell.row, cell.col)).length;
    selectedCellsMetric.textContent = selectedCellsLabel();
    selectedLockMetric.textContent = selectedLock
      ? selectedLock.lockZ === false
        ? `XY @ (${selectedLock.targetX.toFixed(1)}, ${selectedLock.targetY.toFixed(1)}); Z free`
        : `XYZ @ (${selectedLock.targetX.toFixed(1)}, ${selectedLock.targetY.toFixed(1)}, ${selectedLock.targetZ.toFixed(2)})`
      : selectedCellKeys.size > 1
      ? `${lockedSelectedCount}/${selectedCellKeys.size} selected locked`
      : "free";
    activeLocksMetric.textContent = String(xyLocks.size);
    maxLiftMetric.textContent = `${(report?.maxLift || 0).toFixed(3)} mm`;
    maxTiltMetric.textContent = `${(report?.maxTiltDeg || 0).toFixed(2)} deg`;
    maxLockDriftMetric.textContent = `${(report?.maxLockDrift || 0).toFixed(3)} mm`;
  }

  function addAxisLine(points, color) {
    const geometry = new THREE.BufferGeometry().setFromPoints(points.map((p) => new THREE.Vector3(...p)));
    guideRoot.add(new THREE.Line(geometry, new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.75 })));
  }

  addAxisLine(
    [
      [-72, 0, -8],
      [72, 0, -8],
    ],
    0xb24d3a
  );
  addAxisLine(
    [
      [0, -72, -8],
      [0, 72, -8],
    ],
    0x1d6fd6
  );
  addAxisLine(
    [
      [0, 0, -12],
      [0, 0, 22],
    ],
    0x22a36f
  );

  const floorGrid = new THREE.GridHelper(260, 26, 0x8b97a5, 0xc4ccd4);
  floorGrid.rotation.x = Math.PI / 2;
  floorGrid.position.z = -8.1;
  guideRoot.add(floorGrid);

  scene.add(new THREE.HemisphereLight(0xffffff, 0x5f6974, 1.65));
  const key = new THREE.DirectionalLight(0xffffff, 2.15);
  key.position.set(36, -48, 82);
  key.castShadow = true;
  key.shadow.mapSize.width = 1024;
  key.shadow.mapSize.height = 1024;
  scene.add(key);
  const fill = new THREE.DirectionalLight(0xb8d4ff, 0.65);
  fill.position.set(-46, 54, 42);
  scene.add(fill);

  function updateCamera() {
    cameraState.radius = clamp(cameraState.radius, minCameraRadius(), maxCameraRadius());
    const r = cameraState.radius;
    const ce = Math.cos(cameraState.elevation);
    camera.position.set(
      cameraState.target.x + r * ce * Math.cos(cameraState.azimuth),
      cameraState.target.y + r * ce * Math.sin(cameraState.azimuth),
      cameraState.target.z + r * Math.sin(cameraState.elevation)
    );
    updateCameraClipPlanes();
    camera.lookAt(cameraState.target);
  }

  function framedCameraRadius() {
    const verticalHalfFov = degToRad(camera.fov) * 0.5;
    const horizontalHalfFov = Math.atan(Math.tan(verticalHalfFov) * Math.max(0.01, camera.aspect || 1));
    return Math.max(180, latestLatticeRadius * 1.08 / Math.sin(Math.min(verticalHalfFov, horizontalHalfFov)));
  }

  function minCameraRadius() {
    return Math.max(52, Math.min(180, latestLatticeRadius * 0.16));
  }

  function maxCameraRadius() {
    return Math.max(280, framedCameraRadius() * 1.75);
  }

  function updateCameraClipPlanes() {
    const far = Math.max(3000, maxCameraRadius() + latestLatticeRadius * 2.8 + CAD.cellWidthMm * 4);
    if (Math.abs(camera.near - 1) > 1e-9 || Math.abs(camera.far - far) > 1) {
      camera.near = 1;
      camera.far = far;
      camera.updateProjectionMatrix();
    }
  }

  function fitCurrentView() {
    cameraState.radius = framedCameraRadius();
    cameraState.target.set(0, 0, 0);
    updateCamera();
  }

  function setView(view) {
    currentView = view || "iso";
    const radius = framedCameraRadius();
    if (view === "top") {
      cameraState.azimuth = -Math.PI / 2;
      cameraState.elevation = Math.PI / 2 - 0.001;
      cameraState.radius = radius;
    } else if (view === "front") {
      cameraState.azimuth = -Math.PI / 2;
      cameraState.elevation = 0.03;
      cameraState.radius = radius;
    } else if (view === "side") {
      cameraState.azimuth = 0;
      cameraState.elevation = 0.05;
      cameraState.radius = radius;
    } else {
      cameraState.azimuth = -0.82;
      cameraState.elevation = 0.58;
      cameraState.radius = radius;
    }
    cameraState.target.set(0, 0, 0);
    updateCamera();
    document.querySelectorAll("[data-view]").forEach((button) => {
      button.setAttribute("aria-pressed", String(button.dataset.view === view));
    });
  }

  function angleLabel(requestedDeg, effectiveRad) {
    const effectiveDeg = radToDeg(effectiveRad);
    if (Math.abs(requestedDeg - effectiveDeg) > 0.25) {
      return `${requestedDeg.toFixed(0)} -> ${effectiveDeg.toFixed(0)} deg`;
    }
    const travel = Number(aTopInput.max) - Number(aTopInput.min);
    return `${effectiveDeg.toFixed(travel < 0.01 ? 6 : 3)} deg`;
  }

  function updateCollisionReadout(solution) {
    const report = solution.report || { maxPenetration: 0, minClearance: 0 };
    collisionState.classList.remove("status-ok", "status-adjusted", "status-blocked");
    if (solution.status === "blocked") {
      collisionState.classList.add("status-blocked");
      collisionState.textContent = "blocked";
    } else if (solution.status === "adjusted") {
      collisionState.classList.add("status-adjusted");
      collisionState.textContent = "adjusted";
    } else {
      collisionState.classList.add("status-ok");
      collisionState.textContent = "clear";
    }
    penetrationMetric.textContent = `${report.maxPenetration.toFixed(3)} mm`;
    clearanceMetric.textContent = `${Math.max(0, report.minClearance).toFixed(3)} mm`;
    effectiveBLowerMetric.textContent = `${radToDeg(solution.pose.bLowerRad).toFixed(1)} deg`;
    effectiveBTopMetric.textContent = `${radToDeg(solution.pose.bTopRad).toFixed(1)} deg`;
    effectiveCLowerMetric.textContent = `${radToDeg(solution.pose.cLowerRad).toFixed(1)} deg`;
    effectiveCTopMetric.textContent = `${radToDeg(solution.pose.cTopWorldRad - solution.pose.cLowerRad).toFixed(1)} deg`;
  }

  function centerPitchCoordinates(xRows, yCols) {
    // Center every strand separately: choosing its first node as a reference
    // otherwise introduces a spurious shear when the same sheet is rotated.
    const rows = xRows.length, cols = xRows[0].length;
    for (let row = 0; row < rows; row += 1) {
      const mean = xRows[row].reduce((sum, x) => sum + x, 0) / cols;
      for (let col = 0; col < cols; col += 1) xRows[row][col] -= mean;
    }
    for (let col = 0; col < cols; col += 1) {
      const mean = yCols.reduce((sum, row) => sum + row[col], 0) / rows;
      for (let row = 0; row < rows; row += 1) yCols[row][col] -= mean;
    }
  }

  function buildLatticeSnapshot(driveField = null, requestedBLower = degToRad(B_LOWER_REFERENCE_DEG), requestedBTop = degToRad(B_TOP_REFERENCE_DEG)) {
    const rows = cellViews.length;
    const cols = rows ? cellViews[0].length : 0;
    const snapshot = Array.from({ length: rows }, () => Array.from({ length: cols }, () => null));
    const pitch = Array.from({ length: rows }, () => Array.from({ length: cols }, () => CAD.cellWidthMm));
    const xRows = Array.from({ length: rows }, () => Array.from({ length: cols }, () => 0));
    const yCols = Array.from({ length: rows }, () => Array.from({ length: cols }, () => 0));

    for (let row = 0; row < rows; row += 1) {
      for (let col = 0; col < cols; col += 1) {
        const commandDeg = driveField ? driveField[row][col] : 0;
        const solution = solveNonPenetratingPose(degToRad(commandDeg), requestedBLower, requestedBTop);
        const localFrame = displayFrameForPose(solution.pose);
        const localStep = localFrame.mapPoint(solution.pose.bCenter).sub(localFrame.mapPoint(new THREE.Vector2(0, 0)));
        snapshot[row][col] = { pose: solution.pose, frame: localFrame, center: new THREE.Vector2(), bottomRotation: 0, topRotation: 0 };
        pitch[row][col] = Math.max(1, localStep.length());
      }
    }

    for (let row = 0; row < rows; row += 1) {
      for (let col = 1; col < cols; col += 1) {
        xRows[row][col] = xRows[row][col - 1] + 0.5 * (pitch[row][col - 1] + pitch[row][col]);
      }
    }
    for (let col = 0; col < cols; col += 1) {
      for (let row = 1; row < rows; row += 1) {
        yCols[row][col] = yCols[row - 1][col] + 0.5 * (pitch[row - 1][col] + pitch[row][col]);
      }
    }

    centerPitchCoordinates(xRows, yCols);
    let meanX = 0;
    let meanY = 0;
    for (let row = 0; row < rows; row += 1) {
      for (let col = 0; col < cols; col += 1) {
        meanX += xRows[row][col];
        meanY += yCols[row][col];
      }
    }
    meanX /= Math.max(1, rows * cols);
    meanY /= Math.max(1, rows * cols);

    for (let row = 0; row < rows; row += 1) {
      for (let col = 0; col < cols; col += 1) {
        const cell = cellViews[row][col];
        const entry = snapshot[row][col];
        const evenCell = cell.parity === 0;
        entry.center = new THREE.Vector2(xRows[row][col] - meanX, yCols[row][col] - meanY);
        entry.bottomRotation = (evenCell ? 0 : entry.pose.bLowerRad) + entry.frame.rotation;
        entry.topRotation = (evenCell ? entry.pose.aTopRad : entry.pose.bTopWorldRad) + entry.frame.rotation;
      }
    }

    return snapshot;
  }

  function uniformDriveField(angleDeg) {
    const rows = cellViews.length;
    const cols = rows ? cellViews[0].length : 0;
    return Array.from({ length: rows }, () => Array.from({ length: cols }, () => angleDeg));
  }

  function latticeBodiesFromSnapshot(snapshot) {
    const bodies = [];
    snapshot.forEach((rowEntries, row) => {
      rowEntries.forEach((entry, col) => {
        if (!entry) return;
        const cellLabel = `${row}|${col}`;
        const partLabel = `R${row + 1}C${col + 1}`;
        bodies.push({
          cell: cellLabel,
          part: `${partLabel} lower`,
          origin: entry.center,
          rotation: entry.bottomRotation,
          zMin: -CAD.bodyThicknessMm,
          zMax: 0,
        });
        bodies.push({
          cell: cellLabel,
          part: `${partLabel} upper`,
          origin: entry.center,
          rotation: entry.topRotation,
          zMin: 0,
          zMax: CAD.bodyThicknessMm,
        });
      });
    });
    return bodies;
  }

  function latticeBodiesFromRenderedCells() {
    const bodies = [];
    cellViews.forEach((rowViews) => {
      rowViews.forEach((cell) => {
        const cellLabel = `${cell.row}|${cell.col}`;
        const partLabel = `R${cell.row + 1}C${cell.col + 1}`;
        bodies.push({
          cell: cellLabel,
          part: `${partLabel} lower`,
          origin: cell.center,
          rotation: cell.bottomRotation,
          zMin: -CAD.bodyThicknessMm,
          zMax: 0,
        });
        bodies.push({
          cell: cellLabel,
          part: `${partLabel} upper`,
          origin: cell.center,
          rotation: cell.topRotation,
          zMin: 0,
          zMax: CAD.bodyThicknessMm,
        });
      });
    });
    return bodies;
  }

  function latticeCollisionReportForSnapshot(snapshot) {
    return collisionReportForBodies(latticeBodiesFromSnapshot(snapshot));
  }

  function latticeCollisionReportForDrive(angleDeg, requestedBLower, requestedBTop) {
    if (cellViews.length <= 1 && (cellViews[0]?.length || 0) <= 1) {
      return { maxPenetration: 0, minClearance: 0, collisionCount: 0, collisions: [], clear: true };
    }
    return latticeCollisionReportForSnapshot(buildLatticeSnapshot(uniformDriveField(angleDeg), requestedBLower, requestedBTop));
  }

  function latticeCollisionReportFromRenderedCells() {
    return collisionReportForBodies(latticeBodiesFromRenderedCells());
  }

  function mergeCollisionReports(first, second) {
    return {
      maxPenetration: Math.max(first.maxPenetration || 0, second.maxPenetration || 0),
      minClearance: Math.min(first.minClearance || 0, second.minClearance || 0),
      collisionCount: (first.collisionCount || 0) + (second.collisionCount || 0),
      collisions: [...(first.collisions || []), ...(second.collisions || [])],
      clear: Boolean(first.clear && second.clear),
    };
  }

  function updateLatticeFromPose(pose, displayFrame, driveField = null, requestedBLower = degToRad(B_LOWER_REFERENCE_DEG), requestedBTop = degToRad(B_TOP_REFERENCE_DEG)) {
    const rows = cellViews.length;
    const cols = rows ? cellViews[0].length : 0;
    const local = Array.from({ length: rows }, () => Array.from({ length: cols }, () => ({ pose, frame: displayFrame })));
    const pitch = Array.from({ length: rows }, () => Array.from({ length: cols }, () => Math.max(1, displayFrame.mapPoint(pose.bCenter).sub(displayFrame.mapPoint(new THREE.Vector2(0, 0))).length())));
    const xRows = Array.from({ length: rows }, () => Array.from({ length: cols }, () => 0));
    const yCols = Array.from({ length: rows }, () => Array.from({ length: cols }, () => 0));

    if (driveField) {
      for (let row = 0; row < rows; row += 1) {
        for (let col = 0; col < cols; col += 1) {
          const localSolution = solveNonPenetratingPose(degToRad(driveField[row][col]), requestedBLower, requestedBTop);
          const localFrame = displayFrameForPose(localSolution.pose);
          const localStep = localFrame.mapPoint(localSolution.pose.bCenter).sub(localFrame.mapPoint(new THREE.Vector2(0, 0)));
          local[row][col] = { pose: localSolution.pose, frame: localFrame };
          pitch[row][col] = Math.max(1, localStep.length());
        }
      }
    }

    for (let row = 0; row < rows; row += 1) {
      for (let col = 1; col < cols; col += 1) {
        xRows[row][col] = xRows[row][col - 1] + 0.5 * (pitch[row][col - 1] + pitch[row][col]);
      }
    }
    for (let col = 0; col < cols; col += 1) {
      for (let row = 1; row < rows; row += 1) {
        yCols[row][col] = yCols[row - 1][col] + 0.5 * (pitch[row - 1][col] + pitch[row][col]);
      }
    }

    centerPitchCoordinates(xRows, yCols);
    let meanX = 0;
    let meanY = 0;
    for (let row = 0; row < rows; row += 1) {
      for (let col = 0; col < cols; col += 1) {
        meanX += xRows[row][col];
        meanY += yCols[row][col];
      }
    }
    meanX /= Math.max(1, rows * cols);
    meanY /= Math.max(1, rows * cols);

    const baseCenters = Array.from({ length: rows }, () => Array.from({ length: cols }, () => new THREE.Vector2()));
    for (let row = 0; row < rows; row += 1) {
      for (let col = 0; col < cols; col += 1) {
        baseCenters[row][col].set(xRows[row][col] - meanX, yCols[row][col] - meanY);
      }
    }
    const constrained = constrainedLayoutSolve(baseCenters, pitch);
    latestBoundaryReport = constrained.report;

    let maxRadius = CAD.siteRadiusMm * 2;
    const centerSum = new THREE.Vector2(0, 0);
    let centerCount = 0;

    for (let row = 0; row < rows; row += 1) {
      for (let col = 0; col < cols; col += 1) {
        const cell = cellViews[row][col];
        const center = constrained.centers[row][col];
        const z = constrained.z[row][col];
        const evenCell = cell.parity === 0;
        const localPose = local[row][col].pose;
        const localFrame = local[row][col].frame;
        const commandDeg = driveField ? driveField[row][col] : radToDeg(localPose.aTopRad);
        const bottomRotation = (evenCell ? 0 : localPose.bLowerRad) + localFrame.rotation;
        const topRotation = (evenCell ? localPose.aTopRad : localPose.bTopWorldRad) + localFrame.rotation;
        applyCellColor(cell, commandDeg, cachedDriveEnvelope(requestedBLower, requestedBTop));
        cell.center.copy(center);
        cell.z = z;
        cell.tiltX = constrained.tiltX[row][col];
        cell.tiltY = constrained.tiltY[row][col];
        cell.bottomRotation = bottomRotation;
        cell.topRotation = topRotation;
        cell.group.position.set(center.x, center.y, z);
        const normal = new THREE.Vector3(Math.sin(cell.tiltY),
          -Math.sin(cell.tiltX) * Math.cos(cell.tiltY), Math.cos(cell.tiltX) * Math.cos(cell.tiltY));
        cell.group.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), normal);
        cell.bottom.rotation.z = bottomRotation;
        cell.top.rotation.z = topRotation;
        cell.bottom.scale.z = cell.top.scale.z = plateThicknessMm() / CAD.bodyThicknessMm;
        cell.bottom.position.z = -plateThicknessMm() * 0.5;
        cell.top.position.z = plateThicknessMm() * 0.5;
        maxRadius = Math.max(maxRadius, Math.hypot(center.length(), z) + CAD.cellWidthMm);
        centerSum.add(center);
        centerCount += 1;
      }
    }

    if (driveField && xyLocks.size <= 1) {
      const envelope = cachedDriveEnvelope(requestedBLower, requestedBTop);
      const previousCenters = contactMemory?.freeCenters;
      const previousBase = contactMemory?.freeBaseCenters;
      if (pinHoleClearanceMm() > 0 && previousCenters?.length === rows * cols && previousBase?.length === rows * cols) {
        cellViews.flat().forEach((cell, index) => {
          cell.center.x += previousCenters[index][0] - previousBase[index][0];
          cell.center.y += previousCenters[index][1] - previousBase[index][1];
        });
      }
      const bodies = new Map(cellViews.flat().map((cell, index) => {
        const angle = degToRad(driveField[cell.row][cell.col]);
        const epsilon = 1e-4;
        const next = solveNonPenetratingPose(angle + epsilon, requestedBLower, requestedBTop).pose;
        const frame = displayFrameForPose(next);
        const top = (cell.parity === 0 ? next.aTopRad : next.bTopWorldRad) + frame.rotation;
        const bottom = (cell.parity === 0 ? 0 : next.bLowerRad) + frame.rotation;
        return [cell, { angle, initial: angle, center: cell.center,
          yaw: contactMemory?.freeYaw?.[index] || 0, top: cell.topRotation, bottom: cell.bottomRotation,
          topSlope: angleDelta(top, cell.topRotation) / epsilon,
          bottomSlope: angleDelta(bottom, cell.bottomRotation) / epsilon,
          weight: cell.row === selectedDriveCell.row && cell.col === selectedDriveCell.col ? 0 : 1 / CAD.siteRadiusMm ** 2 }];
      }));
      const pairs = explicitConnectionPinPairs().map(([first, second]) => {
        const a = bodies.get(first.cell), b = bodies.get(second.cell);
        const sa = a[first.layer + "Slope"], sb = b[second.layer + "Slope"];
        const va = SITE_VECTORS[first.site], vb = SITE_VECTORS[second.site];
        return { a, b, sa, sb,
          baseA: Math.atan2(va.y, va.x) + a[first.layer] - sa * a.initial,
          baseB: Math.atan2(vb.y, vb.x) + b[second.layer] - sb * b.initial };
      });
      const clearance = 2 * pinHoleClearanceMm();
      // Free sheets must close every pin pair together, rather than independently
      // laying out rows and columns and treating their mismatch as an actuator stop.
      for (let iteration = 0; iteration < 4000; iteration += 1) {
        let residual = 0;
        for (const { a, b, sa, sb, baseA, baseB } of pairs) {
          const ta = baseA + sa * a.angle + a.yaw, tb = baseB + sb * b.angle + b.yaw;
          const ax = CAD.siteRadiusMm * Math.cos(ta), ay = CAD.siteRadiusMm * Math.sin(ta);
          const bx = CAD.siteRadiusMm * Math.cos(tb), by = CAD.siteRadiusMm * Math.sin(tb);
          const dx = b.center.x + bx - a.center.x - ax, dy = b.center.y + by - a.center.y - ay;
          const distance = Math.hypot(dx, dy);
          const excess = Math.max(0, distance - clearance);
          residual = Math.max(residual, excess);
          if (excess <= 1e-8) continue;
          const nx = dx / distance, ny = dy / distance;
          const ya = nx * ay - ny * ax, yb = -nx * by + ny * bx;
          const ga = ya * sa, gb = yb * sb;
          const step = 1.5 * excess / (2 + a.weight * (ga * ga + ya * ya) + b.weight * (gb * gb + yb * yb));
          a.center.x += nx * step; a.center.y += ny * step;
          b.center.x -= nx * step; b.center.y -= ny * step;
          a.angle = clamp(a.angle - step * a.weight * ga, degToRad(envelope.min), degToRad(envelope.max));
          b.angle = clamp(b.angle - step * b.weight * gb, degToRad(envelope.min), degToRad(envelope.max));
          a.yaw -= step * a.weight * ya;
          b.yaw -= step * b.weight * yb;
        }
        if (residual < 1e-5) break;
      }
      const fitsPairs = () => pairs.every(({ a, b, sa, sb, baseA, baseB }) => {
        const ta = baseA + sa * a.angle + a.yaw, tb = baseB + sb * b.angle + b.yaw;
        return Math.hypot(b.center.x - a.center.x + CAD.siteRadiusMm * (Math.cos(tb) - Math.cos(ta)),
          b.center.y - a.center.y + CAD.siteRadiusMm * (Math.sin(tb) - Math.sin(ta))) <= clearance + 1e-5;
      });
      if (!fitsPairs()) {
        // Numerical recovery changes passive coordinates, never the commanded
        // angle. Uniform motion is a known feasible point for an unlocked sheet.
        const driver = bodies.get(cellAt(selectedDriveCell.row, selectedDriveCell.col));
        const uniform = buildLatticeSnapshot(uniformDriveField(radToDeg(driver.angle)), requestedBLower, requestedBTop);
        const targets = [...bodies.entries()].map(([cell, body]) => ({ body,
          center: body.center.clone(), angle: body.angle, yaw: body.yaw,
          uniformCenter: rotate2(uniform[cell.row][cell.col].center, driver.yaw) }));
        const angle = driver.angle, yaw = driver.yaw;
        const mix = (amount) => targets.forEach((target) => {
          target.body.center.copy(target.uniformCenter).lerp(target.center, amount);
          target.body.angle = angle + (target.angle - angle) * amount;
          target.body.yaw = yaw + (target.yaw - yaw) * amount;
        });
        let accepted = 0, rejected = 1;
        for (let iteration = 0; iteration < 32; iteration += 1) {
          const trial = (accepted + rejected) * 0.5;
          mix(trial);
          if (fitsPairs()) accepted = trial;
          else rejected = trial;
        }
        mix(accepted);
      }
      bodies.forEach((body, cell) => {
        driveField[cell.row][cell.col] = radToDeg(body.angle);
        cell.topRotation = body.top + body.topSlope * (body.angle - body.initial) + body.yaw;
        cell.bottomRotation = body.bottom + body.bottomSlope * (body.angle - body.initial) + body.yaw;
        cell.top.rotation.z = cell.topRotation;
        cell.bottom.rotation.z = cell.bottomRotation;
        cell.group.position.set(cell.center.x, cell.center.y, 0);
        maxRadius = Math.max(maxRadius, cell.center.length() + CAD.cellWidthMm);
        applyCellColor(cell, driveField[cell.row][cell.col], envelope);
      });
      if (latestContactField?.field === driveField) {
        const field = latestContactField;
        field.freeYaw = [...bodies.values()].map((body) => body.yaw);
        field.freeBaseCenters = buildLatticeSnapshot(driveField, requestedBLower, requestedBTop)
          .flat().map((entry) => entry.center.toArray());
        field.maxAbs = field.reach = field.dieOffRadius = 0;
        driveField.forEach((row, r) => row.forEach((value, c) => {
          const delta = value - field.restDrive;
          field.deltaField[r][c] = delta;
          field.maxAbs = Math.max(field.maxAbs, Math.abs(delta));
          if (Math.abs(delta) > 1e-6) {
            field.reach += 1;
            field.dieOffRadius = Math.max(field.dieOffRadius, field.distance[r][c]);
          }
        }));
      }
    }
    if (xyLocks.size <= 1) {
      // A single center anchor changes the world frame, not the internal mechanism.
      if (latestContactField?.field === driveField) {
        latestContactField.freeCenters = cellViews.flat().map((cell) => cell.center.toArray());
      }
      const lock = [...xyLocks.values()][0];
      if (lock) {
        const anchor = cellAt(lock.row, lock.col);
        freeTranslation = [lock.targetX - anchor.center.x, lock.targetY - anchor.center.y,
          lock.lockZ ? lock.targetZ - anchor.z : freeTranslation[2]];
      }
      cellViews.flat().forEach((cell) => {
        cell.center.x += freeTranslation[0];
        cell.center.y += freeTranslation[1];
        cell.z += freeTranslation[2];
        cell.group.position.set(cell.center.x, cell.center.y, cell.z);
      });
      if (latestContactField?.field === driveField) latestContactField.freeTranslation = freeTranslation.slice();
      constrained.report.maxLift = Math.max(...cellViews.flat().map((cell) => Math.abs(cell.z)));
    }
    centerSum.set(0, 0);
    cellViews.flat().forEach((cell) => centerSum.add(cell.center));
    latestLatticeRadius = Math.max(...cellViews.flat().map((cell) => cell.center.length() + CAD.cellWidthMm));
    const systemCenter = centerCount ? centerSum.multiplyScalar(1 / centerCount) : new THREE.Vector2(0, 0);
    buildCenterLineGeometry();
    latticeLine.visible = showCenterLinesInput.checked;
    updateConnectionPins();
    updateSelectionMarker();
    updateLockMarkers();
    return {
      stepX: new THREE.Vector2(pitch[0]?.[0] || 0, 0),
      stepY: new THREE.Vector2(0, pitch[0]?.[0] || 0),
      systemCenter,
      boundaryReport: constrained.report,
    };
  }

  let lockedDriveRangeCache = null;

  function lockedDriveRange(envelope, requestedBLower, requestedBTop) {
    const selectedKey = cellKey(selectedDriveCell.row, selectedDriveCell.col);
    const key = JSON.stringify([rowCountInput.value, colCountInput.value, pinRadiusRatio(), plateThicknessMm(),
      [...xyLocks.entries()], selectedKey, [...cellDriveCommands.entries()].filter(([key]) => key !== selectedKey), contactMemory, envelope]);
    if (lockedDriveRangeCache?.key === key) return lockedDriveRangeCache.range;
    const savedCommand = cellDriveCommands.get(selectedKey);
    const hadCommand = cellDriveCommands.has(selectedKey);
    const rest = contactRestDriveDeg(envelope);
    let anchor = selectedDriveCommand(envelope);
    const fits = (angle) => {
      setDriveCommandForCell(selectedDriveCell.row, selectedDriveCell.col, angle, envelope);
      const field = contactDriveField(envelope);
      latestContactField = field;
      latestRestSnapshot = buildLatticeSnapshot(uniformDriveField(field.restDrive), requestedBLower, requestedBTop);
      const { pose } = solveNonPenetratingPose(degToRad(angle), requestedBLower, requestedBTop);
      updateLatticeFromPose(pose, displayFrameForPose(pose), field.field, requestedBLower, requestedBTop);
      return latestPinFitReport.clear;
    };
    try {
      if (!fits(anchor)) anchor = rest;
      let range = { min: anchor, max: anchor };
      if (fits(anchor)) {
        const bound = (target) => {
          let accepted = anchor;
          const direction = Math.sign(target - anchor);
          while (Math.abs(target - accepted) > 1e-8) {
            const trial = accepted + direction * Math.min(1, Math.abs(target - accepted));
            if (fits(trial)) { accepted = trial; continue; }
            let rejected = trial;
            for (let iteration = 0; iteration < 16; iteration += 1) {
              const mid = (accepted + rejected) * 0.5;
              if (fits(mid)) accepted = mid;
              else rejected = mid;
            }
            break;
          }
          return accepted;
        };
        const min = bound(envelope.min);
        const max = bound(envelope.max);
        if (max >= min) range = { min, max };
      }
      lockedDriveRangeCache = { key, range };
      return range;
    } finally {
      if (hadCommand) cellDriveCommands.set(selectedKey, savedCommand);
      else cellDriveCommands.delete(selectedKey);
      // Trial poses must never be left on screen by the unchanged-state fast path.
      lastMechanismUpdateKey = "";
      // Trial poses are overwritten by the real update before the next render.
      lastMechanismUpdateKey = "";
    }
  }

  function controlStateKey() {
    return JSON.stringify([aTopInput.value, actuationModeInput.value, rowCountInput.value, colCountInput.value,
      pinRadiusRatioInput.value, plateThicknessInput.value, showPinsInput.checked, showCenterLinesInput.checked,
      selectedDriveCell, [...selectedCellKeys], [...xyLocks.entries()], [...cellDriveCommands.entries()]]);
  }

  function updateMechanism(timeMs = 0) {
    if (!animateInput.checked && !selectedDriveCommandDirty && lastMechanismUpdateKey && controlStateKey() === lastControlStateKey) return;
    const contactMode = actuationModeInput.value === "contact";
    pinRadiusRatioInput.disabled = !contactMode;
    if (!contactMode) pinRadiusRatioInput.value = "1";
    const requestedBLower = degToRad(B_LOWER_REFERENCE_DEG);
    const requestedBTop = degToRad(B_TOP_REFERENCE_DEG);
    const driveEnvelope = cachedDriveEnvelope(requestedBLower, requestedBTop);
    aTopInput.step = contactMode ? "any" : "0.1";

    if (driveEnvelope) {
      aTopInput.min = String(driveEnvelope.min);
      aTopInput.max = String(driveEnvelope.max);
    } else {
      aTopInput.min = String(A_TOP_FULL_MIN_DEG);
      aTopInput.max = String(A_TOP_FULL_MAX_DEG);
    }

    const jammedDriveDeg = jammedLockDriveDeg(driveEnvelope);
    const jammed = jammedDriveDeg !== null;
    aTopInput.disabled = jammed;
    animateInput.disabled = jammed;
    animationHalfPeriodSecondsInput.disabled = jammed;
    aTopInput.title = jammed ? "Jammed: position locks and zero pin clearance allow 0 degrees of movement." : "";
    if (jammed) animateInput.checked = false;

    let requestedATopDeg = Number(aTopInput.value);
    if (animateInput.checked && driveEnvelope) {
      const range = widestRange(driveEnvelope);
      const halfPeriodSeconds = Number(animationHalfPeriodSecondsInput.value) || 2.5;
      requestedATopDeg = mildlyEasedPingPong(timeMs, range[0], range[1], halfPeriodSeconds, ANIMATION_EASE_BLEND);
      aTopInput.value = requestedATopDeg.toFixed(2);
    }

    let aTopDeg = requestedATopDeg;
    let clamped = false;
    if (driveEnvelope && !envelopeContainsDrive(aTopDeg, driveEnvelope)) {
      aTopDeg = nearestFeasibleDrive(aTopDeg, driveEnvelope);
      aTopInput.value = aTopDeg.toFixed(0);
      clamped = true;
    }
    if (jammedDriveDeg !== null) {
      aTopDeg = jammedDriveDeg;
      aTopInput.value = aTopDeg.toFixed(2);
      clamped = true;
    }
    const enteringContactMode = contactMode && lastActuationMode !== "contact";
    if (enteringContactMode) {
      seedContactCommandsFromUniformDrive(aTopDeg, driveEnvelope);
    }
    const range = contactMode && xyLocks.size > 1 && !jammed && driveEnvelope
      ? lockedDriveRange(driveEnvelope, requestedBLower, requestedBTop) : null;
    if (range) {
      aTopInput.step = "any";
      aTopInput.min = String(range.min);
      aTopInput.max = String(range.max);
      if (animateInput.checked) {
        requestedATopDeg = mildlyEasedPingPong(timeMs, range.min, range.max, Number(animationHalfPeriodSecondsInput.value) || 2.5, ANIMATION_EASE_BLEND);
      }
      aTopDeg = clamp(requestedATopDeg, range.min, range.max);
      clamped = clamped || Math.abs(aTopDeg - requestedATopDeg) > 1e-9;
      aTopInput.value = String(aTopDeg);
      aTopInput.disabled = animateInput.disabled = animationHalfPeriodSecondsInput.disabled = range.max - range.min < 1e-9;
      aTopInput.title = aTopInput.disabled
        ? "No feasible motion found on the current geometric branch. Check locks and pin clearance."
        : "Admissible range on the current geometric branch, checked against pin-hole clearance.";
      if (aTopInput.disabled) animateInput.checked = false;
    }
    if (contactMode && (selectedDriveCommandDirty || animateInput.checked || (range && clamped))) {
      setDriveCommandForCell(selectedDriveCell.row, selectedDriveCell.col, aTopDeg, driveEnvelope);
    }
    selectedDriveCommandDirty = false;
    lastActuationMode = actuationModeInput.value;

    const mechanismUpdateKey = [
      aTopDeg.toFixed(6),
      rowCountInput.value,
      colCountInput.value,
      actuationModeInput.value,
      selectedDriveCell.row,
      selectedDriveCell.col,
      JSON.stringify([...selectedCellKeys].sort()),
      showPinsInput.checked ? 1 : 0,
      showCenterLinesInput.checked ? 1 : 0,
      pinRadiusRatioInput.value,
      plateThicknessInput.value,
      JSON.stringify([...xyLocks.entries()]),
      JSON.stringify([...cellDriveCommands.entries()]),
      driveEnvelope ? driveEnvelope.min : "none",
      driveEnvelope ? driveEnvelope.max : "none",
    ].join("|");
    if (mechanismUpdateKey === lastMechanismUpdateKey) return;
    lastMechanismUpdateKey = mechanismUpdateKey;

    allowedDriveMetric.textContent = jammed ? "0 deg travel (jammed)" : range ? `${range.min.toFixed(3)} to ${range.max.toFixed(3)} deg (geometry / clearance)` : formatDriveRanges(driveEnvelope);
    driveClampMetric.textContent = jammed ? "position locks / zero clearance" : clamped ? `${requestedATopDeg.toFixed(0)} -> ${aTopDeg.toFixed(0)} deg` : driveEnvelope ? "none" : "no feasible range";

    let solution = solveNonPenetratingPose(degToRad(aTopDeg), requestedBLower, requestedBTop);
    if (clamped && solution.status === "clear") {
      solution = { ...solution, status: "adjusted" };
    }

    const pose = solution.pose;
    const displayFrame = displayFrameForPose(pose);
    let contactField = contactMode ? contactDriveField(driveEnvelope) : null;
    latestContactField = contactField;
    latestRestSnapshot = contactField ? buildLatticeSnapshot(uniformDriveField(contactField.restDrive), requestedBLower, requestedBTop) : null;
    let lattice = updateLatticeFromPose(pose, displayFrame, contactField?.field || null, requestedBLower, requestedBTop);
    if (contactField && latestPinFitReport && !latestPinFitReport.clear) {
      let low = 0;
      let high = 1;
      let best = scaledContactField(contactField, 0);
      latestContactField = best;
      lattice = updateLatticeFromPose(pose, displayFrame, best.field, requestedBLower, requestedBTop);
      for (let iteration = 0; iteration < 10; iteration += 1) {
        const mid = (low + high) * 0.5;
        const candidate = scaledContactField(contactField, mid);
        latestContactField = candidate;
        const candidateLattice = updateLatticeFromPose(pose, displayFrame, candidate.field, requestedBLower, requestedBTop);
        if (latestPinFitReport?.clear) {
          low = mid;
          best = candidate;
          lattice = candidateLattice;
        } else {
          high = mid;
        }
      }
      contactField = best;
      latestContactField = contactField;
      lattice = updateLatticeFromPose(pose, displayFrame, contactField.field, requestedBLower, requestedBTop);
      driveClampMetric.textContent = "pin-hole clearance limit";
    }
    const latticeReport = latticeCollisionReportFromRenderedCells();
    const combinedReport = mergeCollisionReports(solution.report, latticeReport);
    const displayedSolution = {
      ...solution,
      report: combinedReport,
      status: combinedReport.clear ? solution.status : "blocked",
    };
    const aCell = cellAt(0, 0);
    const bCell = cellAt(0, 1);
    const cCell = cellAt(0, 2);
    const displayB = displayFrame.mapPoint(pose.bCenter);
    const displayC = displayFrame.mapPoint(pose.cCenter);
    const aAttachWorld = displayFrame.mapPoint(pose.aAttachWorld);
    const bAttachWorld = displayFrame.mapPoint(pose.bCenter.clone().add(pose.bAttachOffset));
    const aLowerAttachWorld = displayFrame.mapPoint(pose.aLowerAttachWorld);
    const bTopAttachWorld = displayFrame.mapPoint(pose.bTopAttachWorld);
    const residuals = pinResiduals(pose);

    const effectiveDriveDeg = contactField ? contactField.field[selectedDriveCell.row][selectedDriveCell.col] : radToDeg(pose.aTopRad);
    aTopOut.textContent = angleLabel(requestedATopDeg, degToRad(effectiveDriveDeg));
    residualMetric.textContent = `${aAttachWorld.distanceTo(bAttachWorld).toFixed(3)} mm`;
    pitchMetric.textContent = `${lattice.stepX.length().toFixed(1)} mm`;
    aCenterMetric.textContent = aCell ? formatPair(aCell.center) : "n/a";
    bCenterMetric.textContent = bCell ? formatPair(bCell.center) : "n/a";
    cCenterMetric.textContent = cCell ? formatPair(cCell.center) : "n/a";
    systemCenterMetric.textContent = formatPair(lattice.systemCenter);
    secondResidualMetric.textContent = `${aLowerAttachWorld.distanceTo(bTopAttachWorld).toFixed(3)} mm`;
    secondTargetMetric.textContent = formatPair(aLowerAttachWorld);
    secondBMetric.textContent = formatPair(bTopAttachWorld);
    bcPrimaryResidualMetric.textContent = `${residuals.bcPrimary.toFixed(3)} mm`;
    bcSecondResidualMetric.textContent = `${residuals.bcSecond.toFixed(3)} mm`;
    bcPitchMetric.textContent = `${displayC.clone().sub(displayB).length().toFixed(1)} mm`;
    selectedCellMetric.textContent = selectedCellLabel();
    dieOffRadiusMetric.textContent = contactField ? `${contactField.dieOffRadius} steps (${contactField.reach}/${cellViews.length * (cellViews[0]?.length || 0)} cells)` : "uniform";
    maxPropagatedMetric.textContent = contactField ? `${contactField.maxAbs.toFixed(1)} deg` : `${Math.abs(aTopDeg).toFixed(1)} deg`;
    updateCollisionReadout(displayedSolution);
    document.getElementById("pinFitMetric").textContent = latestPinFitReport?.clear ? "within tolerance" : "invalid configuration";
    document.getElementById("bodyThicknessMetric").textContent = `${plateThicknessMm().toFixed(2)} mm each`;
    updateContactReadout(contactField);
    if (contactField && latestPinFitReport?.clear) {
      contactMemory = { restDrive: contactField.restDrive, field: contactField.field.map((row) => row.slice()),
        freeCenters: xyLocks.size <= 1 ? contactField.freeCenters : null,
        freeTranslation: contactField.freeTranslation,
        freeBaseCenters: contactField.freeBaseCenters,
        freeYaw: contactField.freeYaw,
        pins: connectionPinMeshes.map((pin) => pin.position.toArray()) };
      contactField.field.forEach((row, r) => row.forEach((value, c) => cellDriveCommands.set(cellKey(r, c), value)));
    }
    lastControlStateKey = controlStateKey();
    updateBoundaryReadout(lattice.boundaryReport || latestBoundaryReport);
  }

  function continuitySnapshot() {
    return cellViews.flatMap((row) =>
      row.map((cell) => ({
        row: cell.row,
        col: cell.col,
        center: cell.center.clone(),
        topRotation: cell.topRotation,
        bottomRotation: cell.bottomRotation,
      }))
    );
  }

  function maxSnapshotJump(previous, current) {
    return maxSnapshotEntryJump(previous, current);
  }

  function scanStringContinuity(options = {}) {
    const ratios = options.ratios || [0.35, 0.5, 0.72, 0.9, 1.0];
    const stepDeg = Math.max(0.1, Number(options.stepDeg) || 0.5);
    const maxCenterJumpMm = Number(options.maxCenterJumpMm) || MAX_DRIVE_STEP_CENTER_JUMP_MM;
    const maxRotationJumpDeg = Number(options.maxRotationJumpDeg) || MAX_DRIVE_STEP_ROTATION_JUMP_DEG;
    const saved = {
      rows: rowCountInput.value,
      cols: colCountInput.value,
      drive: aTopInput.value,
      mode: actuationModeInput.value,
      ratio: pinRadiusRatioInput.value,
      selected: { ...selectedDriveCell },
      selectedKeys: [...selectedCellKeys],
      commands: [...cellDriveCommands.entries()],
      memory: contactMemory, translation: freeTranslation.slice(),
    };
    const failures = [];

    try {
      [
        [1, 15],
        [15, 1],
      ].forEach(([rows, cols]) => {
        rowCountInput.value = String(rows);
        colCountInput.value = String(cols);
        selectedDriveCell = { row: 0, col: 0 };
        selectedCellKeys.clear();
        selectedCellKeys.add("0,0");
        cellDriveCommands.clear();
        rebuildLattice();
        actuationModeInput.value = "contact";
        ratios.forEach((ratio) => {
          pinRadiusRatioInput.value = String(ratio);
          lastMechanismUpdateKey = "";
          updateMechanism();
          const minDrive = Number(aTopInput.min);
          const maxDrive = Number(aTopInput.max);
          ["descending", "ascending"].forEach((direction) => {
            let previous = null;
            const drives = [];
            if (direction === "descending") {
              for (let drive = maxDrive; drive >= minDrive - 1e-9; drive -= stepDeg) drives.push(Math.max(minDrive, drive));
            } else {
              for (let drive = minDrive; drive <= maxDrive + 1e-9; drive += stepDeg) drives.push(Math.min(maxDrive, drive));
            }
            drives.forEach((drive) => {
              aTopInput.value = String(drive);
              selectedDriveCommandDirty = true;
              lastMechanismUpdateKey = "";
              updateMechanism();
              const current = continuitySnapshot();
              if (previous) {
                const jump = maxSnapshotJump(previous.snapshot, current);
                if (jump.centerJump > maxCenterJumpMm || jump.rotationJumpDeg > maxRotationJumpDeg) {
                  failures.push({
                    rows,
                    cols,
                    ratio,
                    direction,
                    fromDrive: previous.drive,
                    toDrive: drive,
                    centerJump: Number(jump.centerJump.toFixed(3)),
                    rotationJumpDeg: Number(jump.rotationJumpDeg.toFixed(3)),
                  });
                }
              }
              previous = { drive, snapshot: current };
            });
          });
        });
      });
    } finally {
      rowCountInput.value = saved.rows;
      colCountInput.value = saved.cols;
      selectedDriveCell = saved.selected;
      selectedCellKeys.clear();
      saved.selectedKeys.forEach((key) => selectedCellKeys.add(key));
      cellDriveCommands.clear();
      saved.commands.forEach(([key, value]) => cellDriveCommands.set(key, value));
      rebuildLattice();
      contactMemory = saved.memory;
      freeTranslation = saved.translation;
      actuationModeInput.value = saved.mode;
      pinRadiusRatioInput.value = saved.ratio;
      aTopInput.value = saved.drive;
      lastMechanismUpdateKey = "";
      updateMechanism();
    }

    return {
      ok: failures.length === 0,
      failures,
      checked: { strings: ["1x15", "15x1"], ratios, stepDeg, maxCenterJumpMm, maxRotationJumpDeg },
    };
  }

  function scanContinuityScenario(options = {}) {
    const rows = clampInt(options.rows ?? 1, Number(rowCountInput.min), Number(rowCountInput.max));
    const cols = clampInt(options.cols ?? 15, Number(colCountInput.min), Number(colCountInput.max));
    const sourceRow = clampInt(options.sourceRow ?? 0, 0, Math.max(0, rows - 1));
    const sourceCol = clampInt(options.sourceCol ?? 0, 0, Math.max(0, cols - 1));
    const ratio = clamp(Number(options.ratio ?? pinRadiusRatio()), Number(pinRadiusRatioInput.min), Number(pinRadiusRatioInput.max));
    const stepDeg = Math.max(0.1, Number(options.stepDeg) || 0.5);
    const directions = options.direction ? [options.direction] : ["descending", "ascending"];
    const maxCenterJumpMm = Number(options.maxCenterJumpMm) || MAX_DRIVE_STEP_CENTER_JUMP_MM;
    const maxRotationJumpDeg = Number(options.maxRotationJumpDeg) || MAX_DRIVE_STEP_ROTATION_JUMP_DEG;
    const saved = {
      rows: rowCountInput.value,
      cols: colCountInput.value,
      drive: aTopInput.value,
      mode: actuationModeInput.value,
      ratio: pinRadiusRatioInput.value,
      selected: { ...selectedDriveCell },
      selectedKeys: [...selectedCellKeys],
      commands: [...cellDriveCommands.entries()],
      memory: contactMemory, translation: freeTranslation.slice(),
    };
    const failures = [];

    try {
      rowCountInput.value = String(rows);
      colCountInput.value = String(cols);
      selectedDriveCell = { row: sourceRow, col: sourceCol };
      selectedCellKeys.clear();
      selectedCellKeys.add(cellKey(sourceRow, sourceCol));
      cellDriveCommands.clear();
      rebuildLattice();
      actuationModeInput.value = "contact";
      pinRadiusRatioInput.value = String(ratio);
      lastActuationMode = "uniform";
      lastMechanismUpdateKey = "";
      updateMechanism();
      const minDrive = Number(aTopInput.min);
      const maxDrive = Number(aTopInput.max);

      directions.forEach((direction) => {
        let previous = null;
        const drives = [];
        if (direction === "descending") {
          for (let drive = maxDrive; drive >= minDrive - 1e-9; drive -= stepDeg) drives.push(Math.max(minDrive, drive));
        } else {
          for (let drive = minDrive; drive <= maxDrive + 1e-9; drive += stepDeg) drives.push(Math.min(maxDrive, drive));
        }
        drives.forEach((drive) => {
          aTopInput.value = String(drive);
          selectedDriveCommandDirty = true;
          lastMechanismUpdateKey = "";
          updateMechanism();
          const current = continuitySnapshot();
          if (previous) {
            const jump = maxSnapshotJump(previous.snapshot, current);
            if (jump.centerJump > maxCenterJumpMm || jump.rotationJumpDeg > maxRotationJumpDeg) {
              failures.push({
                direction,
                fromDrive: previous.drive,
                toDrive: drive,
                centerJump: Number(jump.centerJump.toFixed(3)),
                rotationJumpDeg: Number(jump.rotationJumpDeg.toFixed(3)),
              });
            }
          }
          previous = { drive, snapshot: current };
        });
      });
    } finally {
      rowCountInput.value = saved.rows;
      colCountInput.value = saved.cols;
      selectedDriveCell = saved.selected;
      selectedCellKeys.clear();
      saved.selectedKeys.forEach((key) => selectedCellKeys.add(key));
      cellDriveCommands.clear();
      saved.commands.forEach(([key, value]) => cellDriveCommands.set(key, value));
      rebuildLattice();
      contactMemory = saved.memory;
      freeTranslation = saved.translation;
      actuationModeInput.value = saved.mode;
      pinRadiusRatioInput.value = saved.ratio;
      aTopInput.value = saved.drive;
      lastActuationMode = saved.mode;
      lastMechanismUpdateKey = "";
      updateMechanism();
    }

    return {
      ok: failures.length === 0,
      failures,
      checked: { rows, cols, sourceRow, sourceCol, ratio, stepDeg, directions, maxCenterJumpMm, maxRotationJumpDeg },
    };
  }

  function scanContinuityScenarios(scenarios = []) {
    const results = scenarios.map((scenario) => ({
      name: scenario.name || `${scenario.rows}x${scenario.cols}`,
      ...scanContinuityScenario(scenario),
    }));
    return {
      ok: results.every((result) => result.ok),
      results,
      failures: results.flatMap((result) => result.failures.map((failure) => ({ scenario: result.name, ...failure }))),
    };
  }

  function lockedStringBacklashScenario(options = {}) {
    const ratio = clamp(Number(options.ratio ?? 1), Number(pinRadiusRatioInput.min), Number(pinRadiusRatioInput.max));
    const saved = {
      rows: rowCountInput.value,
      cols: colCountInput.value,
      drive: aTopInput.value,
      mode: actuationModeInput.value,
      ratio: pinRadiusRatioInput.value,
      selected: { ...selectedDriveCell },
      selectedKeys: [...selectedCellKeys],
      commands: [...cellDriveCommands.entries()],
      memory: contactMemory, translation: freeTranslation.slice(),
      locks: [...xyLocks.entries()],
      lastMode: lastActuationMode,
      showPins: showPinsInput.checked,
    };

    try {
      const count = clampInt(options.count ?? 10, 1, 30);
      const rows = clampInt(options.rows ?? (options.vertical ? count : 1), 1, 30);
      const cols = clampInt(options.cols ?? (options.vertical ? 1 : count), 1, 30);
      rowCountInput.value = String(rows);
      colCountInput.value = String(cols);
      showPinsInput.checked = options.showPins !== false;
      pinRadiusRatioInput.value = String(ratio);
      selectedDriveCell = { row: options.sourceRow ?? Math.floor(rows / 2), col: options.sourceCol ?? Math.floor(cols / 2) };
      selectedCellKeys.clear();
      selectedCellKeys.add(cellKey(selectedDriveCell.row, selectedDriveCell.col));
      cellDriveCommands.clear();
      xyLocks.clear();
      rebuildLattice();
      actuationModeInput.value = "uniform";
      lastActuationMode = "uniform";
      lastMechanismUpdateKey = "";
      updateMechanism();
      aTopInput.value = aTopInput.min;
      lastMechanismUpdateKey = "";
      updateMechanism();

      const lockSites = options.locks || (rows > 1 && cols > 1
        ? [[0, 0], [0, cols - 1], [rows - 1, 0], [rows - 1, cols - 1]]
        : [[0, 0], [rows - 1, cols - 1]]);
      lockSites.map(([row, col]) => cellAt(row, col)).forEach((cell) => {
        xyLocks.set(cellKey(cell.row, cell.col), {
          row: cell.row,
          col: cell.col,
          targetX: cell.center.x + (options.lockOffset?.[0] || 0),
          targetY: cell.center.y + (options.lockOffset?.[1] || 0),
          targetZ: (cell.z || 0) + (options.lockOffset?.[2] || 0),
          targetDriveDeg: Number(aTopInput.value),
          lockXY: true,
          lockZ: options.lockModes?.[cellKey(cell.row, cell.col)] !== "xy" && options.lockMode !== "xy",
        });
      });

      let uniformEffectiveDriveDeltaDeg = null;
      if (exactPinTiltLimitRad() <= 1e-9) {
        actuationModeInput.value = "uniform";
        aTopInput.value = aTopInput.max;
        lastMechanismUpdateKey = "";
        updateMechanism();
        uniformEffectiveDriveDeltaDeg = Math.abs(Number(aTopInput.value) - Number(aTopInput.min));
      }

      actuationModeInput.value = "contact";
      pinRadiusRatioInput.value = String(ratio);
      aTopInput.value = aTopInput.min;
      lastActuationMode = "uniform";
      lastMechanismUpdateKey = "";
      updateMechanism();
      const samples = [];
      for (const fraction of options.driveOffsets || options.driveFractions || [1]) {
        aTopInput.value = String(Number(aTopInput.min) + (options.driveOffsets ? fraction : fraction * (Number(aTopInput.max) - Number(aTopInput.min))));
        selectedDriveCommandDirty = true;
        lastMechanismUpdateKey = "";
        updateMechanism();
        if (options.driveFractions || options.driveOffsets) samples.push({
          fraction,
          pinFit: { ...latestPinFitReport },
          maxLift: latestBoundaryReport.maxLift,
          maxLockDrift: latestBoundaryReport.maxLockDrift,
          maxEdgeLengthError: latestBoundaryReport.maxEdgeLengthError,
          slider: { min: Number(aTopInput.min), max: Number(aTopInput.max), value: Number(aTopInput.value) },
          effectiveDrive: latestContactField.field[selectedDriveCell.row][selectedDriveCell.col],
          cells: cellViews.flat().map((cell) => ({ row: cell.row, col: cell.col,
            center: [cell.center.x, cell.center.y, cell.z], normal: cellNormal(cell).toArray(),
            delta: latestContactField.deltaField[cell.row][cell.col],
            lock: lockForCell(cell.row, cell.col) || null })),
          pins: explicitConnectionPinPairs().map((pair, index) => {
            const pin = connectionPinMeshes[index];
            return {
              origin: pin.position.toArray(),
              axis: new THREE.Vector3(0, 1, 0).applyQuaternion(pin.quaternion).toArray(),
              radius: pinRadiusMm(),
              halfLength: pin.geometry.parameters.height * 0.5,
              holes: pair.map((ref) => {
                const site = SITE_VECTORS[ref.site];
                ref.cell.group.updateMatrixWorld(true);
                const layer = ref.layer === "top" ? ref.cell.top : ref.cell.bottom;
                return {
                  center: layer.localToWorld(new THREE.Vector3(site.x, site.y, 0)).toArray(),
                  normal: cellNormal(ref.cell).toArray(),
                  radius: holeRadiusMm(),
                  thickness: plateThicknessMm(),
                };
              }),
            };
          }),
        });
      }

      return {
        ratio,
        rows,
        cols,
        diametralBacklashMm: diametralBacklashMm(),
        exactTiltDeg: exactPinTiltLimitDeg(),
        boundaryActuationScale: boundaryActuationScale(),
        maxLift: latestBoundaryReport?.maxLift || 0,
        maxTiltDeg: latestBoundaryReport?.maxTiltDeg || 0,
        maxAdmissibleEdgeRise: latestBoundaryReport?.maxAdmissibleEdgeRise || 0,
        backlashLimitedLift: Boolean(latestBoundaryReport?.backlashLimitedLift),
        maxEffectiveDriveDeg: latestContactField?.maxAbs || 0,
        contactConstraintScale: latestContactField?.constraintScale ?? 1,
        uniformEffectiveDriveDeltaDeg,
        driveDisabled: aTopInput.disabled,
        animationDisabled: animateInput.disabled,
        slider: { min: Number(aTopInput.min), max: Number(aTopInput.max), value: Number(aTopInput.value) },
        pinFit: { ...(latestPinFitReport || {}) },
        samples,
      };
    } finally {
      rowCountInput.value = saved.rows;
      colCountInput.value = saved.cols;
      selectedDriveCell = saved.selected;
      selectedCellKeys.clear();
      saved.selectedKeys.forEach((key) => selectedCellKeys.add(key));
      cellDriveCommands.clear();
      saved.commands.forEach(([key, value]) => cellDriveCommands.set(key, value));
      xyLocks.clear();
      saved.locks.forEach(([key, value]) => xyLocks.set(key, value));
      rebuildLattice();
      contactMemory = saved.memory;
      freeTranslation = saved.translation;
      actuationModeInput.value = saved.mode;
      pinRadiusRatioInput.value = saved.ratio;
      aTopInput.value = saved.drive;
      lastActuationMode = saved.lastMode;
      showPinsInput.checked = saved.showPins;
      lastMechanismUpdateKey = "";
      updateMechanism();
    }
  }

  function pinFitSnapshot(options = {}) {
    const rows = clampInt(options.rows ?? rowCountInput.value, Number(rowCountInput.min), Number(rowCountInput.max));
    const cols = clampInt(options.cols ?? colCountInput.value, Number(colCountInput.min), Number(colCountInput.max));
    const ratio = clamp(Number(options.ratio ?? pinRadiusRatio()), Number(pinRadiusRatioInput.min), Number(pinRadiusRatioInput.max));
    const drive = options.drive;
    const saved = {
      rows: rowCountInput.value,
      cols: colCountInput.value,
      drive: aTopInput.value,
      mode: actuationModeInput.value,
      ratio: pinRadiusRatioInput.value,
      selected: { ...selectedDriveCell },
      selectedKeys: [...selectedCellKeys],
      commands: [...cellDriveCommands.entries()],
      locks: [...xyLocks.entries()],
      memory: contactMemory, translation: freeTranslation.slice(),
      lastMode: lastActuationMode,
    };
    try {
      rowCountInput.value = String(rows);
      colCountInput.value = String(cols);
      pinRadiusRatioInput.value = String(ratio);
      rebuildLattice();
      if (drive !== undefined) aTopInput.value = String(drive);
      lastMechanismUpdateKey = "";
      updateMechanism();
      return {
        rows,
        cols,
        ratio: pinRadiusRatio(),
        ...(latestPinFitReport || {}),
      };
    } finally {
      rowCountInput.value = saved.rows;
      colCountInput.value = saved.cols;
      pinRadiusRatioInput.value = saved.ratio;
      aTopInput.value = saved.drive;
      actuationModeInput.value = saved.mode;
      selectedDriveCell = saved.selected;
      selectedCellKeys.clear();
      saved.selectedKeys.forEach((key) => selectedCellKeys.add(key));
      cellDriveCommands.clear();
      saved.commands.forEach(([key, value]) => cellDriveCommands.set(key, value));
      xyLocks.clear();
      saved.locks.forEach(([key, value]) => xyLocks.set(key, value));
      rebuildLattice();
      contactMemory = saved.memory;
      freeTranslation = saved.translation;
      lastActuationMode = saved.lastMode;
      lastMechanismUpdateKey = "";
      updateMechanism();
    }
  }

  window.RADAttachmentDiagnostics = {
    currentState() {
      const rect = renderer.domElement.getBoundingClientRect();
      return {
        selected: { ...selectedDriveCell },
        field: contactMemory?.field.map((row) => row.slice()) || null,
        pinFit: { ...latestPinFitReport },
        pins: connectionPinMeshes.map((pin) => pin.position.toArray()),
        cells: cellViews.flat().map((cell) => {
          const projected = cell.group.position.clone().project(camera);
          return { row: cell.row, col: cell.col, center: cell.group.position.toArray(),
            normal: cellNormal(cell).toArray(), thickness: cell.top.scale.z * CAD.bodyThicknessMm,
            topZ: cell.top.position.z, lock: lockForCell(cell.row, cell.col) || null,
            screen: [rect.left + (projected.x + 1) * rect.width / 2, rect.top + (1 - projected.y) * rect.height / 2] };
        }),
      };
    },
    scanStringContinuity,
    scanContinuityScenario,
    scanContinuityScenarios,
    lockedStringBacklashScenario,
    pinFitSnapshot,
  };

  function lockSelectedCellXY() {
    selectedCells().forEach((cell) => {
      xyLocks.set(cellKey(cell.row, cell.col), {
        row: cell.row,
        col: cell.col,
        targetX: cell.center.x,
        targetY: cell.center.y,
        targetZ: cell.z || 0,
        targetDriveDeg: Number(aTopInput.value),
        lockXY: true,
        lockZ: positionLockModeInput.value !== "xy",
      });
    });
    lastMechanismUpdateKey = "";
    updateMechanism();
  }

  function unlockSelectedCellXY() {
    selectedCells().forEach((cell) => {
      xyLocks.delete(cellKey(cell.row, cell.col));
    });
    lastMechanismUpdateKey = "";
    updateMechanism();
  }

  function clearXYLocks() {
    xyLocks.clear();
    lastMechanismUpdateKey = "";
    updateMechanism();
  }

  function resize() {
    const rect = mount.getBoundingClientRect();
    const width = Math.max(1, Math.floor(rect.width));
    const height = Math.max(1, Math.floor(rect.height));
    renderer.setSize(width, height, false);
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
    fitCurrentView();
  }

  let dragging = false;
  let dragMode = "rotate";
  let lastX = 0;
  let lastY = 0;
  let pointerDownX = 0;
  let pointerDownY = 0;
  let pointerDragDistance = 0;

  function selectCellFromPointer(event) {
    const rect = renderer.domElement.getBoundingClientRect();
    pointer.x = ((event.clientX - rect.left) / Math.max(1, rect.width)) * 2 - 1;
    pointer.y = -(((event.clientY - rect.top) / Math.max(1, rect.height)) * 2 - 1);
    raycaster.setFromCamera(pointer, camera);
    const hit = raycaster.intersectObjects(selectableMeshes, false)[0];
    if (!hit || hit.object.userData.cellSelectable !== true) return;
    const row = hit.object.userData.cellRow;
    const col = hit.object.userData.cellCol;
    const key = cellKey(row, col);
    selectedDriveCell = {
      row,
      col,
    };
    if (multiSelectModeInput.checked || event.ctrlKey || event.metaKey || event.shiftKey) {
      if (selectedCellKeys.has(key) && selectedCellKeys.size > 1) {
        selectedCellKeys.delete(key);
      } else {
        selectedCellKeys.add(key);
      }
      if (!selectedCellKeys.size) selectedCellKeys.add(key);
    } else {
      selectedCellKeys.clear();
      selectedCellKeys.add(key);
    }
    if (actuationModeInput.value === "contact") {
      const envelope = cachedDriveEnvelope(degToRad(B_LOWER_REFERENCE_DEG), degToRad(B_TOP_REFERENCE_DEG));
      aTopInput.min = String(envelope?.min ?? A_TOP_FULL_MIN_DEG);
      aTopInput.max = String(envelope?.max ?? A_TOP_FULL_MAX_DEG);
      aTopInput.value = String(selectedDriveCommand(envelope));
    }
    lastMechanismUpdateKey = "";
    updateSelectionMarker();
    selectedCellMetric.textContent = selectedCellLabel();
    updateBoundaryReadout();
  }

  renderer.domElement.addEventListener("pointerdown", (event) => {
    event.preventDefault();
    dragging = true;
    dragMode = event.button === 2 || event.shiftKey ? "pan" : "rotate";
    lastX = event.clientX;
    lastY = event.clientY;
    pointerDownX = event.clientX;
    pointerDownY = event.clientY;
    pointerDragDistance = 0;
    renderer.domElement.setPointerCapture(event.pointerId);
  });

  renderer.domElement.addEventListener("pointermove", (event) => {
    if (!dragging) return;
    const dx = event.clientX - lastX;
    const dy = event.clientY - lastY;
    pointerDragDistance += Math.hypot(dx, dy);
    lastX = event.clientX;
    lastY = event.clientY;
    if (dragMode === "pan") {
      const forward = new THREE.Vector3().subVectors(cameraState.target, camera.position).normalize();
      const right = new THREE.Vector3().crossVectors(forward, camera.up).normalize();
      const up = new THREE.Vector3().crossVectors(right, forward).normalize();
      const scale = cameraState.radius * 0.0018;
      cameraState.target.add(right.multiplyScalar(-dx * scale)).add(up.multiplyScalar(dy * scale));
    } else {
      cameraState.azimuth -= dx * 0.006;
      cameraState.elevation = Math.max(-1.15, Math.min(1.45, cameraState.elevation + dy * 0.006));
    }
    updateCamera();
  });

  renderer.domElement.addEventListener("pointerup", (event) => {
    dragging = false;
    if (event.button === 0 && pointerDragDistance < 4 && Math.hypot(event.clientX - pointerDownX, event.clientY - pointerDownY) < 4) {
      selectCellFromPointer(event);
    }
    try {
      renderer.domElement.releasePointerCapture(event.pointerId);
    } catch (_error) {
      // Pointer capture can already be released by the browser.
    }
  });

  renderer.domElement.addEventListener(
    "wheel",
    (event) => {
      event.preventDefault();
      cameraState.radius = clamp(cameraState.radius * (event.deltaY > 0 ? 1.08 : 0.92), minCameraRadius(), maxCameraRadius());
      updateCamera();
    },
    { passive: false }
  );
  renderer.domElement.addEventListener("contextmenu", (event) => event.preventDefault());

  document.querySelectorAll("[data-view]").forEach((button) => {
    button.addEventListener("click", () => setView(button.dataset.view));
  });
  document.getElementById("fitView").addEventListener("click", fitCurrentView);
  document.getElementById("resetView").addEventListener("click", () => setView("iso"));
  lockSelectedXYButton.addEventListener("click", lockSelectedCellXY);
  unlockSelectedXYButton.addEventListener("click", unlockSelectedCellXY);
  clearXYLocksButton.addEventListener("click", clearXYLocks);
  [
    aTopInput,
    actuationModeInput,
    animateInput,
    animationHalfPeriodSecondsInput,
    showPinsInput,
    showCenterLinesInput,
    pinRadiusRatioInput,
    plateThicknessInput,
    multiSelectModeInput,
  ].forEach((input) => {
    input.addEventListener("input", () => {
      if (input === aTopInput) selectedDriveCommandDirty = true;
      lastMechanismUpdateKey = "";
      updateMechanism();
    });
    input.addEventListener("change", () => {
      if (input === aTopInput) selectedDriveCommandDirty = true;
      lastMechanismUpdateKey = "";
      updateMechanism();
    });
  });
  applyGridButton.addEventListener("click", () => {
    rebuildLattice();
    updateMechanism();
    setView(currentView);
  });
  [rowCountInput, colCountInput].forEach((input) => {
    input.addEventListener("change", () => {
      rebuildLattice();
      updateMechanism();
      setView(currentView);
    });
  });

  new ResizeObserver(resize).observe(mount);
  resize();
  rebuildLattice();
  updateMechanism();
  setView("iso");

  function render(timeMs) {
    updateMechanism(timeMs);
    renderer.render(scene, camera);
    requestAnimationFrame(render);
  }

  requestAnimationFrame(render);
})();
