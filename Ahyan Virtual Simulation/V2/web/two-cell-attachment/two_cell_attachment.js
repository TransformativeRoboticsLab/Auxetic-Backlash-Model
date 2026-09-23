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
  const contactTestCommandInput = document.getElementById("contactTestCommand");
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
  const pinRadiusRatioOut = document.getElementById("pinRadiusRatioOut");
  const contactTestCommandOut = document.getElementById("contactTestCommandOut");
  const holeRadiusMetric = document.getElementById("holeRadiusMetric");
  const pinRadiusMetric = document.getElementById("pinRadiusMetric");
  const radialClearanceMetric = document.getElementById("radialClearanceMetric");
  const freeAngleMetric = document.getElementById("freeAngleMetric");
  const relativeBacklashMetric = document.getElementById("relativeBacklashMetric");
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
  let renderedPinRadius = null;
  let latestContactField = null;
  let latestRestSnapshot = null;
  let selectedDriveCell = { row: 0, col: 0 };
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
    return `r${selectedDriveCell.row + 1} c${selectedDriveCell.col + 1}`;
  }

  function holeRadiusMm() {
    return CAD.nominalHoleDiameterMm * 0.5;
  }

  function pinRadiusRatio() {
    const value = Number(pinRadiusRatioInput.value);
    return clamp(Number.isFinite(value) ? value : 0.72, Number(pinRadiusRatioInput.min), Number(pinRadiusRatioInput.max));
  }

  function pinRadiusMm() {
    return holeRadiusMm() * pinRadiusRatio();
  }

  function pinHoleClearanceMm() {
    return Math.max(0, holeRadiusMm() - pinRadiusMm());
  }

  function normalizedBacklash() {
    return clamp(pinHoleClearanceMm() / CAD.siteRadiusMm, 0, 1);
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

  function contactRestDriveDeg(envelope) {
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

  function contactDriveField(requestedDeg, envelope) {
    const rows = cellViews.length;
    const cols = rows ? cellViews[0].length : 0;
    const restDrive = contactRestDriveDeg(envelope);
    const field = Array.from({ length: rows }, () => Array.from({ length: cols }, () => restDrive));
    const deltaField = Array.from({ length: rows }, () => Array.from({ length: cols }, () => 0));
    const distance = Array.from({ length: rows }, () => Array.from({ length: cols }, () => Infinity));
    const sourceRow = clampInt(selectedDriveCell.row, 0, Math.max(0, rows - 1));
    const sourceCol = clampInt(selectedDriveCell.col, 0, Math.max(0, cols - 1));
    const clampDrive = (value) => (envelope && !envelopeContainsDrive(value, envelope) ? nearestFeasibleDrive(value, envelope) : clamp(value, A_TOP_FULL_MIN_DEG, A_TOP_FULL_MAX_DEG));
    const sourceDrive = clampDrive(requestedDeg);
    const queue = [[sourceRow, sourceCol, sourceDrive - restDrive, 0]];

    while (queue.length) {
      const [row, col, delta, depth] = queue.shift();
      if (row < 0 || row >= rows || col < 0 || col >= cols) continue;
      if (Math.abs(delta) <= Math.abs(deltaField[row][col]) + 1e-9 && depth >= distance[row][col]) continue;
      deltaField[row][col] = delta;
      field[row][col] = clampDrive(restDrive + delta);
      distance[row][col] = Math.min(distance[row][col], depth);
      const residual = contactTransmit(delta);
      if (Math.abs(residual) <= 1e-6) continue;
      neighborCells(row, col, rows, cols).forEach(([nextRow, nextCol]) => {
        if (Math.abs(residual) > Math.abs(deltaField[nextRow][nextCol]) + 1e-9 || depth + 1 < distance[nextRow][nextCol]) {
          queue.push([nextRow, nextCol, residual, depth + 1]);
        }
      });
    }

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
    const pad = annularDisk(CAD.padRadiusMm, holeRadiusMm(), CAD.bodyThicknessMm + 0.05, material, 44);
    pad.position.set(x, y, 0);
    parent.add(pad);
    addEdges(parent, pad);
  }

  function createCrossPart(name, material) {
    const group = new THREE.Group();
    group.name = name;
    addArm(group, 0, material);
    addArm(group, Math.PI / 2, material);
    const hub = annularDisk(CAD.hubRadiusMm, holeRadiusMm(), CAD.bodyThicknessMm + 0.08, material, 48);
    group.add(hub);
    addEdges(group, hub);
    addPad(group, CAD.siteRadiusMm, 0, material);
    addPad(group, -CAD.siteRadiusMm, 0, material);
    addPad(group, 0, CAD.siteRadiusMm, material);
    addPad(group, 0, -CAD.siteRadiusMm, material);
    return group;
  }

  function materialPairForCell() {
    return { top: materials.upper, bottom: materials.lower };
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
      center: new THREE.Vector2(),
      topRotation: 0,
      bottomRotation: 0,
    };
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
    const rows = clampInt(rowCountInput.value, Number(rowCountInput.min), Number(rowCountInput.max));
    const cols = clampInt(colCountInput.value, Number(colCountInput.min), Number(colCountInput.max));
    rowCountInput.value = String(rows);
    colCountInput.value = String(cols);
    resizeFloorGridForDimensions(rows, cols);
    lastMechanismUpdateKey = "";
    removeChildren(latticeRoot);
    removeChildren(pinRoot);
    cellViews = [];
    connections = [];
    connectionPinMeshes = [];
    selectableMeshes.length = 0;
    selectedDriveCell.row = Math.min(selectedDriveCell.row, rows - 1);
    selectedDriveCell.col = Math.min(selectedDriveCell.col, cols - 1);

    for (let row = 0; row < rows; row += 1) {
      const rowViews = [];
      for (let col = 0; col < cols; col += 1) {
        const cell = createCell(row, col);
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
        [makeHoleRef(from, "top", "north"), makeHoleRef(to, "bottom", "west")],
        [makeHoleRef(from, "bottom", "west"), makeHoleRef(to, "top", "west")],
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

  function buildCenterLineGeometry() {
    const points = [];
    for (let row = 0; row < cellViews.length; row += 1) {
      for (let col = 0; col < cellViews[row].length; col += 1) {
        const cell = cellViews[row][col];
        if (col + 1 < cellViews[row].length) {
          const next = cellViews[row][col + 1];
          points.push(new THREE.Vector3(cell.center.x, cell.center.y, 0.3), new THREE.Vector3(next.center.x, next.center.y, 0.3));
        }
        if (row + 1 < cellViews.length) {
          const next = cellViews[row + 1][col];
          points.push(new THREE.Vector3(cell.center.x, cell.center.y, 0.3), new THREE.Vector3(next.center.x, next.center.y, 0.3));
        }
      }
    }
    latticeLine.geometry.dispose();
    latticeLine.geometry = new THREE.BufferGeometry().setFromPoints(points);
  }

  function updateConnectionPins() {
    const visible = showPinsInput.checked;
    const contactMode = actuationModeInput.value === "contact" && latestContactField;
    const drivenPoints = contactMode ? contactPinPoints(latestContactField) : explicitConnectionPinPoints();
    const pinPoints = visible ? dedupePinPoints(contactMode ? drivenPoints : [...drivenPoints, ...overlappingHolePinPoints()]) : [];
    ensureConnectionPinMeshCount(pinPoints.length);
    updatePinGeometry();
    connectionPinMeshes.forEach((pin, index) => {
      const point = pinPoints[index];
      pin.visible = Boolean(visible && point);
      if (point) pin.position.set(point.x, point.y, 0);
    });
  }

  function updateContactReadout() {
    const ratio = pinRadiusRatio();
    const holeRadius = holeRadiusMm();
    const pinRadius = pinRadiusMm();
    const clearance = pinHoleClearanceMm();
    const normalized = normalizedBacklash();
    const backlash = researchAngularBacklashDeg();
    const command = Number(contactTestCommandInput.value) || 0;
    const transmitted = Math.sign(command) * Math.max(0, Math.abs(command) - backlash);
    const floating = Math.abs(command) <= backlash + 1e-9;

    pinRadiusRatioInput.value = ratio.toFixed(2);
    pinRadiusRatioOut.textContent = ratio.toFixed(2);
    contactTestCommandOut.textContent = `${command.toFixed(1)} deg`;
    holeRadiusMetric.textContent = `${holeRadius.toFixed(3)} mm`;
    pinRadiusMetric.textContent = `${pinRadius.toFixed(3)} mm`;
    radialClearanceMetric.textContent = `${clearance.toFixed(3)} mm`;
    freeAngleMetric.textContent = `${backlash.toFixed(2)} deg`;
    relativeBacklashMetric.textContent = normalized.toFixed(4);
    contactStateMetric.textContent = floating ? "floating in clearance" : "contact engaged";
    transmittedCommandMetric.textContent = `${transmitted.toFixed(1)} deg`;
  }

  function updateSelectionMarker() {
    const cell = cellAt(selectedDriveCell.row, selectedDriveCell.col);
    if (!cell || actuationModeInput.value !== "contact") {
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
    selectionRing.position.set(cell.center.x, cell.center.y, CAD.bodyThicknessMm * 1.9);
    selectionRing.visible = true;
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
    const aspect = Math.max(0.45, camera.aspect || 1);
    const aspectPenalty = aspect < 1 ? 1 / aspect : 1;
    return Math.max(180, latestLatticeRadius * (1.18 + aspectPenalty * 0.32) + 80);
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
    return `${effectiveDeg.toFixed(0)} deg`;
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

    let maxRadius = CAD.siteRadiusMm * 2;
    const centerSum = new THREE.Vector2(0, 0);
    let centerCount = 0;

    for (let row = 0; row < rows; row += 1) {
      for (let col = 0; col < cols; col += 1) {
        const cell = cellViews[row][col];
        const center = new THREE.Vector2(xRows[row][col] - meanX, yCols[row][col] - meanY);
        const evenCell = cell.parity === 0;
        const localPose = local[row][col].pose;
        const localFrame = local[row][col].frame;
        const bottomRotation = (evenCell ? 0 : localPose.bLowerRad) + localFrame.rotation;
        const topRotation = (evenCell ? localPose.aTopRad : localPose.bTopWorldRad) + localFrame.rotation;
        cell.center.copy(center);
        cell.bottomRotation = bottomRotation;
        cell.topRotation = topRotation;
        cell.group.position.set(center.x, center.y, 0);
        cell.bottom.rotation.z = bottomRotation;
        cell.top.rotation.z = topRotation;
        maxRadius = Math.max(maxRadius, center.length() + CAD.cellWidthMm);
        centerSum.add(center);
        centerCount += 1;
      }
    }

    latestLatticeRadius = maxRadius;
    const systemCenter = centerCount ? centerSum.multiplyScalar(1 / centerCount) : new THREE.Vector2(0, 0);
    buildCenterLineGeometry();
    latticeLine.visible = showCenterLinesInput.checked;
    updateConnectionPins();
    updateSelectionMarker();
    return {
      stepX: new THREE.Vector2(pitch[0]?.[0] || 0, 0),
      stepY: new THREE.Vector2(0, pitch[0]?.[0] || 0),
      systemCenter,
    };
  }

  function updateMechanism(timeMs = 0) {
    const requestedBLower = degToRad(B_LOWER_REFERENCE_DEG);
    const requestedBTop = degToRad(B_TOP_REFERENCE_DEG);
    const driveEnvelope = cachedDriveEnvelope(requestedBLower, requestedBTop);
    const contactMode = actuationModeInput.value === "contact";

    if (driveEnvelope) {
      aTopInput.min = String(driveEnvelope.min);
      aTopInput.max = String(driveEnvelope.max);
    } else {
      aTopInput.min = String(A_TOP_FULL_MIN_DEG);
      aTopInput.max = String(A_TOP_FULL_MAX_DEG);
    }

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

    const mechanismUpdateKey = [
      aTopDeg.toFixed(2),
      rowCountInput.value,
      colCountInput.value,
      actuationModeInput.value,
      selectedDriveCell.row,
      selectedDriveCell.col,
      showPinsInput.checked ? 1 : 0,
      showCenterLinesInput.checked ? 1 : 0,
      pinRadiusRatioInput.value,
      contactTestCommandInput.value,
      driveEnvelope ? driveEnvelope.min : "none",
      driveEnvelope ? driveEnvelope.max : "none",
    ].join("|");
    if (mechanismUpdateKey === lastMechanismUpdateKey) return;
    lastMechanismUpdateKey = mechanismUpdateKey;

    allowedDriveMetric.textContent = formatDriveRanges(driveEnvelope);
    driveClampMetric.textContent = clamped ? `${requestedATopDeg.toFixed(0)} -> ${aTopDeg.toFixed(0)} deg` : driveEnvelope ? "none" : "no feasible range";

    let solution = solveNonPenetratingPose(degToRad(aTopDeg), requestedBLower, requestedBTop);
    if (clamped && solution.status === "clear") {
      solution = { ...solution, status: "adjusted" };
    }

    const pose = solution.pose;
    const displayFrame = displayFrameForPose(pose);
    const contactField = contactMode ? contactDriveField(aTopDeg, driveEnvelope) : null;
    latestContactField = contactField;
    latestRestSnapshot = contactField ? buildLatticeSnapshot(uniformDriveField(contactField.restDrive), requestedBLower, requestedBTop) : null;
    const lattice = updateLatticeFromPose(pose, displayFrame, contactField?.field || null, requestedBLower, requestedBTop);
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

    aTopOut.textContent = angleLabel(requestedATopDeg, pose.aTopRad);
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
    updateContactReadout();
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
      rebuildLattice();
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

  window.RADAttachmentDiagnostics = {
    scanStringContinuity,
  };

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
    selectedDriveCell = {
      row: hit.object.userData.cellRow,
      col: hit.object.userData.cellCol,
    };
    actuationModeInput.value = "contact";
    lastMechanismUpdateKey = "";
    updateMechanism();
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
  [
    aTopInput,
    actuationModeInput,
    animateInput,
    animationHalfPeriodSecondsInput,
    showPinsInput,
    showCenterLinesInput,
    pinRadiusRatioInput,
    contactTestCommandInput,
  ].forEach((input) => {
    input.addEventListener("input", () => {
      lastMechanismUpdateKey = "";
      updateMechanism();
    });
    input.addEventListener("change", () => {
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
