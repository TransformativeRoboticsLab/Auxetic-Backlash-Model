(function () {
  "use strict";

  const THREE = window.THREE;
  const mount = document.getElementById("threeMount");
  const aTopInput = document.getElementById("aTopAngle");
  const aSiteInput = document.getElementById("aSite");
  const bSiteInput = document.getElementById("bSite");
  const bTopSiteInput = document.getElementById("bTopSite");
  const aLowerSiteInput = document.getElementById("aLowerSite");
  const animateInput = document.getElementById("animate");
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

  const SITE_VECTORS = Object.freeze({
    east: new THREE.Vector2(CAD.siteRadiusMm, 0),
    north: new THREE.Vector2(0, CAD.siteRadiusMm),
    west: new THREE.Vector2(-CAD.siteRadiusMm, 0),
    south: new THREE.Vector2(0, -CAD.siteRadiusMm),
  });

  if (!THREE || !mount) {
    if (mount) mount.textContent = "Three.js did not load.";
    return;
  }

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(42, 1, 0.1, 800);
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
  let centerlineFrame = { key: "", angle: null };

  const materials = {
    aTop: new THREE.MeshStandardMaterial({ color: 0xda7b27, roughness: 0.48, metalness: 0.15 }),
    aBottom: new THREE.MeshStandardMaterial({ color: 0x44515f, roughness: 0.55, metalness: 0.12 }),
    bTop: new THREE.MeshStandardMaterial({ color: 0x2f77bd, roughness: 0.5, metalness: 0.12 }),
    bBottom: new THREE.MeshStandardMaterial({ color: 0x7b5cc8, roughness: 0.55, metalness: 0.1 }),
    cTop: new THREE.MeshStandardMaterial({ color: 0x2f9c7c, roughness: 0.5, metalness: 0.12 }),
    cBottom: new THREE.MeshStandardMaterial({ color: 0xa2567c, roughness: 0.55, metalness: 0.1 }),
    edge: new THREE.LineBasicMaterial({ color: 0x121820, transparent: true, opacity: 0.52 }),
    hole: new THREE.MeshStandardMaterial({ color: 0x15191f, roughness: 0.72, metalness: 0.03 }),
    pin: new THREE.MeshStandardMaterial({ color: 0x22a36f, roughness: 0.32, metalness: 0.22 }),
    trace: new THREE.LineBasicMaterial({ color: 0x1d6fd6, transparent: true, opacity: 0.42 }),
    centerLine: new THREE.LineBasicMaterial({ color: 0x6b7480, transparent: true, opacity: 0.9 }),
    constraint: new THREE.LineBasicMaterial({ color: 0x22a36f, transparent: true, opacity: 0.9 }),
  };

  function rotate2(vector, angle) {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return new THREE.Vector2(c * vector.x - s * vector.y, s * vector.x + c * vector.y);
  }

  function attachmentSelectionKey() {
    return [aSiteInput.value, bSiteInput.value, bTopSiteInput.value, aLowerSiteInput.value].join("|");
  }

  function displayFrameForPose(pose) {
    const key = attachmentSelectionKey();
    const solvedAngle = Math.atan2(pose.bCenter.y, pose.bCenter.x);
    const safeSolvedAngle = Number.isFinite(solvedAngle) ? solvedAngle : 0;
    if (centerlineFrame.key !== key || centerlineFrame.angle === null) {
      centerlineFrame = { key, angle: safeSolvedAngle };
    }
    const displayRotation = centerlineFrame.angle - safeSolvedAngle;
    const displayedBVector = rotate2(pose.bCenter, displayRotation);
    const shift = displayedBVector.clone().multiplyScalar(-1);
    return {
      rotation: displayRotation,
      centerlineAngle: centerlineFrame.angle,
      mapPoint(point) {
        return rotate2(point, displayRotation).add(shift);
      },
    };
  }

  function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
  }

  function degToRad(value) {
    return (value * Math.PI) / 180;
  }

  function radToDeg(value) {
    return (value * 180) / Math.PI;
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
    const aAttachLocal = SITE_VECTORS[aSiteInput.value] || SITE_VECTORS.east;
    const bAttachLocal = SITE_VECTORS[bSiteInput.value] || SITE_VECTORS.south;
    const aLowerAttachLocal = SITE_VECTORS[aLowerSiteInput.value] || SITE_VECTORS.north;
    const bTopAttachLocal = SITE_VECTORS[bTopSiteInput.value] || SITE_VECTORS.south;
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
      aAttachLocal,
      bAttachLocal,
      aLowerAttachLocal,
      bTopAttachLocal,
      bcPrimaryBLocal,
      bcPrimaryCLocal,
      bcSecondBLocal,
      bcSecondCLocal,
      aAttachWorld,
      bAttachOffset,
      bCenter,
      cCenter,
      aLowerAttachWorld,
      bTopAttachOffset,
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
    const bodies = pose.bodies.map((body) => ({ ...body, primitives: bodyPrimitives(body) }));
    let maxPenetration = 0;
    let minClearance = Infinity;
    const collisions = [];
    for (let i = 0; i < bodies.length; i += 1) {
      for (let j = i + 1; j < bodies.length; j += 1) {
        const firstBody = bodies[i];
        const secondBody = bodies[j];
        if (firstBody.cell === secondBody.cell || !zBandsOverlap(firstBody, secondBody)) continue;
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

  function secondPinResidual(pose) {
    return pinResiduals(pose).max;
  }

  function exactSecondPinCandidates(aTopRad, requestedBLowerRad, requestedBTopRad) {
    const aAttachLocal = SITE_VECTORS[aSiteInput.value] || SITE_VECTORS.east;
    const bAttachLocal = SITE_VECTORS[bSiteInput.value] || SITE_VECTORS.south;
    const aLowerAttachLocal = SITE_VECTORS[aLowerSiteInput.value] || SITE_VECTORS.north;
    const bTopAttachLocal = SITE_VECTORS[bTopSiteInput.value] || SITE_VECTORS.south;
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

    function considerPose(pose) {
      const report = collisionReport(pose);
      const residual = secondPinResidual(pose);
      const angleScore = Math.abs(angleDelta(pose.bLowerRad, requestedBLowerRad)) + 0.65 * Math.abs(angleDelta(pose.bTopRad, requestedBTopRad));
      const score = angleScore + report.maxPenetration * 10;
      const candidate = { pose, report, residual, score, angleScore };
      const candidateFeasible = report.clear && residual <= FEASIBLE_PIN_RESIDUAL_MM;
      const bestFeasible = best && best.report.clear && best.residual <= FEASIBLE_PIN_RESIDUAL_MM;
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
    }

    exactSecondPinCandidates(aTopRad, requestedBLowerRad, requestedBTopRad).forEach(considerPose);

    if (!best) {
      const fallbackPose = buildPose(aTopRad, requestedBLowerRad, requestedBTopRad);
      const fallbackReport = collisionReport(fallbackPose);
      return {
        pose: fallbackPose,
        report: fallbackReport,
        residual: secondPinResidual(fallbackPose),
        score: Infinity,
        angleScore: Infinity,
        status: "blocked",
      };
    }

    const attached = best.residual <= FEASIBLE_PIN_RESIDUAL_MM;
    const clear = best.report.clear;
    return {
      ...best,
      status: clear && attached ? (best.angleScore <= 1e-9 ? "clear" : "adjusted") : "blocked",
    };
  }

  function solveNonPenetratingPose(aTopRad, requestedBLowerRad, requestedBTopRad) {
    return evaluateExactLoopClosure(aTopRad, requestedBLowerRad, requestedBTopRad);
  }

  function isDriveAngleFeasible(angleDeg, requestedBLowerRad, requestedBTopRad) {
    const solution = evaluateExactLoopClosure(degToRad(angleDeg), requestedBLowerRad, requestedBTopRad);
    return solution.status !== "blocked";
  }

  function computeDriveEnvelope(requestedBLowerRad, requestedBTopRad) {
    const feasible = [];
    for (let angleDeg = A_TOP_FULL_MIN_DEG; angleDeg <= A_TOP_FULL_MAX_DEG; angleDeg += DRIVE_SAMPLE_STEP_DEG) {
      if (isDriveAngleFeasible(angleDeg, requestedBLowerRad, requestedBTopRad)) feasible.push(angleDeg);
    }
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
    const key = [
      aSiteInput.value,
      bSiteInput.value,
      bTopSiteInput.value,
      aLowerSiteInput.value,
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

  function widestRange(envelope) {
    return envelope.ranges.reduce((best, range) => (range[1] - range[0] > best[1] - best[0] ? range : best), envelope.ranges[0]);
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

  function addEdges(parent, mesh) {
    const edges = new THREE.LineSegments(new THREE.EdgesGeometry(mesh.geometry, 24), materials.edge);
    edges.position.copy(mesh.position);
    edges.rotation.copy(mesh.rotation);
    edges.scale.copy(mesh.scale);
    parent.add(edges);
  }

  function addArm(parent, rotationZ, material) {
    const arm = new THREE.Mesh(
      new THREE.BoxGeometry(CAD.siteRadiusMm * 2, CAD.armWidthMm, CAD.bodyThicknessMm),
      material
    );
    arm.rotation.z = rotationZ;
    arm.castShadow = true;
    arm.receiveShadow = true;
    parent.add(arm);
    addEdges(parent, arm);
  }

  function addPad(parent, x, y, material) {
    const pad = cylinderZ(CAD.padRadiusMm, CAD.bodyThicknessMm + 0.05, material, 56);
    pad.position.set(x, y, 0);
    parent.add(pad);
    addEdges(parent, pad);

    const hole = cylinderZ(CAD.nominalHoleDiameterMm * 0.5, CAD.bodyThicknessMm + 0.24, materials.hole, 40);
    hole.position.set(x, y, 0.04);
    parent.add(hole);
  }

  function createCrossPart(name, material) {
    const group = new THREE.Group();
    group.name = name;
    addArm(group, 0, material);
    addArm(group, Math.PI / 2, material);
    const hub = cylinderZ(CAD.hubRadiusMm, CAD.bodyThicknessMm + 0.08, material, 64);
    group.add(hub);
    addEdges(group, hub);
    const centerHole = cylinderZ(CAD.nominalHoleDiameterMm * 0.5, CAD.bodyThicknessMm + 0.26, materials.hole, 48);
    centerHole.position.z = 0.05;
    group.add(centerHole);
    addPad(group, CAD.siteRadiusMm, 0, material);
    addPad(group, -CAD.siteRadiusMm, 0, material);
    addPad(group, 0, CAD.siteRadiusMm, material);
    addPad(group, 0, -CAD.siteRadiusMm, material);
    return group;
  }

  function createCell(name, topMaterial, bottomMaterial) {
    const group = new THREE.Group();
    group.name = name;
    const bottom = createCrossPart(`${name} lower cross`, bottomMaterial);
    const top = createCrossPart(`${name} upper cross`, topMaterial);
    bottom.position.z = -CAD.bodyThicknessMm * 0.5;
    top.position.z = CAD.bodyThicknessMm * 0.5;
    group.add(bottom);
    group.add(top);
    return { group, top, bottom };
  }

  const cellA = createCell("Cell A", materials.aTop, materials.aBottom);
  const cellB = createCell("Cell B", materials.bTop, materials.bBottom);
  const cellC = createCell("Cell C", materials.cTop, materials.cBottom);
  scene.add(cellA.group);
  scene.add(cellB.group);
  scene.add(cellC.group);

  const sharedPinAB = cylinderZ(CAD.nominalHoleDiameterMm * 0.42, CAD.bodyThicknessMm * 2.8, materials.pin, 48);
  const sharedPinSecond = cylinderZ(CAD.nominalHoleDiameterMm * 0.42, CAD.bodyThicknessMm * 2.8, materials.pin, 48);
  const sharedPinBCPrimary = cylinderZ(CAD.nominalHoleDiameterMm * 0.42, CAD.bodyThicknessMm * 2.8, materials.pin, 48);
  const sharedPinBCSecond = cylinderZ(CAD.nominalHoleDiameterMm * 0.42, CAD.bodyThicknessMm * 2.8, materials.pin, 48);
  scene.add(sharedPinAB);
  scene.add(sharedPinSecond);
  scene.add(sharedPinBCPrimary);
  scene.add(sharedPinBCSecond);

  const aMarker = cylinderZ(1.35, 1.0, materials.aTop, 24);
  const bMarker = cylinderZ(1.35, 1.0, materials.bBottom, 24);
  const bTopMarker = cylinderZ(1.35, 1.0, materials.bTop, 24);
  const aLowerMarker = cylinderZ(1.35, 1.0, materials.aBottom, 24);
  const bcPrimaryBMarker = cylinderZ(1.35, 1.0, materials.bBottom, 24);
  const bcPrimaryCMarker = cylinderZ(1.35, 1.0, materials.cTop, 24);
  const bcSecondBMarker = cylinderZ(1.35, 1.0, materials.bTop, 24);
  const bcSecondCMarker = cylinderZ(1.35, 1.0, materials.cBottom, 24);
  scene.add(aMarker);
  scene.add(bMarker);
  scene.add(bTopMarker);
  scene.add(aLowerMarker);
  scene.add(bcPrimaryBMarker);
  scene.add(bcPrimaryCMarker);
  scene.add(bcSecondBMarker);
  scene.add(bcSecondCMarker);

  const centerLine = new THREE.Line(new THREE.BufferGeometry(), materials.centerLine);
  const constraintLine = new THREE.Line(new THREE.BufferGeometry(), materials.constraint);
  const constraintLineSecond = new THREE.Line(new THREE.BufferGeometry(), materials.constraint);
  const constraintLineBCPrimary = new THREE.Line(new THREE.BufferGeometry(), materials.constraint);
  const constraintLineBCSecond = new THREE.Line(new THREE.BufferGeometry(), materials.constraint);
  const traceLine = new THREE.Line(new THREE.BufferGeometry(), materials.trace);
  scene.add(centerLine);
  scene.add(constraintLine);
  scene.add(constraintLineSecond);
  scene.add(constraintLineBCPrimary);
  scene.add(constraintLineBCSecond);
  scene.add(traceLine);

  function addAxisLine(points, color) {
    const geometry = new THREE.BufferGeometry().setFromPoints(points.map((p) => new THREE.Vector3(...p)));
    scene.add(new THREE.Line(geometry, new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.75 })));
  }

  addAxisLine(
    [
      [-52, 0, -8],
      [112, 0, -8],
    ],
    0xb24d3a
  );
  addAxisLine(
    [
      [0, -42, -8],
      [0, 122, -8],
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

  const grid = new THREE.GridHelper(180, 18, 0x8b97a5, 0xc4ccd4);
  grid.rotation.x = Math.PI / 2;
  grid.position.z = -8.1;
  scene.add(grid);

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

  function setLine(line, points) {
    line.geometry.dispose();
    line.geometry = new THREE.BufferGeometry().setFromPoints(points);
  }

  function updateCamera() {
    const r = cameraState.radius;
    const ce = Math.cos(cameraState.elevation);
    camera.position.set(
      cameraState.target.x + r * ce * Math.cos(cameraState.azimuth),
      cameraState.target.y + r * ce * Math.sin(cameraState.azimuth),
      cameraState.target.z + r * Math.sin(cameraState.elevation)
    );
    camera.lookAt(cameraState.target);
  }

  function setView(view) {
    if (view === "top") {
      cameraState.azimuth = -Math.PI / 2;
      cameraState.elevation = Math.PI / 2 - 0.001;
      cameraState.radius = 330;
    } else if (view === "front") {
      cameraState.azimuth = -Math.PI / 2;
      cameraState.elevation = 0.03;
      cameraState.radius = 360;
    } else if (view === "side") {
      cameraState.azimuth = 0;
      cameraState.elevation = 0.05;
      cameraState.radius = 360;
    } else {
      cameraState.azimuth = -0.82;
      cameraState.elevation = 0.58;
      cameraState.radius = 360;
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

  function updateMechanism(timeMs = 0) {
    const bLowerDeg = B_LOWER_REFERENCE_DEG;
    const bTopDeg = B_TOP_REFERENCE_DEG;
    const requestedBLower = degToRad(bLowerDeg);
    const requestedBTop = degToRad(bTopDeg);
    const driveEnvelope = cachedDriveEnvelope(requestedBLower, requestedBTop);
    allowedDriveMetric.textContent = formatDriveRanges(driveEnvelope);

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
      const mid = (range[0] + range[1]) * 0.5;
      const amplitude = Math.max(0, (range[1] - range[0]) * 0.5);
      requestedATopDeg = Math.round(mid + amplitude * Math.sin(timeMs * 0.001));
      aTopInput.value = requestedATopDeg.toFixed(0);
    }

    let aTopDeg = requestedATopDeg;
    let clamped = false;
    if (driveEnvelope && !isDriveAngleFeasible(aTopDeg, requestedBLower, requestedBTop)) {
      aTopDeg = nearestFeasibleDrive(aTopDeg, driveEnvelope);
      aTopInput.value = aTopDeg.toFixed(0);
      clamped = true;
    }

    driveClampMetric.textContent = clamped ? `${requestedATopDeg.toFixed(0)} -> ${aTopDeg.toFixed(0)} deg` : driveEnvelope ? "none" : "no feasible range";

    const requestedA = degToRad(aTopDeg);
    let solution = solveNonPenetratingPose(requestedA, requestedBLower, requestedBTop);
    if (clamped && solution.status === "clear") {
      solution = { ...solution, status: "adjusted" };
    }
    const pose = solution.pose;
    const displayFrame = displayFrameForPose(pose);
    const aCenterWorld = displayFrame.mapPoint(new THREE.Vector2(0, 0));
    const bCenterWorld = displayFrame.mapPoint(pose.bCenter);
    const cCenterWorld = displayFrame.mapPoint(pose.cCenter);
    const systemCenterWorld = aCenterWorld.clone().add(bCenterWorld).add(cCenterWorld).multiplyScalar(1 / 3);
    const aAttachWorld = displayFrame.mapPoint(pose.aAttachWorld);
    const bAttachOffset = pose.bAttachOffset;
    const bCenter = pose.bCenter;
    const aLowerAttachWorld = displayFrame.mapPoint(pose.aLowerAttachWorld);
    const bTopAttachWorld = displayFrame.mapPoint(pose.bTopAttachWorld);
    const bcPrimaryWorld = displayFrame.mapPoint(pose.bcPrimaryBWorld);
    const bcPrimaryCWorld = displayFrame.mapPoint(pose.bcPrimaryCWorld);
    const bcSecondWorld = displayFrame.mapPoint(pose.bcSecondBWorld);
    const bcSecondCWorld = displayFrame.mapPoint(pose.bcSecondCWorld);

    cellA.group.position.set(aCenterWorld.x, aCenterWorld.y, 0);
    cellA.bottom.rotation.z = displayFrame.rotation;
    cellA.top.rotation.z = pose.aTopRad + displayFrame.rotation;
    cellB.group.position.set(bCenterWorld.x, bCenterWorld.y, 0);
    cellB.bottom.rotation.z = pose.bLowerRad + displayFrame.rotation;
    cellB.top.rotation.z = pose.bTopWorldRad + displayFrame.rotation;
    cellC.group.position.set(cCenterWorld.x, cCenterWorld.y, 0);
    cellC.bottom.rotation.z = pose.cLowerRad + displayFrame.rotation;
    cellC.top.rotation.z = pose.cTopWorldRad + displayFrame.rotation;

    sharedPinAB.position.set(aAttachWorld.x, aAttachWorld.y, 0);
    sharedPinSecond.position.set(aLowerAttachWorld.x, aLowerAttachWorld.y, 0);
    sharedPinBCPrimary.position.set(bcPrimaryWorld.x, bcPrimaryWorld.y, 0);
    sharedPinBCSecond.position.set(bcSecondWorld.x, bcSecondWorld.y, 0);
    aMarker.position.set(aAttachWorld.x, aAttachWorld.y, CAD.bodyThicknessMm + 1.2);
    bMarker.position.set(aAttachWorld.x, aAttachWorld.y, -CAD.bodyThicknessMm - 1.2);
    bTopMarker.position.set(bTopAttachWorld.x, bTopAttachWorld.y, CAD.bodyThicknessMm + 1.2);
    aLowerMarker.position.set(aLowerAttachWorld.x, aLowerAttachWorld.y, -CAD.bodyThicknessMm - 1.2);
    bcPrimaryBMarker.position.set(bcPrimaryWorld.x, bcPrimaryWorld.y, -CAD.bodyThicknessMm - 1.2);
    bcPrimaryCMarker.position.set(bcPrimaryCWorld.x, bcPrimaryCWorld.y, CAD.bodyThicknessMm + 1.2);
    bcSecondBMarker.position.set(bcSecondWorld.x, bcSecondWorld.y, CAD.bodyThicknessMm + 1.2);
    bcSecondCMarker.position.set(bcSecondCWorld.x, bcSecondCWorld.y, -CAD.bodyThicknessMm - 1.2);

    setLine(centerLine, [
      new THREE.Vector3(aCenterWorld.x, aCenterWorld.y, 0.3),
      new THREE.Vector3(bCenterWorld.x, bCenterWorld.y, 0.3),
      new THREE.Vector3(cCenterWorld.x, cCenterWorld.y, 0.3),
    ]);
    setLine(constraintLine, [
      new THREE.Vector3(aAttachWorld.x, aAttachWorld.y, -CAD.bodyThicknessMm - 2),
      new THREE.Vector3(aAttachWorld.x, aAttachWorld.y, CAD.bodyThicknessMm + 2),
    ]);
    setLine(constraintLineSecond, [
      new THREE.Vector3(aLowerAttachWorld.x, aLowerAttachWorld.y, -CAD.bodyThicknessMm - 2),
      new THREE.Vector3(aLowerAttachWorld.x, aLowerAttachWorld.y, CAD.bodyThicknessMm + 2),
    ]);
    setLine(constraintLineBCPrimary, [
      new THREE.Vector3(bcPrimaryWorld.x, bcPrimaryWorld.y, -CAD.bodyThicknessMm - 2),
      new THREE.Vector3(bcPrimaryWorld.x, bcPrimaryWorld.y, CAD.bodyThicknessMm + 2),
    ]);
    setLine(constraintLineBCSecond, [
      new THREE.Vector3(bcSecondWorld.x, bcSecondWorld.y, -CAD.bodyThicknessMm - 2),
      new THREE.Vector3(bcSecondWorld.x, bcSecondWorld.y, CAD.bodyThicknessMm + 2),
    ]);

    const centerlineDirection = new THREE.Vector2(Math.cos(displayFrame.centerlineAngle), Math.sin(displayFrame.centerlineAngle));
    const halfPath = Math.max(CAD.siteRadiusMm * 3.2, bCenter.length() * 1.35);
    setLine(traceLine, [
      new THREE.Vector3(-centerlineDirection.x * halfPath, -centerlineDirection.y * halfPath, CAD.bodyThicknessMm + 2.6),
      new THREE.Vector3(centerlineDirection.x * halfPath, centerlineDirection.y * halfPath, CAD.bodyThicknessMm + 2.6),
    ]);

    const bAttachWorld = displayFrame.mapPoint(bCenter.clone().add(bAttachOffset));
    const residual = aAttachWorld.distanceTo(bAttachWorld);
    const pitch = bCenter.length();
    const residuals = pinResiduals(pose);
    const secondResidual = aLowerAttachWorld.distanceTo(bTopAttachWorld);
    const bcPitch = cCenterWorld.clone().sub(bCenterWorld).length();

    aTopOut.textContent = angleLabel(requestedATopDeg, pose.aTopRad);
    residualMetric.textContent = `${residual.toFixed(3)} mm`;
    pitchMetric.textContent = `${pitch.toFixed(1)} mm`;
    aCenterMetric.textContent = `(${aCenterWorld.x.toFixed(1)}, ${aCenterWorld.y.toFixed(1)})`;
    bCenterMetric.textContent = `(${bCenterWorld.x.toFixed(1)}, ${bCenterWorld.y.toFixed(1)})`;
    cCenterMetric.textContent = `(${cCenterWorld.x.toFixed(1)}, ${cCenterWorld.y.toFixed(1)})`;
    systemCenterMetric.textContent = `(${systemCenterWorld.x.toFixed(1)}, ${systemCenterWorld.y.toFixed(1)})`;
    secondResidualMetric.textContent = `${secondResidual.toFixed(3)} mm`;
    secondTargetMetric.textContent = `(${aLowerAttachWorld.x.toFixed(1)}, ${aLowerAttachWorld.y.toFixed(1)})`;
    secondBMetric.textContent = `(${bTopAttachWorld.x.toFixed(1)}, ${bTopAttachWorld.y.toFixed(1)})`;
    bcPrimaryResidualMetric.textContent = `${residuals.bcPrimary.toFixed(3)} mm`;
    bcSecondResidualMetric.textContent = `${residuals.bcSecond.toFixed(3)} mm`;
    bcPitchMetric.textContent = `${bcPitch.toFixed(1)} mm`;
    updateCollisionReadout(solution);
  }

  function resize() {
    const rect = mount.getBoundingClientRect();
    const width = Math.max(1, Math.floor(rect.width));
    const height = Math.max(1, Math.floor(rect.height));
    renderer.setSize(width, height, false);
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
  }

  let dragging = false;
  let lastX = 0;
  let lastY = 0;

  renderer.domElement.addEventListener("pointerdown", (event) => {
    dragging = true;
    lastX = event.clientX;
    lastY = event.clientY;
    renderer.domElement.setPointerCapture(event.pointerId);
  });

  renderer.domElement.addEventListener("pointermove", (event) => {
    if (!dragging) return;
    const dx = event.clientX - lastX;
    const dy = event.clientY - lastY;
    lastX = event.clientX;
    lastY = event.clientY;
    cameraState.azimuth -= dx * 0.006;
    cameraState.elevation = Math.max(-1.15, Math.min(1.45, cameraState.elevation + dy * 0.006));
    updateCamera();
  });

  renderer.domElement.addEventListener("pointerup", (event) => {
    dragging = false;
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
      cameraState.radius = Math.max(54, Math.min(460, cameraState.radius * (event.deltaY > 0 ? 1.08 : 0.92)));
      updateCamera();
    },
    { passive: false }
  );

  document.querySelectorAll("[data-view]").forEach((button) => {
    button.addEventListener("click", () => setView(button.dataset.view));
  });
  document.getElementById("resetView").addEventListener("click", () => setView("iso"));
  [aTopInput, aSiteInput, bSiteInput, bTopSiteInput, aLowerSiteInput, animateInput].forEach((input) => {
    input.addEventListener("input", () => updateMechanism());
    input.addEventListener("change", () => updateMechanism());
  });

  new ResizeObserver(resize).observe(mount);
  resize();
  setView("iso");

  function render(timeMs) {
    updateMechanism(timeMs);
    renderer.render(scene, camera);
    requestAnimationFrame(render);
  }

  requestAnimationFrame(render);
})();
