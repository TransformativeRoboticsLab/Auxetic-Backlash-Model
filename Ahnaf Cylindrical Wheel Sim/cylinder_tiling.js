(function () {
  "use strict";

  const THREE = window.THREE;
  const mount = document.getElementById("threeMount");

  const alphaInput = document.getElementById("alpha");
  const alphaOut = document.getElementById("alphaOut");
  const backlashInput = document.getElementById("backlash");
  const backlashOut = document.getElementById("backlashOut");
  const animateInput = document.getElementById("animate");
  const ringCountInput = document.getElementById("ringCount");
  const ringCountOut = document.getElementById("ringCountOut");
  const rowCountInput = document.getElementById("rowCount");
  const rowCountOut = document.getElementById("rowCountOut");
  const axialPitchInput = document.getElementById("axialPitch");
  const axialPitchOut = document.getElementById("axialPitchOut");
  const aSiteInput = document.getElementById("aSite");
  const bSiteInput = document.getElementById("bSite");
  const pitchMetric = document.getElementById("pitchMetric");
  const turnMetric = document.getElementById("turnMetric");
  const diameterMetric = document.getElementById("diameterMetric");
  const radiusMetric = document.getElementById("radiusMetric");
  const closureMetric = document.getElementById("closureMetric");
  const totalCellsMetric = document.getElementById("totalCellsMetric");
  const heightMetric = document.getElementById("heightMetric");
  const colorLegend = document.getElementById("colorLegend");
  const alphaCommandMetric = document.getElementById("alphaCommandMetric");
  const alphaEffectiveMetric = document.getElementById("alphaEffectiveMetric");
  const thetaMetric = document.getElementById("thetaMetric");
  const deadzoneState = document.getElementById("deadzoneState");
  const deadzoneRangeMetric = document.getElementById("deadzoneRangeMetric");
  const deadzoneBand = document.getElementById("deadzoneBand");
  const deadzoneMarker = document.getElementById("deadzoneMarker");
  const collisionEnabledInput = document.getElementById("collisionEnabled");
  const collisionState = document.getElementById("collisionState");
  const penetrationMetric = document.getElementById("penetrationMetric");
  const clearanceMetric = document.getElementById("clearanceMetric");
  const pairsCheckedMetric = document.getElementById("pairsCheckedMetric");
  const selectedCellStatus = document.getElementById("selectedCellStatus");
  const selectedIndexMetric = document.getElementById("selectedIndexMetric");
  const selectedCenterMetric = document.getElementById("selectedCenterMetric");
  const selectedBottomRotMetric = document.getElementById("selectedBottomRotMetric");
  const selectedTopRotMetric = document.getElementById("selectedTopRotMetric");
  const saveJsonBtn = document.getElementById("saveJsonBtn");
  const loadJsonBtn = document.getElementById("loadJsonBtn");
  const loadJsonFile = document.getElementById("loadJsonFile");
  const focusToggle = document.getElementById("focusToggle");
  const frameCellBtn = document.getElementById("frameCellBtn");
  const isolateCellBtn = document.getElementById("isolateCellBtn");
  const cellRoleSelect = document.getElementById("cellRoleSelect");
  const cellAlphaInput = document.getElementById("cellAlpha");
  const cellAlphaOut = document.getElementById("cellAlphaOut");
  const selectedCellAlphaMetric = document.getElementById("selectedCellAlphaMetric");
  const actuatorCountMetric = document.getElementById("actuatorCountMetric");
  const lockedCountMetric = document.getElementById("lockedCountMetric");
  const clearRolesBtn = document.getElementById("clearRolesBtn");
  const batchSelectedCountMetric = document.getElementById("batchSelectedCountMetric");
  const batchRoleSelect = document.getElementById("batchRoleSelect");
  const batchAlphaInput = document.getElementById("batchAlpha");
  const batchAlphaOut = document.getElementById("batchAlphaOut");
  const applyBatchBtn = document.getElementById("applyBatchBtn");
  const clearBatchSelectionBtn = document.getElementById("clearBatchSelectionBtn");
  const presetBarrelBtn = document.getElementById("presetBarrelBtn");
  const presetConeBtn = document.getElementById("presetConeBtn");
  const presetSaddleBtn = document.getElementById("presetSaddleBtn");
  const minDiameterMetric = document.getElementById("minDiameterMetric");
  const maxDiameterMetric = document.getElementById("maxDiameterMetric");
  const captureKeyframeBtn = document.getElementById("captureKeyframeBtn");
  const playSequenceBtn = document.getElementById("playSequenceBtn");
  const keyframeList = document.getElementById("keyframeList");
  const keyframeCountMetric = document.getElementById("keyframeCountMetric");
  const saveSequenceBtn = document.getElementById("saveSequenceBtn");
  const loadSequenceBtn = document.getElementById("loadSequenceBtn");
  const loadSequenceFile = document.getElementById("loadSequenceFile");
  const exportObjBtn = document.getElementById("exportObjBtn");
  const targetDiameterInput = document.getElementById("targetDiameter");
  const targetDiameterOut = document.getElementById("targetDiameterOut");
  const fitDiameterBtn = document.getElementById("fitDiameterBtn");
  const achievedDiameterMetric = document.getElementById("achievedDiameterMetric");
  const fitResidualMetric = document.getElementById("fitResidualMetric");
  const heatmapEnabledInput = document.getElementById("heatmapEnabled");
  const showMeasurementsInput = document.getElementById("showMeasurements");

  const CAD = Object.freeze({
    cellWidthMm: 55.604331,
    nominalHoleDiameterMm: 3.4,
    bodyThicknessMm: 4.0,
    padRadiusMm: 4.9,
    hubRadiusMm: 4.6,
    armWidthMm: 5.4,
    siteRadiusMm: 22.1,
  });

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

  const ALPHA_MIN = Number(alphaInput.min);
  const ALPHA_MAX = Number(alphaInput.max);
  const ALPHA_REFERENCE = 1.0;
  const THETA_SAFETY_MIN_DEG = -80;
  const THETA_SAFETY_MAX_DEG = 80;
  const RING_COUNT_MIN = Number(ringCountInput.min);
  const RING_COUNT_MAX = Number(ringCountInput.max);
  const ROW_COUNT_MIN = Number(rowCountInput.min);
  const ROW_COUNT_MAX = Number(rowCountInput.max);

  // Paper-supported RAD angle law (docs/research_grounding.md): theta = 70*alpha - 60.
  function alphaToThetaDeg(alphaValue) {
    return 70 * alphaValue - 60;
  }

  // Bidirectional dead-zone / ReLU coupling law from the same source:
  // f(x) = max(0, x - b) + min(x + b, 0). Applied here around alpha=1 (the
  // reference/undilated cell state) as this primitive's own single global
  // actuator's backlash gap, rather than the paper's original inter-cell
  // coupling use - see the in-page note for that caveat.
  function reluDeadzone(x, b) {
    return Math.max(0, x - b) + Math.min(x + b, 0);
  }

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(42, 1, 0.1, 1600);
  camera.up.set(0, 0, 1);

  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.shadowMap.enabled = true;
  // Pin the canvas's on-screen CSS size to its container explicitly. Without
  // this, renderer.setSize(w, h, false) (false = don't let three.js manage
  // the CSS size) leaves the canvas's displayed size to fall back to its
  // width/height attributes, which are the backing-buffer resolution
  // (w*devicePixelRatio) rather than CSS pixels whenever the display isn't
  // at 100% scaling. That made the canvas render larger than its container
  // on every resize pass, which (combined with resize() running every
  // frame, below) compounded into a runaway feedback loop.
  renderer.domElement.style.display = "block";
  renderer.domElement.style.width = "100%";
  renderer.domElement.style.height = "100%";
  mount.appendChild(renderer.domElement);

  const cameraState = {
    radius: 420,
    azimuth: -0.82,
    elevation: 0.58,
    target: new THREE.Vector3(0, 0, 0),
  };

  const sharedMaterials = {
    edge: new THREE.LineBasicMaterial({ color: 0x121820, transparent: true, opacity: 0.52 }),
    hole: new THREE.MeshStandardMaterial({ color: 0x15191f, roughness: 0.72, metalness: 0.03 }),
    actuatorMarker: new THREE.MeshStandardMaterial({ color: 0xe14b4b, emissive: 0x5a1414, roughness: 0.3 }),
    lockedMarker: new THREE.MeshStandardMaterial({ color: 0x4b7de1, emissive: 0x142a5a, roughness: 0.3 }),
    batchMarker: new THREE.MeshStandardMaterial({ color: 0x39e07a, emissive: 0x0f5a2c, roughness: 0.3 }),
  };

  // Row cross-color palette, cycled if there are more rows than palette entries.
  const ROW_PALETTE = [
    { top: 0xda7b27, bottom: 0x44515f },
    { top: 0x2f77bd, bottom: 0x7b5cc8 },
    { top: 0x2f9c7c, bottom: 0xa2567c },
    { top: 0xc2a83e, bottom: 0x5c7a99 },
    { top: 0xb2503a, bottom: 0x3c8f8f },
    { top: 0x7c8f3c, bottom: 0x8f3c7c },
  ];
  const rowMaterials = [];

  function materialsForRow(row) {
    if (!rowMaterials[row]) {
      const entry = ROW_PALETTE[row % ROW_PALETTE.length];
      rowMaterials[row] = {
        top: new THREE.MeshStandardMaterial({ color: entry.top, roughness: 0.5, metalness: 0.12 }),
        bottom: new THREE.MeshStandardMaterial({ color: entry.bottom, roughness: 0.55, metalness: 0.1 }),
      };
    }
    return rowMaterials[row];
  }

  function rotate2(vector, angle) {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return new THREE.Vector2(c * vector.x - s * vector.y, s * vector.x + c * vector.y);
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

  // --- Material clearance / collision (adapted from two_cell_attachment.js's
  // disk/capsule primitive checker, generalized to arbitrary neighbor pairs
  // instead of three hardcoded named cells). ---
  const BODY_EPSILON_MM = 1e-4;

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
      { type: "capsule", start: new THREE.Vector2(-r, 0), end: new THREE.Vector2(r, 0), radius: armRadius },
      { type: "capsule", start: new THREE.Vector2(0, -r), end: new THREE.Vector2(0, r), radius: armRadius },
      { type: "disk", center: new THREE.Vector2(0, 0), radius: CAD.hubRadiusMm },
      { type: "disk", center: SITE_VECTORS.east, radius: CAD.padRadiusMm },
      { type: "disk", center: SITE_VECTORS.north, radius: CAD.padRadiusMm },
      { type: "disk", center: SITE_VECTORS.west, radius: CAD.padRadiusMm },
      { type: "disk", center: SITE_VECTORS.south, radius: CAD.padRadiusMm },
    ];
  }

  const LOCAL_BODY_PRIMITIVES = localBodyPrimitives();

  function transformPoint(point, origin, angle) {
    return rotate2(point, angle).add(origin);
  }

  function bodyPrimitives(body) {
    return LOCAL_BODY_PRIMITIVES.map((primitive) => {
      if (primitive.type === "disk") {
        return {
          type: "disk",
          cell: body.cell,
          center: transformPoint(primitive.center, body.origin, body.rotation),
          radius: primitive.radius,
        };
      }
      return {
        type: "capsule",
        cell: body.cell,
        start: transformPoint(primitive.start, body.origin, body.rotation),
        end: transformPoint(primitive.end, body.origin, body.rotation),
        radius: primitive.radius,
      };
    });
  }

  // Checks material clearance only between physically adjacent cells
  // (circumferential ring neighbors and axial row neighbors), not every
  // n*m x n*m pair - non-adjacent cells sit far apart at typical ring
  // radii, and this keeps the check O(n*m) instead of O((n*m)^2).
  function checkNeighborClearance(n, m, rings, axialPitch, thetaPerRow) {
    function bodiesForCell(row, i) {
      const center = rings[row].centers[i];
      const z = row * axialPitch;
      const origin2 = new THREE.Vector2(center.x, center.y);
      return [
        { cell: `${row}-${i}-lower`, origin: origin2, rotation: rings[row].bottomRot[i], zMin: z - CAD.bodyThicknessMm, zMax: z },
        {
          cell: `${row}-${i}-upper`,
          origin: origin2,
          rotation: rings[row].bottomRot[i] + thetaPerRow[row][i],
          zMin: z,
          zMax: z + CAD.bodyThicknessMm,
        },
      ];
    }

    let maxPenetration = 0;
    let minClearance = Infinity;
    let pairsChecked = 0;

    function checkCellPair(rowA, iA, rowB, iB) {
      const bodiesA = bodiesForCell(rowA, iA).map((body) => ({ ...body, primitives: bodyPrimitives(body) }));
      const bodiesB = bodiesForCell(rowB, iB).map((body) => ({ ...body, primitives: bodyPrimitives(body) }));
      bodiesA.forEach((bodyA) => {
        bodiesB.forEach((bodyB) => {
          if (!zBandsOverlap(bodyA, bodyB)) return;
          bodyA.primitives.forEach((first) => {
            bodyB.primitives.forEach((second) => {
              const clearance = primitiveDistance(first, second);
              minClearance = Math.min(minClearance, clearance);
              if (clearance < 0) maxPenetration = Math.max(maxPenetration, -clearance);
            });
          });
        });
      });
      pairsChecked += 1;
    }

    for (let row = 0; row < m; row += 1) {
      for (let i = 0; i < n; i += 1) {
        checkCellPair(row, i, row, (i + 1) % n);
      }
    }
    for (let row = 0; row < m - 1; row += 1) {
      for (let i = 0; i < n; i += 1) {
        checkCellPair(row, i, row + 1, i);
      }
    }

    return {
      maxPenetration,
      minClearance: Number.isFinite(minClearance) ? minClearance : 0,
      pairsChecked,
      clear: maxPenetration <= BODY_EPSILON_MM,
    };
  }

  function cylinderZ(radius, depth, material, segments = 48) {
    const mesh = new THREE.Mesh(new THREE.CylinderGeometry(radius, radius, depth, segments), material);
    mesh.rotation.x = Math.PI / 2;
    mesh.castShadow = true;
    mesh.receiveShadow = true;
    return mesh;
  }

  function addEdges(parent, mesh) {
    const edges = new THREE.LineSegments(new THREE.EdgesGeometry(mesh.geometry, 24), sharedMaterials.edge);
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
    const pad = cylinderZ(CAD.padRadiusMm, CAD.bodyThicknessMm + 0.05, material, 40);
    pad.position.set(x, y, 0);
    parent.add(pad);
    addEdges(parent, pad);

    const hole = cylinderZ(CAD.nominalHoleDiameterMm * 0.5, CAD.bodyThicknessMm + 0.24, sharedMaterials.hole, 28);
    hole.position.set(x, y, 0.04);
    parent.add(hole);
  }

  function createCrossPart(name, material) {
    const group = new THREE.Group();
    group.name = name;
    addArm(group, 0, material);
    addArm(group, Math.PI / 2, material);
    const hub = cylinderZ(CAD.hubRadiusMm, CAD.bodyThicknessMm + 0.08, material, 48);
    group.add(hub);
    addEdges(group, hub);
    const centerHole = cylinderZ(CAD.nominalHoleDiameterMm * 0.5, CAD.bodyThicknessMm + 0.26, sharedMaterials.hole, 32);
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

  function disposeGroup(group) {
    group.traverse((object) => {
      if (object.geometry) object.geometry.dispose();
    });
  }

  let cellPool = [];
  let poolN = -1;
  let poolM = -1;

  // Each cell gets its own material instances (not shared per row) so the
  // alpha heatmap toggle can recolor individual cells without disturbing
  // their row-mates.
  function createCellMaterials(row) {
    const entry = ROW_PALETTE[row % ROW_PALETTE.length];
    return {
      top: new THREE.MeshStandardMaterial({ color: entry.top, roughness: 0.5, metalness: 0.12 }),
      bottom: new THREE.MeshStandardMaterial({ color: entry.bottom, roughness: 0.55, metalness: 0.1 }),
    };
  }

  function rebuildPool(n, m) {
    cellPool.forEach((cell) => {
      scene.remove(cell.group);
      disposeGroup(cell.group);
      if (cell.topMaterial) cell.topMaterial.dispose();
      if (cell.bottomMaterial) cell.bottomMaterial.dispose();
      if (cell.roleMarker) scene.remove(cell.roleMarker);
      if (cell.batchMarker) scene.remove(cell.batchMarker);
    });
    cellPool = [];
    multiSelected.clear();
    for (let row = 0; row < m; row += 1) {
      for (let i = 0; i < n; i += 1) {
        const mats = createCellMaterials(row);
        const cell = createCell(`ring${row}-cell${i}`, mats.top, mats.bottom);
        cell.topMaterial = mats.top;
        cell.bottomMaterial = mats.bottom;
        cell.baseRow = row;
        scene.add(cell.group);
        const roleMarker = cylinderZ(2.6, CAD.bodyThicknessMm * 2.6, sharedMaterials.actuatorMarker, 20);
        roleMarker.visible = false;
        scene.add(roleMarker);
        cell.roleMarker = roleMarker;
        const batchMarker = cylinderZ(4.4, CAD.bodyThicknessMm * 3.6, sharedMaterials.batchMarker, 24);
        batchMarker.visible = false;
        scene.add(batchMarker);
        cell.batchMarker = batchMarker;
        cellPool.push(cell);
      }
    }
    poolN = n;
    poolM = m;
  }

  function updateLegend(m) {
    colorLegend.innerHTML = "";
    for (let row = 0; row < m; row += 1) {
      const mats = materialsForRow(row);
      const topColor = `#${mats.top.color.getHexString()}`;
      const bottomColor = `#${mats.bottom.color.getHexString()}`;
      const upper = document.createElement("span");
      upper.innerHTML = `<i class="swatch" style="background:${topColor}"></i>Row ${row + 1} upper cross`;
      const lower = document.createElement("span");
      lower.innerHTML = `<i class="swatch" style="background:${bottomColor}"></i>Row ${row + 1} lower cross`;
      colorLegend.appendChild(upper);
      colorLegend.appendChild(lower);
    }
  }

  // Close n single-pin joints into a ring: cell i's upper-cross hole (aSite,
  // at absolute rotation bottomRot[i] + thetaPerCell[i]) meets cell i+1's
  // lower-cross hole (bSite, at absolute rotation bottomRot[i+1]), where each
  // step turns by a fixed 2*pi/n regardless of thetaPerCell. This is the same
  // forward relation two_cell_attachment.js uses for its A-B pin
  // (aAttachWorld - bAttachOffset), just walked n times with a fixed turn
  // increment instead of laid flat.
  //
  // When every thetaPerCell[i] is equal this is an exact regular n-gon (the
  // classic turtle-graphics closure: identical local steps + a turn summing
  // to 2*pi always return to the start) - closureResidual is ~0. When cells
  // differ, the fixed-turn construction still walks all the way around
  // (heading always returns to 0 after n steps), but the *chord lengths*
  // differ per joint, so the loop generally will not land back exactly on
  // its own start - closureResidual becomes a genuine, meaningful measure of
  // how inconsistent the per-cell commands are with a physically closed
  // ring, the same way real backlash/compliance would have to absorb that
  // mismatch in hardware.
  function buildRing(n, thetaPerCell) {
    const aAttachLocal = SITE_VECTORS[aSiteInput.value] || SITE_VECTORS.east;
    const bAttachLocal = SITE_VECTORS[bSiteInput.value] || SITE_VECTORS.west;
    const turn = (2 * Math.PI) / n;
    const centers = [new THREE.Vector2(0, 0)];
    const bottomRot = [0];
    let heading = 0;
    let pitchSum = 0;
    let closureResidual = 0;
    for (let i = 0; i < n; i += 1) {
      const aAttachWorld = rotate2(aAttachLocal, heading + thetaPerCell[i]);
      heading += turn;
      const bAttachWorld = rotate2(bAttachLocal, heading);
      const step = aAttachWorld.clone().sub(bAttachWorld);
      const nextCenter = centers[i].clone().add(step);
      pitchSum += step.length();
      if (i < n - 1) {
        centers.push(nextCenter);
        bottomRot.push(heading);
      } else {
        closureResidual = nextCenter.distanceTo(centers[0]);
      }
    }
    const centroid = centers.reduce((acc, c) => acc.add(c), new THREE.Vector2(0, 0)).multiplyScalar(1 / n);
    let radius = 0;
    centers.forEach((c) => {
      radius += c.distanceTo(centroid);
    });
    radius /= n;
    // bottomRot[i] is the cell's own chain-walking heading, not its polar
    // position angle around the ring - for a turtle-graphics polygon walk
    // those differ by a constant phase (heading leads the true radial
    // direction, exact amount depends on n and the chosen attachment sites).
    // Radial cell orientation needs the true outward direction from the
    // ring's own centroid, so compute it directly from each cell's actual
    // position rather than reusing bottomRot.
    const radialAngle = centers.map((c) => Math.atan2(c.y - centroid.y, c.x - centroid.x));
    return { centers, bottomRot, radialAngle, pitch: pitchSum / n, closureResidual, radius, diameter: 2 * radius, centroid };
  }

  // --- Per-cell actuation: each cell is "free" (tracks a backlash-gated
  // relaxation toward its neighbors' alpha), "actuator" (a directly
  // commanded alpha), or "locked" (frozen at whatever alpha it had when
  // locked). This is the discrete analogue of the paper's coupling law
  // (docs/research_grounding.md) applied over the ring/cylinder's actual
  // neighbor graph (circumferential + axial) instead of only around a
  // single global reference - actuating one cell should tug its neighbors
  // by die-off distance, gated by the same backlash gap, rather than every
  // cell being an independent free input.
  let cellRoles = [];
  let rolesN = -1;
  let rolesM = -1;

  function ensureRoles(n, m) {
    if (n === rolesN && m === rolesM) return;
    cellRoles = [];
    for (let row = 0; row < m; row += 1) {
      const r = [];
      for (let i = 0; i < n; i += 1) r.push({ role: "free", alpha: ALPHA_REFERENCE });
      cellRoles.push(r);
    }
    rolesN = n;
    rolesM = m;
  }

  function isFixedCell(row, i) {
    return cellRoles[row][i].role !== "free";
  }

  const PROPAGATION_ITERATIONS = 24;

  function computeCellAlphas(n, m, baselineAlpha, backlash) {
    let alphas = [];
    for (let row = 0; row < m; row += 1) alphas.push(new Array(n).fill(baselineAlpha));
    for (let row = 0; row < m; row += 1) {
      for (let i = 0; i < n; i += 1) {
        if (isFixedCell(row, i)) alphas[row][i] = cellRoles[row][i].alpha;
      }
    }
    for (let iter = 0; iter < PROPAGATION_ITERATIONS; iter += 1) {
      const next = alphas.map((r) => r.slice());
      for (let row = 0; row < m; row += 1) {
        for (let i = 0; i < n; i += 1) {
          if (isFixedCell(row, i)) continue;
          const neighbors = [alphas[row][(i - 1 + n) % n], alphas[row][(i + 1) % n]];
          if (row > 0) neighbors.push(alphas[row - 1][i]);
          if (row < m - 1) neighbors.push(alphas[row + 1][i]);
          const avg = neighbors.reduce((a, b) => a + b, 0) / neighbors.length;
          next[row][i] = alphas[row][i] + reluDeadzone(avg - alphas[row][i], backlash);
        }
      }
      alphas = next;
    }
    return alphas;
  }

  function cellThetaDeg(alphaValue) {
    return clamp(alphaToThetaDeg(alphaValue), THETA_SAFETY_MIN_DEG, THETA_SAFETY_MAX_DEG);
  }

  // Draws the backlash dead-zone as a shaded band under the Drive slider
  // (a separate ruler element, not a background painted onto the range
  // input itself - modern Chrome's native "filled track" range-input
  // styling draws over any CSS background set directly on the element, so
  // that approach was tried first and silently did nothing visible), plus
  // a marker for the current commanded alpha. Makes it visually obvious -
  // not just readable in a separate metric - which alpha values are "free"
  // (inside the gap, dragging does nothing) versus "engaged" (outside it,
  // dragging moves the wheel). This is the direct fix for alpha and
  // backlash otherwise reading as "doing the same thing": from an
  // already-engaged alpha, both controls do shift the same effective value
  // (by design - that's the real coupling law), but nothing on screen
  // explained why until now.
  function updateAlphaSliderTrack(backlash, alphaCommand) {
    const range = ALPHA_MAX - ALPHA_MIN;
    const lowPct = clamp(((ALPHA_REFERENCE - backlash - ALPHA_MIN) / range) * 100, 0, 100);
    const highPct = clamp(((ALPHA_REFERENCE + backlash - ALPHA_MIN) / range) * 100, 0, 100);
    deadzoneBand.style.left = `${lowPct}%`;
    deadzoneBand.style.width = `${Math.max(0, highPct - lowPct)}%`;
    const markerPct = clamp(((alphaCommand - ALPHA_MIN) / range) * 100, 0, 100);
    deadzoneMarker.style.left = `${markerPct}%`;
  }

  // Alpha heatmap: blue (contracted, near ALPHA_MIN) -> light gray (neutral,
  // mid-range) -> red (expanded, near ALPHA_MAX), so coupling propagation
  // and per-row diameter differences are visible at a glance instead of
  // only readable by clicking through cells one at a time.
  const HEATMAP_LOW = new THREE.Color(0x2f6fd6);
  const HEATMAP_MID = new THREE.Color(0xd8d8d8);
  const HEATMAP_HIGH = new THREE.Color(0xd63a2f);
  const heatmapColorScratch = new THREE.Color();

  // Per Jacob's feedback: cells were rendering face-on to the axle (their
  // flat cross plane normal to world Z), stacking rings like coins on a
  // rod. The reference structure instead has each cell's face tangent to
  // the cylinder's own surface - normal pointing radially outward, like
  // scales on the tube - with its arms sweeping in the tangential/axial
  // plane. This only changes each cell GROUP's base orientation; the
  // existing bottom/top .rotation.z calls are left untouched and keep
  // working exactly as before, now just measured within this reoriented
  // local frame instead of world XY - matching "equations unchanged, only
  // the visualization needs adjusting". Local X (where the east/west
  // SITE_VECTORS point) maps to the tangential direction, so the existing
  // circumferential attach math is unaffected; local Y maps to axial.
  const radialScratchX = new THREE.Vector3();
  const radialScratchY = new THREE.Vector3(0, 0, 1);
  const radialScratchZ = new THREE.Vector3();
  const radialScratchMatrix = new THREE.Matrix4();

  function setRadialOrientation(object3d, theta) {
    const cosT = Math.cos(theta);
    const sinT = Math.sin(theta);
    radialScratchX.set(-sinT, cosT, 0); // tangential
    radialScratchZ.set(cosT, sinT, 0); // radial (outward)
    radialScratchMatrix.makeBasis(radialScratchX, radialScratchY, radialScratchZ);
    object3d.quaternion.setFromRotationMatrix(radialScratchMatrix);
  }

  function heatColor(alphaValue) {
    const t = clamp((alphaValue - ALPHA_MIN) / (ALPHA_MAX - ALPHA_MIN), 0, 1);
    if (t < 0.5) heatmapColorScratch.lerpColors(HEATMAP_LOW, HEATMAP_MID, t / 0.5);
    else heatmapColorScratch.lerpColors(HEATMAP_MID, HEATMAP_HIGH, (t - 0.5) / 0.5);
    return heatmapColorScratch;
  }

  // Diameter vs alpha is NOT monotonic (it rises then falls as cells
  // over-rotate past their most-open pose - verified numerically before
  // building this: for n=10, diameter peaks around alpha=1.4 then drops
  // back down toward alpha=2.0), so this can't be solved with a bisection
  // search or a closed-form inverse of theta=70*alpha-60. Instead it
  // exhaustively evaluates every alpha the slider can reach (a uniform
  // ring, all cells at the same theta) and keeps whichever is closest to
  // the target - simple, robust to non-monotonicity, and cheap since a
  // whole ring build is just O(n).
  function fitAlphaToDiameter(targetDiameter, n, backlash) {
    let bestAlpha = ALPHA_REFERENCE;
    let bestDiameter = 0;
    let bestDiff = Infinity;
    const step = 0.01;
    const steps = Math.round((ALPHA_MAX - ALPHA_MIN) / step);
    for (let s = 0; s <= steps; s += 1) {
      const a = Math.round((ALPHA_MIN + s * step) * 100) / 100;
      // Must match updateMechanism's pipeline exactly (commanded alpha ->
      // backlash dead-zone -> theta), otherwise the alpha this picks would
      // get shifted by the dead-zone once actually applied and render at a
      // different diameter than what was just fit.
      const effectiveAlpha = ALPHA_REFERENCE + reluDeadzone(a - ALPHA_REFERENCE, backlash);
      const thetaRad = degToRad(cellThetaDeg(effectiveAlpha));
      const diameter = buildRing(n, new Array(n).fill(thetaRad)).diameter;
      const diff = Math.abs(diameter - targetDiameter);
      if (diff < bestDiff) {
        bestDiff = diff;
        bestAlpha = a;
        bestDiameter = diameter;
      }
    }
    return { alpha: bestAlpha, achievedDiameter: bestDiameter };
  }

  function addAxisLine(points, color) {
    const geometry = new THREE.BufferGeometry().setFromPoints(points.map((p) => new THREE.Vector3(...p)));
    scene.add(new THREE.Line(geometry, new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.75 })));
  }

  addAxisLine(
    [
      [-140, 0, -14],
      [140, 0, -14],
    ],
    0xb24d3a
  );
  addAxisLine(
    [
      [0, -140, -14],
      [0, 140, -14],
    ],
    0x1d6fd6
  );
  addAxisLine(
    [
      [0, 0, -14],
      [0, 0, 280],
    ],
    0x22a36f
  );

  const grid = new THREE.GridHelper(320, 32, 0x8b97a5, 0xc4ccd4);
  grid.rotation.x = Math.PI / 2;
  grid.position.z = -14.1;
  scene.add(grid);

  scene.add(new THREE.HemisphereLight(0xffffff, 0x5f6974, 1.65));
  const key = new THREE.DirectionalLight(0xffffff, 2.15);
  key.position.set(60, -90, 160);
  key.castShadow = true;
  key.shadow.mapSize.width = 1024;
  key.shadow.mapSize.height = 1024;
  scene.add(key);
  const fill = new THREE.DirectionalLight(0xb8d4ff, 0.65);
  fill.position.set(-90, 100, 90);
  scene.add(fill);

  const selectionMarker = cylinderZ(3.2, CAD.bodyThicknessMm * 3.2, new THREE.MeshStandardMaterial({ color: 0xffd54a, emissive: 0x7a5c00, roughness: 0.3 }), 32);
  selectionMarker.visible = false;
  scene.add(selectionMarker);

  // CAD-style dimension line across row 0's diameter: a straight line
  // through the ring's own centroid plus small perpendicular end ticks,
  // for a visual scale reference alongside the numeric "wheel diameter"
  // readout - most useful together with Export OBJ, when eyeballing a
  // configuration against real-world size before pulling it into CAD.
  const measurementMaterial = new THREE.LineBasicMaterial({ color: 0x121820 });
  const measurementLine = new THREE.Line(new THREE.BufferGeometry(), measurementMaterial);
  const measurementTickA = new THREE.Line(new THREE.BufferGeometry(), measurementMaterial);
  const measurementTickB = new THREE.Line(new THREE.BufferGeometry(), measurementMaterial);
  measurementLine.visible = false;
  measurementTickA.visible = false;
  measurementTickB.visible = false;
  scene.add(measurementLine, measurementTickA, measurementTickB);

  function setLineGeometry(line, points) {
    line.geometry.dispose();
    line.geometry = new THREE.BufferGeometry().setFromPoints(points);
  }

  // A dimension line drawn straight through the ring's own center (at the
  // cells' own z) gets occluded by the cell geometry itself from most
  // camera angles - confirmed by screenshot before landing on this
  // approach. Real CAD/engineering drawings offset the dimension line
  // clear of the part instead, connected back to the actual measured
  // points by short extension lines, which is what this draws: extension
  // lines from the true diameter endpoints (at row 0's own z) down to a
  // dimension line comfortably below the whole assembly.
  const MEASUREMENT_DROP_MM = 30;

  function updateMeasurementLine(ring0, z) {
    const c = ring0.centroid;
    const r = ring0.radius;
    const dimZ = z - MEASUREMENT_DROP_MM;
    const pointA = new THREE.Vector3(c.x - r, c.y, z);
    const pointB = new THREE.Vector3(c.x + r, c.y, z);
    const dimA = new THREE.Vector3(c.x - r, c.y, dimZ);
    const dimB = new THREE.Vector3(c.x + r, c.y, dimZ);
    setLineGeometry(measurementLine, [dimA, dimB]);
    setLineGeometry(measurementTickA, [pointA, dimA]);
    setLineGeometry(measurementTickB, [pointB, dimB]);
  }

  const raycaster = new THREE.Raycaster();
  const pointerNdc = new THREE.Vector2();
  let selected = null; // { row, i }
  const multiSelected = new Set(); // "row,i" keys - batch selection for differential dilation
  function cellKey(row, i) {
    return `${row},${i}`;
  }
  let isolateActive = false;
  const selectedWorld = new THREE.Vector3();
  let lastCellAlphas = null;
  let lastCellRolesSnapshot = null;
  let lastRingsSnapshot = null;

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
      cameraState.radius = 380;
    } else if (view === "front") {
      cameraState.azimuth = -Math.PI / 2;
      cameraState.elevation = 0.03;
      cameraState.radius = 420;
    } else if (view === "side") {
      cameraState.azimuth = 0;
      cameraState.elevation = 0.05;
      cameraState.radius = 420;
    } else {
      cameraState.azimuth = -0.82;
      cameraState.elevation = 0.58;
      cameraState.radius = 420;
    }
    cameraState.target.set(0, 0, 60);
    updateCamera();
    document.querySelectorAll("[data-view]").forEach((button) => {
      button.setAttribute("aria-pressed", String(button.dataset.view === view));
    });
  }

  function updateMechanism(timeMs = 0) {
    const n = clamp(Math.round(Number(ringCountInput.value)), RING_COUNT_MIN, RING_COUNT_MAX);
    const m = clamp(Math.round(Number(rowCountInput.value)), ROW_COUNT_MIN, ROW_COUNT_MAX);
    const axialPitch = Number(axialPitchInput.value);

    let alphaCommand = Number(alphaInput.value);
    if (animateInput.checked) {
      const mid = (ALPHA_MIN + ALPHA_MAX) * 0.5;
      const amplitude = (ALPHA_MAX - ALPHA_MIN) * 0.5;
      alphaCommand = Math.round((mid + amplitude * Math.sin(timeMs * 0.0005)) * 100) / 100;
      alphaInput.value = String(alphaCommand);
    }
    const backlash = Number(backlashInput.value);
    const baselineAlpha = ALPHA_REFERENCE + reluDeadzone(alphaCommand - ALPHA_REFERENCE, backlash);
    const baselineThetaDeg = cellThetaDeg(baselineAlpha);
    const inDeadzone = Math.abs(alphaCommand - ALPHA_REFERENCE) <= backlash;

    if (n !== poolN || m !== poolM) {
      rebuildPool(n, m);
      updateLegend(m);
      if (selected && (selected.row >= m || selected.i >= n)) {
        selected = null;
      }
    }
    ensureRoles(n, m);

    const cellAlphas = computeCellAlphas(n, m, baselineAlpha, backlash);
    lastCellAlphas = cellAlphas;
    lastCellRolesSnapshot = cellRoles;
    const thetaPerRow = cellAlphas.map((row) => row.map((a) => degToRad(cellThetaDeg(a))));
    const rings = [];
    for (let row = 0; row < m; row += 1) rings.push(buildRing(n, thetaPerRow[row]));
    lastRingsSnapshot = rings;

    const showMeasurements = showMeasurementsInput.checked && !isolateActive;
    measurementLine.visible = showMeasurements;
    measurementTickA.visible = showMeasurements;
    measurementTickB.visible = showMeasurements;
    if (showMeasurements) updateMeasurementLine(rings[0], 0);

    for (let row = 0; row < m; row += 1) {
      const z = row * axialPitch;
      for (let i = 0; i < n; i += 1) {
        const cell = cellPool[row * n + i];
        const center = rings[row].centers[i];
        cell.group.position.set(center.x, center.y, z);
        setRadialOrientation(cell.group, rings[row].radialAngle[i]);
        // The group's own orientation now carries the cell's true outward
        // direction (radialAngle[i], measured from the ring's own centroid -
        // not bottomRot[i], which is the chain-walking heading used for pin
        // closure and leads the true radial direction by a construction-
        // dependent phase, not a fixed 90deg, so reusing it here pointed the
        // wide cross face radially instead of tangentially). The child
        // crosses' local z-rotation is only the *relative* dilation twist
        // between layers on top of that heading, not heading + twist -
        // otherwise heading would be double-applied. The bottom cross's
        // 4-fold symmetry makes 0 an arbitrary but equally valid reference;
        // top stays exactly theta ahead of bottom, preserving the same
        // "opens up by theta" dilation visual as before.
        cell.bottom.rotation.z = 0;
        cell.top.rotation.z = thetaPerRow[row][i];
        const role = cellRoles[row][i].role;
        if (role === "free") {
          cell.roleMarker.visible = false;
        } else {
          cell.roleMarker.visible = true;
          cell.roleMarker.material = role === "actuator" ? sharedMaterials.actuatorMarker : sharedMaterials.lockedMarker;
          cell.roleMarker.position.set(center.x, center.y, z);
        }
        if (multiSelected.has(cellKey(row, i))) {
          cell.batchMarker.visible = true;
          cell.batchMarker.position.set(center.x, center.y, z - 0.01);
        } else {
          cell.batchMarker.visible = false;
        }
        if (heatmapEnabledInput.checked) {
          const color = heatColor(cellAlphas[row][i]);
          cell.topMaterial.color.copy(color);
          cell.bottomMaterial.color.copy(color);
        } else {
          const entry = ROW_PALETTE[row % ROW_PALETTE.length];
          cell.topMaterial.color.setHex(entry.top);
          cell.bottomMaterial.color.setHex(entry.bottom);
        }
      }
    }

    alphaOut.textContent = alphaCommand.toFixed(2);
    backlashOut.textContent = backlash.toFixed(2);
    alphaCommandMetric.textContent = alphaCommand.toFixed(2);
    alphaEffectiveMetric.textContent = baselineAlpha.toFixed(2);
    thetaMetric.textContent = `${baselineThetaDeg.toFixed(1)} deg`;
    deadzoneState.textContent = inDeadzone ? "free (dead zone)" : "engaged";
    deadzoneState.classList.toggle("status-adjusted", inDeadzone);
    deadzoneState.classList.toggle("status-ok", !inDeadzone);
    deadzoneRangeMetric.textContent = `${(ALPHA_REFERENCE - backlash).toFixed(2)} - ${(ALPHA_REFERENCE + backlash).toFixed(2)}`;
    updateAlphaSliderTrack(backlash, alphaCommand);

    ringCountOut.textContent = String(n);
    rowCountOut.textContent = String(m);
    axialPitchOut.textContent = `${axialPitch.toFixed(0)} mm`;
    pitchMetric.textContent = `${rings[0].pitch.toFixed(1)} mm`;
    turnMetric.textContent = `${radToDeg((2 * Math.PI) / n).toFixed(1)} deg`;
    diameterMetric.textContent = `${rings[0].diameter.toFixed(1)} mm`;
    radiusMetric.textContent = `${rings[0].radius.toFixed(1)} mm`;
    closureMetric.textContent = `${rings[0].closureResidual.toFixed(4)} mm`;
    totalCellsMetric.textContent = String(n * m);
    heightMetric.textContent = `${((m - 1) * axialPitch).toFixed(1)} mm`;

    let minDiameter = Infinity;
    let maxDiameter = -Infinity;
    rings.forEach((r) => {
      minDiameter = Math.min(minDiameter, r.diameter);
      maxDiameter = Math.max(maxDiameter, r.diameter);
    });
    minDiameterMetric.textContent = minDiameter.toFixed(1);
    maxDiameterMetric.textContent = maxDiameter.toFixed(1);

    let actuatorCount = 0;
    let lockedCount = 0;
    for (let row = 0; row < m; row += 1) {
      for (let i = 0; i < n; i += 1) {
        if (cellRoles[row][i].role === "actuator") actuatorCount += 1;
        else if (cellRoles[row][i].role === "locked") lockedCount += 1;
      }
    }
    actuatorCountMetric.textContent = String(actuatorCount);
    lockedCountMetric.textContent = String(lockedCount);

    if (collisionEnabledInput.checked) {
      const report = checkNeighborClearance(n, m, rings, axialPitch, thetaPerRow);
      collisionState.textContent = report.clear ? "clear" : "blocked";
      collisionState.classList.toggle("status-ok", report.clear);
      collisionState.classList.toggle("status-blocked", !report.clear);
      penetrationMetric.textContent = `${report.maxPenetration.toFixed(3)} mm`;
      clearanceMetric.textContent = `${Math.max(0, report.minClearance).toFixed(3)} mm`;
      pairsCheckedMetric.textContent = String(report.pairsChecked);
    } else {
      collisionState.textContent = "not checked";
      collisionState.classList.remove("status-ok", "status-blocked");
      penetrationMetric.textContent = "-";
      clearanceMetric.textContent = "-";
      pairsCheckedMetric.textContent = "0";
    }

    if (selected && selected.row < m && selected.i < n) {
      const row = selected.row;
      const i = selected.i;
      const center = rings[row].centers[i];
      const z = row * axialPitch;
      selectionMarker.visible = true;
      selectionMarker.position.set(center.x, center.y, z + 0.01);
      selectedWorld.set(center.x, center.y, z);
      selectedCellStatus.textContent = `row ${row}, cell ${i}`;
      selectedCellStatus.classList.add("status-ok");
      selectedIndexMetric.textContent = `${row}, ${i}`;
      selectedCenterMetric.textContent = `(${center.x.toFixed(1)}, ${center.y.toFixed(1)}, ${z.toFixed(1)})`;
      selectedBottomRotMetric.textContent = `${radToDeg(rings[row].bottomRot[i]).toFixed(1)} deg`;
      selectedTopRotMetric.textContent = `${radToDeg(rings[row].bottomRot[i] + thetaPerRow[row][i]).toFixed(1)} deg`;
      selectedCellAlphaMetric.textContent = cellAlphas[row][i].toFixed(2);
      frameCellBtn.disabled = false;
      isolateCellBtn.disabled = false;
      cellRoleSelect.disabled = false;
      const role = cellRoles[row][i].role;
      if (document.activeElement !== cellRoleSelect) cellRoleSelect.value = role;
      cellAlphaInput.disabled = role === "free";
      if (document.activeElement !== cellAlphaInput) {
        cellAlphaInput.value = String(cellAlphas[row][i]);
      }
      cellAlphaOut.textContent = cellAlphas[row][i].toFixed(2);
    } else {
      selectionMarker.visible = false;
      selectedCellStatus.textContent = "no cell selected";
      selectedCellStatus.classList.remove("status-ok");
      selectedIndexMetric.textContent = "-";
      selectedCenterMetric.textContent = "-";
      selectedBottomRotMetric.textContent = "-";
      selectedTopRotMetric.textContent = "-";
      selectedCellAlphaMetric.textContent = "-";
      frameCellBtn.disabled = true;
      isolateCellBtn.disabled = true;
      cellRoleSelect.disabled = true;
      cellAlphaInput.disabled = true;
      cellAlphaOut.textContent = "-";
      if (isolateActive) {
        isolateActive = false;
        isolateCellBtn.textContent = "Isolate Cell";
      }
    }

    for (let row = 0; row < m; row += 1) {
      for (let i = 0; i < n; i += 1) {
        const cell = cellPool[row * n + i];
        const visible = !isolateActive || (selected && row === selected.row && i === selected.i);
        cell.group.visible = visible;
        cell.roleMarker.visible = visible && cellRoles[row][i].role !== "free";
        cell.batchMarker.visible = visible && multiSelected.has(cellKey(row, i));
      }
    }
    grid.visible = !isolateActive;
  }

  let lastWidth = -1;
  let lastHeight = -1;

  function resize() {
    const rect = mount.getBoundingClientRect();
    const width = Math.max(1, Math.floor(rect.width));
    const height = Math.max(1, Math.floor(rect.height));
    // Skip degenerate/unsettled layout sizes (e.g. a page that reports 0x0
    // for a frame or two while its CSS grid is still resolving) rather than
    // baking a corrupted aspect ratio into the camera - render() retries
    // this every frame, so a real size is picked up as soon as it appears.
    if (width <= 2 || height <= 2) return;
    if (width === lastWidth && height === lastHeight) return;
    lastWidth = width;
    lastHeight = height;
    renderer.setSize(width, height, false);
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
  }

  let dragging = false;
  let lastX = 0;
  let lastY = 0;
  let pointerDownX = 0;
  let pointerDownY = 0;
  const CLICK_MOVE_THRESHOLD_PX = 4;

  function pickCellAt(clientX, clientY, additive) {
    const rect = renderer.domElement.getBoundingClientRect();
    pointerNdc.x = ((clientX - rect.left) / rect.width) * 2 - 1;
    pointerNdc.y = -((clientY - rect.top) / rect.height) * 2 + 1;
    raycaster.setFromCamera(pointerNdc, camera);
    const hits = raycaster.intersectObjects(
      cellPool.map((cell) => cell.group),
      true
    );
    if (!hits.length) {
      if (!additive) selected = null;
      return;
    }
    let node = hits[0].object;
    while (node && node.parent && node.parent !== scene) node = node.parent;
    const index = cellPool.findIndex((cell) => cell.group === node);
    if (index === -1) {
      if (!additive) selected = null;
      return;
    }
    const picked = { row: Math.floor(index / poolN), i: index % poolN };
    if (additive) {
      // Shift-click builds a batch selection for differential dilation
      // (apply a role/alpha to many cells at once) without disturbing the
      // single-cell inspector in the Selected Cell panel.
      const key = cellKey(picked.row, picked.i);
      if (multiSelected.has(key)) multiSelected.delete(key);
      else multiSelected.add(key);
      updateBatchSelectionUi();
    } else {
      selected = picked;
    }
  }

  renderer.domElement.addEventListener("pointerdown", (event) => {
    dragging = true;
    lastX = event.clientX;
    lastY = event.clientY;
    pointerDownX = event.clientX;
    pointerDownY = event.clientY;
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
    const moved = Math.hypot(event.clientX - pointerDownX, event.clientY - pointerDownY);
    if (moved <= CLICK_MOVE_THRESHOLD_PX) {
      pickCellAt(event.clientX, event.clientY, event.shiftKey);
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
      cameraState.radius = Math.max(60, Math.min(900, cameraState.radius * (event.deltaY > 0 ? 1.08 : 0.92)));
      updateCamera();
    },
    { passive: false }
  );

  document.querySelectorAll("[data-view]").forEach((button) => {
    button.addEventListener("click", () => setView(button.dataset.view));
  });
  document.getElementById("resetView").addEventListener("click", () => setView("iso"));

  focusToggle.addEventListener("click", () => {
    const active = document.body.classList.toggle("focus-mode");
    focusToggle.textContent = active ? "Controls" : "Focus";
  });

  document.getElementById("resetAllBtn").addEventListener("click", () => {
    // .defaultValue/.defaultChecked reflect each input's original HTML
    // attribute, so this always matches what's actually declared in
    // index.html rather than a second, driftable copy of the same numbers.
    [alphaInput, backlashInput, ringCountInput, rowCountInput, axialPitchInput, targetDiameterInput].forEach((input) => {
      input.value = input.defaultValue;
    });
    [aSiteInput, bSiteInput].forEach((select) => {
      select.value = Array.from(select.options).find((option) => option.defaultSelected).value;
    });
    [animateInput, collisionEnabledInput, heatmapEnabledInput].forEach((input) => {
      input.checked = input.defaultChecked;
    });
    clearAllRoles();
    selected = null;
    multiSelected.clear();
    updateBatchSelectionUi();
    if (isolateActive) {
      isolateActive = false;
      isolateCellBtn.textContent = "Isolate Cell";
    }
    if (document.body.classList.contains("focus-mode")) {
      document.body.classList.remove("focus-mode");
      focusToggle.textContent = "Focus";
    }
    stopSequence();
    setView("iso");
    updateMechanism();
  });

  frameCellBtn.addEventListener("click", () => {
    if (!selected) return;
    cameraState.target.copy(selectedWorld);
    updateCamera();
  });

  isolateCellBtn.addEventListener("click", () => {
    if (!selected) return;
    isolateActive = !isolateActive;
    isolateCellBtn.textContent = isolateActive ? "Show Lattice" : "Isolate Cell";
  });

  cellRoleSelect.addEventListener("change", () => {
    if (!selected) return;
    const { row, i } = selected;
    const newRole = cellRoleSelect.value;
    if (newRole !== "free" && cellRoles[row][i].role === "free") {
      // Seed the actuator/lock value from the cell's last computed alpha so
      // switching roles doesn't cause a visible jump in the mechanism.
      cellRoles[row][i].alpha = clamp(Number(cellAlphaInput.value) || ALPHA_REFERENCE, ALPHA_MIN, ALPHA_MAX);
    }
    cellRoles[row][i].role = newRole;
    updateMechanism();
  });

  cellAlphaInput.addEventListener("input", () => {
    if (!selected) return;
    const { row, i } = selected;
    if (cellRoles[row][i].role === "free") return;
    cellRoles[row][i].alpha = clamp(Number(cellAlphaInput.value), ALPHA_MIN, ALPHA_MAX);
    updateMechanism();
  });

  function clearAllRoles() {
    for (let row = 0; row < cellRoles.length; row += 1) {
      for (let i = 0; i < cellRoles[row].length; i += 1) {
        cellRoles[row][i] = { role: "free", alpha: ALPHA_REFERENCE };
      }
    }
  }

  clearRolesBtn.addEventListener("click", () => {
    clearAllRoles();
    updateMechanism();
  });

  function updateBatchSelectionUi() {
    const count = multiSelected.size;
    batchSelectedCountMetric.textContent = String(count);
    applyBatchBtn.disabled = count === 0;
    clearBatchSelectionBtn.disabled = count === 0;
  }

  batchAlphaInput.addEventListener("input", () => {
    batchAlphaOut.textContent = Number(batchAlphaInput.value).toFixed(2);
  });

  applyBatchBtn.addEventListener("click", () => {
    if (!multiSelected.size) return;
    const role = batchRoleSelect.value;
    const alpha = clamp(Number(batchAlphaInput.value), ALPHA_MIN, ALPHA_MAX);
    multiSelected.forEach((key) => {
      const [row, i] = key.split(",").map(Number);
      if (!cellRoles[row] || !cellRoles[row][i]) return;
      cellRoles[row][i] = { role, alpha: role === "free" ? ALPHA_REFERENCE : alpha };
    });
    updateMechanism();
  });

  clearBatchSelectionBtn.addEventListener("click", () => {
    multiSelected.clear();
    updateBatchSelectionUi();
    updateMechanism();
  });

  // Shape presets: command every cell in a row to the same alpha, but vary
  // that alpha row-to-row, so each ring settles at a different diameter
  // (rows already solve independently - see computeCellAlphas) producing a
  // visibly non-uniform wheel profile in one click, as a fast demo of
  // differential dilation on top of the shift-click batch tool above.
  function applyRowAlphaPreset(rowAlphaFn) {
    const m = cellRoles.length;
    for (let row = 0; row < m; row += 1) {
      const n = cellRoles[row].length;
      const alpha = clamp(rowAlphaFn(row, m), ALPHA_MIN, ALPHA_MAX);
      for (let i = 0; i < n; i += 1) {
        cellRoles[row][i] = { role: "actuator", alpha };
      }
    }
    updateMechanism();
  }

  presetBarrelBtn.addEventListener("click", () => {
    applyRowAlphaPreset((row, m) => {
      const mid = (m - 1) / 2;
      const t = mid === 0 ? 0 : 1 - Math.abs(row - mid) / mid;
      return ALPHA_REFERENCE + t * 0.5;
    });
  });

  presetConeBtn.addEventListener("click", () => {
    applyRowAlphaPreset((row, m) => {
      const t = m <= 1 ? 0 : row / (m - 1);
      return 0.75 + t * 0.7;
    });
  });

  presetSaddleBtn.addEventListener("click", () => {
    applyRowAlphaPreset((row, m) => {
      const mid = (m - 1) / 2;
      const t = mid === 0 ? 0 : 1 - Math.abs(row - mid) / mid;
      return ALPHA_REFERENCE - t * 0.4;
    });
  });

  function applyDiameterFit() {
    const n = clamp(Math.round(Number(ringCountInput.value)), RING_COUNT_MIN, RING_COUNT_MAX);
    const target = Number(targetDiameterInput.value);
    const backlash = Number(backlashInput.value);
    const fit = fitAlphaToDiameter(target, n, backlash);
    clearAllRoles();
    alphaInput.value = String(fit.alpha);
    achievedDiameterMetric.textContent = `${fit.achievedDiameter.toFixed(1)} mm`;
    fitResidualMetric.textContent = `${(fit.achievedDiameter - target).toFixed(1)} mm`;
    updateMechanism();
  }

  targetDiameterInput.addEventListener("input", () => {
    targetDiameterOut.textContent = `${targetDiameterInput.value} mm`;
    // Live-apply like every other slider in the app - previously this only
    // updated its own label, and the wheel would not actually change until
    // the separate Fit Alpha button was clicked, which read as "target
    // diameter doesn't do anything" since nothing else here works that way.
    applyDiameterFit();
  });

  fitDiameterBtn.addEventListener("click", () => {
    applyDiameterFit();
  });

  [
    alphaInput,
    backlashInput,
    ringCountInput,
    rowCountInput,
    axialPitchInput,
    aSiteInput,
    bSiteInput,
    animateInput,
    collisionEnabledInput,
    heatmapEnabledInput,
    showMeasurementsInput,
  ].forEach((input) => {
    input.addEventListener("input", () => updateMechanism());
    input.addEventListener("change", () => updateMechanism());
  });

  const SAVE_FORMAT = "rad-cylinder-tiling.v1";

  function currentStateSnapshot() {
    return {
      format: SAVE_FORMAT,
      alpha: Number(alphaInput.value),
      backlash: Number(backlashInput.value),
      ringCount: Number(ringCountInput.value),
      rowCount: Number(rowCountInput.value),
      axialPitch: Number(axialPitchInput.value),
      aSite: aSiteInput.value,
      bSite: bSiteInput.value,
      cellRoles,
    };
  }

  function applyStateSnapshot(state) {
    if (!state || typeof state !== "object") return;
    if (Number.isFinite(state.alpha)) alphaInput.value = String(clamp(state.alpha, ALPHA_MIN, ALPHA_MAX));
    if (Number.isFinite(state.backlash)) backlashInput.value = String(clamp(state.backlash, Number(backlashInput.min), Number(backlashInput.max)));
    if (Number.isFinite(state.ringCount)) ringCountInput.value = String(clamp(Math.round(state.ringCount), RING_COUNT_MIN, RING_COUNT_MAX));
    if (Number.isFinite(state.rowCount)) rowCountInput.value = String(clamp(Math.round(state.rowCount), ROW_COUNT_MIN, ROW_COUNT_MAX));
    if (Number.isFinite(state.axialPitch)) axialPitchInput.value = String(clamp(state.axialPitch, Number(axialPitchInput.min), Number(axialPitchInput.max)));
    if (state.aSite && SITE_VECTORS[state.aSite]) aSiteInput.value = state.aSite;
    if (state.bSite && SITE_VECTORS[state.bSite]) bSiteInput.value = state.bSite;

    const n = clamp(Math.round(Number(ringCountInput.value)), RING_COUNT_MIN, RING_COUNT_MAX);
    const m = clamp(Math.round(Number(rowCountInput.value)), ROW_COUNT_MIN, ROW_COUNT_MAX);
    ensureRoles(n, m);
    if (Array.isArray(state.cellRoles)) {
      for (let row = 0; row < m && row < state.cellRoles.length; row += 1) {
        const savedRow = state.cellRoles[row];
        if (!Array.isArray(savedRow)) continue;
        for (let i = 0; i < n && i < savedRow.length; i += 1) {
          const savedCell = savedRow[i];
          if (!savedCell || typeof savedCell !== "object") continue;
          const role = ["free", "actuator", "locked"].includes(savedCell.role) ? savedCell.role : "free";
          const alphaValue = Number.isFinite(savedCell.alpha) ? clamp(savedCell.alpha, ALPHA_MIN, ALPHA_MAX) : ALPHA_REFERENCE;
          cellRoles[row][i] = { role, alpha: alphaValue };
        }
      }
    }
    updateMechanism();
  }

  // Serializes the actual rendered geometry (not a re-derivation of it) of
  // every currently-visible cell into a Wavefront OBJ: walks each cell's
  // Three.js meshes (arms, hubs, pads, hole placeholders - not the role
  // marker or selection marker, since those live outside cell.group and so
  // are naturally skipped by traversing only cell.group), transforms their
  // vertices into world-space millimeters via the already-computed
  // matrixWorld, and writes triangle faces. This does not perform boolean
  // subtraction - hole locations remain solid placeholder cylinders, same
  // as they are on screen - so it's a geometry/measurement reference for
  // CAD, not a print-ready model.
  function buildObjText() {
    scene.updateMatrixWorld(true);
    const lines = ["# Cylindrical Wheel Tiling Simulator - geometry export", `# ${new Date().toISOString()}`];
    let vertexOffset = 0;
    const vertex = new THREE.Vector3();
    cellPool.forEach((cell, cellIndex) => {
      if (!cell.group.visible) return;
      lines.push(`o cell_${cellIndex}`);
      cell.group.traverse((object) => {
        if (!object.isMesh) return;
        const geometry = object.geometry;
        const position = geometry.attributes.position;
        if (!position) return;
        for (let vi = 0; vi < position.count; vi += 1) {
          vertex.fromBufferAttribute(position, vi);
          vertex.applyMatrix4(object.matrixWorld);
          lines.push(`v ${vertex.x.toFixed(4)} ${vertex.y.toFixed(4)} ${vertex.z.toFixed(4)}`);
        }
        const index = geometry.index;
        if (index) {
          for (let fi = 0; fi < index.count; fi += 3) {
            const a = index.getX(fi) + 1 + vertexOffset;
            const b = index.getX(fi + 1) + 1 + vertexOffset;
            const c = index.getX(fi + 2) + 1 + vertexOffset;
            lines.push(`f ${a} ${b} ${c}`);
          }
        } else {
          for (let fi = 0; fi + 2 < position.count; fi += 3) {
            lines.push(`f ${fi + 1 + vertexOffset} ${fi + 2 + vertexOffset} ${fi + 3 + vertexOffset}`);
          }
        }
        vertexOffset += position.count;
      });
    });
    return lines.join("\n");
  }

  exportObjBtn.addEventListener("click", () => {
    const blob = new Blob([buildObjText()], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = "rad-cylinder-tiling.obj";
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  });

  saveJsonBtn.addEventListener("click", () => {
    const blob = new Blob([JSON.stringify(currentStateSnapshot(), null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = "rad-cylinder-tiling.json";
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  });

  loadJsonBtn.addEventListener("click", () => loadJsonFile.click());
  loadJsonFile.addEventListener("change", () => {
    const file = loadJsonFile.files && loadJsonFile.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      try {
        applyStateSnapshot(JSON.parse(String(reader.result)));
      } catch (_error) {
        window.alert("Could not parse that JSON file.");
      }
    };
    reader.readAsText(file);
    loadJsonFile.value = "";
  });

  // --- Timeline: named keyframes are full state snapshots (same shape as
  // Save JSON), applied as discrete jumps - no interpolation between poses.
  const SEQUENCE_FORMAT = "rad-cylinder-tiling-sequence.v1";
  let keyframes = [];
  let activeKeyframeIndex = -1;
  let playbackTimer = null;

  function stopSequence() {
    if (playbackTimer !== null) {
      clearInterval(playbackTimer);
      playbackTimer = null;
    }
    playSequenceBtn.textContent = "Play Sequence";
  }

  function goToKeyframe(index) {
    if (index < 0 || index >= keyframes.length) return;
    activeKeyframeIndex = index;
    applyStateSnapshot(keyframes[index].snapshot);
    renderKeyframeList();
  }

  function renderKeyframeList() {
    keyframeList.innerHTML = "";
    keyframes.forEach((keyframe, index) => {
      const row = document.createElement("div");
      row.className = "keyframe-row" + (index === activeKeyframeIndex ? " active" : "");
      const label = document.createElement("span");
      label.textContent = keyframe.label;
      const goBtn = document.createElement("button");
      goBtn.type = "button";
      goBtn.textContent = "Go";
      goBtn.addEventListener("click", () => {
        stopSequence();
        goToKeyframe(index);
      });
      const deleteBtn = document.createElement("button");
      deleteBtn.type = "button";
      deleteBtn.textContent = "×";
      deleteBtn.setAttribute("aria-label", `Delete ${keyframe.label}`);
      deleteBtn.addEventListener("click", () => {
        stopSequence();
        keyframes.splice(index, 1);
        if (activeKeyframeIndex === index) activeKeyframeIndex = -1;
        else if (activeKeyframeIndex > index) activeKeyframeIndex -= 1;
        renderKeyframeList();
      });
      row.appendChild(label);
      row.appendChild(goBtn);
      row.appendChild(deleteBtn);
      keyframeList.appendChild(row);
    });
    keyframeCountMetric.textContent = String(keyframes.length);
    playSequenceBtn.disabled = keyframes.length === 0;
    saveSequenceBtn.disabled = keyframes.length === 0;
  }

  captureKeyframeBtn.addEventListener("click", () => {
    keyframes.push({ label: `Keyframe ${keyframes.length + 1}`, snapshot: currentStateSnapshot() });
    activeKeyframeIndex = keyframes.length - 1;
    renderKeyframeList();
  });

  playSequenceBtn.addEventListener("click", () => {
    if (playbackTimer !== null) {
      stopSequence();
      return;
    }
    if (keyframes.length === 0) return;
    let index = 0;
    goToKeyframe(index);
    playSequenceBtn.textContent = "Stop";
    playbackTimer = setInterval(() => {
      index += 1;
      if (index >= keyframes.length) {
        stopSequence();
        return;
      }
      goToKeyframe(index);
    }, 1200);
  });

  saveSequenceBtn.addEventListener("click", () => {
    const payload = { format: SEQUENCE_FORMAT, keyframes };
    const blob = new Blob([JSON.stringify(payload, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = "rad-cylinder-tiling-sequence.json";
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  });

  loadSequenceBtn.addEventListener("click", () => loadSequenceFile.click());
  loadSequenceFile.addEventListener("change", () => {
    const file = loadSequenceFile.files && loadSequenceFile.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      try {
        const payload = JSON.parse(String(reader.result));
        if (!Array.isArray(payload.keyframes)) throw new Error("missing keyframes array");
        stopSequence();
        keyframes = payload.keyframes
          .filter((k) => k && typeof k === "object" && k.snapshot)
          .map((k, index) => ({ label: typeof k.label === "string" ? k.label : `Keyframe ${index + 1}`, snapshot: k.snapshot }));
        activeKeyframeIndex = -1;
        renderKeyframeList();
      } catch (_error) {
        window.alert("Could not parse that sequence JSON file.");
      }
    };
    reader.readAsText(file);
    loadSequenceFile.value = "";
  });

  renderKeyframeList();

  new ResizeObserver(resize).observe(mount);
  resize();
  setView("iso");

  // Minimal read-only hook for automated testing: the per-cell alpha grid
  // isn't otherwise DOM-observable, and guessing screen coordinates to
  // click a specific (row, i) cell is fragile across camera angles/
  // viewports. Exposes no way to mutate simulator state.
  window.__cylinderTilingDebug = {
    getCellAlphas: () => lastCellAlphas,
    getCellRoles: () => lastCellRolesSnapshot,
    getSelected: () => selected,
    getCellTopColor: (row, i) => {
      const n = poolN;
      const cell = cellPool[row * n + i];
      return cell ? `#${cell.topMaterial.color.getHexString()}` : null;
    },
    getCellGroupNormal: (row, i) => {
      // World-space direction of the cell group's local Z axis (its
      // "thickness"/face-normal direction) - radially outward once
      // oriented correctly, (0,0,1)-ish when still axial-facing.
      const n = poolN;
      const cell = cellPool[row * n + i];
      if (!cell) return null;
      const normal = new THREE.Vector3(0, 0, 1).applyQuaternion(cell.group.quaternion);
      return { x: normal.x, y: normal.y, z: normal.z };
    },
    getRingGeometry: (row) => {
      // Exposes each ring's actual center positions and centroid so tests
      // can verify cell orientation against the ring's true geometric
      // outward direction, not just against the orientation code's own
      // inputs (bottomRot is the chain-walking heading, not the polar
      // position angle, and trusting it directly was the root cause of a
      // real bug where cells ended up rotated 90-ish degrees off).
      const r = lastRingsSnapshot && lastRingsSnapshot[row];
      if (!r) return null;
      return {
        centroid: { x: r.centroid.x, y: r.centroid.y },
        centers: r.centers.map((c) => ({ x: c.x, y: c.y })),
      };
    },
  };

  function render(timeMs) {
    resize();
    updateMechanism(timeMs);
    renderer.render(scene, camera);
    requestAnimationFrame(render);
  }

  requestAnimationFrame(render);
})();
