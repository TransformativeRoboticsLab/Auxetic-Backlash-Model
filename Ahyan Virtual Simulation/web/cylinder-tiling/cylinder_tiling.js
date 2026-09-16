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
  function checkNeighborClearance(n, m, ring, axialPitch, aTopRad) {
    function bodiesForCell(row, i) {
      const center = ring.centers[i];
      const z = row * axialPitch;
      const origin2 = new THREE.Vector2(center.x, center.y);
      return [
        { cell: `${row}-${i}-lower`, origin: origin2, rotation: ring.bottomRot[i], zMin: z - CAD.bodyThicknessMm, zMax: z },
        { cell: `${row}-${i}-upper`, origin: origin2, rotation: ring.bottomRot[i] + aTopRad, zMin: z, zMax: z + CAD.bodyThicknessMm },
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

  function rebuildPool(n, m) {
    cellPool.forEach((cell) => {
      scene.remove(cell.group);
      disposeGroup(cell.group);
    });
    cellPool = [];
    for (let row = 0; row < m; row += 1) {
      const mats = materialsForRow(row);
      for (let i = 0; i < n; i += 1) {
        const cell = createCell(`ring${row}-cell${i}`, mats.top, mats.bottom);
        scene.add(cell.group);
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

  // Close n identical single-pin joints into a ring: cell i's upper-cross hole
  // (aSite, at absolute rotation bottomRot[i] + aTopRad) meets cell i+1's
  // lower-cross hole (bSite, at absolute rotation bottomRot[i+1]), where each
  // step turns by 2*pi/n. This is the same forward relation
  // two_cell_attachment.js uses for its A-B pin (aAttachWorld - bAttachOffset),
  // just walked n times with a fixed turn increment instead of laid flat -
  // a regular n-gon closes exactly by construction, which the residual below verifies.
  function buildRing(n, aTopRad) {
    const aAttachLocal = SITE_VECTORS[aSiteInput.value] || SITE_VECTORS.east;
    const bAttachLocal = SITE_VECTORS[bSiteInput.value] || SITE_VECTORS.west;
    const turn = (2 * Math.PI) / n;
    const centers = [new THREE.Vector2(0, 0)];
    const bottomRot = [0];
    let heading = 0;
    let pitch = 0;
    let closureResidual = 0;
    for (let i = 0; i < n; i += 1) {
      const aAttachWorld = rotate2(aAttachLocal, heading + aTopRad);
      heading += turn;
      const bAttachWorld = rotate2(bAttachLocal, heading);
      const step = aAttachWorld.clone().sub(bAttachWorld);
      const nextCenter = centers[i].clone().add(step);
      if (i === 0) pitch = step.length();
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
    return { centers, bottomRot, pitch, closureResidual, radius, diameter: 2 * radius };
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

  const raycaster = new THREE.Raycaster();
  const pointerNdc = new THREE.Vector2();
  let selected = null; // { row, i }

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
    const alphaEffective = ALPHA_REFERENCE + reluDeadzone(alphaCommand - ALPHA_REFERENCE, backlash);
    const thetaDeg = clamp(alphaToThetaDeg(alphaEffective), THETA_SAFETY_MIN_DEG, THETA_SAFETY_MAX_DEG);
    const aTopRad = degToRad(thetaDeg);
    const inDeadzone = Math.abs(alphaCommand - ALPHA_REFERENCE) <= backlash;

    if (n !== poolN || m !== poolM) {
      rebuildPool(n, m);
      updateLegend(m);
      if (selected && (selected.row >= m || selected.i >= n)) {
        selected = null;
      }
    }

    const ring = buildRing(n, aTopRad);

    for (let row = 0; row < m; row += 1) {
      const z = row * axialPitch;
      for (let i = 0; i < n; i += 1) {
        const cell = cellPool[row * n + i];
        const center = ring.centers[i];
        cell.group.position.set(center.x, center.y, z);
        cell.bottom.rotation.z = ring.bottomRot[i];
        cell.top.rotation.z = ring.bottomRot[i] + aTopRad;
      }
    }

    alphaOut.textContent = alphaCommand.toFixed(2);
    backlashOut.textContent = backlash.toFixed(2);
    alphaCommandMetric.textContent = alphaCommand.toFixed(2);
    alphaEffectiveMetric.textContent = alphaEffective.toFixed(2);
    thetaMetric.textContent = `${thetaDeg.toFixed(1)} deg`;
    deadzoneState.textContent = inDeadzone ? "free (dead zone)" : "engaged";
    deadzoneState.classList.toggle("status-adjusted", inDeadzone);
    deadzoneState.classList.toggle("status-ok", !inDeadzone);

    ringCountOut.textContent = String(n);
    rowCountOut.textContent = String(m);
    axialPitchOut.textContent = `${axialPitch.toFixed(0)} mm`;
    pitchMetric.textContent = `${ring.pitch.toFixed(1)} mm`;
    turnMetric.textContent = `${radToDeg((2 * Math.PI) / n).toFixed(1)} deg`;
    diameterMetric.textContent = `${ring.diameter.toFixed(1)} mm`;
    radiusMetric.textContent = `${ring.radius.toFixed(1)} mm`;
    closureMetric.textContent = `${ring.closureResidual.toFixed(4)} mm`;
    totalCellsMetric.textContent = String(n * m);
    heightMetric.textContent = `${((m - 1) * axialPitch).toFixed(1)} mm`;

    if (collisionEnabledInput.checked) {
      const report = checkNeighborClearance(n, m, ring, axialPitch, aTopRad);
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
      const center = ring.centers[selected.i];
      const z = selected.row * axialPitch;
      selectionMarker.visible = true;
      selectionMarker.position.set(center.x, center.y, z);
      selectedCellStatus.textContent = `row ${selected.row}, cell ${selected.i}`;
      selectedCellStatus.classList.add("status-ok");
      selectedIndexMetric.textContent = `${selected.row}, ${selected.i}`;
      selectedCenterMetric.textContent = `(${center.x.toFixed(1)}, ${center.y.toFixed(1)}, ${z.toFixed(1)})`;
      selectedBottomRotMetric.textContent = `${radToDeg(ring.bottomRot[selected.i]).toFixed(1)} deg`;
      selectedTopRotMetric.textContent = `${radToDeg(ring.bottomRot[selected.i] + aTopRad).toFixed(1)} deg`;
    } else {
      selectionMarker.visible = false;
      selectedCellStatus.textContent = "no cell selected";
      selectedCellStatus.classList.remove("status-ok");
      selectedIndexMetric.textContent = "-";
      selectedCenterMetric.textContent = "-";
      selectedBottomRotMetric.textContent = "-";
      selectedTopRotMetric.textContent = "-";
    }
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

  function pickCellAt(clientX, clientY) {
    const rect = renderer.domElement.getBoundingClientRect();
    pointerNdc.x = ((clientX - rect.left) / rect.width) * 2 - 1;
    pointerNdc.y = -((clientY - rect.top) / rect.height) * 2 + 1;
    raycaster.setFromCamera(pointerNdc, camera);
    const hits = raycaster.intersectObjects(
      cellPool.map((cell) => cell.group),
      true
    );
    if (!hits.length) {
      selected = null;
      return;
    }
    let node = hits[0].object;
    while (node && node.parent && node.parent !== scene) node = node.parent;
    const index = cellPool.findIndex((cell) => cell.group === node);
    if (index === -1) {
      selected = null;
      return;
    }
    selected = { row: Math.floor(index / poolN), i: index % poolN };
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
      pickCellAt(event.clientX, event.clientY);
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
    updateMechanism();
  }

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

  new ResizeObserver(resize).observe(mount);
  resize();
  setView("iso");

  function render(timeMs) {
    resize();
    updateMechanism(timeMs);
    renderer.render(scene, camera);
    requestAnimationFrame(render);
  }

  requestAnimationFrame(render);
})();
