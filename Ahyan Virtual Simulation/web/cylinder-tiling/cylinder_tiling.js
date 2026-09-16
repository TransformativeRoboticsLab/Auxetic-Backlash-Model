(function () {
  "use strict";

  const THREE = window.THREE;
  const mount = document.getElementById("threeMount");

  const aTopInput = document.getElementById("aTopAngle");
  const aTopOut = document.getElementById("aTopOut");
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

  const A_TOP_MIN_DEG = Number(aTopInput.min);
  const A_TOP_MAX_DEG = Number(aTopInput.max);
  const RING_COUNT_MIN = Number(ringCountInput.min);
  const RING_COUNT_MAX = Number(ringCountInput.max);
  const ROW_COUNT_MIN = Number(rowCountInput.min);
  const ROW_COUNT_MAX = Number(rowCountInput.max);

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

    let aTopDeg = Number(aTopInput.value);
    if (animateInput.checked) {
      const mid = (A_TOP_MIN_DEG + A_TOP_MAX_DEG) * 0.5;
      const amplitude = (A_TOP_MAX_DEG - A_TOP_MIN_DEG) * 0.5;
      aTopDeg = Math.round(mid + amplitude * Math.sin(timeMs * 0.0005));
      aTopInput.value = String(aTopDeg);
    }
    const aTopRad = degToRad(aTopDeg);

    if (n !== poolN || m !== poolM) {
      rebuildPool(n, m);
      updateLegend(m);
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

    aTopOut.textContent = `${aTopDeg.toFixed(0)} deg`;
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
      cameraState.radius = Math.max(60, Math.min(900, cameraState.radius * (event.deltaY > 0 ? 1.08 : 0.92)));
      updateCamera();
    },
    { passive: false }
  );

  document.querySelectorAll("[data-view]").forEach((button) => {
    button.addEventListener("click", () => setView(button.dataset.view));
  });
  document.getElementById("resetView").addEventListener("click", () => setView("iso"));
  [aTopInput, ringCountInput, rowCountInput, axialPitchInput, aSiteInput, bSiteInput, animateInput].forEach((input) => {
    input.addEventListener("input", () => updateMechanism());
    input.addEventListener("change", () => updateMechanism());
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
