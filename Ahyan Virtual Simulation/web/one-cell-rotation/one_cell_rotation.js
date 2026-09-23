(function () {
  "use strict";

  const THREE = window.THREE;
  const mount = document.getElementById("threeMount");
  const angleInput = document.getElementById("angle");
  const animateInput = document.getElementById("animate");
  const angleOut = document.getElementById("angleOut");
  const thetaMetric = document.getElementById("thetaMetric");
  const arcMetric = document.getElementById("arcMetric");
  const padMetric = document.getElementById("padMetric");

  const CAD = Object.freeze({
    cellWidthMm: 55.604331,
    nominalHoleDiameterMm: 3.4,
    bodyThicknessMm: 4.0,
    padRadiusMm: 4.9,
    hubRadiusMm: 4.6,
    armWidthMm: 5.4,
    siteRadiusMm: 22.1,
  });

  if (!THREE || !mount) {
    if (mount) mount.textContent = "Three.js did not load.";
    return;
  }

  const scene = new THREE.Scene();
  scene.background = null;

  const camera = new THREE.PerspectiveCamera(42, 1, 0.1, 600);
  camera.up.set(0, 0, 1);

  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.shadowMap.enabled = true;
  mount.appendChild(renderer.domElement);

  const cameraState = {
    radius: 92,
    azimuth: -0.78,
    elevation: 0.55,
    target: new THREE.Vector3(0, 0, 0),
  };

  const materials = {
    top: new THREE.MeshStandardMaterial({ color: 0xda7b27, roughness: 0.48, metalness: 0.15 }),
    bottom: new THREE.MeshStandardMaterial({ color: 0x44515f, roughness: 0.55, metalness: 0.12 }),
    edge: new THREE.LineBasicMaterial({ color: 0x121820, transparent: true, opacity: 0.56 }),
    hole: new THREE.MeshStandardMaterial({ color: 0x15191f, roughness: 0.7, metalness: 0.05 }),
    axis: new THREE.MeshStandardMaterial({ color: 0x25a36f, roughness: 0.35, metalness: 0.2 }),
    trace: new THREE.LineBasicMaterial({ color: 0x1d6fd6, transparent: true, opacity: 0.9 }),
    ghost: new THREE.MeshBasicMaterial({ color: 0x1d6fd6, transparent: true, opacity: 0.36 }),
  };

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

    const hole = cylinderZ(CAD.nominalHoleDiameterMm * 0.5, CAD.bodyThicknessMm + 0.22, materials.hole, 40);
    hole.position.set(x, y, 0.03);
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
    centerHole.position.z = 0.04;
    group.add(centerHole);

    addPad(group, CAD.siteRadiusMm, 0, material);
    addPad(group, -CAD.siteRadiusMm, 0, material);
    addPad(group, 0, CAD.siteRadiusMm, material);
    addPad(group, 0, -CAD.siteRadiusMm, material);
    return group;
  }

  const bottomCross = createCrossPart("fixed lower cross", materials.bottom);
  const topCross = createCrossPart("driven upper cross", materials.top);
  scene.add(bottomCross);
  scene.add(topCross);

  const centerAxis = cylinderZ(CAD.nominalHoleDiameterMm * 0.34, 34, materials.axis, 40);
  centerAxis.position.z = 0;
  scene.add(centerAxis);

  const startMarker = cylinderZ(1.35, 1.0, materials.ghost, 24);
  startMarker.position.set(CAD.siteRadiusMm, 0, 0);
  scene.add(startMarker);

  const movingMarker = cylinderZ(1.6, 1.2, materials.axis, 24);
  scene.add(movingMarker);

  const arcLine = new THREE.Line(new THREE.BufferGeometry(), materials.trace);
  scene.add(arcLine);

  function addAxisLine(points, color) {
    const geometry = new THREE.BufferGeometry().setFromPoints(points.map((p) => new THREE.Vector3(...p)));
    scene.add(new THREE.Line(geometry, new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.9 })));
  }

  addAxisLine(
    [
      [-34, 0, -8],
      [34, 0, -8],
    ],
    0xb24d3a
  );
  addAxisLine(
    [
      [0, -34, -8],
      [0, 34, -8],
    ],
    0x1d6fd6
  );
  addAxisLine(
    [
      [0, 0, -12],
      [0, 0, 24],
    ],
    0x25a36f
  );

  const grid = new THREE.GridHelper(82, 10, 0x8b97a5, 0xc4ccd4);
  grid.rotation.x = Math.PI / 2;
  grid.position.z = -8.1;
  scene.add(grid);

  const hemi = new THREE.HemisphereLight(0xffffff, 0x5f6974, 1.7);
  scene.add(hemi);

  const key = new THREE.DirectionalLight(0xffffff, 2.1);
  key.position.set(28, -38, 70);
  key.castShadow = true;
  key.shadow.mapSize.width = 1024;
  key.shadow.mapSize.height = 1024;
  scene.add(key);

  const fill = new THREE.DirectionalLight(0xb8d4ff, 0.65);
  fill.position.set(-35, 40, 35);
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
      cameraState.radius = 88;
    } else if (view === "front") {
      cameraState.azimuth = -Math.PI / 2;
      cameraState.elevation = 0.02;
      cameraState.radius = 92;
    } else if (view === "side") {
      cameraState.azimuth = 0;
      cameraState.elevation = 0.05;
      cameraState.radius = 92;
    } else {
      cameraState.azimuth = -0.78;
      cameraState.elevation = 0.55;
      cameraState.radius = 92;
    }
    updateCamera();
    document.querySelectorAll("[data-view]").forEach((button) => {
      button.setAttribute("aria-pressed", String(button.dataset.view === view));
    });
  }

  function updateMechanism(timeMs = 0) {
    if (animateInput.checked) {
      const animatedAngle = 52 * Math.sin(timeMs * 0.001);
      angleInput.value = animatedAngle.toFixed(0);
    }
    const angleDeg = Number(angleInput.value);
    const theta = (angleDeg * Math.PI) / 180;
    const topZ = CAD.bodyThicknessMm * 0.5;
    const bottomZ = -CAD.bodyThicknessMm * 0.5;

    topCross.position.z = topZ;
    bottomCross.position.z = bottomZ;
    topCross.rotation.z = theta;
    bottomCross.rotation.z = 0;
    movingMarker.position.set(CAD.siteRadiusMm * Math.cos(theta), CAD.siteRadiusMm * Math.sin(theta), topZ + 2.8);
    startMarker.position.z = topZ + 2.4;

    const steps = Math.max(8, Math.ceil(Math.abs(angleDeg) / 3));
    const points = [];
    for (let i = 0; i <= steps; i += 1) {
      const a = (theta * i) / steps;
      points.push(new THREE.Vector3(CAD.siteRadiusMm * Math.cos(a), CAD.siteRadiusMm * Math.sin(a), topZ + 2.7));
    }
    arcLine.geometry.dispose();
    arcLine.geometry = new THREE.BufferGeometry().setFromPoints(points);

    const rimArc = Math.abs(theta) * CAD.siteRadiusMm;
    const padOffset = 2 * CAD.siteRadiusMm * Math.abs(Math.sin(theta / 2));
    angleOut.textContent = `${angleDeg.toFixed(0)} deg`;
    thetaMetric.textContent = `${angleDeg.toFixed(1)} deg`;
    arcMetric.textContent = `${rimArc.toFixed(1)} mm`;
    padMetric.textContent = `${padOffset.toFixed(1)} mm`;
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
      // Pointer capture may already be released by the browser.
    }
  });

  renderer.domElement.addEventListener(
    "wheel",
    (event) => {
      event.preventDefault();
      cameraState.radius = Math.max(42, Math.min(180, cameraState.radius * (event.deltaY > 0 ? 1.08 : 0.92)));
      updateCamera();
    },
    { passive: false }
  );

  document.querySelectorAll("[data-view]").forEach((button) => {
    button.addEventListener("click", () => setView(button.dataset.view));
  });
  document.getElementById("resetView").addEventListener("click", () => setView("iso"));
  [angleInput, animateInput].forEach((input) => {
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
