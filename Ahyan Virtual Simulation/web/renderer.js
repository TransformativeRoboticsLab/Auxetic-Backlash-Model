(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  class RadRenderer {
    constructor(mount, onSelect, onViewChange = null, onHover = null) {
      if (!window.THREE) throw new Error("Three.js did not load.");
      this.THREE = window.THREE;
      this.mount = mount;
      this.onSelect = onSelect;
      this.onViewChange = onViewChange;
      this.onHover = onHover;
      this.raycaster = new this.THREE.Raycaster();
      this.pointer = new this.THREE.Vector2();
      this.selectables = [];
      this.hoveredCell = null;
      this.cellGroups = new Map();
      this.latticeConnectors = [];
      this.referenceShape = null;
      this.builtShape = null;
      this.drag = null;
      this.targetState = null;
      this.targetSim = null;
      this.displaySim = null;
      this.lastRebuild = 0;
      this.frameCount = 0;
      this.lastFpsTime = performance.now();
      this.diagnostics = { fps: 0, drawCalls: 0, triangles: 0, geometries: 0, textures: 0, cells: 0, surfaceVertices: 0 };

      this.scene = new this.THREE.Scene();
      this.scene.background = new this.THREE.Color(0xeef3f8);
      this.perspectiveCamera = new this.THREE.PerspectiveCamera(45, 1, 0.1, 200);
      this.orthographicCamera = new this.THREE.OrthographicCamera(-8, 8, 8, -8, 0.1, 200);
      this.camera = this.perspectiveCamera;
      this.perspectiveCamera.up.set(0, 0, 1);
      this.orthographicCamera.up.set(0, 0, 1);
      this.cameraTarget = new this.THREE.Vector3(0, 0, 0);
      this.spherical = { radius: 15.6, theta: -Math.PI / 4, phi: Math.acos(1 / Math.sqrt(3)) };
      this.viewMode = "iso";
      this.projectionMode = "perspective";
      this.viewportAspect = 1;

      this.renderer = new this.THREE.WebGLRenderer({ antialias: true });
      this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
      this.renderer.shadowMap.enabled = true;
      mount.appendChild(this.renderer.domElement);

      this.root = new this.THREE.Group();
      this.referenceRoot = new this.THREE.Group();
      this.displacementVectorRoot = new this.THREE.Group();
      this.cellRoot = new this.THREE.Group();
      this.linkageRoot = new this.THREE.Group();
      this.membraneRoot = new this.THREE.Group();
      this.targetRoot = new this.THREE.Group();
      this.surfaceContourRoot = new this.THREE.Group();
      this.errorRoot = new this.THREE.Group();
      this.surfaceNormalRoot = new this.THREE.Group();
      this.gapRoot = new this.THREE.Group();
      this.actuatorRoot = new this.THREE.Group();
      this.influenceFootprintRoot = new this.THREE.Group();
      this.paintBrushPreviewRoot = new this.THREE.Group();
      this.measurementRoot = new this.THREE.Group();
      this.partCalloutRoot = new this.THREE.Group();
      this.root.add(this.referenceRoot, this.displacementVectorRoot, this.linkageRoot, this.cellRoot, this.membraneRoot, this.targetRoot, this.surfaceContourRoot, this.errorRoot, this.surfaceNormalRoot, this.gapRoot, this.actuatorRoot, this.influenceFootprintRoot, this.paintBrushPreviewRoot, this.measurementRoot, this.partCalloutRoot);
      this.scene.add(this.root);

      this.materials = this.createMaterials();
      this.overlayMaterials = new Map();
      this.geometries = this.createGeometries();
      this.surfaceShape = null;
      this.membraneMesh = null;
      this.targetMesh = null;
      this.buildScene();
      this.bind();
      this.resize();
      this.setView("iso");
      this.animate();
    }

    createMaterials() {
      const T = this.THREE;
      return {
        plate: new T.MeshStandardMaterial({ color: 0x14799d, roughness: 0.42, metalness: 0.08 }),
        selected: new T.MeshStandardMaterial({ color: 0x2d6cdf, roughness: 0.36, metalness: 0.12 }),
        hovered: new T.MeshStandardMaterial({ color: 0xe6b64a, roughness: 0.36, metalness: 0.1 }),
        actuated: new T.MeshStandardMaterial({ color: 0xb44838, roughness: 0.38, metalness: 0.1 }),
        locked: new T.MeshStandardMaterial({ color: 0x667c2d, roughness: 0.5, metalness: 0.06 }),
        abstractBody: new T.MeshStandardMaterial({ color: 0xd9e3ec, roughness: 0.5, metalness: 0.04, transparent: true, opacity: 0.52 }),
        abstractInner: new T.MeshStandardMaterial({ color: 0x1e6b88, roughness: 0.44, metalness: 0.08, transparent: true, opacity: 0.78 }),
        abstractEdge: new T.LineBasicMaterial({ color: 0x172330, transparent: true, opacity: 0.88 }),
        abstractGuide: new T.MeshBasicMaterial({ color: 0x172330, transparent: true, opacity: 0.68 }),
        abstractDatum: new T.MeshBasicMaterial({ color: 0x5d6b78, transparent: true, opacity: 0.32 }),
        paperOuter: new T.MeshStandardMaterial({ color: 0xcbd8e3, roughness: 0.48, metalness: 0.06, transparent: true, opacity: 0.74 }),
        paperInner: new T.MeshStandardMaterial({ color: 0x14799d, roughness: 0.42, metalness: 0.08, transparent: true, opacity: 0.88 }),
        paperClearance: new T.MeshBasicMaterial({ color: 0xe6b64a, transparent: true, opacity: 0.7 }),
        hinge: new T.MeshStandardMaterial({ color: 0x222b36, roughness: 0.32, metalness: 0.35 }),
        brace: new T.MeshStandardMaterial({ color: 0x314152, roughness: 0.38, metalness: 0.22 }),
        linkage: new T.MeshStandardMaterial({ color: 0x253545, roughness: 0.34, metalness: 0.28 }),
        reference: new T.LineBasicMaterial({ color: 0x5d6b78, transparent: true, opacity: 0.42 }),
        displacementVector: new T.LineBasicMaterial({ color: 0xb44838, transparent: true, opacity: 0.72 }),
        displacementHead: new T.MeshBasicMaterial({ color: 0xb44838, transparent: true, opacity: 0.78 }),
        plateEdge: new T.LineBasicMaterial({ color: 0x10202d, transparent: true, opacity: 0.45 }),
        fastener: new T.MeshStandardMaterial({ color: 0x111820, roughness: 0.3, metalness: 0.55 }),
        gap: new T.MeshBasicMaterial({ color: 0xe6b64a, transparent: true, opacity: 0.72 }),
        gapLine: new T.LineBasicMaterial({ color: 0xe6b64a, transparent: true, opacity: 0.85 }),
        stop: new T.MeshStandardMaterial({ color: 0xd2a235, roughness: 0.44, metalness: 0.18 }),
        membrane: new T.MeshStandardMaterial({
          color: 0x8cc8d8,
          transparent: true,
          opacity: 0.32,
          side: T.DoubleSide,
          roughness: 0.25,
          metalness: 0.02,
        }),
        actuator: new T.MeshStandardMaterial({ color: 0xb44838, roughness: 0.28, metalness: 0.25 }),
        actuatorPositive: new T.MeshStandardMaterial({ color: 0xb44838, roughness: 0.26, metalness: 0.22 }),
        actuatorNegative: new T.MeshStandardMaterial({ color: 0x2d6cdf, roughness: 0.26, metalness: 0.22 }),
        target: new T.MeshStandardMaterial({
          color: 0x6f5cc2,
          transparent: true,
          opacity: 0.22,
          side: T.DoubleSide,
          roughness: 0.35,
        }),
        measurement: new T.MeshBasicMaterial({ color: 0x202936 }),
        contour: new T.LineBasicMaterial({ color: 0x202936, transparent: true, opacity: 0.42 }),
        targetContour: new T.LineBasicMaterial({ color: 0x6f5cc2, transparent: true, opacity: 0.72 }),
        selection: new T.LineBasicMaterial({ color: 0x2d6cdf, linewidth: 2 }),
        actuatorEnvelope: new T.MeshBasicMaterial({ color: 0xb44838, transparent: true, opacity: 0.16 }),
        actuatorGaugeTick: new T.MeshBasicMaterial({ color: 0x202936, transparent: true, opacity: 0.62 }),
        errorPositive: new T.LineBasicMaterial({ color: 0xb44838, transparent: true, opacity: 0.85 }),
        errorNegative: new T.LineBasicMaterial({ color: 0x2d6cdf, transparent: true, opacity: 0.85 }),
        surfaceNormal: new T.LineBasicMaterial({ color: 0x2e7d55, transparent: true, opacity: 0.78 }),
        surfaceNormalHead: new T.MeshBasicMaterial({ color: 0x2e7d55, transparent: true, opacity: 0.82 }),
        footprintAlpha: new T.LineBasicMaterial({ color: 0x2d6cdf, transparent: true, opacity: 0.72 }),
        footprintZ: new T.LineBasicMaterial({ color: 0xb44838, transparent: true, opacity: 0.76 }),
        footprintSource: new T.LineBasicMaterial({ color: 0x202936, transparent: true, opacity: 0.9 }),
        paintBrushPreview: new T.LineBasicMaterial({ color: 0x0f6b8f, transparent: true, opacity: 0.9 }),
        partCallout: new T.LineBasicMaterial({ color: 0x1f2a36, transparent: true, opacity: 0.72 }),
        scaleMarker: new T.LineBasicMaterial({ color: 0x202936, transparent: true, opacity: 0.58 }),
      };
    }

    createGeometries() {
      const T = this.THREE;
      const hinge = new T.CylinderGeometry(0.035, 0.035, 0.12, 16);
      hinge.rotateX(Math.PI / 2);
      const screwHead = new T.CylinderGeometry(0.026, 0.026, 0.012, 18);
      screwHead.rotateX(Math.PI / 2);
      return {
        plate: this.makeBeveledPlateGeometry(),
        plateEdge: new T.EdgesGeometry(this.makeBeveledPlateGeometry(), 24),
        abstractBody: new T.BoxGeometry(1, 1, 0.045),
        abstractInner: new T.BoxGeometry(1, 1, 0.065),
        abstractCorner: new T.CylinderGeometry(0.045, 0.045, 0.05, 16),
        abstractGuide: new T.BoxGeometry(1, 0.018, 0.018),
        abstractDatum: new T.BoxGeometry(1, 0.01, 0.01),
        paperPin: new T.CylinderGeometry(1, 1, 0.055, 20),
        paperClearanceRing: new T.TorusGeometry(1, 0.045, 8, 30),
        screwHead,
        screwSlot: new T.BoxGeometry(0.04, 0.006, 0.006),
        hinge,
        link: new T.BoxGeometry(0.52, 0.035, 0.035),
        interCellLink: new T.BoxGeometry(1, 0.035, 0.045),
        pivotBoss: new T.CylinderGeometry(0.058, 0.058, 0.08, 18),
        diagonalBrace: new T.BoxGeometry(1, 0.026, 0.034),
        actuatorRail: new T.CylinderGeometry(0.028, 0.028, 1, 20),
        actuatorSleeve: new T.CylinderGeometry(0.044, 0.044, 0.32, 20),
        actuatorTip: new T.SphereGeometry(0.065, 16, 12),
        actuatorStop: new T.TorusGeometry(0.085, 0.01, 8, 24),
        actuatorTravelGauge: new T.BoxGeometry(0.055, 0.018, 1),
        actuatorTravelNeedle: new T.BoxGeometry(0.14, 0.028, 0.028),
        actuatorGaugeTick: new T.BoxGeometry(0.18, 0.018, 0.018),
        alphaActuatorRail: new T.BoxGeometry(1, 0.032, 0.032),
        alphaActuatorSleeve: new T.BoxGeometry(0.24, 0.06, 0.05),
        alphaActuatorStop: new T.BoxGeometry(0.032, 0.13, 0.08),
        backlashStop: new T.BoxGeometry(0.11, 0.028, 0.075),
        displacementHead: new T.ConeGeometry(0.045, 0.11, 12),
      };
    }

    buildScene() {
      const T = this.THREE;
      const hemi = new T.HemisphereLight(0xffffff, 0x6a7280, 1.2);
      const key = new T.DirectionalLight(0xffffff, 1.4);
      key.position.set(6, -8, 10);
      key.castShadow = true;
      key.shadow.mapSize.set(2048, 2048);
      this.scene.add(hemi, key);
      const grid = new T.GridHelper(18, 18, 0x6f7b89, 0xd0d6df);
      grid.rotation.x = Math.PI / 2;
      grid.position.z = -0.55;
      this.scene.add(grid);
      const axes = new T.AxesHelper(2.2);
      axes.position.set(-7.8, -7.8, -0.48);
      this.scene.add(axes);
      this.addAxisLabel("X", new T.Vector3(-5.35, -7.8, -0.48), 0xb44838);
      this.addAxisLabel("Y", new T.Vector3(-7.8, -5.35, -0.48), 0x2e7d55);
      this.addAxisLabel("Z", new T.Vector3(-7.8, -7.8, 1.95), 0x2d6cdf);
      this.addScaleMarkers();
    }

    addAxisLabel(text, position, color) {
      const T = this.THREE;
      const canvas = document.createElement("canvas");
      canvas.width = 96;
      canvas.height = 96;
      const ctx = canvas.getContext("2d");
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.fillStyle = `#${color.toString(16).padStart(6, "0")}`;
      ctx.font = "700 54px Inter, Arial, sans-serif";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(text, 48, 48);
      const texture = new T.CanvasTexture(canvas);
      const material = new T.SpriteMaterial({ map: texture, transparent: true });
      const sprite = new T.Sprite(material);
      sprite.position.copy(position);
      sprite.scale.set(0.42, 0.42, 0.42);
      sprite.userData.disposeMaterial = true;
      this.scene.add(sprite);
      return sprite;
    }

    addScaleMarkers() {
      const T = this.THREE;
      const root = new T.Group();
      root.position.set(5.15, -7.85, -0.48);
      const z = 0.02;
      const points = [new T.Vector3(0, 0, z), new T.Vector3(3, 0, z)];
      const base = new T.Line(new T.BufferGeometry().setFromPoints(points), this.materials.scaleMarker);
      base.userData.disposeGeometry = true;
      root.add(base);
      for (let i = 0; i <= 3; i += 1) {
        const tickHeight = i === 0 || i === 3 ? 0.24 : 0.16;
        const tick = new T.Line(
          new T.BufferGeometry().setFromPoints([
            new T.Vector3(i, -tickHeight / 2, z),
            new T.Vector3(i, tickHeight / 2, z),
          ]),
          this.materials.scaleMarker
        );
        tick.userData.disposeGeometry = true;
        root.add(tick);
        const label = this.makeSceneLabel(`${i}`, new T.Vector3(i, -0.34, z + 0.02), 0x202936, 0.22);
        root.add(label);
      }
      root.add(this.makeSceneLabel("normalized cell units", new T.Vector3(1.5, -0.62, z + 0.02), 0x4d5866, 0.28));
      this.scene.add(root);
    }

    makeSceneLabel(text, position, color, scale = 0.42) {
      const T = this.THREE;
      const canvas = document.createElement("canvas");
      canvas.width = 256;
      canvas.height = 96;
      const ctx = canvas.getContext("2d");
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.fillStyle = `#${color.toString(16).padStart(6, "0")}`;
      ctx.font = "700 34px Inter, Arial, sans-serif";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(text, canvas.width / 2, canvas.height / 2);
      const texture = new T.CanvasTexture(canvas);
      const material = new T.SpriteMaterial({ map: texture, transparent: true });
      const sprite = new T.Sprite(material);
      sprite.position.copy(position);
      sprite.scale.set(scale * (canvas.width / canvas.height), scale, scale);
      sprite.userData.disposeMaterial = true;
      return sprite;
    }

    bind() {
      window.addEventListener("resize", () => this.resize());
      const canvas = this.renderer.domElement;
      canvas.addEventListener("contextmenu", (event) => event.preventDefault());
      canvas.addEventListener("pointerdown", (event) => {
        const mode = event.shiftKey || event.button === 1 || event.button === 2 ? "pan" : "orbit";
        this.drag = {
          x: event.clientX,
          y: event.clientY,
          moved: false,
          mode,
          theta: this.spherical.theta,
          phi: this.spherical.phi,
          target: this.cameraTarget.clone(),
        };
        canvas.setPointerCapture(event.pointerId);
      });
      canvas.addEventListener("pointermove", (event) => {
        if (!this.drag) {
          this.hoverPick(event);
          return;
        }
        const dx = event.clientX - this.drag.x;
        const dy = event.clientY - this.drag.y;
        if (Math.abs(dx) + Math.abs(dy) > 3) this.drag.moved = true;
        if (this.drag.mode === "pan") {
          this.panCamera(dx, dy);
        } else {
          this.spherical.theta = this.drag.theta - dx * 0.006;
          this.spherical.phi = Math.max(0.18, Math.min(Math.PI - 0.18, this.drag.phi + dy * 0.006));
        }
        this.updateCamera();
        if (this.drag.moved) this.setCustomView(true);
      });
      canvas.addEventListener("pointerup", (event) => {
        if (this.drag && this.drag.mode === "orbit" && !this.drag.moved) this.pick(event);
        this.drag = null;
        this.hoverPick(event);
      });
      canvas.addEventListener("pointerleave", () => this.setHoveredCell(null));
      canvas.addEventListener("wheel", (event) => {
        event.preventDefault();
        this.spherical.radius = Math.max(4, Math.min(34, this.spherical.radius + event.deltaY * 0.012));
        this.updateCamera();
        this.setCustomView(true);
      });
    }

    panCamera(dx, dy) {
      const T = this.THREE;
      const cameraDirection = new T.Vector3();
      this.camera.getWorldDirection(cameraDirection);
      const right = new T.Vector3().crossVectors(cameraDirection, this.camera.up).normalize();
      const up = new T.Vector3().crossVectors(right, cameraDirection).normalize();
      const scale = this.spherical.radius * 0.0018;
      this.cameraTarget.copy(this.drag.target).addScaledVector(right, -dx * scale).addScaledVector(up, dy * scale);
    }

    resize() {
      const rect = this.mount.getBoundingClientRect();
      const width = Math.max(320, rect.width);
      const height = Math.max(320, rect.height);
      this.renderer.setSize(width, height, false);
      this.viewportAspect = width / height;
      this.updateProjectionMatrices();
    }

    updateCamera() {
      const { radius, theta, phi } = this.spherical;
      this.applyCameraUpForCurrentView();
      const position = new this.THREE.Vector3(
        radius * Math.sin(phi) * Math.cos(theta),
        radius * Math.sin(phi) * Math.sin(theta),
        radius * Math.cos(phi)
      );
      for (const camera of [this.perspectiveCamera, this.orthographicCamera]) {
        camera.position.copy(position);
        camera.lookAt(this.cameraTarget);
      }
      this.camera = this.projectionMode === "orthographic" ? this.orthographicCamera : this.perspectiveCamera;
      this.updateProjectionMatrices();
    }

    updateProjectionMatrices() {
      const aspect = Math.max(0.1, this.viewportAspect || 1);
      this.perspectiveCamera.aspect = aspect;
      this.perspectiveCamera.updateProjectionMatrix();
      const viewSize = Math.max(3.2, this.spherical.radius * 0.82);
      this.orthographicCamera.left = (-viewSize * aspect) / 2;
      this.orthographicCamera.right = (viewSize * aspect) / 2;
      this.orthographicCamera.top = viewSize / 2;
      this.orthographicCamera.bottom = -viewSize / 2;
      this.orthographicCamera.updateProjectionMatrix();
    }

    applyCameraUpForCurrentView() {
      if (this.viewMode === "top") {
        this.perspectiveCamera.up.set(0, 1, 0);
        this.orthographicCamera.up.set(0, 1, 0);
      } else {
        this.perspectiveCamera.up.set(0, 0, 1);
        this.orthographicCamera.up.set(0, 0, 1);
      }
    }

    setProjectionMode(mode) {
      this.projectionMode = mode === "orthographic" ? "orthographic" : "perspective";
      this.camera = this.projectionMode === "orthographic" ? this.orthographicCamera : this.perspectiveCamera;
      this.updateCamera();
      return this.projectionMode;
    }

    setView(mode) {
      const normalized = ["iso", "top", "front", "side"].includes(mode) ? mode : "iso";
      this.viewMode = normalized;
      if (normalized === "top") this.setCameraPosition(0, -0.001, 13);
      else if (normalized === "front") this.setCameraPosition(0, -13, 0.001);
      else if (normalized === "side") this.setCameraPosition(13, 0, 0.001);
      else this.setCameraPosition(9, -9, 9);
      this.updateCamera();
      this.notifyViewChange();
    }

    resetView() {
      this.cameraTarget.set(0, 0, 0);
      this.setView("iso");
    }

    setCustomView(forceNotify = false) {
      if (this.viewMode === "custom") {
        if (forceNotify) this.notifyViewChange();
        return;
      }
      this.viewMode = "custom";
      this.notifyViewChange();
    }

    notifyViewChange() {
      if (typeof this.onViewChange === "function") this.onViewChange(this.viewMode);
    }

    getCameraState() {
      return {
        mode: this.viewMode || "custom",
        projection: this.projectionMode,
        radius: this.spherical.radius,
        theta: this.spherical.theta,
        phi: this.spherical.phi,
        target: {
          x: this.cameraTarget.x,
          y: this.cameraTarget.y,
          z: this.cameraTarget.z,
        },
      };
    }

    getOrientationState() {
      const T = this.THREE;
      const fromTarget = new T.Vector3().copy(this.camera.position).sub(this.cameraTarget);
      if (fromTarget.lengthSq() < 1e-9) fromTarget.set(0, 0, 1);
      fromTarget.normalize();
      const up = new T.Vector3().copy(this.camera.up).normalize();
      return {
        viewAxis: this.formatDominantAxis(fromTarget),
        upAxis: this.formatDominantAxis(up),
        projection: this.projectionMode === "orthographic" ? "ortho" : "persp",
        mode: this.viewMode || "custom",
      };
    }

    formatDominantAxis(vector) {
      const entries = [
        ["X", vector.x],
        ["Y", vector.y],
        ["Z", vector.z],
      ];
      let dominant = entries[0];
      for (const entry of entries.slice(1)) {
        if (Math.abs(entry[1]) > Math.abs(dominant[1])) dominant = entry;
      }
      return `${dominant[1] >= 0 ? "+" : "-"}${dominant[0]}`;
    }

    applyCameraState(cameraState) {
      if (!cameraState) return;
      this.setProjectionMode(cameraState.projection || this.projectionMode || "perspective");
      if (typeof cameraState.radius === "number" && typeof cameraState.theta === "number" && typeof cameraState.phi === "number") {
        this.viewMode = cameraState.mode || "custom";
        this.spherical = {
          radius: Math.max(4, Math.min(34, cameraState.radius)),
          theta: cameraState.theta,
          phi: Math.max(0.02, Math.min(Math.PI - 0.02, cameraState.phi)),
        };
        if (cameraState.target) {
          this.cameraTarget.set(
            Number(cameraState.target.x) || 0,
            Number(cameraState.target.y) || 0,
            Number(cameraState.target.z) || 0
          );
        }
        this.updateCamera();
        this.notifyViewChange();
      } else {
        this.setView(cameraState.mode || "iso");
      }
    }

    frameCell(state, sim) {
      const { r, c } = state.selection || { r: 0, c: 0 };
      const center = sim?.centers?.[r]?.[c];
      if (!center) return;
      this.cameraTarget.set(center.x, center.y, center.z);
      this.spherical.radius = Math.max(4, Math.min(9, Math.max(state.grid.rows, state.grid.cols) * 0.82));
      this.viewMode = "custom";
      this.updateCamera();
      this.notifyViewChange();
    }

    setCameraPosition(x, y, z) {
      const radius = Math.max(0.001, Math.hypot(x, y, z));
      this.spherical = {
        radius,
        theta: Math.atan2(y, x),
        phi: Math.acos(Math.max(-1, Math.min(1, z / radius))),
      };
    }

    clearGroup(group) {
      while (group.children.length) {
        const child = group.children.pop();
        child.traverse((obj) => {
          if (obj.geometry && obj.userData.disposeGeometry) obj.geometry.dispose();
          if (obj.material?.map && obj.userData.disposeMaterial) obj.material.map.dispose();
          if (obj.material && obj.userData.disposeMaterial) obj.material.dispose();
        });
      }
    }

    overlayColor(state, sim, r, c) {
      const mode = state.view.overlayMode;
      if (state.selection.r === r && state.selection.c === c) return this.materials.selected;
      if (this.hoveredCell?.r === r && this.hoveredCell?.c === c) return this.materials.hovered;
      if (state.cells.locked[r][c]) return this.materials.locked;
      if (mode === "state" && (Math.abs(state.cells.commandAlpha[r][c]) > 1e-9 || Math.abs(state.cells.commandZ[r][c]) > 1e-9)) {
        return this.materials.actuated;
      }
      if (mode !== "state") {
        let t = 0.5;
        if (mode === "height") t = (sim.height[r][c] + 0.8) / 1.6;
        else if (mode === "zresidual") t = ((sim.zResidual?.[r]?.[c] || 0) + 0.8) / 1.6;
        else if (mode === "influence") t = (sim.influence[r][c] + 0.75) / 1.5;
        else if (mode === "theta") t = (sim.theta[r][c] + 45) / 105;
        else if (mode === "error") t = (sim.targetError[r][c] + 0.8) / 1.6;
        else if (mode === "travel") t = (Math.abs(state.cells.commandAlpha[r][c]) + Math.abs(state.cells.commandZ[r][c])) / 1.55;
        else if (mode === "saturation") t = sim.saturation?.[r]?.[c] || RAD.commandSaturation(state, state.cells.commandAlpha[r][c], state.cells.commandZ[r][c]);
        else if (mode === "strain") t = this.localLinkStrainStrength(sim, r, c);
        else if (mode === "displacement") t = (sim.displacement?.[r]?.[c] || 0) / Math.max(1e-9, sim.metrics?.maxReferenceDisplacement || 1);
        else if (mode === "slope") t = (sim.slope?.magnitude?.[r]?.[c] || 0) / Math.max(1e-9, sim.metrics?.maxSurfaceSlope || 1);
        else if (mode === "inverse") t = this.inversePlanStrength(state, r, c);
        else if (mode === "sensitivity") t = this.sensitivityStrength(state, r, c);
        else if (mode === "reachability") t = this.reachabilityStrength(state, r, c);
        else t = (sim.alpha[r][c] - state.grid.alphaMin) / (state.grid.alphaMax - state.grid.alphaMin);
        return this.overlayMaterial(mode, t);
      }
      return this.materials.plate;
    }

    sensitivityStrength(state, r, c) {
      const sensitivity = state.inverse?.sensitivity;
      const value = sensitivity?.map?.[r]?.[c] || 0;
      const scale = Math.max(1e-9, sensitivity?.maxGain || 0);
      return value / scale;
    }

    reachabilityStrength(state, r, c) {
      const jacobian = state.inverse?.jacobian;
      const value = jacobian?.coverageMap?.[r]?.[c] || 0;
      const scale = Math.max(1e-9, jacobian?.maxCoverage || 0);
      return value / scale;
    }

    localLinkStrainStrength(sim, r, c) {
      const values = [];
      if (sim.linkStrain?.horizontal?.[r]?.[c] !== undefined) values.push(Math.abs(sim.linkStrain.horizontal[r][c]));
      if (sim.linkStrain?.horizontal?.[r]?.[c - 1] !== undefined) values.push(Math.abs(sim.linkStrain.horizontal[r][c - 1]));
      if (sim.linkStrain?.vertical?.[r]?.[c] !== undefined) values.push(Math.abs(sim.linkStrain.vertical[r][c]));
      if (sim.linkStrain?.vertical?.[r - 1]?.[c] !== undefined) values.push(Math.abs(sim.linkStrain.vertical[r - 1][c]));
      const local = values.length ? Math.max(...values) : 0;
      const scale = Math.max(0.025, sim.metrics?.maxAbsLinkStrain || 0);
      return local / scale;
    }

    inversePlanStrength(state, r, c) {
      const preview = state.inverse?.preview;
      if (preview?.contribution?.[r]?.[c] !== undefined) return this.inversePreviewContributionStrength(preview, r, c);
      const plan = state.inverse?.plan || {};
      const linear = state.inverse?.linearSolution || {};
      const candidates = [...(plan.commands || []), ...(linear.commands || []), ...(plan.candidates || [])];
      const best = candidates.reduce((max, candidate) => Math.max(max, Math.max(0, candidate.improvement || 0)), 0);
      if (best <= 0) return 0;
      const found = candidates.find((candidate) => candidate.r === r && candidate.c === c);
      return found ? Math.max(0, found.improvement || 0) / best : 0;
    }

    inversePreviewContributionStrength(preview, r, c) {
      const scale = Math.max(1e-9, preview.maxContribution || 0);
      const signed = preview.contribution[r][c] / scale;
      return 0.5 + 0.5 * Math.max(-1, Math.min(1, signed));
    }

    overlayMaterial(mode, value) {
      const T = this.THREE;
      const bucket = Math.max(0, Math.min(48, Math.round(value * 48)));
      const key = `${mode}:${bucket}`;
      if (!this.overlayMaterials.has(key)) {
        const material = this.materials.plate.clone();
        const t = bucket / 48;
        material.color = new T.Color().setHSL(0.58 - 0.46 * t, 0.62, 0.48);
        this.overlayMaterials.set(key, material);
      }
      return this.overlayMaterials.get(key);
    }

    isSelectedCell(state, r, c) {
      return state.selection?.r === r && state.selection?.c === c;
    }

    isCellVisible(state, r, c) {
      return state.view.isolateSelected !== true || this.isSelectedCell(state, r, c);
    }

    explodedCellAmount(state, r, c) {
      return state.view.explodedSelected === true && this.isSelectedCell(state, r, c) && (state.view.cellVisualMode || "abstract") !== "abstract" ? 0.22 : 0;
    }

    explodedPlatePosition(position, amount, z = 0) {
      if (amount <= 0) return [position[0], position[1], z];
      const length = Math.max(1e-6, Math.hypot(position[0], position[1]));
      return [
        position[0] + (position[0] / length) * amount,
        position[1] + (position[1] / length) * amount,
        z + amount * 0.55,
      ];
    }

    renderState(state, sim) {
      this.targetState = state;
      this.targetSim = sim;
      if (
        !this.displaySim ||
        this.displaySim.alpha.length !== sim.alpha.length ||
        this.displaySim.alpha[0].length !== sim.alpha[0].length
      ) {
        this.displaySim = this.cloneSim(sim);
      }
      this.ensureCellObjects(state);
      this.updateDynamicObjects(this.displaySim);
    }

    cloneSim(sim) {
      return JSON.parse(JSON.stringify(sim));
    }

    stepDisplaySim() {
      if (!this.targetSim || !this.displaySim) return false;
      let maxDelta = 0;
      const lerp = Math.max(0.04, Math.min(0.55, Number(this.targetState?.view?.animationResponse || 0.18)));
      const keys = ["alpha", "theta", "height", "zResidual", "influence", "target", "saturation", "slope"];
      for (const key of keys) {
        if (!this.targetSim[key]) continue;
        if (!this.displaySim[key]) this.displaySim[key] = this.cloneSim(this.targetSim[key]);
        for (let r = 0; r < this.targetSim[key].length; r += 1) {
          for (let c = 0; c < this.targetSim[key][r].length; c += 1) {
            const before = this.displaySim[key][r][c];
            const after = before + (this.targetSim[key][r][c] - before) * lerp;
            this.displaySim[key][r][c] = after;
            maxDelta = Math.max(maxDelta, Math.abs(after - this.targetSim[key][r][c]));
          }
        }
      }
      for (let r = 0; r < this.targetSim.centers.length; r += 1) {
        for (let c = 0; c < this.targetSim.centers[r].length; c += 1) {
          for (const axis of ["x", "y", "z"]) {
            const before = this.displaySim.centers[r][c][axis];
            const after = before + (this.targetSim.centers[r][c][axis] - before) * lerp;
            this.displaySim.centers[r][c][axis] = after;
            maxDelta = Math.max(maxDelta, Math.abs(after - this.targetSim.centers[r][c][axis]));
          }
        }
      }
      this.displaySim.metrics = this.targetSim.metrics;
      return maxDelta > 0.002;
    }

    makeBeveledPlateGeometry() {
      const T = this.THREE;
      const shape = new T.Shape();
      const half = 0.15;
      const radius = 0.035;
      shape.moveTo(-half + radius, -half);
      shape.lineTo(half - radius, -half);
      shape.quadraticCurveTo(half, -half, half, -half + radius);
      shape.lineTo(half, half - radius);
      shape.quadraticCurveTo(half, half, half - radius, half);
      shape.lineTo(-half + radius, half);
      shape.quadraticCurveTo(-half, half, -half, half - radius);
      shape.lineTo(-half, -half + radius);
      shape.quadraticCurveTo(-half, -half, -half + radius, -half);
      const geometry = new T.ExtrudeGeometry(shape, {
        depth: 0.075,
        bevelEnabled: true,
        bevelSegments: 2,
        bevelSize: 0.012,
        bevelThickness: 0.012,
      });
      geometry.center();
      return geometry;
    }

    createPlateHardware() {
      const T = this.THREE;
      const parts = [];
      const edge = new T.LineSegments(this.geometries.plateEdge, this.materials.plateEdge);
      edge.position.z = 0.002;
      parts.push(edge);
      for (const x of [-0.075, 0.075]) {
        const screw = new T.Group();
        const head = new T.Mesh(this.geometries.screwHead, this.materials.fastener);
        const slot = new T.Mesh(this.geometries.screwSlot, this.materials.hinge);
        head.castShadow = true;
        slot.castShadow = true;
        screw.position.set(x, 0, 0.052);
        slot.position.z = 0.008;
        screw.add(head, slot);
        parts.push(screw);
      }
      return parts;
    }

    rebuild(sim) {
      this.ensureCellObjects(this.targetState, true);
      this.updateDynamicObjects(sim);
    }

    ensureCellObjects(state, force = false) {
      const shapeKey = `${state.grid.rows}x${state.grid.cols}:${state.view.surfaceSubdivisions || 5}:${state.view.surfaceInterpolation || "smooth"}`;
      if (!force && this.builtShape === shapeKey) return;
      this.clearGroup(this.cellRoot);
      this.clearGroup(this.linkageRoot);
      this.clearGroup(this.gapRoot);
      this.clearGroup(this.actuatorRoot);
      this.selectables = [];
      const { rows, cols, cellSize } = state.grid;
      const T = this.THREE;
      const plateGeometry = this.geometries.plate;
      const hingeGeometry = this.geometries.hinge;
      const linkGeometry = this.geometries.link;
      const actuatorGeometry = this.geometries.actuatorRail;
      const actuatorSleeveGeometry = this.geometries.actuatorSleeve;

      this.cellGroups.clear();
      this.latticeConnectors = [];
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          const group = new T.Group();
          group.userData = { r, c };
          const record = { group, plates: [], plateHardware: [], hinges: [], links: [], pivotBosses: [], braces: [], abstract: null, paperRad: null, gap: null, stops: [], actuator: null };
          for (let i = 0; i < 4; i += 1) {
            const plate = new T.Mesh(plateGeometry, this.materials.plate);
            plate.castShadow = true;
            plate.receiveShadow = true;
            plate.userData = { r, c, selectable: true };
            const hardware = this.createPlateHardware();
            for (const part of hardware) plate.add(part);
            group.add(plate);
            record.plates.push(plate);
            record.plateHardware.push(...hardware);
            this.selectables.push(plate);
            const hinge = new T.Mesh(hingeGeometry, this.materials.hinge);
            hinge.castShadow = true;
            group.add(hinge);
            record.hinges.push(hinge);
            const pivotBoss = new T.Mesh(this.geometries.pivotBoss, this.materials.brace);
            pivotBoss.rotation.x = Math.PI / 2;
            pivotBoss.castShadow = true;
            group.add(pivotBoss);
            record.pivotBosses.push(pivotBoss);
          }
          for (let i = 0; i < 4; i += 1) {
            const link = new T.Mesh(linkGeometry, this.materials.hinge);
            group.add(link);
            record.links.push(link);
          }
          for (let i = 0; i < 2; i += 1) {
            const brace = new T.Mesh(this.geometries.diagonalBrace, this.materials.brace);
            brace.castShadow = true;
            brace.receiveShadow = true;
            group.add(brace);
            record.braces.push(brace);
          }
          const abstract = new T.Group();
          const outer = new T.Mesh(this.geometries.abstractBody, this.materials.abstractBody);
          const inner = new T.Mesh(this.geometries.abstractInner, this.materials.abstractInner);
          const outerEdge = new T.LineSegments(new T.EdgesGeometry(this.geometries.abstractBody, 15), this.materials.abstractEdge);
          const innerEdge = new T.LineSegments(new T.EdgesGeometry(this.geometries.abstractInner, 15), this.materials.abstractEdge);
          const guideX = new T.Mesh(this.geometries.abstractGuide, this.materials.abstractGuide);
          const guideY = new T.Mesh(this.geometries.abstractGuide, this.materials.abstractGuide);
          const datumX = new T.Mesh(this.geometries.abstractDatum, this.materials.abstractDatum);
          const datumY = new T.Mesh(this.geometries.abstractDatum, this.materials.abstractDatum);
          const corners = [];
          outer.userData = { r, c, selectable: true };
          outer.castShadow = true;
          outer.receiveShadow = true;
          inner.castShadow = true;
          inner.receiveShadow = true;
          inner.position.z = 0.06;
          outerEdge.position.z = 0.03;
          innerEdge.position.z = 0.095;
          guideX.position.z = 0.14;
          guideY.position.z = 0.14;
          datumX.position.z = 0.12;
          datumY.position.z = 0.12;
          datumY.rotation.z = Math.PI / 2;
          abstract.add(outer, inner, outerEdge, innerEdge, datumX, datumY, guideX, guideY);
          this.selectables.push(outer);
          for (let i = 0; i < 4; i += 1) {
            const corner = new T.Mesh(this.geometries.abstractCorner, this.materials.hinge);
            corner.rotation.x = Math.PI / 2;
            corner.castShadow = true;
            abstract.add(corner);
            corners.push(corner);
          }
          group.add(abstract);
          record.abstract = { group: abstract, outer, inner, outerEdge, innerEdge, guideX, guideY, datumX, datumY, corners };

          const paperRad = new T.Group();
          const paperOuter = new T.Mesh(this.geometries.abstractBody, this.materials.paperOuter);
          const paperInner = new T.Mesh(this.geometries.abstractInner, this.materials.paperInner);
          const paperOuterEdge = new T.LineSegments(new T.EdgesGeometry(this.geometries.abstractBody, 15), this.materials.abstractEdge);
          const paperInnerEdge = new T.LineSegments(new T.EdgesGeometry(this.geometries.abstractInner, 15), this.materials.abstractEdge);
          const paperDatumX = new T.Mesh(this.geometries.abstractDatum, this.materials.abstractDatum);
          const paperDatumY = new T.Mesh(this.geometries.abstractDatum, this.materials.abstractDatum);
          const paperPins = [];
          const paperClearanceRings = [];
          paperOuter.userData = { r, c, selectable: true };
          paperOuter.castShadow = true;
          paperOuter.receiveShadow = true;
          paperInner.castShadow = true;
          paperInner.receiveShadow = true;
          paperInner.position.z = 0.075;
          paperOuterEdge.position.z = 0.035;
          paperInnerEdge.position.z = 0.115;
          paperDatumX.position.z = 0.135;
          paperDatumY.position.z = 0.135;
          paperDatumY.rotation.z = Math.PI / 2;
          paperRad.add(paperOuter, paperInner, paperOuterEdge, paperInnerEdge, paperDatumX, paperDatumY);
          this.selectables.push(paperOuter);
          for (let i = 0; i < 4; i += 1) {
            const pin = new T.Mesh(this.geometries.paperPin, this.materials.hinge);
            const clearanceRing = new T.Mesh(this.geometries.paperClearanceRing, this.materials.paperClearance);
            pin.rotation.x = Math.PI / 2;
            clearanceRing.position.z = 0.165;
            pin.castShadow = true;
            paperRad.add(clearanceRing, pin);
            paperPins.push(pin);
            paperClearanceRings.push(clearanceRing);
          }
          group.add(paperRad);
          record.paperRad = { group: paperRad, outer: paperOuter, inner: paperInner, outerEdge: paperOuterEdge, innerEdge: paperInnerEdge, datumX: paperDatumX, datumY: paperDatumY, pins: paperPins, clearanceRings: paperClearanceRings };
          this.cellRoot.add(group);

          const gapGeometry = new T.TorusGeometry(1, 0.006, 6, 48);
          const gap = new T.Mesh(gapGeometry, this.materials.gap);
          gap.userData.disposeGeometry = true;
          this.gapRoot.add(gap);
          record.gap = gap;
          for (let i = 0; i < 8; i += 1) {
            const stop = new T.Mesh(this.geometries.backlashStop, this.materials.stop);
            stop.castShadow = true;
            stop.receiveShadow = true;
            this.gapRoot.add(stop);
            record.stops.push(stop);
          }

          const actuator = new T.Group();
          const rail = new T.Mesh(actuatorGeometry, this.materials.hinge);
          const slider = new T.Mesh(actuatorSleeveGeometry, this.materials.actuator);
          const tip = new T.Mesh(this.geometries.actuatorTip, this.materials.actuator);
          const envelope = new T.Mesh(new T.CylinderGeometry(0.07, 0.07, 1, 20), this.materials.actuatorEnvelope);
          const baseStop = new T.Mesh(this.geometries.actuatorStop, this.materials.hinge);
          const topStop = new T.Mesh(this.geometries.actuatorStop, this.materials.actuator);
          const travelGauge = new T.Mesh(this.geometries.actuatorTravelGauge, this.materials.actuatorEnvelope);
          const travelNeedle = new T.Mesh(this.geometries.actuatorTravelNeedle, this.materials.actuator);
          const alphaRail = new T.Mesh(this.geometries.alphaActuatorRail, this.materials.hinge);
          const alphaSleeve = new T.Mesh(this.geometries.alphaActuatorSleeve, this.materials.actuator);
          const alphaNeedle = new T.Mesh(this.geometries.actuatorTravelNeedle, this.materials.actuator);
          const alphaMinStop = new T.Mesh(this.geometries.alphaActuatorStop, this.materials.hinge);
          const alphaMaxStop = new T.Mesh(this.geometries.alphaActuatorStop, this.materials.actuator);
          const zTravelTicks = [];
          const alphaTravelTicks = [];
          for (let i = 0; i < 5; i += 1) {
            const zTick = new T.Mesh(this.geometries.actuatorGaugeTick, this.materials.actuatorGaugeTick);
            const alphaTick = new T.Mesh(this.geometries.actuatorGaugeTick, this.materials.actuatorGaugeTick);
            zTravelTicks.push(zTick);
            alphaTravelTicks.push(alphaTick);
          }
          envelope.userData.disposeGeometry = true;
          envelope.rotation.x = Math.PI / 2;
          rail.rotation.x = Math.PI / 2;
          slider.rotation.x = Math.PI / 2;
          rail.castShadow = true;
          slider.castShadow = true;
          baseStop.castShadow = true;
          topStop.castShadow = true;
          travelNeedle.castShadow = true;
          alphaRail.castShadow = true;
          alphaSleeve.castShadow = true;
          alphaNeedle.castShadow = true;
          alphaMinStop.castShadow = true;
          alphaMaxStop.castShadow = true;
          actuator.add(envelope, rail, slider, tip, baseStop, topStop, travelGauge, travelNeedle, alphaRail, alphaSleeve, alphaNeedle, alphaMinStop, alphaMaxStop, ...zTravelTicks, ...alphaTravelTicks);
          this.actuatorRoot.add(actuator);
          record.actuator = { group: actuator, envelope, rail, slider, tip, baseStop, topStop, travelGauge, travelNeedle, alphaRail, alphaSleeve, alphaNeedle, alphaMinStop, alphaMaxStop, zTravelTicks, alphaTravelTicks };
          this.cellGroups.set(`${r},${c}`, record);
        }
      }
      this.createLatticeConnectors(state);
      this.syncReferenceLattice(state);
      this.builtShape = shapeKey;
    }

    syncReferenceLattice(state) {
      const shapeKey = `${state.grid.rows}x${state.grid.cols}:${state.grid.cellSize}`;
      if (this.referenceShape === shapeKey) return;
      this.clearGroup(this.referenceRoot);
      const T = this.THREE;
      const { rows, cols, cellSize } = state.grid;
      const points = [];
      const x0 = -((cols - 1) * cellSize) / 2;
      const y0 = -((rows - 1) * cellSize) / 2;
      const z = -0.075;
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c + 1 < cols; c += 1) {
          points.push(new T.Vector3(x0 + c * cellSize, y0 + r * cellSize, z));
          points.push(new T.Vector3(x0 + (c + 1) * cellSize, y0 + r * cellSize, z));
        }
      }
      for (let r = 0; r + 1 < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          points.push(new T.Vector3(x0 + c * cellSize, y0 + r * cellSize, z));
          points.push(new T.Vector3(x0 + c * cellSize, y0 + (r + 1) * cellSize, z));
        }
      }
      const geometry = new T.BufferGeometry().setFromPoints(points);
      geometry.userData.disposeGeometry = true;
      const lines = new T.LineSegments(geometry, this.materials.reference);
      lines.userData.disposeGeometry = true;
      this.referenceRoot.add(lines);
      this.referenceShape = shapeKey;
    }

    createLatticeConnectors(state) {
      const T = this.THREE;
      const { rows, cols } = state.grid;
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          if (c + 1 < cols) this.addLatticeConnector(r, c, r, c + 1, "x");
          if (r + 1 < rows) this.addLatticeConnector(r, c, r + 1, c, "y");
        }
      }
    }

    addLatticeConnector(r0, c0, r1, c1, axis) {
      const T = this.THREE;
      const group = new T.Group();
      const rod = new T.Mesh(this.geometries.interCellLink, this.materials.linkage);
      const pinA = new T.Mesh(this.geometries.hinge, this.materials.hinge);
      const pinB = new T.Mesh(this.geometries.hinge, this.materials.hinge);
      rod.castShadow = true;
      rod.receiveShadow = true;
      pinA.castShadow = true;
      pinB.castShadow = true;
      group.add(rod, pinA, pinB);
      this.linkageRoot.add(group);
      this.latticeConnectors.push({ r0, c0, r1, c1, axis, group, rod, pinA, pinB });
    }

    updateDynamicObjects(sim) {
      this.clearGroup(this.errorRoot);
      this.clearGroup(this.surfaceNormalRoot);
      this.clearGroup(this.measurementRoot);
      this.clearGroup(this.displacementVectorRoot);
      this.clearGroup(this.partCalloutRoot);
      this.clearGroup(this.influenceFootprintRoot);
      this.clearGroup(this.paintBrushPreviewRoot);
      const state = this.targetState;
      const { rows, cols, cellSize } = state.grid;
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          const record = this.cellGroups.get(`${r},${c}`);
          const center = sim.centers[r][c];
          const theta = (sim.theta[r][c] * Math.PI) / 180;
          const radius = 0.31 * Math.sqrt(sim.alpha[r][c]) * cellSize;
          const positions = [
            [-radius, -radius],
            [radius, -radius],
            [radius, radius],
            [-radius, radius],
          ];
          const cellVisible = this.isCellVisible(state, r, c);
          const explodeAmount = this.explodedCellAmount(state, r, c);
          const explodedPositions = positions.map((position) => this.explodedPlatePosition(position, explodeAmount));
          record.group.visible = cellVisible;
          record.group.position.set(center.x, center.y, center.z);
          const material = this.overlayColor(state, sim, r, c);
          const visualMode = state.view.cellVisualMode || "abstract";
          const abstractMode = visualMode === "abstract";
          const paperRadMode = visualMode === "paperRad";
          const mechanismMode = visualMode === "mechanism";
          record.abstract.group.visible = abstractMode;
          record.paperRad.group.visible = paperRadMode;
          for (let i = 0; i < 4; i += 1) {
            const [px, py] = positions[i];
            const [ex, ey, ez] = explodedPositions[i];
            const plate = record.plates[i];
            plate.visible = mechanismMode;
            plate.material = material;
            plate.position.set(ex, ey, ez);
            plate.rotation.z = theta;
            for (const part of record.plateHardware) part.visible = mechanismMode && state.view.fastenersVisible !== false;
            record.hinges[i].visible = mechanismMode;
            record.hinges[i].position.set(ex, ey, 0.08 + ez + explodeAmount * 0.35);
            record.pivotBosses[i].visible = mechanismMode && state.view.pivotsVisible !== false;
            record.pivotBosses[i].position.set(ex, ey, 0.13 + ez + explodeAmount * 0.45);
            record.pivotBosses[i].scale.setScalar(0.85 + Math.max(0, state.grid.backlash) * 0.65);
          }
          for (let i = 0; i < 4; i += 1) {
            const a = explodedPositions[i];
            const b = explodedPositions[(i + 1) % explodedPositions.length];
            const link = record.links[i];
            link.visible = mechanismMode;
            link.position.set((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2 - 0.02);
            link.scale.x = Math.hypot(a[0] - b[0], a[1] - b[1]);
            link.rotation.z = Math.atan2(b[1] - a[1], b[0] - a[0]);
          }
          this.updateDiagonalBraces(record, state, explodedPositions, explodeAmount);
          for (const brace of record.braces) brace.visible = mechanismMode && state.view.pivotsVisible !== false;
          const outerSize = Math.max(0.42, radius * 2 + state.grid.backlash * 0.45);
          const innerSize = Math.max(0.18, radius * 1.16);
          record.abstract.outer.material = material;
          record.abstract.outer.scale.set(outerSize, outerSize, 1);
          record.abstract.inner.scale.set(innerSize, innerSize, 1);
          record.abstract.inner.rotation.z = theta;
          record.abstract.outerEdge.scale.set(outerSize, outerSize, 1);
          record.abstract.innerEdge.scale.set(innerSize, innerSize, 1);
          record.abstract.innerEdge.rotation.z = theta;
          record.abstract.datumX.scale.set(outerSize * 0.86, 1, 1);
          record.abstract.datumY.scale.set(outerSize * 0.86, 1, 1);
          record.abstract.guideX.scale.set(innerSize * 0.9, 1, 1);
          record.abstract.guideY.scale.set(innerSize * 0.9, 1, 1);
          record.abstract.guideX.rotation.z = theta;
          record.abstract.guideY.rotation.z = theta + Math.PI / 2;
          for (let i = 0; i < record.abstract.corners.length; i += 1) {
            const [px, py] = positions[i];
            record.abstract.corners[i].position.set(px, py, 0.13);
            record.abstract.corners[i].scale.setScalar(0.9 + state.grid.backlash * 0.9);
          }
          const paperOuterSize = Math.max(0.44, radius * 2.05 + state.grid.backlash * 0.5);
          const paperInnerSize = Math.max(0.22, radius * 1.18);
          const pinVisualRadius = Math.max(0.028, Math.min(0.14, Number(state.grid.pinRadius ?? 0.18) * cellSize * 0.32));
          const holeVisualRadius = Math.max(pinVisualRadius + 0.006, Math.min(0.17, Number(state.grid.holeRadius ?? 0.225) * cellSize * 0.32));
          record.paperRad.outer.material = material;
          record.paperRad.outer.scale.set(paperOuterSize, paperOuterSize, 1);
          record.paperRad.inner.scale.set(paperInnerSize, paperInnerSize, 1);
          record.paperRad.inner.rotation.z = theta;
          record.paperRad.outerEdge.scale.set(paperOuterSize, paperOuterSize, 1);
          record.paperRad.innerEdge.scale.set(paperInnerSize, paperInnerSize, 1);
          record.paperRad.innerEdge.rotation.z = theta;
          record.paperRad.datumX.scale.set(paperOuterSize * 0.92, 1, 1);
          record.paperRad.datumY.scale.set(paperOuterSize * 0.92, 1, 1);
          for (let i = 0; i < record.paperRad.pins.length; i += 1) {
            const [px, py] = positions[i];
            record.paperRad.pins[i].position.set(px, py, 0.15);
            record.paperRad.pins[i].scale.set(pinVisualRadius, pinVisualRadius, 1);
            record.paperRad.clearanceRings[i].position.set(px, py, 0.17);
            record.paperRad.clearanceRings[i].scale.set(holeVisualRadius, holeVisualRadius, holeVisualRadius);
          }
          record.gap.visible = cellVisible && state.view.gapsVisible;
          record.gap.position.set(center.x, center.y, center.z + 0.11 + explodeAmount * 0.75);
          record.gap.scale.setScalar(radius + state.grid.backlash * 0.35 + explodeAmount * 0.15);
          this.updateBacklashStops(record, state, center, radius, theta, cellVisible, explodeAmount);

          const active = Math.abs(state.cells.commandAlpha[r][c]) > 1e-9 || Math.abs(state.cells.commandZ[r][c]) > 1e-9;
          const actuator = record.actuator;
          actuator.group.visible = this.shouldShowActuator(state, r, c, active);
          if (actuator.group.visible) {
            const limits = RAD.commandLimits(state);
            const travel = Math.max(0.42, 0.48 + limits.z * 1.05 + Math.max(limits.alphaContract, limits.alphaExpand) * 0.18);
            const commandZ = state.cells.commandZ[r][c];
            const saturation = RAD.commandSaturation(state, state.cells.commandAlpha[r][c], commandZ);
            actuator.group.position.set(center.x, center.y, center.z + explodeAmount * 0.5);
            actuator.envelope.position.set(0, 0, 0.12 + travel / 2);
            actuator.envelope.scale.y = travel;
            actuator.rail.position.set(0, 0, 0.12 + travel / 2);
            actuator.rail.scale.y = travel;
            const normalizedZ = Math.max(-1, Math.min(1, commandZ / Math.max(1e-9, limits.z)));
            const sliderZ = Math.max(0.18, Math.min(0.12 + travel - 0.1, 0.12 + travel / 2 + normalizedZ * travel * 0.42));
            actuator.slider.position.set(0, 0, sliderZ);
            actuator.slider.material = saturation >= 0.98 ? this.materials.actuated : commandZ >= 0 ? this.materials.actuatorPositive : this.materials.actuatorNegative;
            actuator.tip.position.set(0, 0, 0.16 + travel);
            actuator.tip.material = actuator.slider.material;
            actuator.baseStop.position.set(0, 0, 0.12);
            actuator.topStop.position.set(0, 0, 0.12 + travel);
            actuator.travelGauge.position.set(0.16, 0, 0.12 + travel / 2);
            actuator.travelGauge.scale.z = travel;
            actuator.travelNeedle.position.set(0.16, 0, sliderZ);
            actuator.travelNeedle.material = actuator.slider.material;
            actuator.travelNeedle.scale.x = 0.7 + Math.min(1.2, saturation) * 0.75;
            this.updateZTravelGaugeTicks(actuator, travel);

            const commandAlpha = state.cells.commandAlpha[r][c];
            const alphaLimit = commandAlpha < 0 ? limits.alphaContract : limits.alphaExpand;
            const normalizedAlpha = Math.max(-1, Math.min(1, commandAlpha / Math.max(1e-9, alphaLimit)));
            const alphaTravel = Math.max(0.5, 0.56 + Math.max(limits.alphaContract, limits.alphaExpand) * 0.32);
            const alphaX = normalizedAlpha * alphaTravel * 0.38;
            const alphaMaterial = saturation >= 0.98 ? this.materials.actuated : commandAlpha >= 0 ? this.materials.actuatorPositive : this.materials.actuatorNegative;
            actuator.alphaRail.position.set(0, -0.28, 0.22);
            actuator.alphaRail.scale.x = alphaTravel;
            actuator.alphaSleeve.position.set(alphaX, -0.28, 0.22);
            actuator.alphaSleeve.material = alphaMaterial;
            actuator.alphaNeedle.position.set(alphaX, -0.37, 0.22);
            actuator.alphaNeedle.rotation.z = Math.PI / 2;
            actuator.alphaNeedle.scale.x = 0.8 + Math.abs(normalizedAlpha) * 0.8;
            actuator.alphaNeedle.material = alphaMaterial;
            actuator.alphaMinStop.position.set(-alphaTravel * 0.42, -0.28, 0.22);
            actuator.alphaMaxStop.position.set(alphaTravel * 0.42, -0.28, 0.22);
            this.updateAlphaTravelGaugeTicks(actuator, alphaTravel);
          }
        }
      }

      this.syncSurfaceMeshes(state, sim);
      this.updateLatticeConnectors(state, sim);
      const isolated = state.view.isolateSelected === true;
      this.referenceRoot.visible = !isolated && state.view.referenceVisible !== false;
      if (!isolated) this.renderDisplacementVectors(state, sim);
      if (!isolated && (state.view.overlayMode === "error" || state.view.targetErrorVectorsVisible !== false) && state.target.type !== "none") this.renderErrorRods(state, sim);
      if (!isolated) this.renderSurfaceNormalVectors(state, sim);
      this.renderSelectedInfluenceFootprint(state, sim);
      this.renderPaintBrushPreview(state, sim);
      if (state.view.measurementsVisible) this.renderMeasurements(state, sim);
      this.renderExplodedPartCallouts(state, sim);
      this.fitRoot(rows, cols);
    }

    updateZTravelGaugeTicks(actuator, travel) {
      const fractions = [-1, -0.5, 0, 0.5, 1];
      for (let i = 0; i < actuator.zTravelTicks.length; i += 1) {
        const tick = actuator.zTravelTicks[i];
        const normalized = fractions[i] || 0;
        tick.position.set(0.16, 0, 0.12 + travel / 2 + normalized * travel * 0.42);
        tick.rotation.set(0, 0, 0);
        tick.scale.x = normalized === 0 ? 1.15 : 0.78;
      }
    }

    updateAlphaTravelGaugeTicks(actuator, alphaTravel) {
      const fractions = [-1, -0.5, 0, 0.5, 1];
      for (let i = 0; i < actuator.alphaTravelTicks.length; i += 1) {
        const tick = actuator.alphaTravelTicks[i];
        const normalized = fractions[i] || 0;
        tick.position.set(normalized * alphaTravel * 0.38, -0.45, 0.22);
        tick.rotation.z = Math.PI / 2;
        tick.scale.x = normalized === 0 ? 1.15 : 0.78;
      }
    }

    renderSelectedInfluenceFootprint(state, sim) {
      if (state.view.influenceFootprintVisible === false || typeof RAD.selectedCellFootprint !== "function") return;
      const T = this.THREE;
      const { rows, cols, cellSize } = state.grid;
      const { r: sr, c: sc } = state.selection || { r: 0, c: 0 };
      const source = sim.centers?.[sr]?.[sc];
      if (!source) return;
      const footprint = RAD.selectedCellFootprint(state, sr, sc);
      const alphaSource = Math.abs(state.cells.commandAlpha[sr][sc]);
      const zSource = Math.abs(state.cells.commandZ[sr][sc]);
      const maxAlpha = Math.max(1e-9, alphaSource);
      const maxZ = Math.max(1e-9, zSource);
      const sourcePoint = new T.Vector3(source.x, source.y, source.z + 0.34);
      this.addFootprintRing(sourcePoint, cellSize * 0.42, this.materials.footprintSource);
      if (alphaSource < 1e-9 && zSource < 1e-9) return;

      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          if (r === sr && c === sc) continue;
          const alphaValue = footprint.alpha[r][c] || 0;
          const zValue = footprint.zResidual[r][c] || 0;
          const alphaStrength = Math.abs(alphaValue) / maxAlpha;
          const zStrength = Math.abs(zValue) / maxZ;
          if (alphaStrength < 0.035 && zStrength < 0.035) continue;
          const center = sim.centers[r][c];
          const z = center.z + 0.28;
          if (alphaStrength >= 0.035) {
            const point = new T.Vector3(center.x, center.y, z + 0.02);
            const radius = cellSize * (0.19 + 0.2 * Math.min(1, alphaStrength));
            this.addFootprintRing(point, radius, this.materials.footprintAlpha);
            this.addFootprintLine(sourcePoint, point, this.materials.footprintAlpha, alphaStrength);
          }
          if (zStrength >= 0.035) {
            const point = new T.Vector3(center.x, center.y, z + 0.1);
            const radius = cellSize * (0.12 + 0.18 * Math.min(1, zStrength));
            this.addFootprintRing(point, radius, this.materials.footprintZ);
            this.addFootprintLine(sourcePoint, point, this.materials.footprintZ, zStrength);
          }
        }
      }
    }

    addFootprintRing(center, radius, material) {
      const T = this.THREE;
      const points = [];
      const segments = 48;
      for (let i = 0; i <= segments; i += 1) {
        const a = (i / segments) * Math.PI * 2;
        points.push(new T.Vector3(center.x + Math.cos(a) * radius, center.y + Math.sin(a) * radius, center.z));
      }
      const geometry = new T.BufferGeometry().setFromPoints(points);
      geometry.userData.disposeGeometry = true;
      const ring = new T.Line(geometry, material);
      ring.userData.disposeGeometry = true;
      this.influenceFootprintRoot.add(ring);
    }

    addFootprintLine(start, end, material, strength) {
      if (new this.THREE.Vector3().subVectors(end, start).lengthSq() < 1e-4) return;
      const T = this.THREE;
      const midpoint = new T.Vector3().lerpVectors(start, end, 0.5);
      midpoint.z += 0.08 + 0.16 * Math.min(1, strength);
      const geometry = new T.BufferGeometry().setFromPoints([start, midpoint, end]);
      geometry.userData.disposeGeometry = true;
      const line = new T.Line(geometry, material);
      line.userData.disposeGeometry = true;
      this.influenceFootprintRoot.add(line);
    }

    renderPaintBrushPreview(state, sim) {
      if (state.view.paintMode !== true || typeof RAD.brushCells !== "function") return;
      const hover = this.hoveredCell || state.selection;
      if (!hover) return;
      const cells = RAD.brushCells(state, hover.r, hover.c, state.view.paintRadius);
      const radius = Math.max(0.18, state.grid.cellSize * 0.36);
      for (const cell of cells) {
        const center = sim.centers?.[cell.r]?.[cell.c];
        if (!center) continue;
        const selected = cell.r === hover.r && cell.c === hover.c;
        this.addBrushPreviewRing(
          new this.THREE.Vector3(center.x, center.y, center.z + (selected ? 0.43 : 0.39)),
          radius * (selected ? 1.12 : 0.92)
        );
      }
    }

    addBrushPreviewRing(center, radius) {
      const T = this.THREE;
      const points = [];
      const segments = 40;
      for (let i = 0; i <= segments; i += 1) {
        const a = (i / segments) * Math.PI * 2;
        points.push(new T.Vector3(center.x + Math.cos(a) * radius, center.y + Math.sin(a) * radius, center.z));
      }
      const geometry = new T.BufferGeometry().setFromPoints(points);
      geometry.userData.disposeGeometry = true;
      const ring = new T.Line(geometry, this.materials.paintBrushPreview);
      ring.userData.disposeGeometry = true;
      this.paintBrushPreviewRoot.add(ring);
    }

    renderSurfaceNormalVectors(state, sim) {
      if (state.view.surfaceNormalsVisible !== true) return;
      const T = this.THREE;
      const { rows, cols } = state.grid;
      const step = Math.max(1, Math.ceil(Math.max(rows, cols) / 8));
      const maxTilt = Math.max(1e-9, sim.metrics?.maxNormalTilt || 1);
      for (let r = 0; r < rows; r += step) {
        for (let c = 0; c < cols; c += step) {
          const center = sim.centers[r][c];
          const normal = sim.slope?.normal?.[r]?.[c] || { x: 0, y: 0, z: 1 };
          const tilt = sim.slope?.tilt?.[r]?.[c] || 0;
          const length = 0.28 + 0.28 * Math.min(1, tilt / maxTilt);
          const start = new T.Vector3(center.x, center.y, center.z + 0.24);
          const end = new T.Vector3(center.x + normal.x * length, center.y + normal.y * length, center.z + 0.24 + normal.z * length);
          this.addSurfaceNormalVector(start, end);
        }
      }
    }

    addSurfaceNormalVector(start, end) {
      const T = this.THREE;
      const direction = new T.Vector3().subVectors(end, start);
      const length = direction.length();
      if (length < 1e-5) return;
      const headLength = Math.min(0.11, Math.max(0.055, length * 0.3));
      const shaftEnd = new T.Vector3().copy(end).addScaledVector(direction.clone().normalize(), -headLength);
      const lineGeometry = new T.BufferGeometry().setFromPoints([start, shaftEnd]);
      lineGeometry.userData.disposeGeometry = true;
      const line = new T.Line(lineGeometry, this.materials.surfaceNormal);
      line.userData.disposeGeometry = true;
      const head = new T.Mesh(this.geometries.displacementHead, this.materials.surfaceNormalHead);
      head.position.copy(end);
      head.quaternion.setFromUnitVectors(new T.Vector3(0, 1, 0), direction.clone().normalize());
      head.scale.setScalar(0.72);
      this.surfaceNormalRoot.add(line, head);
    }

    renderDisplacementVectors(state, sim) {
      if (state.view.displacementVectorsVisible === false) return;
      const T = this.THREE;
      const { rows, cols } = state.grid;
      const cells = [];
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) cells.push({ r, c, d: sim.displacement?.[r]?.[c] || 0 });
      }
      cells.sort((a, b) => b.d - a.d);
      const keep = new Map();
      for (const cell of cells.slice(0, Math.min(24, cells.length))) keep.set(`${cell.r},${cell.c}`, cell);
      keep.set(`${state.selection.r},${state.selection.c}`, { r: state.selection.r, c: state.selection.c, d: sim.displacement?.[state.selection.r]?.[state.selection.c] || 0 });
      const minVisible = Math.max(0.015, (sim.metrics.maxReferenceDisplacement || 0) * 0.04);
      for (const cell of keep.values()) {
        if (cell.d < minVisible && !(cell.r === state.selection.r && cell.c === state.selection.c)) continue;
        const ref = RAD.referenceCenter(state, cell.r, cell.c);
        const current = sim.centers[cell.r][cell.c];
        this.addDisplacementVector(new T.Vector3(ref.x, ref.y, ref.z + 0.03), new T.Vector3(current.x, current.y, current.z + 0.03), cell.d, sim.metrics.maxReferenceDisplacement || 1);
      }
    }

    addDisplacementVector(start, end, magnitude, maxMagnitude) {
      const T = this.THREE;
      const direction = new T.Vector3().subVectors(end, start);
      const length = direction.length();
      if (length < 1e-5) return;
      const headLength = Math.min(0.16, Math.max(0.07, length * 0.28));
      const shaftEnd = new T.Vector3().copy(end).addScaledVector(direction.clone().normalize(), -headLength);
      const lineGeometry = new T.BufferGeometry().setFromPoints([start, shaftEnd]);
      lineGeometry.userData.disposeGeometry = true;
      const line = new T.Line(lineGeometry, this.materials.displacementVector);
      line.userData.disposeGeometry = true;
      const head = new T.Mesh(this.geometries.displacementHead, this.materials.displacementHead);
      head.position.copy(end);
      head.quaternion.setFromUnitVectors(new T.Vector3(0, 1, 0), direction.clone().normalize());
      const scale = 0.75 + Math.min(1.5, magnitude / Math.max(1e-9, maxMagnitude)) * 0.85;
      head.scale.setScalar(scale);
      this.displacementVectorRoot.add(line, head);
    }

    updateLatticeConnectors(state, sim) {
      const visible = state.view.isolateSelected !== true && state.view.linkagesVisible !== false && (state.view.cellVisualMode || "abstract") !== "abstract";
      for (const connector of this.latticeConnectors) {
        connector.group.visible = visible;
        if (!visible) continue;
        const a = this.connectorEndpoint(state, sim, connector.r0, connector.c0, connector.axis, 1);
        const b = this.connectorEndpoint(state, sim, connector.r1, connector.c1, connector.axis, -1);
        this.placeConnector(connector, a, b);
        const strain = this.connectorStrain(sim, connector);
        connector.rod.material = state.view.overlayMode === "strain" ? this.overlayMaterial("strain", 0.5 + 0.5 * Math.min(1, Math.abs(strain) / Math.max(0.025, sim.metrics.maxAbsLinkStrain || 0.025))) : this.materials.linkage;
      }
    }

    connectorStrain(sim, connector) {
      if (connector.axis === "x") return sim.linkStrain?.horizontal?.[connector.r0]?.[connector.c0] || 0;
      return sim.linkStrain?.vertical?.[connector.r0]?.[connector.c0] || 0;
    }

    connectorEndpoint(state, sim, r, c, axis, side) {
      const center = sim.centers[r][c];
      const radius = 0.31 * Math.sqrt(sim.alpha[r][c]) * state.grid.cellSize;
      const z = center.z + 0.075;
      if (axis === "x") return { x: center.x + side * radius, y: center.y, z };
      return { x: center.x, y: center.y + side * radius, z };
    }

    placeConnector(connector, a, b) {
      const dx = b.x - a.x;
      const dy = b.y - a.y;
      const dz = b.z - a.z;
      const length = Math.max(0.04, Math.hypot(dx, dy, dz));
      connector.group.position.set((a.x + b.x) / 2, (a.y + b.y) / 2, (a.z + b.z) / 2);
      connector.group.rotation.z = Math.atan2(dy, dx);
      connector.group.rotation.y = -Math.atan2(dz, Math.hypot(dx, dy));
      connector.rod.position.set(0, 0, 0);
      connector.rod.scale.x = length;
      connector.pinA.position.set(-length / 2, 0, 0);
      connector.pinB.position.set(length / 2, 0, 0);
      const pinScale = 0.8 + Math.min(0.8, length) * 0.18;
      connector.pinA.scale.setScalar(pinScale);
      connector.pinB.scale.setScalar(pinScale);
    }

    shouldShowActuator(state, r, c, active) {
      if (!this.isCellVisible(state, r, c)) return false;
      if ((state.view.cellVisualMode || "abstract") === "abstract") return false;
      if (state.view.actuatorsVisible === false || !active) return false;
      const mode = state.view.actuatorDisplayMode || "selected";
      const selected = state.selection?.r === r && state.selection?.c === c;
      const planned = this.isPlannedActuator(state, r, c);
      if (mode === "active") return true;
      if (mode === "planned") return planned || selected;
      return selected || planned;
    }

    isPlannedActuator(state, r, c) {
      const plan = state.inverse?.plan || {};
      const items = [...(plan.commands || []), ...(plan.history || [])];
      return items.some((item) => item && item.r === r && item.c === c);
    }

    updateDiagonalBraces(record, state, positions, explodeAmount = 0) {
      const visible = state.view.pivotsVisible !== false;
      const pairs = [
        [positions[0], positions[2]],
        [positions[1], positions[3]],
      ];
      for (let i = 0; i < record.braces.length; i += 1) {
        const [a, b] = pairs[i];
        const brace = record.braces[i];
        brace.visible = visible;
        brace.position.set((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, 0.045 + explodeAmount * 0.28);
        brace.scale.x = Math.hypot(a[0] - b[0], a[1] - b[1]);
        brace.rotation.z = Math.atan2(b[1] - a[1], b[0] - a[0]);
      }
    }

    updateBacklashStops(record, state, center, radius, theta, cellVisible = true, explodeAmount = 0) {
      const visible = cellVisible && state.view.stopsVisible !== false && (state.view.cellVisualMode || "abstract") !== "abstract";
      const stopRadius = radius + 0.12 + state.grid.backlash * 0.42 + explodeAmount * 0.65;
      const clearance = 0.12 + state.grid.backlash * 0.7;
      for (let i = 0; i < record.stops.length; i += 1) {
        const quadrant = Math.floor(i / 2);
        const side = i % 2 === 0 ? -1 : 1;
        const angle = theta + Math.PI / 4 + quadrant * (Math.PI / 2) + side * clearance;
        const stop = record.stops[i];
        stop.visible = visible;
        stop.position.set(center.x + Math.cos(angle) * stopRadius, center.y + Math.sin(angle) * stopRadius, center.z + 0.165 + explodeAmount * 0.8);
        stop.rotation.z = angle + Math.PI / 2;
        stop.scale.x = 1 + state.grid.backlash * 1.8;
      }
    }

    syncSurfaceMeshes(state, sim) {
      const shapeKey = `${state.grid.rows}x${state.grid.cols}`;
      if (this.surfaceShape !== shapeKey) {
        this.clearGroup(this.membraneRoot);
        this.clearGroup(this.targetRoot);
        this.clearGroup(this.surfaceContourRoot);
        this.membraneMesh = null;
        this.targetMesh = null;
        this.surfaceShape = shapeKey;
      }
      if (!this.membraneMesh) {
        this.membraneMesh = this.createSurfaceMesh(state, sim, "membrane");
        this.membraneRoot.add(this.membraneMesh);
      }
      if (!this.targetMesh) {
        this.targetMesh = this.createSurfaceMesh(state, sim, "target");
        this.targetRoot.add(this.targetMesh);
      }
      const isolated = state.view.isolateSelected === true;
      this.membraneMesh.visible = !isolated && state.view.membraneVisible;
      this.targetMesh.visible = !isolated && state.view.targetVisible && state.target.type !== "none";
      this.updateSurfaceMesh(this.membraneMesh, state, sim, "membrane");
      this.updateSurfaceMesh(this.targetMesh, state, sim, "target");
      if (isolated) this.clearGroup(this.surfaceContourRoot);
      else this.renderSurfaceContours(state, sim);
    }

    createSurfaceMesh(state, sim, type) {
      const T = this.THREE;
      const geometry = this.buildSurfaceGeometry(state, sim, type);
      const material = type === "target" ? this.materials.target : this.materials.membrane.clone();
      if (type !== "target") material.vertexColors = true;
      const mesh = new T.Mesh(geometry, material);
      mesh.receiveShadow = type !== "target";
      mesh.userData.disposeGeometry = true;
      if (type !== "target") mesh.userData.disposeMaterial = true;
      return mesh;
    }

    updateSurfaceMesh(mesh, state, sim, type) {
      const next = this.surfaceArrays(state, sim, type);
      const position = mesh.geometry.getAttribute("position");
      if (!position || position.count !== next.vertices.length / 3) {
        mesh.geometry.dispose();
        mesh.geometry = this.buildSurfaceGeometry(state, sim, type);
        return;
      }
      for (let i = 0; i < next.vertices.length / 3; i += 1) {
        position.setXYZ(i, next.vertices[i * 3], next.vertices[i * 3 + 1], next.vertices[i * 3 + 2]);
      }
      position.needsUpdate = true;
      const colorAttr = mesh.geometry.getAttribute("color");
      if (colorAttr && next.colors.length) {
        for (let i = 0; i < next.colors.length / 3; i += 1) {
          colorAttr.setXYZ(i, next.colors[i * 3], next.colors[i * 3 + 1], next.colors[i * 3 + 2]);
        }
        colorAttr.needsUpdate = true;
      }
      mesh.geometry.computeVertexNormals();
    }

    buildSurfaceGeometry(state, sim, type) {
      const T = this.THREE;
      const geometry = new T.BufferGeometry();
      const { vertices, colors, indices } = this.surfaceArrays(state, sim, type);
      geometry.setAttribute("position", new T.Float32BufferAttribute(vertices, 3));
      if (colors.length) geometry.setAttribute("color", new T.Float32BufferAttribute(colors, 3));
      geometry.setIndex(indices);
      geometry.computeVertexNormals();
      return geometry;
    }

    surfaceArrays(state, sim, type) {
      const T = this.THREE;
      const { rows, cols } = state.grid;
      const vertices = [];
      const colors = [];
      const indices = [];
      const color = new T.Color();
      const subdivisions = Math.max(2, Math.min(8, Number(state.view.surfaceSubdivisions || 5)));
      const sampleRows = (rows - 1) * subdivisions + 1;
      const sampleCols = (cols - 1) * subdivisions + 1;
      const smooth = state.view.surfaceInterpolation !== "linear";
      for (let sr = 0; sr < sampleRows; sr += 1) {
        const gr = sr / subdivisions;
        for (let sc = 0; sc < sampleCols; sc += 1) {
          const gc = sc / subdivisions;
          const p = this.samplePoint(sim.centers, gr, gc, smooth);
          const alpha = this.sampleScalarGrid(sim.alpha, gr, gc, smooth);
          const target = this.sampleScalarGrid(sim.target, gr, gc, smooth);
          vertices.push(p.x, p.y, (type === "target" ? target + 0.2 : p.z + 0.14));
          const t = (alpha - state.grid.alphaMin) / (state.grid.alphaMax - state.grid.alphaMin);
          color.setHSL(0.56 - 0.44 * t, 0.55, 0.56);
          if (type !== "target") colors.push(color.r, color.g, color.b);
        }
      }
      for (let r = 0; r < sampleRows - 1; r += 1) {
        for (let c = 0; c < sampleCols - 1; c += 1) {
          const a = r * sampleCols + c;
          indices.push(a, a + 1, a + sampleCols, a + 1, a + sampleCols + 1, a + sampleCols);
        }
      }
      return { vertices, colors, indices };
    }

    renderSurfaceContours(state, sim) {
      this.clearGroup(this.surfaceContourRoot);
      if (state.view.surfaceContoursVisible !== true) return;
      const T = this.THREE;
      const { rows, cols } = state.grid;
      const subdivisions = Math.max(2, Math.min(8, Number(state.view.surfaceSubdivisions || 5)));
      const mode = state.view.contourMode || "membrane";
      const smooth = state.view.surfaceInterpolation !== "linear";
      const rowSamples = (cols - 1) * subdivisions + 1;
      const colSamples = (rows - 1) * subdivisions + 1;
      for (let r = 0; r < rows; r += 1) {
        const points = [];
        let errorSum = 0;
        for (let i = 0; i < rowSamples; i += 1) {
          const gc = i / subdivisions;
          const p = this.surfaceContourPoint(state, sim, r, gc, smooth, mode);
          points.push(new T.Vector3(p.x, p.y, p.z));
          errorSum += p.error;
        }
        this.addSurfaceContourLine(points, mode, errorSum / rowSamples);
      }
      for (let c = 0; c < cols; c += 1) {
        const points = [];
        let errorSum = 0;
        for (let i = 0; i < colSamples; i += 1) {
          const gr = i / subdivisions;
          const p = this.surfaceContourPoint(state, sim, gr, c, smooth, mode);
          points.push(new T.Vector3(p.x, p.y, p.z));
          errorSum += p.error;
        }
        this.addSurfaceContourLine(points, mode, errorSum / colSamples);
      }
    }

    surfaceContourPoint(state, sim, gr, gc, smooth, mode) {
      const p = this.samplePoint(sim.centers, gr, gc, smooth);
      const target = this.sampleScalarGrid(sim.target, gr, gc, smooth);
      const error = this.sampleScalarGrid(sim.targetError, gr, gc, smooth);
      if (mode === "target") return { x: p.x, y: p.y, z: target + 0.24, error };
      return { x: p.x, y: p.y, z: p.z + 0.2, error };
    }

    addSurfaceContourLine(points, mode, meanError) {
      const T = this.THREE;
      const geometry = new T.BufferGeometry().setFromPoints(points);
      geometry.userData.disposeGeometry = true;
      const material = mode === "target" ? this.materials.targetContour : mode === "error" && meanError >= 0 ? this.materials.errorPositive : mode === "error" ? this.materials.errorNegative : this.materials.contour;
      const line = new T.Line(geometry, material);
      line.userData.disposeGeometry = true;
      this.surfaceContourRoot.add(line);
    }

    sampleScalarGrid(grid, r, c, smooth) {
      return this.sampleScalarField(grid, r, c, smooth, (value) => value);
    }

    sampleScalarField(grid, r, c, smooth, read) {
      if (!smooth || grid.length < 3 || grid[0].length < 3) return this.bilinearField(grid, r, c, read);
      const r1 = Math.floor(r);
      const c1 = Math.floor(c);
      const tr = r - r1;
      const tc = c - c1;
      const rows = grid.length;
      const cols = grid[0].length;
      const rowSamples = [];
      for (let dr = -1; dr <= 2; dr += 1) {
        const rr = Math.max(0, Math.min(rows - 1, r1 + dr));
        const p0 = read(grid[rr][Math.max(0, Math.min(cols - 1, c1 - 1))]);
        const p1 = read(grid[rr][Math.max(0, Math.min(cols - 1, c1))]);
        const p2 = read(grid[rr][Math.max(0, Math.min(cols - 1, c1 + 1))]);
        const p3 = read(grid[rr][Math.max(0, Math.min(cols - 1, c1 + 2))]);
        rowSamples.push(this.catmullRom(p0, p1, p2, p3, tc));
      }
      return this.catmullRom(rowSamples[0], rowSamples[1], rowSamples[2], rowSamples[3], tr);
    }

    bilinearField(grid, r, c, read) {
      const rows = grid.length;
      const cols = grid[0].length;
      const r0 = Math.min(rows - 2, Math.max(0, Math.floor(r)));
      const c0 = Math.min(cols - 2, Math.max(0, Math.floor(c)));
      const r1 = Math.min(rows - 1, r0 + 1);
      const c1 = Math.min(cols - 1, c0 + 1);
      const tr = rows === 1 ? 0 : r - r0;
      const tc = cols === 1 ? 0 : c - c0;
      return this.bilerpScalar(read(grid[r0][c0]), read(grid[r0][c1]), read(grid[r1][c0]), read(grid[r1][c1]), tr, tc);
    }

    samplePoint(grid, r, c, smooth) {
      return {
        x: this.sampleScalarField(grid, r, c, smooth, (p) => p.x),
        y: this.sampleScalarField(grid, r, c, smooth, (p) => p.y),
        z: this.sampleScalarField(grid, r, c, smooth, (p) => p.z),
      };
    }

    catmullRom(p0, p1, p2, p3, t) {
      const t2 = t * t;
      const t3 = t2 * t;
      return 0.5 * (2 * p1 + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2 + (-p0 + 3 * p1 - 3 * p2 + p3) * t3);
    }

    renderErrorRods(state, sim) {
      const T = this.THREE;
      const { rows, cols } = state.grid;
      const errors = [];
      for (let r = 0; r < rows; r += 1) {
        for (let c = 0; c < cols; c += 1) {
          errors.push({ r, c, error: sim.targetError[r][c], magnitude: Math.abs(sim.targetError[r][c]) });
        }
      }
      errors.sort((a, b) => b.magnitude - a.magnitude);
      const keep = new Set(errors.slice(0, Math.min(18, errors.length)).map((item) => `${item.r},${item.c}`));
      keep.add(`${state.selection.r},${state.selection.c}`);
      for (const key of keep) {
        const [r, c] = key.split(",").map(Number);
        const p = sim.centers[r][c];
        const actual = new T.Vector3(p.x, p.y, p.z + 0.18);
        const target = new T.Vector3(p.x, p.y, sim.target[r][c] + 0.2);
        const geometry = new T.BufferGeometry().setFromPoints([actual, target]);
        geometry.userData.disposeGeometry = true;
        const line = new T.Line(geometry, sim.targetError[r][c] >= 0 ? this.materials.errorPositive : this.materials.errorNegative);
        line.userData.disposeGeometry = true;
        this.errorRoot.add(line);
      }
    }

    bilerpScalar(a00, a01, a10, a11, tr, tc) {
      return a00 * (1 - tr) * (1 - tc) + a01 * (1 - tr) * tc + a10 * tr * (1 - tc) + a11 * tr * tc;
    }

    bilerpPoint(p00, p01, p10, p11, tr, tc) {
      return {
        x: this.bilerpScalar(p00.x, p01.x, p10.x, p11.x, tr, tc),
        y: this.bilerpScalar(p00.y, p01.y, p10.y, p11.y, tr, tc),
        z: this.bilerpScalar(p00.z, p01.z, p10.z, p11.z, tr, tc),
      };
    }

    makeLabel(text) {
      const T = this.THREE;
      const canvas = document.createElement("canvas");
      canvas.width = 300;
      canvas.height = 168;
      const ctx = canvas.getContext("2d");
      ctx.fillStyle = "rgba(255,255,255,0.84)";
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.strokeStyle = "rgba(45,60,78,0.34)";
      ctx.strokeRect(1, 1, canvas.width - 2, canvas.height - 2);
      ctx.fillStyle = "#18202b";
      ctx.font = "22px system-ui, sans-serif";
      for (const [i, line] of text.split("\n").slice(0, 5).entries()) ctx.fillText(line, 14, 30 + i * 27);
      const texture = new T.CanvasTexture(canvas);
      const material = new T.SpriteMaterial({ map: texture, transparent: true });
      const sprite = new T.Sprite(material);
      sprite.userData.disposeMaterial = true;
      sprite.scale.set(1.82, 1.02, 1);
      return sprite;
    }

    makePartCalloutLabel(text) {
      const T = this.THREE;
      const canvas = document.createElement("canvas");
      canvas.width = 256;
      canvas.height = 72;
      const ctx = canvas.getContext("2d");
      ctx.fillStyle = "rgba(255,255,255,0.92)";
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.strokeStyle = "rgba(31,42,54,0.35)";
      ctx.strokeRect(1, 1, canvas.width - 2, canvas.height - 2);
      ctx.fillStyle = "#1f2a36";
      ctx.font = "700 24px Inter, Arial, sans-serif";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(text, canvas.width / 2, canvas.height / 2);
      const texture = new T.CanvasTexture(canvas);
      const material = new T.SpriteMaterial({ map: texture, transparent: true });
      const sprite = new T.Sprite(material);
      sprite.userData.disposeMaterial = true;
      sprite.scale.set(0.96, 0.27, 1);
      return sprite;
    }

    renderExplodedPartCallouts(state, sim) {
      if (state.view.explodedSelected !== true || (state.view.cellVisualMode || "abstract") === "abstract") return;
      const { r, c } = state.selection;
      const record = this.cellGroups.get(`${r},${c}`);
      if (!record || !record.group.visible) return;
      const T = this.THREE;
      const center = sim.centers[r][c];
      const groupOrigin = new T.Vector3(center.x, center.y, center.z);
      const plate = record.plates[0].position.clone().add(groupOrigin);
      const hinge = record.hinges[0].position.clone().add(groupOrigin);
      const gap = record.gap.position.clone();
      const actuator = record.actuator.group.visible ? record.actuator.slider.position.clone().add(record.actuator.group.position) : null;
      const alphaActuator = record.actuator.group.visible ? record.actuator.alphaSleeve.position.clone().add(record.actuator.group.position) : null;
      const callouts = [
        { text: "rotating plate", anchor: plate, label: plate.clone().add(new T.Vector3(-0.88, -0.46, 0.46)) },
        { text: "hinge pin", anchor: hinge, label: hinge.clone().add(new T.Vector3(-0.78, 0.36, 0.54)) },
        { text: "backlash gap", anchor: gap, label: gap.clone().add(new T.Vector3(0.92, 0.22, 0.42)) },
      ];
      if (actuator) callouts.push({ text: "z actuator", anchor: actuator, label: actuator.clone().add(new T.Vector3(0.84, -0.24, 0.3)) });
      if (alphaActuator) callouts.push({ text: "alpha actuator", anchor: alphaActuator, label: alphaActuator.clone().add(new T.Vector3(0.78, -0.56, 0.18)) });
      for (const callout of callouts) this.addPartCallout(callout.anchor, callout.label, callout.text);
    }

    addPartCallout(anchor, labelPosition, text) {
      const T = this.THREE;
      const leader = new T.Line(new T.BufferGeometry().setFromPoints([anchor, labelPosition]), this.materials.partCallout);
      leader.userData.disposeGeometry = true;
      const label = this.makePartCalloutLabel(text);
      label.position.copy(labelPosition);
      this.partCalloutRoot.add(leader, label);
    }

    renderMeasurements(state, sim) {
      const { r, c } = state.selection;
      const center = sim.centers[r][c];
      const measurementMode = state.view.measurementMode || "all";
      const T = this.THREE;
      if (state.view.measurementLabelsVisible === true) {
        const label = this.makeLabel(this.measurementLabel(state, sim, r, c, measurementMode));
        const labelPosition = this.measurementLabelPosition(state, center);
        label.position.copy(labelPosition);
        const labelLeader = new T.Line(
          new T.BufferGeometry().setFromPoints([
            new T.Vector3(center.x, center.y, center.z + 0.32),
            new T.Vector3(labelPosition.x, labelPosition.y, labelPosition.z - 0.18),
          ]),
          this.materials.measurement
        );
        labelLeader.userData.disposeGeometry = true;
        this.measurementRoot.add(labelLeader, label);
      }
      const ringRadius = this.measurementRingRadius(state, sim, r, c, measurementMode);
      const lineGeometry = new T.BufferGeometry().setFromPoints([
        new T.Vector3(center.x, center.y, 0),
        new T.Vector3(center.x, center.y, center.z),
      ]);
      const heightLine = new T.Line(lineGeometry, this.materials.measurement);
      heightLine.userData.disposeGeometry = true;
      this.measurementRoot.add(heightLine);
      const ring = new T.Line(
        new T.BufferGeometry().setFromPoints(this.selectionRingPoints(center, ringRadius)),
        this.materials.selection
      );
      ring.userData.disposeGeometry = true;
      this.measurementRoot.add(ring);
      if (measurementMode === "theta" || measurementMode === "all") this.renderThetaGuide(state, sim, r, c, center);
      if (measurementMode === "backlash" || measurementMode === "all") this.renderBacklashGuide(state, sim, r, c, center);
      if (measurementMode === "height" || measurementMode === "alpha" || measurementMode === "all") this.renderActuatorTravelGuide(state, sim, r, c, center);
      const xAxis = new T.Line(
        new T.BufferGeometry().setFromPoints([
          new T.Vector3(center.x, center.y, center.z + 0.18),
          new T.Vector3(center.x + 0.62, center.y, center.z + 0.18),
        ]),
        this.materials.errorPositive
      );
      xAxis.userData.disposeGeometry = true;
      const yAxis = new T.Line(
        new T.BufferGeometry().setFromPoints([
          new T.Vector3(center.x, center.y, center.z + 0.18),
          new T.Vector3(center.x, center.y + 0.62, center.z + 0.18),
        ]),
        this.materials.errorNegative
      );
      yAxis.userData.disposeGeometry = true;
      this.measurementRoot.add(xAxis, yAxis);
    }

    measurementLabelPosition(state, center) {
      const T = this.THREE;
      const { rows, cols, cellSize } = state.grid;
      const sideX = center.x >= 0 ? -1 : 1;
      const sideY = center.y >= 0 ? -1 : 1;
      const halfX = Math.max(0.5, (cols - 1) * cellSize * 0.5);
      const halfY = Math.max(0.5, (rows - 1) * cellSize * 0.5);
      const margin = Math.max(1.05, cellSize * 1.2);
      const x = sideX * (halfX + margin);
      const y = sideY * (halfY + margin * 0.78);
      const z = Math.max(center.z + 0.84, 0.86);
      return new T.Vector3(x, y, z);
    }

    renderActuatorTravelGuide(state, sim, r, c, center) {
      const T = this.THREE;
      const limits = RAD.commandLimits(state);
      const commandZ = state.cells.commandZ[r][c];
      const commandAlpha = state.cells.commandAlpha[r][c];
      const zScale = 0.62;
      const zBase = center.z + 0.18;
      const zX = center.x + 0.78;
      const zActual = zBase + (commandZ / Math.max(1e-9, limits.z)) * zScale;
      this.addMeasurementLine(
        [new T.Vector3(zX, center.y, zBase - zScale), new T.Vector3(zX, center.y, zBase + zScale)],
        this.materials.actuatorEnvelope
      );
      for (const fraction of [-1, 0, 1]) {
        const z = zBase + fraction * zScale;
        this.addMeasurementLine(
          [new T.Vector3(zX - 0.12, center.y, z), new T.Vector3(zX + 0.12, center.y, z)],
          fraction === 0 ? this.materials.measurement : this.materials.actuatorEnvelope
        );
      }
      this.addMeasurementLine(
        [new T.Vector3(zX - 0.18, center.y, zActual), new T.Vector3(zX + 0.18, center.y, zActual)],
        commandZ >= 0 ? this.materials.errorPositive : this.materials.errorNegative
      );

      const alphaScale = 0.68;
      const alphaY = center.y - 0.82;
      const alphaZ = center.z + 0.3;
      const alphaLimit = commandAlpha < 0 ? limits.alphaContract : limits.alphaExpand;
      const alphaActual = center.x + (commandAlpha / Math.max(1e-9, alphaLimit)) * alphaScale;
      this.addMeasurementLine(
        [new T.Vector3(center.x - alphaScale, alphaY, alphaZ), new T.Vector3(center.x + alphaScale, alphaY, alphaZ)],
        this.materials.actuatorEnvelope
      );
      for (const fraction of [-1, 0, 1]) {
        const x = center.x + fraction * alphaScale;
        this.addMeasurementLine(
          [new T.Vector3(x, alphaY - 0.12, alphaZ), new T.Vector3(x, alphaY + 0.12, alphaZ)],
          fraction === 0 ? this.materials.measurement : this.materials.actuatorEnvelope
        );
      }
      this.addMeasurementLine(
        [new T.Vector3(alphaActual, alphaY - 0.16, alphaZ), new T.Vector3(alphaActual, alphaY + 0.16, alphaZ)],
        commandAlpha >= 0 ? this.materials.errorPositive : this.materials.errorNegative
      );
    }

    addMeasurementLine(points, material) {
      const T = this.THREE;
      const line = new T.Line(new T.BufferGeometry().setFromPoints(points), material);
      line.userData.disposeGeometry = true;
      this.measurementRoot.add(line);
      return line;
    }

    measurementLabel(state, sim, r, c, mode) {
      const residual = sim.target[r][c] - sim.height[r][c];
      const alpha = sim.alpha[r][c];
      const theta = sim.theta[r][c];
      const height = sim.height[r][c];
      const zResidual = sim.zResidual?.[r]?.[c] || 0;
      const backlash = state.grid.backlash;
      const influence = sim.influence[r][c];
      const linkStrain = this.localLinkStrainStrength(sim, r, c) * Math.max(0.025, sim.metrics?.maxAbsLinkStrain || 0);
      const referenceDisplacement = sim.displacement?.[r]?.[c] || 0;
      const slope = sim.slope?.magnitude?.[r]?.[c] || 0;
      const normalTilt = sim.slope?.tilt?.[r]?.[c] || 0;
      if (mode === "alpha") {
        const limits = RAD.commandLimits(state);
        const alphaLimit = state.cells.commandAlpha[r][c] < 0 ? limits.alphaContract : limits.alphaExpand;
        const alphaTravel = Math.abs(state.cells.commandAlpha[r][c]) / Math.max(1e-9, alphaLimit);
        return `alpha ${alpha.toFixed(3)}\ncommand ${state.cells.commandAlpha[r][c].toFixed(2)}\nalpha travel ${alphaTravel.toFixed(2)}\nrange ${state.grid.alphaMin.toFixed(2)}-${state.grid.alphaMax.toFixed(2)}\ntheta law ${theta.toFixed(1)} deg`;
      }
      if (mode === "theta") {
        return `theta ${theta.toFixed(1)} deg\nlaw: theta = 70 alpha - 60\nalpha ${alpha.toFixed(3)}\ncmd alpha ${state.cells.commandAlpha[r][c].toFixed(2)}`;
      }
      if (mode === "height") {
        return `z ${height.toFixed(3)}\ncmd z ${state.cells.commandZ[r][c].toFixed(2)}\nz residual ${zResidual.toFixed(3)}\ntarget ${sim.target[r][c].toFixed(3)}\nresidual ${residual.toFixed(3)}\nmean signed ${sim.metrics.meanSignedTargetError.toFixed(3)}`;
      }
      if (mode === "backlash") {
        return `backlash b ${backlash.toFixed(3)}\ninfluence ${influence.toFixed(3)}\ndie-off ${Number.isFinite(sim.dieOff[r][c]) ? sim.dieOff[r][c] : "locked"}\ndead-zone +/-${backlash.toFixed(2)}`;
      }
      return `alpha ${alpha.toFixed(3)}  theta ${theta.toFixed(1)}\nz ${height.toFixed(3)}  z residual ${zResidual.toFixed(3)}\nnormal tilt ${normalTilt.toFixed(1)}  residual ${residual.toFixed(3)}\ninfluence ${influence.toFixed(3)}  link strain ${linkStrain.toFixed(3)}\ncmd a ${state.cells.commandAlpha[r][c].toFixed(2)}  cmd z ${state.cells.commandZ[r][c].toFixed(2)}\nb ${backlash.toFixed(2)}  ref disp ${referenceDisplacement.toFixed(3)}`;
    }

    measurementRingRadius(state, sim, r, c, mode) {
      if (mode === "alpha") return 0.32 + Math.sqrt(Math.max(0.01, sim.alpha[r][c])) * 0.18;
      if (mode === "backlash") return 0.45 + state.grid.backlash * 0.55;
      return 0.48;
    }

    renderThetaGuide(state, sim, r, c, center) {
      const T = this.THREE;
      const theta = (sim.theta[r][c] * Math.PI) / 180;
      const radius = 0.72;
      const points = [new T.Vector3(center.x, center.y, center.z + 0.24)];
      const steps = 16;
      for (let i = 0; i <= steps; i += 1) {
        const a = (theta * i) / steps;
        points.push(new T.Vector3(center.x + Math.cos(a) * radius, center.y + Math.sin(a) * radius, center.z + 0.24));
      }
      const arc = new T.Line(new T.BufferGeometry().setFromPoints(points), theta >= 0 ? this.materials.errorPositive : this.materials.errorNegative);
      arc.userData.disposeGeometry = true;
      this.measurementRoot.add(arc);
    }

    renderBacklashGuide(state, sim, r, c, center) {
      const T = this.THREE;
      const radius = 0.45 + state.grid.backlash * 0.55;
      const clearance = Math.max(0.08, state.grid.backlash * 0.8);
      const base = (sim.theta[r][c] * Math.PI) / 180 + Math.PI / 4;
      for (const side of [-1, 1]) {
        const angle = base + side * clearance;
        const line = new T.Line(
          new T.BufferGeometry().setFromPoints([
            new T.Vector3(center.x, center.y, center.z + 0.32),
            new T.Vector3(center.x + Math.cos(angle) * radius, center.y + Math.sin(angle) * radius, center.z + 0.32),
          ]),
          this.materials.gapLine
        );
        line.userData.disposeGeometry = true;
        this.measurementRoot.add(line);
      }
    }

    selectionRingPoints(center, radius) {
      const T = this.THREE;
      const points = [];
      for (let i = 0; i <= 48; i += 1) {
        const a = (i / 48) * Math.PI * 2;
        points.push(new T.Vector3(center.x + Math.cos(a) * radius, center.y + Math.sin(a) * radius, center.z + 0.16));
      }
      return points;
    }

    fitRoot(rows, cols) {
      const desired = Math.max(rows, cols) * 1.75;
      if (Math.abs(this.spherical.radius - desired) > 4) {
        this.spherical.radius = Math.max(7, desired);
        this.updateCamera();
      }
    }

    pick(event) {
      const rect = this.renderer.domElement.getBoundingClientRect();
      this.pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
      this.pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
      this.raycaster.setFromCamera(this.pointer, this.camera);
      const hits = this.raycaster.intersectObjects(this.selectables, false);
      if (hits.length && hits[0].object.userData.selectable) {
        this.onSelect(hits[0].object.userData.r, hits[0].object.userData.c);
      }
    }

    hoverPick(event) {
      const rect = this.renderer.domElement.getBoundingClientRect();
      this.pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
      this.pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
      this.raycaster.setFromCamera(this.pointer, this.camera);
      const hits = this.raycaster.intersectObjects(this.selectables, false);
      if (hits.length && hits[0].object.userData.selectable) {
        this.setHoveredCell({ r: hits[0].object.userData.r, c: hits[0].object.userData.c });
      } else {
        this.setHoveredCell(null);
      }
    }

    setHoveredCell(cell) {
      const changed = (this.hoveredCell?.r ?? null) !== (cell?.r ?? null) || (this.hoveredCell?.c ?? null) !== (cell?.c ?? null);
      if (!changed) return;
      this.hoveredCell = cell;
      this.renderer.domElement.style.cursor = cell ? "pointer" : "";
      if (this.displaySim) this.updateDynamicObjects(this.displaySim);
      if (typeof this.onHover === "function") this.onHover(cell);
    }

    animate() {
      requestAnimationFrame(() => this.animate());
      if (this.stepDisplaySim()) {
        const now = performance.now();
        if (now - this.lastRebuild > 32) {
          this.rebuild(this.displaySim);
          this.lastRebuild = now;
        }
      }
      this.renderer.render(this.scene, this.camera);
      this.updateDiagnostics();
    }

    updateDiagnostics() {
      this.frameCount += 1;
      const now = performance.now();
      if (now - this.lastFpsTime < 500) return;
      const info = this.renderer.info;
      this.diagnostics = {
        fps: Math.round((this.frameCount * 1000) / (now - this.lastFpsTime)),
        drawCalls: info.render.calls,
        triangles: info.render.triangles,
        geometries: info.memory.geometries,
        textures: info.memory.textures,
        cells: this.cellGroups.size,
        surfaceVertices:
          (this.membraneMesh?.geometry?.getAttribute("position")?.count || 0) +
          (this.targetMesh?.geometry?.getAttribute("position")?.count || 0),
      };
      this.frameCount = 0;
      this.lastFpsTime = now;
    }

    getDiagnostics() {
      return { ...this.diagnostics };
    }
  }

  RAD.RadRenderer = RadRenderer;
})();
