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
  const pinRowsInput = document.getElementById("pinRows");
  const reachBands = document.getElementById("reachBands");
  const realizedMarker = document.getElementById("realizedMarker");
  const driveRealizedMetric = document.getElementById("driveRealizedMetric");
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
  const circumferentialPinMetric = document.getElementById("circumferentialPinMetric");
  const axialPinMetric = document.getElementById("axialPinMetric");
  const jointTiltLimitMetric = document.getElementById("jointTiltLimitMetric");
  const circumferentialBendMetric = document.getElementById("circumferentialBendMetric");
  const axialBendMetric = document.getElementById("axialBendMetric");
  const minCellsClosureMetric = document.getElementById("minCellsClosureMetric");
  const constraintsEnabledInput = document.getElementById("constraintsEnabled");
  const constraintStateMetric = document.getElementById("constraintStateMetric");
  const realizedAlphaMetric = document.getElementById("realizedAlphaMetric");
  const feasibleRangeMetric = document.getElementById("feasibleRangeMetric");
  const pinDiameterInput = document.getElementById("pinDiameter");
  const pinDiameterOut = document.getElementById("pinDiameterOut");
  const showPinsInput = document.getElementById("showPins");
  const pinClearanceMetric = document.getElementById("pinClearanceMetric");
  const pinBOverLMetric = document.getElementById("pinBOverLMetric");
  const pinDeltaPhiMetric = document.getElementById("pinDeltaPhiMetric");
  const pinAlphaDeadzoneMetric = document.getElementById("pinAlphaDeadzoneMetric");
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

  // Max tilt of one plate (thickness t) on a pin (diameter d) in a hole
  // (diameter D), exact contact form from the pin/backlash one-pager:
  // d + t*sin(theta) = D*cos(theta)  =>  theta = acos(d/sqrt(D^2+t^2)) - atan(t/D).
  // ~= b/t for b = D - d << t.
  function pinTiltLimitRad(d, D, t) {
    return Math.acos(d / Math.sqrt(D * D + t * t)) - Math.atan(t / D);
  }
  // Pin sized as a fraction of the hole, with backlash derived from the
  // clearance - concept from Ahyan's two-cell attachment V2 (pin radius
  // ratio -> radial clearance b -> b/L -> dphi = asin(b/L)). Default 3.0 mm
  // is an assumed M3 bolt (hardware photos show bolts through the pads);
  // not measured - confirm against the real parts.
  function pinDiameterMm() {
    const value = Number(pinDiameterInput.value);
    return clamp(Number.isFinite(value) ? value : 3.0, Number(pinDiameterInput.min), CAD.nominalHoleDiameterMm);
  }

  function pinRadialClearanceMm() {
    return Math.max(0, (CAD.nominalHoleDiameterMm - pinDiameterMm()) / 2);
  }

  // In-plane angular dead zone of one arm swinging about its neighbor's pin.
  function pinDeltaPhiRad() {
    return Math.asin(clamp(pinRadialClearanceMm() / CAD.siteRadiusMm, 0, 1));
  }

  // A joint is two plates (one per cell) on one floating bolt, so their
  // relative tilt limit is the sum of each plate's own limit.
  function jointTiltLimitRad() {
    return 2 * Math.max(0, pinTiltLimitRad(pinDiameterMm(), CAD.nominalHoleDiameterMm, CAD.bodyThicknessMm));
  }

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

  // --- 3D material collision ("no fusing through other solids") ---
  // Each cross is modeled at its real thickness as a union of flat parts:
  // the hub disc, four pad discs and four half-arms (hub -> pad), each a 2D
  // shape (disc, or rounded segment for an arm) extruded through the
  // plate's thickness along the cell's own normal, in the cell's actual 3D
  // pose. Penetration between two parts is estimated by sampling points
  // across each part (rim, interior, three depths) and evaluating the other
  // part's exact signed distance. The two pads one joint pin passes
  // through, and the arms leading to them, are exempt from each other:
  // whether that joint can bend far enough is the backlash tilt check.
  const COLLISION_TOLERANCE_MM = 0.05;
  const PLATE_HALF_THICKNESS = CAD.bodyThicknessMm / 2;

  function partSamples(part) {
    const inPlane = [];
    if (part.kind === "disc") {
      inPlane.push([part.a.x, part.a.y]);
      for (let k = 0; k < 8; k += 1) {
        const angle = (k * Math.PI) / 4;
        inPlane.push([part.a.x + 0.995 * part.r * Math.cos(angle), part.a.y + 0.995 * part.r * Math.sin(angle)]);
      }
    } else {
      const dx = part.b.x - part.a.x;
      const dy = part.b.y - part.a.y;
      const length = Math.hypot(dx, dy);
      const px = -dy / length;
      const py = dx / length;
      [0.25, 0.5, 0.75].forEach((f) => {
        const cx = part.a.x + f * dx;
        const cy = part.a.y + f * dy;
        [0, 0.995, -0.995].forEach((o) => inPlane.push([cx + o * part.r * px, cy + o * part.r * py]));
      });
    }
    const samples = [];
    [-0.95, 0, 0.95].forEach((depth) => inPlane.forEach(([u, v]) => samples.push([u, v, depth * PLATE_HALF_THICKNESS])));
    return samples;
  }

  const CROSS_PARTS = (() => {
    const parts = [{ kind: "disc", site: null, a: new THREE.Vector2(0, 0), b: new THREE.Vector2(0, 0), r: CAD.hubRadiusMm }];
    Object.entries(SITE_VECTORS).forEach(([site, vector]) => {
      parts.push({ kind: "disc", site, a: vector.clone(), b: vector.clone(), r: CAD.padRadiusMm });
      parts.push({ kind: "arm", site, a: new THREE.Vector2(0, 0), b: vector.clone(), r: CAD.armWidthMm / 2 });
    });
    return parts.map((part) => {
      const halfLength = part.a.distanceTo(part.b) / 2;
      return {
        ...part,
        samples: partSamples(part),
        mid: part.a.clone().add(part.b).multiplyScalar(0.5),
        bound: Math.hypot(halfLength + part.r, PLATE_HALF_THICKNESS),
      };
    });
  })();
  const CELL_BOUND = Math.hypot(CAD.siteRadiusMm + CAD.padRadiusMm, CAD.bodyThicknessMm);

  function segmentDistance2(px, py, ax, ay, bx, by) {
    const dx = bx - ax;
    const dy = by - ay;
    const len2 = dx * dx + dy * dy;
    let t = len2 > 1e-12 ? ((px - ax) * dx + (py - ay) * dy) / len2 : 0;
    t = t < 0 ? 0 : t > 1 ? 1 : t;
    const qx = ax + t * dx - px;
    const qy = ay + t * dy - py;
    return Math.sqrt(qx * qx + qy * qy);
  }

  // Exact signed distance to a flat part in its own layer coordinates
  // (u, v in the plate, w along the normal from the plate's mid-plane).
  function partSignedDistance(part, u, v, w) {
    const inPlane = segmentDistance2(u, v, part.a.x, part.a.y, part.b.x, part.b.y) - part.r;
    const through = Math.abs(w) - PLATE_HALF_THICKNESS;
    if (inPlane <= 0 && through <= 0) return Math.max(inPlane, through);
    const a = Math.max(inPlane, 0);
    const b = Math.max(through, 0);
    return Math.sqrt(a * a + b * b);
  }

  // Both crosses of one cell as world-space layer frames: origin on the
  // layer's mid-plane at the hub, e1/e2 the cross's own in-plane axes
  // (twisted by crossRot about the normal), e3 the outward normal.
  function cellLayers(ring, i) {
    const t = ring.tangents[i];
    const nrm = ring.normals[i];
    const hub = new THREE.Vector3(ring.centers[i].x, ring.centers[i].y, ring.baseZ + ring.hubZ[i]);
    return ["lower", "upper"].map((layer) => {
      const phi = layer === "upper" ? ring.crossRotUpper[i] : ring.crossRotLower[i];
      const offset = layer === "upper" ? PLATE_HALF_THICKNESS : -PLATE_HALF_THICKNESS;
      const c = Math.cos(phi);
      const s = Math.sin(phi);
      const e1 = new THREE.Vector3(t.x * c, t.y * c, s);
      const e2 = new THREE.Vector3(-t.x * s, -t.y * s, c);
      const e3 = new THREE.Vector3(nrm.x, nrm.y, 0);
      const origin = hub.clone().addScaledVector(e3, offset);
      const parts = CROSS_PARTS.map((part) => {
        const center = origin.clone().addScaledVector(e1, part.mid.x).addScaledVector(e2, part.mid.y);
        const worldSamples = part.samples.map(([u, v, w]) =>
          origin.clone().addScaledVector(e1, u).addScaledVector(e2, v).addScaledVector(e3, w)
        );
        return { part, center, worldSamples };
      });
      return { layer, origin, e1, e2, e3, parts };
    });
  }

  function toLayer(layer, point) {
    const dx = point.x - layer.origin.x;
    const dy = point.y - layer.origin.y;
    const dz = point.z - layer.origin.z;
    return [
      dx * layer.e1.x + dy * layer.e1.y + dz * layer.e1.z,
      dx * layer.e2.x + dy * layer.e2.y + dz * layer.e2.z,
      dx * layer.e3.x + dy * layer.e3.y + dz * layer.e3.z,
    ];
  }

  // Deepest sampled penetration between two parts (positive = overlap),
  // and the sample point where it occurs.
  function partPairPenetration(layerA, entryA, layerB, entryB) {
    let depth = -Infinity;
    let where = null;
    entryA.worldSamples.forEach((point) => {
      const [u, v, w] = toLayer(layerB, point);
      const d = -partSignedDistance(entryB.part, u, v, w);
      if (d > depth) {
        depth = d;
        where = point;
      }
    });
    entryB.worldSamples.forEach((point) => {
      const [u, v, w] = toLayer(layerA, point);
      const d = -partSignedDistance(entryA.part, u, v, w);
      if (d > depth) {
        depth = d;
        where = point;
      }
    });
    return { depth, where };
  }

  // Checks every cell against its ring neighbors (one and two over) and the
  // three nearest cells in the next row, culled by bounding spheres.
  // `pinnedRows`: rows are joined, so north/south pads of stacked cells are
  // an axial joint. `regionRows`: only check pairs touching these rows.
  function checkCollisions(rings, n, m, pinnedRows, regionRows = null) {
    const aSite = aSiteInput.value;
    const bSite = bSiteInput.value;
    const layerCache = rings.map(() => new Array(n));
    const layersOf = (row, i) => layerCache[row][i] || (layerCache[row][i] = cellLayers(rings[row], i));
    const inRegion = (row) => !regionRows || regionRows.includes(row);
    const hubs = rings.map((ring) => ring.centers.map((c, i) => new THREE.Vector3(c.x, c.y, ring.baseZ + ring.hubZ[i])));
    let maxPenetration = 0;
    let minClearance = Infinity;
    let pairsChecked = 0;
    const contacts = [];

    const secondPin = hasSecondJointPin();
    // Joint pads (and the arms leading to them) are exempt from each other.
    // Circumferential joints: cell i's aSite pads with cell i+1's bSite pads.
    // Axial joints: a cell's north pads with the south pads of the cell
    // stacked above it - north/south are opposite sites, so both layer
    // pairings line up, like the second circumferential pin.
    function isJointPair(joint, layerA, entryA, layerB, entryB) {
      if (entryA.part.site !== joint.a || entryB.part.site !== joint.b) return false;
      if (layerA.layer === "upper" && layerB.layer === "lower") return true;
      return joint.both && layerA.layer === "lower" && layerB.layer === "upper";
    }
    const ringJoint = { a: aSite, b: bSite, both: secondPin };
    const axialJoint = { a: "north", b: "south", both: true };

    function checkPair(rowA, iA, rowB, iB, joint) {
      if (!inRegion(rowA) && !inRegion(rowB)) return;
      if (hubs[rowA][iA].distanceTo(hubs[rowB][iB]) > 2 * CELL_BOUND) return;
      pairsChecked += 1;
      layersOf(rowA, iA).forEach((layerA) => {
        layersOf(rowB, iB).forEach((layerB) => {
          layerA.parts.forEach((entryA) => {
            layerB.parts.forEach((entryB) => {
              if (joint && isJointPair(joint, layerA, entryA, layerB, entryB)) return;
              const gap = entryA.center.distanceTo(entryB.center) - entryA.part.bound - entryB.part.bound;
              if (gap > 0) {
                minClearance = Math.min(minClearance, gap);
                return;
              }
              const { depth, where } = partPairPenetration(layerA, entryA, layerB, entryB);
              minClearance = Math.min(minClearance, -depth);
              if (depth > COLLISION_TOLERANCE_MM) {
                maxPenetration = Math.max(maxPenetration, depth);
                contacts.push({
                  depth,
                  point: where,
                  cellA: [rowA, iA],
                  cellB: [rowB, iB],
                  partA: `${layerA.layer} ${entryA.part.kind} ${entryA.part.site || "hub"}`,
                  partB: `${layerB.layer} ${entryB.part.kind} ${entryB.part.site || "hub"}`,
                });
              }
            });
          });
        });
      });
    }

    for (let row = 0; row < m; row += 1) {
      for (let i = 0; i < n; i += 1) {
        checkPair(row, i, row, (i + 1) % n, ringJoint);
        if (n > 4) checkPair(row, i, row, (i + 2) % n, null);
        if (row + 1 < m) {
          for (let di = -1; di <= 1; di += 1) {
            checkPair(row, i, row + 1, (i + di + n) % n, di === 0 && pinnedRows ? axialJoint : null);
          }
        }
      }
    }
    contacts.sort((a, b) => b.depth - a.depth);
    return {
      maxPenetration,
      minClearance: Number.isFinite(minClearance) ? Math.max(0, minClearance) : 0,
      pairsChecked,
      contacts: contacts.slice(0, 60),
      clear: maxPenetration <= COLLISION_TOLERANCE_MM,
    };
  }

  // Real pin-hole alignment: each cell's east/west (circumferential) and
  // north/south (axial) pad centers are where a physical pin would pass
  // through into the matching neighbor. Computed directly from the actual
  // rendered per-cell frame (tangent/normal from computeCellFrame, the same
  // basis setRadialOrientation used), not the flat 2D approximation
  // checkNeighborClearance above still uses - this is what actually answers
  // "do the holes line up", which body-clearance distance alone can't.
  // Circumferential gaps are a genuine, honest residual (like ring closure)
  // rather than exactly zero: a straight cross arm can't perfectly face two
  // curved-ring neighbors at once, even with the bisector heading tangent
  // chosen to minimize it. Axial gaps are exactly zero unless rows differ
  // in diameter (differential dilation/barrel/cone shapes), in which case
  // they honestly reflect that the rows no longer stack as a true cylinder.
  // Physical pins (idea from Ahyan's two-cell attachment V2: a solid pin in
  // each hole, radius from the pin control, hideable). One through each
  // cell's hub joining its two layers, one at each circumferential joint and
  // one at each axial joint, each along the local radial direction (the
  // bolt axis for radially facing cells). Joint pins sit midway between the
  // two cells' lower-layer pads, the same pads the pin-gap readout uses.
  const pinMaterial = new THREE.MeshStandardMaterial({ color: 0x1f2933, roughness: 0.36, metalness: 0.35 });
  const pinGeometry = new THREE.CylinderGeometry(1, 1, 1, 20);
  const pinMeshes = [];
  const PIN_UP = new THREE.Vector3(0, 1, 0);
  const pinAxisScratch = new THREE.Vector3();

  function placePin(index, x, y, z, axisX, axisY, axisZ, length, radius) {
    while (pinMeshes.length <= index) {
      const mesh = new THREE.Mesh(pinGeometry, pinMaterial);
      mesh.castShadow = true;
      scene.add(mesh);
      pinMeshes.push(mesh);
    }
    const pin = pinMeshes[index];
    pinAxisScratch.set(axisX, axisY, axisZ).normalize();
    pin.quaternion.setFromUnitVectors(PIN_UP, pinAxisScratch);
    pin.position.set(x, y, z);
    pin.scale.set(radius, length, radius);
    pin.visible = true;
  }

  function updatePins(n, m, rings, cellFramesByRow, pinnedRows) {
    let count = 0;
    if (showPinsInput.checked && !isolateActive) {
      const radius = pinDiameterMm() / 2;
      const aSite = aSiteInput.value;
      const bSite = bSiteInput.value;
      const pinLayers = jointPinLayers();
      for (let row = 0; row < m; row += 1) {
        for (let i = 0; i < n; i += 1) {
          const frame = cellFramesByRow[row][i];
          const center = rings[row].centers[i];
          const hubZ = rings[row].baseZ + rings[row].hubZ[i];
          placePin(count++, center.x, center.y, hubZ, frame.normal.x, frame.normal.y, 0, CAD.bodyThicknessMm * 2.55, radius);

          // Joint pins pass through the actual pinned pad pairs.
          const next = (i + 1) % n;
          const frameB = cellFramesByRow[row][next];
          pinLayers.forEach(([layerA, layerB]) => {
            const padA = padWorld(rings[row], i, layerA, aSite);
            const padB = padWorld(rings[row], next, layerB, bSite);
            placePin(
              count++,
              (padA.x + padB.x) / 2,
              (padA.y + padB.y) / 2,
              (padA.z + padB.z) / 2,
              frame.normal.x + frameB.normal.x,
              frame.normal.y + frameB.normal.y,
              0,
              CAD.bodyThicknessMm * 2.95,
              radius
            );
          });

          // Axial joint pins between a cell's north pads and the south pads
          // of the cell above (both layer pairings, as north/south are
          // opposite sites) - only when the rows are pinned together.
          if (pinnedRows && row + 1 < m) {
            const frameUp = cellFramesByRow[row + 1][i];
            [
              ["upper", "lower"],
              ["lower", "upper"],
            ].forEach(([layerA, layerB]) => {
              const padA = padWorld(rings[row], i, layerA, "north");
              const padB = padWorld(rings[row + 1], i, layerB, "south");
              placePin(
                count++,
                (padA.x + padB.x) / 2,
                (padA.y + padB.y) / 2,
                (padA.z + padB.z) / 2,
                frame.normal.x + frameUp.normal.x,
                frame.normal.y + frameUp.normal.y,
                0,
                CAD.bodyThicknessMm * 2.95,
                radius
              );
            });
          }
        }
      }
    }
    for (let index = count; index < pinMeshes.length; index += 1) pinMeshes[index].visible = false;
    lastPinCount = count;
  }

  // Red dots where parts would overlap (the deepest sampled points). Those
  // points lie inside solid parts, so they're drawn on top of everything.
  const contactMaterial = new THREE.MeshBasicMaterial({ color: 0xff2a36, depthTest: false, transparent: true, opacity: 0.9 });
  const contactGeometry = new THREE.SphereGeometry(1.8, 12, 8);
  const contactMeshes = [];

  function updateContactMarkers(contacts) {
    contacts.forEach((contact, index) => {
      if (!contactMeshes[index]) {
        const mesh = new THREE.Mesh(contactGeometry, contactMaterial);
        mesh.renderOrder = 10;
        scene.add(mesh);
        contactMeshes.push(mesh);
      }
      contactMeshes[index].position.copy(contact.point);
      contactMeshes[index].visible = true;
    });
    for (let index = contacts.length; index < contactMeshes.length; index += 1) contactMeshes[index].visible = false;
  }

  // How far a pinned pad pair sits off one common pin axis (the bisector of
  // the two cells' normals); their spacing along the axis is just the
  // stacked plates, so this is what decides whether one pin fits both.
  function padPairOffAxis(padA, padB, normalA, normalB) {
    const axis = new THREE.Vector3(normalA.x + normalB.x, normalA.y + normalB.y, 0).normalize();
    const offset = padA.sub(padB);
    return offset.addScaledVector(axis, -offset.dot(axis)).length();
  }

  function computePinAlignment(n, m, rings, cellFramesByRow) {
    let maxCircumferential = 0;
    let maxAxial = 0;
    let maxCircumferentialBend = 0;
    let maxAxialBend = 0;
    let jointsChecked = 0;
    const aSite = aSiteInput.value;
    const bSite = bSiteInput.value;
    const pinLayers = jointPinLayers();
    for (let row = 0; row < m; row += 1) {
      for (let i = 0; i < n; i += 1) {
        const next = (i + 1) % n;
        const frameA = cellFramesByRow[row][i];
        const frameB = cellFramesByRow[row][next];
        pinLayers.forEach(([layerA, layerB]) => {
          const lateral = padPairOffAxis(
            padWorld(rings[row], i, layerA, aSite),
            padWorld(rings[row], next, layerB, bSite),
            frameA.normal,
            frameB.normal
          );
          maxCircumferential = Math.max(maxCircumferential, lateral);
        });
        // With radially oriented cells the joint bolt is radial, so the
        // angle between neighbors' normals is a tilt of the plates on that
        // bolt - the quantity the backlash tilt limit bounds.
        const bend = Math.acos(clamp(frameA.normal.dot(frameB.normal), -1, 1));
        maxCircumferentialBend = Math.max(maxCircumferentialBend, bend);
        jointsChecked += 1;
      }
    }
    for (let row = 0; row < m - 1; row += 1) {
      for (let i = 0; i < n; i += 1) {
        // A cell's north pads against the south pads of the cell above:
        // ~0 when rows are pinned and share a diameter; with a set
        // (unpinned) pitch it is the gap between the rows.
        [
          ["upper", "lower"],
          ["lower", "upper"],
        ].forEach(([layerA, layerB]) => {
          const lateral = padPairOffAxis(
            padWorld(rings[row], i, layerA, "north"),
            padWorld(rings[row + 1], i, layerB, "south"),
            rings[row].normals[i],
            rings[row + 1].normals[i]
          );
          maxAxial = Math.max(maxAxial, lateral);
        });
        // Cells stay axis-aligned, so a change in radius between rows has to
        // be taken up as tilt at the axial joint: the profile's slope angle.
        const radiusNorth = rings[row].centers[i].length();
        const radiusSouth = rings[row + 1].centers[i].length();
        const rise = rings[row + 1].baseZ + rings[row + 1].hubZ[i] - rings[row].baseZ - rings[row].hubZ[i];
        const bend = Math.atan2(Math.abs(radiusSouth - radiusNorth), rise);
        maxAxialBend = Math.max(maxAxialBend, bend);
        jointsChecked += 1;
      }
    }
    const tiltLimit = jointTiltLimitRad();
    return {
      maxCircumferential,
      maxAxial,
      maxCircumferentialBend,
      maxAxialBend,
      tiltLimit,
      circumferentialOk: maxCircumferentialBend <= tiltLimit + 1e-9,
      axialOk: maxAxialBend <= tiltLimit + 1e-9,
      // Infinity when the pin fills the hole (no clearance to tilt in).
      minCellsForBacklashClosure: tiltLimit > 1e-9 ? Math.ceil((2 * Math.PI) / tiltLimit) : Infinity,
      jointsChecked,
    };
  }

  function cylinderZ(radius, depth, material, segments = 48) {
    const mesh = new THREE.Mesh(new THREE.CylinderGeometry(radius, radius, depth, segments), material);
    mesh.rotation.x = Math.PI / 2;
    mesh.castShadow = true;
    mesh.receiveShadow = true;
    return mesh;
  }

  // Selection / batch / role markers: flat rings lying on a cell's outer
  // face around its hub, so the hub pin passes through the ring's middle.
  // Concentric radii keep all three readable on one cell at once.
  function faceRingMarker(radius, material) {
    const mesh = new THREE.Mesh(new THREE.TorusGeometry(radius, 0.9, 10, 40), material);
    mesh.renderOrder = 2;
    return mesh;
  }

  const MARKER_Z_AXIS = new THREE.Vector3(0, 0, 1);
  const markerNormalScratch = new THREE.Vector3();
  function placeFaceMarker(mesh, center, z, normal) {
    const offset = CAD.bodyThicknessMm + 0.9;
    markerNormalScratch.set(normal.x, normal.y, 0).normalize();
    mesh.position.set(center.x + normal.x * offset, center.y + normal.y * offset, z);
    mesh.quaternion.setFromUnitVectors(MARKER_Z_AXIS, markerNormalScratch);
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
        const roleMarker = faceRingMarker(8.8, sharedMaterials.actuatorMarker);
        roleMarker.visible = false;
        scene.add(roleMarker);
        cell.roleMarker = roleMarker;
        const batchMarker = faceRingMarker(7.4, sharedMaterials.batchMarker);
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

  // One cell's two circumferential joint pins, in the plane of the cell.
  // The dilation twist is split symmetrically - lower cross at -theta/2,
  // upper at +theta/2 - so the upper aSite pin (joint with the next cell)
  // and the lower bSite pin (joint with the previous cell) end up level:
  // the chord between them is rotated by psi so it lies along the cell's
  // own tangential axis. Twisting only the upper cross (lower fixed) would
  // tilt that chord by theta/2 and send the row off axially in a helix.
  // For east/west sites: pitch = 2*L*cos(theta/2), psi = 0.
  function cellJointGeometry(theta) {
    const aLocal = SITE_VECTORS[aSiteInput.value] || SITE_VECTORS.east;
    const bLocal = SITE_VECTORS[bSiteInput.value] || SITE_VECTORS.west;
    const aPin = rotate2(aLocal, theta / 2);
    const bPin = rotate2(bLocal, -theta / 2);
    const chord = aPin.clone().sub(bPin);
    const pitch = Math.max(chord.length(), 1e-6);
    const psi = chord.length() > 1e-9 ? Math.atan2(chord.y, chord.x) : 0;
    const aInFace = rotate2(aPin, -psi);
    const bInFace = rotate2(bPin, -psi);
    return {
      pitch,
      crossRotLower: -theta / 2 - psi,
      crossRotUpper: theta / 2 - psi,
      // Pin positions relative to the hub, in face coordinates
      // (x along the tangent, y along the axle). aInFace.y === bInFace.y.
      bPinAlong: bInFace.x,
      pinAxial: aInFace.y,
    };
  }

  // With opposite attachment sites (east/west or north/south) the symmetric
  // twist also brings a second pad pair exactly together at every joint:
  // cell i's LOWER aSite pad and cell i+1's UPPER bSite pad, offset from the
  // first pin by 2*L*sin(theta/2) along the axle (rot(v, t/2) - rot(v, -t/2)
  // = 2 sin(t/2) perp(v), so the offsets match only when a = -b). That is a
  // second pin per joint - the mirrored double-pin closure.
  function hasSecondJointPin() {
    const a = SITE_VECTORS[aSiteInput.value] || SITE_VECTORS.east;
    const b = SITE_VECTORS[bSiteInput.value] || SITE_VECTORS.west;
    return a.clone().add(b).lengthSq() < 1e-9;
  }

  // Circumradius of a polygon with sides `sides` inscribed in a circle
  // (every vertex on the circle), and the central angle of each side.
  // Exact for any mix of side lengths, so the ring always closes; the
  // physical question becomes whether each joint can bend that far
  // (backlash tilt check) and whether bodies collide.
  function solveCyclicPolygon(sides) {
    const n = sides.length;
    const maxSide = Math.max(...sides);
    const maxIndex = sides.indexOf(maxSide);
    const halfMin = maxSide / 2;
    const angleSum = (R) => sides.reduce((sum, p) => sum + 2 * Math.asin(Math.min(1, p / (2 * R))), 0);
    // If even the tightest circle (longest side as a diameter) leaves the
    // arcs short of 2*pi, the circumcenter lies outside the polygon and the
    // longest side takes the major arc instead.
    const centerInside = angleSum(halfMin) >= 2 * Math.PI;
    const residual = centerInside
      ? (R) => angleSum(R) - 2 * Math.PI
      : (R) =>
          sides.reduce((sum, p, i) => (i === maxIndex ? sum : sum + Math.asin(Math.min(1, p / (2 * R)))), 0) -
          Math.asin(Math.min(1, maxSide / (2 * R)));
    let lo = halfMin;
    let hi = halfMin * 2;
    // centerInside: residual falls with R; otherwise it rises. Bracket, then bisect.
    const sign = centerInside ? 1 : -1;
    while (sign * residual(hi) > 0) hi *= 2;
    for (let iter = 0; iter < 80; iter += 1) {
      const mid = (lo + hi) / 2;
      if (sign * residual(mid) > 0) lo = mid;
      else hi = mid;
    }
    const R = (lo + hi) / 2;
    const angles = sides.map((p, i) => {
      const minor = 2 * Math.asin(Math.min(1, p / (2 * R)));
      return !centerInside && i === maxIndex ? 2 * Math.PI - minor : minor;
    });
    const closure = Math.abs(angles.reduce((a, b) => a + b, 0) - 2 * Math.PI) * R;
    return { R, angles, closure, centerInside, n };
  }

  // Build one ring: cell i sits flat on the polygon face between joint pin
  // i (its lower bSite pin, shared with cell i-1) and joint pin i+1 (its
  // upper aSite pin, shared with cell i+1). Pins are the polygon vertices
  // on a circle about the axle, so neighboring pins coincide exactly.
  // Each cell's frame: tangent along its face (vertex i -> i+1), normal
  // pointing outward, axial along the axle. The whole ring is then rotated
  // so cell 0's hub sits at polar angle 0, so rows stack in columns.
  function buildRing(n, thetaPerCell) {
    const joints = thetaPerCell.map((theta) => cellJointGeometry(theta));
    const sides = joints.map((j) => j.pitch);
    const cyclic = solveCyclicPolygon(sides);
    const vertices = [];
    let angle = 0;
    for (let i = 0; i < n; i += 1) {
      vertices.push(new THREE.Vector2(cyclic.R * Math.cos(angle), cyclic.R * Math.sin(angle)));
      angle += cyclic.angles[i];
    }
    const meanPinAxial = joints.reduce((sum, j) => sum + j.pinAxial, 0) / n;
    const centers = [];
    const tangents = [];
    const normals = [];
    const hubZ = [];
    for (let i = 0; i < n; i += 1) {
      const start = vertices[i];
      const end = vertices[(i + 1) % n];
      const tangent = end.clone().sub(start).normalize();
      const normal = new THREE.Vector2(tangent.y, -tangent.x);
      tangents.push(tangent);
      normals.push(normal);
      centers.push(start.clone().addScaledVector(tangent, -joints[i].bPinAlong));
      // Pins of a ring share one height; each hub sits below its pins by
      // its own pinAxial, re-centered so the ring's mean hub height is 0.
      hubZ.push(meanPinAxial - joints[i].pinAxial);
    }
    const phase = -Math.atan2(centers[0].y, centers[0].x);
    const origin = new THREE.Vector2(0, 0);
    [vertices, centers, tangents, normals].forEach((list) => list.forEach((v) => v.rotateAround(origin, phase)));
    const radius = centers.reduce((sum, c) => sum + c.length(), 0) / n;
    const heading = tangents.map((t) => Math.atan2(t.y, t.x));
    return {
      centers,
      tangents,
      normals,
      vertices,
      hubZ,
      pinZ: meanPinAxial,
      crossRotLower: joints.map((j) => j.crossRotLower),
      crossRotUpper: joints.map((j) => j.crossRotUpper),
      // Absolute in-plane rotation of each cross about the cell normal,
      // measured from the ring's +x axis (for the Selected Cell readout).
      bottomRot: heading.map((h, i) => h + joints[i].crossRotLower),
      topRot: heading.map((h, i) => h + joints[i].crossRotUpper),
      radialAngle: centers.map((c) => Math.atan2(c.y, c.x)),
      pitch: sides.reduce((a, b) => a + b, 0) / n,
      pinCircleRadius: cyclic.R,
      closureResidual: cyclic.closure,
      radius,
      diameter: 2 * radius,
      centroid: origin.clone(),
    };
  }

  // World position of one pad's center (on its layer's mid-plane).
  function padWorld(ring, i, layer, site) {
    const face = rotate2(SITE_VECTORS[site], layer === "upper" ? ring.crossRotUpper[i] : ring.crossRotLower[i]);
    const t = ring.tangents[i];
    const nrm = ring.normals[i];
    const offset = layer === "upper" ? PLATE_HALF_THICKNESS : -PLATE_HALF_THICKNESS;
    return new THREE.Vector3(
      ring.centers[i].x + face.x * t.x + offset * nrm.x,
      ring.centers[i].y + face.x * t.y + offset * nrm.y,
      ring.baseZ + ring.hubZ[i] + face.y
    );
  }

  // The pinned pad pairs at the joint between cell i and cell i+1:
  // [cell i's layer, cell i+1's layer] for each pin.
  function jointPinLayers() {
    return hasSecondJointPin()
      ? [
          ["upper", "lower"],
          ["lower", "upper"],
        ]
      : [["upper", "lower"]];
  }

  // --- Constraint solver: no fusing through, no discontinuous motion ---
  // The realized pose never teleports to a requested one. It walks there
  // from the last valid pose in small steps (1 deg of cell twist per step,
  // the step Ahyan's V2 samples its drive at), collision-checking each step.
  // A step that would push parts into each other, or that makes any cell
  // jump more than 4 mm or turn more than 5 deg at once (Ahyan's
  // discontinuity thresholds), is refused; the pose then settles at the
  // contact point found by bisection and is "held" there.
  const CONSTRAINT_STEP_RAD = (1 * Math.PI) / 180;
  const MAX_STEP_CENTER_JUMP_MM = 4;
  const MAX_STEP_ROTATION_JUMP_RAD = (5 * Math.PI) / 180;
  const MAX_CONSTRAINT_STEPS_PER_FRAME = 24;
  let acceptedPose = null;
  let snapPosePending = true;
  let lastSolve = null;

  // Row spacing: `spacing === null` means the rows are pinned together, so
  // the spacing is whatever makes each row's north pads meet the south pads
  // of the row above (2L*cos(theta/2) for a uniform ring - the same as the
  // spacing around the ring). A number is a set pitch with the rows not
  // joined. Each ring gets its row's base height as ring.baseZ.
  function pinnedRowPitch(lower, upper) {
    let sum = 0;
    for (let i = 0; i < lower.centers.length; i += 1) {
      const northTop = lower.hubZ[i] + rotate2(SITE_VECTORS.north, lower.crossRotUpper[i]).y;
      const southBottom = upper.hubZ[i] + rotate2(SITE_VECTORS.south, upper.crossRotLower[i]).y;
      sum += northTop - southBottom;
    }
    return sum / lower.centers.length;
  }

  function assignRowHeights(rings, spacing) {
    if (!rings.length) return;
    rings[0].baseZ = 0;
    for (let row = 1; row < rings.length; row += 1) {
      rings[row].baseZ = rings[row - 1].baseZ + (spacing === null ? pinnedRowPitch(rings[row - 1], rings[row]) : spacing);
    }
  }

  function currentRowSpacing() {
    return pinRowsInput.checked ? null : Number(axialPitchInput.value);
  }

  // `regionRows`: only check cell pairs touching these rows (null = all).
  function buildPose(n, m, spacing, grid, withCollisions, regionRows = null) {
    const rings = grid.map((row) => buildRing(n, row));
    assignRowHeights(rings, spacing);
    return {
      grid,
      rings,
      report: withCollisions ? checkCollisions(rings, n, m, spacing === null, regionRows) : null,
    };
  }

  function angleGap(a, b) {
    return Math.abs(Math.atan2(Math.sin(a - b), Math.cos(a - b)));
  }

  function stepIsContinuous(from, to) {
    for (let row = 0; row < from.rings.length; row += 1) {
      const a = from.rings[row];
      const b = to.rings[row];
      for (let i = 0; i < a.centers.length; i += 1) {
        const jump = Math.hypot(
          a.centers[i].x - b.centers[i].x,
          a.centers[i].y - b.centers[i].y,
          a.baseZ + a.hubZ[i] - b.baseZ - b.hubZ[i]
        );
        if (jump > MAX_STEP_CENTER_JUMP_MM) return false;
        if (angleGap(a.bottomRot[i], b.bottomRot[i]) > MAX_STEP_ROTATION_JUMP_RAD) return false;
        if (angleGap(a.topRot[i], b.topRot[i]) > MAX_STEP_ROTATION_JUMP_RAD) return false;
      }
    }
    return true;
  }

  function gridSignature(grid) {
    return grid.map((row) => row.map((v) => v.toFixed(7)).join(",")).join(";");
  }

  // One constrained step toward `target`: first every cell together; if
  // that is blocked, each row on its own, and within a row whose cells
  // want different things, each cell on its own - so parts that can move
  // keep moving while blocked parts hold at contact. Moves shrink from one
  // step to 1/8 step so a blocked part creeps up to its contact point.
  // Only rows a move touches are re-checked for collisions.
  const STEP_SCALES = [1, 0.5, 0.25, 0.125];

  function advanceGrid(grid, target, maxStep, onlyRow = null, onlyCell = null) {
    return grid.map((row, r) =>
      row.map((v, i) => {
        if ((onlyRow !== null && r !== onlyRow) || (onlyCell !== null && i !== onlyCell)) return v;
        const delta = target[r][i] - v;
        return Math.abs(delta) <= maxStep ? target[r][i] : v + Math.sign(delta) * maxStep;
      })
    );
  }

  function constrainedStep(n, m, spacing, start, target) {
    const pinned = spacing === null;
    let pose = start;
    let moved = false;
    let discontinuous = false;
    // Try `grid` from the current pose, checking only `regionRows` (null =
    // all). Allowed if continuous and it doesn't add overlap in that region
    // (leaving an already-overlapping pose is fine if overlap doesn't grow).
    const attempt = (grid, regionRows) => {
      const candidate = buildPose(n, m, spacing, grid, true, regionRows);
      if (!stepIsContinuous(pose, candidate)) {
        discontinuous = true;
        return false;
      }
      const before = regionRows ? checkCollisions(pose.rings, n, m, pinned, regionRows).maxPenetration : pose.report.maxPenetration;
      if (candidate.report.maxPenetration > Math.max(COLLISION_TOLERANCE_MM, before + 1e-6)) return false;
      pose = regionRows ? buildPose(n, m, spacing, grid, true) : candidate;
      moved = true;
      return true;
    };
    const tryScaled = (onlyRow, onlyCell, region) =>
      STEP_SCALES.some((scale) => attempt(advanceGrid(pose.grid, target, CONSTRAINT_STEP_RAD * scale, onlyRow, onlyCell), region));

    if (attempt(advanceGrid(pose.grid, target, CONSTRAINT_STEP_RAD), null)) return { pose, moved, discontinuous };
    for (let r = 0; r < m; r += 1) {
      if (pose.grid[r].every((v, i) => Math.abs(v - target[r][i]) < 1e-9)) continue;
      const region = [r - 1, r, r + 1].filter((x) => x >= 0 && x < m);
      if (tryScaled(r, null, region)) continue;
      // A row with one shared command holds as a unit; only rows whose
      // cells want different things get split up cell by cell.
      if (target[r].every((v) => Math.abs(v - target[r][0]) < 1e-9)) continue;
      for (let i = 0; i < n; i += 1) {
        if (Math.abs(pose.grid[r][i] - target[r][i]) >= 1e-9) tryScaled(r, i, region);
      }
    }
    return { pose, moved, discontinuous };
  }

  function realizePose(n, m, spacing, requested) {
    const enforce = constraintsEnabledInput.checked;
    const withCollisions = enforce || collisionEnabledInput.checked;
    const topoKey = [n, m, spacing, aSiteInput.value, bSiteInput.value].join("|");
    const key = `${topoKey}|${enforce}|${withCollisions}`;
    const requestSig = gridSignature(requested);
    if (
      lastSolve &&
      !snapPosePending &&
      lastSolve.key === key &&
      lastSolve.requestSig === requestSig &&
      lastSolve.status !== "moving"
    ) {
      return lastSolve;
    }

    let pose;
    let status;
    if (!enforce || snapPosePending || !acceptedPose || acceptedPose.topoKey !== topoKey || acceptedPose.key !== key) {
      pose = buildPose(n, m, spacing, requested, withCollisions);
      status = !enforce ? "off" : pose.report.clear ? "free" : "start-collides";
      snapPosePending = false;
    } else {
      let current = acceptedPose.pose;
      const reached = () => current.grid.every((row, r) => row.every((v, i) => Math.abs(v - requested[r][i]) < 1e-9));
      status = current.report.clear ? "free" : "start-collides";
      // Work in bounded slices per frame so the page stays responsive;
      // "moving" means the next frame continues.
      const deadline = performance.now() + 40;
      let steps = 0;
      while (!reached()) {
        const result = constrainedStep(n, m, spacing, current, requested);
        current = result.pose;
        if (!result.moved) {
          status = result.discontinuous ? "discontinuous" : "held";
          break;
        }
        steps += 1;
        if (steps >= MAX_CONSTRAINT_STEPS_PER_FRAME || performance.now() > deadline) {
          if (!reached()) status = "moving";
          break;
        }
      }
      if (reached()) status = current.report.clear ? "free" : "start-collides";
      pose = current;
    }
    acceptedPose = { topoKey, key, pose };
    lastSolve = { key, requestSig, status, pose };
    return lastSolve;
  }

  // Collision-free ranges of one shared effective alpha for the current
  // ring/row/site settings: sweep alpha in 0.01 steps, keep poses with no
  // overlap, and split a range wherever consecutive poses jump
  // discontinuously. Cached per setting; a uniform ring's rows are
  // identical, so two rows capture every row-to-row interaction.
  let envelopeCache = { key: "", value: null };

  function uniformAlphaEnvelope(n, m, spacing) {
    const key = [n, Math.min(m, 2), spacing, aSiteInput.value, bSiteInput.value].join("|");
    if (envelopeCache.key === key) return envelopeCache.value;
    const rows = Math.min(m, 2);
    const ranges = [];
    let open = null;
    let previous = null;
    for (let k = 0; k <= Math.round((ALPHA_MAX - ALPHA_MIN) / 0.01); k += 1) {
      const alpha = Math.round((ALPHA_MIN + k * 0.01) * 100) / 100;
      const theta = degToRad(cellThetaDeg(alpha));
      const pose = buildPose(n, rows, spacing, Array.from({ length: rows }, () => new Array(n).fill(theta)), true);
      const ok = pose.report.clear && (!previous || !open || stepIsContinuous(previous, pose));
      if (ok) {
        if (!open) open = [alpha, alpha];
        else open[1] = alpha;
      } else if (open) {
        ranges.push(open);
        open = pose.report.clear ? [alpha, alpha] : null;
      }
      previous = pose;
    }
    if (open) ranges.push(open);
    envelopeCache = { key, value: ranges };
    return ranges;
  }

  function envelopeRangeContaining(ranges, alpha) {
    return ranges.find(([lo, hi]) => alpha >= lo - 1e-9 && alpha <= hi + 1e-9) || null;
  }

  // Limit a shared effective alpha to the collision-free range the pose is
  // in now (so it never jumps to a separate range), or, before there is a
  // pose, the range nearest the command. Keeps REACH_MARGIN from each edge.
  const REACH_MARGIN = 0.02;
  function limitToReachable(alpha, ranges) {
    if (!ranges.length) return alpha;
    const realized = lastRealizedAlphas ? lastRealizedAlphas.flat() : [];
    const current = realized.length ? realized.reduce((a, b) => a + b, 0) / realized.length : null;
    const distance = (range, a) => (a < range[0] ? range[0] - a : a > range[1] ? a - range[1] : 0);
    const nearest = (a) => ranges.reduce((best, r) => (distance(r, a) < distance(best, a) ? r : best), ranges[0]);
    const range = current === null ? nearest(alpha) : nearest(current);
    const lo = range[1] - range[0] > 2 * REACH_MARGIN ? range[0] + REACH_MARGIN : (range[0] + range[1]) / 2;
    const hi = range[1] - range[0] > 2 * REACH_MARGIN ? range[1] - REACH_MARGIN : (range[0] + range[1]) / 2;
    return clamp(alpha, lo, hi);
  }

  function thetaToAlpha(theta) {
    return (radToDeg(theta) + 60) / 70;
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
  // Effective alpha -> the commanded alpha that produces it (the inverse of
  // the backlash dead zone around alpha = 1).
  function effectiveToCommand(effective, backlash, side) {
    if (effective > ALPHA_REFERENCE) return effective + backlash;
    if (effective < ALPHA_REFERENCE) return effective - backlash;
    return side === "high" ? ALPHA_REFERENCE + backlash : ALPHA_REFERENCE - backlash;
  }

  // Green segments: commands whose uniform pose is collision-free. Hollow
  // marker: where the realized pose is, on the same command scale.
  function updateReachRuler(envelope, backlash, realizedMean, enforced) {
    const range = ALPHA_MAX - ALPHA_MIN;
    const pct = (a) => clamp(((a - ALPHA_MIN) / range) * 100, 0, 100);
    reachBands.innerHTML = "";
    envelope.forEach(([lo, hi]) => {
      const commandLo = effectiveToCommand(lo, backlash, "low");
      const commandHi = effectiveToCommand(hi, backlash, "high");
      // Skip ranges no slider command can produce (e.g. pushed past the
      // slider's end by the dead-zone shift).
      if (commandHi < ALPHA_MIN || commandLo > ALPHA_MAX) return;
      const left = pct(commandLo);
      const right = pct(commandHi);
      const band = document.createElement("div");
      band.className = "reach-band";
      band.style.left = `${left}%`;
      band.style.width = `${Math.max(0.6, right - left)}%`;
      reachBands.appendChild(band);
    });
    realizedMarker.style.display = enforced ? "" : "none";
    realizedMarker.style.left = `${pct(effectiveToCommand(realizedMean, backlash, "high"))}%`;
  }

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

  function setRadialOrientation(object3d, tangent2, normal2) {
    radialScratchX.set(tangent2.x, tangent2.y, 0); // tangential (hole-aligned)
    radialScratchZ.set(normal2.x, normal2.y, 0); // radial (outward)
    radialScratchMatrix.makeBasis(radialScratchX, radialScratchY, radialScratchZ);
    object3d.quaternion.setFromRotationMatrix(radialScratchMatrix);
  }

  function heatColor(alphaValue) {
    const t = clamp((alphaValue - ALPHA_MIN) / (ALPHA_MAX - ALPHA_MIN), 0, 1);
    if (t < 0.5) heatmapColorScratch.lerpColors(HEATMAP_LOW, HEATMAP_MID, t / 0.5);
    else heatmapColorScratch.lerpColors(HEATMAP_MID, HEATMAP_HIGH, (t - 0.5) / 0.5);
    return heatmapColorScratch;
  }

  // Diameter vs alpha is NOT monotonic over the slider's range: the pin
  // chord 2L*cos(theta/2) is longest at zero twist (theta = 0, alpha ~0.86)
  // and shrinks as the twist grows either way, so the ring is widest there
  // and narrower toward both ends. That rules out a bisection search or a
  // closed-form inverse of theta=70*alpha-60 over the full range. Instead it
  // exhaustively evaluates every alpha the slider can reach (a uniform
  // ring, all cells at the same theta) and keeps whichever is closest to
  // the target - simple, robust to non-monotonicity, and cheap since a
  // whole ring build is just O(n).
  // `allowedRanges` (effective-alpha ranges) restricts the search to poses
  // the constraint solver can actually reach; null means unrestricted.
  function fitAlphaToDiameter(targetDiameter, n, backlash, allowedRanges = null) {
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
      if (allowedRanges && !envelopeRangeContaining(allowedRanges, Math.round(effectiveAlpha * 100) / 100)) continue;
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

  const selectionMarker = faceRingMarker(6.2, new THREE.MeshStandardMaterial({ color: 0xffd54a, emissive: 0x7a5c00, roughness: 0.3 }));
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
  let lastPinAlignment = null;
  let lastPinCount = 0;
  let lastCollisionReport = null;
  let lastRealizedAlphas = null;
  let lastConstraintStatus = "off";
  let lastEnvelope = [];

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
    const spacing = currentRowSpacing();

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

    // With constraints on, the shared drive can't go past the reachable
    // range - like Ahyan's V2 clamping the drive to the nearest feasible
    // value - with a small margin kept from contact so the rest of the
    // structure isn't jammed against its stops.
    const envelope = uniformAlphaEnvelope(n, m, spacing);
    lastEnvelope = envelope;
    const drivenBaseline = constraintsEnabledInput.checked ? limitToReachable(baselineAlpha, envelope) : baselineAlpha;
    const driveLimited = Math.abs(drivenBaseline - baselineAlpha) > 1e-9;

    const cellAlphas = computeCellAlphas(n, m, drivenBaseline, backlash);
    lastCellAlphas = cellAlphas;
    lastCellRolesSnapshot = cellRoles;
    const requestedTheta = cellAlphas.map((row) => row.map((a) => degToRad(cellThetaDeg(a))));
    const solve = realizePose(n, m, spacing, requestedTheta);
    const rings = solve.pose.rings;
    const thetaPerRow = solve.pose.grid;
    const realizedAlphas = thetaPerRow.map((row) => row.map(thetaToAlpha));
    lastRealizedAlphas = realizedAlphas;
    lastRingsSnapshot = rings;
    const cellFramesByRow = [];
    for (let row = 0; row < m; row += 1) cellFramesByRow.push(new Array(n));

    const showMeasurements = showMeasurementsInput.checked && !isolateActive;
    measurementLine.visible = showMeasurements;
    measurementTickA.visible = showMeasurements;
    measurementTickB.visible = showMeasurements;
    if (showMeasurements) updateMeasurementLine(rings[0], 0);

    for (let row = 0; row < m; row += 1) {
      for (let i = 0; i < n; i += 1) {
        const cell = cellPool[row * n + i];
        const center = rings[row].centers[i];
        const z = rings[row].baseZ + rings[row].hubZ[i];
        cell.group.position.set(center.x, center.y, z);
        const frame = { tangent: rings[row].tangents[i], normal: rings[row].normals[i] };
        cellFramesByRow[row][i] = frame;
        setRadialOrientation(cell.group, frame.tangent, frame.normal);
        // Cross twists are measured in the cell's own face frame (group
        // local x = tangent, y = axle, z = outward normal): lower at
        // -theta/2, upper at +theta/2 (see cellJointGeometry).
        cell.bottom.rotation.z = rings[row].crossRotLower[i];
        cell.top.rotation.z = rings[row].crossRotUpper[i];
        const role = cellRoles[row][i].role;
        if (role === "free") {
          cell.roleMarker.visible = false;
        } else {
          cell.roleMarker.visible = true;
          cell.roleMarker.material = role === "actuator" ? sharedMaterials.actuatorMarker : sharedMaterials.lockedMarker;
          placeFaceMarker(cell.roleMarker, center, z, frame.normal);
        }
        if (multiSelected.has(cellKey(row, i))) {
          cell.batchMarker.visible = true;
          placeFaceMarker(cell.batchMarker, center, z, frame.normal);
        } else {
          cell.batchMarker.visible = false;
        }
        if (heatmapEnabledInput.checked) {
          const color = heatColor(realizedAlphas[row][i]);
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
    const meanRowPitch = m > 1 ? rings[m - 1].baseZ / (m - 1) : 0;
    axialPitchInput.disabled = spacing === null;
    axialPitchOut.textContent =
      spacing === null ? (m > 1 ? `${meanRowPitch.toFixed(1)} mm (set by the pins)` : "set by the pins") : `${spacing.toFixed(0)} mm`;
    pitchMetric.textContent = `${rings[0].pitch.toFixed(1)} mm`;
    turnMetric.textContent = `${radToDeg((2 * Math.PI) / n).toFixed(1)} deg`;
    diameterMetric.textContent = `${rings[0].diameter.toFixed(1)} mm`;
    radiusMetric.textContent = `${rings[0].radius.toFixed(1)} mm`;
    closureMetric.textContent = `${rings[0].closureResidual.toFixed(4)} mm`;
    totalCellsMetric.textContent = String(n * m);
    heightMetric.textContent = `${rings[m - 1].baseZ.toFixed(1)} mm`;

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

    const report = solve.pose.report;
    lastCollisionReport = report;
    updateContactMarkers(collisionEnabledInput.checked && report && !isolateActive ? report.contacts : []);
    if (collisionEnabledInput.checked && report) {
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

    const CONSTRAINT_LABELS = {
      off: "not enforced",
      free: "free - pose reached",
      moving: "moving toward the command",
      held: "held at contact (would fuse through)",
      discontinuous: "held (step would jump discontinuously)",
      "start-collides": "start pose overlaps - can only move out",
      limited: "drive limited to the reachable range",
    };
    const status = driveLimited && (solve.status === "free" || solve.status === "moving") ? "limited" : solve.status;
    lastConstraintStatus = status;
    constraintStateMetric.textContent = CONSTRAINT_LABELS[status] || status;
    constraintStateMetric.classList.toggle("status-ok", status === "free" || status === "moving");
    constraintStateMetric.classList.toggle(
      "status-blocked",
      status === "held" || status === "discontinuous" || status === "start-collides" || status === "limited"
    );
    let meanRealized = 0;
    let maxLag = 0;
    realizedAlphas.forEach((row, r) =>
      row.forEach((a, i) => {
        meanRealized += a / (n * m);
        maxLag = Math.max(maxLag, Math.abs(a - thetaToAlpha(requestedTheta[r][i])));
      })
    );
    realizedAlphaMetric.textContent =
      maxLag > 0.005 ? `${meanRealized.toFixed(2)} mean (lags command by up to ${maxLag.toFixed(2)})` : `${meanRealized.toFixed(2)} mean`;
    feasibleRangeMetric.textContent = envelope.length
      ? envelope.map(([lo, hi]) => (lo === hi ? lo.toFixed(2) : `${lo.toFixed(2)} - ${hi.toFixed(2)}`)).join(", ")
      : "none";
    updateReachRuler(envelope, backlash, meanRealized, constraintsEnabledInput.checked);
    const stoppedAtContact = status === "held" || status === "discontinuous";
    driveRealizedMetric.textContent =
      status === "limited"
        ? `${meanRealized.toFixed(2)} - command ${alphaCommand.toFixed(2)} is out of reach`
        : stoppedAtContact
          ? `${meanRealized.toFixed(2)} - stopped at contact`
          : status === "moving"
            ? `${meanRealized.toFixed(2)} - moving`
            : meanRealized.toFixed(2);
    driveRealizedMetric.classList.toggle("status-blocked", stoppedAtContact || status === "limited");

    const pinAlignment = computePinAlignment(n, m, rings, cellFramesByRow);
    lastPinAlignment = pinAlignment;
    circumferentialPinMetric.textContent = `${pinAlignment.maxCircumferential.toFixed(3)} mm`;
    axialPinMetric.textContent = `${pinAlignment.maxAxial.toFixed(3)} mm`;
    jointTiltLimitMetric.textContent = `${radToDeg(pinAlignment.tiltLimit).toFixed(2)} deg`;
    const setBendStatus = (element, bend, ok) => {
      element.textContent = `${radToDeg(bend).toFixed(1)} deg - ${ok ? "within backlash" : "exceeds backlash"}`;
      element.classList.toggle("status-ok", ok);
      element.classList.toggle("status-blocked", !ok);
    };
    setBendStatus(circumferentialBendMetric, pinAlignment.maxCircumferentialBend, pinAlignment.circumferentialOk);
    setBendStatus(axialBendMetric, pinAlignment.maxAxialBend, pinAlignment.axialOk);
    minCellsClosureMetric.textContent = Number.isFinite(pinAlignment.minCellsForBacklashClosure)
      ? String(pinAlignment.minCellsForBacklashClosure)
      : "none (pin fills hole)";

    const pinClearance = pinRadialClearanceMm();
    const deltaPhi = pinDeltaPhiRad();
    pinDiameterOut.textContent = `${pinDiameterMm().toFixed(2)} mm`;
    pinClearanceMetric.textContent = `${pinClearance.toFixed(3)} mm`;
    pinBOverLMetric.textContent = (pinClearance / CAD.siteRadiusMm).toFixed(4);
    pinDeltaPhiMetric.textContent = `${radToDeg(deltaPhi).toFixed(2)} deg`;
    // theta = 70*alpha - 60, so one degree of cell rotation is 1/70 alpha.
    pinAlphaDeadzoneMetric.textContent = `+/- ${(radToDeg(deltaPhi) / 70).toFixed(4)}`;
    updatePins(n, m, rings, cellFramesByRow, spacing === null);

    if (selected && selected.row < m && selected.i < n) {
      const row = selected.row;
      const i = selected.i;
      const center = rings[row].centers[i];
      const z = rings[row].baseZ + rings[row].hubZ[i];
      selectionMarker.visible = true;
      placeFaceMarker(selectionMarker, center, z, rings[row].normals[i]);
      selectedWorld.set(center.x, center.y, z);
      selectedCellStatus.textContent = `row ${row}, cell ${i}`;
      selectedCellStatus.classList.add("status-ok");
      selectedIndexMetric.textContent = `${row}, ${i}`;
      selectedCenterMetric.textContent = `(${center.x.toFixed(1)}, ${center.y.toFixed(1)}, ${z.toFixed(1)})`;
      selectedBottomRotMetric.textContent = `${radToDeg(rings[row].bottomRot[i]).toFixed(1)} deg`;
      selectedTopRotMetric.textContent = `${radToDeg(rings[row].topRot[i]).toFixed(1)} deg`;
      const commandedCellAlpha = cellAlphas[row][i];
      const realizedCellAlpha = realizedAlphas[row][i];
      selectedCellAlphaMetric.textContent =
        Math.abs(realizedCellAlpha - commandedCellAlpha) > 0.005
          ? `${realizedCellAlpha.toFixed(2)} (commanded ${commandedCellAlpha.toFixed(2)})`
          : realizedCellAlpha.toFixed(2);
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
    [alphaInput, backlashInput, ringCountInput, rowCountInput, axialPitchInput, targetDiameterInput, pinDiameterInput].forEach((input) => {
      input.value = input.defaultValue;
    });
    [aSiteInput, bSiteInput].forEach((select) => {
      select.value = Array.from(select.options).find((option) => option.defaultSelected).value;
    });
    [animateInput, collisionEnabledInput, heatmapEnabledInput, showPinsInput, constraintsEnabledInput, pinRowsInput].forEach((input) => {
      input.checked = input.defaultChecked;
    });
    clearAllRoles();
    // A reset is a fresh setup, not a motion - jump straight to the new pose.
    snapPosePending = true;
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
  // (rows solve independently) - a one-click non-uniform wheel profile.
  // `rowWidthFn(row, m)` gives each row's wanted width in [0, 1] (1 = widest).
  // Alphas are picked from the collision-free window the pose is currently
  // in (or the widest window), so the constraint solver can actually reach
  // them. Diameter falls as alpha rises in that window - more twist
  // shortens the pin-to-pin span 2L*cos(theta/2) - so the widest row gets
  // the window's lowest alpha. Without a window, fall back to 1.25-1.75.
  function applyRowWidthPreset(rowWidthFn) {
    const margin = 0.02;
    let lo = 1.25;
    let hi = 1.75;
    if (lastEnvelope.length) {
      const flat = lastRealizedAlphas ? lastRealizedAlphas.flat() : [];
      const current = flat.length ? Math.round((flat.reduce((a, b) => a + b, 0) / flat.length) * 100) / 100 : null;
      const range =
        (current !== null && envelopeRangeContaining(lastEnvelope, current)) ||
        lastEnvelope.reduce((best, r) => (r[1] - r[0] > best[1] - best[0] ? r : best), lastEnvelope[0]);
      if (range[1] - range[0] > 2 * margin) {
        lo = range[0] + margin;
        hi = range[1] - margin;
      }
    }
    const m = cellRoles.length;
    for (let row = 0; row < m; row += 1) {
      const n = cellRoles[row].length;
      const width = clamp(rowWidthFn(row, m), 0, 1);
      const alpha = Math.round(clamp(hi - width * (hi - lo), ALPHA_MIN, ALPHA_MAX) * 100) / 100;
      for (let i = 0; i < n; i += 1) {
        cellRoles[row][i] = { role: "actuator", alpha };
      }
    }
    updateMechanism();
  }

  // 1 at the middle row(s), 0 at the ends.
  function middleWeight(row, m) {
    const mid = (m - 1) / 2;
    return mid === 0 ? 1 : 1 - Math.abs(row - mid) / mid;
  }

  presetBarrelBtn.addEventListener("click", () => applyRowWidthPreset((row, m) => middleWeight(row, m)));
  presetConeBtn.addEventListener("click", () => applyRowWidthPreset((row, m) => (m <= 1 ? 0.5 : row / (m - 1))));
  presetSaddleBtn.addEventListener("click", () => applyRowWidthPreset((row, m) => 1 - middleWeight(row, m)));

  function applyDiameterFit() {
    const n = clamp(Math.round(Number(ringCountInput.value)), RING_COUNT_MIN, RING_COUNT_MAX);
    const target = Number(targetDiameterInput.value);
    const backlash = Number(backlashInput.value);
    // With constraints on, only fit within the collision-free range the
    // current pose is in (the pose can't cross a collision to reach another).
    let allowedRanges = null;
    if (constraintsEnabledInput.checked && lastEnvelope.length) {
      const flat = lastRealizedAlphas ? lastRealizedAlphas.flat() : [];
      const current = flat.length ? Math.round((flat.reduce((a, b) => a + b, 0) / flat.length) * 100) / 100 : null;
      const range = current === null ? null : envelopeRangeContaining(lastEnvelope, current);
      // Same margin from contact the drive limit keeps, so the fitted alpha
      // is exactly what gets applied.
      allowedRanges = (range ? [range] : lastEnvelope).map(([lo, hi]) =>
        hi - lo > 2 * REACH_MARGIN ? [lo + REACH_MARGIN, hi - REACH_MARGIN] : [(lo + hi) / 2, (lo + hi) / 2]
      );
    }
    const fit = fitAlphaToDiameter(target, n, backlash, allowedRanges);
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
    pinDiameterInput,
    showPinsInput,
    constraintsEnabledInput,
    pinRowsInput,
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
      pinDiameter: Number(pinDiameterInput.value),
      ringCount: Number(ringCountInput.value),
      rowCount: Number(rowCountInput.value),
      axialPitch: Number(axialPitchInput.value),
      pinRows: pinRowsInput.checked,
      aSite: aSiteInput.value,
      bSite: bSiteInput.value,
      cellRoles,
    };
  }

  function applyStateSnapshot(state) {
    if (!state || typeof state !== "object") return;
    if (Number.isFinite(state.alpha)) alphaInput.value = String(clamp(state.alpha, ALPHA_MIN, ALPHA_MAX));
    if (Number.isFinite(state.backlash)) backlashInput.value = String(clamp(state.backlash, Number(backlashInput.min), Number(backlashInput.max)));
    if (Number.isFinite(state.pinDiameter)) {
      pinDiameterInput.value = String(clamp(state.pinDiameter, Number(pinDiameterInput.min), Number(pinDiameterInput.max)));
    }
    if (Number.isFinite(state.ringCount)) ringCountInput.value = String(clamp(Math.round(state.ringCount), RING_COUNT_MIN, RING_COUNT_MAX));
    if (Number.isFinite(state.rowCount)) rowCountInput.value = String(clamp(Math.round(state.rowCount), ROW_COUNT_MIN, ROW_COUNT_MAX));
    if (Number.isFinite(state.axialPitch)) axialPitchInput.value = String(clamp(state.axialPitch, Number(axialPitchInput.min), Number(axialPitchInput.max)));
    // Files saved before rows could be pinned used a set pitch.
    pinRowsInput.checked = state.pinRows === true;
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
        // Loading a file sets up a pose; it isn't a motion to walk through.
        snapPosePending = true;
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
    // World position of one pad (site) on one layer of a rendered cell, read
    // from the actual meshes - lets tests check that joined pads coincide.
    getSiteWorld: (row, i, layer, site) => {
      const cell = cellPool[row * poolN + i];
      if (!cell || !SITE_VECTORS[site]) return null;
      const cross = layer === "upper" ? cell.top : cell.bottom;
      cell.group.updateMatrixWorld(true);
      const local = SITE_VECTORS[site];
      const world = cross.localToWorld(new THREE.Vector3(local.x, local.y, 0));
      return { x: world.x, y: world.y, z: world.z };
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
    getPinAlignment: () => lastPinAlignment,
    getRealizedAlphas: () => lastRealizedAlphas,
    getConstraintState: () => lastConstraintStatus,
    getEnvelope: () => lastEnvelope,
    getCollisionReport: () =>
      lastCollisionReport && {
        maxPenetration: lastCollisionReport.maxPenetration,
        clear: lastCollisionReport.clear,
        pairsChecked: lastCollisionReport.pairsChecked,
        contactCount: lastCollisionReport.contacts.length,
        deepest: lastCollisionReport.contacts[0]
          ? {
              cellA: lastCollisionReport.contacts[0].cellA,
              cellB: lastCollisionReport.contacts[0].cellB,
              partA: lastCollisionReport.contacts[0].partA,
              partB: lastCollisionReport.contacts[0].partB,
              depth: lastCollisionReport.contacts[0].depth,
            }
          : null,
      },
    // Collision-check rows at the given effective alphas (one per row),
    // without touching the live pose.
    probeRowAlphas: (rowAlphas) => {
      const n = clamp(Math.round(Number(ringCountInput.value)), RING_COUNT_MIN, RING_COUNT_MAX);
      const grid = rowAlphas.map((a) => new Array(n).fill(degToRad(cellThetaDeg(a))));
      const pose = buildPose(n, grid.length, currentRowSpacing(), grid, true);
      const c = pose.report.contacts[0];
      return { maxPenetration: pose.report.maxPenetration, deepest: c ? `${c.cellA}:${c.partA} x ${c.cellB}:${c.partB}` : null };
    },
    // Collision-check a uniform ring at a given effective alpha without
    // touching the live pose (used to map the collision-free range).
    probeUniformAlpha: (alphaEffective) => {
      const n = clamp(Math.round(Number(ringCountInput.value)), RING_COUNT_MIN, RING_COUNT_MAX);
      const m = Math.min(2, clamp(Math.round(Number(rowCountInput.value)), ROW_COUNT_MIN, ROW_COUNT_MAX));
      const theta = degToRad(cellThetaDeg(alphaEffective));
      const report = buildPose(n, m, currentRowSpacing(), Array.from({ length: m }, () => new Array(n).fill(theta)), true).report;
      return {
        maxPenetration: report.maxPenetration,
        clear: report.clear,
        contacts: report.contacts.slice(0, 8).map((c) => `${c.cellA}:${c.partA} x ${c.cellB}:${c.partB} ${c.depth.toFixed(2)}`),
      };
    },
    getPinInfo: () => ({
      visibleCount: lastPinCount,
      renderedRadius: pinMeshes.find((pin) => pin.visible)?.scale.x ?? null,
      diameterMm: pinDiameterMm(),
      deltaPhiDeg: radToDeg(pinDeltaPhiRad()),
    }),
    // Client-space point over the visible cell hub nearest the camera, so
    // tests can click a real cell instead of guessing (the ring is hollow
    // and centered on the axle, so the canvas center looks into empty space).
    getNearestCellClientPoint: (onlyRow) => {
      const rect = renderer.domElement.getBoundingClientRect();
      let best = null;
      cellPool.forEach((cell, index) => {
        if (!cell.group.visible) return;
        if (onlyRow !== undefined && Math.floor(index / poolN) !== onlyRow) return;
        const world = cell.group.position.clone();
        const dist = world.distanceTo(camera.position);
        const ndc = world.clone().project(camera);
        if (Math.abs(ndc.x) > 1 || Math.abs(ndc.y) > 1 || ndc.z > 1) return;
        if (!best || dist < best.dist) {
          best = {
            dist,
            x: rect.left + ((ndc.x + 1) / 2) * rect.width,
            y: rect.top + ((1 - ndc.y) / 2) * rect.height,
            row: Math.floor(index / poolN),
            i: index % poolN,
          };
        }
      });
      return best;
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
