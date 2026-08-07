(function () {
  "use strict";

  const RAD = (window.RAD = window.RAD || {});

  function fmt(value) {
    const number = Number(value || 0);
    if (!Number.isFinite(number) || Math.abs(number) < 1e-12) return "0";
    return Number(number.toFixed(9)).toString();
  }

  function squareVertices(center, side, thetaDegrees, z) {
    const half = side / 2;
    const theta = (thetaDegrees * Math.PI) / 180;
    const cos = Math.cos(theta);
    const sin = Math.sin(theta);
    return [
      [-half, -half],
      [half, -half],
      [half, half],
      [-half, half],
    ].map(([x, y]) => ({
      x: center.x + x * cos - y * sin,
      y: center.y + x * sin + y * cos,
      z,
    }));
  }

  function extrudedPolygon(vertices, thickness) {
    if (!Array.isArray(vertices) || vertices.length < 3) {
      throw new Error("plate polygon needs at least three vertices");
    }
    const lower = vertices.map((v) => ({ x: v.x, y: v.y, z: v.z - thickness / 2 }));
    const upper = vertices.map((v) => ({ x: v.x, y: v.y, z: v.z + thickness / 2 }));
    const n = vertices.length;
    const faces = [];
    for (let i = 1; i < n - 1; i += 1) {
      faces.push([0, i + 1, i]);
      faces.push([n, n + i, n + i + 1]);
    }
    for (let i = 0; i < n; i += 1) {
      const j = (i + 1) % n;
      faces.push([i, j, n + j]);
      faces.push([i, n + j, n + i]);
    }
    return { vertices: lower.concat(upper), faces };
  }

  function cylinder(center, radius, height, segments) {
    if (segments < 6) throw new Error("pinSegments must be at least 6");
    const vertices = [];
    for (let layer = 0; layer < 2; layer += 1) {
      const z = center.z + (layer === 0 ? -height / 2 : height / 2);
      for (let i = 0; i < segments; i += 1) {
        const angle = (i / segments) * Math.PI * 2;
        vertices.push({
          x: center.x + Math.cos(angle) * radius,
          y: center.y + Math.sin(angle) * radius,
          z,
        });
      }
    }
    vertices.push({ x: center.x, y: center.y, z: center.z - height / 2 });
    vertices.push({ x: center.x, y: center.y, z: center.z + height / 2 });
    const bottomCenter = segments * 2;
    const topCenter = bottomCenter + 1;
    const faces = [];
    for (let i = 0; i < segments; i += 1) {
      const j = (i + 1) % segments;
      faces.push([i, j, segments + j]);
      faces.push([i, segments + j, segments + i]);
      faces.push([bottomCenter, j, i]);
      faces.push([topCenter, segments + i, segments + j]);
    }
    return { vertices, faces };
  }

  function midpoint(a, b) {
    return { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, z: (a.z + b.z) / 2 };
  }

  function sideMidpoint(cell, side) {
    const pairs = {
      bottom: [0, 1],
      right: [1, 2],
      top: [2, 3],
      left: [3, 0],
    };
    const [i, j] = pairs[side];
    return midpoint(cell.outer[i], cell.outer[j]);
  }

  function barBetween(start, end, width) {
    const dx = end.x - start.x;
    const dy = end.y - start.y;
    const dz = end.z - start.z;
    const length = Math.hypot(dx, dy, dz);
    if (length < 1e-9) throw new Error("connector endpoints must be distinct");
    const ux = dx / length;
    const uy = dy / length;
    let vx = -uy;
    let vy = ux;
    const vl = Math.hypot(vx, vy) || 1;
    vx = (vx / vl) * (width / 2);
    vy = (vy / vl) * (width / 2);
    const wz = width / 2;
    const vertices = [
      { x: start.x - vx, y: start.y - vy, z: start.z - wz },
      { x: start.x + vx, y: start.y + vy, z: start.z - wz },
      { x: start.x + vx, y: start.y + vy, z: start.z + wz },
      { x: start.x - vx, y: start.y - vy, z: start.z + wz },
      { x: end.x - vx, y: end.y - vy, z: end.z - wz },
      { x: end.x + vx, y: end.y + vy, z: end.z - wz },
      { x: end.x + vx, y: end.y + vy, z: end.z + wz },
      { x: end.x - vx, y: end.y - vy, z: end.z + wz },
    ];
    return {
      vertices,
      faces: [
        [0, 1, 2],
        [0, 2, 3],
        [4, 6, 5],
        [4, 7, 6],
        [0, 4, 5],
        [0, 5, 1],
        [1, 5, 6],
        [1, 6, 2],
        [2, 6, 7],
        [2, 7, 3],
        [3, 7, 4],
        [3, 4, 0],
      ],
    };
  }

  function component(name, kind, cell, data) {
    return { name, kind, cell, vertices: data.vertices, faces: data.faces };
  }

  function buildPaperRadMesh(state, options = {}) {
    const sim = options.sim || RAD.simulate(state);
    const cellSize = Number(state.grid.cellSize || 1);
    const backlash = Math.max(0, Number(state.grid.backlash || 0));
    const pinRadius = Math.max(
      0.02 * cellSize,
      Number(state.grid.pinRadius ?? 0.18) * cellSize * 0.32
    );
    const plateThickness = Math.max(1e-4, Number(options.plateThickness ?? 0.035));
    const pinHeight = Math.max(1e-4, Number(options.pinHeight ?? 0.075));
    const connectorWidth = Math.max(1e-4, Number(options.connectorWidth ?? 0.045));
    const pinSegments = Math.max(6, Math.floor(Number(options.pinSegments ?? 12)));
    const includePins = options.includePins !== false;
    const includeConnectors = options.includeConnectors !== false;
    const cells = [];
    const components = [];

    for (let r = 0; r < state.grid.rows; r += 1) {
      cells[r] = [];
      for (let c = 0; c < state.grid.cols; c += 1) {
        const center = sim.centers[r][c];
        const alpha = sim.alpha[r][c];
        const theta = sim.theta[r][c];
        const activeSide = 0.62 * cellSize * Math.sqrt(Math.max(0.001, alpha));
        const outerSide = activeSide + backlash * cellSize * 0.5;
        const innerSide = Math.max(0.18 * cellSize, activeSide * 0.58);
        const outer = squareVertices(center, outerSide, 0, center.z);
        const inner = squareVertices(center, innerSide, theta, center.z + 0.08 * cellSize);
        const record = { outer, inner, center, alpha, theta };
        cells[r][c] = record;
        const prefix = `cell_${r}_${c}`;
        components.push(
          component(
            `${prefix}_outer_plate`,
            "outer_plate",
            [r, c],
            extrudedPolygon(outer, plateThickness)
          )
        );
        components.push(
          component(
            `${prefix}_inner_plate`,
            "inner_plate",
            [r, c],
            extrudedPolygon(inner, plateThickness)
          )
        );
        if (includePins) {
          outer.concat(inner).forEach((joint, index) => {
            const part = index < 4 ? "outer" : "inner";
            components.push(
              component(
                `${prefix}_${part}_pin_${index % 4}`,
                "pin",
                [r, c],
                cylinder(joint, pinRadius, pinHeight, pinSegments)
              )
            );
          });
        }
      }
    }

    if (includeConnectors) {
      for (let r = 0; r < state.grid.rows; r += 1) {
        for (let c = 0; c < state.grid.cols; c += 1) {
          if (c + 1 < state.grid.cols) {
            components.push(
              component(
                `connector_${r}_${c}_to_${r}_${c + 1}`,
                "connector",
                null,
                barBetween(
                  sideMidpoint(cells[r][c], "right"),
                  sideMidpoint(cells[r][c + 1], "left"),
                  connectorWidth
                )
              )
            );
          }
          if (r + 1 < state.grid.rows) {
            components.push(
              component(
                `connector_${r}_${c}_to_${r + 1}_${c}`,
                "connector",
                null,
                barBetween(
                  sideMidpoint(cells[r][c], "top"),
                  sideMidpoint(cells[r + 1][c], "bottom"),
                  connectorWidth
                )
              )
            );
          }
        }
      }
    }

    const vertexCount = components.reduce((sum, item) => sum + item.vertices.length, 0);
    const faceCount = components.reduce((sum, item) => sum + item.faces.length, 0);
    return { components, source: { rows: state.grid.rows, cols: state.grid.cols }, vertexCount, faceCount };
  }

  function exportPaperRadMeshObj(stateOrMesh, options = {}) {
    const mesh = stateOrMesh?.components ? stateOrMesh : buildPaperRadMesh(stateOrMesh, options);
    const lines = [
      "# RAD paper lattice normalized browser OBJ export",
      `# components ${mesh.components.length}`,
      `# vertices ${mesh.vertexCount}`,
      `# faces ${mesh.faceCount}`,
    ];
    let offset = 1;
    for (const part of mesh.components) {
      lines.push(`o ${part.name}`);
      lines.push(`# kind ${part.kind}`);
      if (part.cell) lines.push(`# cell ${part.cell[0]} ${part.cell[1]}`);
      for (const vertex of part.vertices) {
        lines.push(`v ${fmt(vertex.x)} ${fmt(vertex.y)} ${fmt(vertex.z)}`);
      }
      for (const face of part.faces) {
        lines.push(`f ${face[0] + offset} ${face[1] + offset} ${face[2] + offset}`);
      }
      offset += part.vertices.length;
    }
    return `${lines.join("\n")}\n`;
  }

  RAD.buildPaperRadMesh = buildPaperRadMesh;
  RAD.exportPaperRadMeshObj = exportPaperRadMeshObj;
})();
