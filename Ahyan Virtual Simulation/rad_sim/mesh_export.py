from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .cell_geometry import (
    PaperRADLatticeGeometry,
    build_paper_rad_lattice_geometry,
)
from .models import LatticeConfig, LatticeState, SimulationResult


@dataclass(frozen=True)
class MeshComponent:
    name: str
    kind: str
    vertices: np.ndarray
    faces: np.ndarray
    cell: tuple[int, int] | None = None

    @property
    def vertex_count(self) -> int:
        return int(self.vertices.shape[0])

    @property
    def face_count(self) -> int:
        return int(self.faces.shape[0])


@dataclass(frozen=True)
class PaperRADMesh:
    components: tuple[MeshComponent, ...]
    source: PaperRADLatticeGeometry

    @property
    def vertex_count(self) -> int:
        return sum(component.vertex_count for component in self.components)

    @property
    def face_count(self) -> int:
        return sum(component.face_count for component in self.components)

    @property
    def component_count(self) -> int:
        return len(self.components)

    @property
    def bounds(self) -> tuple[np.ndarray, np.ndarray]:
        if not self.components:
            zeros = np.zeros(3, dtype=float)
            return zeros, zeros
        vertices = np.vstack([component.vertices for component in self.components])
        return vertices.min(axis=0), vertices.max(axis=0)

    def kind_counts(self) -> dict[str, int]:
        return dict(Counter(component.kind for component in self.components))


def _validate_positive(name: str, value: float) -> float:
    value = float(value)
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _extruded_polygon(
    vertices: np.ndarray, thickness: float
) -> tuple[np.ndarray, np.ndarray]:
    if vertices.ndim != 2 or vertices.shape[1] != 3 or vertices.shape[0] < 3:
        raise ValueError("vertices must have shape (n, 3) with n >= 3")
    lower = vertices.copy()
    upper = vertices.copy()
    lower[:, 2] -= thickness / 2.0
    upper[:, 2] += thickness / 2.0
    n = vertices.shape[0]
    faces: list[tuple[int, int, int]] = []
    for i in range(1, n - 1):
        faces.append((0, i + 1, i))
        faces.append((n, n + i, n + i + 1))
    for i in range(n):
        j = (i + 1) % n
        faces.append((i, j, n + j))
        faces.append((i, n + j, n + i))
    return np.vstack([lower, upper]), np.asarray(faces, dtype=int)


def _cylinder(
    center: np.ndarray,
    radius: float,
    height: float,
    segments: int,
) -> tuple[np.ndarray, np.ndarray]:
    if segments < 6:
        raise ValueError("pin_segments must be at least 6")
    angles = np.linspace(0.0, 2.0 * np.pi, segments, endpoint=False)
    ring = np.column_stack([np.cos(angles) * radius, np.sin(angles) * radius])
    bottom = np.column_stack(
        [
            center[0] + ring[:, 0],
            center[1] + ring[:, 1],
            np.full(segments, center[2] - height / 2.0),
        ]
    )
    top = bottom.copy()
    top[:, 2] = center[2] + height / 2.0
    bottom_center = np.array([[center[0], center[1], center[2] - height / 2.0]])
    top_center = np.array([[center[0], center[1], center[2] + height / 2.0]])
    faces: list[tuple[int, int, int]] = []
    bottom_center_index = segments * 2
    top_center_index = bottom_center_index + 1
    for i in range(segments):
        j = (i + 1) % segments
        faces.append((i, j, segments + j))
        faces.append((i, segments + j, segments + i))
        faces.append((bottom_center_index, j, i))
        faces.append((top_center_index, segments + i, segments + j))
    vertices = np.vstack([bottom, top, bottom_center, top_center])
    return vertices, np.asarray(faces, dtype=int)


def _bar_between(
    start: np.ndarray, end: np.ndarray, width: float
) -> tuple[np.ndarray, np.ndarray]:
    axis = end - start
    length = np.linalg.norm(axis)
    if length < 1e-12:
        raise ValueError("connector endpoints must be distinct")
    u = axis / length
    reference = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(u, reference)) > 0.95:
        reference = np.array([0.0, 1.0, 0.0])
    v = np.cross(u, reference)
    v = v / np.linalg.norm(v) * (width / 2.0)
    w = np.cross(u, v)
    w = w / np.linalg.norm(w) * (width / 2.0)
    vertices = np.array(
        [
            start - v - w,
            start + v - w,
            start + v + w,
            start - v + w,
            end - v - w,
            end + v - w,
            end + v + w,
            end - v + w,
        ],
        dtype=float,
    )
    faces = np.asarray(
        [
            (0, 1, 2),
            (0, 2, 3),
            (4, 6, 5),
            (4, 7, 6),
            (0, 4, 5),
            (0, 5, 1),
            (1, 5, 6),
            (1, 6, 2),
            (2, 6, 7),
            (2, 7, 3),
            (3, 7, 4),
            (3, 4, 0),
        ],
        dtype=int,
    )
    return vertices, faces


def _component(
    name: str,
    kind: str,
    cell: tuple[int, int] | None,
    data: tuple[np.ndarray, np.ndarray],
) -> MeshComponent:
    vertices, faces = data
    return MeshComponent(name=name, kind=kind, cell=cell, vertices=vertices, faces=faces)


def build_paper_rad_lattice_mesh(
    config: LatticeConfig,
    state_or_result: LatticeState | SimulationResult | None = None,
    *,
    plate_thickness: float = 0.035,
    pin_height: float = 0.075,
    pin_segments: int = 12,
    connector_width: float = 0.045,
    include_pins: bool = True,
    include_connectors: bool = True,
) -> PaperRADMesh:
    """Build a CAD-style normalized triangle mesh for the paper RAD lattice.

    This is an inspectable mesh approximation of the paper-grounded topology,
    not a manufacturing model. Thicknesses and pin visual radii are configurable
    because the extracted RAD papers do not provide exact CAD part dimensions.
    """
    plate_thickness = _validate_positive("plate_thickness", plate_thickness)
    pin_height = _validate_positive("pin_height", pin_height)
    connector_width = _validate_positive("connector_width", connector_width)
    lattice = build_paper_rad_lattice_geometry(config, state_or_result)
    components: list[MeshComponent] = []

    for record in lattice.cells:
        r, c = record.index
        geometry = record.geometry
        prefix = f"cell_{r}_{c}"
        components.append(
            _component(
                f"{prefix}_outer_plate",
                "outer_plate",
                record.index,
                _extruded_polygon(geometry.outer_part, plate_thickness),
            )
        )
        components.append(
            _component(
                f"{prefix}_inner_plate",
                "inner_plate",
                record.index,
                _extruded_polygon(geometry.inner_part, plate_thickness),
            )
        )
        if include_pins:
            pin_radius = max(
                0.02 * config.cell_size,
                config.pin_radius * config.cell_size * 0.32,
            )
            for joint in (*geometry.outer_joints, *geometry.inner_joints):
                components.append(
                    _component(
                        f"{prefix}_{joint.part}_pin_{joint.index}",
                        "pin",
                        record.index,
                        _cylinder(joint.position, pin_radius, pin_height, pin_segments),
                    )
                )

    if include_connectors:
        for index, connector in enumerate(lattice.connectors):
            name = (
                f"connector_{index}_{connector.first[0]}_{connector.first[1]}"
                f"_to_{connector.second[0]}_{connector.second[1]}"
            )
            components.append(
                _component(
                    name,
                    "connector",
                    None,
                    _bar_between(connector.start, connector.end, connector_width),
                )
            )

    return PaperRADMesh(components=tuple(components), source=lattice)


def export_paper_rad_mesh_obj(mesh: PaperRADMesh) -> str:
    lines = [
        "# RAD paper lattice normalized OBJ export",
        f"# components {mesh.component_count}",
        f"# vertices {mesh.vertex_count}",
        f"# faces {mesh.face_count}",
    ]
    vertex_offset = 1
    for component in mesh.components:
        lines.append(f"o {component.name}")
        lines.append(f"# kind {component.kind}")
        if component.cell is not None:
            lines.append(f"# cell {component.cell[0]} {component.cell[1]}")
        for vertex in component.vertices:
            lines.append(f"v {vertex[0]:.9g} {vertex[1]:.9g} {vertex[2]:.9g}")
        for face in component.faces:
            a, b, c = face + vertex_offset
            lines.append(f"f {a} {b} {c}")
        vertex_offset += component.vertex_count
    return "\n".join(lines) + "\n"


def iter_obj_vertices(obj_text: str) -> Iterable[tuple[float, float, float]]:
    for line in obj_text.splitlines():
        if not line.startswith("v "):
            continue
        _, x, y, z = line.split()
        yield float(x), float(y), float(z)
