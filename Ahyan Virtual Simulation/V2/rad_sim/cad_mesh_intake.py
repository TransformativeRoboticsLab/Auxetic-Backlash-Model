from __future__ import annotations

import argparse
import json
import math
import re
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from .cell_geometry import RADHardwareProfile, calibration_readiness


CAD_MESH_AUDIT_SCHEMA = "rad-sim.cad-mesh-audit.v1"
SUPPORTED_CAD_BOUND_FORMATS = (".obj", ".stl", ".step", ".stp")


@dataclass(frozen=True)
class MeshBounds:
    minimum: np.ndarray
    maximum: np.ndarray
    vertex_count: int
    face_count: int

    @property
    def size(self) -> np.ndarray:
        return self.maximum - self.minimum

    @property
    def planar_span_mm(self) -> float:
        size = self.size
        return float(max(size[0], size[1]))

    @property
    def thickness_mm(self) -> float:
        return float(self.size[2])


def _iter_obj_vertex_faces(text: str) -> tuple[list[tuple[float, float, float]], int]:
    vertices: list[tuple[float, float, float]] = []
    face_count = 0
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("v "):
            parts = stripped.split()
            if len(parts) >= 4:
                vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
        elif stripped.startswith("f "):
            face_count += 1
    return vertices, face_count


def _is_binary_stl(data: bytes) -> bool:
    if len(data) < 84:
        return False
    face_count = struct.unpack_from("<I", data, 80)[0]
    return 84 + face_count * 50 == len(data)


def _iter_binary_stl(data: bytes) -> tuple[list[tuple[float, float, float]], int]:
    face_count = struct.unpack_from("<I", data, 80)[0]
    vertices: list[tuple[float, float, float]] = []
    offset = 84
    for _ in range(face_count):
        offset += 12
        for _vertex_index in range(3):
            vertices.append(struct.unpack_from("<fff", data, offset))
            offset += 12
        offset += 2
    return vertices, int(face_count)


def _iter_ascii_stl(text: str) -> tuple[list[tuple[float, float, float]], int]:
    vertices: list[tuple[float, float, float]] = []
    face_count = 0
    for line in text.splitlines():
        parts = line.strip().split()
        if len(parts) == 4 and parts[0].lower() == "vertex":
            vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
        elif len(parts) >= 2 and parts[0].lower() == "facet" and parts[1].lower() == "normal":
            face_count += 1
    return vertices, face_count


_STEP_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
_STEP_POINT_RE = re.compile(
    rf"CARTESIAN_POINT\s*\(\s*(?:'[^']*'|\$)?\s*,\s*\(\s*({_STEP_NUMBER})\s*,\s*({_STEP_NUMBER})\s*,\s*({_STEP_NUMBER})\s*\)\s*\)",
    re.IGNORECASE,
)


def _iter_step_cartesian_points(text: str) -> tuple[list[tuple[float, float, float]], int]:
    vertices = [
        (float(match.group(1)), float(match.group(2)), float(match.group(3)))
        for match in _STEP_POINT_RE.finditer(text)
    ]
    face_count = len(re.findall(r"\bADVANCED_FACE\s*\(", text, flags=re.IGNORECASE))
    return vertices, face_count


def mesh_bounds_from_file(path: str | Path) -> MeshBounds:
    source = Path(path)
    suffix = source.suffix.lower()
    data = source.read_bytes()
    if suffix == ".obj":
        vertices, face_count = _iter_obj_vertex_faces(data.decode("utf-8", errors="ignore"))
    elif suffix == ".stl":
        if _is_binary_stl(data):
            vertices, face_count = _iter_binary_stl(data)
        else:
            vertices, face_count = _iter_ascii_stl(data.decode("utf-8", errors="ignore"))
    elif suffix in {".step", ".stp"}:
        vertices, face_count = _iter_step_cartesian_points(data.decode("utf-8", errors="ignore"))
    else:
        raise ValueError("CAD mesh audit supports .obj, .stl, .step, and .stp exports")
    if not vertices:
        raise ValueError(f"no mesh vertices found in {source}")
    array = np.asarray(vertices, dtype=float)
    return MeshBounds(
        minimum=array.min(axis=0),
        maximum=array.max(axis=0),
        vertex_count=int(array.shape[0]),
        face_count=int(face_count),
    )


def hardware_profile_from_mesh_bounds(
    bounds: MeshBounds,
    *,
    name: str = "cad-mesh-derived",
    source: str = "exported CAD mesh bounds",
    side_length_mm: float | None = None,
    pin_radius_mm: float | None = None,
    hole_radius_mm: float | None = None,
    backlash_mm: float | None = None,
    boss_radius_mm: float | None = None,
    joint_stack_height_mm: float | None = None,
    fabrication_hole_tolerance_mm: float = 0.1,
    notes: str = "",
) -> RADHardwareProfile:
    inferred_side = float(side_length_mm) if side_length_mm is not None else bounds.planar_span_mm
    inferred_stack = float(joint_stack_height_mm) if joint_stack_height_mm is not None else bounds.thickness_mm
    return RADHardwareProfile(
        name=name,
        source=source,
        side_length_mm=inferred_side,
        fabrication_hole_tolerance_mm=fabrication_hole_tolerance_mm,
        backlash_mm=backlash_mm,
        pin_radius_mm=pin_radius_mm,
        hole_radius_mm=hole_radius_mm,
        plate_thickness_mm=bounds.thickness_mm,
        joint_stack_height_mm=inferred_stack,
        boss_radius_mm=boss_radius_mm,
        notes=notes,
    )


def cad_mesh_audit(
    path: str | Path,
    *,
    name: str = "cad-mesh-derived",
    source: str | None = None,
    side_length_mm: float | None = None,
    pin_radius_mm: float | None = None,
    hole_radius_mm: float | None = None,
    backlash_mm: float | None = None,
    boss_radius_mm: float | None = None,
    joint_stack_height_mm: float | None = None,
    fabrication_hole_tolerance_mm: float = 0.1,
) -> dict[str, object]:
    source_path = Path(path)
    bounds = mesh_bounds_from_file(source_path)
    profile = hardware_profile_from_mesh_bounds(
        bounds,
        name=name,
        source=source or str(source_path),
        side_length_mm=side_length_mm,
        pin_radius_mm=pin_radius_mm,
        hole_radius_mm=hole_radius_mm,
        backlash_mm=backlash_mm,
        boss_radius_mm=boss_radius_mm,
        joint_stack_height_mm=joint_stack_height_mm,
        fabrication_hole_tolerance_mm=fabrication_hole_tolerance_mm,
        notes="Mesh bounds derived from exported CAD; pin/hole/boss values require segmentation or explicit measurements.",
    )
    readiness = calibration_readiness(profile)
    size = bounds.size
    finite_bounds = all(math.isfinite(float(value)) for value in (*bounds.minimum, *bounds.maximum))
    return {
        "schema": CAD_MESH_AUDIT_SCHEMA,
        "sourceFile": str(source_path).replace("\\", "/"),
        "supportedFormats": list(SUPPORTED_CAD_BOUND_FORMATS),
        "mesh": {
            "format": source_path.suffix.lower(),
            "vertexCount": bounds.vertex_count,
            "faceCount": bounds.face_count,
            "boundsMm": {
                "min": bounds.minimum.tolist(),
                "max": bounds.maximum.tolist(),
                "size": size.tolist(),
            },
            "planarSpanMm": bounds.planar_span_mm,
            "thicknessMm": bounds.thickness_mm,
            "finiteBounds": finite_bounds,
        },
        "hardwareProfile": profile.to_dict(),
        "calibrationReadiness": {
            "level": readiness.level,
            "visualReady": readiness.visual_ready,
            "meshReady": readiness.mesh_ready,
            "solverReady": readiness.solver_ready,
            "missingFields": list(readiness.missing_fields),
            "meshMissingFields": list(readiness.mesh_missing_fields),
            "solverGaps": list(readiness.solver_gaps),
            "summary": readiness.summary,
        },
        "claimBoundary": {
            "physicalAccuracyValidated": False,
            "reason": "Mesh bounds calibrate visible geometry scale only; exact contact requires segmented pin-hole surfaces and bench holdouts.",
        },
    }


def export_cad_mesh_audit_json(path: str | Path, **kwargs: object) -> str:
    return json.dumps(cad_mesh_audit(path, **kwargs), indent=2)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit an exported RAD cell OBJ/STL mesh and derive a hardware profile.")
    parser.add_argument("mesh", help="Path to exported one-cell OBJ or STL")
    parser.add_argument("--out", default="", help="Optional JSON output path")
    parser.add_argument("--name", default="cad-mesh-derived")
    parser.add_argument("--side-length-mm", type=float, default=None)
    parser.add_argument("--pin-radius-mm", type=float, default=None)
    parser.add_argument("--hole-radius-mm", type=float, default=None)
    parser.add_argument("--backlash-mm", type=float, default=None)
    parser.add_argument("--boss-radius-mm", type=float, default=None)
    parser.add_argument("--joint-stack-height-mm", type=float, default=None)
    parser.add_argument("--fabrication-hole-tolerance-mm", type=float, default=0.1)
    args = parser.parse_args(list(argv) if argv is not None else None)
    payload = export_cad_mesh_audit_json(
        args.mesh,
        name=args.name,
        side_length_mm=args.side_length_mm,
        pin_radius_mm=args.pin_radius_mm,
        hole_radius_mm=args.hole_radius_mm,
        backlash_mm=args.backlash_mm,
        boss_radius_mm=args.boss_radius_mm,
        joint_stack_height_mm=args.joint_stack_height_mm,
        fabrication_hole_tolerance_mm=args.fabrication_hole_tolerance_mm,
    )
    if args.out:
        Path(args.out).write_text(payload, encoding="utf-8")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
