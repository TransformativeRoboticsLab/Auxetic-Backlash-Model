from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from .models import SimulationResult


LOCK_FILENAME_RE = re.compile(
    r"^AOM_L(?P<count>\d+)"
    r"(?:_(?P<cells>\d+(?:_\d+)*?)_(?P<angles>\d+))?"
    r"_(?P<state>\d+)\.csv$"
)


@dataclass(frozen=True)
class LockConfigurationRecord:
    filename: str
    lock_count: int
    lock_cells: tuple[int, ...]
    lock_angles: tuple[int, ...]
    state_index: int
    stable_state_count: int
    realized_lock_cells: tuple[int, ...]
    failed_lock_cells: tuple[int, ...]
    marker_count: int
    centered_markers: np.ndarray
    raw_markers: np.ndarray
    claim_label: str = "measured single-strand lock configuration"


@dataclass(frozen=True)
class LockDataset:
    records: tuple[LockConfigurationRecord, ...]
    source_path: Path
    marker_summary_file: Path
    strand_cells: int = 12

    @property
    def configuration_count(self) -> int:
        return len(self.records)

    @property
    def marker_count(self) -> int:
        return max((record.marker_count for record in self.records), default=0)


def parse_lock_filename(filename: str) -> tuple[int, tuple[int, ...], tuple[int, ...], int]:
    match = LOCK_FILENAME_RE.match(filename)
    if not match:
        raise ValueError(f"unsupported lock dataset filename: {filename}")
    lock_count = int(match.group("count"))
    state_index = int(match.group("state"))
    cell_text = match.group("cells")
    angle_text = match.group("angles")
    lock_cells = tuple(int(value) for value in cell_text.split("_")) if cell_text else ()
    if angle_text:
        if len(angle_text) % 2 != 0:
            raise ValueError(f"lock angle block must be two-digit chunks: {filename}")
        lock_angles = tuple(
            int(angle_text[index : index + 2])
            for index in range(0, len(angle_text), 2)
        )
    else:
        lock_angles = ()
    if lock_count != len(lock_cells):
        raise ValueError(f"lock count does not match lock cells in {filename}")
    if lock_angles and len(lock_angles) != lock_count:
        raise ValueError(f"lock count does not match lock angles in {filename}")
    return lock_count, lock_cells, lock_angles, state_index


def _consecutive_failed_lock_cells(lock_cells: tuple[int, ...]) -> tuple[int, ...]:
    known_failures = {
        (1, 2, 3): (2,),
        (4, 5, 6): (5,),
        (7, 8, 9): (8,),
        (10, 11, 12): (10, 11),
    }
    ordered_tuple = tuple(sorted(lock_cells))
    if ordered_tuple in known_failures:
        return known_failures[ordered_tuple]
    failed: list[int] = []
    ordered = list(ordered_tuple)
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end] == ordered[end - 1] + 1:
            end += 1
        run = ordered[start:end]
        if len(run) >= 3:
            failed.extend(run[1:-1])
        start = end
    return tuple(failed)


def realized_lock_cells(
    lock_cells: Iterable[int],
    failed_lock_cells: Iterable[int] | None = None,
) -> tuple[int, ...]:
    failed = set(failed_lock_cells or ())
    return tuple(cell for cell in lock_cells if cell not in failed)


def _load_marker_summary(path: Path) -> dict[str, list[tuple[str, float, float, float]]]:
    rows: dict[str, list[tuple[str, float, float, float]]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"file", "marker", "mean_x", "mean_y", "mean_z"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError(f"marker summary missing columns: {sorted(required)}")
        for raw in reader:
            filename = str(raw["file"])
            rows.setdefault(filename, []).append(
                (
                    str(raw["marker"]),
                    float(raw["mean_x"]),
                    float(raw["mean_y"]),
                    float(raw["mean_z"]),
                )
            )
    for values in rows.values():
        values.sort(key=lambda item: item[1])
    return rows


def _stable_counts(summary: dict[str, list[tuple[str, float, float, float]]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for filename in summary:
        try:
            lock_count, lock_cells, lock_angles, _state_index = parse_lock_filename(filename)
        except ValueError:
            continue
        key = _configuration_key(lock_count, lock_cells, lock_angles)
        counts[key] = counts.get(key, 0) + 1
    return counts


def _configuration_key(
    lock_count: int,
    lock_cells: tuple[int, ...],
    lock_angles: tuple[int, ...],
) -> str:
    cells = "_".join(str(cell) for cell in lock_cells)
    angles = "".join(f"{angle:02d}" for angle in lock_angles)
    return f"L{lock_count}:{cells}:{angles}"


def _center_markers(raw: np.ndarray) -> np.ndarray:
    if raw.size == 0:
        return raw.reshape(0, 3)
    centered = raw - np.mean(raw, axis=0, keepdims=True)
    # Use the long strand direction as a canonical positive axis.
    order = np.argsort(centered[:, 2])
    centered = centered[order]
    if centered[-1, 2] < centered[0, 2]:
        centered = centered[::-1]
    return centered


def load_lock_dataset(
    logs_path: str | Path = r"C:\Users\ahyan\OneDrive\Desktop\logs",
    *,
    strand_cells: int = 12,
) -> LockDataset:
    source = Path(logs_path)
    summary_path = source / "marker_xyz_means_all_37_files.csv"
    if not summary_path.exists():
        raise FileNotFoundError(f"missing marker summary file: {summary_path}")
    summary = _load_marker_summary(summary_path)
    counts = _stable_counts(summary)
    records: list[LockConfigurationRecord] = []
    for filename, markers in sorted(summary.items()):
        try:
            lock_count, lock_cells, lock_angles, state_index = parse_lock_filename(filename)
        except ValueError:
            continue
        failed = _consecutive_failed_lock_cells(lock_cells)
        raw = np.array([[x, y, z] for _marker, x, y, z in markers], dtype=float)
        centered = _center_markers(raw)
        records.append(
            LockConfigurationRecord(
                filename=filename,
                lock_count=lock_count,
                lock_cells=lock_cells,
                lock_angles=lock_angles,
                state_index=state_index,
                stable_state_count=counts[_configuration_key(lock_count, lock_cells, lock_angles)],
                realized_lock_cells=realized_lock_cells(lock_cells, failed),
                failed_lock_cells=failed,
                marker_count=int(centered.shape[0]),
                centered_markers=centered,
                raw_markers=raw,
            )
        )
    return LockDataset(tuple(records), source, summary_path, strand_cells=strand_cells)


def lock_mask(lock_cells: Iterable[int], strand_cells: int = 12) -> np.ndarray:
    mask = np.zeros(strand_cells, dtype=float)
    for cell in lock_cells:
        index = int(cell) - 1
        if 0 <= index < strand_cells:
            mask[index] = 1.0
    return mask


def lock_angle_vector(
    lock_cells: Iterable[int],
    lock_angles: Iterable[int],
    strand_cells: int = 12,
) -> np.ndarray:
    values = np.zeros(strand_cells, dtype=float)
    for cell, angle in zip(lock_cells, lock_angles, strict=False):
        index = int(cell) - 1
        if 0 <= index < strand_cells:
            values[index] = float(angle)
    return values


def lock_feature_vector(
    lock_cells: Iterable[int],
    lock_angles: Iterable[int],
    *,
    strand_cells: int = 12,
    state_index: int = 1,
    realized_cells: Iterable[int] | None = None,
) -> np.ndarray:
    nominal = lock_mask(lock_cells, strand_cells)
    realized = lock_mask(realized_cells if realized_cells is not None else lock_cells, strand_cells)
    angles = lock_angle_vector(lock_cells, lock_angles, strand_cells) / 40.0
    return np.concatenate(
        [
            nominal,
            realized,
            nominal - realized,
            angles,
            np.array([float(state_index)], dtype=float),
        ]
    )


def _resample_shape(markers: np.ndarray, count: int) -> np.ndarray:
    if markers.shape[0] == count:
        return markers.copy()
    if markers.shape[0] == 0:
        return np.zeros((count, 3), dtype=float)
    t_old = np.linspace(0.0, 1.0, markers.shape[0])
    t_new = np.linspace(0.0, 1.0, count)
    out = np.zeros((count, 3), dtype=float)
    for axis in range(3):
        out[:, axis] = np.interp(t_new, t_old, markers[:, axis])
    return out


def _markers_to_sim_axes(markers: np.ndarray, strand_axis: str) -> np.ndarray:
    """Map mocap axes into simulator axes.

    The lab plotting script treats X/Z as the ground plane and Y as vertical.
    The simulator uses X/Y as the ground plane and Z as vertical.
    """

    if strand_axis == "row":
        return np.column_stack([markers[:, 2], markers[:, 0], markers[:, 1]])
    if strand_axis == "column":
        return np.column_stack([markers[:, 0], markers[:, 2], markers[:, 1]])
    raise ValueError("strand_axis must be 'row' or 'column'")


def _rotation_between_vectors(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    source_norm = np.linalg.norm(source)
    target_norm = np.linalg.norm(target)
    if source_norm < 1e-12 or target_norm < 1e-12:
        return np.eye(3)
    source_unit = source / source_norm
    target_unit = target / target_norm
    dot = float(np.clip(np.dot(source_unit, target_unit), -1.0, 1.0))
    if np.isclose(dot, 1.0):
        return np.eye(3)
    if np.isclose(dot, -1.0):
        basis = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(source_unit, basis))) > 0.9:
            basis = np.array([0.0, 1.0, 0.0])
        axis = np.cross(source_unit, basis)
        axis = axis / np.linalg.norm(axis)
        return -np.eye(3) + 2.0 * np.outer(axis, axis)
    axis = np.cross(source_unit, target_unit)
    skew = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ]
    )
    axis_norm_sq = float(np.dot(axis, axis))
    return np.eye(3) + skew + skew @ skew * ((1.0 - dot) / axis_norm_sq)


def _position_anchor_indices(
    position_locked_mask: np.ndarray | None,
    count: int,
) -> tuple[int, int] | tuple[()]:
    if position_locked_mask is None:
        return ()
    mask = np.asarray(position_locked_mask, dtype=bool).reshape(-1)
    if mask.shape[0] != count:
        raise ValueError(f"position lock mask must contain {count} values")
    anchors = np.flatnonzero(mask)
    if anchors.size < 2:
        return ()
    return int(anchors[0]), int(anchors[-1])


def _align_coordinates_to_position_locks(
    measured_coordinates: np.ndarray,
    simulated_coordinates: np.ndarray,
    position_locked_mask: np.ndarray | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    measured = np.asarray(measured_coordinates, dtype=float)
    simulated = np.asarray(simulated_coordinates, dtype=float)
    if measured.shape != simulated.shape or measured.ndim != 2 or measured.shape[1] != 3:
        raise ValueError("measured and simulated coordinates must both have shape (N, 3)")
    anchors = _position_anchor_indices(position_locked_mask, measured.shape[0])
    if anchors:
        first, last = anchors
        source_vec = measured[last] - measured[first]
        target_vec = simulated[last] - simulated[first]
        source_length = float(np.linalg.norm(source_vec))
        target_length = float(np.linalg.norm(target_vec))
        if source_length > 1e-12 and target_length > 1e-12:
            rotation = _rotation_between_vectors(source_vec, target_vec)
            scale = target_length / source_length
            aligned = (measured - measured[first]) @ rotation.T * scale + simulated[first]
            aligned[first] = simulated[first]
            aligned[last] = simulated[last]
            return aligned, {
                "endpointAnchored": True,
                "anchors": anchors,
                "scaleAppliedToMeasuredData": float(scale),
            }
    measured_centered = measured - np.mean(measured, axis=0, keepdims=True)
    simulated_center = np.mean(simulated, axis=0, keepdims=True)
    measured_extent = float(np.ptp(measured[:, 0]))
    simulated_extent = float(np.ptp(simulated[:, 0]))
    if measured_extent > 1e-12 and simulated_extent > 1e-12:
        scale = simulated_extent / measured_extent
    else:
        scale = np.linalg.norm(simulated - simulated_center) / max(
            np.linalg.norm(measured_centered),
            1e-12,
        )
    aligned = measured_centered * scale + simulated_center
    return aligned, {
        "endpointAnchored": False,
        "anchors": (),
        "scaleAppliedToMeasuredData": float(scale),
    }


def anchor_coordinates_to_position_locks(
    measured_coordinates: np.ndarray,
    simulated_coordinates: np.ndarray,
    position_locked_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Align measured strand coordinates to simulated positional lock anchors."""

    aligned, _metadata = _align_coordinates_to_position_locks(
        measured_coordinates,
        simulated_coordinates,
        position_locked_mask,
    )
    return aligned


def _dataset_arrays(
    dataset: LockDataset,
) -> tuple[np.ndarray, np.ndarray, int]:
    marker_count = dataset.marker_count
    features = np.vstack(
        [
            lock_feature_vector(
                record.lock_cells,
                record.lock_angles,
                strand_cells=dataset.strand_cells,
                state_index=record.state_index,
                realized_cells=record.realized_lock_cells,
            )
            for record in dataset.records
        ]
    )
    targets = np.vstack(
        [_resample_shape(record.centered_markers, marker_count).reshape(1, -1) for record in dataset.records]
    )
    return features, targets, marker_count


def lock_record_cell_coordinates(
    record: LockConfigurationRecord,
    *,
    strand_cells: int = 12,
    strand_axis: str = "row",
) -> np.ndarray:
    """Return measured marker geometry resampled as per-cell simulator coordinates."""

    markers = _resample_shape(record.centered_markers, strand_cells)
    markers = _markers_to_sim_axes(markers, strand_axis)
    return markers - np.mean(markers, axis=0, keepdims=True)


def lock_coordinate_training_arrays(
    dataset: LockDataset,
    *,
    strand_axis: str = "row",
) -> tuple[np.ndarray, np.ndarray]:
    """Build inputs and per-cell 3D-coordinate targets for learned surrogates."""

    features = []
    targets = []
    for record in dataset.records:
        features.append(
            lock_feature_vector(
                record.lock_cells,
                record.lock_angles,
                strand_cells=dataset.strand_cells,
                state_index=record.state_index,
                realized_cells=record.realized_lock_cells,
            )
        )
        targets.append(
            lock_record_cell_coordinates(
                record,
                strand_cells=dataset.strand_cells,
                strand_axis=strand_axis,
            ).reshape(-1)
        )
    return np.vstack(features), np.vstack(targets)


def lock_coordinate_dataset(
    dataset: LockDataset,
    *,
    strand_axis: str = "row",
) -> dict[str, object]:
    records = []
    for record in dataset.records:
        coordinates = lock_record_cell_coordinates(
            record,
            strand_cells=dataset.strand_cells,
            strand_axis=strand_axis,
        )
        records.append(
            {
                "filename": record.filename,
                "lockCount": record.lock_count,
                "lockCells": record.lock_cells,
                "lockAngles": record.lock_angles,
                "stateIndex": record.state_index,
                "stableStateCount": record.stable_state_count,
                "failedLockCells": record.failed_lock_cells,
                "realizedLockCells": record.realized_lock_cells,
                "featureVector": lock_feature_vector(
                    record.lock_cells,
                    record.lock_angles,
                    strand_cells=dataset.strand_cells,
                    state_index=record.state_index,
                    realized_cells=record.realized_lock_cells,
                ),
                "cellCoordinates": coordinates,
                "rawMarkerCount": record.marker_count,
            }
        )
    return {
        "schema": "rad-sim.lock-coordinate-dataset.v1",
        "sourcePath": str(dataset.source_path),
        "markerSummaryFile": str(dataset.marker_summary_file),
        "strandCells": dataset.strand_cells,
        "strandAxis": strand_axis,
        "configurationCount": dataset.configuration_count,
        "inputDimension": 4 * dataset.strand_cells + 1,
        "outputDimension": 3 * dataset.strand_cells,
        "records": records,
        "claimLabels": {
            "target": "per-cell 3D coordinate targets for empirical surrogates",
            "axes": "lab X/Z ground plane and Y vertical remapped to simulator X/Y/Z",
        },
    }


def export_lock_coordinate_dataset_json(
    dataset: LockDataset,
    *,
    strand_axis: str = "row",
) -> str:
    return json.dumps(
        lock_coordinate_dataset(dataset, strand_axis=strand_axis),
        indent=2,
        default=_json_default,
    )


def predict_lock_strand_shape(
    dataset: LockDataset,
    lock_cells: Iterable[int],
    lock_angles: Iterable[int] | None = None,
    *,
    state_index: int = 1,
    realized_cells: Iterable[int] | None = None,
    method: str = "rbf",
    length_scale: float = 1.5,
) -> dict[str, object]:
    """Predict a measured 12-cell strand shape from the lock dataset.

    The default is an RBF/nearest-neighbor surrogate. With only 37 measured
    shapes this is deliberately conservative; a large neural network would
    mostly memorize the data.
    """

    cells = tuple(int(cell) for cell in lock_cells)
    angles = tuple(int(angle) for angle in (lock_angles or (30,) * len(cells)))
    failed = _consecutive_failed_lock_cells(cells)
    realized = tuple(realized_cells) if realized_cells is not None else realized_lock_cells(cells, failed)
    features, targets, marker_count = _dataset_arrays(dataset)
    query = lock_feature_vector(
        cells,
        angles,
        strand_cells=dataset.strand_cells,
        state_index=state_index,
        realized_cells=realized,
    )
    distances = np.linalg.norm(features - query[None, :], axis=1)
    nearest_order = np.argsort(distances)
    if method == "nearest" or np.isclose(distances[nearest_order[0]], 0.0):
        weights = np.zeros_like(distances)
        weights[nearest_order[0]] = 1.0
    elif method == "rbf":
        scale = max(float(length_scale), 1e-9)
        weights = np.exp(-(distances**2) / (2.0 * scale * scale))
        if float(np.sum(weights)) <= 0:
            weights = np.zeros_like(distances)
            weights[nearest_order[0]] = 1.0
        else:
            weights = weights / np.sum(weights)
    else:
        raise ValueError("method must be 'rbf' or 'nearest'")
    predicted = (weights[:, None] * targets).sum(axis=0).reshape(marker_count, 3)
    neighbors = [
        {
            "filename": dataset.records[index].filename,
            "distance": float(distances[index]),
            "weight": float(weights[index]),
            "lockCells": dataset.records[index].lock_cells,
            "realizedLockCells": dataset.records[index].realized_lock_cells,
            "stateIndex": dataset.records[index].state_index,
        }
        for index in nearest_order[: min(5, len(nearest_order))]
    ]
    return {
        "schema": "rad-sim.single-strand-lock-prediction.v1",
        "method": method,
        "lockCells": cells,
        "lockAngles": angles,
        "stateIndex": int(state_index),
        "failedLockCells": failed,
        "realizedLockCells": tuple(int(cell) for cell in realized),
        "markerCount": int(marker_count),
        "predictedMarkers": predicted,
        "nearestNeighbors": neighbors,
        "claimLabels": {
            "surrogate": "empirical measured-shape surrogate",
            "physics": "not a first-principles lock mechanics law",
            "failureMode": "consecutive middle locks may be nominal but unrealized",
        },
    }


def lock_dataset_summary(dataset: LockDataset) -> dict[str, object]:
    stable_counts = {}
    failed_records = []
    for record in dataset.records:
        stable_counts.setdefault(record.stable_state_count, 0)
        stable_counts[record.stable_state_count] += 1
        if record.failed_lock_cells:
            failed_records.append(
                {
                    "filename": record.filename,
                    "lockCells": record.lock_cells,
                    "failedLockCells": record.failed_lock_cells,
                    "realizedLockCells": record.realized_lock_cells,
                }
            )
    return {
        "schema": "rad-sim.single-strand-lock-dataset-summary.v1",
        "sourcePath": str(dataset.source_path),
        "markerSummaryFile": str(dataset.marker_summary_file),
        "configurationCount": dataset.configuration_count,
        "markerCount": dataset.marker_count,
        "strandCells": dataset.strand_cells,
        "stableStateHistogram": stable_counts,
        "failedConsecutiveLockRecords": failed_records,
        "claimLabels": {
            "dataset": "measured 12-cell single-strand lock dataset",
            "lockModel": "calibration evidence for nominal-to-realized lock operators",
        },
    }


def export_lock_dataset_summary_json(dataset: LockDataset) -> str:
    return json.dumps(lock_dataset_summary(dataset), indent=2, default=_json_default)


def _json_default(value: object) -> object:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return list(value)
    raise TypeError(f"object is not JSON serializable: {type(value)!r}")


def calibrate_single_strand_from_lock_dataset(
    result: SimulationResult,
    dataset: LockDataset,
    *,
    state_index: int = 1,
    method: str = "rbf",
) -> dict[str, object]:
    """Attach measured single-strand marker geometry to a simulation result.

    This does not mutate the lattice solver. It provides an empirical measured
    strand overlay that can be compared against the current kinematic centers.
    """

    if result.config.rows == 1:
        strand_axis = "row"
        lock_cells = tuple(int(col + 1) for col in np.nonzero(result.state.locked_mask[0])[0])
        position_locked = result.state.position_locked_mask[0]
    elif result.config.cols == 1:
        strand_axis = "column"
        lock_cells = tuple(int(row + 1) for row in np.nonzero(result.state.locked_mask[:, 0])[0])
        position_locked = result.state.position_locked_mask[:, 0]
    else:
        raise ValueError("single-strand calibration requires a 1 x N or N x 1 lattice")
    prediction = predict_lock_strand_shape(
        dataset,
        lock_cells,
        (30,) * len(lock_cells),
        state_index=state_index,
        method=method,
    )
    simulated = (
        result.deformed_centers_3d[0]
        if result.config.rows == 1
        else result.deformed_centers_3d[:, 0, :]
    )
    measured = _resample_shape(
        np.asarray(prediction["predictedMarkers"], dtype=float),
        simulated.shape[0],
    )
    measured = _markers_to_sim_axes(measured, strand_axis)
    aligned_measured, alignment = _align_coordinates_to_position_locks(
        measured,
        simulated,
        position_locked,
    )
    residual = aligned_measured - simulated
    return {
        "schema": "rad-sim.single-strand-lock-calibration-overlay.v1",
        "strandAxis": strand_axis,
        "lockCells": lock_cells,
        "positionLockedCells": tuple(int(index + 1) for index in np.flatnonzero(position_locked)),
        "prediction": {
            key: value
            for key, value in prediction.items()
            if key != "predictedMarkers"
        },
        "measuredMarkers": aligned_measured,
        "simulatedCenters": simulated.copy(),
        "residual": residual,
        "metrics": {
            "rmsResidual": float(np.sqrt(np.mean(residual**2))),
            "maxResidual": float(np.max(np.linalg.norm(residual, axis=1))),
            "scaleAppliedToMeasuredData": float(alignment["scaleAppliedToMeasuredData"]),
            "endpointAnchored": bool(alignment["endpointAnchored"]),
            "positionLockedAnchors": tuple(int(index + 1) for index in alignment["anchors"]),
        },
        "claimLabels": {
            "overlay": "empirical lock-dataset alignment for single-strand simulation",
            "solver": "diagnostic overlay only; does not prove kinematic solver accuracy",
        },
    }


def calibrate_repeated_row_locks_from_lock_dataset(
    result: SimulationResult,
    dataset: LockDataset,
    *,
    state_index: int = 1,
    method: str = "rbf",
) -> dict[str, object]:
    """Build a measured-sheet overlay for identical lock patterns on each row.

    This is for the experimental case where a 12-cell measured strand is used
    as the row template and the same nominal locks are repeated across rows.
    """

    if result.config.cols < 1:
        raise ValueError("lattice must have at least one column")
    row_patterns = [
        tuple(int(col + 1) for col in np.nonzero(result.state.locked_mask[row])[0])
        for row in range(result.config.rows)
    ]
    first_pattern = row_patterns[0] if row_patterns else ()
    if any(pattern != first_pattern for pattern in row_patterns):
        raise ValueError("repeated-row calibration requires identical lock masks on every row")
    prediction = predict_lock_strand_shape(
        dataset,
        first_pattern,
        (30,) * len(first_pattern),
        state_index=state_index,
        method=method,
    )
    measured_row = _markers_to_sim_axes(
        _resample_shape(np.asarray(prediction["predictedMarkers"], dtype=float), result.config.cols),
        "row",
    )
    sheet = np.zeros_like(result.deformed_centers_3d, dtype=float)
    endpoint_anchored_rows = 0
    position_anchor_cells = 0
    scales: list[float] = []
    for row in range(result.config.rows):
        row_position_locked = result.state.position_locked_mask[row]
        aligned_row, alignment = _align_coordinates_to_position_locks(
            measured_row,
            result.deformed_centers_3d[row],
            row_position_locked,
        )
        sheet[row] = aligned_row
        if alignment["endpointAnchored"]:
            endpoint_anchored_rows += 1
        position_anchor_cells += len(alignment["anchors"])
        scales.append(float(alignment["scaleAppliedToMeasuredData"]))
    residual = sheet - result.deformed_centers_3d
    return {
        "schema": "rad-sim.repeated-row-lock-calibration-overlay.v1",
        "rowLockPattern": first_pattern,
        "rowCount": int(result.config.rows),
        "colCount": int(result.config.cols),
        "prediction": {
            key: value
            for key, value in prediction.items()
            if key != "predictedMarkers"
        },
        "calibratedSheetCenters": sheet,
        "kinematicCenters": result.deformed_centers_3d.copy(),
        "residual": residual,
        "metrics": {
            "rmsResidual": float(np.sqrt(np.mean(residual**2))),
            "maxResidual": float(np.max(np.linalg.norm(residual.reshape(-1, 3), axis=1))),
            "scaleAppliedToMeasuredData": float(np.mean(scales)) if scales else 1.0,
            "endpointAnchoredRows": int(endpoint_anchored_rows),
            "positionLockedAnchorCells": int(position_anchor_cells),
        },
        "claimLabels": {
            "overlay": "measured single-strand shape repeated row-wise across a sheet",
            "boundaryCondition": "valid when each row has the same nominal lock mask",
            "solver": "calibration overlay only; does not replace the full sheet equilibrium solver",
        },
    }
