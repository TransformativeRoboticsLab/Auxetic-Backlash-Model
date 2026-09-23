from __future__ import annotations

import json
from typing import Iterable

import numpy as np

from .lock_dataset import (
    LockDataset,
    _consecutive_failed_lock_cells,
    lock_coordinate_training_arrays,
    lock_feature_vector,
    realized_lock_cells,
)


def _json_default(value: object) -> object:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, tuple):
        return list(value)
    raise TypeError(f"object is not JSON serializable: {type(value)!r}")


def _standardize(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = np.mean(values, axis=0, keepdims=True)
    scale = np.std(values, axis=0, keepdims=True)
    scale = np.where(scale < 1e-9, 1.0, scale)
    return (values - mean) / scale, mean.reshape(-1), scale.reshape(-1)


def _forward(
    x: np.ndarray,
    w1: np.ndarray,
    b1: np.ndarray,
    w2: np.ndarray,
    b2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    hidden = np.tanh(x @ w1 + b1)
    return hidden @ w2 + b2, hidden


def train_lock_coordinate_mlp(
    dataset: LockDataset,
    *,
    hidden_units: int = 24,
    epochs: int = 800,
    learning_rate: float = 0.01,
    weight_decay: float = 1e-4,
    validation_fraction: float = 0.2,
    seed: int = 7,
    strand_axis: str = "row",
) -> dict[str, object]:
    """Train a small empirical MLP from lock state to per-cell 3D coordinates.

    The purpose is to establish the neural-surrogate pipeline. With the current
    small dataset, the model should be treated as an empirical interpolator and
    not as a validated general mechanics law.
    """

    x_raw, y_raw = lock_coordinate_training_arrays(dataset, strand_axis=strand_axis)
    if x_raw.shape[0] < 2:
        raise ValueError("at least two measured configurations are required for MLP training")
    x, x_mean, x_scale = _standardize(x_raw)
    y, y_mean, y_scale = _standardize(y_raw)
    rng = np.random.default_rng(seed)
    indices = rng.permutation(x.shape[0])
    validation_count = int(round(x.shape[0] * max(0.0, min(0.8, validation_fraction))))
    validation_indices = indices[:validation_count]
    train_indices = indices[validation_count:]
    if train_indices.size == 0:
        train_indices = indices
        validation_indices = np.array([], dtype=int)

    x_train = x[train_indices]
    y_train = y[train_indices]
    input_dim = x.shape[1]
    output_dim = y.shape[1]
    hidden = max(1, int(hidden_units))
    w1 = rng.normal(0.0, np.sqrt(2.0 / input_dim), size=(input_dim, hidden))
    b1 = np.zeros(hidden)
    w2 = rng.normal(0.0, np.sqrt(2.0 / hidden), size=(hidden, output_dim))
    b2 = np.zeros(output_dim)
    m = [np.zeros_like(w1), np.zeros_like(b1), np.zeros_like(w2), np.zeros_like(b2)]
    v = [np.zeros_like(w1), np.zeros_like(b1), np.zeros_like(w2), np.zeros_like(b2)]
    beta1 = 0.9
    beta2 = 0.999
    eps = 1e-8

    for step in range(1, int(epochs) + 1):
        pred, hidden_state = _forward(x_train, w1, b1, w2, b2)
        d_pred = 2.0 * (pred - y_train) / max(1, pred.size)
        grad_w2 = hidden_state.T @ d_pred + weight_decay * w2
        grad_b2 = np.sum(d_pred, axis=0)
        d_hidden = (d_pred @ w2.T) * (1.0 - hidden_state**2)
        grad_w1 = x_train.T @ d_hidden + weight_decay * w1
        grad_b1 = np.sum(d_hidden, axis=0)
        grads = [grad_w1, grad_b1, grad_w2, grad_b2]
        params = [w1, b1, w2, b2]
        for index, grad in enumerate(grads):
            m[index] = beta1 * m[index] + (1.0 - beta1) * grad
            v[index] = beta2 * v[index] + (1.0 - beta2) * (grad * grad)
            m_hat = m[index] / (1.0 - beta1**step)
            v_hat = v[index] / (1.0 - beta2**step)
            params[index] -= learning_rate * m_hat / (np.sqrt(v_hat) + eps)

    train_pred, _hidden = _forward(x_train, w1, b1, w2, b2)
    train_mse = float(np.mean((train_pred - y_train) ** 2))
    validation_mse = None
    if validation_indices.size:
        validation_pred, _hidden = _forward(x[validation_indices], w1, b1, w2, b2)
        validation_mse = float(np.mean((validation_pred - y[validation_indices]) ** 2))

    return {
        "schema": "rad-sim.lock-coordinate-mlp.v1",
        "strandCells": dataset.strand_cells,
        "strandAxis": strand_axis,
        "inputDimension": int(input_dim),
        "outputDimension": int(output_dim),
        "hiddenUnits": int(hidden),
        "epochs": int(epochs),
        "learningRate": float(learning_rate),
        "weightDecay": float(weight_decay),
        "seed": int(seed),
        "xMean": x_mean,
        "xScale": x_scale,
        "yMean": y_mean,
        "yScale": y_scale,
        "weights": {
            "w1": w1,
            "b1": b1,
            "w2": w2,
            "b2": b2,
        },
        "metrics": {
            "sampleCount": int(x_raw.shape[0]),
            "trainCount": int(train_indices.size),
            "validationCount": int(validation_indices.size),
            "trainMseStandardized": train_mse,
            "validationMseStandardized": validation_mse,
        },
        "claimLabels": {
            "model": "small empirical MLP from lock operator state to 3D cell coordinates",
            "status": "working surrogate scaffold, not validated mechanics",
            "nextEvidence": "requires larger train/validation dataset across locks, actuators, boundaries, and sheets",
        },
    }


def predict_lock_coordinate_mlp(
    model: dict[str, object],
    lock_cells: Iterable[int],
    lock_angles: Iterable[int] | None = None,
    *,
    state_index: int = 1,
    realized_cells: Iterable[int] | None = None,
) -> dict[str, object]:
    cells = tuple(int(cell) for cell in lock_cells)
    angles = tuple(int(angle) for angle in (lock_angles or (30,) * len(cells)))
    strand_cells = int(model["strandCells"])
    failed = _consecutive_failed_lock_cells(cells)
    realized = tuple(realized_cells) if realized_cells is not None else realized_lock_cells(cells, failed)
    feature = lock_feature_vector(
        cells,
        angles,
        strand_cells=strand_cells,
        state_index=state_index,
        realized_cells=realized,
    )
    x = (feature - np.asarray(model["xMean"], dtype=float)) / np.asarray(model["xScale"], dtype=float)
    weights = model["weights"]
    prediction, _hidden = _forward(
        x.reshape(1, -1),
        np.asarray(weights["w1"], dtype=float),
        np.asarray(weights["b1"], dtype=float),
        np.asarray(weights["w2"], dtype=float),
        np.asarray(weights["b2"], dtype=float),
    )
    y = prediction.reshape(-1) * np.asarray(model["yScale"], dtype=float) + np.asarray(model["yMean"], dtype=float)
    return {
        "schema": "rad-sim.lock-coordinate-mlp-prediction.v1",
        "lockCells": cells,
        "lockAngles": angles,
        "stateIndex": int(state_index),
        "failedLockCells": failed,
        "realizedLockCells": realized,
        "cellCoordinates": y.reshape(strand_cells, 3),
        "claimLabels": {
            "prediction": "neural surrogate coordinate prediction",
            "status": model.get("claimLabels", {}).get("status", "empirical surrogate"),
        },
    }


def export_lock_coordinate_mlp_json(model: dict[str, object]) -> str:
    return json.dumps(model, indent=2, default=_json_default)
