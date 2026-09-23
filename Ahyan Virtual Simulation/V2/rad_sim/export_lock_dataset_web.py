from __future__ import annotations

import json
from pathlib import Path

from .lock_dataset import load_lock_dataset, lock_coordinate_dataset


DEFAULT_LOGS_PATH = Path(r"C:\Users\ahyan\OneDrive\Desktop\logs")


def browser_lock_dataset_payload(logs_path: str | Path = DEFAULT_LOGS_PATH) -> dict[str, object]:
    dataset = load_lock_dataset(logs_path)
    payload = lock_coordinate_dataset(dataset)
    payload["schema"] = "rad-sim.browser-lock-coordinate-dataset.v1"
    payload["runtime"] = "browser"
    payload["sourcePath"] = "processed-local-lock-logs"
    payload["markerSummaryFile"] = "marker_xyz_means_all_37_files.csv"
    return payload


def write_browser_lock_dataset_artifact(
    output_path: str | Path,
    *,
    logs_path: str | Path = DEFAULT_LOGS_PATH,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = browser_lock_dataset_payload(logs_path)
    text = (
        "window.RAD_LOCK_DATASET = "
        + json.dumps(payload, indent=2, default=_json_default)
        + ";\n"
    )
    output.write_text(text, encoding="utf-8")
    return output


def _json_default(value: object) -> object:
    if hasattr(value, "tolist"):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return list(value)
    raise TypeError(f"object is not JSON serializable: {type(value)!r}")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    output = root / "web" / "data" / "lock_dataset.js"
    written = write_browser_lock_dataset_artifact(output)
    print(written)


if __name__ == "__main__":
    main()
