"""Export compact two-cell external-fidelity evidence for the browser UI."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

from .two_cell_bench import CAD_RAD_CELL_REFERENCE


WEB_TWO_CELL_EXTERNAL_FIDELITY_SCHEMA = "rad-sim.web-two-cell-external-fidelity-summary.v1"

DEFAULT_BENCHMARK = Path(
    "outputs/two_cell_external_fidelity_mjcf_full_real_comparison/"
    "two_cell_external_fidelity_benchmark_summary.json"
)
DEFAULT_PROFILE = Path(
    "outputs/two_cell_external_fidelity_mjcf_full_real_correction/"
    "two_cell_external_fidelity_correction_profile.json"
)
DEFAULT_APPLICATION = Path(
    "outputs/two_cell_external_fidelity_mjcf_full_real_corrected_preview/"
    "two_cell_external_fidelity_correction_application.json"
)
DEFAULT_OUT_JSON = Path("web/data/two_cell_external_fidelity_summary.json")
DEFAULT_REPRESENTATIVE_CASE_IDS = (
    "FM_000",
    "FM_028",
    "FM_056",
    "FM_084",
    "FM_112",
    "FM_140",
    "FM_168",
    "FM_196",
    "FM_224",
    "FM_251",
)

REPRESENTATIVE_ROW_FIELDS = (
    "caseId",
    "matrixIndex",
    "connector",
    "leftSite",
    "rightSite",
    "actuationCase",
    "lockMode",
    "backlash",
    "alphaCommand",
    "zCommand",
    "gravityForce",
    "pinRadius",
    "holeRadius",
    "pinHoleClearanceMm",
    "predictedContactMode",
    "observedContactMode",
    "lockHeldObserved",
    "predictedRightCellZ",
    "correctedRightCellZ",
    "observedRightCellZ",
    "baselineResidualRightCellZ",
    "correctedResidualRightCellZ",
    "predictedVerticalSlipMm",
    "correctedVerticalSlipMm",
    "observedVerticalSlipMm",
    "baselineResidualVerticalSlipMm",
    "correctedResidualVerticalSlipMm",
    "predictedTotalSlipMm",
    "correctedTotalSlipMm",
    "observedTotalSlipMm",
)

SUMMARY_CORRECTION_FIELDS = (
    "observedRightCellZ",
    "observedVerticalSlipMm",
    "observedLateralSlipMm",
    "observedTotalSlipMm",
)


def _load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _subset_mapping(mapping: dict[str, Any], keys: Iterable[str]) -> dict[str, Any]:
    return {key: mapping[key] for key in keys if key in mapping}


def _compact_rows(rows: Iterable[dict[str, Any]], case_ids: Iterable[str]) -> list[dict[str, Any]]:
    wanted = set(str(case_id) for case_id in case_ids)
    selected = [row for row in rows if str(row.get("caseId", "")) in wanted]
    if not selected:
        selected = list(rows)[:24]
    return [_subset_mapping(dict(row), REPRESENTATIVE_ROW_FIELDS) for row in selected]


def build_two_cell_external_fidelity_web_summary(
    benchmark: dict[str, Any],
    profile: dict[str, Any],
    application: dict[str, Any],
    *,
    source_artifacts: dict[str, str] | None = None,
    representative_case_ids: Iterable[str] = DEFAULT_REPRESENTATIVE_CASE_IDS,
) -> dict[str, Any]:
    cad_reference = CAD_RAD_CELL_REFERENCE.to_dict()
    artifact_cad_reference = profile.get("cadReference") or application.get("cadReference")
    if isinstance(artifact_cad_reference, dict):
        cad_reference.update(artifact_cad_reference)
    evaluation = dict(application.get("evaluation", {}))
    comparison = dict(benchmark.get("comparisonSummary", {}))
    run_summary = dict(benchmark.get("runSummary", {}))
    application_summary = dict(application.get("summary", {}))
    correction_summary = dict(profile.get("summary", {}))
    remaining_evidence = sorted(
        set(
            [
                *list(benchmark.get("summary", {}).get("remainingEvidence", [])),
                *list(correction_summary.get("remainingEvidence", [])),
                *list(application_summary.get("remainingEvidence", [])),
            ]
        )
    )
    field_corrections = dict(profile.get("fieldCorrections", {}))
    return {
        "schema": WEB_TWO_CELL_EXTERNAL_FIDELITY_SCHEMA,
        "model": "browser-compact-two-cell-external-fidelity-evidence",
        "sourceArtifacts": source_artifacts or {},
        "cadReference": cad_reference,
        "runSummary": run_summary,
        "comparisonSummary": comparison,
        "correctionSummary": correction_summary,
        "applicationSummary": application_summary,
        "applicationEvaluation": evaluation,
        "proposedReducedProxyUpdates": benchmark.get("proposedReducedProxyUpdates", {}),
        "fieldCorrections": _subset_mapping(field_corrections, SUMMARY_CORRECTION_FIELDS),
        "axisSummaries": benchmark.get("axisSummaries", {}),
        "worstAxisResiduals": benchmark.get("worstAxisResiduals", []),
        "representativeRows": _compact_rows(application.get("rows", []), representative_case_ids),
        "summary": {
            "available": True,
            "externalRunComplete": bool(benchmark.get("summary", {}).get("externalRunComplete")),
            "coversFullFidelityMatrix": bool(benchmark.get("summary", {}).get("coversFullFidelityMatrix")),
            "caseCount": run_summary.get("caseCount") or evaluation.get("caseCount"),
            "connectorRowCount": run_summary.get("expectedConnectorMeasurementRowCount")
            or application_summary.get("rowCount"),
            "measurementRowCount": application_summary.get("measurementRowCount"),
            "observedScalarCount": evaluation.get("observedScalarCount")
            or comparison.get("observedScalarCount"),
            "baselineRms": evaluation.get("baselineRms"),
            "correctedRms": evaluation.get("correctedRms"),
            "improvement": evaluation.get("improvement"),
            "holdoutImprovement": correction_summary.get("holdoutImprovement"),
            "contactModeAccuracy": evaluation.get("contactModeAccuracy")
            if evaluation.get("contactModeAccuracy") is not None
            else comparison.get("contactModeAccuracy"),
            "lockHeldAccuracy": evaluation.get("lockHeldAccuracy")
            if evaluation.get("lockHeldAccuracy") is not None
            else comparison.get("lockHeldAccuracy"),
            "readyForProxyPreview": bool(application_summary.get("readyForProxyPreview")),
            "passesReducedModelTolerance": bool(
                benchmark.get("summary", {}).get("passesReducedModelTolerance")
            ),
            "physicalAccuracyValidated": False,
            "remainingEvidence": remaining_evidence,
            "remainingEvidenceCount": len(remaining_evidence),
        },
        "claimLabels": {
            "externalRun": "MuJoCo ran all generated proxy two-cell matrix cases.",
            "correctedPreview": "Affine correction improves the proxy preview but remains fitted to proxy data.",
            "physicalAccuracy": "false until segmented CAD contact and independent bench holdout pass.",
        },
    }


def write_two_cell_external_fidelity_web_summary(
    benchmark_path: str | Path = DEFAULT_BENCHMARK,
    profile_path: str | Path = DEFAULT_PROFILE,
    application_path: str | Path = DEFAULT_APPLICATION,
    out_json: str | Path = DEFAULT_OUT_JSON,
) -> dict[str, Path]:
    benchmark = _load_json(benchmark_path)
    profile = _load_json(profile_path)
    application = _load_json(application_path)
    out_json = Path(out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    payload = build_two_cell_external_fidelity_web_summary(
        benchmark,
        profile,
        application,
        source_artifacts={
            "benchmark": str(benchmark_path),
            "correctionProfile": str(profile_path),
            "correctedPreview": str(application_path),
        },
    )
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    out_js = out_json.with_suffix(".js")
    out_js.write_text(
        "(function () {\n"
        "  const payload = "
        + json.dumps(payload, indent=2)
        + ";\n"
        "  window.RAD_EXTERNAL_FIDELITY_SUMMARY = payload;\n"
        "  window.RAD = window.RAD || {};\n"
        "  window.RAD.EXTERNAL_FIDELITY_SUMMARY = payload;\n"
        "})();\n",
        encoding="utf-8",
    )
    return {"json": out_json, "js": out_js}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export compact two-cell external-fidelity evidence for web/index.html."
    )
    parser.add_argument("--benchmark", default=str(DEFAULT_BENCHMARK))
    parser.add_argument("--profile", default=str(DEFAULT_PROFILE))
    parser.add_argument("--application", default=str(DEFAULT_APPLICATION))
    parser.add_argument("--out", default=str(DEFAULT_OUT_JSON))
    args = parser.parse_args()
    files = write_two_cell_external_fidelity_web_summary(
        args.benchmark,
        args.profile,
        args.application,
        args.out,
    )
    print(f"wrote {files['json']}")
    print(f"wrote {files['js']}")


if __name__ == "__main__":
    main()
