from __future__ import annotations

import argparse
import json
from pathlib import Path

from .models import LatticeConfig
from .two_cell_bench import (
    TwoCellBenchControls,
    export_cad_rad_cell_archive_audit_json,
    export_cad_rad_cell_layout_json,
    export_cad_rad_cell_reference_profile_json,
    export_one_cell_cad_export_audit_json,
    export_two_cell_actuation_sweep_csv,
    export_two_cell_actuation_sweep_json,
    export_two_cell_backlash_sweep_csv,
    export_two_cell_bench_json,
    export_two_cell_contact_phase_map_csv,
    export_two_cell_contact_phase_map_json,
    export_two_cell_connector_contact_csv,
    export_two_cell_connector_contact_json,
    export_two_cell_connector_contact_sweep_csv,
    export_two_cell_connector_contact_sweep_json,
    export_two_cell_connector_measurement_template_csv,
    export_two_cell_connector_measurement_template_json,
    export_two_cell_cad_contact_decomposition_csv,
    export_two_cell_cad_contact_decomposition_json,
    export_two_cell_external_results_template_csv,
    export_two_cell_external_results_template_json,
    export_two_cell_external_fidelity_matrix_manifest_csv,
    export_two_cell_external_fidelity_matrix_manifest_json,
    export_two_cell_external_fidelity_mjcf_run_json,
    export_two_cell_exact_contact_handoff_plan_csv,
    export_two_cell_exact_contact_handoff_plan_json,
    export_two_cell_fidelity_matrix_measurement_template_csv,
    export_two_cell_fidelity_matrix_measurement_template_json,
    export_two_cell_measurement_template_csv,
    export_two_cell_mjcf_proxy_json,
    export_two_cell_mjcf_xml,
    export_two_cell_mujoco_proxy_run_json,
    export_two_cell_physical_fidelity_matrix_csv,
    export_two_cell_physical_fidelity_matrix_json,
    export_two_cell_radius_backlash_phase_diagram_csv,
    export_two_cell_radius_backlash_phase_diagram_json,
    export_two_cell_radius_backlash_transition_report_csv,
    export_two_cell_radius_backlash_transition_report_json,
    export_two_cell_physical_response_atlas_csv,
    export_two_cell_physical_response_atlas_json,
    export_two_cell_physical_test_packet_json,
    export_two_cell_physical_simulation_suite_csv,
    export_two_cell_physical_simulation_suite_json,
    export_two_cell_physics_validation_report_csv,
    export_two_cell_physics_validation_report_json,
    export_two_cell_quasistatic_json,
    export_two_cell_quasistatic_sweep_csv,
    export_two_cell_quasistatic_sweep_json,
    export_two_cell_segmented_cad_intake_templates_csv,
    export_two_cell_segmented_cad_intake_templates_json,
    export_two_cell_segmented_cad_intake_validation_csv,
    export_two_cell_segmented_cad_intake_validation_json,
    export_two_cell_segmented_cad_readiness_csv,
    export_two_cell_segmented_cad_readiness_json,
    two_cell_mjcf_proxy_report,
    two_cell_segmented_cad_intake_templates,
)


def write_two_cell_bench_packet_artifacts(
    out: str | Path,
    *,
    backlash: float = 0.1,
    pin_radius: float = 0.18,
    hole_radius: float = 0.225,
    alpha_command: float = -0.35,
    z_command: float = 0.35,
    hole_sweep_max: float = 0.5,
    hole_sweep_steps: int = 9,
) -> dict[str, Path]:
    output_dir = Path(out)
    output_dir.mkdir(parents=True, exist_ok=True)
    config = LatticeConfig(
        rows=1,
        cols=2,
        backlash=backlash,
        pin_radius=pin_radius,
        hole_radius=hole_radius,
    )
    controls = TwoCellBenchControls(
        alpha_command=alpha_command,
        z_command=z_command,
        hole_sweep_max=hole_sweep_max,
        hole_sweep_steps=hole_sweep_steps,
    )
    external_mjcf_dir = output_dir / "external_fidelity_matrix_mjcf"
    atlas_mjcf_dir = output_dir / "response_atlas_priority_mjcf"
    files = {
        "cad_archive_audit": output_dir / "cad_rad_cell_archive_audit.json",
        "cad_layout": output_dir / "cad_rad_cell_layout.json",
        "cad_reference_profile": output_dir / "cad_rad_cell_reference_profile.json",
        "one_cell_cad_export_audit": output_dir / "one_cell_cad_export_audit.json",
        "segmented_cad_readiness": output_dir / "two_cell_segmented_cad_readiness.json",
        "segmented_cad_readiness_csv": output_dir / "two_cell_segmented_cad_readiness.csv",
        "segmented_cad_intake_templates": output_dir / "two_cell_segmented_cad_intake_templates.json",
        "segmented_cad_intake_templates_csv": output_dir / "two_cell_segmented_cad_intake_templates.csv",
        "segmented_cad_intake_validation": output_dir / "two_cell_segmented_cad_intake_validation.json",
        "segmented_cad_intake_validation_csv": output_dir / "two_cell_segmented_cad_intake_validation.csv",
        "segmented_cad_template_dir": output_dir / "segmented_cad_intake_templates",
        "bench": output_dir / "two_cell_bench.json",
        "sweep": output_dir / "two_cell_backlash_sweep.csv",
        "actuation_sweep": output_dir / "two_cell_actuation_sweep.csv",
        "actuation_report": output_dir / "two_cell_actuation_sweep.json",
        "mjcf_xml": output_dir / "two_cell_cad_proxy.xml",
        "mjcf_report": output_dir / "two_cell_cad_proxy.json",
        "mujoco_run": output_dir / "two_cell_mujoco_proxy_run.json",
        "physical_suite": output_dir / "two_cell_physical_simulation_suite.json",
        "physical_suite_csv": output_dir / "two_cell_physical_simulation_suite.csv",
        "physical_fidelity_matrix": output_dir / "two_cell_physical_fidelity_matrix.json",
        "physical_fidelity_matrix_csv": output_dir / "two_cell_physical_fidelity_matrix.csv",
        "physical_response_atlas": output_dir / "two_cell_physical_response_atlas.json",
        "physical_response_atlas_csv": output_dir / "two_cell_physical_response_atlas.csv",
        "contact_phase_map": output_dir / "two_cell_contact_phase_map.json",
        "contact_phase_map_csv": output_dir / "two_cell_contact_phase_map.csv",
        "radius_backlash_phase_diagram": output_dir / "two_cell_radius_backlash_phase_diagram.json",
        "radius_backlash_phase_diagram_csv": output_dir / "two_cell_radius_backlash_phase_diagram.csv",
        "radius_backlash_transition_report": output_dir / "two_cell_radius_backlash_transition_report.json",
        "radius_backlash_transition_report_csv": output_dir / "two_cell_radius_backlash_transition_report.csv",
        "cad_contact_decomposition": output_dir / "two_cell_cad_contact_decomposition.json",
        "cad_contact_decomposition_csv": output_dir / "two_cell_cad_contact_decomposition.csv",
        "exact_contact_handoff_plan": output_dir / "two_cell_exact_contact_handoff_plan.json",
        "exact_contact_handoff_plan_csv": output_dir / "two_cell_exact_contact_handoff_plan.csv",
        "fidelity_matrix_measurement_template": output_dir / "two_cell_fidelity_matrix_measurement_template.json",
        "fidelity_matrix_measurement_template_csv": output_dir / "two_cell_fidelity_matrix_measurement_template.csv",
        "external_fidelity_matrix_manifest": output_dir / "two_cell_external_fidelity_matrix_manifest.json",
        "external_fidelity_matrix_manifest_csv": output_dir / "two_cell_external_fidelity_matrix_manifest.csv",
        "external_fidelity_matrix_mjcf_dir": external_mjcf_dir,
        "external_fidelity_matrix_mjcf_index": output_dir / "two_cell_external_fidelity_matrix_mjcf_index.json",
        "external_fidelity_matrix_mjcf_run": output_dir / "two_cell_external_fidelity_matrix_mjcf_run.json",
        "response_atlas_priority_mjcf_dir": atlas_mjcf_dir,
        "response_atlas_priority_mjcf_index": output_dir / "two_cell_response_atlas_priority_mjcf_index.json",
        "response_atlas_priority_mjcf_run": output_dir / "two_cell_response_atlas_priority_mjcf_run.json",
        "quasistatic": output_dir / "two_cell_quasistatic.json",
        "quasistatic_sweep": output_dir / "two_cell_quasistatic_sweep.csv",
        "quasistatic_report": output_dir / "two_cell_quasistatic_sweep.json",
        "connector_contact": output_dir / "two_cell_connector_contact.json",
        "connector_contact_csv": output_dir / "two_cell_connector_contact.csv",
        "connector_contact_sweep": output_dir / "two_cell_connector_contact_sweep.json",
        "connector_contact_sweep_csv": output_dir / "two_cell_connector_contact_sweep.csv",
        "connector_measurement_template": output_dir / "two_cell_connector_measurement_template.json",
        "connector_measurement_template_csv": output_dir / "two_cell_connector_measurement_template.csv",
        "external_template": output_dir / "two_cell_external_results_template.csv",
        "external_template_json": output_dir / "two_cell_external_results_template.json",
        "physics_validation": output_dir / "two_cell_physics_validation_report.json",
        "physics_validation_csv": output_dir / "two_cell_physics_validation_report.csv",
        "packet": output_dir / "two_cell_physical_test_packet.json",
        "template": output_dir / "two_cell_measurement_template.csv",
    }
    files["cad_archive_audit"].write_text(export_cad_rad_cell_archive_audit_json(), encoding="utf-8")
    files["cad_layout"].write_text(export_cad_rad_cell_layout_json(config), encoding="utf-8")
    files["cad_reference_profile"].write_text(export_cad_rad_cell_reference_profile_json(config), encoding="utf-8")
    files["one_cell_cad_export_audit"].write_text(export_one_cell_cad_export_audit_json(config), encoding="utf-8")
    files["segmented_cad_readiness"].write_text(
        export_two_cell_segmented_cad_readiness_json(config, controls),
        encoding="utf-8",
    )
    files["segmented_cad_readiness_csv"].write_text(
        export_two_cell_segmented_cad_readiness_csv(config, controls),
        encoding="utf-8",
    )
    files["segmented_cad_intake_templates"].write_text(
        export_two_cell_segmented_cad_intake_templates_json(config, controls),
        encoding="utf-8",
    )
    files["segmented_cad_intake_templates_csv"].write_text(
        export_two_cell_segmented_cad_intake_templates_csv(config, controls),
        encoding="utf-8",
    )
    files["segmented_cad_intake_validation"].write_text(
        export_two_cell_segmented_cad_intake_validation_json(config, controls),
        encoding="utf-8",
    )
    files["segmented_cad_intake_validation_csv"].write_text(
        export_two_cell_segmented_cad_intake_validation_csv(config, controls),
        encoding="utf-8",
    )
    files["segmented_cad_template_dir"].mkdir(parents=True, exist_ok=True)
    for item in two_cell_segmented_cad_intake_templates(config, controls)["templates"]:
        relative = str(item["path"]).replace("assets/cad/segmented/", "", 1)
        target = files["segmented_cad_template_dir"] / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(str(item["content"]), encoding="utf-8")
    files["bench"].write_text(export_two_cell_bench_json(config, controls), encoding="utf-8")
    files["sweep"].write_text(export_two_cell_backlash_sweep_csv(config, controls), encoding="utf-8")
    files["actuation_sweep"].write_text(export_two_cell_actuation_sweep_csv(config, controls), encoding="utf-8")
    files["actuation_report"].write_text(export_two_cell_actuation_sweep_json(config, controls), encoding="utf-8")
    files["mjcf_xml"].write_text(export_two_cell_mjcf_xml(config, controls), encoding="utf-8")
    files["mjcf_report"].write_text(export_two_cell_mjcf_proxy_json(config, controls), encoding="utf-8")
    files["mujoco_run"].write_text(export_two_cell_mujoco_proxy_run_json(config, controls), encoding="utf-8")
    files["physical_suite"].write_text(export_two_cell_physical_simulation_suite_json(config, controls), encoding="utf-8")
    files["physical_suite_csv"].write_text(export_two_cell_physical_simulation_suite_csv(config, controls), encoding="utf-8")
    files["physical_fidelity_matrix"].write_text(
        export_two_cell_physical_fidelity_matrix_json(config, controls),
        encoding="utf-8",
    )
    files["physical_fidelity_matrix_csv"].write_text(
        export_two_cell_physical_fidelity_matrix_csv(config, controls),
        encoding="utf-8",
    )
    files["physical_response_atlas"].write_text(
        export_two_cell_physical_response_atlas_json(config, controls),
        encoding="utf-8",
    )
    files["physical_response_atlas_csv"].write_text(
        export_two_cell_physical_response_atlas_csv(config, controls),
        encoding="utf-8",
    )
    files["contact_phase_map"].write_text(
        export_two_cell_contact_phase_map_json(config, controls),
        encoding="utf-8",
    )
    files["contact_phase_map_csv"].write_text(
        export_two_cell_contact_phase_map_csv(config, controls),
        encoding="utf-8",
    )
    files["radius_backlash_phase_diagram"].write_text(
        export_two_cell_radius_backlash_phase_diagram_json(config, controls),
        encoding="utf-8",
    )
    files["radius_backlash_phase_diagram_csv"].write_text(
        export_two_cell_radius_backlash_phase_diagram_csv(config, controls),
        encoding="utf-8",
    )
    files["radius_backlash_transition_report"].write_text(
        export_two_cell_radius_backlash_transition_report_json(config, controls),
        encoding="utf-8",
    )
    files["radius_backlash_transition_report_csv"].write_text(
        export_two_cell_radius_backlash_transition_report_csv(config, controls),
        encoding="utf-8",
    )
    files["cad_contact_decomposition"].write_text(
        export_two_cell_cad_contact_decomposition_json(config, controls),
        encoding="utf-8",
    )
    files["cad_contact_decomposition_csv"].write_text(
        export_two_cell_cad_contact_decomposition_csv(config, controls),
        encoding="utf-8",
    )
    files["exact_contact_handoff_plan"].write_text(
        export_two_cell_exact_contact_handoff_plan_json(config, controls),
        encoding="utf-8",
    )
    files["exact_contact_handoff_plan_csv"].write_text(
        export_two_cell_exact_contact_handoff_plan_csv(config, controls),
        encoding="utf-8",
    )
    files["fidelity_matrix_measurement_template"].write_text(
        export_two_cell_fidelity_matrix_measurement_template_json(config, controls),
        encoding="utf-8",
    )
    files["fidelity_matrix_measurement_template_csv"].write_text(
        export_two_cell_fidelity_matrix_measurement_template_csv(config, controls),
        encoding="utf-8",
    )
    files["external_fidelity_matrix_manifest"].write_text(
        export_two_cell_external_fidelity_matrix_manifest_json(config, controls),
        encoding="utf-8",
    )
    files["external_fidelity_matrix_manifest_csv"].write_text(
        export_two_cell_external_fidelity_matrix_manifest_csv(config, controls),
        encoding="utf-8",
    )
    manifest = json.loads(files["external_fidelity_matrix_manifest"].read_text(encoding="utf-8"))
    external_mjcf_dir.mkdir(parents=True, exist_ok=True)
    mjcf_index = []
    manifest_by_case_id = {str(case["caseId"]): case for case in manifest["caseRows"]}

    def write_case_proxy(case: dict, directory: Path, *, category: str | None = None) -> dict:
        case_controls = TwoCellBenchControls(
            alpha_command=float(case["alphaCommand"]),
            z_command=float(case["zCommand"]),
            hole_sweep_max=hole_sweep_max,
            hole_sweep_steps=hole_sweep_steps,
            left_position_locked=case["lockMode"] != "left_free",
            right_locked=case["lockMode"] == "right_state_locked",
            right_position_locked=case["lockMode"] == "right_position_locked",
        )
        case_config = LatticeConfig(
            rows=1,
            cols=2,
            backlash=float(case["backlash"]),
            pin_radius=float(case["pinRadius"]),
            hole_radius=float(case["holeRadius"]),
        )
        gravity_force = max(0.0, float(case["gravityForce"]))
        proxy = two_cell_mjcf_proxy_report(
            case_config,
            case_controls,
            hole_radius=float(case["holeRadius"]),
            pin_radius=float(case["pinRadius"]),
            lock_mode=str(case["lockMode"]),
            gravity=(0.0, 0.0, -9.81 * gravity_force),
        )
        xml_path = external_mjcf_dir / f"{case['caseId']}.xml"
        if category is not None:
            xml_path = directory / f"{case['caseId']}.xml"
        xml_path.write_text(str(proxy["xml"]), encoding="utf-8")
        out = {
            "caseId": case["caseId"],
            "matrixIndex": case["matrixIndex"],
            "actuationCase": case["actuationCase"],
            "path": str(xml_path.relative_to(output_dir)).replace("\\", "/"),
            "lockMode": case["lockMode"],
            "backlash": case["backlash"],
            "pinRadius": case["pinRadius"],
            "holeRadius": case["holeRadius"],
            "pinHoleClearanceMm": case["pinHoleClearanceMm"],
            "pinRadiusMm": proxy["scale"]["pinRadiusMm"],
            "holeRadiusMm": proxy["scale"]["holeRadiusMm"],
            "xyMmPerModelUnit": proxy["scale"]["xyMmPerModelUnit"],
            "zMmPerModelUnit": proxy["scale"]["zMmPerModelUnit"],
            "alphaCommand": case["alphaCommand"],
            "zCommand": case["zCommand"],
            "gravityForce": case["gravityForce"],
            "bodyCount": proxy["summary"]["bodyCount"],
            "actuatorCount": proxy["summary"]["actuatorCount"],
            "connectorContactPairCount": proxy["summary"]["connectorContactPairCount"],
        }
        if category is not None:
            out["atlasPriorityCategory"] = category
        return out

    for case in manifest["caseRows"]:
        mjcf_index.append(write_case_proxy(case, external_mjcf_dir))
    files["external_fidelity_matrix_mjcf_index"].write_text(
        json.dumps(
            {
                "schema": "rad-sim.two-cell-external-fidelity-matrix-mjcf-index.v1",
                "manifestSchema": manifest["schema"],
                "caseCount": len(mjcf_index),
                "directory": str(external_mjcf_dir.relative_to(output_dir)).replace("\\", "/"),
                "files": mjcf_index,
                "physicalAccuracyValidated": False,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    files["external_fidelity_matrix_mjcf_run"].write_text(
        export_two_cell_external_fidelity_mjcf_run_json(
            files["external_fidelity_matrix_mjcf_index"],
            engine_availability={"mujoco": False},
        ),
        encoding="utf-8",
    )
    atlas = json.loads(files["physical_response_atlas"].read_text(encoding="utf-8"))
    atlas_mjcf_dir.mkdir(parents=True, exist_ok=True)
    atlas_mjcf_index = []
    seen_atlas_cases = set()
    for priority_case in atlas["benchPriority"]["firstCasesToMeasure"]:
        case_id = str(priority_case.get("caseId") or f"FM_{int(priority_case['index']):03d}")
        if case_id in seen_atlas_cases:
            continue
        seen_atlas_cases.add(case_id)
        case = manifest_by_case_id.get(case_id)
        if case is None:
            case = manifest_by_case_id.get(f"FM_{int(priority_case['index']):03d}")
        if case is None:
            continue
        atlas_mjcf_index.append(
            write_case_proxy(
                case,
                atlas_mjcf_dir,
                category=str(priority_case.get("category") or "firstCasesToMeasure"),
            )
        )
    files["response_atlas_priority_mjcf_index"].write_text(
        json.dumps(
            {
                "schema": "rad-sim.two-cell-response-atlas-priority-mjcf-index.v1",
                "atlasSchema": atlas["schema"],
                "manifestSchema": manifest["schema"],
                "sourcePriority": "physicalResponseAtlas.benchPriority.firstCasesToMeasure",
                "caseCount": len(atlas_mjcf_index),
                "directory": str(atlas_mjcf_dir.relative_to(output_dir)).replace("\\", "/"),
                "files": atlas_mjcf_index,
                "physicalAccuracyValidated": False,
                "claimBoundary": {
                    "allowedClaim": "ranked reduced-proxy MJCF cases for first external or bench tests",
                    "blockedClaim": "exact real-cell mechanics until these cases are run against segmented CAD/contact or bench tracking",
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    files["response_atlas_priority_mjcf_run"].write_text(
        export_two_cell_external_fidelity_mjcf_run_json(
            files["response_atlas_priority_mjcf_index"],
            engine_availability={"mujoco": False},
        ),
        encoding="utf-8",
    )
    files["quasistatic"].write_text(export_two_cell_quasistatic_json(config, controls), encoding="utf-8")
    files["quasistatic_sweep"].write_text(export_two_cell_quasistatic_sweep_csv(config, controls), encoding="utf-8")
    files["quasistatic_report"].write_text(export_two_cell_quasistatic_sweep_json(config, controls), encoding="utf-8")
    files["connector_contact"].write_text(export_two_cell_connector_contact_json(config, controls), encoding="utf-8")
    files["connector_contact_csv"].write_text(export_two_cell_connector_contact_csv(config, controls), encoding="utf-8")
    files["connector_contact_sweep"].write_text(
        export_two_cell_connector_contact_sweep_json(config, controls),
        encoding="utf-8",
    )
    files["connector_contact_sweep_csv"].write_text(
        export_two_cell_connector_contact_sweep_csv(config, controls),
        encoding="utf-8",
    )
    files["connector_measurement_template"].write_text(
        export_two_cell_connector_measurement_template_json(config, controls),
        encoding="utf-8",
    )
    files["connector_measurement_template_csv"].write_text(
        export_two_cell_connector_measurement_template_csv(config, controls),
        encoding="utf-8",
    )
    files["external_template"].write_text(export_two_cell_external_results_template_csv(config, controls), encoding="utf-8")
    files["external_template_json"].write_text(
        export_two_cell_external_results_template_json(config, controls),
        encoding="utf-8",
    )
    files["physics_validation"].write_text(
        export_two_cell_physics_validation_report_json(config, controls),
        encoding="utf-8",
    )
    files["physics_validation_csv"].write_text(
        export_two_cell_physics_validation_report_csv(config, controls),
        encoding="utf-8",
    )
    files["packet"].write_text(
        export_two_cell_physical_test_packet_json(config, controls),
        encoding="utf-8",
    )
    files["template"].write_text(
        export_two_cell_measurement_template_csv(config, controls),
        encoding="utf-8",
    )
    return files


def main() -> None:
    parser = argparse.ArgumentParser(description="Export RAD two-cell physical bench packet artifacts.")
    parser.add_argument("--out", default="outputs/two_cell_bench_packet")
    parser.add_argument("--backlash", type=float, default=0.1)
    parser.add_argument("--pin-radius", type=float, default=0.18)
    parser.add_argument("--hole-radius", type=float, default=0.225)
    parser.add_argument("--alpha-command", type=float, default=-0.35)
    parser.add_argument("--z-command", type=float, default=0.35)
    parser.add_argument("--hole-sweep-max", type=float, default=0.5)
    parser.add_argument("--hole-sweep-steps", type=int, default=9)
    args = parser.parse_args()
    files = write_two_cell_bench_packet_artifacts(
        args.out,
        backlash=args.backlash,
        pin_radius=args.pin_radius,
        hole_radius=args.hole_radius,
        alpha_command=args.alpha_command,
        z_command=args.z_command,
        hole_sweep_max=args.hole_sweep_max,
        hole_sweep_steps=args.hole_sweep_steps,
    )
    for label, path in files.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
