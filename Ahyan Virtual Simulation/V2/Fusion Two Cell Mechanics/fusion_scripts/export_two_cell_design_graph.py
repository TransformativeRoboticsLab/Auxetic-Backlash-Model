"""Fusion script for exporting a RAD two-cell design graph.

Run inside Autodesk Fusion with the two-connected-cell design open:

Scripts and Add-Ins -> Scripts -> create/select this script -> Run

The script writes a JSON design graph and a STEP file into the project exports
folder. It is intentionally conservative: it extracts the CAD graph and raw
mechanism evidence needed by the simulator, but it does not claim that the
simulator equations are physically validated until CAD pose samples and bench
data are compared.
"""

import datetime
import json
import os
import traceback

import adsk.core
import adsk.fusion


SCHEMA_VERSION = "rad.two_cell_fusion_export.v1"
DEFAULT_OUTPUT_DIR = (
    r"C:\Users\ahyan\OneDrive\Documents\TRL\Ahyan Virtual Simulation"
    r"\Fusion Two Cell Mechanics\exports"
)
EXPORT_BASENAME = "two_cell"
ACTUATION_NAME_HINTS = (
    "actuat",
    "alpha",
    "theta",
    "angle",
    "drive",
    "motion",
    "open",
    "dilate",
)


def matrix_to_list(matrix):
    try:
        return list(matrix.asArray())
    except Exception:
        return []


def matrix_translation(matrix_values):
    if len(matrix_values) >= 12:
        return [
            safe_float(matrix_values[3]),
            safe_float(matrix_values[7]),
            safe_float(matrix_values[11]),
        ]
    return []


def point_to_list(point):
    try:
        return [float(point.x), float(point.y), float(point.z)]
    except Exception:
        return []


def vector_to_list(vector):
    try:
        return [float(vector.x), float(vector.y), float(vector.z)]
    except Exception:
        return []


def box_to_dict(box):
    try:
        return {
            "min": point_to_list(box.minPoint),
            "max": point_to_list(box.maxPoint),
        }
    except Exception:
        return {}


def safe_float(value):
    try:
        return float(value)
    except Exception:
        return None


def safe_attr(obj, name, default=None):
    try:
        return getattr(obj, name)
    except Exception:
        return default


def obj_type(obj):
    try:
        return obj.objectType
    except Exception:
        return type(obj).__name__ if obj is not None else ""


def compact_value(value):
    if value is None:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    if hasattr(value, "asArray"):
        return matrix_to_list(value)
    if all(hasattr(value, name) for name in ("x", "y", "z")):
        return point_to_list(value)
    return str(value)


def snapshot_attrs(obj, attrs):
    data = {}
    for attr in attrs:
        value = safe_attr(obj, attr)
        if value is not None:
            data[attr] = compact_value(value)
    return data


def safe_name(name):
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in name)


def get_occurrence_name(occurrence):
    if not occurrence:
        return ""
    try:
        return occurrence.fullPathName
    except Exception:
        return occurrence.name


def joint_motion_type(joint):
    try:
        return joint.jointMotion.objectType
    except Exception:
        return "unknown"


def parameter_record(param):
    if not param:
        return {}
    name = safe_attr(param, "name", "")
    return {
        "name": name,
        "source": obj_type(param),
        "owner": safe_attr(safe_attr(param, "parentComponent"), "name", ""),
        "expression": safe_attr(param, "expression", ""),
        "value": safe_float(safe_attr(param, "value")),
        "unit": safe_attr(param, "unit", ""),
        "comment": safe_attr(param, "comment", ""),
        "role": "actuation_candidate"
        if any(hint in name.lower() for hint in ACTUATION_NAME_HINTS)
        else "",
    }


def geometry_frame(geometry):
    if not geometry:
        return {}
    return {
        "object_type": obj_type(geometry),
        "origin": point_to_list(safe_attr(geometry, "origin")),
        "primary_axis": vector_to_list(safe_attr(geometry, "primaryAxisVector")),
        "secondary_axis": vector_to_list(safe_attr(geometry, "secondaryAxisVector")),
        "third_axis": vector_to_list(safe_attr(geometry, "thirdAxisVector")),
    }


def motion_limits(motion):
    data = snapshot_attrs(
        motion,
        (
            "jointType",
            "rotationAxis",
            "slideDirection",
            "rotationValue",
            "slideValue",
            "pitchValue",
            "pitch",
            "offset",
        ),
    )
    rotation_axis_vector = vector_to_list(safe_attr(motion, "rotationAxisVector"))
    slide_direction_vector = vector_to_list(safe_attr(motion, "slideDirectionVector"))
    if rotation_axis_vector:
        data["rotationAxisVector"] = rotation_axis_vector
    if slide_direction_vector:
        data["slideDirectionVector"] = slide_direction_vector
    for name in (
        "rotationLimits",
        "slideLimits",
        "pitchLimits",
        "customLimits",
    ):
        limits = safe_attr(motion, name)
        if limits:
            data[name] = snapshot_attrs(
                limits,
                (
                    "isMinimumValueEnabled",
                    "isMaximumValueEnabled",
                    "minimumValue",
                    "maximumValue",
                    "restValue",
                ),
            )
    return data


def collect_parameters(design):
    parameters = []
    seen = set()
    collections = (
        ("user", safe_attr(design, "userParameters", [])),
        ("all", safe_attr(design, "allParameters", [])),
    )
    for source, collection in collections:
        try:
            iterable = [collection.item(i) for i in range(collection.count)]
        except Exception:
            iterable = collection
        for param in iterable:
            record = parameter_record(param)
            name = record.get("name", "")
            owner = record.get("owner", "")
            key = (source, owner, name)
            if key in seen:
                continue
            seen.add(key)
            record["source"] = source
            parameters.append(record)
    return parameters


def body_record(body):
    record = {
        "name": safe_attr(body, "name", ""),
        "role": "",
        "object_type": obj_type(body),
        "entity_token": safe_attr(body, "entityToken", ""),
        "bounding_box": box_to_dict(safe_attr(body, "boundingBox")),
        "is_solid": bool(safe_attr(body, "isSolid", False)),
        "face_count": 0,
        "edge_count": 0,
        "cylindrical_faces": [],
    }
    try:
        record["face_count"] = body.faces.count
        record["edge_count"] = body.edges.count
    except Exception:
        pass
    try:
        physical = body.physicalProperties
        record["physical_properties"] = snapshot_attrs(
            physical,
            (
                "area",
                "volume",
                "mass",
                "density",
                "centerOfMass",
            ),
        )
    except Exception:
        pass
    try:
        for face in body.faces:
            geometry = face.geometry
            if "Cylinder" in obj_type(geometry):
                record["cylindrical_faces"].append(
                    {
                        "temp_id": safe_attr(face, "tempId", ""),
                        "radius": safe_float(safe_attr(geometry, "radius")),
                        "origin": point_to_list(safe_attr(geometry, "origin")),
                        "axis": vector_to_list(safe_attr(geometry, "axis")),
                    }
                )
    except Exception:
        pass
    return record


def collect_occurrences(root_component):
    occurrences = []
    for occurrence in root_component.allOccurrences:
        bodies = []
        try:
            for body in occurrence.bRepBodies:
                bodies.append({"name": body.name, "role": ""})
        except Exception:
            pass

        occurrences.append(
            {
                "name": occurrence.name,
                "full_path": get_occurrence_name(occurrence),
                "component": occurrence.component.name if occurrence.component else "",
                "role": "",
                "transform": matrix_to_list(occurrence.transform2),
                "bounding_box": box_to_dict(safe_attr(occurrence, "boundingBox")),
                "bodies": bodies,
            }
        )
    return occurrences


def collect_rigid_groups(root_component):
    groups = []
    for collection_name in ("allRigidGroups", "rigidGroups"):
        collection = safe_attr(root_component, collection_name)
        if not collection:
            continue
        try:
            iterable = [collection.item(i) for i in range(collection.count)]
        except Exception:
            iterable = collection
        for group in iterable:
            occurrences = []
            group_occurrences = safe_attr(group, "occurrences")
            if group_occurrences:
                try:
                    occurrences = [
                        get_occurrence_name(group_occurrences.item(i))
                        for i in range(group_occurrences.count)
                    ]
                except Exception:
                    pass
            groups.append(
                {
                    "name": safe_attr(group, "name", ""),
                    "source": collection_name,
                    "is_suppressed": bool(safe_attr(group, "isSuppressed", False)),
                    "occurrences": occurrences,
                }
            )
    return groups


def collect_joint_collection(collection, kind):
    joints = []
    if not collection:
        return joints
    try:
        iterable = [collection.item(i) for i in range(collection.count)]
    except Exception:
        iterable = collection
    for joint in iterable:
        geometry_one = None
        geometry_two = None
        try:
            geometry_one = joint.geometryOrOriginOne
            geometry_two = joint.geometryOrOriginTwo
        except Exception:
            pass

        motion = safe_attr(joint, "jointMotion")
        frame_one = geometry_frame(geometry_one)
        frame_two = geometry_frame(geometry_two)
        axis = frame_one.get("primary_axis") or frame_two.get("primary_axis") or []
        origin = frame_one.get("origin") or frame_two.get("origin") or []

        joints.append(
            {
                "name": joint.name,
                "kind": kind,
                "object_type": obj_type(joint),
                "joint_type": joint_motion_type(joint),
                "occurrence_one": get_occurrence_name(joint.occurrenceOne),
                "occurrence_two": get_occurrence_name(joint.occurrenceTwo),
                "axis": axis,
                "origin": origin,
                "limits": motion_limits(motion) if motion else {},
                "role": "",
                "geometry_one": frame_one,
                "geometry_two": frame_two,
                "is_suppressed": bool(safe_attr(joint, "isSuppressed", False)),
            }
        )
    return joints


def collect_joints(root_component):
    joints = []
    joints.extend(collect_joint_collection(safe_attr(root_component, "allJoints"), "joint"))
    joints.extend(
        collect_joint_collection(safe_attr(root_component, "allAsBuiltJoints"), "as_built_joint")
    )
    return joints


def collect_motion_links(root_component):
    links = []
    seen = set()
    components = [root_component]
    try:
        for occurrence in root_component.allOccurrences:
            if occurrence.component:
                components.append(occurrence.component)
    except Exception:
        pass

    for component in components:
        collection = safe_attr(component, "motionLinks")
        if not collection:
            continue
        try:
            iterable = [collection.item(i) for i in range(collection.count)]
        except Exception:
            iterable = collection
        for link in iterable:
            token = safe_attr(link, "entityToken", "") or safe_attr(link, "name", "")
            if token in seen:
                continue
            seen.add(token)
            value_one = parameter_record(safe_attr(link, "valueOne"))
            value_two = parameter_record(safe_attr(link, "valueTwo"))
            links.append(
                {
                    "name": safe_attr(link, "name", ""),
                    "object_type": obj_type(link),
                    "parent_component": safe_attr(safe_attr(link, "parentComponent"), "name", ""),
                    "joint_one": safe_attr(safe_attr(link, "jointOne"), "name", ""),
                    "joint_two": safe_attr(safe_attr(link, "jointTwo"), "name", ""),
                    "motion_one": compact_value(safe_attr(link, "motionOne")),
                    "motion_two": compact_value(safe_attr(link, "motionTwo")),
                    "value_one": value_one,
                    "value_two": value_two,
                    "is_reversed": bool(safe_attr(link, "isReversed", False)),
                    "is_suppressed": bool(safe_attr(link, "isSuppressed", False)),
                    "entity_token": token,
                    "error_or_warning": safe_attr(link, "errorOrWarningMessage", ""),
                }
            )
    return links


def collect_contact_candidates(occurrences):
    bodies = []
    for occurrence in occurrences:
        for body in occurrence.get("bodies", []):
            name = (occurrence.get("full_path", "") + "/" + body.get("name", "")).lower()
            cylinders = body.get("cylindrical_faces", [])
            role = ""
            if any(token in name for token in ("pin", "screw", "shaft", "bolt")):
                role = "pin_candidate"
            if any(token in name for token in ("hole", "slot", "bushing", "clearance")):
                role = "hole_candidate" if not role else role + "_and_hole_candidate"
            if role or cylinders:
                bodies.append(
                    {
                        "occurrence": occurrence.get("full_path", ""),
                        "body": body.get("name", ""),
                        "role_hint": role,
                        "cylindrical_faces": cylinders,
                        "bounding_box": body.get("bounding_box", {}),
                    }
                )
    candidates = []
    for i, one in enumerate(bodies):
        for two in bodies[i + 1 :]:
            role_text = one.get("role_hint", "") + " " + two.get("role_hint", "")
            if not role_text.strip() and not (
                one.get("cylindrical_faces") and two.get("cylindrical_faces")
            ):
                continue
            candidates.append(
                {
                    "one": one,
                    "two": two,
                    "reason": "name/cylindrical-face proximity candidate; validate in Fusion",
                }
            )
            if len(candidates) >= 300:
                return candidates
    return candidates


def occurrence_transforms(root_component):
    transforms = {}
    for occurrence in root_component.allOccurrences:
        transforms[get_occurrence_name(occurrence)] = matrix_to_list(occurrence.transform2)
    return transforms


def center_from_box(box):
    if not box or "min" not in box or "max" not in box:
        return []
    try:
        return [
            0.5 * (float(box["min"][0]) + float(box["max"][0])),
            0.5 * (float(box["min"][1]) + float(box["max"][1])),
            0.5 * (float(box["min"][2]) + float(box["max"][2])),
        ]
    except Exception:
        return []


def infer_cell_centers(occurrences):
    cells = []
    top_level = []
    for occurrence in occurrences:
        name = (occurrence.get("full_path") or occurrence.get("name") or "").lower()
        component = (occurrence.get("component") or "").lower()
        is_top_level = "+" not in (occurrence.get("full_path") or "")
        is_unit_cell = "unit cell" in name or "unit cell" in component
        if is_top_level and is_unit_cell:
            center = matrix_translation(occurrence.get("transform", []))
            if center:
                top_level.append(
                    {
                        "name": occurrence.get("full_path", occurrence.get("name", "")),
                        "center": center,
                        "source": "top_level_occurrence_origin",
                    }
                )
    if len(top_level) >= 2:
        return top_level

    for occurrence in occurrences:
        name = (occurrence.get("full_path") or occurrence.get("name") or "").lower()
        component = (occurrence.get("component") or "").lower()
        if "+" not in (occurrence.get("full_path") or "") and (
            "cell" in name or "cell" in component or "rad" in name or "rad" in component
        ):
            center = center_from_box(occurrence.get("bounding_box"))
            if center:
                cells.append(
                    {
                        "name": occurrence.get("full_path", occurrence.get("name", "")),
                        "center": center,
                        "source": "top_level_bounding_box",
                    }
                )
    return cells


def pitch_between(first, second):
    try:
        dx = first[0] - second[0]
        dy = first[1] - second[1]
        dz = first[2] - second[2]
        return (dx * dx + dy * dy + dz * dz) ** 0.5
    except Exception:
        return None


def capture_pose(root_component, occurrences, sample_id, actuation_value, label):
    if occurrences is None:
        occurrences = collect_occurrences(root_component)
    cells = infer_cell_centers(occurrences)
    left_center = cells[0]["center"] if len(cells) > 0 else []
    right_center = cells[1]["center"] if len(cells) > 1 else []
    return {
        "sample_id": sample_id,
        "label": label,
        "actuation": actuation_value,
        "cell_metrics": {
            "left_alpha": None,
            "right_alpha": None,
            "left_theta": None,
            "right_theta": None,
            "left_center": left_center,
            "right_center": right_center,
            "pitch": pitch_between(left_center, right_center),
            "inferred_cell_centers": cells,
        },
        "occurrence_transforms": occurrence_transforms(root_component),
        "tracked_points_world": {},
        "interference": None,
    }


def find_actuation_parameter(design):
    try:
        params = design.userParameters
        for i in range(params.count):
            param = params.item(i)
            if any(hint in param.name.lower() for hint in ACTUATION_NAME_HINTS):
                return param
    except Exception:
        pass
    return None


def find_limited_joint_drive(root_component):
    try:
        joints = [root_component.allJoints.item(i) for i in range(root_component.allJoints.count)]
    except Exception:
        try:
            joints = list(root_component.allJoints)
        except Exception:
            joints = []

    candidates = []
    for joint in joints:
        motion = safe_attr(joint, "jointMotion")
        if not motion:
            continue
        for value_attr, limits_attr, unit in (
            ("rotationValue", "rotationLimits", "rad"),
            ("slideValue", "slideLimits", "cm"),
            ("pitchValue", "pitchLimits", "rad_or_cm_by_motion"),
        ):
            limits = safe_attr(motion, limits_attr)
            minimum = safe_float(safe_attr(limits, "minimumValue"))
            maximum = safe_float(safe_attr(limits, "maximumValue"))
            current = safe_float(safe_attr(motion, value_attr))
            has_min = bool(safe_attr(limits, "isMinimumValueEnabled", False))
            has_max = bool(safe_attr(limits, "isMaximumValueEnabled", False))
            if has_min and has_max and minimum is not None and maximum is not None and maximum > minimum:
                candidates.append(
                    {
                        "joint": joint,
                        "motion": motion,
                        "value_attr": value_attr,
                        "limits_attr": limits_attr,
                        "unit": unit,
                        "minimum": minimum,
                        "maximum": maximum,
                        "current": current,
                        "priority": 0 if value_attr == "rotationValue" else 1,
                    }
                )
    if not candidates:
        return None, []
    candidates.sort(key=lambda item: (item["priority"], item["maximum"] - item["minimum"]))
    return candidates[0], [
        {
            "joint": safe_attr(item["joint"], "name", ""),
            "motion_type": joint_motion_type(item["joint"]),
            "value_attr": item["value_attr"],
            "limits_attr": item["limits_attr"],
            "unit": item["unit"],
            "minimum": item["minimum"],
            "maximum": item["maximum"],
            "current": item["current"],
        }
        for item in candidates
    ]


def expression_for_value(value, unit):
    if unit:
        return str(value) + " " + unit
    return str(value)


def collect_pose_samples(design, root_component, occurrences):
    samples = [capture_pose(root_component, occurrences, 0, 0.0, "current_pose")]
    parameter = find_actuation_parameter(design)
    if not parameter:
        joint_drive, candidates = find_limited_joint_drive(root_component)
        if not joint_drive:
            return samples, {
                "sampled": False,
                "reason": "No user parameter name matched actuation hints and no limited joint drive was found.",
            }

        motion = joint_drive["motion"]
        value_attr = joint_drive["value_attr"]
        original_value = safe_float(safe_attr(motion, value_attr))
        minimum = joint_drive["minimum"]
        maximum = joint_drive["maximum"]
        current = joint_drive["current"]
        midpoint = 0.5 * (minimum + maximum)
        values = [
            minimum,
            0.5 * (minimum + midpoint),
            current if current is not None else midpoint,
            0.5 * (midpoint + maximum),
            maximum,
        ]
        try:
            driven_samples = []
            for index, value in enumerate(values, start=1):
                setattr(motion, value_attr, value)
                adsk.doEvents()
                driven_samples.append(
                    capture_pose(
                        root_component,
                        None,
                        index,
                        value,
                        safe_attr(joint_drive["joint"], "name", "") + "." + value_attr,
                    )
                )
            samples.extend(driven_samples)
            return samples, {
                "sampled": True,
                "source": "limited_joint_drive",
                "joint": safe_attr(joint_drive["joint"], "name", ""),
                "motion_type": joint_motion_type(joint_drive["joint"]),
                "value_attr": value_attr,
                "unit": joint_drive["unit"],
                "values": values,
                "candidate_limited_joints": candidates,
                "restored_value": original_value,
            }
        except Exception as exc:
            return samples, {
                "sampled": False,
                "source": "limited_joint_drive",
                "joint": safe_attr(joint_drive["joint"], "name", ""),
                "reason": "Driving the limited joint failed; current pose was still exported.",
                "error": str(exc),
                "candidate_limited_joints": candidates,
            }
        finally:
            if original_value is not None:
                try:
                    setattr(motion, value_attr, original_value)
                    adsk.doEvents()
                except Exception:
                    pass

    original_expression = parameter.expression
    baseline = safe_float(parameter.value)
    if baseline is None:
        return samples, {
            "sampled": False,
            "parameter": parameter.name,
            "reason": "Matched parameter had no numeric value.",
        }

    unit = parameter.unit
    scale = max(abs(baseline), 1.0)
    values = [baseline - 0.2 * scale, baseline - 0.1 * scale, baseline, baseline + 0.1 * scale, baseline + 0.2 * scale]
    driven_samples = []
    try:
        for index, value in enumerate(values, start=1):
            parameter.expression = expression_for_value(value, unit)
            adsk.doEvents()
            driven_samples.append(
                capture_pose(root_component, None, index, value, parameter.name)
            )
        samples.extend(driven_samples)
        return samples, {
            "sampled": True,
            "parameter": parameter.name,
            "unit": unit,
            "values": values,
            "restored_expression": original_expression,
        }
    except Exception as exc:
        return samples, {
            "sampled": False,
            "parameter": parameter.name,
            "reason": "Pose driving failed; current pose was still exported.",
            "error": str(exc),
        }
    finally:
        try:
            parameter.expression = original_expression
            adsk.doEvents()
        except Exception:
            pass


def export_step(design, output_dir):
    step_path = os.path.join(output_dir, EXPORT_BASENAME + ".step")
    export_manager = design.exportManager
    step_options = export_manager.createSTEPExportOptions(step_path)
    success = export_manager.execute(step_options)
    return step_path, bool(success)


def run(_context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        product = app.activeProduct

        if not isinstance(product, adsk.fusion.Design):
            ui.messageBox("The active document is not a Fusion design.")
            return

        design = product
        root = design.rootComponent
        document = app.activeDocument
        design_name = document.name if document else "two_cell_design"
        output_dir = DEFAULT_OUTPUT_DIR
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        step_path, step_success = export_step(design, output_dir)
        occurrences = collect_occurrences(root)
        pose_samples, pose_sampling = collect_pose_samples(design, root, occurrences)

        data = {
            "schema_version": SCHEMA_VERSION,
            "source": {
                "design_name": design_name,
                "source_file": design_name,
                "export_timestamp": datetime.datetime.utcnow().isoformat() + "Z",
                "notes": "CAD graph export. Validate roles, contact pairs, and pose samples before deriving final equations.",
            },
            "units": {
                "length": design.unitsManager.defaultLengthUnits,
                "angle": "rad_internal",
            },
            "coordinate_frame": {
                "origin": "Fusion root component origin unless overwritten",
                "x_axis": "unset",
                "y_axis": "unset",
                "z_axis": "unset",
            },
            "parameters": collect_parameters(design),
            "occurrences": occurrences,
            "rigid_groups": collect_rigid_groups(root),
            "joints": collect_joints(root),
            "motion_links": collect_motion_links(root),
            "tracked_points": [],
            "pose_samples": pose_samples,
            "contact_candidates": collect_contact_candidates(occurrences),
            "pose_sampling": pose_sampling,
            "exports": {
                "step_path": step_path,
                "step_success": step_success,
            },
        }

        json_path = os.path.join(output_dir, EXPORT_BASENAME + "_design_graph.json")
        with open(json_path, "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2)

        ui.messageBox(
            "RAD two-cell export complete.\n\n"
            + "JSON:\n"
            + json_path
            + "\n\nSTEP:\n"
            + step_path
            + "\n\nSTEP success: "
            + str(step_success)
        )
    except Exception:
        if ui:
            ui.messageBox("Export failed:\n" + traceback.format_exc())


def stop(_context):
    pass
