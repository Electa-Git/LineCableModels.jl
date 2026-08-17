from __future__ import annotations

import itertools
import re
import shutil
import time
import traceback
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import mhi.pscad

from .util import (
    case_id,
    copy_project_tree,
    flat_snapshot,
    jsonable,
    project_namespace,
    sha256_file,
    sha256_json,
    slug,
    snapshot,
    utc_now,
    wait_for_artifacts,
    write_json,
)


@dataclass(frozen=True)
class LineKey:
    iid: int
    definition: str
    component_type: str
    name: str


@dataclass(frozen=True)
class Variant:
    resistivity_ohm_m: float | None
    formulations: dict[str, Any]
    reduction: str

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def all_aclines(project: Any) -> list[Any]:
    return project.find_all("TLine") + project.find_all("Cable")


def local_canvas(project: Any, line: Any) -> Any:
    try:
        return line.canvas()
    except ValueError:
        # Some official archives renamed the .pscx for V5 without rewriting
        # the scope prefix embedded in local definition references.
        return project.canvas(str(line.defn_name).split(":", 1)[-1])


def values(component: Any) -> dict[str, Any]:
    result = component.parameters()
    return result if isinstance(result, dict) else {}


def choice_values(component: Any, parameter: str) -> list[Any]:
    raw = component.range(parameter)
    if isinstance(raw, range):
        return list(raw)
    if isinstance(raw, (set, frozenset, list, tuple)):
        return sorted(raw, key=str)
    raise ValueError(f"{component.defn_name}.{parameter} is not a finite choice: {raw!r}")


def choose(available: Iterable[Any], *tokens: str) -> Any:
    items = list(available)
    for token in tokens:
        for item in items:
            if token.casefold() in str(item).casefold():
                return item
    raise ValueError(f"None of {tokens!r} appears in {items!r}")


def component_snapshot(component: Any, *, include_ranges: bool = False) -> dict[str, Any]:
    params = values(component)
    row: dict[str, Any] = {
        "iid": getattr(component, "iid", None),
        "type": type(component).__name__,
        "definition": getattr(component, "defn_name", None),
        "parameters": jsonable(params),
    }
    if include_ranges:
        ranges = {}
        selector_names = {
            "EarthForm",
            "EarthForm2",
            "EarthForm3",
            "GrRho",
            "Output",
            "Numf",
            "FDIS",
            "CPASS",
            "NFP",
            "DCenab",
            "DCCOR",
            "Interp1",
            "ElimGW",
            "LC",
            "elimp",
        }
        selector_names.update(name for name in params if name.casefold().startswith("elim"))
        for name in params.keys() & selector_names:
            try:
                candidate = component.range(name)
            except Exception:
                continue
            if isinstance(candidate, (range, set, frozenset, list, tuple)):
                ranges[name] = jsonable(candidate)
        row["ranges"] = ranges
    return row


def line_snapshot(project: Any, line: Any, *, include_ranges: bool = False) -> dict[str, Any]:
    canvas = local_canvas(project, line)
    return {
        "outer": component_snapshot(line, include_ranges=include_ranges),
        "internal_canvas": [
            component_snapshot(component, include_ranges=include_ranges)
            for component in canvas.components()
        ],
    }


def unique_lines(project: Any) -> list[LineKey]:
    result: list[LineKey] = []
    seen: set[str] = set()
    for line in all_aclines(project):
        definition = str(line.defn_name)
        key = definition.casefold()
        if key in seen:
            continue
        seen.add(key)
        params = values(line)
        result.append(
            LineKey(
                iid=int(line.iid),
                definition=definition,
                component_type=type(line).__name__,
                name=str(params.get("Name", params.get("name", definition.split(":")[-1]))),
            )
        )
    return result


def find_line(project: Any, key: LineKey) -> Any:
    for line in all_aclines(project):
        if int(line.iid) == key.iid:
            return line
    for line in all_aclines(project):
        if str(line.defn_name).casefold() == key.definition.casefold():
            return line
    raise LookupError(f"Line definition not found: {key.definition}")


def relevant_components(project: Any, line: Any) -> dict[str, Any]:
    components = local_canvas(project, line).components()
    ground = next(
        (
            component
            for component in components
            if str(component.defn_name).casefold().endswith(":line_ground")
            or {"EarthForm", "EarthForm2", "GRRES"}.issubset(values(component))
        ),
        None,
    )
    frequency = next(
        (
            component
            for component in components
            if str(component.defn_name).casefold().endswith(":line_frephase_options")
            or {"FS", "FE", "Numf", "Output"}.issubset(values(component))
        ),
        None,
    )
    defns = [str(component.defn_name).casefold() for component in components]
    parameter_sets = [values(component) for component in components]
    aerial = any("line_tower" in name for name in defns) or any(
        "ElimGW" in params and "NG" in params for params in parameter_sets
    )
    underground = any("cable_" in name or "ground_plane" in name for name in defns) or any(
        "LL" in params or "cnum" in params or any(re.fullmatch(r"LL\d+", key) for key in params)
        for params in parameter_sets
    )
    return {
        "all": components,
        "ground": ground,
        "frequency": frequency,
        "aerial": aerial,
        "underground": underground,
    }


def requested_choices(setting: Any, available: list[Any]) -> list[Any]:
    if setting == "available" or setting is None:
        return available
    requested = setting if isinstance(setting, list) else [setting]
    lookup = {str(item).casefold(): item for item in available}
    result = []
    for item in requested:
        key = str(item).casefold()
        if key not in lookup:
            raise ValueError(f"Requested choice {item!r} is not available: {available!r}")
        result.append(lookup[key])
    return result


def formulation_axes(info: dict[str, Any], config: dict[str, Any]) -> list[tuple[str, list[Any]]]:
    axis_config = config["axes"]["earth_formulation"]
    if not axis_config.get("enabled", True):
        return []
    ground = info["ground"]
    if ground is None:
        return []
    axes: list[tuple[str, list[Any]]] = []
    if info["aerial"]:
        available = choice_values(ground, "EarthForm2")
        axes.append(("EarthForm2", requested_choices(axis_config.get("aerial"), available)))
    if info["underground"]:
        available = choice_values(ground, "EarthForm")
        axes.append(("EarthForm", requested_choices(axis_config.get("underground"), available)))
    if info["aerial"] and info["underground"]:
        available = choice_values(ground, "EarthForm3")
        axes.append(("EarthForm3", requested_choices(axis_config.get("mixed"), available)))
    return axes


def has_reduction_control(info: dict[str, Any]) -> bool:
    for component in info["all"]:
        params = values(component)
        if "ElimGW" in params and int(params.get("NG", 0)) > 0:
            return True
        if "LC" in params or ("elimp" in params and bool(params.get("poins", True))):
            return True
        if any(re.fullmatch(r"LC\d+", name, flags=re.IGNORECASE) for name in params):
            return True
        if any(name.casefold().startswith("elim") for name in params):
            return True
    return False


def variants_for(info: dict[str, Any], config: dict[str, Any]) -> list[Variant]:
    if config["frequency_dependent_phase"].get("enabled", True) and info["frequency"] is None:
        return []
    rho_config = config["axes"]["earth_resistivity"]
    resistivities: list[float | None]
    if rho_config.get("enabled", True) and info["ground"] is not None:
        resistivities = [float(item) for item in rho_config.get("values_ohm_m", [])]
    else:
        resistivities = [None]

    axes = formulation_axes(info, config)
    if axes:
        formulations = [dict(zip((name for name, _ in axes), combination)) for combination in itertools.product(*(items for _, items in axes))]
    else:
        formulations = [{}]

    reduction_config = config["axes"]["kron_reduction"]
    if reduction_config.get("enabled", True) and has_reduction_control(info):
        reductions = [str(item) for item in reduction_config.get("states", ["retained", "reduced"])]
    else:
        reductions = ["not_applicable"]

    return [
        Variant(rho, formulation, reduction)
        for rho, formulation, reduction in itertools.product(resistivities, formulations, reductions)
    ]


def apply_frequency_options(component: Any, config: dict[str, Any]) -> dict[str, Any]:
    settings = config["frequency_dependent_phase"]
    if not settings.get("enabled", True):
        return {}
    params = values(component)
    mutation: dict[str, Any] = {}
    requested = {
        "Output": "YES" if settings.get("detailed_output", True) else "NO",
        "FS": float(settings["start_hz"]),
        "FE": float(settings["end_hz"]),
        "Numf": int(settings["samples"]),
        "FDIS": settings.get("distribution", "LOG_LINEAR"),
        "Interp1": "YES" if settings.get("interpolate_travel_time", True) else "NO",
        "DCenab": "ENABLED" if settings.get("dc_correction", True) else "DISABLED",
        "DCCOR": settings.get("dc_correction_method", "FUNCTIONAL_FORM"),
    }
    passivity = settings.get("passivity", {})
    requested.update(
        {
            "CPASS": passivity.get("mode", "DISABLE"),
            "NFP": int(passivity.get("samples", 1000)),
            "FSP": float(passivity.get("start_hz", settings["start_hz"])),
            "FEP": float(passivity.get("end_hz", settings["end_hz"])),
        }
    )
    for name, value in requested.items():
        if name in params:
            mutation[name] = value
    if mutation:
        component.parameters(**mutation)
    return mutation


def conductor_layers(layer_layout: Any) -> int:
    """Return the number of conducting layers, including the core."""
    if isinstance(layer_layout, (int, float)):
        return int(layer_layout) // 2 + 1
    matches = re.findall(r"(?:^|_)C\d+", str(layer_layout), flags=re.IGNORECASE)
    return max(1, len(matches))


def eliminated_by_lc(lc_value: Any, concentric_layers: int, params: dict[str, Any], suffix: str = "") -> int:
    text = str(lc_value).casefold()
    if "all" in text and "concentric" in text:
        return concentric_layers
    if "outermost" in text:
        return min(1, concentric_layers)
    if "spec" in text:
        count = 0
        for layer in range(1, concentric_layers + 1):
            value = params.get(f"elim{layer}{suffix}")
            if value is True or "elim" in str(value).casefold():
                count += 1
        return count
    return 0


def apply_reduction(line: Any, info: dict[str, Any], state: str) -> list[dict[str, Any]]:
    if state == "not_applicable":
        return []
    if state not in {"retained", "reduced"}:
        raise ValueError(f"Unknown reduction state: {state}")
    changed: list[dict[str, Any]] = []
    dimension_delta = 0
    for component in info["all"]:
        params = values(component)
        mutation: dict[str, Any] = {}
        if "ElimGW" in params and int(params.get("NG", 0)) > 0:
            options = choice_values(component, "ElimGW")
            old_eliminated = int(params["NG"]) if "enable" in str(params["ElimGW"]).casefold() else 0
            new_eliminated = 0 if state == "retained" else int(params["NG"])
            dimension_delta += old_eliminated - new_eliminated
            mutation["ElimGW"] = choose(options, "DISABLED" if state == "retained" else "ENABLED")
        if "LC" in params:
            options = choice_values(component, "LC")
            concentric = conductor_layers(params.get("LL")) - 1
            old_eliminated = eliminated_by_lc(params["LC"], concentric, params)
            new_eliminated = 0 if state == "retained" else concentric
            dimension_delta += old_eliminated - new_eliminated
            if state == "retained":
                mutation["LC"] = choose(options, "NONE", "RETAIN")
            else:
                mutation["LC"] = choose(options, "ALL_CONCENTRIC", "OUTERMOST", "ELIM")
        elif "elimp" in params:
            if bool(params.get("poins", True)):
                old_eliminated = 1 if bool(params["elimp"]) else 0
                new_eliminated = 0 if state == "retained" else 1
                dimension_delta += old_eliminated - new_eliminated
                mutation["elimp"] = state == "reduced"
            cable_count = int(params.get("cnum", 0))
            for cable_index in range(1, cable_count + 1):
                lc_name = f"LC{cable_index}"
                ll_name = f"LL{cable_index}"
                if lc_name not in params:
                    continue
                concentric = conductor_layers(params.get(ll_name)) - 1
                old_eliminated = eliminated_by_lc(
                    params[lc_name], concentric, params, str(cable_index)
                )
                new_eliminated = 0 if state == "retained" else concentric
                dimension_delta += old_eliminated - new_eliminated
                options = choice_values(component, lc_name)
                mutation[lc_name] = (
                    choose(options, "NONE", "RETAIN")
                    if state == "retained"
                    else choose(options, "ALL_CONCENTRIC", "OUTERMOST", "ELIM")
                )
        elif "ElimGW" not in params:
            elimination_names = [name for name in params if name.casefold().startswith("elim")]
            for name in elimination_names:
                options = choice_values(component, name)
                mutation[name] = choose(options, "RETAIN", "DISABLE") if state == "retained" else choose(options, "ELIM", "ENABLE")
        if mutation:
            component.parameters(**mutation)
            changed.append({"component": str(component.defn_name), "iid": int(component.iid), "parameters": jsonable(mutation)})
    if dimension_delta:
        outer = values(line)
        old_dimension = int(outer.get("Dim", 0))
        # Cable uses Dim=0 as an automatic conductor-count sentinel. Overhead
        # TLine objects have a fixed positive port count which must be updated.
        if old_dimension > 0:
            new_dimension = old_dimension + dimension_delta
            if new_dimension < 1:
                raise ValueError(f"Reduction would produce invalid line dimension {new_dimension}")
            line.parameters(Dim=new_dimension)
            changed.append(
                {
                    "component": str(line.defn_name),
                    "iid": int(line.iid),
                    "parameters": {"Dim": new_dimension},
                    "reason": "match PSCAD interface nodes to retained/equivalent conductor count",
                }
            )
    return changed


def apply_variant(project: Any, line: Any, variant: Variant, config: dict[str, Any]) -> dict[str, Any]:
    info = relevant_components(project, line)
    mutations: dict[str, Any] = {"frequency": {}, "ground": {}, "reduction": []}
    if config["frequency_dependent_phase"].get("enabled", True):
        if info["frequency"] is None:
            raise RuntimeError(f"{line.defn_name} is not configured as a frequency-dependent phase model")
        mutations["frequency"] = apply_frequency_options(info["frequency"], config)

    if info["ground"] is not None:
        ground_mutation = dict(variant.formulations)
        if variant.resistivity_ohm_m is not None:
            ground_params = values(info["ground"])
            if "GrRho" in ground_params:
                ground_mutation["GrRho"] = choose(
                    choice_values(info["ground"], "GrRho"), "CONSTANT_RESISTIVITY"
                )
            ground_mutation["GRRES"] = variant.resistivity_ohm_m
        if ground_mutation:
            info["ground"].parameters(**ground_mutation)
        mutations["ground"] = jsonable(ground_mutation)
    mutations["reduction"] = apply_reduction(line, info, variant.reduction)
    return mutations


def changed_files(before: dict[str, tuple[int, int]], after: dict[str, tuple[int, int]]) -> list[str]:
    return sorted(
        [name for name, stamp in after.items() if before.get(name) != stamp],
        key=str.casefold,
    )


def copy_artifacts(temp: Path, names: list[str], destination: Path) -> list[dict[str, Any]]:
    rows = []
    for name in names:
        source = temp / name
        if not source.is_file():
            continue
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        rows.append(
            {
                "path": str(Path("artifacts") / name),
                "size": target.stat().st_size,
                "sha256": sha256_file(target),
                "kind": classify_artifact(target),
            }
        )
    return rows


def copy_flat_artifacts(root: Path, names: list[str], destination: Path) -> list[dict[str, Any]]:
    rows = []
    for name in names:
        source = root / name
        if not source.is_file():
            continue
        target = destination / "process_cwd" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        rows.append(
            {
                "path": str(Path("artifacts") / "process_cwd" / name),
                "size": target.stat().st_size,
                "sha256": sha256_file(target),
                "kind": classify_artifact(target),
            }
        )
    return rows


def classify_artifact(path: Path) -> str:
    suffix = path.suffix.casefold()
    if suffix == ".out" and path.parent.name.casefold() == "process_cwd":
        return "detailed_frequency_output"
    return {
        ".tli": "normalized_lcp_input",
        ".cli": "normalized_lcp_input",
        ".tlo": "emtdc_constants",
        ".clo": "emtdc_constants",
        ".out": "human_readable_matrices",
        ".log": "lcp_diagnostics",
    }.get(suffix, "detailed_or_auxiliary_output")


def solver_rejection(record_dir: Path, artifacts: list[dict[str, Any]]) -> dict[str, Any] | None:
    markers = (
        "tline.exe: error",
        "emttl ending!",
        "forrtl: severe",
        "singularity occurs",
    )
    for artifact in artifacts:
        if artifact["kind"] != "lcp_diagnostics":
            continue
        path = record_dir / artifact["path"]
        try:
            text = path.read_text(encoding="cp1252", errors="replace")
        except OSError:
            continue
        matched = [marker for marker in markers if marker in text.casefold()]
        if matched:
            tail = "\n".join(text.splitlines()[-30:])
            return {"markers": matched, "log": artifact["path"], "tail": tail}
    return None


def acquire_certificate(app: Any, allowed: bool) -> str:
    if app.licensed():
        return str(app.get_current_certificate())
    if not allowed:
        raise RuntimeError("PSCAD has no active certificate and acquisition is disabled")
    certificates = [item for item in app.get_available_certificates(refresh=True).values() if item.available]
    if not certificates:
        raise RuntimeError("No PSCAD certificate is available")
    result = app.get_certificate(certificates[0])
    if result != 0 or not app.licensed():
        raise RuntimeError(f"PSCAD certificate acquisition failed with result {result}")
    return str(app.get_current_certificate())


class PscadSession:
    def __init__(self, config: dict[str, Any]):
        self.config = config
        self.app = None
        self.certificate = None

    def __enter__(self) -> "PscadSession":
        settings = self.config["pscad"]
        self.app = mhi.pscad.launch(
            version=str(settings["version"]),
            x64=bool(settings.get("x64", True)),
            minimize=bool(settings.get("minimize", True)),
            splash=False,
            timeout=int(settings.get("launch_timeout_seconds", 60)),
        )
        self.certificate = acquire_certificate(self.app, bool(settings.get("acquire_certificate", True)))
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self.app is not None:
            # PSCAD certificate retention is user-managed; do not release it.
            try:
                self.app.quit()
            except (ConnectionError, OSError):
                # PSCAD occasionally closes the Automation socket itself after
                # a large final LCP solve. Record manifests are already durable.
                pass

    def load(self, source: Path) -> Any:
        assert self.app is not None
        self.app.load(str(source))
        return self.app.project(project_namespace(source))


def inspect_source(source: Path, config: dict[str, Any]) -> dict[str, Any]:
    source = source.resolve()
    with PscadSession(config) as session:
        project = session.load(source)
        rows = []
        for key in unique_lines(project):
            line = find_line(project, key)
            info = relevant_components(project, line)
            variants = variants_for(info, config)
            rows.append(
                {
                    "key": asdict(key),
                    "geometry_class": {
                        "aerial": info["aerial"],
                        "underground": info["underground"],
                        "mixed": info["aerial"] and info["underground"],
                    },
                    "variant_count": len(variants),
                    "variants": [variant.as_dict() for variant in variants],
                    "model": line_snapshot(project, line, include_ranges=True),
                }
            )
        result = {
            "created_at": utc_now(),
            "source": str(source),
            "source_sha256": sha256_file(source),
            "pscad_version": str(config["pscad"]["version"]),
            "certificate": session.certificate,
            "definitions": rows,
            "total_variants": sum(row["variant_count"] for row in rows),
        }
        project.unload()
        return result


def run_source(
    source: Path,
    config: dict[str, Any],
    *,
    output_root: Path | None = None,
    line_selector: str | None = None,
    dry_run: bool = False,
    limit: int | None = None,
    resume: bool = True,
    session: PscadSession | None = None,
) -> dict[str, Any]:
    source = source.resolve()
    output_root = (output_root or Path(config["output"]["root"])) .resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    owns_session = session is None
    session = session or PscadSession(config)
    if owns_session:
        session.__enter__()
    assert session.app is not None
    summary: dict[str, Any] = {
        "created_at": utc_now(),
        "source": str(source),
        "source_sha256": sha256_file(source),
        "pscad_version": str(config["pscad"]["version"]),
        "certificate": session.certificate,
        "records": [],
        "excluded_definitions": [],
    }
    source_group = output_root / f"{slug(source.stem, 38)}_{summary['source_sha256'][:8]}"
    try:
        donor_project = session.load(source)
        line_keys = unique_lines(donor_project)
        if line_selector:
            needle = line_selector.casefold()
            line_keys = [
                key
                for key in line_keys
                if needle in key.definition.casefold() or needle in key.name.casefold() or needle == str(key.iid)
            ]
            if not line_keys:
                raise LookupError(f"No line/cable definition matched {line_selector!r}")
        planned: list[tuple[LineKey, Variant]] = []
        donor_models: dict[str, Any] = {}
        for key in line_keys:
            line = find_line(donor_project, key)
            info = relevant_components(donor_project, line)
            donor_models[key.definition] = line_snapshot(donor_project, line, include_ranges=True)
            line_variants = variants_for(info, config)
            if not line_variants:
                summary["excluded_definitions"].append(
                    {
                        "line": asdict(key),
                        "reason": "not a frequency-dependent phase model in the donor PSCAD datamodel",
                    }
                )
            planned.extend((key, variant) for variant in line_variants)
        donor_project.unload()

        if limit is not None:
            planned = planned[:limit]
        summary["planned_variants"] = len(planned)
        if dry_run:
            summary["plan"] = [
                {"line": asdict(key), "variant": variant.as_dict()} for key, variant in planned
            ]
            return summary

        for index, (key, variant) in enumerate(planned, start=1):
            specification = {
                "source_sha256": summary["source_sha256"],
                "line_definition": key.definition,
                "variant": variant.as_dict(),
                "frequency_dependent_phase": config["frequency_dependent_phase"],
                "pscad_version": summary["pscad_version"],
            }
            spec_hash = sha256_json(specification)
            record_dir = source_group / slug(key.definition.split(":")[-1], 30) / f"v_{spec_hash[:16]}"
            manifest_path = record_dir / "manifest.json"
            if resume and manifest_path.exists():
                try:
                    import json

                    existing = json.loads(manifest_path.read_text(encoding="utf-8"))
                    if existing.get("specification_sha256") == spec_hash and existing.get("status") in {"success", "solver_rejected"}:
                        resumed_status = f"resumed_{existing['status']}"
                        summary["records"].append({"path": str(record_dir), "status": resumed_status})
                        print(f"[{index}/{len(planned)}] resume {key.definition} {spec_hash[:8]}", flush=True)
                        continue
                except (OSError, ValueError):
                    pass

            print(f"[{index}/{len(planned)}] run {key.definition} {variant.as_dict()}", flush=True)
            record_dir.mkdir(parents=True, exist_ok=True)
            work_dir = record_dir / "project"
            variant_source = copy_project_tree(source, work_dir)
            manifest: dict[str, Any] = {
                "schema_version": 1,
                "created_at": utc_now(),
                "status": "running",
                "source": str(source),
                "source_sha256": summary["source_sha256"],
                "source_label": "PSCAD reference",
                "line": asdict(key),
                "variant": variant.as_dict(),
                "specification_sha256": spec_hash,
                "pscad_version": summary["pscad_version"],
                "certificate": session.certificate,
                "donor_model": donor_models[key.definition],
                "mutations": {},
                "artifacts": [],
            }
            write_json(manifest_path, manifest)
            project = None
            started = time.monotonic()
            try:
                try:
                    project = session.load(variant_source)
                except (ConnectionError, OSError):
                    if not owns_session:
                        raise
                    # A completed high-port LCP fit can close PSCAD's
                    # Automation socket only after the result was harvested.
                    # Reconnect here and retry this same variant instead of
                    # consuming a whole campaign pass to discover the stale
                    # socket on the next non-resumed record.
                    session.__exit__(None, None, None)
                    session = PscadSession(config)
                    session.__enter__()
                    summary.setdefault("automation_restarts", 0)
                    summary["automation_restarts"] += 1
                    project = session.load(variant_source)
                line = find_line(project, key)
                manifest["mutations"] = apply_variant(project, line, variant, config)
                manifest["declared_model"] = line_snapshot(project, line, include_ranges=True)
                project.save()

                temp = Path(project.temp_folder)
                before = snapshot(temp)
                process_cwd = Path.cwd()
                cwd_before = flat_snapshot(process_cwd)
                line.compile()
                try:
                    immediate_messages = project.messages()
                except (ConnectionError, OSError):
                    immediate_messages = []
                has_solver_error = any(
                    str(getattr(message, "status", "")).casefold() == "error"
                    for message in immediate_messages
                )
                if has_solver_error:
                    time.sleep(0.5)
                    after = snapshot(temp)
                else:
                    after = wait_for_artifacts(
                        temp,
                        before,
                        float(config["pscad"].get("artifact_timeout_seconds", 300)),
                        float(config["pscad"].get("artifact_quiet_seconds", 2.0)),
                        process_cwd=process_cwd,
                        cwd_before=cwd_before,
                    )
                try:
                    messages = project.messages()
                    manifest["pscad_messages"] = [
                        jsonable(message._asdict() if hasattr(message, "_asdict") else message)
                        for message in messages
                    ]
                    manifest["pscad_output"] = project.output()
                except (ConnectionError, OSError) as diagnostic_error:
                    manifest["automation_warning"] = (
                        "PSCAD closed its diagnostic socket after the LCP solve: "
                        f"{diagnostic_error}"
                    )
                names = changed_files(before, after)
                cwd_after = flat_snapshot(process_cwd)
                cwd_names = changed_files(cwd_before, cwd_after)
                artifacts_dir = record_dir / "artifacts"
                manifest["artifacts"] = copy_artifacts(temp, names, artifacts_dir)
                manifest["artifacts"].extend(
                    copy_flat_artifacts(process_cwd, cwd_names, artifacts_dir)
                )
                lcp_inputs = [
                    record_dir / row["path"]
                    for row in manifest["artifacts"]
                    if row["kind"] == "normalized_lcp_input"
                ]
                constants = [row for row in manifest["artifacts"] if row["kind"] == "emtdc_constants"]
                reports = [row for row in manifest["artifacts"] if row["kind"] == "human_readable_matrices"]
                detailed = [row for row in manifest["artifacts"] if row["kind"] == "detailed_frequency_output"]
                if lcp_inputs:
                    manifest["case_id"] = case_id(summary["pscad_version"], lcp_inputs)
                rejection = solver_rejection(record_dir, manifest["artifacts"])
                if lcp_inputs and rejection:
                    manifest["status"] = "solver_rejected"
                    manifest["solver_diagnostic"] = rejection
                elif lcp_inputs and constants and (reports or detailed):
                    manifest["status"] = "success"
                    if not reports:
                        manifest.setdefault("warnings", []).append(
                            "PSCAD omitted the ordinary .out summary; detailed frequency outputs and EMTDC constants are complete."
                        )
                else:
                    raise RuntimeError(
                        "Incomplete LCP artifact set: "
                        f"inputs={len(lcp_inputs)}, constants={len(constants)}, reports={len(reports)}"
                    )
            except Exception as exc:
                manifest["status"] = "failed"
                manifest["error"] = {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}
                if owns_session and isinstance(exc, (ConnectionError, OSError)):
                    try:
                        session.__exit__(None, None, None)
                    except Exception:
                        pass
                    session = PscadSession(config)
                    session.__enter__()
                    summary.setdefault("automation_restarts", 0)
                    summary["automation_restarts"] += 1
            finally:
                manifest["elapsed_seconds"] = time.monotonic() - started
                write_json(manifest_path, manifest)
                if project is not None:
                    try:
                        project.unload()
                    except Exception:
                        pass
            summary["records"].append(
                {
                    "path": str(record_dir),
                    "status": manifest["status"],
                    "case_id": manifest.get("case_id"),
                    "error": manifest.get("error", {}).get("message"),
                }
            )
        return summary
    finally:
        summary["completed_at"] = utc_now()
        summary["successful"] = sum(row["status"] in {"success", "resumed_success"} for row in summary["records"])
        summary["solver_rejected"] = sum(
            row["status"] in {"solver_rejected", "resumed_solver_rejected"}
            for row in summary["records"]
        )
        summary["provisioned"] = summary["successful"] + summary["solver_rejected"]
        summary["failed"] = sum(row["status"] == "failed" for row in summary["records"])
        write_json(source_group / "run_summary.json", summary)
        if owns_session:
            session.__exit__(None, None, None)
