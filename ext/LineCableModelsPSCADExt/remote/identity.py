"""Identify PSCAD's selected line-constants implementation without solving a model."""

import hashlib
import platform
from importlib import metadata
from pathlib import Path

from mhi.common import process
import mhi.pscad


def identify(version, app=None):
    if version != "5.1.0":
        raise ValueError("The LCM PSCAD adapter supports version 5.1.0")
    executable = Path(process.find_exe("PSCAD", version=version, x64=True)).resolve()
    root = executable.parents[2]
    automation_version = metadata.version("mhi.pscad")
    if automation_version != "3.1.2":
        raise ValueError("The LCM PSCAD adapter requires mhi.pscad 3.1.2")
    owned = app is None
    try:
        if owned:
            app = mhi.pscad.launch(version=version, x64=True, minimize=True,
                                   splash=False, silence=True, timeout=60)
        if str(app.version) != version:
            raise ValueError("PSCAD launched an unexpected application version")
        # Version-bound workaround: mhi.pscad 3.1.2 settings() initializes its
        # Fortran codec even for a read, raising "Unable to retrieve detected
        # software" on this 5.1 station. This is its own underlying read-only
        # settings call; bypass only that unrelated decoder, not solver settings.
        # Revisit when the required automation version fixes that reader.
        selected = app._settings({})["file_lcp"]
        expanded = str(selected).replace("$(HomeDir)", str(root))
        solver = Path(expanded)
        if "$(" in expanded or not solver.is_absolute():
            raise ValueError("Cannot identify PSCAD's selected LCP executable: " + str(selected))
        files = {
            "pscad": executable,
            "line_constants": solver.resolve(),
            "master_library": root / "master.pslx",
        }
        result = {
            "schema": "1",
            "version": version,
            "python_version": platform.python_version(),
            "automation_version": automation_version,
            "common_version": metadata.version("mhi.common"),
            "profile": "saved profile; numerical inputs supplied explicitly",
        }
        for name, path in files.items():
            result[name + "_path"] = str(path)
            with path.open("rb") as stream:
                result[name + "_sha256"] = hashlib.file_digest(stream, "sha256").hexdigest()
        return result
    finally:
        if owned and app is not None:
            app.quit()
