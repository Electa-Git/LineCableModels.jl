"""Local automation double. Matrices test file transport, never solver accuracy."""
from pathlib import Path
from types import ModuleType, SimpleNamespace
from unittest.mock import patch
import importlib.metadata
import runpy
import sys


class Component:
    def __init__(self, name, parameters, state):
        self.defn_name = name
        self.values = parameters.copy()
        self.state = state

    def parameters(self, **updates):
        for name, value in updates.items():
            if name == self.state.reject_field:
                continue
            if name == "Output" and value == 1:
                value = "YES"
            if name == "EarthForm" and value == 2:
                value = "DIRECT_NUMERICAL_INTEGRATION"
            self.values[name] = value
        return self.values.copy()

    def canvas(self):
        return self.state.canvas

    def compile(self):
        self.state.compiled += 1
        if self.state.compile_error:
            raise RuntimeError("synthetic compile failure")
        output = Path(self.state.raw)
        output.mkdir(exist_ok=True)
        for suffix in ("zm", "zp", "ym", "yp"):
            (output / ("synthetic_" + suffix + ".out")).write_text(
                "LOG10(FN) FN element\n-1 0.1 2\n1 10 3\n"
            )


class Project:
    def __init__(self, state):
        self.state = state
        self.temp_folder = state.raw

    def find_all(self, kind):
        return [self.state.line] * self.state.line_count if kind == "Cable" else []

    def canvas(self, name):
        if self.state.canvas_fallback:
            raise RuntimeError("use component canvas")
        return self.state.canvas

    def messages(self):
        if self.state.diagnostics_error:
            raise RuntimeError("synthetic message failure")
        return ["automation fixture diagnostic"]

    def output(self):
        if self.state.diagnostics_error:
            raise RuntimeError("synthetic output failure")
        return "automation fixture output"

    def save(self):
        self.state.saved += 1

    def unload(self):
        self.state.unloaded += 1
        if self.state.cleanup_error:
            raise RuntimeError("synthetic unload failure")


class Application:
    def __init__(self, state):
        self.state = state

    @property
    def version(self):
        return self.state.version

    def licensed(self):
        return self.state.licensed

    def load(self, path):
        self.state.loaded.append(path)

    def project(self, name):
        return self.state.project

    def quit(self):
        self.state.quit += 1
        if self.state.cleanup_error:
            raise RuntimeError("synthetic quit failure")


def install(raw):
    state = SimpleNamespace(
        raw=raw, version="5.1.0", licensed=True, compiled=0, saved=0,
        unloaded=0, quit=0, loaded=[], line_count=1, reject_field="",
        compile_error=False, canvas_fallback=False, diagnostics_error=False,
        cleanup_error=False, identity_calls=0, change_identity_after=100,
        automation_version="3.1.2",
    )
    state.line = Component("fixture:Cable", {"Name": "old", "Freq": 50.0}, state)
    state.frequency = Component("master:line_frephase_options",
                                {"FS": 1., "FE": 2., "Numf": 100, "Output": "NO"}, state)
    state.ground = Component("master:line_ground", {"EarthForm": "WEDEPOHL"}, state)
    state.cable = Component("master:Cable_Coax", {"elim1": "RETAIN", "elim2": 0}, state)
    state.components = [state.frequency, state.ground, state.cable]
    state.canvas = SimpleNamespace(components=lambda: state.components)
    state.project = Project(state)
    state.app = Application(state)
    module = ModuleType("mhi.pscad")
    module.launch = lambda **kwargs: state.app
    parent = ModuleType("mhi")
    parent.pscad = module

    def identify(version, app):
        state.identity_calls += 1
        return {"version": version, "installation":
                "changed" if state.identity_calls > state.change_identity_after else "fixture"}

    state.patches = [
        patch.dict(sys.modules, {"mhi": parent, "mhi.pscad": module}),
        patch.object(importlib.metadata, "version", lambda name: state.automation_version),
        patch.object(runpy, "run_path", lambda path: {"identify": identify}),
    ]
    for item in state.patches:
        item.start()
    return state


def uninstall(state):
    for item in reversed(state.patches):
        item.stop()
