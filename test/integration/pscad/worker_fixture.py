"""Current automation protocol double; independent channel bytes, no solver oracle."""
from pathlib import Path
from types import ModuleType, SimpleNamespace
from unittest.mock import patch
import importlib.metadata
import runpy
import sys

CHANNELS = ("zm", "zp", "ym", "yp")


def output_text(channel):
    index = CHANNELS.index(channel) + 1
    return "LOG10(FN) FN element\n" + "".join(
        f"{power} {10.0**power:g} {100*index+7*sample}\n"
        for sample, power in enumerate((-1, 1), 1))


class AutomationSession:
    """One isolated session with configurable boundary failures and call counts."""
    def __init__(self, raw):
        self.raw, self.version, self.licensed = raw, "5.1.0", True
        self.automation_version = "3.1.2"
        self.compiled = self.saved = self.unloaded = self.quit = 0
        self.identity_calls, self.change_identity_after = 0, 100
        self.loaded, self.line_count, self.reject_field = [], 1, ""
        self.compile_error = self.canvas_fallback = self.diagnostics_error = self.cleanup_error = False
        self.malformed_channel = ""
        self.line = ParameterBlock(self, "current:Cable", Name="unset", Freq=37.0)
        self.frequency = ParameterBlock(self, "master:line_frephase_options", FS=3., FE=5., Numf=2, Output="NO")
        self.ground = ParameterBlock(self, "master:line_ground", EarthForm="WEDEPOHL")
        self.cable = ParameterBlock(self, "master:Cable_Coax", elim1="RETAIN", elim2=0)
        self.components = [self.frequency, self.ground, self.cable]
        self.canvas = SimpleNamespace(components=lambda: list(self.components))
        self.project = ProjectSession(self)
        self.app = ApplicationSession(self)

    def identify(self, version, app):
        self.identity_calls += 1
        return {"version": version, "installation":
                "changed" if self.identity_calls > self.change_identity_after else "fixture"}


class ParameterBlock:
    def __init__(self, session, definition, **parameters):
        self.session, self.defn_name, self.values = session, definition, parameters

    def parameters(self, **updates):
        translations = {("Output", 1): "YES", ("EarthForm", 2): "DIRECT_NUMERICAL_INTEGRATION"}
        self.values.update({key: translations.get((key, value), value)
                            for key, value in updates.items() if key != self.session.reject_field})
        return dict(self.values)

    def canvas(self):
        return self.session.canvas

    def compile(self):
        self.session.compiled += 1
        if self.session.compile_error:
            raise RuntimeError("synthetic compile failure")
        destination = Path(self.session.raw)
        destination.mkdir()
        for channel in CHANNELS:
            text = "malformed channel\n" if channel == self.session.malformed_channel else output_text(channel)
            (destination / f"current_{channel}.out").write_text(text)


class ProjectSession:
    def __init__(self, session):
        self.session, self.temp_folder = session, session.raw

    def find_all(self, kind):
        return [self.session.line for _ in range(self.session.line_count)] if kind == "Cable" else []

    def canvas(self, name):
        if self.session.canvas_fallback:
            raise RuntimeError("component owns canvas")
        return self.session.canvas

    def messages(self):
        if self.session.diagnostics_error:
            raise RuntimeError("message retrieval failed")
        return ["current transport diagnostic"]

    def output(self):
        if self.session.diagnostics_error:
            raise RuntimeError("output retrieval failed")
        return "current automation output"

    def save(self):
        self.session.saved += 1

    def unload(self):
        self.session.unloaded += 1
        if self.session.cleanup_error:
            raise RuntimeError("synthetic unload failure")


class ApplicationSession:
    def __init__(self, session):
        self.session = session

    @property
    def version(self):
        return self.session.version

    def licensed(self):
        return self.session.licensed

    def load(self, path):
        self.session.loaded.append(path)

    def project(self, name):
        return self.session.project

    def quit(self):
        self.session.quit += 1
        if self.session.cleanup_error:
            raise RuntimeError("synthetic quit failure")


def install(raw):
    state = AutomationSession(raw)
    parent, automation = ModuleType("mhi"), ModuleType("mhi.pscad")
    parent.pscad = automation
    automation.launch = lambda **kwargs: state.app
    state.patches = [patch.dict(sys.modules, {"mhi": parent, "mhi.pscad": automation}),
                     patch.object(importlib.metadata, "version", lambda _: state.automation_version),
                     patch.object(runpy, "run_path", lambda _: {"identify": state.identify})]
    for replacement in state.patches:
        replacement.start()
    return state


def uninstall(state):
    for replacement in reversed(state.patches):
        replacement.stop()
