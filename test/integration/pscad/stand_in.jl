# The PSCAD backend runs against this stand-in for `mhi.pscad` without PSCAD. The stand-in
# application loads the staged project from the shared directory and applies the settings as
# PSCAD reads them back. A compile call of the line writes the detailed files of the solve
# into the working directory, the test's `local_root`. The matrices are protocol fixtures,
# not references. A stand-in for the SSH command forwards the port and starts PSCAD.
@testmodule PSCADStandIn begin
    using LineCableModels
    using PythonCall

    const SOURCE = raw"""
import math
import os
import sys
import types
import xml.etree.ElementTree as ElementTree
from pathlib import PureWindowsPath

# PSCAD's readback of the choice fields that a native run sets.
CHOICES = {
    "EarthForm": {0: "WEDEPOHL", 2: "DIRECT_NUMERICAL_INTEGRATION", 3: "SAAD"},
    "EarthForm2": {0: "DERISEMLYEN", 2: "DIRECT_NUMERICAL_INTEGRATION"},
    "EarthForm3": {0: "AMETANIL", 2: "LUCCA"},
    "enablf": {1: "YES"},
    "Output": {1: "YES"},
}
# The fields that PSCAD's definitions add to the fields of a staged project.
DEFAULTS = {"master:line_frephase_options": {"Output": "NO"}}
# The entry (row, column) at the k-th frequency of each detailed file is scale*(100r+10c+k).
SCALES = {"zm": 1.0, "zp": 0.0, "ym": 1e-6, "yp": 0.0}


def fields(element):
    return {p.get("name"): p.get("value") for p in element.find("paramlist").iter("param")}


class Component:
    def __init__(self, application, definition, parameters):
        self._application = application
        self.defn_name = definition
        self._parameters = {**DEFAULTS.get(definition.lower(), {}), **parameters}

    def parameters(self, **updates):
        for name, value in updates.items():
            value = CHOICES.get(name, {}).get(value, value)
            self._parameters[name] = self._application.readbacks.get(name, value)
        return dict(self._parameters)


class Line(Component):
    def __init__(self, application, definition, parameters, canvas):
        super().__init__(application, definition, parameters)
        self._canvas = canvas

    def compile(self):
        application = self._application
        application.compiles += 1
        application.messages.append("ACLine segment manual solve: started.")
        if application.compile_error:
            raise RuntimeError(application.compile_error)
        if not application.solves:
            return
        components = self._canvas.components()
        frequency = next(c for c in components
                         if c.defn_name.lower() == "master:line_frephase_options")._parameters
        start, end, steps = float(frequency["FS"]), float(frequency["FE"]), int(frequency["Numf"])
        size = sum(1 for c in components if c.defn_name.lower() == "master:cable_coax"
                   for k in range(1, 5) if c._parameters.get(f"CONNAM{k}", "none") != "none")
        name = self._parameters["Name"]
        span = math.log10(end) - math.log10(start)
        for kind, scale in SCALES.items():
            path = os.path.join(application.working_directory, f"{name}_{kind}.out")
            with open(path, "w") as stream:
                stream.write(" LOG10(FN) FN matrix\n")
                for k in range(1, steps + 2):
                    f = 10 ** (math.log10(start) + (k - 1) * span / steps)
                    entries = [scale * (100 * r + 10 * c + k)
                               for r in range(1, size + 1) for c in range(1, size + 1)]
                    entries += [0.0] if application.malformed else []
                    stream.write(" ".join(repr(x) for x in [math.log10(f), f] + entries) + "\n")
        with open(os.path.join(application.working_directory, f"{name}_zs.out"), "w") as stream:
            stream.write("sequence impedance\n")
        with open(os.path.join(application.working_directory, f"{name}.log"), "w") as stream:
            stream.write("line constants log\n")


class Canvas:
    def __init__(self):
        self._components = []

    def components(self):
        return list(self._components)


class Project:
    def __init__(self, application, path):
        root = ElementTree.parse(path).getroot()
        self._application = application
        self.name = os.path.splitext(os.path.basename(path))[0]
        self.namespace = root.get("name")
        self.saves = 0
        row = next(d for d in root.iter("Definition") if d.get("name") == "CableSystem")
        self._canvas = Canvas()
        self._canvas._components = [Component(application, u.get("defn"), fields(u))
                                    for u in row.iter("User")]
        wire = next(w for w in root.iter("Wire") if w.get("classid") == "Cable")
        self.line = Line(application, wire.get("defn"), fields(wire.find("User")), self._canvas)

    def find_all(self, kind):
        return [self.line] if kind == "Cable" else []

    def canvas(self, name):
        return self._canvas

    def messages(self):
        return list(self._application.messages)

    def output(self):
        if self._application.output_error:
            raise RuntimeError(self._application.output_error)
        return "stand-in project output"

    def save(self):
        self.saves += 1

    def unload(self):
        del self._application.loaded[self.name]


class Application:
    version = "5.1.0"

    def __init__(self, working_directory, shared_root, station):
        self.working_directory = working_directory
        self.shared_root = shared_root
        self.station = station
        self.silence = False
        self.connects = 0
        self.compiles = 0
        self.quits = 0
        self.loads = []
        self.loaded = {}
        self.messages = ["PSCAD stand-in message"]
        self.file_lcp = r"$(HomeDir)\bin\win64\tline.exe"
        self.compile_error = ""
        self.output_error = ""
        self.solves = True
        self.malformed = False
        self.readbacks = {}

    def licensed(self):
        return True

    def settings(self):
        raise RuntimeError("Unable to retrieve detected software")

    def _settings(self, settings):
        return {"file_lcp": self.file_lcp}

    # PSCAD exits, and with it the SSH command that started it.
    def quit(self):
        self.quits += 1
        os.remove(os.path.join(self.station, "pscad-listens"))

    def projects(self):
        return [{"name": "master", "type": "Library", "description": "Master Library"}] + \
            [{"name": name, "type": "Case", "description": ""} for name in self.loaded]

    def project(self, name):
        return self.loaded[name]

    def _load(self, path):
        self.loads.append(path)
        parts = PureWindowsPath(path).relative_to(PureWindowsPath(self.shared_root)).parts
        project = Project(self, os.path.join(self.working_directory, *parts))
        self.loaded[project.name] = project


def connect(host=None, port=0, timeout=5):
    pscad.application.connects += 1
    return pscad.application


pscad = types.ModuleType("mhi.pscad")
pscad.connect = connect
pscad.application = None
mhi = types.ModuleType("mhi")
mhi.__path__ = []
mhi.pscad = pscad
sys.modules["mhi"] = mhi
sys.modules["mhi.pscad"] = pscad
"""
    const PYTHON = pydict()
    pyexec(SOURCE, PYTHON)

    # The SSH command's stand-in. Its first argument is the station's directory of the test,
    # and the others are those of an SSH command. With `-N` it forwards: it listens on the
    # caller's port of `-L` and keeps every connection open. `netstat -ano -p TCP` lists the
    # PSCAD that listens on the station, `pscad-listens`, which holds its port. With a command
    # that starts PSCAD, it records the command, then PSCAD listens until it quits or the
    # command is ended. A `never-listens` file makes a start that never listens. Without
    # Julia's signal handlers, SIGTERM ends the stand-in at once. It runs one thread, since
    # another thread would meet a GC safepoint without those handlers.
    const SSH = joinpath(mktempdir(), "ssh.jl")
    write(SSH, raw"""
        # The stand-in, a program of its own, loads the standard library whatever load path
        # its caller passes on.
        push!(LOAD_PATH, "@stdlib")
        using Sockets
        station, arguments = ARGS[1], ARGS[2:end]
        listens = joinpath(station, "pscad-listens")
        if "-N" in arguments
            write(joinpath(station, "forward-$(getpid()).pid"), "")
            caller = parse(Int, split(arguments[findfirst(==("-L"), arguments) + 1], ':')[2])
            server = listen(ip"127.0.0.1", caller)
            held = TCPSocket[]
            while true
                push!(held, accept(server))
            end
        elseif last(arguments) == "netstat -ano -p TCP"
            println("\nActive Connections\n\n  Proto  Local Address          Foreign Address        State           PID")
            isfile(listens) && println("  TCP    127.0.0.1:$(read(listens, String))    0.0.0.0:0    LISTENING    4242")
        else
            write(joinpath(station, "start-$(getpid()).pid"), "")
            write(joinpath(station, "started.txt"), last(arguments))
            while isfile(joinpath(station, "never-listens"))
                sleep(1)
            end
            write(listens, match(r"/port:(\d+)", last(arguments))[1])
            while isfile(listens)
                sleep(0.1)
            end
        end
        """)
    const SHARED_ROOT = raw"Z:\exchange"
    const PORT = 64339
    const PSCAD_EXECUTABLE = raw"C:\Program Files (x86)\PSCAD\5.1.0\bin\win64\Pscad.exe"
    const STATION = Ref("")

    # A fresh stand-in PSCAD whose working directory is `local_root`. It listens on the
    # station when it is `listening`, and the run borrows it. Otherwise the run starts it.
    function install(local_root::AbstractString; listening::Bool = true)
        STATION[] = mktempdir()
        listening && write(joinpath(STATION[], "pscad-listens"), string(PORT))
        application = PYTHON["Application"](local_root, SHARED_ROOT, STATION[])
        PYTHON["pscad"].application = application
        return application
    end

    # The station of the stand-in, the named tuple form of the `remote` computation option.
    remote(local_root::AbstractString; timeout = 60) =
        (host = "station", port = PORT, shared_root = SHARED_ROOT,
            pscad_executable = PSCAD_EXECUTABLE, local_root, timeout,
            ssh = [joinpath(Sys.BINDIR, Base.julia_exename()), "--startup-file=no",
                "--handle-signals=no", "--threads=1,0", SSH, STATION[]])

    # The processes of the SSH command's stand-in that a run started, by their process ids.
    processes(kind) = [parse(Int, chop(file; head = length(kind) + 1, tail = 4))
                       for file in readdir(STATION[]) if startswith(file, kind * "-")]
    # Whether no process has the process id `pid`.
    ended(pid) = ccall(:kill, Cint, (Cint, Cint), pid, 0) != 0

    # The stand-in's counters and states, as Julia values.
    connects(application) = pyconvert(Int, application.connects)
    quits(application) = pyconvert(Int, application.quits)
    compiles(application) = pyconvert(Int, application.compiles)
    silence(application) = pyconvert(Bool, application.silence)
    loaded(application) = pyconvert(Vector{String}, pybuiltins.list(application.loaded))
end
