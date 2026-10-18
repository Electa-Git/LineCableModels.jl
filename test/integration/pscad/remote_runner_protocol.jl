using Test
using PythonCall
using TOML

const REPOSITORY_ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
const REMOTE_ROOT=joinpath(
    REPOSITORY_ROOT,"src","pscad","remote")

include(joinpath(REMOTE_ROOT,"runner.jl"))

const PYTHON_AUTOMATION_FIXTURE=raw"""
import importlib.metadata
import os
import runpy
import sys
import types
from collections import namedtuple
import xml.etree.ElementTree as ET

IDENTITY = {
    "schema": "1",
    "version": "5.1.0",
    "python_version": "fixture",
    "automation_version": "3.1.2",
    "common_version": "fixture",
    "profile": "fixture",
    "pscad_path": "fixture",
    "pscad_sha256": "0" * 64,
    "line_constants_path": "fixture",
    "line_constants_sha256": "1" * 64,
    "master_library_path": "fixture",
    "master_library_sha256": "2" * 64,
}

Message = namedtuple("Message", "text label status name")

class Component:
    def __init__(self, definition, **parameters):
        self.defn_name = definition
        self._parameters = parameters

    def parameters(self, **updates):
        self._parameters.update(updates)
        if updates.get("Output") == 1:
            self._parameters["Output"] = "YES"
        return dict(self._parameters)

class Line(Component):
    def __init__(self):
        super().__init__("master:fixture", Name="unset")

    def compile(self):
        LAST_APP._project.compiled = True
        if os.environ.get("LCM_RUNNER_FIXTURE_FAIL") == "1" or os.environ.get("LCM_RUNNER_MODE") == "build_error_throw":
            raise RuntimeError("synthetic compile failure")
        root = os.environ["LCM_RUNNER_FIXTURE_OUTPUT"]
        for suffix in ("_zm.out", "_zp.out", "_ym.out", "_yp.out"):
            with open(os.path.join(root, "raw" + suffix), "w") as stream:
                stream.write("header\n")
                for index in range(3):
                    stream.write(f"{index} {index + 1}\n")

class Canvas:
    def __init__(self, components):
        self._components = components

    def components(self):
        return self._components

class Project:
    def __init__(self):
        self.line = Line()
        self.frequency = Component(
            "master:line_frephase_options", FS=0.0, FE=0.0, Numf=0, Output="NO")
        self.ground = Component("master:line_ground")
        self.cable = Component("master:cable_coax", elim1=0, elim2="RETAIN")
        self._canvas = Canvas([self.frequency, self.ground, self.cable])
        self.temp_folder = os.environ["LCM_RUNNER_FIXTURE_OUTPUT"]
        self.unloaded = False
        self.filename = None
        self.compiled = False
        self.ran = False

    def find_all(self, kind):
        return [self.line] if kind == "TLine" else []

    def canvas(self, name):
        return self._canvas

    def messages(self):
        status = "error" if os.environ.get("LCM_RUNNER_MODE") in ("build_error", "build_error_throw") else "normal"
        text = "An alert was issued while solving line constants, check the log file 'generated.log'" if status == "error" else "fixture message"
        return [Message(text, "Build", status, "generated")]

    def output(self):
        return "fixture project output"

    def save(self):
        if os.environ.get("LCM_RUNNER_MODE") == "moved":
            document = ET.parse(self.filename)
            for node in document.findall("./definitions/Definition/schematic/User"):
                node.set("y", str(int(node.get("y")) - 18))
            document.write(self.filename)

    def unload(self):
        self.unloaded = True

class Application:
    version = "5.1.0"
    def __init__(self, silence=False):
        self._project = Project()
        self._silence = silence
        self.silence_updates = []
        self.silence_at_load = None
        self.restore_after_unload = None
        self.closed = False
        self.quit_calls = 0
        self.settings_calls = 0
        self.loaded_path = None
        self.other = Project()
        self.other.filename = os.path.join(os.environ["LCM_RUNNER_FIXTURE_OUTPUT"], "user.pscx")
        mode = os.environ.get("LCM_RUNNER_MODE")
        self.existing = None
        if mode in ("stale", "foreign"):
            directory = os.environ["LCM_RUNNER_SCRATCH"] if mode == "stale" else os.environ["LCM_RUNNER_SCRATCH"] + "-foreign"
            self.existing = Project()
            self.existing.filename = os.path.join(directory, "old-run", "generated.pscx")

    @property
    def silence(self):
        return self._silence

    @silence.setter
    def silence(self, value):
        self.silence_updates.append(value)
        if len(self.silence_updates) > 1:
            self.restore_after_unload = self._project.unloaded
            if os.environ.get("LCM_RUNNER_MODE") in ("restore_error", "run_restore_error"):
                raise RuntimeError("synthetic silence restore failure")
        self._silence = value

    def licensed(self):
        return os.environ.get("LCM_RUNNER_MODE") != "license"

    def load(self, path):
        self.silence_at_load = self.silence
        if not self.silence:
            raise RuntimeError("modal dialog would block")
        self.loaded_path = path
        self._project.filename = path
        if os.environ.get("LCM_RUNNER_MODE") == "load_error":
            raise RuntimeError("synthetic partial load failure")

    def project(self, name):
        if name == "user_project":
            return self.other
        if self.existing is not None and not self.existing.unloaded and name == "generated":
            return self.existing
        if self.loaded_path is None:
            raise RuntimeError("project not loaded")
        if name != os.path.splitext(os.path.basename(self.loaded_path))[0]:
            raise RuntimeError("unknown project " + name)
        return self._project

    def projects(self):
        entries = [{"name": "master", "type": "Library", "description": "Master Library"},
                   {"name": "user_project", "type": "Case", "description": "User case"}]
        if self.existing is not None and not self.existing.unloaded:
            entries.append({"name": "generated", "type": "Case", "description": "Existing case"})
        if self.loaded_path is not None and not self._project.unloaded:
            entries.append({"name": "generated", "type": "Case", "description": "Our case"})
        return entries

    def settings(self):
        self.settings_calls += 1
        if os.environ.get("LCM_RUNNER_MODE") == "session0":
            raise SystemError("Unable to retrieve detected software")
        return {"fortran_version": "GFortran 8.1 (64-bit)"}

    def _settings(self, updates):
        return {"fortran_version": "fixture-fortran"}

    def quit(self):
        self.quit_calls += 1
        self.closed = True

pscad = types.ModuleType("mhi.pscad")
def launch(**kwargs):
    global LAST_APP, LAST_LAUNCH
    LAST_LAUNCH = kwargs
    LAST_APP = Application(silence=True)
    return LAST_APP
pscad.launch = launch
def connect(**kwargs):
    global LAST_APP, LAST_CONNECT
    LAST_CONNECT = kwargs
    LAST_APP = Application(silence=os.environ.get("LCM_RUNNER_SILENCE") == "true")
    return LAST_APP
pscad.connect = connect
mhi = types.ModuleType("mhi")
mhi.pscad = pscad
sys.modules["mhi"] = mhi
sys.modules["mhi.pscad"] = pscad
importlib.metadata.version = lambda name: "3.1.2" if name == "mhi.pscad" else "fixture"
runpy._lcm_original_run_path = getattr(runpy, "_lcm_original_run_path", runpy.run_path)
runpy.run_path = lambda path: {"identify": lambda version, app: dict(IDENTITY)}
"""

function runner_arguments(directory,output;verbosity="0", mode="constants", port=0, owned=true)
    project=joinpath(directory,"generated.pscx")
    write(project, "<project><definitions><Definition name=\"Main\"><schematic>" *
        "<User id=\"1\" x=\"72\" y=\"3168\" defn=\"master:voltmetergnd\"/>" *
        "</schematic></Definition></definitions></project>")
    input=Dict(
        "schema_version"=>4,
        "native_settings"=>Dict(
            "ground"=>Dict{String,Any}(),
            "frequency"=>Dict{String,Any}(),
            "configuration"=>Dict{String,Any}(),
        ),
    )
    open(joinpath(directory,"computation.toml"),"w") do io
        TOML.print(io,input;sorted=true)
    end
    return [project,output,"generated","fixture","pscad","10.0","1000.0","2",
        "5.1.0",verbosity,mode,string(port),string(owned),directory]
end

@testset "Remote PSCAD runner protocol without PSCAD" begin
    pyexec(PYTHON_AUTOMATION_FIXTURE,Main)

    @test _same_value(1,1.0)
    @test !_same_value(1,2)
    @test _same_value("YES",pystr("YES"))
    for value in (0,"no","FALSE"," disabled ","retain","none")
        @test _disabled(value)
    end
    @test !_disabled("YES")
    @test_throws ArgumentError main(String[])

    mktempdir(mkpath(joinpath(REPOSITORY_ROOT, ".linecablemodels", "test", "runner-protocol"))) do directory
        output=joinpath(directory,"outputs")
        arguments=runner_arguments(directory,output)
        @test_throws ArgumentError main(setindex!(copy(arguments),"3",10))
        @test_throws ArgumentError main(setindex!(copy(arguments),"5.0.2",9))
        @test_throws ArgumentError main(setindex!(copy(arguments),"bad-stem",4))
        @test_throws ArgumentError main(setindex!(copy(arguments), "1023", 12))
        @test_throws ArgumentError main(setindex!(copy(arguments), "false", 13))
        @test arguments[11:14] == ["constants", "0", "true", directory]

        python_environment=pyimport("os").environ
        python_environment["LCM_RUNNER_FIXTURE_OUTPUT"]=output
        try
            @test main(arguments) === nothing
        finally
            python_environment.pop("LCM_RUNNER_FIXTURE_OUTPUT",nothing)
        end
        for name in ("result_zm.out","result_zp.out","result_ym.out","result_yp.out")
            @test _data_rows(joinpath(output,name)) == 3
        end
        @test TOML.parsefile(joinpath(output,"solver.toml"))["version"] == "5.1.0"
        @test TOML.parsefile(joinpath(output,"native-settings.toml")) == Dict(
            "ground"=>Dict{String,Any}(),
            "frequency"=>Dict{String,Any}(),
            "configuration"=>Dict{String,Any}(),
        )
        console=read(joinpath(output,"pscad-console.txt"),String)
        @test occursin("Starting PSCAD line-constants computation",console)
        @test occursin("Collected detailed PSCAD Z and Y outputs",console)
    end

    mktempdir(mkpath(joinpath(REPOSITORY_ROOT, ".linecablemodels", "test", "runner-protocol"))) do directory
        output=joinpath(directory,"failed")
        arguments=runner_arguments(directory,output;verbosity="2")
        python_environment=pyimport("os").environ
        python_environment["LCM_RUNNER_FIXTURE_OUTPUT"]=output
        python_environment["LCM_RUNNER_FIXTURE_FAIL"]="1"
        caught=try
            main(arguments)
            nothing
        catch error
            error
        finally
            python_environment.pop("LCM_RUNNER_FIXTURE_OUTPUT",nothing)
            python_environment.pop("LCM_RUNNER_FIXTURE_FAIL",nothing)
        end
        @test caught isa ErrorException
        @test occursin("line-constants computation",sprint(showerror,caught))
        @test occursin("synthetic compile failure",read(joinpath(output,"pscad-console.txt"),String))
    end
end

const PYTHON_TRANSIENT_FIXTURE=raw"""
import shutil

def transient_parameters(self, **updates):
    if not hasattr(self, "_settings"):
        self._settings = dict(time_duration="0.019", time_step="1000",
                              sample_step="1000", PlotType="1", output_filename="unset.out")
    self._settings.update(updates)
    return dict(self._settings)

def transient_run(self):
    self.ran = True
    if os.environ.get("LCM_RUNNER_MODE") in ("run_error", "run_restore_error", "build_error_throw"):
        raise RuntimeError("synthetic run failure")
    root = os.environ["LCM_RUNNER_FIXTURE_OUTPUT"]
    stem = self.parameters()["output_filename"][:-4]
    target = os.path.join(root, "native")
    os.makedirs(target, exist_ok=True)
    fixture = os.environ["LCM_RUNNER_FIXTURE_CHANNELS"]
    shutil.copyfile(os.path.join(fixture, "channels.inf"), os.path.join(target, stem + ".inf"))
    shutil.copyfile(os.path.join(fixture, "channels_01.out"), os.path.join(target, stem + "_01.out"))

Project.parameters = transient_parameters
Project.run = transient_run
"""

@testset "Remote PSCAD transient runner protocol without PSCAD" begin
    pyexec(PYTHON_AUTOMATION_FIXTURE,Main)
    pyexec(PYTHON_TRANSIENT_FIXTURE,Main)
    root=mkpath(joinpath(REPOSITORY_ROOT,".linecablemodels","test","transient-runner"))
    mktempdir(root) do directory
        output=joinpath(directory,"outputs")
        arguments=runner_arguments(directory,output; mode="transient")
        @test arguments[4] == "fixture"
        @test arguments[11] == "transient"
        @test arguments[12:14] == ["0", "true", directory]
        @test_throws ArgumentError main(setindex!(copy(arguments), "invalid", 11))
        environment=pyimport("os").environ
        environment["LCM_RUNNER_FIXTURE_OUTPUT"]=output
        environment["LCM_RUNNER_FIXTURE_CHANNELS"]=joinpath(
            REPOSITORY_ROOT,"test","fixtures","pscad","transient")
        environment["LCM_RUNNER_FIXTURE_FAIL"]="1"
        try
            @test main(arguments) === nothing
        finally
            for name in ("LCM_RUNNER_FIXTURE_OUTPUT","LCM_RUNNER_FIXTURE_CHANNELS",
                    "LCM_RUNNER_FIXTURE_FAIL")
                environment.pop(name,nothing)
            end
        end
        @test TOML.parsefile(joinpath(output,"execution.toml")) ==
            Dict("mode"=>"transient","stem"=>"fixture")
        settings=TOML.parsefile(joinpath(output,"native-settings.toml"))["project"]
        @test settings["output_filename"] == "fixture.out"
        @test settings["sample_step"] == "1000"
        @test parse(Float64,read(joinpath(output,"timing.txt"),String)) >= 0
        @test isfile(joinpath(output,"solver.toml"))
        @test !isfile(joinpath(output,"result_zm.out"))
        for name in ("fixture.inf","fixture_01.out")
            @test read(joinpath(output,name)) == read(joinpath(output,"native",name))
        end
        console=read(joinpath(output,"pscad-console.txt"),String)
        @test occursin("Starting PSCAD project.run()",console)
        @test occursin("Collected PSCAD transient outputs",console)
        @test occursin("Project fitting readbacks: FS=0.0 Hz, FE=0.0 Hz, Numf=0", console)
        @test !occursin("Applied frequency range 10.0 Hz to 1000.0 Hz", console)
    end
end

@testset "Remote PSCAD application acquisition and load geometry" begin
    pyexec(PYTHON_AUTOMATION_FIXTURE, Main)
    pyexec(PYTHON_TRANSIENT_FIXTURE, Main)
    root = mkpath(joinpath(REPOSITORY_ROOT, ".linecablemodels", "test", "station-runner"))
    for mode in ("constants", "transient"), (port, owned) in ((0, true), (50000, false), (50000, true)),
            scenario in ("success", "stale", "foreign", "moved", "load_error", "session0",
                "run_error", "build_error", "build_error_throw", "restore_error", "run_restore_error")
        mktempdir(root) do directory
            output = joinpath(directory, "outputs")
            arguments = runner_arguments(directory, output; mode, port, owned)
            if scenario in ("success", "build_error", "build_error_throw")
                write(joinpath(directory, "generated.log"),
                    "omitted first diagnostic\n" * repeat("line constants detail\n", 20) *
                    "synthetic line constants failure\n")
                write(joinpath(mkpath(joinpath(directory, "generated.gf81")), "cable.log"),
                    "synthetic cable log ending\n")
            end
            @test arguments[11:14] == [mode, string(port), string(owned), directory]
            environment = pyimport("os").environ
            environment["LCM_RUNNER_FIXTURE_OUTPUT"] = output
            environment["LCM_RUNNER_FIXTURE_CHANNELS"] = joinpath(REPOSITORY_ROOT, "test", "fixtures", "pscad", "transient")
            environment["LCM_RUNNER_SCRATCH"] = directory
            environment["LCM_RUNNER_MODE"] = scenario
            caught = try
                main(arguments)
                nothing
            catch error
                error
            finally
                for key in ("LCM_RUNNER_FIXTURE_OUTPUT", "LCM_RUNNER_FIXTURE_CHANNELS", "LCM_RUNNER_SCRATCH", "LCM_RUNNER_MODE")
                    environment.pop(key, nothing)
                end
            end
            failed = scenario in ("foreign", "moved", "load_error", "build_error", "build_error_throw") ||
                (mode == "transient" && scenario in ("session0", "run_error", "run_restore_error"))
            @test (caught !== nothing) == failed
            @test pyconvert(Bool, pyeval("LAST_APP.closed", Main)) == owned
            @test !pyconvert(Bool, pyeval("LAST_APP.other.unloaded", Main))
            before_load = scenario == "foreign" || (mode == "transient" && scenario == "session0")
            @test pyconvert(Bool, pyeval("LAST_APP.loaded_path is None", Main)) == before_load
            @test pyconvert(Bool, pyeval("LAST_APP._project.unloaded", Main)) == !before_load
            @test pyconvert(Int, pyeval("LAST_APP.settings_calls", Main)) == Int(mode == "transient")
            restore_fails = !owned && scenario in ("restore_error", "run_restore_error")
            @test pyconvert(Bool, pyeval("LAST_APP.silence", Main)) == (owned || restore_fails)
            @test pyconvert(Vector{Bool}, pyeval("LAST_APP.silence_updates", Main)) ==
                (owned ? [true] : [true, false])
            if !before_load
                @test pyconvert(Bool, pyeval("LAST_APP.silence_at_load", Main))
            end
            if port == 0
                @test pyconvert(Bool, pyeval("LAST_LAUNCH['silence']", Main))
            else
                @test pyconvert(Dict{String, Any}, pyeval("LAST_CONNECT", Main)) ==
                    Dict("host" => "127.0.0.1", "port" => port, "timeout" => 30)
            end
            if scenario == "stale"
                @test pyconvert(Bool, pyeval("LAST_APP.existing.unloaded", Main))
            end
            if scenario == "moved"
                @test occursin("load geometry", sprint(showerror, caught))
                @test occursin("3168", read(joinpath(output, "load-relocations.txt"), String))
                @test !pyconvert(Bool, pyeval("LAST_APP._project.ran or LAST_APP._project.compiled", Main))
            end
            if mode == "transient" && scenario == "session0"
                @test occursin("set pscad_port", sprint(showerror, caught))
                @test occursin("session 0", sprint(showerror, caught))
            end
            if scenario in ("success", "build_error", "build_error_throw")
                @test read(joinpath(output, "generated.log")) == read(joinpath(directory, "generated.log"))
                @test read(joinpath(output, "generated.gf81", "cable.log")) ==
                    read(joinpath(directory, "generated.gf81", "cable.log"))
            end
            if scenario in ("build_error", "build_error_throw")
                @test occursin("error\tBuild\tgenerated", read(joinpath(output, "build-messages.txt"), String))
                failure = TOML.parsefile(joinpath(output, "failure.toml"))
                @test failure["kind"] == "line_constants_error"
                @test failure["errors"] == ["An alert was issued while solving line constants, check the log file 'generated.log'",
                    "No Tline ERROR block found in the referenced log."]
                @test failure["stage"] == ""
                @test failure["log_path"] == "generated.log"
                diagnostic = sprint(showerror, caught)
                @test occursin("An alert was issued while solving line constants", diagnostic)
                @test !occursin("omitted first diagnostic", diagnostic)
                @test !occursin("synthetic cable log ending", diagnostic)
                @test occursin("An alert was issued while solving line constants",
                    read(joinpath(output, "pscad-console.txt"), String))
                if scenario == "build_error_throw"
                    @test occursin(mode == "transient" ? "synthetic run failure" :
                        "synthetic compile failure", diagnostic)
                end
            end
            if mode == "transient" && scenario == "run_restore_error"
                @test occursin("synthetic run failure", sprint(showerror, caught))
                @test !occursin("synthetic silence restore failure", sprint(showerror, caught))
            end
        end
    end
    for silence in (false, true)
        mktempdir(root) do directory
            environment = pyimport("os").environ
            environment["LCM_RUNNER_FIXTURE_OUTPUT"] = joinpath(directory, "outputs")
            environment["LCM_RUNNER_SILENCE"] = string(silence)
            try
                @test main(runner_arguments(directory, joinpath(directory, "outputs"); port=50000, owned=false)) === nothing
                @test pyconvert(Bool, pyeval("LAST_APP.silence", Main)) == silence
                @test pyconvert(Vector{Bool}, pyeval("LAST_APP.silence_updates", Main)) == [true, silence]
                @test pyconvert(Bool, pyeval("LAST_APP.restore_after_unload", Main))
            finally
                environment.pop("LCM_RUNNER_FIXTURE_OUTPUT", nothing)
                environment.pop("LCM_RUNNER_SILENCE", nothing)
            end
        end
    end
end

module PSCADProbeProtocol
using Test, PythonCall

const PYTHON_PROBE_FIXTURE = raw"""
import atexit
import runpy
import sys
import threading
import time
import types

source, host, port, limit, mode, trace = sys.argv[1:]
def shutdown():
    with open(trace + ".shutdown", "w") as stream:
        stream.write("interpreter shutdown")
atexit.register(shutdown)

def record(call):
    with open(trace, "a") as stream:
        stream.write(call + "\n")
    if mode == call:
        while True:
            time.sleep(1)
    if mode in ("error", "shutdown_error"):
        raise RuntimeError("synthetic remote-call failure")

class Application:
    __slots__ = ()
    def is_busy(self):
        record("is_busy")
        return mode == "shutdown_success"
    @property
    def version(self):
        record("version")
        return "5.1.0"

def connect(**kwargs):
    assert kwargs == {"host": host, "port": int(port), "timeout": float(limit)}
    if mode in ("shutdown_success", "shutdown_error"):
        threading.Thread(target=lambda: time.sleep(30), daemon=False).start()
    record("connect")
    if mode == "stderr_success":
        print("synthetic automation diagnostic", file=sys.stderr, flush=True)
    return Application()

mhi = types.ModuleType("mhi")
mhi.pscad = types.ModuleType("mhi.pscad")
mhi.pscad.connect = connect
sys.modules["mhi"] = mhi
sys.modules["mhi.pscad"] = mhi.pscad
sys.argv = [source, host, port, limit]
runpy.run_path(source, run_name="__main__")
"""

@testset "PSCAD responsiveness probe without station" begin
    python = pyconvert(String, pyimport("sys").executable)
    probe = normpath(joinpath(@__DIR__, "../../../src/pscad/remote/probe.py"))
    for mode in ("success", "shutdown_success", "stderr_success",
        "connect", "is_busy", "version", "error",
        "shutdown_error")
        mktempdir(mkpath(joinpath(@__DIR__, "../../../.linecablemodels/test/station-protocol"))) do root
            trace = joinpath(root, "calls.txt")
            console = joinpath(root, "pscad-probe.txt")
            command = `$python -u -c $PYTHON_PROBE_FIXTURE $probe 127.0.0.1 50000 1 $mode $trace`
            process = open(console, "w") do io
                run(pipeline(ignorestatus(command); stdout = io, stderr = io); wait = false)
            end
            try
                @test timedwait(() -> process_exited(process), 10) == :ok
            finally
                process_running(process) && kill(process)
                wait(process)
            end
            output = read(console, String)
            calls = readlines(trace)
            @test !isfile(trace * ".shutdown")
            @test occursin("Probing PSCAD on 127.0.0.1:50000 (limit=1 s)", output)
            if mode in ("success", "shutdown_success", "stderr_success")
                @test process.exitcode == 0
                @test calls == ["connect", "is_busy", "version"]
                busy = mode == "shutdown_success" ? "True" : "False"
                @test occursin("PSCAD remote calls answered: busy=$busy, version=5.1.0", output)
                @test !occursin("probe exceeded", output)
                @test occursin("synthetic automation diagnostic", output) ==
                      (mode == "stderr_success")
            elseif mode in ("error", "shutdown_error")
                @test process.exitcode == 1
                @test calls == ["connect"]
                @test occursin("synthetic remote-call failure", output)
                @test !occursin("probe exceeded", output)
            else
                @test process.exitcode == 124
                @test last(calls) == mode
                @test occursin("PSCAD remote-call probe exceeded 1 s", output)
                @test !occursin("PSCAD remote calls answered", output)
            end
        end
    end
end
end

module PSCADLifecycleProtocol
using Test, TOML

const POWERSHELL_LIFECYCLE_FIXTURE = raw"""
param([string] $Lifecycle, [string] $LocalCase, [string] $Mode, [string] $EntryPoint)
. $Lifecycle
$ownerPath = Join-Path $LocalCase "pscad-owner.txt"
$exe = Join-Path $LocalCase "Pscad.exe"
$created = [DateTime]::Parse("2026-10-08T12:00:00Z").ToUniversalTime()
Set-Content -LiteralPath $ownerPath -Value @(5308, $exe, 54330, $created.ToString('o'))
$application = [PSCustomObject] @{
    ProcessId = 5308
    ExecutablePath = $exe
    CommandLine = 'Pscad.exe /port:54330'
    CreationDate = $created
    ExecutionState = 3
}
if ($EntryPoint -eq 'runner') {
    $ownerPath = Join-Path $LocalCase 'owner.txt'
    $runner = "$LocalCase\invoke-runner.cmd"
    Set-Content -LiteralPath $ownerPath -Value @(5308, $runner)
    $application.CommandLine = $runner
} elseif ($EntryPoint -eq 'simulation') {
    $application.ExecutablePath = "$LocalCase\fixture.exe"
    Remove-Item -LiteralPath $ownerPath
}
$script:frames = @($application)
switch ($Mode) {
    'gone' { $script:frames = @($null) }
    'terminated' { $application.ExecutionState = 7; $script:frames += $null }
    'empty_executable' { $application.ExecutablePath = ''; $script:frames += $null }
    'empty_command' { $application.CommandLine = ' '; $script:frames += $null }
    'empty_creation' { $application.CreationDate = $null; $script:frames += $null }
    'persistent_empty' { $application.ExecutablePath = $null }
    'recovers' {
        $incomplete = [PSCustomObject] @{
            ProcessId = 5308; ExecutablePath = ''; CommandLine = 'Pscad.exe /port:54330'
            CreationDate = $created; ExecutionState = 3
        }
        $script:frames = @($incomplete, $application, $null)
    }
    'different_executable' { $application.ExecutablePath = Join-Path $LocalCase 'other.exe' }
    'different_creation' { $application.CreationDate = $created.AddSeconds(1) }
    'different_port' { $application.CommandLine = 'Pscad.exe /port:54331' }
    'stop_gone' { $script:frames += $null }
    'stop_throw_gone' { $script:frames += $null }
    'stop_exiting' {
        $script:frames += [PSCustomObject] @{
            ProcessId = 5308; ExecutablePath = ''; CommandLine = $null
            CreationDate = $null; ExecutionState = 7
        }
    }
    'stop_different_executable' {
        $script:frames += [PSCustomObject] @{
            ProcessId = 5308; ExecutablePath = 'other.exe'; CommandLine = 'other.exe'
            CreationDate = $created; ExecutionState = 3
        }
    }
    'stop_different_creation' {
        $script:frames += [PSCustomObject] @{
            ProcessId = 5308; ExecutablePath = $exe; CommandLine = 'Pscad.exe /port:54330'
            CreationDate = $created.AddSeconds(1); ExecutionState = 3
        }
    }
}
$script:queries = 0
$script:kills = 0
$script:afterStopQueries = 0
function Get-CimInstance {
    param([string] $ClassName, [string] $Filter)
    if ([string]::IsNullOrEmpty($Filter)) {
        if ($EntryPoint -eq 'simulation') { return $application }
        return
    }
    if ($ClassName -ne 'Win32_Process' -or $Filter -ne 'ProcessId = 5308') {
        throw "Unexpected process query: $ClassName $Filter"
    }
    if ($Mode -eq 'stop_delayed_gone') {
        $script:queries++
        if ($script:kills -eq 0) { return $application }
        $script:afterStopQueries++
        if ($script:afterStopQueries -le 2) { return $application }
        return
    }
    $index = [Math]::Min($script:queries, $script:frames.Count - 1)
    $script:queries++
    return $script:frames[$index]
}
function taskkill.exe {
    if (($args -join ' ') -ne '/PID 5308 /T /F') { throw "Unexpected stop: $args" }
    $script:kills++
    if ($Mode.StartsWith('stop_throw_')) { throw 'synthetic stop command failure' }
    if ($ErrorActionPreference -ne 'Continue') { throw 'Native command requires Continue' }
    Write-Output 'synthetic taskkill output'
    Write-Error 'synthetic taskkill stderr'
    $code = if ($Mode.StartsWith('stop_')) { 128 } else { 0 }
    $global:LASTEXITCODE = $code
}
$exitCode = 0
$caught = ''
try {
    # Represents a runner that has already returned a successful exit code.
    try {
    } finally {
        $grace = if ($Mode -eq 'stop_delayed_gone') { 2 } elseif ($Mode -eq 'matched' -or $Mode.StartsWith('stop_')) { 0 } elseif ($Mode -eq 'persistent_empty') { 1 } else { 2 }
        if ($EntryPoint -eq 'case') {
            Stop-PSCADCase -LocalCase $LocalCase -GraceSeconds $grace
        } elseif ($EntryPoint -eq 'runner') {
            Stop-OwnedRunner -OwnerPath $ownerPath -OwnedRoot $LocalCase -GraceSeconds $grace
        } elseif ($EntryPoint -eq 'simulation') {
            Stop-EMTDC -CaseDirectory $LocalCase -GraceSeconds $grace
        } else {
            Stop-OwnedPSCAD -OwnerPath $ownerPath -GraceSeconds $grace
        }
    }
} catch {
    $exitCode = 1
    $caught = $_.Exception.Message
    [Console]::Error.WriteLine($caught)
}
Set-Content -LiteralPath (Join-Path $LocalCase 'observed.toml') -Value @(
    "queries = $script:queries", "kills = $script:kills", "exit_code = $exitCode",
    "owner_exists = $((Test-Path -LiteralPath $ownerPath).ToString().ToLowerInvariant())",
    "caught = '$caught'"
)
exit $exitCode
"""

@testset "PSCAD owned-process cleanup without station" begin
    powershell = get(ENV, "LINECABLEMODELS_PSCAD_TEST_POWERSHELL", "")
    if isempty(powershell)
        @test_skip "lifecycle protocol requires a local PowerShell interpreter"
    else
        lifecycle = normpath(joinpath(@__DIR__, "../../../src/pscad/remote/lifecycle.ps1"))
        for entry in ("final", "case"),
            mode in ("gone", "terminated", "empty_executable", "empty_command",
                "empty_creation", "persistent_empty", "recovers", "different_executable",
                "different_creation", "different_port", "matched", "stop_gone", "stop_exiting",
                "stop_different_executable", "stop_different_creation", "stop_alive",
                "stop_throw_gone", "stop_throw_alive", "stop_delayed_gone")

            mktempdir(mkpath(joinpath(@__DIR__, "../../../.linecablemodels/test/station-protocol"))) do root
                fixture = joinpath(root, "lifecycle-fixture.ps1")
                write(fixture, POWERSHELL_LIFECYCLE_FIXTURE)
                console = joinpath(root, "console.txt")
                command = `$powershell -NoLogo -NoProfile -NonInteractive -File $fixture -Lifecycle $lifecycle -LocalCase $root -Mode $mode -EntryPoint $entry`
                process = open(console, "w") do io
                    run(pipeline(ignorestatus(command); stdout = io, stderr = io); wait = false)
                end
                try
                    @test timedwait(() -> process_exited(process), 15) == :ok
                finally
                    process_running(process) && kill(process)
                    wait(process)
                end
                stop_failed = mode in ("stop_alive", "stop_throw_alive")
                @test process.exitcode == Int(stop_failed)
                observed = TOML.parsefile(joinpath(root, "observed.toml"))
                @test observed["exit_code"] == Int(stop_failed)
                @test observed["owner_exists"] == stop_failed
                @test observed["kills"] ==
                      Int(mode == "matched" || startswith(mode, "stop_"))
                expected_queries = entry == "final" && mode == "recovers" ? 3 :
                                   startswith(mode, "stop_") ? 2 :
                                   entry == "final" &&
                                   mode in ("terminated", "empty_executable",
                    "empty_command", "empty_creation", "persistent_empty") ? 2 : 1
                @test observed["queries"] >= expected_queries
                console_text = read(console, String)
                killed = mode == "matched" || startswith(mode, "stop_")
                @test occursin("taskkill PID=5308 exit code=", console_text) == killed
                if killed && !startswith(mode, "stop_throw_")
                    @test occursin("synthetic taskkill output", console_text)
                    @test occursin("synthetic taskkill stderr", console_text)
                    @test isfile(joinpath(root, "result", "taskkill-5308-stdout.txt"))
                    @test isfile(joinpath(root, "result", "taskkill-5308-stderr.txt"))
                end
                mode == "stop_delayed_gone" && @test observed["queries"] >= 5
                mismatch = mode in
                           ("different_executable", "different_creation", "different_port",
                    "stop_different_executable", "stop_different_creation")
                @test occursin("Warning: PID 5308", console_text) == mismatch
                @test occursin("no longer belongs to the recorded owned PSCAD", console_text) ==
                      (mode in ("different_executable", "different_creation",
                    "stop_different_executable", "stop_different_creation"))
                @test occursin("Could not stop owned PSCAD process tree 5308", observed["caught"]) ==
                      stop_failed
            end
        end
        for entry in ("runner", "simulation"),
            mode in ("matched", "stop_delayed_gone", "stop_alive")

            mktempdir(mkpath(joinpath(@__DIR__, "../../../.linecablemodels/test/station-protocol"))) do root
                fixture = joinpath(root, "lifecycle-fixture.ps1")
                write(fixture, POWERSHELL_LIFECYCLE_FIXTURE)
                console = joinpath(root, "console.txt")
                command = `$powershell -NoLogo -NoProfile -NonInteractive -File $fixture -Lifecycle $lifecycle -LocalCase $root -Mode $mode -EntryPoint $entry`
                process = open(console, "w") do io
                    run(pipeline(ignorestatus(command); stdout = io, stderr = io); wait = false)
                end
                try
                    @test timedwait(() -> process_exited(process), 15) == :ok
                finally
                    process_running(process) && kill(process)
                    wait(process)
                end
                @test process.exitcode == Int(mode == "stop_alive")
                observed = TOML.parsefile(joinpath(root, "observed.toml"))
                @test observed["kills"] == 1
                mode == "stop_delayed_gone" && @test observed["queries"] >= 3
                @test occursin("synthetic taskkill stderr", read(console, String))
                @test occursin("Could not stop", observed["caught"]) ==
                      (mode == "stop_alive")
            end
        end
    end
end
end

module PSCADSupervisorProtocol
using Test, TOML, PythonCall

const POWERSHELL_SUPERVISOR_FIXTURE = raw"""
param([string] $Supervisor, [string] $Lifecycle, [string] $Root, [string] $Python, [string] $Mode)
$shared = Join-Path $Root 'shared'
$local = Join-Path $Root 'scratch'
$toolkit = Join-Path $shared 'toolkit'
$harness = Join-Path $Root 'harness'
New-Item -ItemType Directory -Path $toolkit, $harness -Force | Out-Null
Copy-Item -LiteralPath $Supervisor -Destination (Join-Path $harness 'supervisor.ps1')
# The actual supervisor runs with local process and Windows service mocks.
Set-Content -LiteralPath (Join-Path $harness 'lifecycle.ps1') -Value @'
. $Lifecycle
$script:cleanupCalls = [Collections.Generic.List[string]]::new()
$script:emtdcCalls = 0
function Get-CimInstance {
    if ($Mode -eq 'desktop_start' -and $script:desktopStarted) {
        return [PSCustomObject] @{
            ProcessId = 17284; ExecutablePath = $Python
            CreationDate = [DateTime]::UtcNow; CommandLine = 'Pscad.exe /port:54330'
        }
    }
}
function Get-NetTCPConnection {
    if ($Mode -eq 'desktop_start' -and -not $script:desktopStarted) { return }
    $script:pscadOwned = $Mode -ne 'stderr_success'
    [PSCustomObject] @{ OwningProcess = 17284 }
}
function New-ScheduledTaskAction { param($Execute, $Argument); return $Execute }
function New-ScheduledTaskPrincipal { param($UserId, $LogonType, $RunLevel); return $UserId }
function New-ScheduledTaskSettingsSet { return }
function Register-ScheduledTask {
    Set-Content -LiteralPath (Join-Path $SharedCase 'task-registered.txt') -Value $true
}
function Start-ScheduledTask { $script:desktopStarted = $true }
function Unregister-ScheduledTask {
    Set-Content -LiteralPath (Join-Path $SharedCase 'task-removed.txt') -Value $true
}
function Invoke-PSCADNative {
    param($Executable, $Arguments, $StdoutPath, $StderrPath)
    if ($Mode -eq 'desktop_start' -and $Arguments[0] -eq '-c') {
        Set-Content -LiteralPath $StdoutPath -Value $Python
        return 0
    }
    & $Executable @Arguments 1> $StdoutPath 2> $StderrPath
    return $LASTEXITCODE
}
function Stop-EMTDC {
    param([string] $CaseDirectory)
    $script:cleanupCalls.Add('EMTDC')
    $script:emtdcCalls++
    if ($script:emtdcCalls -gt 1 -and $Mode -eq 'probe_failure_cleanup_failure') {
        throw 'synthetic EMTDC cleanup failure'
    }
}
function Stop-OwnedPSCAD {
    param([string] $OwnerPath, [int] $GraceSeconds)
    $script:cleanupCalls.Add('PSCAD')
    if ($script:pscadOwned -and $Mode -ne 'desktop_start') { throw 'synthetic owned PSCAD cleanup failure' }
}
function Stop-OwnedRunner {
    param([string] $OwnerPath, [string] $OwnedRoot)
    $script:cleanupCalls.Add('runner')
    if ($Mode -eq 'runner_primary_cleanup_failure' -and (Test-Path -LiteralPath $OwnerPath)) { throw 'synthetic runner cleanup failure' }
}
function Start-Process {
    param($FilePath, $ArgumentList, [switch] $PassThru, $WindowStyle)
    Set-Content -LiteralPath (Join-Path $LocalCase 'result/exit-code.txt') -Value '0'
    Set-Content -LiteralPath (Join-Path $LocalCase 'result/stdout.txt') -Value 'mock runner completed'
    Set-Content -LiteralPath (Join-Path $LocalCase 'result/stderr.txt') -Value ''
    Set-Content -LiteralPath (Join-Path $SharedCase 'runner-started.txt') -Value 'started'
    $process = [PSCustomObject] @{ Id = 42; HasExited = $Mode -ne 'runner_primary_cleanup_failure' }
    $process | Add-Member ScriptMethod WaitForExit { return }
    $process | Add-Member ScriptMethod Refresh { throw 'synthetic primary runner failure' }
    return $process
}
function Copy-Result {
    param([string] $Source, [string] $Destination)
    $script:cleanupCalls.Add('copy')
    if ($Mode -eq 'probe_failure_cleanup_failure') { throw 'synthetic output copy failure' }
    New-Item -ItemType Directory -Path $Destination -Force | Out-Null
    Get-ChildItem -LiteralPath $Source -Force | Copy-Item -Destination $Destination -Recurse -Force
}
function Remove-Item {
    param([string] $LiteralPath, [switch] $Recurse, [switch] $Force, $ErrorAction)
    # Record all cleanup attempts even if earlier steps failed.
    if ($LiteralPath.EndsWith('owner.txt')) {
        Set-Content -LiteralPath (Join-Path $SharedCase 'cleanup.txt') -Value $script:cleanupCalls
    }
    Microsoft.PowerShell.Management\Remove-Item @PSBoundParameters
}
'@
# Override Copy-Result after the supervisor has defined it, at the first port query.
$wrapper = Get-Content -LiteralPath (Join-Path $harness 'lifecycle.ps1') -Raw
$wrapper = $wrapper.Replace('function Get-NetTCPConnection {', @'
function Get-NetTCPConnection {
    Set-Item Function:script:Copy-Result -Value {
        param([string] $Source, [string] $Destination)
        $script:cleanupCalls.Add('copy')
        if ($Mode -eq 'probe_failure_cleanup_failure') { throw 'synthetic output copy failure' }
        New-Item -ItemType Directory -Path $Destination -Force | Out-Null
        Get-ChildItem -LiteralPath $Source -Force | Copy-Item -Destination $Destination -Recurse -Force
    }
'@)
Set-Content -LiteralPath (Join-Path $harness 'lifecycle.ps1') -Value $wrapper
Set-Content -LiteralPath (Join-Path $toolkit 'runner.jl') -Value '# mock runner'
Set-Content -LiteralPath (Join-Path $toolkit 'Project.toml') -Value ''
Set-Content -LiteralPath (Join-Path $shared 'generated.pscx') -Value 'fixture'
$probeExit = if ($Mode -eq 'probe_failure_cleanup_failure') { 1 } else { 0 }
Set-Content -LiteralPath (Join-Path $toolkit 'probe.py') -Value @(
    'import os, sys',
    'print("PSCAD remote calls answered: busy=True, version=5.1.0", flush=True)',
    'print("synthetic native stderr", file=sys.stderr, flush=True)',
    "os._exit($probeExit)"
)
& (Join-Path $harness 'supervisor.ps1') -SharedCase $shared -LocalCase $local -RemoteRoot $Root -Julia $Python -Python $Python -Formulation fixture -FrequencyStart 1 -FrequencyEnd 10 -FrequencyIncrements 1 -Mode constants -ProjectName fixture -OutputStem pi -PSCADVersion 5.1.0 -PSCADPort 54330 -Verbosity 0 -TimeoutSeconds 5
exit $LASTEXITCODE
"""

const POWERSHELL_CASE_FAILURE_FIXTURE = raw"""
param([string] $Lifecycle, [string] $Root, [string] $Python)
. $Lifecycle
$ErrorActionPreference = 'Stop'
$stdout = Join-Path $Root 'native-stdout.txt'
$stderr = Join-Path $Root 'native-stderr.txt'
$code = Invoke-PSCADNative -Executable $Python -Arguments @('-c', 'import sys; print("Pscad.exe"); print("synthetic find_exe stderr", file=sys.stderr)') -StdoutPath $stdout -StderrPath $stderr
if ($code -ne 0 -or $ErrorActionPreference -ne 'Stop') { throw 'Native invocation changed caller error handling' }
$script:stages = [Collections.Generic.List[string]]::new()
function Stop-OwnedRunner { $script:stages.Add('runner'); throw 'primary runner stop failure' }
function Stop-EMTDC { $script:stages.Add('EMTDC'); throw 'secondary EMTDC stop failure' }
function Stop-OwnedPSCAD { $script:stages.Add('PSCAD'); throw 'secondary PSCAD stop failure' }
try {
    Stop-PSCADCase -LocalCase $Root
    throw 'Expected cleanup failure'
} catch {
    [Console]::Error.WriteLine($_.Exception.Message)
    Set-Content -LiteralPath (Join-Path $Root 'stages.txt') -Value $script:stages
}
"""

@testset "PSCAD supervisor primary errors and native stderr without station" begin
    powershell = get(ENV, "LINECABLEMODELS_PSCAD_TEST_POWERSHELL", "")
    if isempty(powershell)
        @test_skip "supervisor protocol requires a local PowerShell interpreter"
    else
        supervisor = normpath(joinpath(@__DIR__, "../../../src/pscad/remote/supervisor.ps1"))
        lifecycle = normpath(joinpath(@__DIR__, "../../../src/pscad/remote/lifecycle.ps1"))
        python = pyconvert(String, pyimport("sys").executable)
        modes = Sys.iswindows() ? ("stderr_success", "desktop_start", "probe_failure_cleanup_failure", "runner_primary_cleanup_failure") :
            ("stderr_success", "probe_failure_cleanup_failure", "runner_primary_cleanup_failure")
        for mode in modes
            mktempdir(mkpath(joinpath(@__DIR__, "../../../.linecablemodels/test/station-protocol"))) do root
                fixture = joinpath(root, "supervisor-fixture.ps1")
                write(fixture, POWERSHELL_SUPERVISOR_FIXTURE)
                console = joinpath(root, "console.txt")
                command = `$powershell -NoLogo -NoProfile -NonInteractive -File $fixture -Supervisor $supervisor -Lifecycle $lifecycle -Root $root -Python $python -Mode $mode`
                process = open(console, "w") do io
                    run(pipeline(ignorestatus(command); stdout = io, stderr = io); wait = false)
                end
                try
                    @test timedwait(() -> process_exited(process), 15) == :ok
                finally
                    process_running(process) && kill(process)
                    wait(process)
                end
                text = read(console, String)
                @test process.exitcode == Int(mode ∉ ("stderr_success", "desktop_start"))
                shared = joinpath(root, "shared")
                @test isfile(joinpath(shared, "runner-started.txt")) ==
                      (mode != "probe_failure_cleanup_failure")
                if mode != "probe_failure_cleanup_failure"
                    output = joinpath(shared, "outputs")
                    @test occursin("PSCAD remote calls answered",
                        read(joinpath(output, "pscad-probe.txt"), String))
                    @test occursin("synthetic native stderr",
                        read(joinpath(output, "pscad-probe-stderr.txt"), String))
                end
                if mode == "desktop_start"
                    @test isfile(joinpath(shared, "task-registered.txt"))
                    @test isfile(joinpath(shared, "task-removed.txt"))
                end
                stages = readlines(joinpath(shared, "cleanup.txt"))
                @test stages[end] == "copy"
                if mode in ("stderr_success", "desktop_start")
                    @test !occursin("Cleanup failure:", text)
                    @test !occursin("not answering remote calls", text)
                else
                    primary = mode == "probe_failure_cleanup_failure" ?
                              "PSCAD on port 54330 is not answering remote calls within 60 s" :
                              "synthetic primary runner failure"
                    primary_at = findfirst(primary, text)
                    cleanup_at = findfirst("Cleanup failure:", text)
                    @test primary_at !== nothing
                    @test cleanup_at !== nothing
                    @test first(primary_at) < first(cleanup_at)
                    @test occursin("synthetic owned PSCAD cleanup failure", text)
                    @test "PSCAD" in stages
                    if mode == "probe_failure_cleanup_failure"
                        @test occursin("synthetic EMTDC cleanup failure", text)
                        @test occursin("synthetic output copy failure", text)
                    else
                        @test occursin("synthetic runner cleanup failure", text)
                        @test "runner" in stages
                    end
                end
            end
        end
        mktempdir(mkpath(joinpath(@__DIR__, "../../../.linecablemodels/test/station-protocol"))) do root
            fixture = joinpath(root, "case-failure-fixture.ps1")
            write(fixture, POWERSHELL_CASE_FAILURE_FIXTURE)
            console = joinpath(root, "console.txt")
            command = `$powershell -NoLogo -NoProfile -NonInteractive -File $fixture -Lifecycle $lifecycle -Root $root -Python $python`
            process = open(console, "w") do io
                run(pipeline(ignorestatus(command); stdout = io, stderr = io); wait = false)
            end
            try
                @test timedwait(() -> process_exited(process), 15) == :ok
            finally
                process_running(process) && kill(process)
                wait(process)
            end
            @test process.exitcode == 0
            @test readlines(joinpath(root, "stages.txt")) == ["runner", "EMTDC", "PSCAD"]
            text = read(console, String)
            @test first(findfirst("primary runner stop failure", text)) <
                  first(findfirst("secondary EMTDC stop failure", text)) <
                  first(findfirst("secondary PSCAD stop failure", text))
            @test occursin("Pscad.exe", read(joinpath(root, "native-stdout.txt"), String))
            @test occursin("synthetic find_exe stderr", read(joinpath(root, "native-stderr.txt"), String))
        end
    end
end
end


module PSCADIdentityProtocol
using Test, PythonCall

@testset "PSCAD raw Fortran identity with mock installation" begin
    directory = mkpath(joinpath(@__DIR__, "../../../.linecablemodels/test/station-identity"))
    mktempdir(directory) do root
        executable = joinpath(mkpath(joinpath(root, "bin", "win64")), "Pscad.exe")
        solver = joinpath(root, "line-constants.exe")
        write(executable, "fixture executable")
        write(solver, "fixture line constants")
        write(joinpath(root, "master.pslx"), "fixture master library")
        globals = pydict()
        globals["EXECUTABLE"] = executable
        globals["SOLVER"] = solver
        globals["IDENTITY_SOURCE"] = normpath(joinpath(@__DIR__, "../../../src/pscad/remote/identity.py"))
        pyexec(raw"""
import runpy
import sys
import types
common = types.ModuleType("mhi.common")
common.process = types.SimpleNamespace(find_exe=lambda *args, **kw: EXECUTABLE)
sys.modules["mhi.common"] = common
import mhi.pscad
identify = runpy._lcm_original_run_path(IDENTITY_SOURCE)["identify"]
closed = []
launches = []
def detected_software_error():
    raise SystemError("Unable to retrieve detected software")
app = types.SimpleNamespace(version="5.1.0",
    _settings=lambda kw: {"fortran_version": "GF81x64", "file_lcp": SOLVER},
    settings=detected_software_error, quit=lambda: closed.append(True))
def launch(**kw):
    launches.append(kw)
    return app
mhi.pscad.launch = launch
attached = identify("5.1.0", app)
launched = identify("5.1.0")
""", globals)
        attached = pyconvert(Dict{String, String}, globals["attached"])
        @test attached == pyconvert(Dict{String, String}, globals["launched"])
        @test attached["fortran_version"] == "GF81x64"
        @test attached["schema"] == "1"
        @test attached["line_constants_path"] == realpath(solver)
        @test all(occursin(r"^[0-9a-f]{64}$", attached[name * "_sha256"])
            for name in ("pscad", "line_constants", "master_library"))
        @test pyconvert(Int, pybuiltins.len(globals["closed"])) == 1
        @test pyconvert(Int, pybuiltins.len(globals["launches"])) == 1
        pyexec(raw"""
app._settings = lambda kw: {"file_lcp": SOLVER}
unknown = identify("5.1.0", app)
""", globals)
        @test pyconvert(String, globals["unknown"]["fortran_version"]) == "unknown"
        @test_throws PyException globals["identify"]("5.0.2", globals["app"])
    end
end
end

@testset "PSCAD stable output sizes and buffered writes" begin
    directory = mkpath(joinpath(REPOSITORY_ROOT, ".linecablemodels", "test", "station-outputs"))
    mktempdir(directory) do root
        write(joinpath(root, "large.inf"), "PGB(1) Desc=\"Voltage\" Group=\"\"\n")
        path = joinpath(root, "large_01.out")
        write(path, "0.0 1.0D+0\n")
        @test _wait_outputs([root], "large", 1; timeout=5, poll_interval=0.001).outs == [path]
        rows = 1_000_000
        open(path, "w") do io
            for _ in 1:rows
                println(io, "0.0 1.0D+0")
            end
        end
        @test _wait_outputs([root], "large", rows; timeout=0.1, poll_interval=0.001).outs == [path]
        @test_throws r"rows per file: \[1000000\]; expected 1000001" _wait_outputs(
            [root], "large", rows + 1; timeout=0.1, poll_interval=0.001)
        write(path, "0.0 1.0D+0\n")
        writer = @async begin
            sleep(0.2)
            open(path, "a") do io
                println(io, "0.1 2.0d+0")
            end
        end
        try
            @test _wait_outputs([root], "large", 2; timeout=5, poll_interval=0.01).outs == [path]
        finally
            wait(writer)
        end
    end
end
