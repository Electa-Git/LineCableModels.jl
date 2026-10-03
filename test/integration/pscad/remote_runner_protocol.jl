using Test
using PythonCall
using TOML

const REPOSITORY_ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
const REMOTE_ROOT=joinpath(
    REPOSITORY_ROOT,"ext","LineCableModelsPSCADExt","remote")

include(joinpath(REMOTE_ROOT,"runner.jl"))

const PYTHON_AUTOMATION_FIXTURE=raw"""
import importlib.metadata
import os
import runpy
import sys
import types

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
        if os.environ.get("LCM_RUNNER_FIXTURE_FAIL") == "1":
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

    def find_all(self, kind):
        return [self.line] if kind == "TLine" else []

    def canvas(self, name):
        return self._canvas

    def messages(self):
        return ["fixture message"]

    def output(self):
        return "fixture project output"

    def save(self):
        pass

    def unload(self):
        self.unloaded = True

class Application:
    version = "5.1.0"

    def __init__(self):
        self._project = Project()
        self.closed = False

    def licensed(self):
        return True

    def load(self, path):
        pass

    def project(self, name):
        return self._project

    def quit(self):
        self.closed = True

pscad = types.ModuleType("mhi.pscad")
pscad.launch = lambda **kwargs: Application()
mhi = types.ModuleType("mhi")
mhi.pscad = pscad
sys.modules["mhi"] = mhi
sys.modules["mhi.pscad"] = pscad
importlib.metadata.version = lambda name: "3.1.2" if name == "mhi.pscad" else "fixture"
runpy.run_path = lambda path: {"identify": lambda version, app: dict(IDENTITY)}
"""

function runner_arguments(directory,output;verbosity="0")
    project=joinpath(directory,"generated.pscx")
    write(project,"fixture")
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
        "5.1.0",verbosity]
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

    mktempdir() do directory
        output=joinpath(directory,"outputs")
        arguments=runner_arguments(directory,output)
        @test_throws ArgumentError main(setindex!(copy(arguments),"3",10))
        @test_throws ArgumentError main(setindex!(copy(arguments),"5.0.2",9))
        @test_throws ArgumentError main(setindex!(copy(arguments),"bad-stem",4))

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
        @test occursin("Starting PSCAD line-constants calculation",console)
        @test occursin("Collected detailed PSCAD Z and Y outputs",console)
    end

    mktempdir() do directory
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
        @test occursin("line-constants calculation",sprint(showerror,caught))
        @test occursin("synthetic compile failure",read(joinpath(output,"pscad-console.txt"),String))
    end
end
