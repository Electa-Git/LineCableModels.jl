using Test
include("fixture_container_engine.jl")
const FixtureEngine = FixtureContainerEngine

@testset "test engine launcher preserves arguments and isolates connection settings" begin
    @test_throws ArgumentError FixtureEngine.launcher(String[], Dict())
    @test_throws ArgumentError FixtureEngine.launcher(["podman"], Dict())
    for command in (["/usr/bin/podman", "--remote=false"],
            ["/usr/bin/docker", "--host", "unix:///run/user/1000/docker.sock"])
        script = FixtureEngine.launcher(command, Dict("LC_ALL" => "C"))
        @test occursin(join(FixtureEngine.shell_word.(command), " "), script)
    end
    mktempdir() do directory
        path = joinpath(directory, "engine")
        # Report only test-owned argument/environment evidence, never host data.
        probe = raw"""
        printf '%s\000' "$@" "${DOCKER_HOST-absent}" "${CONTAINER_HOST-absent}" "${DOCKER_CONTEXT-absent}" "${LCM_FIXTURE_SECRET-absent}" "$LC_ALL"
        """
        write(path, FixtureEngine.launcher(["/bin/sh", "-c", probe, "engine"], Dict("LC_ALL" => "C")))
        chmod(path, 0o700)
        arguments = ["two words", "single'quote", raw"$literal", "", "--dash", "line\nbreak"]
        output = withenv("DOCKER_HOST" => "tcp://foreign.invalid:2375",
                "CONTAINER_HOST" => "ssh://foreign.invalid", "DOCKER_CONTEXT" => "foreign",
                "LCM_FIXTURE_SECRET" => "must-not-be-inherited") do
            read(Cmd([path; arguments]), String)
        end
        @test split(output, '\0')[1:end-1] == [arguments; fill("absent", 4); "C"]
        @test_throws ArgumentError FixtureEngine.main(["docker", path])
        @test_throws ArgumentError FixtureEngine.main(["containerd", joinpath(directory, "new")])
        @test_throws ArgumentError FixtureEngine.main(["podman"])
    end
end

@testset "legacy fixture cleanup preserves failure and owns its targets" begin
    integration = normpath(joinpath(@__DIR__, "..", "..", "test", "integration"))
    for harness in ("run.sh", "run-artifacts.sh")
        source = read(joinpath(integration, harness), String)
        cleanup = match(r"(?ms)^cleanup\(\) \{\n.*?^\}", source)
        @test cleanup !== nothing
        # Execute the actual cleanup body with a fake engine and fake rm. No
        # container or file is touched; unexpected targets fail the fixture.
        setup = raw"""
        set -Eeuo pipefail
        CONTAINER_RUNTIME=fixture_engine
        CONTAINER_NAME=fixture-owned
        broker_started=true
        started=true
        cli_created=true
        worker_data=/tmp/lcm-fixture-scratch
        certificate_parent=/tmp/lcm-fixture-certificates
        fixture_engine() {
            if [[ "$1" == inspect ]]; then printf 'false\n'; return 0; fi
            [[ "$1" == rm && "$2" == -f && $# == 3 ]] || return 95
            [[ "$3" == fixture-owned || "$3" == fixture-owned-cli ]] || return 95
            [[ "$fail_cleanup" == false ]]
        }
        rm() {
            [[ $# == 3 && "$1" == -rf && "$2" == -- ]] || return 95
            [[ "$3" == /tmp/lcm-fixture-scratch || "$3" == /tmp/lcm-fixture-certificates ]]
        }
        """
        for (original, fail, expected) in ((0, false, 0), (7, false, 7),
                (0, true, 1), (7, true, 1))
            script = setup * "\nfail_cleanup=$fail\n" * cleanup.match *
                "\ntrap cleanup EXIT\nexit $original\n"
            process = run(pipeline(ignorestatus(`/bin/bash --noprofile --norc -c $script`);
                stdout=devnull, stderr=devnull))
            @test process.exitcode == expected
        end
    end
end
