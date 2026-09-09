using Test, TOML, UUIDs, LineCableModelsRuntime
const ImageLayoutRT = LineCableModelsRuntime

@testset "image commands and relocated Julia sources follow the shared policy" begin
    root = normpath(joinpath(@__DIR__,"..","..",".."))
    recipe = read(joinpath(root,"playground","worker","containers","Containerfile"),String)
    copies = Dict{String,Vector{Vector{String}}}()
    commands = Dict{String,Vector{String}}()
    entrypoint = String[]
    stage = ""
    for line in split(recipe,'\n')
        if startswith(line,"FROM ")
            stage = last(split(line)); copies[stage] = Vector{String}[]
        elseif startswith(line,"COPY ")
            push!(copies[stage],String.(split(line)[2:end]))
        elseif startswith(line,"CMD ")
            commands[stage] = String.(ImageLayoutRT.JSON3.read(line[5:end]))
        elseif startswith(line,"ENTRYPOINT ")
            entrypoint = String.(ImageLayoutRT.JSON3.read(line[12:end]))
        end
    end
    @test Set(keys(copies)) == Set(("shared","line-parameters","power-flow","julia-terminal"))
    @test !occursin("COPY .",recipe)
    @test !occursin("ln -s",recipe)
    @test entrypoint == ["/usr/local/julia/bin/julia"]
    for name in ("line-parameters","power-flow","julia-terminal")
        kind = name == "julia-terminal" ? :terminal : :scientific
        profile = ProfileDefinition(name,"registry.invalid/lcm@sha256:" * repeat("a",64),repeat("a",64);
            kind,isolation=:container,operations=kind == :scientific ? ("fixture.echo",) : ())
        fence = ImageLayoutRT.Protocol.AssignmentFence(string(uuid4()),string(uuid4()),"alice","main","worker-a",
            string(uuid4()),string(uuid4()),name,"1.0.0",profile.fingerprint,1)
        policy = ContainerPolicy(profile,ResourceReceipt(uuid4(),uuid4(),fence,:podman,repeat("b",64),nothing))
        @test commands[name] == ImageLayoutRT.container_julia_arguments(policy)
        @test "--entrypoint=" * only(entrypoint) in container_create_arguments(policy)
        docker_policy = ContainerPolicy(profile,ResourceReceipt(uuid4(),uuid4(),fence,:docker,repeat("b",64),nothing))
        @test "--entrypoint=" * only(entrypoint) in container_create_arguments(docker_policy)
        mktempdir() do directory
            for fields in [copies["shared"]; copies[name]]
                destination = joinpath(directory,lstrip(last(fields),'/'))
                for relative in fields[1:end-1]
                    source = joinpath(root,relative)
                    @test ispath(source)
                    target = !isdir(source) && endswith(last(fields),"/") ? joinpath(destination,basename(source)) : destination
                    mkpath(dirname(target))
                    cp(source,target)
                end
            end
            project = joinpath(directory,"opt","lcm","source","playground","worker","profiles","active")
            manifest = TOML.parsefile(joinpath(project,"Manifest.toml"))
            for entries in values(manifest["deps"]), entry in entries
                haskey(entry,"path") || continue
                source = normpath(joinpath(project,entry["path"]))
                @test startswith(source,directory * "/")
                @test isfile(joinpath(source,"Project.toml"))
            end
            @test !islink(project)
            julia = joinpath(Sys.BINDIR,Base.julia_exename())
            expression = "using LineCableModelsExecutionCore; @assert !any(m -> nameof(m) in (:Bonito,:NATS,:LineCableModels,:PowerImpedance),values(Base.loaded_modules)); println(\"relocated_core_ok\")"
            command = setenv(`$julia --startup-file=no --compiled-modules=existing --project=$project -e $expression`,
                ["PATH"=>ENV["PATH"],"JULIA_LOAD_PATH"=>"@:@stdlib","JULIA_DEPOT_PATH"=>join(DEPOT_PATH,':'),
                    "JULIA_NUM_THREADS"=>"1","OPENBLAS_NUM_THREADS"=>"1"])
            runner = CommandRunner(;timeout_seconds=30)
            try
                result = run_owned_command!(runner,command)
                @test result.exitcode == 0
                @test strip(result.output) == "relocated_core_ok"
            finally
                close(runner)
            end
        end
    end
end
