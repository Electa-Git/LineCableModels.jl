@testitem "Gmsh FEM / GetDP installs and executes from an empty artifact cache" tags=[:extension,:fem_numerical] begin
    # Resolve the real package artifact in a child process whose artifact depot
    # starts empty. Existing Julia packages load before the depot is replaced.
    mktempdir() do directory
        cache=joinpath(directory,"depot")
        mkpath(cache)
        script=joinpath(directory,"cold_getdp.jl")
        write(script,raw"""
            using LineCableModels, Gmsh
            extension=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
            cache=only(ARGS)
            @assert isempty(readdir(cache))
            empty!(DEPOT_PATH)
            push!(DEPOT_PATH,cache)
            selected=withenv("LINECABLEMODELS_GETDP"=>nothing) do
                extension._getdp_selection(Formulation(:LineCableModelsFEM))
            end
            @assert selected.source === :artifact
            @assert startswith(selected.path,joinpath(cache,"artifacts"))
            @assert isfile(selected.path)
            identity=extension._getdp_identity(selected.path)
            @assert occursin(r"Version\s*:\s*3\.5\.0",identity.info)
            @assert occursin("complex arithmetic",identity.info)
            println("COLD_GETDP_OK ",selected.artifact_hash)
            """)
        command=`$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) --startup-file=no --compiled-modules=existing $script $cache`
        output=read(command,String)
        @test occursin("COLD_GETDP_OK",output)
        @test isdir(joinpath(cache,"artifacts"))
    end
end
