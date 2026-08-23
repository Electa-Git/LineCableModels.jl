@testitem "Core / architecture / removed infrastructure stays removed" tags=[:unit] begin
    root=pkgdir(LineCableModels)
    source_root=joinpath(root, "src")
    source_files=filter(endswith(".jl"),
        [joinpath(directory, file)
         for (directory, _, files) in walkdir(source_root)
         for file in files])
    source=Dict(relpath(path, root)=>read(path, String) for path in source_files)

    @test !isdefined(LineCableModels, :Commons)
    @test !isdefined(LineCableModels, :Utils)
    @test !ispath(joinpath(source_root, "commons"))
    @test !ispath(joinpath(source_root, "utils"))
    @test isempty(filter(path -> basename(path) == "helpers.jl", source_files))

    forbidden=(
        "resolve_T", "coerce_to_T", "BASE_FLOAT", "REALSCALAR", "COMPLEXSCALAR",
        "@construct", "@parameterize", "@measurify", "Meta.parse", "global_logger("
    )
    for token in forbidden
        @test all(!occursin(token, contents) for contents in values(source))
    end
    @test all(!occursin("@assert", contents) for contents in values(source))

    macro_path=joinpath("src", "grammar", "macros.jl")
    macro_facades=Set((
        joinpath("src", "Grammar.jl"),
        joinpath("src", "LineCableModels.jl"),
        joinpath("src", "parametricbuilder", "ParametricBuilder.jl")
    ))
    for (path, contents) in source
        path == macro_path && continue
        if path in macro_facades
            @test occursin("@relax", contents)
            @test !occursin(r"\brecast\b", contents)
            continue
        end
        @test !occursin("@relax", contents)
        @test !occursin(r"\brecast\b", contents)
    end

    maintained=filter(pair -> pair.first != "src/retired.jl", collect(source))
    @test all(!occursin("calc_", contents) for (_, contents) in maintained)

    eager_files=(
        "circstrands.jl", "rectstrands.jl", "strip.jl", "tubular.jl",
        "insulator.jl", "semicon.jl", "conductorgroup.jl", "insulatorgroup.jl",
        "cablecomponent.jl", "cabledesign.jl"
    )
    for file in eager_files
        contents=source[joinpath("src", "datamodel", file)]
        @test !occursin(r"\btemperature\s*::", contents)
        @test !occursin(r";\s*temperature\s*=", contents)
    end
end

@testitem "ParametricBuilder / relax / isolated future provision" tags=[:unit] begin
    using LineCableModels.ParametricBuilder: @relax

    @relax struct FutureValue{T <: Real}
        scalar::T
        values::Vector{T}
        label::Symbol
    end

    value=FutureValue(1.0f0, [2.0, 3.0], :future)
    @test value isa FutureValue{Float64}
    @test value.scalar == 1.0
    @test value.values == [2.0, 3.0]
    converted=convert(FutureValue{Float32}, value)
    @test converted isa FutureValue{Float32}
    @test converted.values == Float32[2, 3]
    @test converted.label === :future
end
