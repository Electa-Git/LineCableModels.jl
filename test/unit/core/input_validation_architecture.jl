@testitem "Core / InputValidation / one direct check per owned input" tags=[:unit] begin
    using RequiredInterfaces

    const IV=LineCableModels.InputValidation
    const Grammar=LineCableModels.Grammar
    const Engine=LineCableModels.Engine

    root=pkgdir(LineCableModels)
    source_root=joinpath(root, "src")
    @test isfile(joinpath(source_root, "inputvalidation", "InputValidation.jl"))
    @test !isfile(joinpath(source_root, "validation", "Validation.jl"))
    @test parentmodule(validate) === IV
    @test !applicable(validate, (answer = 42,))

    source_files=[joinpath(directory, file)
                  for (directory, _, files) in walkdir(source_root)
                  for file in files if endswith(file, ".jl")]
    source=Dict(relpath(path, root)=>read(path, String) for path in source_files)
    for forbidden in (
        "OwnerRule", "Validation.rules", "rules(::Type",
        "_check_material", "_check_earth_material", "_check_earth_layer",
        "_check_earth_model", "_check_line_parameters_problem",
        "_check_cable_constants_problem", "_check_blueprint",
        "_check_workspace"
    )
        @test all(!occursin(forbidden, contents) for contents in values(source))
    end

    expected=Dict(
        Material=>joinpath("src", "materials", "material.jl"),
        LineCableModels.EarthProps.EarthMaterial =>
            joinpath("src", "earthprops", "earthmaterial.jl"),
        EarthLayer=>joinpath("src", "earthprops", "earthlayer.jl"),
        EarthModel=>joinpath("src", "earthprops", "earthmodel.jl"),
        LineCableModels.DataModel.PreviewShape =>
            joinpath("src", "datamodel", "preview", "geometry.jl"),
        CableDesign=>joinpath("src", "datamodel", "design", "cabledesign.jl"),
        LineCableSystem =>
            joinpath("src", "datamodel", "linecablesystem", "linecablesystem.jl"),
        LineParametersProblem=>joinpath("src", "engine", "problems.jl"),
        CableConstantsProblem=>joinpath("src", "engine", "cableconstants.jl"),
        LineParameters => joinpath("src", "engine", "lineparameters", "lineparameters.jl"),
        Engine.CableBlueprint=>joinpath("src", "engine", "blueprint.jl"),
        Engine.LineParametersWorkspace=>joinpath("src", "engine", "input.jl"),
        ModalTransformationProblem=>joinpath("src", "transforms", "problems.jl"),
        ParametricProblem=>joinpath("src", "parametricbuilder", "results.jl")
    )

    function referenced_functions(method::Method)
        found=Set{Function}()
        function visit(value)
            if value isa GlobalRef && isdefined(value.mod, value.name)
                resolved=getproperty(value.mod, value.name)
                resolved isa Function && push!(found, resolved)
            elseif value isa Expr
                foreach(visit, value.args)
            end
            return nothing
        end
        foreach(visit, Base.uncompressed_ast(method).code)
        return found
    end

    function package_owned(value::Union{Function, Module})
        owner=value isa Module ? value : parentmodule(value)
        while owner !== Main && owner !== Base && owner !== Core
            owner === LineCableModels && return true
            owner=parentmodule(owner)
        end
        return false
    end

    admitted_calls=Set{Function}((
        validate,
        LineCableModels.DataModel.area,
        LineCableModels.DataModel.outer_radius,
        nphases
    ))
    for (type, path) in expected
        method=which(validate, Tuple{type})
        @test normpath(abspath(String(method.file))) == normpath(joinpath(root, path))
    end
    # Inspect the dispatch table, not just the historical list above: a newly
    # added input validator must obey the same no-delegation/no-mutation rule.
    for method in methods(validate)
        package_owned(method.module) || continue
        calls=referenced_functions(method)
        private_calls=filter(calls) do call
            package_owned(call) && startswith(String(nameof(call)), "_")
        end
        @test isempty(private_calls)
        unowned_calls=filter(calls) do call
            package_owned(call) &&
                !startswith(String(nameof(call)), "#") &&
                call ∉ admitted_calls
        end
        @test isempty(unowned_calls)
        mutating_calls=filter(calls) do call
            name=String(nameof(call))
            endswith(name, "!") && name != "!"
        end
        @test isempty(mutating_calls)
    end

    function package_type(type)
        owner=parentmodule(Base.unwrap_unionall(type))
        while owner !== Main && owner !== Base && owner !== Core
            owner === LineCableModels && return true
            owner=parentmodule(owner)
        end
        return false
    end

    for abstract_owner in (
        AbstractMaterial,
        AbstractEarthModel,
        AbstractProblemDefinition
    )
        @test RequiredInterfaces.isInterface(abstract_owner)
        @test Tuple(RequiredInterfaces.functions(
            RequiredInterfaces.getInterface(abstract_owner)
        )) == (validate,)
        implementors=filter(
            package_type,
            RequiredInterfaces.nonabstract_subtypes(abstract_owner)
        )
        for implementor in implementors
            @test which(validate, Tuple{implementor}).sig !=
                  which(validate, Tuple{abstract_owner}).sig
        end
    end

    constructor_files=(
        joinpath("src", "materials", "material.jl"),
        joinpath("src", "earthprops", "earthmaterial.jl"),
        joinpath("src", "earthprops", "earthlayer.jl"),
        joinpath("src", "earthprops", "earthmodel.jl"),
        joinpath("src", "datamodel", "preview", "geometry.jl"),
        joinpath("src", "datamodel", "design", "cabledesign.jl"),
        joinpath("src", "datamodel", "linecablesystem", "linecablesystem.jl"),
        joinpath("src", "engine", "problems.jl"),
        joinpath("src", "engine", "cableconstants.jl"),
        joinpath("src", "engine", "blueprint.jl"),
        joinpath("src", "engine", "input.jl"),
        joinpath("src", "engine", "lineparameters", "lineparameters.jl"),
        joinpath("src", "transforms", "problems.jl"),
        joinpath("src", "parametricbuilder", "results.jl")
    )
    for path in constructor_files
        @test occursin(r"validate\s*\(\s*new", source[path])
    end

    gridspace=source[joinpath("src", "gridspace.jl")]
    @test occursin("Target <: AbstractProblemDefinition", gridspace)
    @test occursin("validate(point.build", gridspace)
    lineparameters=source[joinpath("src", "engine", "lineparameters.jl")]
    cableconstants=source[joinpath("src", "engine", "cableconstants.jl")]
    transforms=source[joinpath("src", "transforms", "compute.jl")]
    @test occursin("validate(problem)", lineparameters)
    @test occursin("validate(problem)", cableconstants)
    @test occursin("validate(problem)", transforms)
end
