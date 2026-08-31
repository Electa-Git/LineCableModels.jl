@testitem "Quality / documentation / formulation catalogue ownership" tags = [:quality] begin
    using Base.Docs: DocStr, meta
    using DocStringExtensions: TypedMethodSignatures
    using LineCableModels

    root = pkgdir(LineCableModels)
    page = read(
        joinpath(root, "docs", "src", "transmission-line-parameters.md"),
        String,
    )
    navigation = read(joinpath(root, "docs", "make.jl"), String)

    categories = (
        (
            module_owner = LineCableModels.Engine.InternalImpedance,
            registry = LineCableModels.Engine.InternalImpedance.formulas(),
            path = ("engine", "internalimpedance", "formulas"),
            default = :Schelkunoff1934,
        ),
        (
            module_owner = LineCableModels.Engine.InsulationImpedance,
            registry = LineCableModels.Engine.InsulationImpedance.formulas(),
            path = ("engine", "insulationimpedance", "formulas"),
            default = :Ametani1980,
        ),
        (
            module_owner = LineCableModels.Engine.InsulationAdmittance,
            registry = LineCableModels.Engine.InsulationAdmittance.formulas(),
            path = ("engine", "insulationadmittance", "formulas"),
            default = :Ametani2004,
        ),
        (
            module_owner = LineCableModels.Engine.SemiconAdmittance,
            registry = LineCableModels.Engine.SemiconAdmittance.formulas(),
            path = ("engine", "semiconadmittance", "formulas"),
            default = :Ametani2004,
        ),
        (
            module_owner = LineCableModels.Engine.EarthImpedance,
            registry = LineCableModels.Engine.EarthImpedance.formulas(),
            path = ("engine", "earthimpedance", "formulas"),
            default = :Papadopoulos2010,
        ),
        (
            module_owner = LineCableModels.Engine.EarthAdmittance,
            registry = LineCableModels.Engine.EarthAdmittance.formulas(),
            path = ("engine", "earthadmittance", "formulas"),
            default = :Papadopoulos2010,
        ),
        (
            module_owner = LineCableModels.Transforms,
            registry = LineCableModels.Transforms.formulas(),
            path = ("transforms", "formulas"),
            default = :Chrysochos2014,
        ),
        (
            module_owner = LineCableModels.EarthProps.FD,
            registry = LineCableModels.EarthProps.FD.formulas(),
            path = ("earthprops", "fd", "formulas"),
            default = nothing,
        ),
        (
            module_owner = LineCableModels.EarthProps.EHEM,
            registry = LineCableModels.EarthProps.EHEM.formulas(),
            path = ("earthprops", "ehem", "formulas"),
            default = :Layer,
        ),
    )

    for category in categories
        directory = joinpath(root, "src", category.path...)
        formula_files = sort(filter(
            file -> endswith(file, ".jl"),
            readdir(directory),
        ))

        @test category.module_owner.DEFAULT === category.default
        @test :default ∉ category.registry
        selected_default = category.module_owner.Formula(:default)
        if category.default === nothing
            @test selected_default === nothing
        elseif applicable(formula_id, selected_default)
            @test formula_id(selected_default) === category.default
        else
            @test category.module_owner === LineCableModels.EarthProps.EHEM
            @test selected_default == LineCableModels.EarthProps.EHEM.Layer(-1)
        end

        for file in formula_files
            source = read(joinpath(directory, file), String)
            @test isnothing(match(
                r"(?m)^\s*return\s+(material|model|parameters|input|source)\s*$",
                source,
            ))
        end

        descriptions = Tuple{Any, DocStr}[]
        for (binding, multidoc) in meta(category.module_owner)
            binding.var === :description || continue
            for (typesig, docstring) in multidoc.docs
                dirname(String(docstring.data[:path])) == directory || continue
                push!(descriptions, (typesig, docstring))
            end
        end

        documented_files = sort(map(
            entry -> basename(String(last(entry).data[:path])),
            descriptions,
        ))
        @test documented_files == formula_files

        documented_identifiers = map(descriptions) do (typesig, _)
            formula_type = only(typesig.parameters)
            first(Base.unwrap_unionall(formula_type).parameters)
        end
        @test isempty(setdiff(category.registry, documented_identifiers))

        for (_, docstring) in descriptions
            @test any(
                value -> value isa TypedMethodSignatures,
                collect(docstring.text),
            )
            strings = filter(value -> value isa String, collect(docstring.text))
            scientific_text = join(strings)
            @test occursin("**Identification.**", scientific_text)
            @test occursin("**Expression.**", scientific_text)
            @test occursin("**Reference.**", scientific_text)
            @test isnothing(match(r"(?m)^\d{4}\.", scientific_text))
        end

        module_name = string(category.module_owner)
        @test occursin("Main.$module_name", page)
    end

    @test count("```@eval", page) == length(categories)
    @test !occursin("```@autodocs", page)
    @test !occursin(r"(?m)^### `:", page)
    @test !occursin("**Identification.**", page)
    @test !occursin("**Expression.**", page)
    @test !occursin("**Reference.**", page)
    @test occursin(
        "\"Transmission line parameters\" => \"transmission-line-parameters.md\"",
        navigation,
    )
end
