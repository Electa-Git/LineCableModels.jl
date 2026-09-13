@testitem "Quality / documentation / formula source ownership" tags = [:quality] begin
    using Base.Docs: DocStr, meta
    using DocStringExtensions: TypedMethodSignatures
    using LineCableModels

    root = pkgdir(LineCableModels)

    categories = (
        (
            module_owner = LineCableModels.Materials.TemperatureDependent,
            registry = LineCableModels.Materials.TemperatureDependent.formulas(),
            path = ("materials", "temperaturedependent", "formulas"),
            default = :default
        ),
        (
            module_owner = LineCableModels.Engine.PipeImpedance,
            registry = LineCableModels.Engine.PipeImpedance.formulas(),
            path = ("engine", "pipeimpedance", "formulas"),
            default = :default
        ),
        (
            module_owner = LineCableModels.Engine.InternalImpedance,
            registry = LineCableModels.Engine.InternalImpedance.formulas(),
            path = ("engine", "internalimpedance", "formulas"),
            default = :default
        ),
        (
            module_owner = LineCableModels.Engine.InsulationImpedance,
            registry = LineCableModels.Engine.InsulationImpedance.formulas(),
            path = ("engine", "insulationimpedance", "formulas"),
            default = :default
        ),
        (
            module_owner = LineCableModels.Engine.InsulationAdmittance,
            registry = LineCableModels.Engine.InsulationAdmittance.formulas(),
            path = ("engine", "insulationadmittance", "formulas"),
            default = :default
        ),
        (
            module_owner = LineCableModels.Engine.SemiconAdmittance,
            registry = LineCableModels.Engine.SemiconAdmittance.formulas(),
            path = ("engine", "semiconadmittance", "formulas"),
            default = :default
        ),
        (
            module_owner = LineCableModels.Engine.EarthImpedance,
            registry = LineCableModels.Engine.EarthImpedance.formulas(),
            path = ("engine", "earthimpedance", "formulas"),
            default = :default
        ),
        (
            module_owner = LineCableModels.Engine.EarthAdmittance,
            registry = LineCableModels.Engine.EarthAdmittance.formulas(),
            path = ("engine", "earthadmittance", "formulas"),
            default = :default
        ),
        (
            module_owner = LineCableModels.Transforms,
            registry = LineCableModels.Transforms.formulas(),
            path = ("transforms", "formulas"),
            default = :default
        ),
        (
            module_owner = LineCableModels.Earth.FrequencyDependent,
            registry = LineCableModels.Earth.FrequencyDependent.formulas(),
            path = ("earth", "frequencydependent", "formulas"),
            default = :default
        ),
        (
            module_owner = LineCableModels.Earth.EquivalentHomogeneous,
            registry = LineCableModels.Earth.EquivalentHomogeneous.formulas(),
            path = ("earth", "equivalenthomogeneous", "formulas"),
            default = :default
        )
    )

    for category in categories
        directory = joinpath(root, "src", category.path...)
        formula_files = sort(filter(
            file -> endswith(file, ".jl"),
            readdir(directory)
        ))

        selected_default = category.module_owner.Formula(:default)
        @test formula_id(selected_default) === :default
        @test :default in category.registry

        for file in formula_files
            # Static FrequencyDependent is an intentional identity relation, exercised numerically.
            category.module_owner === LineCableModels.Earth.FrequencyDependent &&
                file == "default.jl" && continue
            source = read(joinpath(directory, file), String)
            @test isnothing(match(
                r"(?m)^\s*return\s+(material|model|parameters|input|source)\s*$",
                source
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
            descriptions
        ))
        @test documented_files == formula_files

        # Every registered choice must reach a real description. Whether its
        # sole authority dispatches on an instance or its type is not a docs
        # invariant; indexing the method-signature representation rejected
        # legitimate type/instance delegation without checking this behavior.
        for identifier in category.registry
            selected=category.module_owner.Formula(identifier)
            @test formula_id(selected)===identifier
            @test description(selected) isa AbstractString
            @test !isempty(description(selected))
            for compact in (false,true)
                @test description(selected;compact)==description(typeof(selected);compact)
                @test description(selected;compact) isa AbstractString
                @test !isempty(description(selected;compact))
            end
        end
        @test !applicable(description,category.module_owner.Formula{:UnregisteredDescriptionTest})

        for (_, docstring) in descriptions
            @test any(
                value -> value isa TypedMethodSignatures,
                collect(docstring.text)
            )
            strings = filter(value -> value isa String, collect(docstring.text))
            scientific_text = join(strings)
            @test occursin("**Identification.**", scientific_text)
            @test occursin("**Expression.**", scientific_text)
            if category.module_owner === LineCableModels.Earth.FrequencyDependent &&
               basename(String(docstring.data[:path])) == "default.jl"
                @test !occursin("**Reference.**", scientific_text)
            else
                @test occursin("**Reference.**", scientific_text)
            end
            @test isnothing(match(r"(?m)^\d{4}\.", scientific_text))
        end
    end
end
