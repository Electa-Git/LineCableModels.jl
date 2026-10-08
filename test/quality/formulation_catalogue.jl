@testitem "Quality / documentation / formula source ownership" tags = [:quality] setup = [FormulaFamilies] begin
    using Base.Docs: Binding, DocStr, meta
    using DocStringExtensions: TypedMethodSignatures
    using LineCableModels

    # Every registered family whose formulas each live in a file of the `formulas`
    # directory beside the module that registers them.
    categories = map(FormulaFamilies.families()) do family
        registration = which(LineCableModels.Commons.formulas, Tuple{Type{family.Formula}})
        (module_owner = family,
            registry = LineCableModels.Commons.formulas(family.Formula),
            directory = joinpath(dirname(String(registration.file)), "formulas"))
    end
    filter!(category -> isdir(category.directory), categories)
    @test length(categories) == 11

    for category in categories
        directory = category.directory
        formula_files = sort(filter(
            file -> endswith(file, ".jl"),
            readdir(directory)
        ))

        selected_default = category.module_owner.Formula(:default)
        @test formula_id(selected_default) !== :default
        @test :default in category.registry

        formula_docs = DocStr[]
        for (binding, multidoc) in meta(category.module_owner)
            binding.var in (:description, :earth_impedance, :earth_potential_coefficient,
                :source_coefficients) || continue
            for docstring in values(multidoc.docs)
                dirname(String(docstring.data[:path])) == directory || continue
                push!(formula_docs, docstring)
            end
        end

        documented_files = sort(unique(map(
            docstring -> basename(String(docstring.data[:path])),
            formula_docs
        )))
        @test documented_files == formula_files

        # Each registered choice must produce a real description. Whether its
        # sole authority dispatches on an instance or its type is not a docs
        # invariant. Indexing the method-signature representation rejected
        # legitimate type or instance delegation without checking this behavior.
        for identifier in category.registry
            selected=category.module_owner.Formula(identifier)
            @test formula_id(selected)===(identifier === :default ?
                                          formula_id(selected_default) : identifier)
            @test description(selected) isa AbstractString
            @test !isempty(description(selected))
            for compact in (false, true)
                @test description(selected; compact)==description(typeof(selected); compact)
                @test description(selected; compact) isa AbstractString
                @test !isempty(description(selected; compact))
            end
        end
        @test !applicable(description, category.module_owner.Formula{:UnregisteredDescriptionTest})

        for docstring in formula_docs
            @test any(
                value -> value isa TypedMethodSignatures,
                collect(docstring.text)
            )
        end
    end

    # Scientific documentation belongs to the implemented equation method.
    # Description methods retain only their compact and full labels.
    for (owner, name, identifier, kind, first_layer, second_layer) in (
        (LineCableModels.Engine.EarthImpedance, :earth_impedance, :saad1996, :self, 2, 2),
        (LineCableModels.Engine.EarthImpedance,
            :earth_impedance, :wedepohl1973, :self, 2, 2),
        (LineCableModels.Engine.EarthImpedance, :earth_impedance, :gary1976, :self, 1, 1),
        (LineCableModels.Engine.EarthImpedance,
            :earth_impedance, :lucca1994, :mutual, 1, 2),
        (LineCableModels.Engine.EarthAdmittance,
            :earth_potential_coefficient, :ideal, :self, 1, 1))
        signature = Tuple{owner.Formula{identifier}, Val{kind},
            Val{first_layer}, Val{second_layer}, Any, Any}
        docstring = meta(owner)[Binding(owner, name)].docs[signature]
        @test basename(String(docstring.data[:path])) == string(identifier, ".jl")
        @test !any(
            doc -> basename(String(doc.data[:path])) == string(identifier, ".jl"),
            values(meta(owner)[Binding(owner, :description)].docs)
        )
    end
end

@testitem "Quality / local shunt formulations" tags = [:quality] begin
    const owner = LineCableModels.Engine.ShuntModel
    @test LineCableModels.Commons.formulas(owner.Formula) == (:default, :equivalent, :boundary)
    for identifier in LineCableModels.Commons.formulas(owner.Formula)
        selected = owner.Formula(identifier)
        @test formula_id(selected) === (identifier === :default ? :equivalent : identifier)
        @test NamedTuple(selected).identifier === formula_id(selected)
        for compact in (false, true)
            @test !isempty(description(selected; compact))
            @test description(selected; compact) == description(typeof(selected); compact)
        end
    end
end
