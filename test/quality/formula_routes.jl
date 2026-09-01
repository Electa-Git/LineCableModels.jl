@testitem "Quality / formula catalogues / bound route ownership" tags=[:quality] begin
    using LineCableModels

    const EN=LineCableModels.Engine
    const EP=LineCableModels.EarthProps
    const FM=LineCableModels.FormulaMethod

    catalogues=(
        (
            owner = EN.InternalImpedance,
            path = ("engine", "internalimpedance", "formulas"),
            template = :internal_impedance,
            methods = (
                :internal_impedance,
                :surface_impedance_state,
                :surface_impedances
            ),
            override = :inner
        ),
        (
            owner = EN.EarthImpedance,
            path = ("engine", "earthimpedance", "formulas"),
            templates = (:earth_impedance, :propagation_constant),
            methods = (
                :earth_impedance,
                :propagation_constant,
                :struve_h1,
                :closed_form_term,
                :downward_coefficients!,
                :upward_coefficients,
                :spectral_kernel,
                :conductor_order,
                :local_layer_depth
            ),
            unregistered = (:DeriSemlyen1981,),
            override = :self
        ),
        (
            owner = EN.EarthAdmittance,
            path = ("engine", "earthadmittance", "formulas"),
            templates = (
                :earth_potential_coefficient,
                :earth_impedance,
                :propagation_constant
            ),
            methods = (
                :earth_potential_coefficient,
                :earth_impedance,
                :propagation_constant,
                :integral_terms
            ),
            description_exceptions = (:IdealGround,),
            override = :self
        ),
        (
            owner = EN.InsulationImpedance,
            path = ("engine", "insulationimpedance", "formulas"),
            template = :insulation_impedance,
            methods = (:insulation_impedance,),
            override = :route
        ),
        (
            owner = EN.InsulationAdmittance,
            path = ("engine", "insulationadmittance", "formulas"),
            template = :insulation_material,
            methods = (:insulation_material,),
            override = :route
        ),
        (
            owner = EN.SemiconAdmittance,
            path = ("engine", "semiconadmittance", "formulas"),
            template = :semicon_material,
            methods = (:semicon_material,),
            override = :route
        ),
        (
            owner = EP.FD,
            path = ("earthprops", "fd", "formulas"),
            template = :earth_material,
            methods = (:earth_material,),
            override = :route
        ),
        (
            owner = EP.EHEM,
            path = ("earthprops", "ehem", "formulas"),
            template = :equivalent_material,
            methods = (
                :equivalent_material,
                :transverse_propagation_constant,
                :equivalent_conductivity_step
            ),
            override = :route
        ),
        (
            owner = LineCableModels.Transforms,
            path = ("transforms", "formulas"),
            template = :modal_operators,
            methods = (
                :modal_operators,
                :modal_basis,
                :newton_workspace,
                :newton_step!,
                :levenberg_marquardt_workspace,
                :levenberg_marquardt_step!,
                :levenberg_marquardt_residual!,
                :levenberg_marquardt_jacobian!
            ),
            override = :route
        )
    )

    binder=(identity, selector, value)->(identity, selector, value)
    bound=FM(Val(:BoundFormula), binder, Val(:interaction))
    @test bound(:runtime) == (
        Val(:BoundFormula), Val(:interaction), :runtime
    )
    @test_throws ArgumentError FM(Val(:BoundFormula), binder, :interaction)

    root=pkgdir(LineCableModels)
    common_methods=(:assumptions, :description, :routes)
    named_method=r"(?ms)^\s*(?:@inline\s+)?function\s+([A-Za-z_][A-Za-z0-9_!]*)\s*\(\s*(.*?)\)"

    for catalogue in catalogues
        owner=catalogue.owner
        identifiers=owner.formulas()
        directory=joinpath(root, "src", catalogue.path...)
        files=sort(filter(endswith(".jl"), readdir(directory)))
        unregistered=haskey(catalogue, :unregistered) ?
                     catalogue.unregistered : ()
        description_exceptions=haskey(catalogue, :description_exceptions) ?
                               catalogue.description_exceptions : ()
        templates=haskey(catalogue, :templates) ?
                  catalogue.templates : (catalogue.template,)
        sources=Dict{Symbol, Tuple{String, String}}()
        for file in files
            source=read(joinpath(directory, file), String)
            matches=collect(eachmatch(
                r"(?m)^:([A-Za-z][A-Za-z0-9]*)\s*$",
                source
            ))
            @test length(matches) == 1
            length(matches) == 1 || continue
            discovered=Symbol(only(matches).captures[1])
            sources[discovered]=(file, source)
        end
        @test Set(keys(sources)) == union(Set(identifiers), Set(unregistered))

        for identifier in identifiers
            formula=owner.Formula(identifier)
            @test formula_id(formula) === identifier
            @test !isempty(description(formula))
            if !(identifier in description_exceptions)
                @test occursin(
                    r"\(\d{4}(?:/\d{4})?\)$",
                    description(formula)
                )
            end

            selected=isdefined(owner, :routes) ?
                     values(owner.routes(formula)) : (formula.route,)
            for route in selected
                @test route isa FM
                @test typeof(route).parameters[1] === identifier
                @test nameof(route.method) in templates
                @test all(argument->argument isa Val, route.arguments)
            end

            @test_throws ArgumentError owner.Formula(identifier; unknown = true)

            _, source=sources[identifier]
            for matched in eachmatch(named_method, source)
                method_name=Symbol(matched.captures[1])
                method_name in common_methods && continue
                @test method_name in catalogue.methods
                arguments=strip(matched.captures[2])
                identity=Regex(
                    "^(?:[A-Za-z_][A-Za-z0-9_]*)?::Val\\{:$identifier\\}"
                )
                @test occursin(identity, arguments)
                @test isnothing(match(r"\d{4}", String(method_name)))
            end
        end

        for identifier in unregistered
            formula=owner.Formula(identifier)
            @test formula_id(formula) === identifier
            @test !isempty(description(formula))
            @test occursin(
                r"\(\d{4}(?:/\d{4})?\)$",
                description(formula)
            )
            @test isempty(owner.routes(formula))
        end

        first_identifier=first(identifiers)
        custom=(arguments...)->arguments
        override=NamedTuple{(catalogue.override,)}((custom,))
        overridden=owner.Formula(first_identifier; override...)
        @test formula_id(overridden) === first_identifier
        if isdefined(owner, :routes)
            @test getproperty(
                owner.routes(overridden), catalogue.override
            ) === custom
        else
            @test overridden.route === custom
        end

        @test_throws ArgumentError owner.Formula(:UnregisteredFormula)
    end

    for owner in (
            EN.InsulationImpedance,
            EN.InsulationAdmittance,
            EN.SemiconAdmittance,
            EP.FD,
            EP.EHEM,
            LineCableModels.Transforms
    )
        @test !isdefined(owner, :Functor)
    end

    insulation=EN.InsulationAdmittance.Formula(:Ametani2004)
    semicon=EN.SemiconAdmittance.Formula(:Ametani2004)
    @test formula_id(insulation) === formula_id(semicon) === :Ametani2004
    @test typeof(insulation) !== typeof(semicon)
end
