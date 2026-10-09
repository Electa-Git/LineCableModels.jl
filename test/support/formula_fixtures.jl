@testmodule FormulaFixtures begin
    using LineCableModels
    const E = LineCableModels.Engine
    const G = LineCableModels.Commons
    const II = E.InternalImpedance
    const EI = E.EarthImpedance
    const EA = E.EarthAdmittance
    const EP = LineCableModels.Earth
    const FD = EP.FrequencyDependent
    const TD = LineCableModels.Materials.TemperatureDependent
    const EH = EP.EquivalentHomogeneous
    const IA = E.InsulationAdmittance
    const SA = E.SemiconAdmittance
    const Expression = LineCableModels.Expression
    const calls = Tuple[]

    struct BufferReplacement <: E.AbstractFormulation end
    G.initialize_buffers(::BufferReplacement, ::Type, input, plan,
        buffers) = merge(buffers, (destination = copy(buffers.destination),))

    # These scientific selections belong to the consumer and leave built-in
    # Formula registrations and private constructors unchanged. The built-in formula lists remain closed.
    # These formulas have expressions up to layer 3 (`Layer…`) or layer 2 (`HalfSpace…`,
    # usable in a recipe). Their expression signatures alone declare the media they handle.
    for (name, parent, owner, operation, layers) in (
        (:LayerImpedance, E.EarthImpedanceFormulation, EI, EI.earth_impedance, 3),
        (:LayerPotential, E.EarthAdmittanceFormulation, EA, EA.earth_potential_coefficient, 3),
        (:HalfSpaceImpedance, E.EarthImpedanceFormulation, EI, EI.earth_impedance, 2),
        (:HalfSpacePotential, E.EarthAdmittanceFormulation, EA, EA.earth_potential_coefficient, 2))
        @eval struct $name{P, O, R} <: $parent
            parameters::P
            options::O
            equivalent_earth::R
            initialized::Vector{Tuple}
        end
        operation_name = GlobalRef(owner, nameof(operation))
        for s in 1:layers, t in 1:layers, kind in (s == t ? (:self, :mutual) : (:mutual,))
            @eval function $operation_name(
                    selected::$name, ::Val{$(QuoteNode(kind))},
                    ::Val{$s}, ::Val{$t}, functor, workspace)
                input = functor.input
                pair = input.pair
                push!(calls,
                    ($(QuoteNode(nameof(owner))), input.jω,
                        pair.row, pair.column, pair.layers, pair.heights,
                        copy(input.rho), input.physical))
                coefficient = 11 * $s + 17 * $t + 3 * pair.row + 5 * pair.column +
                              real(input.jω / (2pi*im)) / 100 +
                              (pair.row == pair.column ? 101 : 0)
                if haskey(input.options.data, :integration)
                    integral = E.SpectralIntegral(λ -> complex(exp(-2λ)))
                    value,
                    _ = E.integrate(integral,
                        input.options.data.integration.method, input.options.data.integration.options,
                        workspace.buffers)
                    coefficient *= 2value
                end
                return selected.parameters.scale *
                       $(owner === EI ? :(coefficient * (1e-4 + 1e-3im)) :
                         :(coefficient * 1e9))
            end
        end
        @eval LineCableModels.formulation_options(::Expression{
            <:$name, typeof($operation)}) = FormulationOptions()
    end
    LineCableModels.formulation_options(::Expression{<:Union{LayerImpedance, HalfSpaceImpedance},
        typeof(EI.earth_impedance),
        A}) where {A <:
                   Tuple{Union{Val{:self}, Val{:mutual}}, Val{1},
        Val{1}}} = FormulationOptions(integration = (method = :quad, options = (;)))

    function G.initialize_buffers(
            selected::Union{LayerImpedance, HalfSpaceImpedance}, ::Type{T}, input, plan,
            buffers) where {T}
        layers = [E.layer_index(entry.pair)
                  for calculation in plan.earth.calculations
                  if something(calculation.impedance, calculation.admittance).formula === selected
                  for entry in calculation.pairs]
        initialized = (1, 1) in layers ?
                      G.initialize_buffers(
                          E.SpectralIntegral, Val(:quad), T, input, plan, buffers) :
                      buffers
        push!(selected.initialized, (layers, get(initialized, :quadrature, nothing)))
        return initialized
    end

    function selection(
            owner; options::Union{NamedTuple, FormulationOptions} = FormulationOptions(),
            scale = 1.0, equivalent_earth = nothing, layers = 3)
        options = options isa NamedTuple ? FormulationOptions(options) : options
        selected_type = layers == 3 ? (owner === EI ? LayerImpedance : LayerPotential) :
                        (owner === EI ? HalfSpaceImpedance : HalfSpacePotential)
        return selected_type((scale = scale,), options, equivalent_earth, Tuple[])
    end

    # These formulas admit the deepest-layer reduction.
    LineCableModels.validate(rule::EH.Formula{:bottommost},
        ::Expression{<:Union{LayerImpedance, LayerPotential, HalfSpaceImpedance, HalfSpacePotential}}) = rule

    # Two more layered earth-impedance formulas. One declares its expressions up to layer 3
    # with typed runtime arguments, the other for every layer with generic `Val{S}`
    # positions. Each records the number of media its expressions receive.
    struct TypedLayerImpedance{P, O} <: E.EarthImpedanceFormulation
        parameters::P
        options::O
        equivalent_earth::Nothing
        media::Vector{Int}
    end
    TypedLayerImpedance() = TypedLayerImpedance((;), FormulationOptions(), nothing, Int[])
    for s in 1:3, t in 1:3, kind in (s == t ? (:self, :mutual) : (:mutual,))
        @eval function EI.earth_impedance(selected::TypedLayerImpedance,
                ::Val{$(QuoteNode(kind))}, ::Val{$s}, ::Val{$t},
                functor::G.Functor, workspace)
            push!(selected.media, length(functor.input.rho))
            pair = functor.input.pair
            return (1e-4 + 1e-3im) * (pair.row == pair.column ? 10 : 1)
        end
    end

    struct GenericLayerImpedance{P, O} <: E.EarthImpedanceFormulation
        parameters::P
        options::O
        equivalent_earth::Nothing
        media::Vector{Int}
    end
    GenericLayerImpedance() = GenericLayerImpedance((;), FormulationOptions(), nothing, Int[])
    function EI.earth_impedance(selected::GenericLayerImpedance,
            ::Union{Val{:self}, Val{:mutual}}, ::Val{S}, ::Val{T},
            functor, workspace) where {S, T}
        push!(selected.media, length(functor.input.rho))
        pair = functor.input.pair
        return (1e-4 + 1e-3im) * (pair.row == pair.column ? 10 : 1)
    end

    for selected_type in (TypedLayerImpedance, GenericLayerImpedance)
        @eval LineCableModels.formulation_options(::Expression{
            <:$selected_type, typeof(EI.earth_impedance)}) = FormulationOptions()
    end

    # These formulas admit nonzero permittivities of either sign.
    function LineCableModels.validate(rho::AbstractVector,
            ::Union{LayerImpedance, LayerPotential, HalfSpaceImpedance, HalfSpacePotential},
            epsilon::AbstractVector, mu::AbstractVector, thickness)
        length(rho) == length(epsilon) == length(mu) ||
            throw(DimensionMismatch("material vectors must align"))
        all(x -> isfinite(x) && !iszero(x), epsilon) || throw(DomainError(epsilon,
            "the layered formula requires nonzero finite permittivities"))
        return rho
    end

    # An algebraic coupled equation uses the same main workspace and explicit
    # calculation action without any formula-owned lifecycle.
    struct CoupledImpedance{P, O} <: E.EarthImpedanceFormulation
        parameters::P
        options::O
        equivalent_earth::Nothing
        solves::Vector{ComplexF64}
    end
    CoupledImpedance() = CoupledImpedance((;), FormulationOptions(), nothing, ComplexF64[])
    LineCableModels.formulation_options(::Expression{
        <:CoupledImpedance, typeof(EI.earth_impedance)}) = FormulationOptions()
    function G.initialize_buffers(
            ::CoupledImpedance, ::Type{T}, input, plan, buffers) where {T}
        geometry=plan.geometry
        return merge(buffers,
            (coupled = Matrix{Complex{T}}(undef,
                length(geometry.radius), length(geometry.radius)),))
    end
    # The parts write into the formula's own matrix, which its conversion completes.
    function G.Functor(selected::CoupledImpedance, input::NamedTuple; workspace)
        return G.Functor(selected, merge(input, (destinations = (workspace.buffers.coupled,),)), (;))
    end
    function E.earth!(selected::CoupledImpedance, functor::G.Functor, workspace)
        values=only(functor.input.destinations)
        # Complete-system coupling adds the sum over all source conductors.
        values .+= functor.input.jω*1e-6*sum(axes(values, 2))
        push!(selected.solves, functor.input.jω)
        return (impedance = values,)
    end
    for source in 1:2, target in 1:2,
        kind in (source==target ? (:self, :mutual) : (:mutual,))
        @eval function EI.earth_impedance(::CoupledImpedance, ::Val{$(QuoteNode(kind))},
                ::Val{$source}, ::Val{$target}, functor, workspace)
            n=length(workspace.plan.geometry.radius)
            pair=functor.input.pair
            return functor.input.jω*1e-6*(pair.row+2pair.column+(pair.row==pair.column ? n :
                                                                 0))
        end
    end

    struct SurfaceLaw{Kinds, P, O} <: E.InternalImpedanceFormulation
        parameters::P
        options::O
        state_inputs::Vector{Tuple}
        evaluations::Vector{Tuple}
    end
    function SurfaceLaw(; kinds = (:inner, :outer, :transfer),
            coefficients = (inner = 2.0+1im, outer = 2.0+1im, transfer = 0.5+0.1im))
        parameters = (coefficients = coefficients,)
        options = FormulationOptions(NamedTuple{kinds}(map(_ -> (;), kinds)))
        SurfaceLaw{kinds, typeof(parameters), typeof(options)}(
            parameters, options, Tuple[], Tuple[])
    end
    # The surface impedances share a serial number of the evaluation point.
    function G.Functor(selected::SurfaceLaw, input::NamedTuple; workspace = nothing)
        (; r_in, r_ex, rho, mu_r, jω) = input
        push!(selected.state_inputs, (r_in, r_ex, rho, mu_r, jω))
        state = (serial = length(selected.state_inputs), rho = rho, radius = r_ex)
        return G.Functor(selected, input, state)
    end
    # Availability is dispatched by actual surface, not inferred from options.
    for kinds in ((:inner, :outer, :transfer), (:outer,), (:transfer,), (:inner,)),
        kind in kinds

        @eval function II.internal_impedance(selected::SurfaceLaw{$kinds},
                ::Val{$(QuoteNode(kind))}, functor, workspace)
            push!(selected.evaluations, (functor.state.serial, Val($(QuoteNode(kind)))))
            return selected.parameters.coefficients[$(QuoteNode(kind))]
        end
    end
    LineCableModels.formulation_options(::Expression{
        <:SurfaceLaw, typeof(II.internal_impedance)}) = FormulationOptions()

    struct SpectralSurface{P, O} <: E.InternalImpedanceFormulation
        parameters::P
        options::O
        seen::Vector{Tuple}
    end
    function SpectralSurface(method = :quad)
        seed = SpectralSurface((;), FormulationOptions(outer = (;)), Tuple[])
        outer = formulation_options(Expression(seed, II.internal_impedance, Val(:outer)),
            FormulationOptions(integration = (method = method,)))
        return SpectralSurface((;), FormulationOptions(outer = outer.data), seed.seen)
    end
    function G.initialize_buffers(
            selected::SpectralSurface, ::Type{T}, input, plan, buffers) where {T}
        push!(selected.seen, (:initialize, buffers))
        return G.initialize_buffers(
            E.SpectralIntegral, Val(:quad), T, input, plan, buffers)
    end
    LineCableModels.formulation_options(::Expression{
        <:SpectralSurface, typeof(II.internal_impedance),
        Tuple{Val{:outer}}}) = FormulationOptions(integration = (
        method = :quad, options = (;)))
    # Its inner and transfer surface impedances are fixed coefficients. Only `outer` integrates.
    II.internal_impedance(::SpectralSurface, ::Val{:inner}, functor, workspace) = 3e-5 + 1e-6im
    II.internal_impedance(::SpectralSurface, ::Val{:transfer}, functor, workspace) = 1e-6 + 0im
    LineCableModels.formulation_options(::Expression{<:SpectralSurface, typeof(II.internal_impedance),
        <:Union{Tuple{Val{:inner}}, Tuple{Val{:transfer}}}}) = FormulationOptions()
    function II.internal_impedance(selected::SpectralSurface, ::Val{:outer}, functor, workspace)
        integration=functor.input.options.data.integration
        push!(selected.seen, (integration.method, workspace))
        integral=E.SpectralIntegral(λ->complex(exp(-2λ)))
        value,
        _=E.integrate(integral, integration.method, integration.options,
            workspace===nothing ? nothing : workspace.buffers)
        return value*1e-4
    end
    LineCableModels.description(::Type{<:SpectralSurface}, ::Val{:method},
        ::Val{:quad}; compact::Bool = false) = "quad"

    struct DispersiveEarth{P, O} <: FD.FrequencyDependentFormulation
        parameters::P
        options::O
        seen::Vector{Tuple}
        events::Vector{Symbol}
    end
    DispersiveEarth(; scale = 100.0,
        exponent = 1,
        events = Symbol[]) = DispersiveEarth(
        (scale = scale, exponent = exponent), FormulationOptions(), Tuple[], events)
    function FD.earth_material(selected::DispersiveEarth, functor, workspace)
        (; material, frequency) = functor.input
        parameters = selected.parameters
        push!(selected.seen, (material.rho, frequency, workspace))
        push!(selected.events, :fd)
        return EP.EarthMaterial(
            material.rho/(1+(frequency/parameters.scale)^parameters.exponent),
            material.eps_r, material.mu_r)
    end
    LineCableModels.formulation_options(::Expression{
        <:DispersiveEarth, typeof(FD.earth_material)}) = FormulationOptions()

    struct InsulationReactance{P, O} <: E.InsulationImpedanceFormulation
        parameters::P
        options::O
    end
    InsulationReactance(inductance = 2) = InsulationReactance((inductance = inductance,), FormulationOptions())
    E.InsulationImpedance.insulation_impedance(selected::InsulationReactance, functor,
        workspace) = selected.parameters.inductance*functor.input.jω
    LineCableModels.formulation_options(::Expression{<:InsulationReactance,
        typeof(E.InsulationImpedance.insulation_impedance)}) = FormulationOptions()

    struct ConstantResistivity{P, O} <: TD.TemperatureDependentFormulation
        parameters::P
        options::O
        seen::Vector{Tuple}
    end
    ConstantResistivity(value) = ConstantResistivity((rho = value,), FormulationOptions(), Tuple[])
    function TD.temperature_resistivity(selected::ConstantResistivity, functor, workspace)
        (; material, temperature) = functor.input
        push!(selected.seen, (material, temperature, workspace))
        return selected.parameters.rho
    end
    LineCableModels.formulation_options(::Expression{
        <:ConstantResistivity, typeof(TD.temperature_resistivity)}) = FormulationOptions()

    struct ScaledResistivity{P, O} <: TD.TemperatureDependentFormulation
        parameters::P
        options::O
        seen::Vector{Tuple}
    end
    ScaledResistivity(scale = 2) = ScaledResistivity((scale = scale,), FormulationOptions(), Tuple[])
    function TD.temperature_resistivity(selected::ScaledResistivity, functor, workspace)
        (; material, temperature) = functor.input
        push!(selected.seen, (material, temperature, workspace))
        return selected.parameters.scale*material.rho
    end
    LineCableModels.formulation_options(::Expression{
        <:ScaledResistivity, typeof(TD.temperature_resistivity)}) = FormulationOptions()

    struct ExponentialResistivity{P, O} <: TD.TemperatureDependentFormulation
        parameters::P
        options::O
    end
    ExponentialResistivity(scale = 1000.0) = ExponentialResistivity((scale = scale,), FormulationOptions())
    function TD.temperature_resistivity(selected::ExponentialResistivity, functor, workspace)
        (; material, temperature) = functor.input
        return material.rho*exp((temperature-material.T0)/selected.parameters.scale)
    end
    LineCableModels.formulation_options(::Expression{
        <:ExponentialResistivity, typeof(TD.temperature_resistivity)}) = FormulationOptions()

    struct DispersiveSoil{P, O} <: FD.FrequencyDependentFormulation
        parameters::P
        options::O
        seen::Vector{Float64}
    end
    DispersiveSoil() = DispersiveSoil(
        (rho_scale = 1000, epsilon_scale = 2000, mu_scale = 10000),
        FormulationOptions(), Float64[])
    function FD.earth_material(selected::DispersiveSoil, functor, workspace)
        (; material, frequency) = functor.input
        p = selected.parameters
        push!(selected.seen, frequency)
        EP.EarthMaterial(material.rho/(1+frequency/p.rho_scale),
            material.eps_r*(1+frequency/p.epsilon_scale), material.mu_r*(1+frequency/p.mu_scale))
    end
    LineCableModels.formulation_options(::Expression{
        <:DispersiveSoil, typeof(FD.earth_material)}) = FormulationOptions()

    struct ScaledSoil{P, O} <: FD.FrequencyDependentFormulation
        parameters::P
        options::O
    end
    ScaledSoil(; rho = 1, epsilon = 1,
        mu = 1) = ScaledSoil((rho = rho, epsilon = epsilon, mu = mu), FormulationOptions())
    function FD.earth_material(selected::ScaledSoil, functor, workspace)
        (; material) = functor.input
        p = selected.parameters
        return EP.EarthMaterial(material.rho*p.rho, material.eps_r*p.epsilon, material.mu_r*p.mu)
    end
    LineCableModels.formulation_options(::Expression{
        <:ScaledSoil, typeof(FD.earth_material)}) = FormulationOptions()

    struct OhmicDielectric{P, O} <: E.InsulationAdmittanceFormulation
        parameters::P
        options::O
        temperatures::Vector{Tuple}
    end
    OhmicDielectric() = OhmicDielectric((;), FormulationOptions(), Tuple[])
    function IA.insulation_material(selected::OhmicDielectric, functor, workspace)
        (; material, frequency, temperature) = functor.input
        push!(selected.temperatures, (material.T0, temperature))
        complex(inv(material.rho), 2pi*frequency*8.8541878128e-12*material.eps_r)
    end
    LineCableModels.formulation_options(::Expression{
        <:OhmicDielectric, typeof(IA.insulation_material)}) = FormulationOptions()

    for (name, parent, operation) in (
        (:InsulationLaw, E.InsulationAdmittanceFormulation, IA.insulation_material),
        (:SemiconLaw, E.SemiconAdmittanceFormulation, SA.semicon_material))
        @eval struct $name{P, O} <: $parent
            parameters::P
            options::O
        end
        @eval $name(; scale = 2) = $name((scale = scale,), FormulationOptions())
        operation_name=GlobalRef(parentmodule(operation), nameof(operation))
        @eval function $operation_name(selected::$name, functor, workspace)
            (; material, frequency, temperature) = functor.input
            return complex(oftype(frequency, selected.parameters.scale)/material.rho,
                frequency*material.eps_r+temperature)
        end
        @eval LineCableModels.formulation_options(::Expression{
            <:$name, typeof($operation)}) = FormulationOptions()
    end

    struct MeanEarth{P, O} <: EH.AbstractRule
        parameters::P
        options::O
        seen::Vector{Tuple}
    end
    MeanEarth() = MeanEarth((;), FormulationOptions(), Tuple[])
    function EH.equivalent_material(selected::MeanEarth, ::Val{Kind}, ::Val{S}, ::Val{T},
            functor, workspace) where {Kind, S, T}
        (; rho, eps_r, mu_r, pair, frequency) = functor.input
        push!(selected.seen, (copy(rho), pair, frequency, workspace))
        soils=2:length(rho)
        return EP.EarthMaterial(sum(rho[soils])/length(soils),
            sum(eps_r[soils])/length(soils), sum(mu_r[soils])/length(soils))
    end
    LineCableModels.formulation_options(::Expression{
        <:MeanEarth, typeof(EH.equivalent_material)}) = FormulationOptions()

    struct SquaredBottomEarth{P, O} <: EH.AbstractRule
        parameters::P
        options::O
        events::Vector{Symbol}
        workspaces::Vector{Any}
    end
    SquaredBottomEarth(events = Symbol[]) = SquaredBottomEarth(
        (scale = 100,), FormulationOptions(), events, Any[])
    function EH.equivalent_material(
            selected::SquaredBottomEarth, ::Val{Kind}, ::Val{S}, ::Val{T},
            functor, workspace) where {Kind, S, T}
        (; rho, eps_r, mu_r) = functor.input
        push!(selected.events, :ehem)
        push!(selected.workspaces, workspace)
        return EP.EarthMaterial(last(rho)^2/selected.parameters.scale, last(eps_r), last(mu_r))
    end
    LineCableModels.formulation_options(::Expression{
        <:SquaredBottomEarth, typeof(EH.equivalent_material)}) = FormulationOptions()
    LineCableModels.validate(reduction::SquaredBottomEarth,
        ::Expression{<:Union{EI.Formula{:unified}, EA.Formula{:unified}}}) = reduction

    struct UserCoaxialShunt{P, O} <: E.ShuntModelFormulation
        parameters::P
        options::O
        response_count::Base.RefValue{Int}
    end
    UserCoaxialShunt() = UserCoaxialShunt((;), FormulationOptions(), Ref(0))
    function E.internal_shunt_response(selected::UserCoaxialShunt, design::CableDesign,
            geometry, T, methods, solutions, design_index)
        selected.response_count[]+=1
        response=E.internal_shunt_response(E.ShuntModel.Formula(:equivalent),
            design, geometry, T, methods, solutions, design_index)
        merge(response, (details = merge(response.details, (requested = :UserCoaxialShunt,)),))
    end

    struct UserCoaxialPipe{P, O} <: E.PipeImpedanceFormulation
        parameters::P
        options::O
    end
    UserCoaxialPipe() = UserCoaxialPipe((;), FormulationOptions())
    # It admits coaxial topology, as the built-in formula without a pipe term does.
    LineCableModels.validate(design::CableDesign, ::UserCoaxialPipe, backend) =
        LineCableModels.validate(design, E.PipeImpedance.Formula(:none), backend)

    # Counter-bearing native types observe choreography without replacing any
    # package method or attaching executable values to a selection record.
    for (name, parent, owner, operation, id, event) in (
        (:CountedInsulationZ, E.InsulationImpedanceFormulation, E.InsulationImpedance,
            E.InsulationImpedance.insulation_impedance, :ametani1980, :local_z),
        (:CountedInsulationY, E.InsulationAdmittanceFormulation,
            IA, IA.insulation_material, :lossy, :local_y),
        (:CountedSemiconY, E.SemiconAdmittanceFormulation,
            SA, SA.semicon_material, :lossy, :local_y))
        @eval struct $name{F, P, O} <: $parent
            base::F
            parameters::P
            options::O
            events::Vector{Symbol}
            workspaces::Vector{Any}
        end
        @eval function $name(events)
            base=$owner.Formula($(QuoteNode(id)))
            $name(base, base.parameters, base.options, events, Any[])
        end
        op=GlobalRef(owner, nameof(operation))
        @eval function $op(selected::$name, functor, workspace)
            push!(selected.events, $(QuoteNode(event)))
            push!(selected.workspaces, workspace)
            $op(selected.base, functor, workspace)
        end
        @eval LineCableModels.formulation_options(::Expression{
            <:$name, typeof($operation)}) = FormulationOptions()
    end
    for (name, parent, owner, operation, event) in (
        (:CountedEarthZ, E.EarthImpedanceFormulation, EI, EI.earth_impedance, :earth_z),
        (:CountedEarthP, E.EarthAdmittanceFormulation,
            EA, EA.earth_potential_coefficient, :earth_y))
        @eval struct $name{P, O} <: $parent
            parameters::P
            options::O
            equivalent_earth::Nothing
            events::Vector{Symbol}
        end
        @eval function $name(events)
            $name((;), FormulationOptions(), nothing, events)
        end
        op=GlobalRef(owner, nameof(operation))
        @eval function $op(selected::$name, kind::Union{Val{:self}, Val{:mutual}},
                s::Val{2}, t::Val{2},
                functor, workspace)
            push!(selected.events, $(QuoteNode(event)))
            pair = functor.input.pair
            # Manufactured coefficients with diagonal dominance test the stage order,
            # independently of any built-in author's numerical implementation.
            return $(owner === EI ? :(1e-4 + 1e-3im) : :(1e9)) *
                   (pair.row == pair.column ? 10 : 1)
        end
        @eval LineCableModels.formulation_options(::Expression{
            <:$name, typeof($operation)}) = FormulationOptions()
    end

    for selected_type in
        (LayerImpedance, LayerPotential, HalfSpaceImpedance, HalfSpacePotential,
        TypedLayerImpedance, GenericLayerImpedance,
        CoupledImpedance, SurfaceLaw, SpectralSurface,
        DispersiveEarth, InsulationReactance, ConstantResistivity, ScaledResistivity,
        ExponentialResistivity, DispersiveSoil, ScaledSoil, OhmicDielectric, InsulationLaw, SemiconLaw,
        MeanEarth, SquaredBottomEarth, UserCoaxialShunt, UserCoaxialPipe,
        CountedInsulationZ, CountedInsulationY,
        CountedSemiconY, CountedEarthZ, CountedEarthP)
        id=Symbol(nameof(selected_type))
        @eval begin
            LineCableModels.formula_id(::$selected_type) = $(QuoteNode(id))
            LineCableModels.formula_id(::Type{<:$selected_type}) = $(QuoteNode(id))
            LineCableModels.description(::$selected_type; compact::Bool = false) = $(string(id))
            LineCableModels.description(::Type{<:$selected_type}; compact::Bool = false) = $(string(id))
            LineCableModels.formulation_options(selected::$selected_type) = selected.options
        end
        if selected_type <: Union{E.EarthImpedanceFormulation, E.EarthAdmittanceFormulation}
            @eval Base.NamedTuple(selected::$selected_type) = (
                identifier = formula_id(selected),
                parameters = selected.parameters, options = selected.options.data,
                equivalent_earth = selected.equivalent_earth === nothing ? nothing :
                                   NamedTuple(selected.equivalent_earth))
        else
            @eval Base.NamedTuple(selected::$selected_type) = (
                identifier = formula_id(selected),
                parameters = selected.parameters, options = selected.options.data)
        end
    end
end
