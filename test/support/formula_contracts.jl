@testmodule FormulaContractModels begin
    using LineCableModels
    const E = LineCableModels.Engine
    const II = E.InternalImpedance
    const EI = E.EarthImpedance
    const EA = E.EarthAdmittance
    const EP = LineCableModels.Earth
    const FD = EP.FrequencyDependent
    const TD = LineCableModels.Materials.TemperatureDependent
    const EH = EP.EquivalentHomogeneous
    const IA = E.InsulationAdmittance
    const SA = E.SemiconAdmittance
    const FM = LineCableModels.FormulaMethod
    const calls = Tuple[]

    struct BufferReplacement <: E.AbstractFormulation end
    E.initialize_buffers(::BufferReplacement, ::Type, input, invariants,
        buffers) = merge(buffers, (destination = copy(buffers.destination),))

    # These are consumer-owned scientific selections, not registrations or
    # private constructors of built-in Formula types. The catalogues stay closed.
    for (name, parent, owner, operation) in (
        (:LayerImpedance, E.EarthImpedanceFormulation, EI, EI.earth_impedance),
        (:LayerPotential, E.EarthAdmittanceFormulation, EA, EA.earth_potential_coefficient))
        @eval struct $name{A, P, O, R} <: $parent
            assumptions::A
            parameters::P
            options::O
            equivalent_earth::R
            initialized::Vector{Tuple}
        end
        operation_name = GlobalRef(owner, nameof(operation))
        for s in 1:3, t in 1:3, kind in (s == t ? (:self, :mutual) : (:mutual,))
            @eval function $operation_name(
                    selected::$name, ::Val{$(QuoteNode(kind))},
                    ::Val{$s}, ::Val{$t}, functor, pair, workspace)
                push!(calls,
                    ($(QuoteNode(nameof(owner))), functor.state.jω,
                        pair.row, pair.column, pair.layers, pair.heights,
                        copy(functor.state.rho), functor.binding.physical_pair))
                coefficient = 11 * $s + 17 * $t + 3 * pair.row + 5 * pair.column +
                              real(functor.state.jω / (2pi*im)) / 100 +
                              (pair.row == pair.column ? 101 : 0)
                if haskey(functor.options.data, :integration)
                    integral = E.SpectralIntegral(λ -> complex(exp(-2λ)))
                    value,
                    _ = E.integrate(functor.options.data.integration.method,
                        integral, functor.options.data.integration.options,
                        workspace.buffers.quadrature)
                    coefficient *= 2value
                end
                return selected.parameters.scale *
                       $(owner === EI ? :(coefficient * (1e-4 + 1e-3im)) :
                         :(coefficient * 1e9))
            end
        end
        @eval LineCableModels.formulation_options(::FM{
            <:$name, typeof($operation)}) = FormulationOptions()
    end
    LineCableModels.formulation_options(::FM{<:LayerImpedance, typeof(EI.earth_impedance),
        A}) where {A <:
                   Tuple{Union{Val{:self}, Val{:mutual}}, Val{1},
        Val{1}}} = FormulationOptions(integration = (method = :quad, options = (;)))

    function E.initialize_buffers(
            selected::LayerImpedance, ::Type{T}, input, invariants, buffers) where {T}
        layers = [E.layer_index(interaction.pair)
                  for call in invariants.earth_bindings.earth_impedance.cases
                  if call.selection === selected for interaction in call.interactions]
        initialized = (1, 1) in layers ?
                      E.initialize_buffers(Val(:quad), T, input, invariants, buffers) :
                      buffers
        push!(selected.initialized, (layers, initialized.quadrature))
        return initialized
    end

    function selection(
            owner; options::Union{NamedTuple, FormulationOptions} = FormulationOptions(),
            layers = 3:3, scale = 1.0)
        options = options isa NamedTuple ? FormulationOptions(options) : options
        physical = (media = Val(:stratified), layers = layers,
            permittivity = :nonzero)
        selected_type = owner === EI ? LayerImpedance : LayerPotential
        return selected_type(physical, (scale = scale,), options, nothing, Tuple[])
    end

    # An algebraic coupled equation uses the same main workspace and explicit
    # calculation action without any formula-owned lifecycle.
    struct CoupledImpedance{A, P, O} <: E.EarthImpedanceFormulation
        assumptions::A
        parameters::P
        options::O
        equivalent_earth::Nothing
        solves::Vector{ComplexF64}
    end
    CoupledImpedance() = CoupledImpedance(
        (media = Val(:homogeneous), layers = 2:2,
            permittivity = :positive),
        (;), FormulationOptions(), nothing, ComplexF64[])
    LineCableModels.formulation_options(::FM{
        <:CoupledImpedance, typeof(EI.earth_impedance)}) = FormulationOptions()
    function E.initialize_buffers(
            ::CoupledImpedance, ::Type{T}, input, invariants, buffers) where {T}
        geometry=invariants.geometry
        return merge(buffers,
            (coupled = Matrix{Complex{T}}(undef,
                length(geometry.radius), length(geometry.radius)),))
    end
    function E.earth!(Z, P, selected::CoupledImpedance, ::Nothing,
            binding, ::Nothing, materials, workspace, frequency)
        values = workspace.buffers.coupled
        s = workspace.input.jω[frequency]
        n = size(values, 1)
        for p in 1:n, q in 1:n

            values[p, q] = s*1e-6*(sum(1:n)+p+2q+(p==q ? n : 0))
        end
        push!(selected.solves, s)
        return E.earth!(Z, selected, binding, materials, s, workspace, nothing)
    end
    for source in 1:2, target in 1:2,
        kind in (source==target ? (:self, :mutual) : (:mutual,))
        @eval EI.earth_impedance(::CoupledImpedance, ::Val{$(QuoteNode(kind))},
            ::Val{$source}, ::Val{$target}, functor, pair,
            workspace) = workspace.buffers.coupled[pair.row, pair.column]
    end

    struct SurfaceLaw{Kinds, P, O} <: E.InternalImpedanceFormulation
        parameters::P
        options::O
        preparations::Vector{Tuple}
        evaluations::Vector{Tuple}
    end
    function SurfaceLaw(; kinds = (:inner, :outer, :transfer),
            coefficients = (inner = 2.0+1im, outer = 2.0+1im, transfer = 0.5+0.1im))
        parameters = (coefficients = coefficients,)
        options = FormulationOptions(NamedTuple{kinds}(map(_ -> (;), kinds)))
        SurfaceLaw{kinds, typeof(parameters), typeof(options)}(
            parameters, options, Tuple[], Tuple[])
    end
    function (selected::SurfaceLaw)(r_in, r_ex, rho, mu_r, jω)
        push!(selected.preparations, (r_in, r_ex, rho, mu_r, jω))
        state = (serial = length(selected.preparations), rho = rho, radius = r_ex)
        return II.Functor(selected, state, selected.options)
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
    LineCableModels.formulation_options(::FM{
        <:SurfaceLaw, typeof(II.internal_impedance)}) = FormulationOptions()

    struct SpectralSurface{P, O} <: E.InternalImpedanceFormulation
        parameters::P
        options::O
        seen::Vector{Tuple}
    end
    function SpectralSurface(method = :quad)
        seed = SpectralSurface((;), FormulationOptions(outer = (;)), Tuple[])
        outer = formulation_options(FM(seed, II.internal_impedance, Val(:outer)),
            FormulationOptions(integration = (method = method,)))
        return SpectralSurface((;), FormulationOptions(outer = outer.data), seed.seen)
    end
    function E.initialize_buffers(
            selected::SpectralSurface, ::Type{T}, input, invariants, buffers) where {T}
        push!(selected.seen, (:initialize, buffers))
        return E.initialize_buffers(Val(:quad), T, input, invariants, buffers)
    end
    (selected::SpectralSurface)(
        r_in, r_ex, rho, mu_r, jω) = II.Functor(selected, (jω = jω,), selected.options)
    LineCableModels.formulation_options(::FM{
        <:SpectralSurface, typeof(II.internal_impedance),
        Tuple{Val{:outer}}}) = FormulationOptions(integration = (
        method = :quad, options = (;)))
    function II.internal_impedance(selected::SpectralSurface, ::Val{:outer}, functor, workspace)
        push!(selected.seen, (functor.options.data.integration.method, workspace))
        integral=E.SpectralIntegral(λ->complex(exp(-2λ)))
        value,
        _=E.integrate(functor.options.data.integration.method, integral,
            functor.options.data.integration.options,
            workspace===nothing ? nothing : workspace.buffers.quadrature)
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
    function FD.earth_material(
            selected::DispersiveEarth, material, frequency, parameters, options, workspace)
        push!(selected.seen, (material.rho, frequency, workspace))
        push!(selected.events, :fd)
        return EP.EarthMaterial(
            material.rho/(1+(frequency/parameters.scale)^parameters.exponent),
            material.eps_r, material.mu_r)
    end
    LineCableModels.formulation_options(::FM{
        <:DispersiveEarth, typeof(FD.earth_material)}) = FormulationOptions()

    struct InsulationReactance{P, O} <: E.InsulationImpedanceFormulation
        parameters::P
        options::O
    end
    InsulationReactance(inductance = 2) = InsulationReactance((inductance = inductance,), FormulationOptions())
    E.InsulationImpedance.insulation_impedance(::InsulationReactance,
        r_in, r_ex, mu_r, s, parameters, options, workspace) = parameters.inductance*s
    LineCableModels.formulation_options(::FM{<:InsulationReactance,
        typeof(E.InsulationImpedance.insulation_impedance)}) = FormulationOptions()

    struct ConstantResistivity{P, O} <: TD.TemperatureDependentFormulation
        parameters::P
        options::O
        seen::Vector{Tuple}
    end
    ConstantResistivity(value) = ConstantResistivity((rho = value,), FormulationOptions(), Tuple[])
    function TD.temperature_resistivity(selected::ConstantResistivity, material,
            temperature, parameters, options, workspace)
        push!(selected.seen, (material, temperature, workspace))
        return parameters.rho
    end
    LineCableModels.formulation_options(::FM{
        <:ConstantResistivity, typeof(TD.temperature_resistivity)}) = FormulationOptions()

    struct ScaledResistivity{P, O} <: TD.TemperatureDependentFormulation
        parameters::P
        options::O
        seen::Vector{Tuple}
    end
    ScaledResistivity(scale = 2) = ScaledResistivity((scale = scale,), FormulationOptions(), Tuple[])
    function TD.temperature_resistivity(selected::ScaledResistivity, material,
            temperature, parameters, options, workspace)
        push!(selected.seen, (material, temperature, workspace))
        return parameters.scale*material.rho
    end
    LineCableModels.formulation_options(::FM{
        <:ScaledResistivity, typeof(TD.temperature_resistivity)}) = FormulationOptions()

    struct ExponentialResistivity{P, O} <: TD.TemperatureDependentFormulation
        parameters::P
        options::O
    end
    ExponentialResistivity(scale = 1000.0) = ExponentialResistivity((scale = scale,), FormulationOptions())
    TD.temperature_resistivity(
        ::ExponentialResistivity, m, t, p, o, w) = m.rho*exp((t-m.T0)/p.scale)
    LineCableModels.formulation_options(::FM{
        <:ExponentialResistivity, typeof(TD.temperature_resistivity)}) = FormulationOptions()

    struct DispersiveSoil{P, O} <: FD.FrequencyDependentFormulation
        parameters::P
        options::O
        seen::Vector{Float64}
    end
    DispersiveSoil() = DispersiveSoil(
        (rho_scale = 1000, epsilon_scale = 2000, mu_scale = 10000),
        FormulationOptions(), Float64[])
    function FD.earth_material(selected::DispersiveSoil, m, f, p, o, w)
        push!(selected.seen, f)
        EP.EarthMaterial(m.rho/(1+f/p.rho_scale), m.eps_r*(1+f/p.epsilon_scale), m.mu_r*(1+f/p.mu_scale))
    end
    LineCableModels.formulation_options(::FM{
        <:DispersiveSoil, typeof(FD.earth_material)}) = FormulationOptions()

    struct ScaledSoil{P, O} <: FD.FrequencyDependentFormulation
        parameters::P
        options::O
    end
    ScaledSoil(; rho = 1, epsilon = 1,
        mu = 1) = ScaledSoil((rho = rho, epsilon = epsilon, mu = mu), FormulationOptions())
    FD.earth_material(::ScaledSoil, m, f, p, o,
        w) = EP.EarthMaterial(m.rho*p.rho, m.eps_r*p.epsilon, m.mu_r*p.mu)
    LineCableModels.formulation_options(::FM{
        <:ScaledSoil, typeof(FD.earth_material)}) = FormulationOptions()

    struct OhmicDielectric{P, O} <: E.InsulationAdmittanceFormulation
        parameters::P
        options::O
        temperatures::Vector{Tuple}
    end
    OhmicDielectric() = OhmicDielectric((;), FormulationOptions(), Tuple[])
    function IA.insulation_material(selected::OhmicDielectric, m, f, t, p, o, w)
        push!(selected.temperatures, (m.T0, t))
        complex(inv(m.rho), 2pi*f*8.8541878128e-12*m.eps_r)
    end
    LineCableModels.formulation_options(::FM{
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
        @eval function $operation_name(
                ::$name, material, frequency, temperature, parameters, options, workspace)
            return complex(
                oftype(frequency, parameters.scale)/material.rho, frequency*material.eps_r+temperature)
        end
        @eval LineCableModels.formulation_options(::FM{
            <:$name, typeof($operation)}) = FormulationOptions()
    end

    struct MeanEarth{P, O} <: EH.AbstractRule
        parameters::P
        options::O
        seen::Vector{Tuple}
    end
    MeanEarth() = MeanEarth((;), FormulationOptions(), Tuple[])
    function EH.equivalent_material(selected::MeanEarth, ::Val{Kind}, ::Val{S}, ::Val{T},
            rho, epsilon, mu, model, pair, frequency, parameters,
            options, workspace) where {Kind, S, T}
        push!(selected.seen, (copy(rho), pair, frequency, workspace))
        soils=2:length(rho)
        return EP.EarthMaterial(sum(rho[soils])/length(soils),
            sum(epsilon[soils])/length(soils), sum(mu[soils])/length(soils))
    end
    LineCableModels.formulation_options(::FM{
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
            rho, epsilon, mu, model, pair, frequency, parameters,
            options, workspace) where {Kind, S, T}
        push!(selected.events, :ehem)
        push!(selected.workspaces, workspace)
        return EP.EarthMaterial(last(rho)^2/parameters.scale, last(epsilon), last(mu))
    end
    LineCableModels.formulation_options(::FM{
        <:SquaredBottomEarth, typeof(EH.equivalent_material)}) = FormulationOptions()
    LineCableModels.validate(
        binding::FM{<:Union{EI.Formula{:unified}, EA.Formula{:unified}}},
        ::SquaredBottomEarth) = binding

    struct FixedModalMaps{P, O} <: LineCableModels.AbstractFormulation
        parameters::P
        options::O
    end
    FixedModalMaps(voltage::AbstractArray{<:Number, 3},
        current::AbstractArray{<:Number, 3}) = FixedModalMaps(
        (voltage = voltage, current = current), FormulationOptions())
    function LineCableModels.Transforms.modal_operators(
            ::FixedModalMaps, source, parameters, options, workspace)
        return LineCableModels.Transforms.ModalOperators(copy(parameters.voltage), copy(parameters.current))
    end
    LineCableModels.formulation_options(::FM{<:FixedModalMaps,
        typeof(LineCableModels.Transforms.modal_operators)}) = FormulationOptions()

    struct UserCoaxialShunt{P, O} <: E.ShuntModelFormulation
        parameters::P
        options::O
        preparations::Base.RefValue{Int}
    end
    UserCoaxialShunt() = UserCoaxialShunt((;), FormulationOptions(), Ref(0))
    function E.internal_shunt_response(selected::UserCoaxialShunt, design::CableDesign,
            geometry, T, methods, solutions, design_index)
        selected.preparations[]+=1
        response=E.internal_shunt_response(E.ShuntModel.Formula(:coaxial),
            design, geometry, T, methods, solutions, design_index)
        merge(response, (details = merge(response.details, (requested = :UserCoaxialShunt,)),))
    end

    struct CoaxialPipePolicy{P, O} <: E.PipeImpedanceFormulation
        parameters::P
        options::O
    end
    CoaxialPipePolicy() = CoaxialPipePolicy((;), FormulationOptions())
    E.Formulation(::LineCableModelsCoaxial, ::CoaxialPipePolicy, ::Val{:coaxial}) = nothing

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
        if owner === E.InsulationImpedance
            @eval function $op(
                    selected::$name, r_in, r_ex, mu_r, s, parameters, options, workspace)
                push!(selected.events, $(QuoteNode(event)))
                push!(selected.workspaces, workspace)
                $op(selected.base, r_in, r_ex, mu_r, s, parameters, options, workspace)
            end
        else
            @eval function $op(selected::$name, material, frequency,
                    temperature, parameters, options, workspace)
                push!(selected.events, $(QuoteNode(event)))
                push!(selected.workspaces, workspace)
                $op(selected.base, material, frequency,
                    temperature, parameters, options, workspace)
            end
        end
        @eval LineCableModels.formulation_options(::FM{
            <:$name, typeof($operation)}) = FormulationOptions()
    end
    for (name, parent, owner, operation, event) in (
        (:CountedEarthZ, E.EarthImpedanceFormulation, EI, EI.earth_impedance, :earth_z),
        (:CountedEarthP, E.EarthAdmittanceFormulation,
        EA, EA.earth_potential_coefficient, :earth_y))
        @eval struct $name{A, P, O} <: $parent
            assumptions::A
            parameters::P
            options::O
            equivalent_earth::Nothing
            events::Vector{Symbol}
        end
        @eval function $name(events)
            $name((media = Val(:homogeneous), layers = 2:2, permittivity = :positive),
                (;), FormulationOptions(), nothing, events)
        end
        op=GlobalRef(owner, nameof(operation))
        @eval function $op(selected::$name, kind::Union{Val{:self}, Val{:mutual}},
                s::Val{2}, t::Val{2},
                functor, pair, workspace)
            push!(selected.events, $(QuoteNode(event)))
            # Manufactured diagonal-dominant coefficients test stage order,
            # independently of any built-in author's numerical implementation.
            return $(owner === EI ? :(1e-4 + 1e-3im) : :(1e9)) *
                   (pair.row == pair.column ? 10 : 1)
        end
        @eval LineCableModels.formulation_options(::FM{
            <:$name, typeof($operation)}) = FormulationOptions()
    end

    for selected_type in (LayerImpedance, LayerPotential, CoupledImpedance, SurfaceLaw, SpectralSurface,
        DispersiveEarth, InsulationReactance, ConstantResistivity, ScaledResistivity,
        ExponentialResistivity, DispersiveSoil, ScaledSoil, OhmicDielectric, InsulationLaw, SemiconLaw,
        MeanEarth, SquaredBottomEarth, FixedModalMaps, UserCoaxialShunt, CoaxialPipePolicy,
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
                identifier = formula_id(selected), assumptions = selected.assumptions,
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
