@testmodule EarthInteractionFixtures begin
    using LineCableModels
    const E=LineCableModels.Engine
    const G=LineCableModels.Commons
    const EI, EA=E.EarthImpedance, E.EarthAdmittance
    const Expression=LineCableModels.Expression

    # The indexed methods calculate the coefficients. The engine supplies reuse.
    for (name, parent, operation, multiplier) in
        ((:PairImpedance, E.EarthImpedanceFormulation, EI.earth_impedance, 1e-4+1e-3im),
        (:PairPotential, E.EarthAdmittanceFormulation, EA.earth_potential_coefficient, 1e9))
        @eval struct $name <: $parent
            parameters::NamedTuple{(), Tuple{}}
            options::FormulationOptions{NamedTuple{(), Tuple{}}}
            equivalent_earth::Nothing
            calls::Vector{Tuple{Int, Int}}
        end
        @eval $name()=$name((;), FormulationOptions(), nothing, Tuple{Int, Int}[])
        @eval LineCableModels.formulation_options(::Expression{
            <:$name, typeof($operation)})=FormulationOptions()
        @eval begin
            LineCableModels.formula_id(::$name)=$(QuoteNode(name))
            LineCableModels.formula_id(::Type{$name})=$(QuoteNode(name))
            LineCableModels.description(::$name; compact::Bool = false)=$(string(name))
            LineCableModels.description(::Type{$name}; compact::Bool = false)=$(string(name))
            LineCableModels.formulation_options(selected::$name)=selected.options
            Base.NamedTuple(selected::$name)=(identifier = formula_id(selected),
                parameters = selected.parameters,
                options = selected.options.data, equivalent_earth = nothing)
        end
        operation_name=GlobalRef(parentmodule(operation), nameof(operation))
        for layer in (1, 2), kind in (:self, :mutual)

            @eval function $operation_name(selected::$name, ::Val{$(QuoteNode(kind))},
                    ::Val{$layer}, ::Val{$layer}, functor, pair, workspace)
                push!(selected.calls, (pair.row, pair.column))
                radii=workspace.plan.geometry.radius
                coefficient=10(pair.row==pair.column) + abs(pair.heights[2]) +
                            2abs(pair.heights[1]) +
                            pair.separation + radii[pair.row] + 3radii[pair.column] +
                            functor.state.rho[2]/100 + imag(functor.state.jω)/(2pi*1000)
                return $multiplier*coefficient
            end
        end
    end

    function E.earth_bindings(::Union{PairImpedance, PairPotential}, binding::NamedTuple, geometry::NamedTuple)
        inputs=map(binding.interactions) do interaction
            pair=interaction.pair
            (pair.row==pair.column, pair.layers, pair.heights, pair.separation,
                geometry.radius[pair.row], geometry.radius[pair.column])
        end
        return merge(binding, (reuse_inputs = inputs,))
    end

    struct IntegralImpedance{P} <: E.EarthImpedanceFormulation
        parameters::P
        options::FormulationOptions{NamedTuple{(), Tuple{}}}
        equivalent_earth::Nothing
        calls::Base.RefValue{Int}
    end
    IntegralImpedance(description;
        positions = false)=IntegralImpedance(
        (; description, positions), FormulationOptions(), nothing, Ref(0))
    LineCableModels.formulation_options(::Expression{
        <:IntegralImpedance, typeof(EI.earth_impedance)})=FormulationOptions()
    function EI.earth_impedance(
            selected::IntegralImpedance, ::Union{Val{:self}, Val{:mutual}},
            ::Val{2}, ::Val{2}, functor, pair, workspace)
        selected.calls[]+=1
        integral=E.SpectralIntegral(x->complex(exp(-x)*cos(20x)))
        context=selected.parameters.positions ?
                (receiver = pair.row, source = pair.column,
            frequency = imag(functor.state.jω)/(2pi), term = :test) :
                selected.parameters.description
        value,
        _=E.integrate(integral, Val(:quad), (rtol = 1e-14, atol = 0.0, maxevals = 15),
            workspace.buffers; observations = workspace.buffers.observations, context)
        return value
    end
    E.earth_bindings(::IntegralImpedance,
        binding::NamedTuple,
        geometry::NamedTuple) = merge(binding, (reuse_inputs = fill((), length(binding.interactions)),))
    G.initialize_buffers(::IntegralImpedance,
        ::Type{T},
        input,
        plan,
        buffers) where {T} =
        G.initialize_buffers(E.SpectralIntegral, Val(:quad), T, input, plan, buffers)

    function workspace(problem, impedance = PairImpedance(); trace = false)
        selected=Formulation(
            earth_impedance = impedance, earth_admittance = PairPotential();
            options = (
                reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
        execution=E.computation_options(LineCableModelsCoaxial, ComputationOptions(; trace))
        blueprints=only(E.flatten(
            LineCableModelsCoaxial(), problem.system.designs, eltype(problem), [selected]))
        work=E.LineParametersWorkspace(problem, selected, execution, blueprints)
        E.materials!(work, selected)
        E.materials!(work, selected, 1)
        return work
    end
end
