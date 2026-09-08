@testmodule FormulaContractModels begin
    using LineCableModels
    const E = LineCableModels.Engine
    const EI = E.EarthImpedance
    const EA = E.EarthAdmittance
    const EP = LineCableModels.Earth
    const FM = LineCableModels.FormulaMethod
    const calls = Tuple[]

    # Finite, deliberately asymmetric manufactured equations. These do not enter
    # the literature inventory. Exact selector methods are the case authority.
    for owner in (EI, EA)
        operation = owner === EI ? EI.earth_impedance : EA.earth_potential_coefficient
        operation_name = GlobalRef(owner, nameof(operation))
        for s in 1:3, t in 1:3, kind in (s == t ? (:self, :mutual) : (:mutual,))
            @eval function $operation_name(
                    ::Val{:ContractLayers}, ::Val{$(QuoteNode(kind))},
                    ::Val{$s}, ::Val{$t}, functor, pair, workspace)
                push!(calls,
                    ($(QuoteNode(nameof(owner))), functor.state.jω,
                        pair.row, pair.column, pair.layers, pair.heights,
                        copy(functor.state.rho), functor.binding.physical_pair))
                coefficient = 10 * $s + $t + (pair.row == pair.column ? 100 : 0)
                if haskey(functor.options, :integration)
                    integral = E.SpectralIntegral(Val(:cosine), λ -> complex(exp(-λ)),
                        (height = 1.0, separation = 0.0), 1.0)
                    coefficient *= 2 * E.integrate(functor.options.integration.method,
                        integral, functor.options.integration.options, workspace)
                end
                return $(owner === EI ? :(coefficient * (1e-4 + 1e-3im)) :
                         :(coefficient * 1e9))
            end
        end
        @eval function E.hooks(::FM{:ContractLayers, typeof($operation)})
            return (configurable = (:Γ, :contribution),
                defaults = (
                    Γ = (jω, materials, layers) -> zero(jω),
                    air = (jω, μ, σ, ε) -> sqrt(jω * μ * (σ + jω * ε)),
                    earth = (jω, μ, σ, ε) -> sqrt(jω * μ * (σ + jω * ε)),
                    permeability = identity, contribution = nothing))
        end
        @eval LineCableModels.computation_options(::FM{
            :ContractLayers, typeof($operation)}) = (;)
        @eval LineCableModels.Formulation(::LineCableModelsCoaxial,
            selected::$owner.Formula{:ContractLayers}) = selected
    end
    LineCableModels.computation_options(::FM{:ContractLayers,
        typeof(EI.earth_impedance),
        A}) where {
        A <:
        Tuple{Union{Val{:self}, Val{:mutual}}, Val{1},
        Val{1}}} = (integration = (method = :quad, options = (;)),)

    function selection(owner; options = (;), hooks = (;))
        physical = (media = Val(:stratified), layers = 3:3,
            longitudinal = :zero, permittivity = :nonzero)
        return owner.Formula{:ContractLayers, typeof(physical), typeof((;)),
            typeof(hooks), typeof(options), Nothing}(physical, (;), hooks, options, nothing)
    end
end
