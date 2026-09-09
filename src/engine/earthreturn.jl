function earth!(
        destination::AbstractMatrix, bindings::NamedTuple{(
            :selection, :cases)},
        earth, jω, formula, Γ, workspace,
        thickness
)
    bindings.selection === formula ||
        throw(ArgumentError("workspace is bound to a different earth formula selection"))
    foreach(bindings.cases) do binding
        earth!(destination, binding, earth, jω,
            binding.selection, Γ, workspace, thickness)
    end
    return destination
end

function earth!(
        destination::AbstractMatrix{Complex{T}}, binding::NamedTuple{(
            :selection, :declaration, :interactions, :reductions)},
        earth, jω, formula, Γ, workspace,
        thickness
) where {T <: Real}
    thickness = media(formula) === Val(:stratified) ? thickness : nothing
    for interaction in binding.interactions
        index, pair=interaction.index, interaction.pair
        rho=@view earth.rho[:, index]
        epsilon=@view earth.epsilon[:, index]
        mu=@view earth.mu[:, index]
        functor=formula(
            rho, epsilon, mu, jω, pair, binding.declaration; Γ,
            thickness, physical_pair = interaction.physical_pair)
        resources = haskey(binding.declaration.options, :integration) ? workspace : nothing
        destination[pair.row, pair.column]=functor(resources)
    end
    return destination
end

@inline _gamma(::Nothing, frequency::Int) = nothing
@inline _gamma(values::AbstractVector, frequency::Int) = values[frequency]
