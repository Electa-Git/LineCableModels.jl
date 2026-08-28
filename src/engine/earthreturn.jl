@inline function compute_earth_return_matrix!(
        destination::AbstractMatrix{Complex{T}},
        cables::AbstractVector{Int},
        input::NamedTuple,
        earth,
        frequency::Int,
        formulation
) where {T <: Real}
    rho = @view earth.rho[:, frequency]
    epsilon = @view earth.epsilon[:, frequency]
    mu = @view earth.mu[:, frequency]
    s = input.jω[frequency]
    @inbounds for column in eachindex(cables), row in eachindex(cables)

        left = cables[row]
        right = cables[column]
        separation = input.horz_sep[left, right]
        heights = (input.vert[left], input.vert[right])
        kind = row == column ? Val(:self) : Val(:mutual)
        destination[row, column] = formulation(
            kind, heights, separation, rho, epsilon, mu, s
        )
    end
    return destination
end
