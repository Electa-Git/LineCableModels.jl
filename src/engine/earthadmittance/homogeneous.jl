@inline propagation(::Val{:full}, jω, μ, σ, ε) = sqrt(jω * μ * (σ + jω * ε))
@inline propagation(::Val{:lossless}, jω, μ, σ, ε) = jω * sqrt(μ * ε)
@inline propagation(::Val{:conductive}, jω, μ, σ, ε) = sqrt(jω * μ * σ)
@inline propagation(::Val{:vacuum}, jω, μ, σ, ε) = jω * sqrt(μ * vacuum_permittivity(ε))

@inline function _geometry(pair)
    h_i = abs(pair.heights[1])
    h_j = abs(pair.heights[2])
    # These retained self expressions use the mutual expression at the radius.
    y_ij = pair.row == pair.column ? something(pair.radius) : pair.separation
    H = h_i + h_j
    d_ij = hypot(y_ij, h_i - h_j)
    D_ij = hypot(y_ij, H)
    return (; h_i, h_j, y_ij, H, d_ij, D_ij)
end
