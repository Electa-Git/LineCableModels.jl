"""
$(TYPEDSIGNATURES)

Inspects all numerical data within a `LineParametersProblem` and determines the
common floating-point type. If any value (frequencies, geometric properties,
material properties, or earth properties) is a `Measurement`, the function
returns `Measurement{Float64}`. Otherwise, it returns `Float64`.
"""
_find_common_type(problem::LineParametersProblem) = resolve_T(problem)

function _get_earth_data(
        functor::AbstractEHEMFormulation,
        earth_model::EarthModel,
        freq::Vector{<:REALSCALAR},
        T::DataType
)
    return functor(earth_model, freq, T)
end

"""
Default method for when no EHEM formulation is provided.
"""
function _get_earth_data(::Nothing,
        earth_model::EarthModel,
        freq::AbstractVector{<:REALSCALAR},
        ::Type{T}) where {T <: REALSCALAR}
    nL = length(earth_model.layers)
    nF = length(freq)

    ρ = Matrix{T}(undef, nL, nF)
    ε = Matrix{T}(undef, nL, nF)
    μ = Matrix{T}(undef, nL, nF)

    @inbounds for i in 1:nL
        L = earth_model.layers[i]
        @assert length(L.rho_g) == nF && length(L.eps_g) == nF && length(L.mu_g) == nF
        # Fill elementwise to avoid temp vectors
        for j in 1:nF
            ρ[i, j] = coerce_to_T(L.rho_g[j], T)
            ε[i, j] = coerce_to_T(L.eps_g[j], T)
            μ[i, j] = coerce_to_T(L.mu_g[j], T)
        end
    end

    return (rho_g = ρ, eps_g = ε, mu_g = μ)
end

@inline function _get_outer_radii(cable_map::AbstractVector{Int},
        r_ext::AbstractVector{T},
        r_ins_ext::AbstractVector{T}) where {T <: Real}
    @assert length(cable_map) == length(r_ext) == length(r_ins_ext)
    n = length(cable_map)
    G = maximum(cable_map)
    gmax = fill(zero(T), G)
    @inbounds for i in 1:n
        g = cable_map[i]
        r = max(r_ext[i], r_ins_ext[i])
        if r > gmax[g]
            gmax[g] = r
        end
    end
    return gmax
end

@inline function _calc_horz_sep!(dest::AbstractMatrix{T},
        horz::AbstractVector{T},
        r_ext::AbstractVector{T},
        r_ins_ext::AbstractVector{T},
        cable_map::AbstractVector{Int}) where {T <: Real}
    @assert size(dest, 1) == size(dest, 2) == length(horz) ==
            length(r_ext) == length(r_ins_ext) == length(cable_map)
    n = length(horz)
    gmax = _get_outer_radii(cable_map, r_ext, r_ins_ext)
    @inbounds for j in 1:n, i in 1:n

        if cable_map[i] == cable_map[j]
            dest[i, j] = gmax[cable_map[i]]
        else
            dest[i, j] = abs(horz[i] - horz[j])
        end
    end
    return dest
end

@inline function _get_cable_indices(ws)
    Nc = ws.n_cables
    idxs_by_cable = [Int[] for _ in 1:Nc]
    @inbounds for i in 1:ws.n_phases
        push!(idxs_by_cable[ws.cable_map[i]], i)
    end
    heads = similar(collect(1:Nc))
    @inbounds for c in 1:Nc
        heads[c] = idxs_by_cable[c][1]   # representative (any member) per cable
    end
    return idxs_by_cable, heads
end

@inline function _to_phase!(A::AbstractMatrix{Complex{T}}) where {T <: REALSCALAR}
    m, n = size(A)

    # Right-multiply by T_I (lower-triangular ones): cumulative sum of columns, right→left
    @inbounds for j in (n - 1):-1:1
        @views A[:, j] .+= A[:, j + 1]
    end

    # Left-multiply by T_V^{-1} (bidiagonal solve): cumulative sum of rows, bottom→top
    @inbounds for i in (m - 1):-1:1
        @views A[i, :] .+= A[i + 1, :]
    end

    return A
end

# function _to_phase!(
# 	M::Matrix{Complex{T}},
# ) where {T <: REALSCALAR}
# 	# Check the size of the M matrix (assuming M is NxN)
# 	N = size(M, 1)

# 	# Build the voltage transformation matrix T_V
# 	T_V = Matrix{T}(I, N, N + 1)  # Start with an identity matrix
# 	for i ∈ 1:N
# 		T_V[i, i+1] = -1  # Set the -1 in the next column
# 	end
# 	T_V = T_V[:, 1:N]  # Remove the last column

# 	# Build the current transformation matrix T_I
# 	T_I = tril(ones(T, N, N))  # Lower triangular matrix of ones

# 	# Compute the new impedance matrix M_prime
# 	M = T_V \ M * T_I

# 	return M
# end
