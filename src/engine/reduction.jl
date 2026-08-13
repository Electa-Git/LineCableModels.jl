using LinearAlgebra: BLAS, BlasFloat

function reorder_indices(map::AbstractVector{<:Integer})
	n = length(map)
	phases = Int[]                     # encounter order of phases > 0
	firsts = Int[]
	sizehint!(firsts, n)
	zeros = Int[]
	sizehint!(zeros, n)
	tails = Dict{Int, Vector{Int}}()   # phase => remaining indices

	seen = Set{Int}()
	@inbounds for (i, p) in pairs(map)
		if p > 0
			if !(p in seen)
				push!(seen, p)
				push!(phases, p)
				push!(firsts, i)
			else
				push!(get!(tails, p, Int[]), i)
			end
		else
			push!(zeros, i)
		end
	end

	perm = Vector{Int}(undef, n)
	k = 1
	@inbounds begin
		for i in firsts
			perm[k] = i
			k += 1
		end
		for p in phases
			if haskey(tails, p)
				for i in tails[p]
					perm[k] = i
					k += 1
				end
			end
		end
		for i in zeros
			perm[k] = i
			k += 1
		end
	end
	return perm
end

# Non-mutating reorder (2D)
function reorder_M(M::AbstractMatrix, map::AbstractVector{<:Integer})
	n = size(M, 1)
	n == size(M, 2) == length(map) || throw(ArgumentError("shape mismatch"))
	perm = reorder_indices(map)
	return M[perm, perm], map[perm]
end


"""
	kronify(M, phase_map)
	Kron elimination
"""
function kronify(
	M::Matrix{Complex{T}},
	phase_map::Vector{Int},
) where {T <: REALSCALAR}
	keep = findall(!=(0), phase_map)
	eliminate = findall(==(0), phase_map)

	M11 = M[keep, keep]
	M12 = M[keep, eliminate]
	M21 = M[eliminate, keep]
	M22 = M[eliminate, eliminate]

	return M11 - (M12 * inv(M22)) * M21
end

"""
	kronify!(M, phase_map, Mred)
	Kron's angry little brother (in-place)
"""
function kronify!(
	M::Matrix{Complex{T}},
	phase_map::Vector{Int},
	Mred::Matrix{Complex{T}},
) where {T <: REALSCALAR}
	keep = findall(!=(0), phase_map)
	eliminate = findall(==(0), phase_map)

	M11 = M[keep, keep]
	M12 = M[keep, eliminate]
	M21 = M[eliminate, keep]
	M22 = M[eliminate, eliminate]
	@views @inbounds Mred .= M11 - (M12 * inv(M22)) * M21
	return nothing
end

# In-place: columns tail -= first (from original), then rows tail -= first (after col pass).
function merge_bundles!(M::AbstractMatrix{T}, ph::AbstractVector{<:Integer}) where {T}
	n = size(M, 1)
	(size(M, 2) == n && length(ph) == n) || throw(ArgumentError("shape mismatch"))

	# Encounter-ordered groups (include phase 0)
	groups   = Vector{Vector{Int}}()
	index_of = Dict{Int, Int}()
	@inbounds for (i, p) in pairs(ph)
		gi = get(index_of, p, 0)
		if gi == 0
			push!(groups, Int[])
			gi = length(groups)
			index_of[p] = gi
		end
		push!(groups[gi], i)
	end

	# -------- Pass 1: columns --------
	@inbounds for grp in groups
		length(grp) > 1 || continue
		i1 = grp[1]
		base_col = @view M[:, i1]          # original column kept intact in this pass
		for t in Iterators.drop(eachindex(grp), 1)
			j   = grp[t]
			col = @view M[:, j]
			if (M isa StridedMatrix{T}) && (T <: BlasFloat)
				BLAS.axpy!(-one(T), base_col, col)  # col -= base_col
			else
				col .-= base_col
			end
		end
	end

	# -------- Pass 2: rows --------
	newmap = copy(ph)
	@inbounds for grp in groups
		length(grp) > 1 || continue
		i1 = grp[1]
		base_row = @view M[i1, :]          # uses row after pass 1 (matches Z1→Z2)
		for t in Iterators.drop(eachindex(grp), 1)
			i   = grp[t]
			row = @view M[i, :]
			if (M isa StridedMatrix{T}) && (T <: BlasFloat)
				BLAS.axpy!(-one(T), base_row, row)  # row -= base_row
			else
				row .-= base_row
			end
			newmap[i] = 0
		end
	end
	return M, newmap
end
