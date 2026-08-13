@inline function Base.getproperty(f::Homogeneous, name::Symbol)
	if name === :s || name === :t || name === :Γx || name === :γ1 || name === :γ2 ||
	   name === :μ2
		return getproperty(from_kernel(f), name)
	end
	return getfield(f, name)  # subtype-specific fields (if any)
end