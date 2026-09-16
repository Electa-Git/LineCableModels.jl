# Records have structural value semantics, not a NamedTuple collection interface.
Base.:(==)(a::FormulationOptions, b::FormulationOptions) = a.data == b.data
Base.:(==)(a::ComputationOptions, b::ComputationOptions) = a.data == b.data
Base.:(==)(a::ComputationDetails, b::ComputationDetails) = a.data == b.data

Base.isequal(a::FormulationOptions, b::FormulationOptions) = isequal(a.data, b.data)
Base.isequal(a::ComputationOptions, b::ComputationOptions) = isequal(a.data, b.data)
Base.isequal(a::ComputationDetails, b::ComputationDetails) = isequal(a.data, b.data)

Base.hash(value::FormulationOptions, seed::UInt) = hash(value.data, hash(:FormulationOptions, seed))
Base.hash(value::ComputationOptions, seed::UInt) = hash(value.data, hash(:ComputationOptions, seed))
Base.hash(value::ComputationDetails, seed::UInt) = hash(value.data, hash(:ComputationDetails, seed))

function Base.show(io::IO, value::FormulationOptions)
    print(io, "FormulationOptions(")
    show(io, value.data)
    print(io, ')')
end

function Base.show(io::IO, value::ComputationOptions)
    print(io, "ComputationOptions(")
    show(io, value.data)
    print(io, ')')
end

function Base.show(io::IO, value::ComputationDetails)
    print(io, "ComputationDetails(")
    show(io, value.data)
    print(io, ')')
end
