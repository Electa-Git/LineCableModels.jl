@testmodule SavedResults begin
    using LineCableModels
    import Measurements

    "Save `value` with `save` and read it back with `import_data(path, :result)`."
    function round_trip(value)
        return mktempdir() do directory
            import_data(save(value, joinpath(directory, "result.jls")), :result)
        end
    end

    # Whether `a` and `b` have the same type at every node and the same leaves. Loading gives
    # uncertainty sources new identities. An uncertain leaf compares its nominal value, its
    # uncertainty and its derivative for each source, with each source named by its nominal
    # value and standard deviation.
    same(a, b) = typeof(a) === typeof(b) && same_contents(a, b)
    same(a::Measurements.Measurement, b::Measurements.Measurement) =
        typeof(a) === typeof(b) && isequal(Measurements.value(a), Measurements.value(b)) &&
        isequal(Measurements.uncertainty(a), Measurements.uncertainty(b)) &&
        derivatives(a) == derivatives(b)
    derivatives(x) = sort!([(source[1], source[2], Measurements.derivative(x, source))
        for source in keys(Measurements.uncertainty_components(x))])

    same_contents(a::Complex, b) = same(real(a), real(b)) && same(imag(a), imag(b))
    same_contents(a::Union{Number, Symbol, AbstractString, Nothing, Missing, Type, Module},
        b) = isequal(a, b)
    same_contents(a::AbstractArray, b) = size(a) == size(b) && all(splat(same), zip(a, b))
    same_contents(a::AbstractDict, b) = length(a) == length(b) &&
        all(haskey(b, key) && same(value, b[key]) for (key, value) in a)
    same_contents(a::Union{Tuple, NamedTuple}, b) = all(splat(same), zip(values(a), values(b)))
    function same_contents(a, b)
        return all(1:nfields(a)) do i
            isdefined(a, i) == isdefined(b, i) &&
                (!isdefined(a, i) || same(getfield(a, i), getfield(b, i)))
        end
    end
end
