"""
$(TYPEDEF)

Store the line parameters of a problem: the series impedance and the shunt admittance per
unit length of the retained conductors at each frequency. Each row has the label of its
conductor, and the record keeps the formulation that computed it and the length of the line.

$(TYPEDFIELDS)
"""
struct LineParameters{T <: Real, F <: AbstractFormulation, N <: NamedTuple}
    "Series impedance of the line \\[Ω/m\\], with one matrix per frequency."
    Z::Array{Complex{T}, 3}
    "Shunt admittance of the line \\[S/m\\], with one matrix per frequency."
    Y::Array{Complex{T}, 3}
    "Frequencies of the scan \\[Hz\\]."
    frequencies::Vector{T}
    """
    Label of each retained conductor: `cable:<i>:<terminal>`, or the labels of a bundle
    joined by `+`.
    """
    terminals::Vector{String}
    "Length of the line \\[m\\]."
    line_length::T
    "Formulation that computed the result."
    formulation::F
    """
    The backend's matrices that build `Z` and `Y`, each an `Array{Complex{T}, 3}` whose last
    index is the frequency, or empty when the computation did not keep them. A backend's
    trace has one type whether or not the computation keeps them. For [`LCMBackend`](@ref),
    `Zin` \\[Ω/m\\] and `Pin` \\[m/F\\] are the matrices inside the cables, by conductor.
    `Zg` \\[Ω/m\\] and `Pg` \\[m/F\\] are the matrices of the earth return, by cable.
    """
    trace::N
    function LineParameters{T}(Z::Array{Complex{T}, 3}, Y::Array{Complex{T}, 3},
            frequencies::Vector{T}, terminals::Vector{String}, line_length::T,
            formulation::F,
            trace::N
    ) where {T <: Real, F <: AbstractFormulation, names,
            N <: NamedTuple{names, <:Tuple{Vararg{Array{Complex{T}, 3}}}}}
        return validate(new{T, F, N}(
            Z, Y, frequencies, terminals, line_length, formulation, trace))
    end
end

# Z and Y are square, of one size, with one slice per frequency and one label per row, and
# finite.
function validate(result::LineParameters)
    n = size(result.Z, 1)
    size(result.Z) == (n, n, length(result.frequencies)) &&
    size(result.Y) == size(result.Z) ||
        throw(DimensionMismatch("LineParameters.Z and Y must be n × n × " *
                                "$(length(result.frequencies)), one slice per frequency; " *
                                "received $(size(result.Z)) and $(size(result.Y))"))
    length(result.terminals) == n || throw(DimensionMismatch(
        "LineParameters.terminals must label the $n rows of Z; received " *
        "$(length(result.terminals)) labels"))
    for matrix in result.trace
        isempty(matrix) || size(matrix, 3) == length(result.frequencies) ||
            throw(DimensionMismatch("LineParameters.trace arrays must have " *
                                    "$(length(result.frequencies)) slices, one per frequency; " *
                                    "received $(size(matrix))"))
    end
    all(isfinite, result.Z) ||
        throw(DomainError(result.Z, "LineParameters.Z must be finite"))
    all(isfinite, result.Y) ||
        throw(DomainError(result.Y, "LineParameters.Y must be finite"))
    return result
end
