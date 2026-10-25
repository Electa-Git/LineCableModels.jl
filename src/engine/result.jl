"""
$(TYPEDEF)

Store the line parameters of a problem: the series impedance and the shunt admittance per
unit length of the retained conductors at each frequency. Each row has the label of its
conductor, and the record keeps the formulation that computed it and the length of the line.

$(TYPEDFIELDS)
"""
struct LineParameters{T <: Real, F <: AbstractFormulation}
    "Series impedance of the line \\[Ω/m\\], by row, column and frequency."
    Z::Array{Complex{T}, 3}
    "Shunt admittance of the line \\[S/m\\], by row, column and frequency."
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
    The matrices that build `Z` and `Y`, or empty arrays when the computation did not keep
    them. `Zin` \\[Ω/m\\] and `Pin` \\[m/F\\] are the matrices inside the cables, by
    conductor. `Zg` \\[Ω/m\\] and `Pg` \\[m/F\\] are the matrices of the earth return, by
    cable. The last index of each is the frequency.
    """
    trace::NamedTuple{(:Zin, :Pin, :Zg, :Pg), NTuple{4, Array{Complex{T}, 3}}}
    function LineParameters{T}(Z::Array{Complex{T}, 3}, Y::Array{Complex{T}, 3},
            frequencies::Vector{T}, terminals::Vector{String}, line_length::T,
            formulation::F,
            trace::NamedTuple{(:Zin, :Pin, :Zg, :Pg), NTuple{4, Array{Complex{T}, 3}}}
    ) where {T <: Real, F <: AbstractFormulation}
        return validate(new{T, F}(
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
    all(isfinite, result.Z) ||
        throw(DomainError(result.Z, "LineParameters.Z must be finite"))
    all(isfinite, result.Y) ||
        throw(DomainError(result.Y, "LineParameters.Y must be finite"))
    return result
end
