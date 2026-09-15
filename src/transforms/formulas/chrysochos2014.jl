"""
$(TYPEDSIGNATURES)

**Identification.** Levenberg–Marquardt tracking of complex modal eigenpairs,
after A. I. Chrysochos, T. A. Papadopoulos, and G. K. Papagiannis (2014).

**Reference.** A. I. Chrysochos, T. A. Papadopoulos, and G. K. Papagiannis,
“Robust Calculation of Frequency-Dependent Transmission-Line Transformation
Matrices Using the Levenberg–Marquardt Method,” *IEEE Transactions on Power
Delivery*, 29(4), 1621–1629, 2014. DOI: 10.1109/TPWRD.2013.2284504.
"""
function description(::Type{<:Formula{:chrysochos2014}}; compact::Bool=false)
    compact ? "Chrysochos" : "Chrysochos Levenberg–Marquardt modal transformation (2014)"
end

"Delegate the normalized Chrysochos route to the package default implementation."
function modal_operators(
        ::Val{:chrysochos2014}, args...
)
    modal_operators(Val(:default), args...)
end

function computation_options(::FormulaMethod{:chrysochos2014, typeof(modal_operators)})
    computation_options(FormulaMethod(Val(:default), modal_operators))
end

function computation_options(
        ::FormulaMethod{:chrysochos2014, typeof(modal_operators)},
        section::Val,
        defaults::NamedTuple,
        supplied::NamedTuple
)
    computation_options(
        FormulaMethod(Val(:default), modal_operators),
        section,
        defaults,
        supplied
    )
end

:chrysochos2014
