"Add one owned value to a mutable collection."
function add! end

"Validate a value and return it unchanged, or throw a native exception."
validate(value) = value

"Return a short scientific description of a formulation or model."
function description end

"Return the maximum number of cable elements admitted by a packing geometry."
function maxfill end

"Return the physical storage basis of a result."
function basis end

"Return the series-impedance values of a result."
function Z end

"Return the shunt-admittance values of a result."
function Y end

function R end
function X end
function L end
function G end
function B end
function C end
function series_impedance end
function shunt_admittance end
function resistance end
function reactance end
function inductance end
function conductance end
function susceptance end
function capacitance end
function frequencies end
function nconductors end
function nfrequencies end
function ncables end
function nphases end

"Return the nominal value of a deterministic or uncertain quantity."
nominal(value) = value
nominal(value::Complex) = complex(nominal(real(value)), nominal(imag(value)))
nominal(values::AbstractArray) = nominal.(values)

"Return the standard uncertainty of a quantity; deterministic numbers return zero."
standard_uncertainty(value::Number) = zero(nominal(value))
standard_uncertainty(::Any) = 0.0

abstract type LineParamsDomain end
struct PhaseDomain <: LineParamsDomain end
struct ModalDomain <: LineParamsDomain end

"Return the domain tag type of a value, or `nothing` when it has no domain."
@inline domain(::Type) = nothing
@inline domain(value) = domain(typeof(value))
