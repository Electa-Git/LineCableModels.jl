"Return a short scientific description of an Engine formulation."
function description end

"Return the series-impedance values of a line-parameter result."
function Z end

"Return the shunt-admittance values of a line-parameter result."
function Y end

function X end
function G end
function B end
function series_impedance end
function shunt_admittance end
function reactance end
function conductance end
function susceptance end
function frequencies end
function nconductors end
function nfrequencies end
"Return absolute numerical errors from an owned comparison result."
function absolute_error end
"Return reference-normalised numerical errors from an owned comparison result."
function relative_error end

"Abstract tag for the physical domain represented by line-parameter matrices."
abstract type LineParamsDomain end

"Return formulation-owned constitutive data for one material class."
function constitutive end
"Tag line parameters expressed in the physical phase domain."
struct PhaseDomain <: LineParamsDomain end
"Tag line parameters obtained from a calculated modal transformation."
struct ModalDomain <: LineParamsDomain end

"Return the domain tag type of a value, or `nothing` when it has no domain."
@inline domain(::Type) = nothing
@inline domain(value) = domain(typeof(value))
