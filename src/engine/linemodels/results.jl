"""
$(TYPEDEF)

Store rational admittance in the phase domain and groups for delayed propagation.
Residues have axes phase × phase × pole.

$(TYPEDFIELDS)
"""
struct WidebandModel{T}
    "Poles for admittance \\[1/s\\], residues \\[dimensionless\\] and direct matrix \\[S\\]."
    admittance::NamedTuple{(:poles, :residues, :direct),
        Tuple{Vector{Complex{T}}, Array{Complex{T}, 3}, Matrix{T}}}
    "Propagation groups with delays \\[s\\], poles \\[1/s\\] and residues \\[dimensionless\\]."
    propagation::Vector{NamedTuple{(:delay, :poles, :residues),
        Tuple{T, Vector{Complex{T}}, Array{Complex{T}, 3}}}}
end

"""
$(TYPEDEF)

Retain a real, constant modal-to-phase transformation and real-pole sections.
The transformation maps modal quantities to phase quantities.

$(TYPEDFIELDS)
"""
struct ModalModel{T}
    "Transformation matrix with real entries, phase × mode \\[dimensionless\\]."
    transformation::Matrix{T}
    "Frequency used to select the transformation \\[Hz\\]."
    transformation_frequency::T
    "Modal impedance sections. A section adds real-pole terms to a constant and delays the sum by `delay` \\[s\\]."
    impedance::Vector{NamedTuple{(:constant, :poles, :residues, :delay),
        Tuple{T, Vector{T}, Vector{T}, T}}}
    "Modal propagation sections. A section adds real-pole terms to a constant and delays the sum by `delay` \\[s\\]."
    propagation::Vector{NamedTuple{(:constant, :poles, :residues, :delay),
        Tuple{T, Vector{T}, Vector{T}, T}}}
end

"""
$(TYPEDEF)

Retain one fitted line model and its measured fitting diagnostics.
The errors are relative fractions. The `f_max` field records frequency \\[Hz\\].

$(TYPEDFIELDS)
"""
struct LineModel{M, T, D} <: AbstractCoreResult
    "Fitted wideband or modal payload."
    model::M
    "Phase labels in physical order."
    coordinates::Vector{String}
    "Segment length \\[m\\]."
    line_length::T
    "Frequencies for the fit \\[Hz\\]."
    frequencies::Vector{T}
    "Admittance and propagation errors: per-entry RMS and maximum matrices, with `f_max` giving the frequency of each maximum."
    fitting_error::NamedTuple{(:admittance, :propagation), Tuple{
        NamedTuple{(:rms, :max, :f_max), Tuple{Matrix{T}, Matrix{T}, Matrix{T}}},
        NamedTuple{(:rms, :max, :f_max), Tuple{Matrix{T}, Matrix{T}, Matrix{T}}}}}
    "Passivity evidence and perturbation."
    passivity::NamedTuple{(:checked, :passive, :bands, :worst, :f_worst, :perturbation)}
    "Direct-current targets and achieved values, or nothing."
    dc::NamedTuple{(:targets, :achieved)}
    "Supplemental computation record."
    details::D
end

"""Select poles \\[1/s\\] by admittance or impedance and propagation group."""
function poles end
"""Select residues \\[dimensionless\\] by admittance or impedance and propagation group."""
function residues end
"""Select propagation delays \\[s\\]."""
function delays end
"""Select relative errors of the fit and the frequencies of their maxima \\[Hz\\]."""
function fitting_error end
"""Select the retained passivity evidence."""
function passivity end

poles(model::WidebandModel) = [model.admittance.poles,
    [group.poles for group in model.propagation]...]
poles(model::ModalModel) = [[section.poles for section in model.impedance]...,
    [section.poles for section in model.propagation]...]
residues(model::WidebandModel) = [model.admittance.residues,
    [group.residues for group in model.propagation]...]
residues(model::ModalModel) = [[section.residues for section in model.impedance]...,
    [section.residues for section in model.propagation]...]
delays(model::WidebandModel) = [group.delay for group in model.propagation]
delays(model::ModalModel) = [section.delay for section in model.propagation]
poles(result::LineModel) = poles(result.model)
residues(result::LineModel) = residues(result.model)
delays(result::LineModel) = delays(result.model)
fitting_error(result::LineModel) = result.fitting_error
passivity(result::LineModel) = result.passivity
details(result::LineModel) = result.details

function completed_inputs(model::LineModel)
    return get(details(model).data,:inputs,(;))
end

"""
$(TYPEDSIGNATURES)

Evaluate fitted characteristic admittance at complex frequencies and DC.
Wideband and modal models use, respectively,

```math
Y_c(s) = D + \\sum_k \\frac{R_k}{s-p_k}, \\qquad
Y_c(s) = T\\,\\mathrm{diag}(1/Z_{c,m}(s))\\,T^\\mathsf{T},
```

where ``T`` is the retained current transformation and ``Z_{c,m}`` is the
fitted modal impedance.

# Arguments

- `result`: Completed fitted line model.
- `s`: Vector of complex frequencies \\[rad/s\\]. Zero selects DC.

# Returns

- Phase-domain admittance \\[S\\], with axes `terminal × terminal × frequency`
  in `result.coordinates` and `s` order.
"""
function Yc(result::LineModel{<:WidebandModel{T}}, s::AbstractVector{<:Complex}) where {T}
    admittance = result.model.admittance
    n = size(admittance.direct, 1)
    values = Array{promote_type(Complex{T}, eltype(s)), 3}(undef, n, n, length(s))
    for (k, point) in enumerate(s)
        @views values[:, :, k] .= admittance.direct
        for j in eachindex(admittance.poles)
            @views values[:, :, k] .+= admittance.residues[:, :, j] ./ (point - admittance.poles[j])
        end
    end
    return values
end

function Yc(result::LineModel{<:ModalModel{T}}, s::AbstractVector{<:Complex}) where {T}
    model = result.model
    transformation = model.transformation
    n = size(transformation, 1)
    S = promote_type(Complex{T}, eltype(s))
    values = Array{S, 3}(undef, n, n, length(s))
    modal = Vector{S}(undef, n)
    for (k, point) in enumerate(s)
        for (mode, section) in enumerate(model.impedance)
            impedance = section.constant + sum(r / (point - p)
                for (p, r) in zip(section.poles, section.residues); init=zero(S))
            modal[mode] = inv(impedance * exp(-point * section.delay))
        end
        @views values[:, :, k] .= transformation * Diagonal(modal) * transpose(transformation)
    end
    return values
end

"""
$(TYPEDSIGNATURES)

Evaluate fitted current-wave propagation at complex frequencies and DC.
Wideband and modal models use, respectively,

```math
H(s) = \\sum_g e^{-s\\tau_g}\\sum_k \\frac{R_{g,k}}{s-p_{g,k}}, \\qquad
H(s) = T\\,\\mathrm{diag}(h_m(s))\\,T^{-1},
```

where ``T`` is the retained current transformation and each fitted modal section
is ``h_m(s)=e^{-s\\tau_m}(c_m+\\sum_k r_{m,k}/(s-p_{m,k}))``.
At DC each delay factor is one and each pole term is ``-r_{m,k}/p_{m,k}``.

# Arguments

- `result`: Completed fitted line model.
- `s`: Vector of complex frequencies \\[rad/s\\]. Zero selects DC.

# Returns

- Phase-domain propagation \\[dimensionless\\], with axes `terminal × terminal × frequency`
  in `result.coordinates` and `s` order.
"""
function H(result::LineModel{<:WidebandModel{T}}, s::AbstractVector{<:Complex}) where {T}
    model = result.model
    n = size(model.admittance.direct, 1)
    values = zeros(promote_type(Complex{T}, eltype(s)), n, n, length(s))
    for (k, point) in enumerate(s), group in model.propagation
        delay = exp(-point * group.delay)
        for j in eachindex(group.poles)
            @views values[:, :, k] .+= group.residues[:, :, j] .* (delay / (point - group.poles[j]))
        end
    end
    return values
end

function H(result::LineModel{<:ModalModel{T}}, s::AbstractVector{<:Complex}) where {T}
    model = result.model
    transformation = model.transformation
    n = size(transformation, 1)
    S = promote_type(Complex{T}, eltype(s))
    values = Array{S, 3}(undef, n, n, length(s))
    modal = Vector{S}(undef, n)
    for (k, point) in enumerate(s)
        for (mode, section) in enumerate(model.propagation)
            modal[mode] = (section.constant + sum(r / (point - p)
                for (p, r) in zip(section.poles, section.residues); init=zero(S))) *
                exp(-point * section.delay)
        end
        @views values[:, :, k] .= (transformation * Diagonal(modal)) / transformation
    end
    return values
end
