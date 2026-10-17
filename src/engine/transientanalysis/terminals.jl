"""
$(TYPEDEF)

Connect one terminal to ground through a series resistance and inductance.
Zero resistance and inductance give a solid connection to ground.

$(TYPEDFIELDS)
"""
struct TerminalImpedance{T}
    "Terminal label."
    terminal::String
    "Resistance \\[Ω\\]."
    R::T
    "Inductance \\[H\\]."
    L::T

    function TerminalImpedance(terminal::AbstractString, R::Real, L::Real)
        r, l = float.(promote(R, L))
        return validate(new{typeof(r)}(String(terminal), r, l))
    end
end

"""
$(TYPEDSIGNATURES)

Assign `R=0` \\[Ω\\] and `L=0` \\[H\\] to a terminal, lifting explicit grids.
"""
function TerminalImpedance(terminal; R=0, L=0)
    return parameterize(TerminalImpedance, TerminalImpedance, (terminal, R, L))
end

function validate(assignment::TerminalImpedance)
    isfinite(assignment.R) && assignment.R >= 0 ||
        throw(DomainError(assignment.R, "resistance must be finite and nonnegative"))
    isfinite(assignment.L) && assignment.L >= 0 ||
        throw(DomainError(assignment.L, "inductance must be finite and nonnegative"))
    return assignment
end

"""
Represent an ideal voltage source between its terminal and ground.
"""
abstract type Waveform end

"""
$(TYPEDEF)

Apply a voltage step of `amplitude` at `t0`.

$(TYPEDFIELDS)
"""
struct Step{T} <: Waveform
    "Terminal label."
    terminal::String
    "Voltage amplitude \\[V\\]."
    amplitude::T
    "Step time \\[s\\]."
    t0::T

    function Step(terminal::AbstractString, amplitude::Real, t0::Real)
        a, t = float.(promote(amplitude, t0))
        return validate(new{typeof(a)}(String(terminal), a, t))
    end
end

"""
$(TYPEDSIGNATURES)

Assign a voltage step \\[V\\] with `t0=0` \\[s\\], lifting explicit grids.
"""
function Step(terminal, amplitude; t0=0)
    return parameterize(Step, Step, (terminal, amplitude, t0))
end

function validate(source::Step)
    isfinite(source.amplitude) || throw(DomainError(source.amplitude, "amplitude must be finite"))
    isfinite(source.t0) || throw(DomainError(source.t0, "step time must be finite"))
    return source
end

(source::Step)(t) = t >= source.t0 ? source.amplitude : zero(source.amplitude)

"""
$(TYPEDEF)

Apply the voltage ``A\\cos(2πft+φ)``.

$(TYPEDFIELDS)
"""
struct Cosine{T} <: Waveform
    "Terminal label."
    terminal::String
    "Voltage amplitude \\[V\\]."
    amplitude::T
    "Frequency \\[Hz\\]."
    frequency::T
    "Phase \\[rad\\]."
    phase::T

    function Cosine(terminal::AbstractString, amplitude::Real, frequency::Real, phase::Real)
        a, f, p = float.(promote(amplitude, frequency, phase))
        return validate(new{typeof(a)}(String(terminal), a, f, p))
    end
end

"""
$(TYPEDSIGNATURES)

Assign a cosine voltage \\[V\\] at `frequency` \\[Hz\\] with `phase=0` \\[rad\\].
Explicit grids lift the assignment.
"""
function Cosine(terminal, amplitude, frequency; phase=0)
    return parameterize(Cosine, Cosine, (terminal, amplitude, frequency, phase))
end

function validate(source::Cosine)
    isfinite(source.amplitude) || throw(DomainError(source.amplitude, "amplitude must be finite"))
    isfinite(source.frequency) && source.frequency > 0 ||
        throw(DomainError(source.frequency, "frequency must be finite and positive"))
    isfinite(source.phase) || throw(DomainError(source.phase, "phase must be finite"))
    return source
end

(source::Cosine)(t) = source.amplitude * cos(2π * source.frequency * t + source.phase)

"""
$(TYPEDEF)

Apply the voltage ``A(e^{-αt}-e^{-βt})``, with ``β>α>0``.

$(TYPEDFIELDS)
"""
struct DoubleExponential{T} <: Waveform
    "Terminal label."
    terminal::String
    "Voltage amplitude \\[V\\]."
    amplitude::T
    "Slow decay rate \\[1/s\\]."
    alpha::T
    "Fast decay rate \\[1/s\\]."
    beta::T

    function DoubleExponential(terminal::AbstractString, amplitude::Real, alpha::Real, beta::Real)
        a, α, β = float.(promote(amplitude, alpha, beta))
        return validate(new{typeof(a)}(String(terminal), a, α, β))
    end
end

"""
$(TYPEDSIGNATURES)

Assign a double exponential voltage \\[V\\] with decay rates \\[1/s\\], lifting explicit grids.
"""
function DoubleExponential(terminal, amplitude, alpha, beta)
    return parameterize(DoubleExponential, DoubleExponential, (terminal, amplitude, alpha, beta))
end

function validate(source::DoubleExponential)
    isfinite(source.amplitude) || throw(DomainError(source.amplitude, "amplitude must be finite"))
    isfinite(source.alpha) && source.alpha > 0 ||
        throw(DomainError(source.alpha, "alpha must be finite and positive"))
    isfinite(source.beta) && source.beta > source.alpha ||
        throw(DomainError(source.beta, "beta must be finite and greater than alpha"))
    return source
end

(source::DoubleExponential)(t) =
    source.amplitude * (exp(-source.alpha * t) - exp(-source.beta * t))

"""
$(TYPEDEF)

Terminate a whole end with the line's characteristic admittance matrix.
"""
struct Matched end
