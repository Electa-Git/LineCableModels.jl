"""Passive inputs and display projections for the curated study; no solver or transport."""
module StudyCases

export AbstractStudyCase, LineParameters, CorridorImpedance, inputs, preparation_inputs,
    operation, role, case_title, result_series

abstract type AbstractStudyCase end
"""Two identical coaxial cables, with grounded sheaths, in homogeneous earth."""
struct LineParameters <: AbstractStudyCase end
"""Existing OHL/UGC reference network; passive length sensitivity, not statistical UQ."""
struct CorridorImpedance <: AbstractStudyCase end

operation(::LineParameters) = "line.frequency_scan"
operation(::CorridorImpedance) = "impedance.evaluate"
role(::LineParameters) = "parameters"
role(::CorridorImpedance) = "power-flow"
case_title(::LineParameters) = "Frequency-dependent line parameters"
case_title(::CorridorImpedance) = "OHL / UGC impedance sensitivity"

function bounded(value, low, high, label)
    value isa Real && !(value isa Bool) && isfinite(value) && low <= value <= high ||
        throw(ArgumentError("$label must be finite and in [$low, $high]"))
    return Float64(value)
end

function frequency_inputs(low, high, points)
    lo = bounded(low, 1, 1e6, "minimum frequency (Hz)")
    hi = bounded(high, 1, 1e6, "maximum frequency (Hz)")
    lo < hi || throw(ArgumentError("maximum frequency must exceed minimum frequency"))
    n = bounded(points, 2, 200, "display sample count")
    isinteger(n) || throw(ArgumentError("display sample count must be an integer"))
    return lo, hi, Int(n)
end

"""
    inputs(::LineParameters; separation_m=0.5, depth_m=1.0,
        earth_resistivity_ohm_m=100.0, minimum_frequency_hz=1.0,
        maximum_frequency_hz=2500.0, frequency_points=40)

Return complete passive `line.frequency_scan` inputs. Distances are in m, earth
resistivity in Ω m, and frequencies in Hz. The specimen fixes a 10 mm core, 5 mm
insulation and 1 mm sheath at 20 °C; the sheath uses the same conductor material.
The log-spaced display sweep has 2–200 samples. No model is solved or prepared.
The worker remains the validation authority.
"""
function inputs(::LineParameters; separation_m=0.5, depth_m=1.0,
        earth_resistivity_ohm_m=100.0, minimum_frequency_hz=1.0,
        maximum_frequency_hz=2500.0, frequency_points=40)
    lo, hi, n = frequency_inputs(minimum_frequency_hz, maximum_frequency_hz, frequency_points)
    frequencies = exp10.(range(log10(lo), log10(hi); length=n))
    frequencies[1], frequencies[end] = lo, hi
    return Dict{String,Any}(
        "core_radius_m"=>0.01, "insulation_thickness_m"=>0.005, "sheath_thickness_m"=>0.001,
        "conductor_resistivity_ohm_m"=>1.7241e-8, "conductor_relative_permeability"=>1.0,
        "insulation_resistivity_ohm_m"=>1e14, "insulation_relative_permittivity"=>2.3,
        "temperature_celsius"=>20.0, "line_length_m"=>1000.0,
        "separation_m"=>bounded(separation_m, 0.05, 100, "separation (m)"),
        "depth_m"=>bounded(depth_m, 0.05, 100, "burial depth (m)"),
        "earth_resistivity_ohm_m"=>bounded(earth_resistivity_ohm_m, 0.01, 1e6, "earth resistivity (Ω m)"),
        "frequencies_hz"=>frequencies)
end

reference_specification() = Dict{String,Any}("case_id"=>"ohl_ugc_transition_v1", "earth_resistivity_ohm_m"=>100.0)

"""
    inputs(::CorridorImpedance; ugc_share=0.5, corridor_length_m=100000.0,
        length_error_percent=5.0, minimum_frequency_hz=1.0,
        maximum_frequency_hz=2500.0, frequency_points=80)

Return complete `impedance.evaluate` inputs for the prepared reference network.
Length is in m, frequency in Hz, `ugc_share` is dimensionless and length error
is in percent. Share is limited to 0.001–0.999 to match the worker's nonzero-branch
regularization. Changes affect passive corridor elements, not active setpoints.
"""
function inputs(::CorridorImpedance; ugc_share=0.5, corridor_length_m=100000.0,
        length_error_percent=5.0, minimum_frequency_hz=1.0,
        maximum_frequency_hz=2500.0, frequency_points=80)
    lo, hi, n = frequency_inputs(minimum_frequency_hz, maximum_frequency_hz, frequency_points)
    return Dict{String,Any}("specification"=>reference_specification(), "prepared_resource_key"=>"",
        "ugc_share"=>bounded(ugc_share, 0.001, 0.999, "UGC share"),
        "corridor_length_m"=>bounded(corridor_length_m, 1000, 1e6, "corridor length (m)"),
        "length_error_percent"=>bounded(length_error_percent, 0, 50, "length error (%)"),
        "minimum_frequency_hz"=>lo, "maximum_frequency_hz"=>hi, "frequency_points"=>n)
end

"""Return explicit representative warm-up inputs; constructing them allocates no worker."""
preparation_inputs(case::LineParameters) = inputs(case)
preparation_inputs(::CorridorImpedance) = Dict{String,Any}("specification"=>reference_specification())

function finite_vector(value; positive=false)
    value isa AbstractVector && 2 <= length(value) <= 200 || throw(ArgumentError("invalid display series length"))
    all(x -> x isa Real && !(x isa Bool) && isfinite(x) && (!positive || x > 0), value) ||
        throw(ArgumentError("invalid display series value"))
    return Float64.(value)
end

function frequency_vector(value)
    frequencies = finite_vector(value; positive=true)
    all(diff(frequencies) .> 0) || throw(ArgumentError("display frequencies must be strictly increasing"))
    return frequencies
end

"""
    result_series(case, result, quantity)

Validate a bounded worker result and return `(frequency, curves, unit)` for display.
Line quantities are `resistance`, `reactance` (Ω/m), `conductance` and `susceptance`
(S/m): real/imaginary matrix entries, not derived L or C. Curves are self (1,1)
and mutual (1,2). Corridor curves retain worker values in dB re 1 Ω for nominal
and ±length error. No solver, interpolation or uncertainty estimator is added;
malformed/nonfinite results are rejected before plotting.
"""
function result_series(::LineParameters, value, quantity="resistance")
    tensor_key, part, unit = if quantity == "resistance"
        ("series_impedance_ohm_per_m", "real", "Ω/m")
    elseif quantity == "reactance"
        ("series_impedance_ohm_per_m", "imag", "Ω/m")
    elseif quantity == "conductance"
        ("shunt_admittance_s_per_m", "real", "S/m")
    elseif quantity == "susceptance"
        ("shunt_admittance_s_per_m", "imag", "S/m")
    else
        throw(ArgumentError("unknown line display quantity"))
    end
    frequencies = frequency_vector(value["frequencies_hz"])
    for key in ("series_impedance_ohm_per_m", "shunt_admittance_s_per_m")
        tensor = value[key]
        tensor isa AbstractVector && length(tensor) == 2 || throw(ArgumentError("expected two matrix rows"))
        for row in tensor
            row isa AbstractVector && length(row) == 2 || throw(ArgumentError("expected two matrix columns"))
            for column in row
                column isa AbstractVector && length(column) == length(frequencies) || throw(ArgumentError("matrix sample count mismatch"))
                for component in ("real", "imag")
                    finite_vector([sample[component] for sample in column])
                end
            end
        end
    end
    row = value[tensor_key][1]
    curves = [(label="Self (1,1)", values=Float64[sample[part] for sample in row[1]]),
        (label="Mutual (1,2)", values=Float64[sample[part] for sample in row[2]])]
    return (frequency=frequencies, curves=curves, unit=unit)
end

function result_series(::CorridorImpedance, value, quantity="impedance")
    quantity == "impedance" || throw(ArgumentError("unknown corridor display quantity"))
    curves = value["curves"]
    keys = ("base_minus_error", "base", "base_plus_error")
    frequencies = frequency_vector(curves[first(keys)]["frequency_hz"])
    series = map(zip(keys, ("Shorter corridor", "Nominal corridor", "Longer corridor"))) do (key, label)
        curve = curves[key]
        frequency_vector(curve["frequency_hz"]) == frequencies || throw(ArgumentError("curve frequencies differ"))
        values = finite_vector(curve["magnitude_db_ohm"])
        length(values) == length(frequencies) || throw(ArgumentError("curve sample count mismatch"))
        (label=label, values=values)
    end
    return (frequency=frequencies, curves=series, unit="dB re 1 Ω")
end

end
