module ICHQP2026FundamentalsDeck

using Bonito
using Makie
using Markdown
using Observables: async_latest, off
using Printf
using QuadGK
using SpecialFunctions
using WGLMakie
using ..PageAuthoring

const MU0 = 4π * 1.0e-7
const EPSILON0 = 8.8541878128e-12
const PLOT_BACKGROUND = RGBf(0.9, 0.9, 0.9)
const PLOT_TEXT = RGBf(0.16, 0.18, 0.18)
const PLOT_GRID = (RGBf(0.45, 0.47, 0.47), 0.35)
const PRESENTATION_ASSET_DIRECTORY = normpath(joinpath(@__DIR__, "..", "assets"))
const SKIN_EFFECT_IMAGE = joinpath(PRESENTATION_ASSET_DIRECTORY, "skeffect.png")
const EARTH_RETURN_IMAGE = joinpath(PRESENTATION_ASSET_DIRECTORY, "earthreturn.png")
const CABLE_DARK_MODE_IMAGE = joinpath(
    PRESENTATION_ASSET_DIRECTORY,
    "cable_dark_mode.svg"
)
const PREFLIGHT_STATE = preflight_state(
    label = "Fundamentals dashboards not activated"
)
const WARMUP_TASK = Ref{Union{Nothing, Task}}(nothing)
const WARMUP_RESULT = Ref{Any}(nothing)
const WARMUP_LOCK = ReentrantLock()
const SI_PREFIX = Dict(
    -12 => "p",
    -9 => "n",
    -6 => "μ",
    -3 => "m",
    0 => "",
    3 => "k",
    6 => "M",
    9 => "G",
    12 => "T"
)

function engineering(value::Real, unit::AbstractString)
    x = Float64(value)
    isfinite(x) || return "$x $unit"
    iszero(x) && return "0 $unit"
    exponent = clamp(3floor(Int, log10(abs(x)) / 3), -12, 12)
    scaled = x / 10.0^exponent
    return @sprintf("%.4g %s%s", scaled, SI_PREFIX[exponent], unit)
end

function complex_engineering(value::Complex, unit::AbstractString)
    sign = imag(value) >= 0 ? "+" : "−"
    return "$(engineering(real(value), unit)) $sign j $(engineering(abs(imag(value)), unit))"
end

function light_axis(figure_position; attributes...)
    return Axis(
        figure_position;
        backgroundcolor = PLOT_BACKGROUND,
        titlecolor = PLOT_TEXT,
        xlabelcolor = PLOT_TEXT,
        ylabelcolor = PLOT_TEXT,
        xticklabelcolor = PLOT_TEXT,
        yticklabelcolor = PLOT_TEXT,
        xgridcolor = PLOT_GRID,
        ygridcolor = PLOT_GRID,
        attributes...
    )
end

Base.@kwdef struct SkinEffectModel
    conductor_radius::Float64 = 5.0e-3
    insulation_radius::Float64 = 10.0e-3
    conductor_relative_permeability::Float64 = 1.0
    insulation_relative_permeability::Float64 = 1.0
    grid_points::Int = 151
    radial_points::Int = 801
    db_min::Float64 = -120.0
    db_max::Float64 = 80.0
end

struct SkinPlotCache
    x_mm::Vector{Float32}
    y_mm::Vector{Float32}
    radius_m::Vector{Float64}
    radius_mm::Vector{Float32}
    profile_index::Matrix{Int32}
    inside_insulation::BitMatrix
end

function validate(model::SkinEffectModel)
    model.conductor_radius > 0 ||
        throw(ArgumentError("conductor radius must be positive"))
    model.insulation_radius > model.conductor_radius ||
        throw(ArgumentError("insulation radius must exceed conductor radius"))
    model.grid_points >= 32 ||
        throw(ArgumentError("grid_points must be at least 32"))
    model.radial_points >= 128 ||
        throw(ArgumentError("radial_points must be at least 128"))
    model.db_min < model.db_max ||
        throw(ArgumentError("db_min must be smaller than db_max"))
    return model
end

function build_skin_cache(model::SkinEffectModel)
    validate(model)
    radius = model.insulation_radius
    x = collect(range(-radius, radius; length = model.grid_points))
    y = collect(range(-radius, radius; length = model.grid_points))
    radial = collect(range(0.0, radius; length = model.radial_points))
    profile_index = zeros(Int32, model.grid_points, model.grid_points)
    inside = falses(model.grid_points, model.grid_points)

    @inbounds for j in eachindex(y), i in eachindex(x)

        rho = hypot(x[i], y[j])
        if rho <= radius
            inside[i, j] = true
            index = round(Int, rho / radius * (model.radial_points - 1)) + 1
            profile_index[i, j] = Int32(clamp(index, 1, model.radial_points))
        end
    end

    return SkinPlotCache(
        Float32.(1.0e3 .* x),
        Float32.(1.0e3 .* y),
        radial,
        Float32.(1.0e3 .* radial),
        profile_index,
        inside
    )
end

# besselix(ν, z) = exp(-abs(real(z))) Iν(z). Reconstructing the scale only
# in this ratio keeps the MHz copper case finite.
@inline function scaled_i0_ratio(z::ComplexF64, reference::ComplexF64)
    scale = exp(abs(real(z)) - abs(real(reference)))
    return scale * besselix(0.0, z) / besselix(0.0, reference)
end

function internal_impedance_per_length(
        omega::Float64,
        conductivity::Float64,
        permeability::Float64,
        radius::Float64
)
    gamma = sqrt(complex(0.0, omega * permeability * conductivity))
    argument = gamma * radius
    if abs(argument) < 1.0e-3
        return inv(π * radius^2 * conductivity) +
               im * omega * permeability / (8π)
    end
    ratio = besselix(0.0, argument) / besselix(1.0, argument)
    return gamma * ratio / (2π * radius * conductivity)
end

@inline function clipped_db20(value::Float64, db_min::Float64, db_max::Float64)
    floor_amplitude = 10.0^(db_min / 20)
    return clamp(20log10(max(value, floor_amplitude)), db_min, db_max)
end

function solve_skin_state(
        frequency::Float64,
        conductivity::Float64,
        model::SkinEffectModel,
        cache::SkinPlotCache
)
    frequency > 0 || throw(ArgumentError("frequency must be positive"))
    conductivity > 0 || throw(ArgumentError("conductivity must be positive"))

    radius = model.conductor_radius
    conductor_permeability = MU0 * model.conductor_relative_permeability
    insulation_permeability = MU0 * model.insulation_relative_permeability
    omega = 2π * frequency
    gamma = sqrt(complex(0.0, omega * conductor_permeability * conductivity))
    surface_argument = gamma * radius
    impedance = internal_impedance_per_length(
        omega,
        conductivity,
        conductor_permeability,
        radius
    )
    dc_resistance = inv(π * radius^2 * conductivity)
    skin_depth = sqrt(2 / (omega * conductor_permeability * conductivity))
    reference = max(abs(impedance), floatmin(Float64))
    profile_db = Vector{Float32}(undef, length(cache.radius_m))

    @inbounds for index in eachindex(cache.radius_m)
        rho = cache.radius_m[index]
        field_per_ampere = if rho <= radius
            impedance * scaled_i0_ratio(
                ComplexF64(gamma * rho),
                ComplexF64(surface_argument)
            )
        else
            impedance + im * omega * insulation_permeability *
                        log(rho / radius) / (2π)
        end
        profile_db[index] = Float32(clipped_db20(
            abs(field_per_ampere) / reference,
            model.db_min,
            model.db_max
        ))
    end

    field_db = fill(Float32(NaN), model.grid_points, model.grid_points)
    @inbounds for index in eachindex(field_db)
        if cache.inside_insulation[index]
            field_db[index] = profile_db[Int(cache.profile_index[index])]
        end
    end

    return (;
        frequency,
        conductivity,
        skin_depth,
        radius_over_skin_depth = radius / skin_depth,
        impedance,
        dc_resistance,
        ac_resistance = real(impedance),
        internal_inductance = imag(impedance) / omega,
        ac_over_dc = real(impedance) / dc_resistance,
        center_db = Float64(first(profile_db)),
        profile_db,
        field_db
    )
end

function make_skin_figure(
        session::Session,
        model::SkinEffectModel,
        cache::SkinPlotCache,
        state::Observable,
        color_range::Observable
)
    conductor_radius_mm = Float32(1.0e3 * model.conductor_radius)
    insulation_radius_mm = Float32(1.0e3 * model.insulation_radius)
    field_db = map(value -> value.field_db, session, state)
    profile_db = map(value -> value.profile_db, session, state)
    title = map(session, state) do value
        "Axial field · $(engineering(value.frequency, "Hz")) · σ = $(engineering(value.conductivity, "S/m"))"
    end
    skin_position = map(
        value -> Float32[(model.conductor_radius - value.skin_depth) * 1.0e3],
        session,
        state
    )
    skin_visible = map(
        value -> value.skin_depth < model.conductor_radius,
        session,
        state
    )

    figure = Figure(size = (1100, 610), backgroundcolor = PLOT_BACKGROUND)
    field_axis = light_axis(
        figure[1, 1];
        title,
        xlabel = "lateral coordinate (mm)",
        ylabel = "vertical coordinate (mm)",
        aspect = DataAspect()
    )
    field_plot = heatmap!(
        field_axis,
        cache.x_mm,
        cache.y_mm,
        field_db;
        colormap = :inferno,
        colorrange = color_range,
        nan_color = :transparent,
        interpolate = true
    )
    theta = range(0.0, 2π; length = 361)
    lines!(
        field_axis,
        conductor_radius_mm .* cos.(theta),
        conductor_radius_mm .* sin.(theta);
        color = :white,
        linewidth = 2.2
    )
    lines!(
        field_axis,
        insulation_radius_mm .* cos.(theta),
        insulation_radius_mm .* sin.(theta);
        color = (:white, 0.7),
        linewidth = 1.6,
        linestyle = :dash
    )
    limits!(
        field_axis,
        -insulation_radius_mm,
        insulation_radius_mm,
        -insulation_radius_mm,
        insulation_radius_mm
    )
    Colorbar(
        figure[1, 2],
        field_plot;
        label = "20 log₁₀ |Eₓ(r) / Eₓ(a)| (dB)",
        labelcolor = PLOT_TEXT,
        ticklabelcolor = PLOT_TEXT,
        width = 16
    )

    profile_axis = light_axis(
        figure[1, 3];
        title = "Radial field profile",
        xlabel = "radius (mm)",
        ylabel = "relative amplitude (dB)",
        limits = map(session, color_range) do range
            (0.0f0, insulation_radius_mm, Float32(first(range)), Float32(last(range)))
        end
    )
    profile_plot = lines!(
        profile_axis,
        cache.radius_mm,
        profile_db;
        linewidth = 3.0,
        color = RGBf(0.19, 0.47, 0.91),
        label = "|Eₓ(r) / Eₓ(a)|"
    )
    hlines!(profile_axis, Float32[0]; color = (:black, 0.35), linewidth = 1)
    vlines!(profile_axis, Float32[conductor_radius_mm]; color = :black, linewidth = 1.8)
    skin_marker = vlines!(
        profile_axis,
        skin_position;
        color = RGBf(0.93, 0.45, 0.08),
        linewidth = 2.0,
        linestyle = :dash,
        visible = skin_visible,
        label = "one skin depth"
    )
    axislegend(
        profile_axis;
        position = :lb,
        backgroundcolor = (:white, 0.85),
        labelcolor = PLOT_TEXT,
        framevisible = false
    )
    DataInspector(figure)
    return (;
        figure,
        field_axis,
        profile_axis,
        field_plot,
        profile_plot,
        skin_marker
    )
end

Base.@kwdef struct EarthReturnModel
    conductor_radius::Float64 = 12.0e-3
    conductor_height::Float64 = 10.0
    conductor_spacing::Float64 = 1.0
    earth_relative_permeability::Float64 = 1.0
    local_y_extent::Float64 = 4.0
    local_z_min::Float64 = -4.0
    local_z_max::Float64 = 14.0
    local_ny::Int = 91
    local_nz::Int = 91
    earth_normalized_extent::Float64 = 4.0
    earth_ny::Int = 91
    earth_nz::Int = 81
    spectrum_points::Int = 256
    spectrum_decay_orders::Float64 = 45.0
    field_db_min::Float64 = -120.0
    current_db_min::Float64 = -120.0
end

struct EarthFieldCache
    y_local::Vector{Float64}
    z_local::Vector{Float64}
    lambda::Vector{Float64}
    weights::Vector{Float64}
end

function validate(model::EarthReturnModel)
    model.conductor_radius > 0 ||
        throw(ArgumentError("conductor radius must be positive"))
    model.conductor_height > model.conductor_radius ||
        throw(ArgumentError("conductor height must exceed conductor radius"))
    model.conductor_spacing > 2model.conductor_radius ||
        throw(ArgumentError("conductor spacing must exceed conductor diameter"))
    model.local_ny >= 32 && model.local_nz >= 32 ||
        throw(ArgumentError("local field grids must have at least 32 points"))
    model.earth_ny >= 32 && model.earth_nz >= 32 ||
        throw(ArgumentError("earth field grids must have at least 32 points"))
    model.spectrum_points >= 128 ||
        throw(ArgumentError("spectrum_points must be at least 128"))
    return model
end

function build_earth_cache(model::EarthReturnModel)
    validate(model)
    y = collect(range(-model.local_y_extent, model.local_y_extent; length = model.local_ny))
    z = collect(range(model.local_z_min, model.local_z_max; length = model.local_nz))
    lambda_max = model.spectrum_decay_orders / model.conductor_height
    lambda_range = range(0.0, lambda_max; length = model.spectrum_points)
    lambda = collect(lambda_range)
    weights = fill(step(lambda_range), model.spectrum_points)
    firstindex(weights) == lastindex(weights) || begin
        first_index = firstindex(weights)
        last_index = lastindex(weights)
        weights[first_index] *= 0.5
        weights[last_index] *= 0.5
    end
    return EarthFieldCache(y, z, lambda, weights)
end

@inline function passive_sqrt(value::ComplexF64)
    root = sqrt(value)
    if real(root) < 0 || (iszero(real(root)) && imag(root) < 0)
        return -root
    end
    return root
end

@inline function db20(value::Real, floor_db::Real)
    floor_amplitude = 10.0^(floor_db / 20)
    return max(20log10(max(Float64(value), floor_amplitude)), floor_db)
end

@inline function source_positions(model::EarthReturnModel)
    half_spacing = model.conductor_spacing / 2
    return (-half_spacing, half_spacing)
end

function conductor_currents(current::Float64, return_ratio::Float64)
    ratio = clamp(return_ratio, 0.0, 1.0)
    return (ComplexF64(current), ComplexF64(-ratio * current))
end

function medium_state(
        frequency::Float64,
        resistivity::Float64,
        relative_permittivity::Float64,
        model::EarthReturnModel
)
    omega = 2π * frequency
    conductivity = inv(resistivity)
    permittivity = EPSILON0 * relative_permittivity
    permeability = MU0 * model.earth_relative_permeability
    gamma_squared = ComplexF64(
        im * omega * permeability *
        (conductivity + im * omega * permittivity)
    )
    skin_depth = sqrt(2resistivity / (omega * permeability))
    propagation = passive_sqrt(gamma_squared)
    attenuation_length = real(propagation) > 0 ? inv(real(propagation)) : Inf
    conduction_ratio = conductivity / (omega * permittivity)
    return (;
        omega,
        conductivity,
        permittivity,
        permeability,
        gamma_squared,
        skin_depth,
        attenuation_length,
        conduction_ratio
    )
end

function spectral_terms(cache::EarthFieldCache, medium)
    propagation = Vector{ComplexF64}(undef, length(cache.lambda))
    reflection = similar(propagation)
    transmission_denominator = similar(propagation)
    @inbounds for index in eachindex(cache.lambda)
        lambda = cache.lambda[index]
        root = passive_sqrt(ComplexF64(lambda^2) + medium.gamma_squared)
        propagation[index] = root
        reflection[index] = iszero(lambda) ? -1.0 + 0im :
                            (lambda - root) / (lambda + root)
        transmission_denominator[index] = lambda + root
    end
    return propagation, reflection, transmission_denominator
end

function external_field_map(
        current::Float64,
        return_ratio::Float64,
        medium,
        model::EarthReturnModel,
        cache::EarthFieldCache
)
    source_y = source_positions(model)
    source_current = conductor_currents(current, return_ratio)
    propagation, reflection, transmission_denominator = spectral_terms(cache, medium)
    magnitude = Matrix{Float32}(undef, length(cache.y_local), length(cache.z_local))
    height = model.conductor_height
    radius = model.conductor_radius

    @inbounds for z_index in eachindex(cache.z_local), y_index in eachindex(cache.y_local)

        y = cache.y_local[y_index]
        z = cache.z_local[z_index]
        hy = 0.0 + 0im
        hz = 0.0 + 0im
        for source in 1:2
            source_current[source] == 0 && continue
            dy = y - source_y[source]
            if z >= 0
                dz = z - height
                radius_squared = dy^2 + dz^2
                if radius_squared > radius^2
                    hy += -source_current[source] * dz / (2π * radius_squared)
                    hz += source_current[source] * dy / (2π * radius_squared)
                end
                reflected_hy = 0.0 + 0im
                reflected_hz = 0.0 + 0im
                for index in eachindex(cache.lambda)
                    lambda = cache.lambda[index]
                    exponential = exp(-lambda * (z + height))
                    weight = cache.weights[index]
                    reflected_hy += weight * reflection[index] * exponential *
                                    cos(lambda * dy)
                    reflected_hz += weight * reflection[index] * exponential *
                                    sin(lambda * dy)
                end
                hy += -source_current[source] * reflected_hy / (2π)
                hz += source_current[source] * reflected_hz / (2π)
            else
                transmitted_hy = 0.0 + 0im
                transmitted_hz = 0.0 + 0im
                for index in eachindex(cache.lambda)
                    lambda = cache.lambda[index]
                    exponential = exp(-lambda * height + propagation[index] * z)
                    weight = cache.weights[index]
                    transmitted_hy += weight *
                                      (propagation[index] /
                                       transmission_denominator[index]) *
                                      exponential * cos(lambda * dy)
                    transmitted_hz += weight *
                                      (lambda / transmission_denominator[index]) *
                                      exponential * sin(lambda * dy)
                end
                hy += source_current[source] * transmitted_hy / π
                hz += source_current[source] * transmitted_hz / π
            end
        end
        magnitude[y_index, z_index] = Float32(hypot(abs(hy), abs(hz)))
    end

    for source in 1:2
        source_current[source] == 0 && source == 2 && continue
        @inbounds for z_index in eachindex(cache.z_local),
            y_index in eachindex(cache.y_local)

            if hypot(
                cache.y_local[y_index] - source_y[source],
                cache.z_local[z_index] - height
            ) <= radius
                magnitude[y_index, z_index] = NaN32
            end
        end
    end

    maximum_field = maximum(filter(isfinite, vec(magnitude)))
    field_db = map(magnitude) do value
        isfinite(value) ?
        Float32(db20(value / maximum_field, model.field_db_min)) : NaN32
    end
    return field_db, maximum_field
end

function earth_current_map(
        current::Float64,
        return_ratio::Float64,
        medium,
        model::EarthReturnModel,
        cache::EarthFieldCache
)
    source_y = source_positions(model)
    source_current = conductor_currents(current, return_ratio)
    propagation, _, transmission_denominator = spectral_terms(cache, medium)
    height = model.conductor_height
    plot_skin_depth = clamp(medium.skin_depth, 0.05, 5.0e3)
    normalized_y = collect(range(
        -model.earth_normalized_extent,
        model.earth_normalized_extent;
        length = model.earth_ny
    ))
    normalized_z = collect(range(
        0.0,
        model.earth_normalized_extent;
        length = model.earth_nz
    ))
    y = plot_skin_depth .* normalized_y
    z = -plot_skin_depth .* normalized_z
    current_density = Matrix{ComplexF64}(undef, length(y), length(z))

    @inbounds for z_index in eachindex(z), y_index in eachindex(y)

        vector_potential = 0.0 + 0im
        for source in 1:2
            source_current[source] == 0 && continue
            dy = y[y_index] - source_y[source]
            integral = 0.0 + 0im
            for index in eachindex(cache.lambda)
                lambda = cache.lambda[index]
                integral += cache.weights[index] *
                            exp(-lambda * height + propagation[index] * z[z_index]) *
                            cos(lambda * dy) / transmission_denominator[index]
            end
            vector_potential += medium.permeability * source_current[source] * integral / π
        end
        current_density[y_index, z_index] = -im * medium.omega * medium.conductivity *
                                            vector_potential
    end

    maximum_current = maximum(abs, current_density)
    current_db = map(current_density) do value
        Float32(db20(abs(value) / max(maximum_current, eps()), model.current_db_min))
    end
    return normalized_y, normalized_z, current_db, maximum_current
end

function generalized_earth_correction(
        separation::Float64,
        first_height::Float64,
        second_height::Float64,
        medium
)
    integrand(lambda) = begin
        root = passive_sqrt(ComplexF64(lambda^2) + medium.gamma_squared)
        exp(-(first_height + second_height) * lambda) *
        cos(separation * lambda) / (lambda + root)
    end
    integral, _ = quadgk(integrand, 0.0, Inf; rtol = 2.0e-7)
    return im * medium.omega * MU0 * integral / π
end

function earth_impedance_state(medium, model::EarthReturnModel)
    height = model.conductor_height
    spacing = model.conductor_spacing
    radius = model.conductor_radius
    earth_self = generalized_earth_correction(0.0, height, height, medium)
    earth_mutual = generalized_earth_correction(spacing, height, height, medium)
    geometric_self = im * medium.omega * MU0 / (2π) * log(2height / radius)
    image_distance = hypot(spacing, 2height)
    geometric_mutual = im * medium.omega * MU0 / (2π) *
                       log(image_distance / spacing)
    external_self = geometric_self + earth_self
    external_mutual = geometric_mutual + earth_mutual
    return (;
        earth_self,
        earth_mutual,
        external_self,
        external_mutual,
        differential = external_self - external_mutual
    )
end

function solve_earth_state(
        frequency::Float64,
        current::Float64,
        resistivity::Float64,
        relative_permittivity::Float64,
        return_ratio::Float64,
        model::EarthReturnModel,
        cache::EarthFieldCache
)
    frequency > 0 || throw(ArgumentError("frequency must be positive"))
    current > 0 || throw(ArgumentError("current must be positive"))
    resistivity > 0 || throw(ArgumentError("earth resistivity must be positive"))
    relative_permittivity >= 1 ||
        throw(ArgumentError("relative permittivity must be at least one"))
    medium = medium_state(frequency, resistivity, relative_permittivity, model)
    field_db, maximum_field = external_field_map(
        current,
        return_ratio,
        medium,
        model,
        cache
    )
    normalized_y, normalized_z, current_db, maximum_current = earth_current_map(
        current,
        return_ratio,
        medium,
        model,
        cache
    )
    impedance = earth_impedance_state(medium, model)
    return (;
        frequency,
        current,
        resistivity,
        relative_permittivity,
        return_ratio,
        medium...,
        field_db,
        maximum_field,
        normalized_y,
        normalized_z,
        current_db,
        maximum_current,
        impedance...
    )
end

function make_earth_figure(
        session::Session,
        model::EarthReturnModel,
        cache::EarthFieldCache,
        state::Observable,
        field_color_range::Observable,
        current_color_range::Observable
)
    field_db = map(value -> value.field_db, session, state)
    current_db = map(value -> value.current_db, session, state)
    normalized_y = map(value -> Float32.(value.normalized_y), session, state)
    normalized_z = map(value -> Float32.(value.normalized_z), session, state)
    field_title = map(session, state) do value
        value.return_ratio > 0.5 ?
        "External |H| · two-wire ±I excitation" :
        "External |H| · single-wire excitation"
    end

    figure = Figure(size = (1120, 610), backgroundcolor = PLOT_BACKGROUND)
    field_axis = light_axis(
        figure[1, 1];
        title = field_title,
        xlabel = "lateral position y (m)",
        ylabel = "height z (m)",
        aspect = DataAspect()
    )
    field_plot = heatmap!(
        field_axis,
        Float32.(cache.y_local),
        Float32.(cache.z_local),
        field_db;
        colorrange = field_color_range,
        colormap = :turbo
    )
    hlines!(field_axis, [0.0f0]; linewidth = 2, color = :black)
    Colorbar(
        figure[1, 2],
        field_plot;
        label = "20 log₁₀(|H| / max |H|) (dB)",
        labelcolor = PLOT_TEXT,
        ticklabelcolor = PLOT_TEXT,
        width = 16
    )
    first_y, second_y = source_positions(model)
    scatter!(
        field_axis,
        Float32[first_y, second_y],
        Float32[model.conductor_height, model.conductor_height];
        marker = :circle,
        markersize = 14,
        color = RGBf(0.19, 0.47, 0.91)
    )

    current_axis = light_axis(
        figure[1, 3];
        title = "Induced earth conduction current |Jₓ|",
        xlabel = "lateral coordinate y / δ",
        ylabel = "depth / δ"
    )
    current_plot = heatmap!(
        current_axis,
        normalized_y,
        normalized_z,
        current_db;
        colorrange = current_color_range,
        colormap = :inferno
    )
    current_axis.yreversed = true
    Colorbar(
        figure[1, 4],
        current_plot;
        label = "20 log₁₀(|Jₓ| / max |Jₓ|) (dB)",
        labelcolor = PLOT_TEXT,
        ticklabelcolor = PLOT_TEXT,
        width = 16
    )
    DataInspector(figure)
    return (;
        figure,
        field_axis,
        current_axis,
        field_plot,
        current_plot
    )
end

function setup(session::Session)
    skin_model = validate(SkinEffectModel())
    skin_cache = build_skin_cache(skin_model)
    frequency_values = collect(0.0:0.025:7.0)
    conductivity_values = collect(4.0:0.02:8.0)
    conductivity_default = conductivity_values[
        argmin(abs.(conductivity_values .- log10(5.8e7)))
    ]
    skin_frequency_slider = Bonito.Slider(
        frequency_values;
        value = 3.0,
        id = "skin-frequency-input",
        ariaLabel = "Skin effect frequency on a logarithmic scale"
    )
    conductivity_slider = Bonito.Slider(
        conductivity_values;
        value = conductivity_default,
        id = "conductor-conductivity-input",
        ariaLabel = "Conductor conductivity on a logarithmic scale"
    )
    skin_db_floor_slider = Bonito.Slider(
        collect(-120.0:5.0:-10.0);
        value = -80.0,
        id = "skin-db-floor-input",
        ariaLabel = "Skin effect plot lower decibel bound"
    )
    skin_db_ceiling_slider = Bonito.Slider(
        collect(0.0:5.0:80.0);
        value = 60.0,
        id = "skin-db-ceiling-input",
        ariaLabel = "Skin effect plot upper decibel bound"
    )
    skin_color_range = map(
        (floor, ceiling) -> (Float32(floor), Float32(ceiling)),
        session,
        skin_db_floor_slider.value,
        skin_db_ceiling_slider.value
    )
    skin_reset_button = Bonito.Button(
        "Reset plot view";
        style = nothing,
        id = "skin-reset-view-control",
        class = "lc-action-button",
        type = "button",
        ariaLabel = "Reset both skin effect plot axes"
    )
    skin_inputs = map(
        (log_frequency, log_conductivity) -> (;
            frequency = 10.0^Float64(log_frequency),
            conductivity = 10.0^Float64(log_conductivity)
        ),
        session,
        skin_frequency_slider.value,
        conductivity_slider.value
    )
    skin_state = Observable(solve_skin_state(
        skin_inputs[].frequency,
        skin_inputs[].conductivity,
        skin_model,
        skin_cache
    ))
    latest_skin_inputs = async_latest(skin_inputs, 1)
    skin_observer = on(latest_skin_inputs) do values
        skin_state[] = solve_skin_state(
            values.frequency,
            values.conductivity,
            skin_model,
            skin_cache
        )
        return nothing
    end
    skin_plot = make_skin_figure(
        session,
        skin_model,
        skin_cache,
        skin_state,
        skin_color_range
    )
    skin_reset_observer = on(skin_reset_button.value) do clicked
        clicked || return nothing
        reset_limits!(skin_plot.field_axis)
        reset_limits!(skin_plot.profile_axis)
        return nothing
    end

    earth_model = validate(EarthReturnModel())
    earth_cache = build_earth_cache(earth_model)
    earth_frequency_values = collect(-1.0:0.025:6.0)
    earth_current_values = collect(0.0:0.025:4.0)
    earth_frequency_default = earth_frequency_values[
        argmin(abs.(earth_frequency_values .- log10(50.0)))
    ]
    earth_current_default = earth_current_values[
        argmin(abs.(earth_current_values .- log10(500.0)))
    ]
    earth_frequency_slider = Bonito.Slider(
        earth_frequency_values;
        value = earth_frequency_default,
        id = "earth-frequency-input",
        ariaLabel = "Earth return frequency on a logarithmic scale"
    )
    current_slider = Bonito.Slider(
        earth_current_values;
        value = earth_current_default,
        id = "earth-current-input",
        ariaLabel = "Conductor current on a logarithmic scale"
    )
    resistivity_slider = Bonito.Slider(
        collect(-1.0:0.025:4.0);
        value = 2.0,
        id = "earth-resistivity-input",
        ariaLabel = "Earth resistivity on a logarithmic scale"
    )
    permittivity_slider = Bonito.Slider(
        collect(0.0:0.02:2.0);
        value = 1.0,
        id = "earth-permittivity-input",
        ariaLabel = "Earth relative permittivity on a logarithmic scale"
    )
    return_ratio_slider = Bonito.Slider(
        collect(0.0:0.01:1.0);
        value = 0.0,
        id = "earth-return-ratio-input",
        ariaLabel = "Second conductor return current ratio"
    )
    earth_field_db_floor_slider = Bonito.Slider(
        collect(-120.0:5.0:-20.0);
        value = -70.0,
        id = "earth-field-db-floor-input",
        ariaLabel = "External magnetic field lower decibel bound"
    )
    earth_current_db_floor_slider = Bonito.Slider(
        collect(-120.0:5.0:-20.0);
        value = -80.0,
        id = "earth-current-db-floor-input",
        ariaLabel = "Earth current density lower decibel bound"
    )
    earth_field_color_range = map(
        floor -> (Float32(floor), 0.0f0),
        session,
        earth_field_db_floor_slider.value
    )
    earth_current_color_range = map(
        floor -> (Float32(floor), 0.0f0),
        session,
        earth_current_db_floor_slider.value
    )
    earth_reset_button = Bonito.Button(
        "Reset plot view";
        style = nothing,
        id = "earth-reset-view-control",
        class = "lc-action-button",
        type = "button",
        ariaLabel = "Reset both earth return plot axes"
    )
    earth_inputs = map(
        (
            log_frequency,
            log_current,
            log_resistivity,
            log_permittivity,
            return_ratio
        ) -> (;
            frequency = 10.0^Float64(log_frequency),
            current = 10.0^Float64(log_current),
            resistivity = 10.0^Float64(log_resistivity),
            relative_permittivity = 10.0^Float64(log_permittivity),
            return_ratio = Float64(return_ratio)
        ),
        session,
        earth_frequency_slider.value,
        current_slider.value,
        resistivity_slider.value,
        permittivity_slider.value,
        return_ratio_slider.value
    )
    earth_state = Observable(solve_earth_state(
        earth_inputs[].frequency,
        earth_inputs[].current,
        earth_inputs[].resistivity,
        earth_inputs[].relative_permittivity,
        earth_inputs[].return_ratio,
        earth_model,
        earth_cache
    ))
    latest_earth_inputs = async_latest(earth_inputs, 1)
    earth_observer = on(latest_earth_inputs) do values
        earth_state[] = solve_earth_state(
            values.frequency,
            values.current,
            values.resistivity,
            values.relative_permittivity,
            values.return_ratio,
            earth_model,
            earth_cache
        )
        return nothing
    end
    earth_plot = make_earth_figure(
        session,
        earth_model,
        earth_cache,
        earth_state,
        earth_field_color_range,
        earth_current_color_range
    )
    earth_reset_observer = on(earth_reset_button.value) do clicked
        clicked || return nothing
        reset_limits!(earth_plot.field_axis)
        reset_limits!(earth_plot.current_axis)
        return nothing
    end

    on(session.on_close) do _
        off(skin_observer)
        off(earth_observer)
        off(skin_reset_observer)
        off(earth_reset_observer)
        return nothing
    end
    return (;
        skin_model,
        skin_cache,
        skin_inputs,
        skin_state,
        skin_plot,
        skin_frequency_slider,
        conductivity_slider,
        skin_db_floor_slider,
        skin_db_ceiling_slider,
        skin_color_range,
        skin_reset_button,
        skin_reset_observer,
        earth_model,
        earth_cache,
        earth_inputs,
        earth_state,
        earth_plot,
        earth_frequency_slider,
        current_slider,
        resistivity_slider,
        permittivity_slider,
        return_ratio_slider,
        earth_field_db_floor_slider,
        earth_current_db_floor_slider,
        earth_field_color_range,
        earth_current_color_range,
        earth_reset_button,
        earth_reset_observer
    )
end

function sources_of_uncertainties_page(::Session, ::Any)
    origins = webpart(
        prose(md"""
        ## Internal and external origins

        - **Geometrical and material properties**
        - **Real field data**, including soil resistivity and the actual conductor layout
        - **Presence of interferences**
        - **Modelling procedure**, from parameter calculation to the EMT representation

        ## Modeling complexity

        - Composite structures: solid, tubular or stranded cores, semiconductors, screens, armors, sheaths, tapes, and water-blocking materials etc

        - Complex interaction with the power system and surrounding media (earth return, thermal effects)

        - High-fidelity models (e.g. FEM) are computationally expensive and/or do not interface well with EMT-type simulation software (e.g. PSCAD)

        - Simplified models (e.g. closed-form approximations) facilitate computations, but introduce additional uncertainties
        """);
        kind = :panel
    )
    skin_effect = webpart(
        local_image(
            SKIN_EFFECT_IMAGE;
            alt = "Cross-sectional visualization of conductor skin effect",
            class = "lc-content-illustration"
        );
        kind = :panel,
        title = "Internal origin · skin effect"
    )
    earth_return = webpart(
        local_image(
            EARTH_RETURN_IMAGE;
            alt = "Visualization of current returning through the earth",
            class = "lc-content-illustration"
        );
        kind = :panel,
        title = "External origin · earth return"
    )
    canvas = webgrid(
        [:origins :skin_effect; :origins :earth_return];
        columns = "minmax(24rem, 1fr) minmax(0, 1fr)",
        compact_columns = "minmax(0, 1fr)",
        rows = "repeat(2, minmax(0, 1fr))",
        stack_rows = "auto minmax(16rem, 1fr) minmax(16rem, 1fr)",
        height = "100%",
        class = "lc-master-layout lc-layout-uncertainty-sources",
        origins,
        skin_effect,
        earth_return
    )
    return (body = slide("Sources of uncertainties in cable parameters", canvas),)
end

function internal_impedances_page(session::Session, state)
    frequency = map(value -> engineering(value.frequency, "Hz"), session, state.skin_inputs)
    conductivity = map(
        value -> engineering(value.conductivity, "S/m"),
        session,
        state.skin_inputs
    )
    skin_depth = map(
        value -> engineering(value.skin_depth, "m"),
        session,
        state.skin_state
    )
    penetration = map(
        value -> @sprintf("%.4g", value.radius_over_skin_depth),
        session,
        state.skin_state
    )
    resistance_ratio = map(
        value -> @sprintf("%.4g", value.ac_over_dc),
        session,
        state.skin_state
    )
    impedance = map(
        value -> complex_engineering(value.impedance, "Ω/m"),
        session,
        state.skin_state
    )
    inductance = map(
        value -> engineering(value.internal_inductance, "H/m"),
        session,
        state.skin_state
    )
    db_floor = map(
        value -> @sprintf("%.0f dB", value),
        session,
        state.skin_db_floor_slider.value
    )
    db_ceiling = map(
        value -> @sprintf("%.0f dB", value),
        session,
        state.skin_db_ceiling_slider.value
    )

    controls = webpart(
        control(
            "Frequency",
            state.skin_frequency_slider;
            value = frequency,
            id = "skin-frequency-control"
        ),
        control(
            "Conductor conductivity",
            state.conductivity_slider;
            value = conductivity,
            id = "conductor-conductivity-control"
        ),
        control(
            "Display floor",
            state.skin_db_floor_slider;
            value = db_floor,
            id = "skin-db-floor-control"
        ),
        control(
            "Display ceiling",
            state.skin_db_ceiling_slider;
            value = db_ceiling,
            id = "skin-db-ceiling-control"
        ),
        value_list(
            "Skin depth δ" => skin_depth,
            "Penetration a/δ" => penetration,
            "Rac / Rdc" => resistance_ratio,
            "Zint" => impedance,
            "Lint" => inductance
        ),
        state.skin_reset_button,
        diagnostic(
            "I₀(γr) / I₀(γa)";
            suffix = " · scaled complex Bessel evaluation"
        );
        kind = :controls,
        class = "lc-fundamentals-controls"
    )
    plot = webpart(
        wgl_figure(state.skin_plot.figure);
        kind = :plot,
        title = "Cross-sectional field and radial profile",
        meta = "adjustable dB scale",
        id = "skin-effect-plot",
        body_class = "lc-transition-plot-host"
    )
    equations = webpart(
        prose(
            md"""
      ``Z_{\mathrm{int}}=\frac{\gamma_c I_0(\gamma_c a)}{2\pi a\sigma_c I_1(\gamma_c a)}``

      ``\gamma_c=\sqrt{j\omega\mu_c\sigma_c}``

      ``L_{\mathrm{int}}=\frac{\operatorname{Im}(Z_{\mathrm{int}})}{\omega}``
      """;
            class = "lc-equation");
        kind = :panel,
        title = "Internal-impedance model",
        class = "lc-equation-panel"
    )
    canvas = layout_controls_plot(
        controls,
        plot;
        footer = equations,
        footer_placement = :under_plot,
        footer_rows = "minmax(0, 1fr) minmax(11.5rem, 0.34fr)",
        height = "100%",
        class = "lc-interactive-grid lc-transition-grid lc-fundamentals-dashboard-grid"
    )
    return (body = slide(
        "Internal impedances",
        canvas;
        lede = md"""
        Drag frequency or conductivity to inspect cylindrical current crowding. Inside the metal, ``J_x = \sigma_c E_x``; the insulation continuation stops at its physical boundary rather than inventing a far external field.
        """
    ),)
end

function earth_return_page(session::Session, state)
    frequency = map(value -> engineering(value.frequency, "Hz"), session, state.earth_inputs)
    current = map(value -> engineering(value.current, "A"), session, state.earth_inputs)
    resistivity = map(
        value -> engineering(value.resistivity, "Ω·m"),
        session,
        state.earth_inputs
    )
    permittivity = map(
        value -> @sprintf("%.3g", value.relative_permittivity),
        session,
        state.earth_inputs
    )
    return_ratio = map(
        value -> @sprintf("%.2f", value.return_ratio),
        session,
        state.earth_inputs
    )
    skin_depth = map(
        value -> engineering(value.skin_depth, "m"),
        session,
        state.earth_state
    )
    conduction_ratio = map(
        value -> @sprintf("%.3g", value.conduction_ratio),
        session,
        state.earth_state
    )
    self_correction = map(
        value -> complex_engineering(value.earth_self, "Ω/m"),
        session,
        state.earth_state
    )
    mutual_correction = map(
        value -> complex_engineering(value.earth_mutual, "Ω/m"),
        session,
        state.earth_state
    )
    differential = map(
        value -> complex_engineering(value.differential, "Ω/m"),
        session,
        state.earth_state
    )
    peak_current = map(
        value -> engineering(value.maximum_current, "A/m²"),
        session,
        state.earth_state
    )
    field_db_floor = map(
        value -> @sprintf("%.0f dB", value),
        session,
        state.earth_field_db_floor_slider.value
    )
    current_db_floor = map(
        value -> @sprintf("%.0f dB", value),
        session,
        state.earth_current_db_floor_slider.value
    )

    controls = webpart(
        control(
            "Frequency",
            state.earth_frequency_slider;
            value = frequency,
            id = "earth-frequency-control"
        ),
        control(
            "Current amplitude",
            state.current_slider;
            value = current,
            id = "earth-current-control"
        ),
        control(
            "Earth resistivity",
            state.resistivity_slider;
            value = resistivity,
            id = "earth-resistivity-control"
        ),
        control(
            "Earth εᵣ",
            state.permittivity_slider;
            value = permittivity,
            id = "earth-permittivity-control"
        ),
        control(
            "Second-conductor return ratio α",
            state.return_ratio_slider;
            value = return_ratio,
            id = "earth-return-ratio-control"
        ),
        control(
            "External |H| display floor",
            state.earth_field_db_floor_slider;
            value = field_db_floor,
            id = "earth-field-db-floor-control"
        ),
        control(
            "Earth |Jₓ| display floor",
            state.earth_current_db_floor_slider;
            value = current_db_floor,
            id = "earth-current-db-floor-control"
        ),
        value_list(
            "Earth skin depth" => skin_depth,
            "σ / (ωε)" => conduction_ratio,
            "Earth self Zᵍ₁₁" => self_correction,
            "Earth mutual Zᵍ₁₂" => mutual_correction,
            "Zext,₁₁ − Zext,₁₂" => differential,
            "Peak induced |Jₓ|" => peak_current
        ),
        state.earth_reset_button;
        kind = :controls,
        class = "lc-fundamentals-controls"
    )
    plot = webpart(
        wgl_figure(state.earth_plot.figure);
        kind = :plot,
        title = "Lossy half-space field solution",
        meta = "independent adjustable dB floors",
        id = "earth-return-plot",
        body_class = "lc-transition-plot-host"
    )
    equations = webpart(
        prose(
            md"""
      ``Z^g_{ij}=\frac{j\omega\mu_0}{\pi}\int_0^\infty\frac{e^{-(h_i+h_j)\lambda}\cos(y_{ij}\lambda)}{\lambda+\sqrt{\lambda^2+\gamma_g^2}}\,\mathrm{d}\lambda``

      ``Z^{\mathrm{ext}}_{ij}=Z^{\mathrm{perf}}_{ij}+Z^g_{ij}.``
      """;
            class = "lc-equation");
        kind = :panel,
        title = "Generalized earth-return impedance",
        class = "lc-equation-panel"
    )
    canvas = layout_controls_plot(
        controls,
        plot;
        footer = equations,
        footer_placement = :under_plot,
        footer_rows = "minmax(0, 1fr) minmax(11.5rem, 0.34fr)",
        height = "100%",
        class = "lc-interactive-grid lc-transition-grid lc-fundamentals-dashboard-grid"
    )
    return (body = slide(
        "Earth return",
        canvas;
        lede = md"""
        Explore the external magnetic field and induced earth-current density for conductors above a homogeneous lossy half-space. The spectral field solution and  earth-return impedance retain both ``\sigma`` and ``\varepsilon``.
        """
    ),)
end

function uncertainty_quantification_page(::Session, ::Any)
    cable_model = webpart(
        local_image(
            CABLE_DARK_MODE_IMAGE;
            alt = "Cable model with uncertain geometrical and material parameters",
            class = "lc-content-illustration"
        );
        kind = :panel,
        title = "Uncertain cable model"
    )
    propagation = webpart(
        prose(md"""
        **Addition and subtraction**

        ``\bar{z}=\bar{x}\pm\bar{y}=(x\pm y)\pm\sqrt{(\delta x)^2+(\delta y)^2}``

        **Multiplication and division**

        ``\bar{z}=(xy\ \text{or}\ x/y)\pm\delta z``

        ``\frac{\delta z}{|z|}=\sqrt{\left(\frac{\delta x}{x}\right)^2+\left(\frac{\delta y}{y}\right)^2}``

        **An arbitrary function**

        ``\delta f=\sqrt{\sum_i\left(\frac{\partial f}{\partial x_i}\delta x_i\right)^2}``

        > **Important**
        > Subtracting nominal values does not subtract their uncertainties: the independent > contributions still combine in quadrature.
        """);
        kind = :panel,
        title = "Linear error propagation"
    )
    canvas = webgrid(
        reshape([:cable_model, :propagation], 1, 2);
        columns = "minmax(0, 13fr) minmax(21rem, 7fr)",
        compact_columns = "minmax(0, 1fr)",
        rows = "minmax(0, 1fr)",
        stack_rows = "minmax(18rem, 1fr) auto",
        height = "100%",
        class = "lc-master-layout lc-layout-uncertainty-quantification",
        cable_model,
        propagation
    )
    return (body = slide(
        "Uncertainty quantification",
        canvas;
        lede = md"""
        Physical quantities are represented by nominal values with associated uncertainties, namely: ``\bar{x}=x\pm\delta x`` and ``\bar{y}=y\pm\delta y``, and propagated with linear error theory using `Measurements.jl`.
        """
    ),)
end

function warm_numerical_defaults!()
    set_preflight!(
        PREFLIGHT_STATE,
        :skin_effect,
        0.12,
        "Compiling the internal-impedance solution…"
    )
    yield()
    skin_model = validate(SkinEffectModel())
    skin_cache = build_skin_cache(skin_model)
    skin_state = solve_skin_state(1.0e3, 5.8e7, skin_model, skin_cache)

    set_preflight!(
        PREFLIGHT_STATE,
        :earth_return,
        0.42,
        "Compiling the earth-return field solution…"
    )
    yield()
    earth_model = validate(EarthReturnModel())
    earth_cache = build_earth_cache(earth_model)
    earth_state = solve_earth_state(
        50.0,
        500.0,
        100.0,
        10.0,
        0.0,
        earth_model,
        earth_cache
    )
    WARMUP_RESULT[] = (;
        skin_model,
        skin_cache,
        skin_state,
        earth_model,
        earth_cache,
        earth_state
    )
    return WARMUP_RESULT[]
end

function warm_dashboard_render!()
    set_preflight!(
        PREFLIGHT_STATE,
        :dashboard_render,
        0.72,
        "Compiling the interactive dashboard render…"
    )
    yield()
    WGLMakie.activate!()
    warm_app = App() do session
        state = setup(session)
        sources = sources_of_uncertainties_page(session, state)
        internal = internal_impedances_page(session, state)
        earth = earth_return_page(session, state)
        uncertainty = uncertainty_quantification_page(session, state)
        return DOM.div(
            sources.body...,
            internal.body...,
            earth.body...,
            uncertainty.body...
        )
    end
    parent = Session(NoConnection(); asset_server = NoServer())
    rendered = nothing
    try
        rendered = Bonito.show_html(IOBuffer(), warm_app; parent)
    finally
        isnothing(rendered) || close(rendered)
        close(parent)
    end
    return nothing
end

function warm_fundamentals!()
    try
        warm_numerical_defaults!()
        warm_dashboard_render!()
        set_preflight!(
            PREFLIGHT_STATE,
            :hot,
            1.0,
            "Fundamentals dashboards hot"
        )
    catch error
        set_preflight!(
            PREFLIGHT_STATE,
            :error,
            PREFLIGHT_STATE[].progress,
            "Fundamentals warmup failed: $(sprint(showerror, error))"
        )
        @error "ICHQP2026 Fundamentals warmup failed" exception = (
            error,
            catch_backtrace()
        )
    end
    return nothing
end

function start_warmup!()
    return lock(WARMUP_LOCK) do
        if preflight_ready(PREFLIGHT_STATE)
            return WARMUP_TASK[]
        end
        task = WARMUP_TASK[]
        !isnothing(task) && !istaskdone(task) && return task
        set_preflight!(
            PREFLIGHT_STATE,
            :queued,
            0.02,
            "Starting Fundamentals warmup…"
        )
        WARMUP_TASK[] = @async warm_fundamentals!()
        errormonitor(WARMUP_TASK[])
        return WARMUP_TASK[]
    end
end

const PREFLIGHT_RESOURCE = preflight_resource(
    id = "fundamentals-dashboards",
    title = "Fundamentals interactive dashboards",
    activate = start_warmup!,
    state = PREFLIGHT_STATE
)

const DECK = deck_descriptor(
    id = "fundamentals",
    group = "ICHQP2026",
    title = "Fundamentals",
    order = 20,
    render = true,
    setup = setup,
    resources = (PREFLIGHT_RESOURCE,),
    pages = (
        deck_page(
            "Sources of uncertainties";
            id = "sources-of-uncertainties",
            class = "lc-fill-page lc-fluid-type",
            build = sources_of_uncertainties_page
        ),
        deck_page(
            "Internal impedances";
            id = "internal-impedances",
            class = "lc-application-slide lc-fill-page",
            build = internal_impedances_page
        ),
        deck_page(
            "Earth return";
            id = "earth-return",
            class = "lc-application-slide lc-fill-page",
            build = earth_return_page
        ),
        deck_page(
            "Uncertainty quantification";
            id = "uncertainty-quantification",
            class = "lc-fill-page lc-fluid-type",
            build = uncertainty_quantification_page
        )
    )
)

end

ICHQP2026FundamentalsDeck.DECK
