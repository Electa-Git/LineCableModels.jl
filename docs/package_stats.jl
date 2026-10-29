using CairoMakie
using CSV
using Downloads
using GeoInterface
using GeoJSON

const PACKAGE_UUID = "dffb7669-8cd4-4628-99dd-7694e356ba37"
const LOG_URL = "https://julialang-logs.s3.amazonaws.com/public_outputs/current/package_requests_by_region.csv.gz"
const LAND_URL = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/v5.1.2/geojson/ne_110m_land.geojson"
const DEFAULT_OUTPUT = joinpath(@__DIR__, "src", "assets", "user-statistics.svg")

# Julia package-server regions are operational regions, not user countries.
const REGION_LOCATION = Dict(
    "us-east" => (label = "US East", longitude = -77.0, latitude = 38.0,
        offset = (12, 14), align = (:left, :bottom)),
    "us-west" => (label = "US West", longitude = -122.0, latitude = 38.0,
        offset = (-12, 14), align = (:right, :bottom)),
    "eu-central" => (label = "EU Central", longitude = 10.0, latitude = 50.0,
        offset = (0, 18), align = (:center, :bottom)),
    "eu-north" => (label = "EU North", longitude = 18.0, latitude = 60.0,
        offset = (0, 18), align = (:center, :bottom)),
    "au" => (label = "Australia", longitude = 151.0, latitude = -33.0,
        offset = (0, 18), align = (:center, :bottom)),
    "jp" => (label = "Japan", longitude = 139.7, latitude = 35.7,
        offset = (12, 14), align = (:left, :bottom)),
    "in" => (label = "India", longitude = 77.2, latitude = 28.6,
        offset = (0, 18), align = (:center, :bottom)),
    "kr" => (label = "Korea", longitude = 127.0, latitude = 37.5,
        offset = (-12, 14), align = (:right, :bottom)),
    "sa" => (label = "South America", longitude = -46.6, latitude = -23.5,
        offset = (0, 18), align = (:center, :bottom)),
    "sg" => (label = "Singapore", longitude = 103.8, latitude = 1.3,
        offset = (0, -18), align = (:center, :top)),
    "cn-east" => (label = "China East", longitude = 121.5, latitude = 31.2,
        offset = (-12, 14), align = (:right, :bottom)),
    "cn-northeast" => (label = "China Northeast", longitude = 123.4, latitude = 41.8,
        offset = (12, 14), align = (:left, :bottom)),
    "cn-southeast" => (label = "China Southeast", longitude = 113.3, latitude = 23.1,
        offset = (0, -18), align = (:center, :top))
)

function package_regions(path::AbstractString)
    counts = Dict{String, Float64}()
    minimum_date = nothing
    maximum_date = nothing

    for row in CSV.File(path)
        ismissing(row.package_uuid) && continue
        String(row.package_uuid) == PACKAGE_UUID || continue
        ismissing(row.client_type) && continue
        String(row.client_type) == "user" || continue
        ismissing(row.region) && continue
        region = String(row.region)
        haskey(REGION_LOCATION, region) || continue
        ismissing(row.request_addrs) && continue

        counts[region] = get(counts, region, 0.0) + Float64(row.request_addrs)
        if !ismissing(row.date_min)
            minimum_date = isnothing(minimum_date) ? row.date_min :
                           min(minimum_date, row.date_min)
        end
        if !ismissing(row.date_max)
            maximum_date = isnothing(maximum_date) ? row.date_max :
                           max(maximum_date, row.date_max)
        end
    end

    rows = [(; region, request_addrs) for (region, request_addrs) in counts]
    sort!(rows; by = row -> row.request_addrs, rev = true)
    dates = (; minimum = minimum_date, maximum = maximum_date)
    return first(rows, min(10, length(rows))), dates
end

function land_geometries(path::AbstractString)
    collection = GeoJSON.read(read(path, String))
    return GeoInterface.geometry.(collect(collection))
end

function draw_statistics(rows, dates, land_path::AbstractString, output::AbstractString)
    figure = Figure(
        size = (1050, 480),
        figure_padding = 8,
        backgroundcolor = RGBf(0.95, 0.97, 0.99)
    )
    axis = Axis(
        figure[2, 1];
        aspect = DataAspect(),
        limits = (-180, 180, -60, 88),
        width = 1000,
        height = 395,
        backgroundcolor = RGBf(0.68, 0.86, 0.97)
    )
    hidedecorations!(axis)
    hidespines!(axis)
    poly!(
        axis,
        land_geometries(land_path);
        color = RGBf(0.96, 0.96, 0.93),
        strokecolor = RGBf(0.68, 0.72, 0.73),
        strokewidth = 0.5
    )

    if isempty(rows)
        text!(axis, 0, 5;
            text = "No public package-server usage recorded yet",
            align = (:center, :center),
            fontsize = 24,
            color = RGBf(0.12, 0.18, 0.24))
        period = "awaiting the first registered-package requests"
    else
        maximum_count = maximum(row.request_addrs for row in rows)
        for row in rows
            location = REGION_LOCATION[String(row.region)]
            size = 14 + 32sqrt(row.request_addrs / maximum_count)
            scatter!(axis, [location.longitude], [location.latitude];
                markersize = size,
                color = RGBAf(0.05, 0.29, 0.53, 0.78),
                strokecolor = :white,
                strokewidth = 1.5)
            text!(axis, location.longitude, location.latitude;
                text = "$(location.label) · $(round(Int, row.request_addrs))",
                offset = location.offset,
                align = location.align,
                fontsize = 13,
                color = RGBf(0.05, 0.10, 0.15))
        end
        period = "$(dates.minimum) to $(dates.maximum)"
    end

    Label(figure[1, 1], "LineCableModels.jl package users";
        fontsize = 25,
        font = :bold,
        color = RGBf(0.05, 0.10, 0.15))
    Label(figure[3, 1],
        "Top Julia package-server regions · summed request addresses · $period";
        fontsize = 13,
        color = RGBf(0.22, 0.28, 0.34))
    rowgap!(figure.layout, 3)
    mkpath(dirname(output))
    save(output, figure)
    return output
end

function replace_statistics(rows, dates, land_path::AbstractString, output::AbstractString)
    mkpath(dirname(output))
    mktempdir(dirname(output)) do directory
        staged = joinpath(directory, basename(output))
        draw_statistics(rows, dates, land_path, staged)
        mv(staged, output; force = true)
    end
    return output
end

function main(args = ARGS)
    source = isempty(args) ? Downloads.download(LOG_URL) : abspath(args[1])
    output = length(args) < 2 ? DEFAULT_OUTPUT : abspath(args[2])
    land = length(args) < 3 ? Downloads.download(LAND_URL) : abspath(args[3])
    rows, dates = package_regions(source)
    result = replace_statistics(rows, dates, land, output)
    @info "Generated package-server usage map" result regions = length(rows)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
