using CairoMakie
using CSV
using DataFrames
using Downloads

const PACKAGE_UUID = "dffb7669-8cd4-4628-99dd-7694e356ba37"
const LOG_URL = "https://julialang-logs.s3.amazonaws.com/public_outputs/current/package_requests_by_region.csv.gz"
const DEFAULT_OUTPUT = joinpath(@__DIR__, "src", "assets", "user-statistics.svg")

# Julia package-server regions are operational regions, not user countries.
const REGION_LOCATION = Dict(
    "us-east" => (label = "US East", longitude = -77.0, latitude = 38.0),
    "us-west" => (label = "US West", longitude = -122.0, latitude = 38.0),
    "eu-central" => (label = "EU Central", longitude = 10.0, latitude = 50.0),
    "eu-north" => (label = "EU North", longitude = 18.0, latitude = 60.0),
    "au" => (label = "Australia", longitude = 151.0, latitude = -33.0),
    "jp" => (label = "Japan", longitude = 139.7, latitude = 35.7),
    "in" => (label = "India", longitude = 77.2, latitude = 28.6),
    "kr" => (label = "Korea", longitude = 127.0, latitude = 37.5),
    "sa" => (label = "South America", longitude = -46.6, latitude = -23.5),
    "sg" => (label = "Singapore", longitude = 103.8, latitude = 1.3),
    "cn-east" => (label = "China East", longitude = 121.5, latitude = 31.2),
    "cn-northeast" => (label = "China Northeast", longitude = 123.4, latitude = 41.8),
    "cn-southeast" => (label = "China Southeast", longitude = 113.3, latitude = 23.1)
)

# A deliberately low-detail silhouette keeps the generated asset small and
# avoids adding a geospatial plotting dependency to the documentation project.
const LAND = [
    [(-168, 15), (-165, 70), (-120, 74), (-82, 50), (-53, 52), (-82, 10),
        (-102, 8), (-116, 28)],
    [(-82, 12), (-52, 8), (-35, -10), (-54, -56), (-76, -50), (-81, -5)],
    [(-18, 35), (2, 72), (45, 70), (65, 55), (105, 78), (180, 65),
        (145, 40), (122, 22), (104, 5), (52, 12), (35, 30)],
    [(-18, 35), (15, 37), (51, 12), (40, -35), (18, -35), (-5, 5)],
    [(112, -10), (154, -10), (153, -42), (115, -36)],
    [(-52, 60), (-25, 82), (-12, 62)],
    [(45, -13), (51, -12), (50, -25), (43, -24)]
]

function package_regions(path::AbstractString)
    frame = CSV.read(path, DataFrame)
    selected = filter(frame) do row
        row.package_uuid == PACKAGE_UUID && !ismissing(row.client_type) &&
            row.client_type == "user" && !ismissing(row.region) &&
            haskey(REGION_LOCATION, String(row.region))
    end
    if isempty(selected)
        return (
            DataFrame(region = String[], request_addrs = Float64[]),
            (minimum = nothing, maximum = nothing)
        )
    end
    grouped = combine(
        groupby(selected, :region),
        :request_addrs => (values -> sum(skipmissing(values))) => :request_addrs
    )
    sort!(grouped, :request_addrs; rev = true)
    dates = (
        minimum = minimum(skipmissing(selected.date_min)),
        maximum = maximum(skipmissing(selected.date_max))
    )
    return first(grouped, min(10, nrow(grouped))), dates
end

function draw_statistics(rows, dates, output::AbstractString)
    figure = Figure(size = (980, 540), backgroundcolor = RGBf(0.95, 0.97, 0.99))
    axis = Axis(
        figure[1, 1];
        aspect = DataAspect(),
        backgroundcolor = RGBf(0.68, 0.86, 0.97),
        limits = (-180, 180, -60, 88)
    )
    hidedecorations!(axis)
    hidespines!(axis)
    for polygon in LAND
        poly!(
            axis,
            Point2f.(first.(polygon), last.(polygon));
            color = RGBf(0.96, 0.96, 0.93),
            strokecolor = RGBf(0.78, 0.80, 0.79),
            strokewidth = 0.8
        )
    end

    if isempty(rows)
        text!(
            axis,
            0,
            5;
            text = "No public package-server usage recorded yet",
            align = (:center, :center),
            fontsize = 24,
            color = RGBf(0.12, 0.18, 0.24)
        )
        period = "awaiting the first registered-package requests"
    else
        maximum_count = maximum(rows.request_addrs)
        for row in eachrow(rows)
            location = REGION_LOCATION[String(row.region)]
            size = 14 + 32sqrt(row.request_addrs / maximum_count)
            scatter!(
                axis,
                [location.longitude],
                [location.latitude];
                markersize = size,
                color = RGBAf(0.05, 0.29, 0.53, 0.78),
                strokecolor = :white,
                strokewidth = 1.5
            )
            text!(
                axis,
                location.longitude,
                location.latitude + 5;
                text = "$(location.label) · $(round(Int, row.request_addrs))",
                align = (:center, :bottom),
                fontsize = 13,
                color = RGBf(0.05, 0.10, 0.15)
            )
        end
        period = "$(dates.minimum) to $(dates.maximum)"
    end

    Label(
        figure[0, 1],
        "LineCableModels.jl package clients";
        fontsize = 25,
        font = :bold,
        color = RGBf(0.05, 0.10, 0.15)
    )
    Label(
        figure[2, 1],
        "Top Julia package-server regions · summed request addresses · $period";
        fontsize = 13,
        color = RGBf(0.22, 0.28, 0.34)
    )
    mkpath(dirname(output))
    save(output, figure)
    return output
end

function main(args = ARGS)
    source = isempty(args) ? Downloads.download(LOG_URL) : abspath(args[1])
    output = length(args) < 2 ? DEFAULT_OUTPUT : abspath(args[2])
    rows, dates = package_regions(source)
    result = draw_statistics(rows, dates, output)
    @info "Generated package-server usage map" result regions = nrow(rows)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
