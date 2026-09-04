using CairoMakie
using JLD2
using SHA

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const FEM_ROOT = joinpath(
    REPOSITORY_ROOT,
    ".linecablemodels",
    "fem",
    "gauntlet",
    "corrected_fullband"
)
const CORPUS_ROOT = joinpath(
    REPOSITORY_ROOT,
    ".linecablemodels",
    "fem",
    "gauntlet",
    "analytical_fullband_corpus"
)
const CASE_ROOT = joinpath(CORPUS_ROOT, "cases")
const DEFAULT_ADMITTANCE_OUTPUT_ROOT = joinpath(
    CORPUS_ROOT,
    "plots",
    "fem_papadopoulos2010_xue2018_infinite"
)
const DEFAULT_IMPEDANCE_OUTPUT_ROOT = joinpath(
    CORPUS_ROOT,
    "plots",
    "fem_papadopoulos2010_xue2018_impedance"
)
const DEFAULT_COMBINED_OUTPUT_ROOT = joinpath(
    CORPUS_ROOT,
    "plots",
    "fem_papadopoulos2010_xue2018"
)
const VALID_FORMATS = Set(("png", "pdf", "svg"))

function usage(io = stdout)
    println(io, """
    Plot persisted FEM, Papadopoulos-2010, and Xue-2018 matrix results.

    Usage:
      julia --project=docs test/gauntlet/plot_fem_matrix_comparison.jl [options] [CASE ...]

    Options:
      --output=DIR    Output directory (default: quantity-specific corpus directory)
      --format=FORMAT png, pdf, or svg (default: png)
      --quantity=NAME admittance, impedance, or both (default: admittance)
      --help          Show this message

    With no CASE arguments, every case having all three persisted artifacts is
    plotted. Each case produces one full matrix grid for the real part and one
    for the imaginary part. Every matrix-element axis overlays FEM,
    Papadopoulos 2010, and Xue 2018. Admittance uses the default Xue 2018
    infinite-depth self/mutual route.
    """)
end

function parse_arguments(arguments)
    output_root = nothing
    format = "png"
    quantity = :admittance
    cases = String[]
    for argument in arguments
        if argument == "--help" || argument == "-h"
            usage()
            exit()
        elseif startswith(argument, "--output=")
            output_root = abspath(split(argument, '='; limit = 2)[2])
        elseif startswith(argument, "--format=")
            format = lowercase(split(argument, '='; limit = 2)[2])
        elseif startswith(argument, "--quantity=")
            value = lowercase(split(argument, '='; limit = 2)[2])
            value in ("admittance", "impedance", "both") || throw(ArgumentError(
                "unsupported quantity '$value'; choose admittance, impedance, or both"
            ))
            quantity = Symbol(value)
        elseif startswith(argument, '-')
            throw(ArgumentError("unknown option: $argument"))
        else
            push!(cases, argument)
        end
    end
    format in VALID_FORMATS || throw(ArgumentError(
        "unsupported format '$format'; choose png, pdf, or svg"
    ))
    return (; output_root, format, quantity, cases)
end

function quantity_config(quantity)
    if quantity === :admittance
        return (
            matrix_key = "Y",
            matrix_label = "Y",
            unit = "S/m",
            formula_field = "earth_admittance",
            papadopoulos_file = "earth_admittance_papadopoulos2010.jld2",
            xue_file = "earth_admittance_xue2018.jld2",
            xue_legend = "Xue 2018 infinite depth"
        )
    elseif quantity === :impedance
        return (
            matrix_key = "Z",
            matrix_label = "Z",
            unit = "Ω/m",
            formula_field = "earth_impedance",
            papadopoulos_file = "earth_impedance_papadopoulos2010.jld2",
            xue_file = "earth_impedance_xue2018.jld2",
            xue_legend = "Xue 2018"
        )
    end
    throw(ArgumentError("unsupported plot quantity: $quantity"))
end

function default_output_root(quantity)
    quantity === :admittance && return DEFAULT_ADMITTANCE_OUTPUT_ROOT
    quantity === :impedance && return DEFAULT_IMPEDANCE_OUTPUT_ROOT
    quantity === :both && return DEFAULT_COMBINED_OUTPUT_ROOT
    throw(ArgumentError("unsupported plot quantity: $quantity"))
end

function artifact_paths(case_id, quantity)
    config = quantity_config(quantity)
    return (
        fem = joinpath(FEM_ROOT, case_id, "reference.jld2"),
        papadopoulos = joinpath(CASE_ROOT, case_id, config.papadopoulos_file),
        xue = joinpath(CASE_ROOT, case_id, config.xue_file)
    )
end

function available_cases(quantity)
    isdir(FEM_ROOT) || error("persisted FEM root is missing: $FEM_ROOT")
    return sort!(filter(readdir(FEM_ROOT)) do case_id
        paths = artifact_paths(case_id, quantity)
        isfile(paths.fem) && isfile(paths.papadopoulos) && isfile(paths.xue)
    end)
end

file_sha256(path) = bytes2hex(sha256(read(path)))

function require_key(document, key, path)
    haskey(document, key) || error("$path does not contain '$key'")
    return document[key]
end

function require_equal(label, actual, expected, path)
    actual == expected || error(
        "$path has incompatible $label: got $(repr(actual)), " *
        "expected $(repr(expected))"
    )
end

function load_case(case_id, quantity)
    config = quantity_config(quantity)
    paths = artifact_paths(case_id, quantity)
    for path in paths
        isfile(path) || error("persisted comparison input is missing: $path")
    end

    fem = JLD2.load(paths.fem)
    papadopoulos = JLD2.load(paths.papadopoulos)
    xue = JLD2.load(paths.xue)

    require_equal("FEM schema", require_key(fem, "schema_version", paths.fem), 3, paths.fem)
    require_equal(
        "Papadopoulos formula",
        require_key(papadopoulos, config.formula_field, paths.papadopoulos),
        :Papadopoulos2010,
        paths.papadopoulos
    )
    require_equal(
        "Xue formula",
        require_key(xue, config.formula_field, paths.xue),
        :Xue2018,
        paths.xue
    )
    require_equal(
        "Xue candidate",
        require_key(xue, "formula", paths.xue),
        :Xue2018,
        paths.xue
    )

    fem_digest = file_sha256(paths.fem)
    frequencies = require_key(fem, "frequencies", paths.fem)
    port_order = require_key(fem, "port_order", paths.fem)
    fem_values = require_key(fem, config.matrix_key, paths.fem)

    for (label, document, path) in (
        ("Papadopoulos", papadopoulos, paths.papadopoulos),
        ("Xue", xue, paths.xue)
    )
        status = require_key(document, "status", path)
        status in (:execution_failure, :not_applicable) && error(
            "$label artifact for $case_id cannot be plotted: status=$status"
        )
        require_equal(
            "$label FEM-reference digest",
            require_key(document, "fem_reference_sha256", path),
            fem_digest,
            path
        )
        require_equal(
            "$label frequencies",
            require_key(document, "frequencies", path),
            frequencies,
            path
        )
        require_equal(
            "$label port order",
            require_key(document, "port_order", path),
            port_order,
            path
        )
    end

    papadopoulos_values = require_key(
        papadopoulos, config.matrix_key, paths.papadopoulos
    )
    xue_values = require_key(xue, config.matrix_key, paths.xue)
    expected_size = (length(port_order), length(port_order), length(frequencies))
    for (label, values) in (
        ("FEM $(config.matrix_label)", fem_values),
        ("Papadopoulos $(config.matrix_label)", papadopoulos_values),
        ("Xue $(config.matrix_label)", xue_values)
    )
        size(values) == expected_size || error(
            "$case_id $label has size $(size(values)); expected $expected_size"
        )
        all(isfinite, values) || error("$case_id $label contains non-finite values")
    end

    return (;
        case_id,
        quantity,
        matrix_label = config.matrix_label,
        unit = config.unit,
        xue_legend = config.xue_legend,
        frequencies,
        port_order = string.(port_order),
        fem_values,
        papadopoulos_values,
        xue_values,
        paths,
        fem_digest
    )
end

function terminal_map(port_order)
    entries = ["$index=$(port_order[index])" for index in eachindex(port_order)]
    lines = String[]
    for first_index in 1:3:length(entries)
        last_index = min(first_index + 2, length(entries))
        push!(lines, join(entries[first_index:last_index], "    "))
    end
    return join(lines, '\n')
end

function matrix_plot(data, component::Symbol)
    transform, component_label = if component === :real
        real, "real($(data.matrix_label))"
    elseif component === :imaginary
        imag, "imaginary($(data.matrix_label))"
    else
        throw(ArgumentError("unsupported component: $component"))
    end

    terminal_count = length(data.port_order)
    figure = Figure(size = (
        max(1150, 290 * terminal_count),
        max(950, 245 * terminal_count + 180)
    ))
    Label(
        figure[0, 1:terminal_count],
        "$(data.case_id) — $component_label per matrix element";
        fontsize = 24,
        font = :bold,
        tellwidth = false
    )

    colours = (
        fem = :black,
        papadopoulos = :darkorange2,
        xue = :dodgerblue3
    )
    styles = (fem = :solid, papadopoulos = :dash, xue = :dot)

    for row in 1:terminal_count, column in 1:terminal_count

        axis = Axis(
            figure[row, column];
            xscale = log10,
            title = "$(data.matrix_label)[$row,$column]",
            titlesize = 13,
            xlabel = row == terminal_count ? "f (Hz)" : "",
            ylabel = column == 1 ? "$component_label ($(data.unit))" : "",
            xlabelsize = 11,
            ylabelsize = 11,
            xticklabelsize = 9,
            yticklabelsize = 9
        )
        hlines!(axis, [0.0]; color = (:grey45, 0.45), linewidth = 0.8)
        lines!(
            axis,
            data.frequencies,
            transform.(@view data.fem_values[row, column, :]);
            color = colours.fem,
            linestyle = styles.fem,
            linewidth = 2.5
        )
        lines!(
            axis,
            data.frequencies,
            transform.(@view data.papadopoulos_values[row, column, :]);
            color = colours.papadopoulos,
            linestyle = styles.papadopoulos,
            linewidth = 2.0
        )
        lines!(
            axis,
            data.frequencies,
            transform.(@view data.xue_values[row, column, :]);
            color = colours.xue,
            linestyle = styles.xue,
            linewidth = 2.0
        )
        row == terminal_count || hidexdecorations!(axis; grid = false)
    end

    legend_elements = [
        LineElement(color = colours.fem, linestyle = styles.fem, linewidth = 3),
        LineElement(
            color = colours.papadopoulos,
            linestyle = styles.papadopoulos,
            linewidth = 3
        ),
        LineElement(color = colours.xue, linestyle = styles.xue, linewidth = 3)
    ]
    Legend(
        figure[terminal_count + 1, 1:terminal_count],
        legend_elements,
        ["FEM", "Papadopoulos 2010", data.xue_legend];
        orientation = :horizontal,
        tellwidth = false,
        framevisible = false
    )
    Label(
        figure[terminal_count + 2, 1:terminal_count],
        terminal_map(data.port_order);
        fontsize = 11,
        halign = :left,
        justification = :left,
        tellwidth = false
    )
    rowgap!(figure.layout, 7)
    colgap!(figure.layout, 7)
    return figure
end

function plot_case(case_id, quantity, output_root, format)
    data = load_case(case_id, quantity)
    case_output = joinpath(output_root, case_id)
    mkpath(case_output)
    outputs = String[]
    for component in (:real, :imaginary)
        output = joinpath(
            case_output,
            "$(data.matrix_label)_$(component).$format"
        )
        save(output, matrix_plot(data, component))
        push!(outputs, output)
    end
    println(
        "PLOTTED\t", case_id,
        "\tquantity=", quantity,
        "\tterminals=", length(data.port_order),
        "\tfrequencies=", length(data.frequencies),
        "\tfem_sha256=", data.fem_digest,
        "\toutputs=", join(outputs, ',')
    )
    flush(stdout)
    return outputs
end

function main(arguments)
    options = parse_arguments(arguments)
    quantities = options.quantity === :both ?
                 (:admittance, :impedance) : (options.quantity,)
    discovered = intersect((available_cases(quantity) for quantity in quantities)...)
    cases = isempty(options.cases) ? discovered : unique(options.cases)
    isempty(cases) && error("no complete persisted FEM/Papadopoulos/Xue cases found")
    output_root = something(options.output_root, default_output_root(options.quantity))
    mkpath(output_root)
    println(
        "PLAN\tcases=", length(cases),
        "\tquantities=", join(quantities, ','),
        "\tformat=", options.format,
        "\toutput=", output_root
    )
    for case_id in cases
        for quantity in quantities
            plot_case(case_id, quantity, output_root, options.format)
        end
    end
    println("COMPLETE\tcases=", length(cases), "\toutput=", output_root)
end

main(ARGS)
