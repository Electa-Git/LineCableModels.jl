const THEORY_SOURCE = joinpath(@__DIR__, "theory")
const THEORY_OUTPUT = joinpath(@__DIR__, "src", "theory")
const THEORY_PRIVATE_FILES = Set(("verification.md",))

function theory_files(source = THEORY_SOURCE)
    files = String[]
    for (directory, _, names) in walkdir(source), name in names
        relative = relpath(joinpath(directory, name), source)
        relative in THEORY_PRIVATE_FILES && continue
        push!(files, relative)
    end
    return sort!(files)
end

function documenter_matrix_math(content)
    lines = String[]
    display_math = false
    for line in split(content, '\n')
        if strip(line) == "```math"
            push!(lines, line)
            display_math = true
        elseif display_math
            push!(lines, line)
            strip(line) == "```" && (display_math = false)
        else
            push!(lines, replace(line, r"(?<!\\)\$([^$\n]+)(?<!\\)\$" =>
                matched -> "``" * matched[2:prevind(matched,lastindex(matched))] * "``"))
        end
    end
    display_math && error("unclosed display-math block in the matrix formulation")
    return join(lines, '\n')
end

function build_theory_pages!(source = THEORY_SOURCE, output = THEORY_OUTPUT)
    # This directory contains only generated copies; authored pages stay in docs/theory.
    rm(output; recursive = true, force = true)
    mkpath(output)
    pages = String[]
    for relative in theory_files(source)
        origin = joinpath(source, relative)
        destination = joinpath(output, relative)
        mkpath(dirname(destination))
        if endswith(relative, ".md")
            # Preserve GitHub fragments in the authored corpus. Documenter retains
            # the capitalization of its automatic heading identifiers.
            content = replace(
                read(origin, String),
                "#identification-and-source" => "#Identification-and-source"
            )
            relative == "matrix_formulation.md" &&
                (content = documenter_matrix_math(content))
            edit_path = relpath(origin, dirname(destination))
            metadata = "```@meta\nEditURL = $(repr(edit_path))\n```\n\n"
            write(destination, metadata, content)
            relative == "contents.md" || push!(pages, joinpath("theory", relative))
        else
            cp(origin, destination)
        end
    end
    return pages
end
