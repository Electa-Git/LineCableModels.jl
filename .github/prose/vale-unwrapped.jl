# Documenter renders each paragraph of a page or a docstring on one line. Several Vale rules
# match a phrase only within one line, so a hard wrap in the source hides it. This script
# gives Vale a copy of each input with every prose paragraph joined onto its first line. It
# then reports each alert at its line and column in the input.
#
# Usage: `julia --startup-file=no vale-unwrapped.jl --config=<vale.ini> [flag...] file...`.
# An argument that starts with a hyphen goes to Vale. The other arguments are the inputs,
# under the current directory. The environment variable `VALE` gives the path of the binary,
# `vale` by default. The output is one line per alert,
# `path:line:column:severity:check:message`, and the exit status is the status of Vale.
# The last line on stderr is the summary sentence of Vale, which counts the inputs.
#
# The copy of a Markdown file joins the whole file. The copy of a Julia file joins its
# docstrings, its `#=` blocks and its runs of `#` comments. Vale reads an input of another
# kind as it is. A copy has the line count of its input, since each joined line leaves a blank
# line in its place. The code stays on its lines for the tree-sitter view of Vale.
#
# The copy leaves the body of an admonition as it is, since its indentation under `!!!`
# makes it code. It also leaves as they are the docstrings that have text on the line of
# their opening quotes. The tokens come from `Base.JuliaSyntax`, which Base vendors but does
# not declare public. The tokenizer is the same from Julia 1.10 to 1.12, the version that CI
# runs. Where Vale places an alert at column 1 of a joined line, the script reports column 1
# of the first line of its paragraph.

using Base.JuliaSyntax: tokenize, kind, @K_str

const FENCE = r"^\s*(```|~~~)"
# The first line of a block other than a paragraph.
const BLOCK = r"^\s*(#{1,6}(\s|$)|[-*+]\s|\d+[.)]\s|\||>|!!!|<|\[[^\]]+\]:|\$\$|[=-]+\s*$|\$\(\w+\)\s*$)"
# The first line of a block that takes no continuation line. List items and quotes take one.
const CLOSED = r"^\s*(#{1,6}(\s|$)|\||!!!|<|\[[^\]]+\]:|\$\$|[=-]+\s*$|\$\(\w+\)\s*$)"
# The text before the opening quotes of a docstring, on their line.
const DOCOPEN = r"^\s*(@doc\s+)?(raw)?$"

indentation(line) = length(line) - length(lstrip(line))
hardbreak(line) = endswith(line, "  ") || endswith(line, "\\")

"""
Index of the paragraph line that each Markdown line continues. A line that starts a block
gets its own index. So does a blank line or a line of code.
"""
function paragraph_heads(lines)
    heads = collect(eachindex(lines))
    fenced = false
    head = 0 # first line of the open paragraph or list item, 0 when none is open
    for (i, line) in enumerate(lines)
        if occursin(FENCE, line)
            fenced = !fenced
            head = 0
        elseif fenced || isempty(strip(line))
            head = 0
        elseif head > 0 && !occursin(BLOCK, line) && !hardbreak(lines[i - 1])
            heads[i] = head
        elseif head == 0 && indentation(line) >= 4 && !occursin(BLOCK, line)
            head = 0 # indented code
        else
            head = occursin(CLOSED, line) ? 0 : i
        end
    end
    return heads
end

"""
Lines of prose in a source. `markdown` is the Markdown text of each line in `rows`. A line
becomes `emptied` once its words move to the first line of their paragraph.
"""
struct ProseRegion
    rows::Vector{Int}
    markdown::Vector{String}
    emptied::Vector{String}
end

"""
The lines of a docstring or of a `#=` block, less the indentation they share.
"""
function body_region(lines, rows)
    shared = minimum((indentation(lines[k]) for k in rows if !isempty(strip(lines[k]))); init = 0x7fff)
    markdown = [String(chop(lines[k]; head = min(shared, indentation(lines[k])), tail = 0))
                for k in rows]
    return ProseRegion(collect(rows), markdown, fill("", length(rows)))
end

"""
The bodies of the docstrings and `#=` blocks of a Julia source, and its runs of whole-line
`#` comments.
"""
function prose_regions(src::String, lines::Vector{String})
    bytes = codeunits(src)
    starts = [1; findall(==(UInt8('\n')), bytes) .+ 1]
    row(byte) = searchsortedlast(starts, byte)
    before(byte) = String(bytes[starts[row(byte)]:(byte - 1)])
    regions = ProseRegion[]
    comments = Tuple{Int, String, String}[] # row, indentation, text after the hash
    opening = nothing
    for t in tokenize(src)
        b = first(t.range)
        if kind(t) == K"\"\"\""
            if opening === nothing
                opening = t
                continue
            end
            first_row, last_row = row(first(opening.range)) + 1, row(b) - 1
            if occursin(DOCOPEN, before(first(opening.range))) &&
               bytes[last(opening.range) + 1] == UInt8('\n') && last_row >= first_row
                push!(regions, body_region(lines, first_row:last_row))
            end
            opening = nothing
        elseif kind(t) == K"Comment" && opening === nothing && isempty(strip(before(b)))
            text = String(bytes[t.range])
            if startswith(text, "#=")
                first_row, last_row = row(b) + 1, row(last(t.range)) - 1
                if last_row >= first_row && strip(lines[first_row - 1]) == "#=" &&
                   strip(lines[last_row + 1]) == "=#"
                    push!(regions, body_region(lines, first_row:last_row))
                end
            elseif text == "#" || startswith(text, "# ")
                words = chop(text; head = min(2, length(text)), tail = 0)
                push!(comments, (row(b), before(b), words))
            end
        end
    end
    # A run of comments is a sequence of consecutive rows at one indentation.
    run = Tuple{Int, String, String}[]
    function close_run()
        isempty(run) && return
        markdown, emptied = [c[3] for c in run], [c[2] * "#" for c in run]
        push!(regions, ProseRegion(first.(run), markdown, emptied))
        empty!(run)
    end
    for c in comments
        !isempty(run) && (c[1] != run[end][1] + 1 || c[2] != run[end][2]) && close_run()
        push!(run, c)
    end
    close_run()
    return regions
end

"""
Join each paragraph of `region` onto its first line in `lines`. Each entry of `pieces[row]`
gives a column of the copied line and the source row and column of the text there.
"""
function unwrap!(lines, pieces, region::ProseRegion)
    heads = paragraph_heads(region.markdown)
    for (k, row) in enumerate(region.rows)
        heads[k] == k && continue
        head = region.rows[heads[k]]
        words = lstrip(region.markdown[k])
        lines[head] = rstrip(lines[head])
        column = length(lines[row]) - length(words) + 1
        push!(pieces[head], (length(lines[head]) + 2, row, column))
        lines[head] *= " " * words
        lines[row] = region.emptied[k]
    end
end

"""
Write the unwrapped copy of `path` to `copied`. Return the column map of its lines.
"""
function write_unwrapped(path, copied)
    src = read(path, String)
    lines = String.(split(src, '\n'))
    pieces = [[(1, row, 1)] for row in eachindex(lines)]
    regions = if endswith(path, ".md")
        [ProseRegion(collect(eachindex(lines)), copy(lines), fill("", length(lines)))]
    else
        prose_regions(src, copy(lines))
    end
    foreach(region -> unwrap!(lines, pieces, region), regions)
    mkpath(dirname(copied))
    write(copied, join(lines, '\n'))
    return pieces
end

"""
The row and column in the source of the text at `column` of the copied line `row`.
"""
function source_location(pieces, row, column)
    piece = findlast(p -> p[1] <= column, pieces[row])
    start, source_row, source_column = pieces[row][piece]
    return source_row, source_column + column - start
end

# Vale writes one record per alert, and the character 0x1E ends each one. A message quotes
# its match, and a match can span two lines.
const TEMPLATE = """
{{- range .Files}}{{- \$path := .Path -}}{{- range .Alerts -}}
{{ \$path }}\t{{ .Line }}\t{{ index .Span 0 }}\t{{ .Severity }}\t{{ .Check }}\t{{ .Message }}\x1e
{{- end -}}{{- end -}}
"""

"""
The count `n` and its noun, singular for one, as Vale writes them.
"""
counted(n, noun) = n == 1 ? "1 $noun" : "$n $(noun)s"

function main(args)
    flags = filter(startswith("-"), args)
    inputs = filter(!startswith("-"), args)
    any(startswith("--output"), flags) && error("vale-unwrapped.jl sets the output of Vale")
    return mktempdir() do temporary
        mirror = realpath(temporary)
        maps = Dict{String, Tuple{String, Vector{Vector{NTuple{3, Int}}}}}()
        targets = String[]
        for input in inputs
            relative = relpath(abspath(input), pwd())
            startswith(relative, "..") && error("$input is outside the current directory")
            if endswith(input, ".md") || endswith(input, ".jl")
                copied = joinpath(mirror, relative)
                maps[copied] = (relative, write_unwrapped(input, copied))
                push!(targets, copied)
            else
                push!(targets, input)
            end
        end
        template = joinpath(mirror, "alerts.tmpl")
        write(template, TEMPLATE)
        vale = get(ENV, "VALE", "vale")
        process = open(pipeline(ignorestatus(`$vale $flags --output=$template $targets`); stderr))
        records = split(read(process, String), '\x1e'; keepempty = false)
        wait(process)
        counts = Dict("error" => 0, "warning" => 0, "suggestion" => 0)
        for record in records
            path, row, column, severity, check, message = split(record, '\t'; limit = 6)
            row, column = parse(Int, row), parse(Int, column)
            if haskey(maps, path)
                path, pieces = maps[path]
                row, column = source_location(pieces, row, column)
            end
            counts[severity] = get(counts, severity, 0) + 1
            println("$path:$row:$column:$severity:$check:", replace(message, '\n' => ' '))
        end
        # Flush the alerts first, so that a stream that joins the two outputs ends with the
        # summary.
        flush(stdout)
        severities = ("error", "warning", "suggestion")
        alerts = join((counted(counts[k], k) for k in severities), ", ", " and ")
        println(stderr, alerts, " in ", counted(length(targets), "file"), ".")
        return process.exitcode
    end
end

exit(main(ARGS))
