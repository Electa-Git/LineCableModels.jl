# Label check. It reads the comments, docstrings and test item names under
# `test/quality`, `test/tools` and `test/support`, and the prose of `test/README.md` and
# `docs/src/developers.md` outside code. A label made of a capital letter and a number, a
# pull request number or a reference to a planning document fails the check. Code, the
# architecture baseline, `src/`, `ext/` and CI are out of scope, because their letters and
# numbers are domain notation.
@testmodule LabelCheck begin
    using JuliaSyntax: JuliaSyntax, @K_str

    const REPOSITORY = dirname(dirname(@__DIR__))
    const DIRECTORIES = ("test/quality", "test/tools", "test/support")
    const DOCUMENTS = ("test/README.md", "docs/src/developers.md")
    const EXCLUDED = ("test/quality/architecture_baseline.toml",)
    const PATTERNS = (r"\b[A-Z][0-9]{1,2}[a-z]?\b", r"\bPR\s*#?[0-9]+",
        r"(?i)\b(?:amendment|roadmap)s?\b")

    line_at(source, index) = count(==('\n'), SubString(source, 1, prevind(source, index))) + 1

    # The prose of a Julia source: comments, docstrings and test item names, each with
    # its line and kind. Code and other strings are left out.
    function julia_prose(source)
        found = Tuple{Int, String, String}[]
        for token in JuliaSyntax.tokenize(source)
            JuliaSyntax.kind(token) == K"Comment" || continue
            push!(found, (line_at(source, first(token.range)), "comment", source[token.range]))
        end
        line = Ref(0)
        visit(x) = if x isa LineNumberNode
            line[] = x.line
        elseif x isa Expr
            if Meta.isexpr(x, :macrocall) && length(x.args) >= 3 && x.args[2] isa LineNumberNode
                line[] = x.args[2].line
                name = x.args[1]
                if name == GlobalRef(Core, Symbol("@doc")) || name === Symbol("@doc")
                    push!(found, (line[], "docstring", text(x.args[3])))
                elseif name === Symbol("@testitem") && x.args[3] isa String
                    push!(found, (line[], "test item name", x.args[3]))
                end
            end
            foreach(visit, x.args)
        end
        visit(Meta.parseall(source))
        return found
    end

    # The literal parts of a string or an interpolated string.
    text(x::AbstractString) = String(x)
    text(x::Expr) = join(text(part) for part in x.args if part isa Union{AbstractString, Expr})
    text(_) = ""

    # The comments of a TOML source.
    function toml_prose(source)
        found = Tuple{Int, String, String}[]
        for (line, raw) in enumerate(split(source, '\n'))
            unquoted = replace(raw, r"\"(?:[^\"\\]|\\.)*\"" => "\"\"")
            index = findfirst('#', unquoted)
            index === nothing || push!(found, (line, "comment", unquoted[index:end]))
        end
        return found
    end

    # The prose of a Markdown source, outside fenced code blocks, code spans and link
    # targets.
    function markdown_prose(source)
        found = Tuple{Int, String, String}[]
        fenced = false
        for (line, raw) in enumerate(split(source, '\n'))
            startswith(lstrip(raw), "```") && (fenced = !fenced; continue)
            fenced && continue
            prose = replace(raw, r"(`+)(?:(?!\1).)*\1" => " ", r"\]\([^)]*\)" => "]")
            isempty(strip(prose)) || push!(found, (line, "prose", prose))
        end
        return found
    end

    # Each label, pull request number or planning word in `entries`, as a report line.
    function findings(path, entries)
        found = String[]
        for (line, kind, prose) in entries, pattern in PATTERNS, m in eachmatch(pattern, prose)
            push!(found, "$path:$line: `$(m.match)` in a $kind")
        end
        return found
    end

    function scope(root = REPOSITORY)
        files = String[]
        for directory in DIRECTORIES, (path, _, names) in walkdir(joinpath(root, directory))
            for name in names
                file = relpath(joinpath(path, name), root)
                endswith(name, ".jl") || endswith(name, ".toml") || continue
                endswith(name, "Manifest.toml") || file in EXCLUDED || push!(files, file)
            end
        end
        return sort!(append!(files, DOCUMENTS))
    end

    function live(root = REPOSITORY)
        found = String[]
        for file in scope(root)
            source = read(joinpath(root, file), String)
            entries = endswith(file, ".jl") ? julia_prose(source) :
                endswith(file, ".toml") ? toml_prose(source) : markdown_prose(source)
            append!(found, findings(file, entries))
        end
        return found
    end
end

@testitem "Quality / labels / comments, docstrings and documentation name things" tags=[:quality] setup=[LabelCheck] begin
    found = LabelCheck.live()
    foreach(println, found)
    @test found == String[]
end

@testitem "Quality / labels / negative controls" tags=[:quality] setup=[LabelCheck] begin
    L = LabelCheck
    # Each probe source has labels in its prose. The check reads only comments,
    # docstrings, test item names and documentation prose, never the code of this item.
    label, other = string('A', 1), string('C', 4, 'b')
    pull, planning = "PR 54", "roadmap"
    julia = """
        # Planted comment $(label).
        "Planted docstring $(pull)."
        f(x) = x
        $(label) = 1
        text = "a string $(other) in code"
        @testitem "planted $(planning)" begin
            #= block $(other) =#
        end
        """
    @test L.findings("probe.jl", L.julia_prose(julia)) == [
        "probe.jl:1: `$(label)` in a comment",
        "probe.jl:7: `$(other)` in a comment",
        "probe.jl:2: `$(pull)` in a docstring",
        "probe.jl:6: `$(planning)` in a test item name"]
    markdown = """
        Prose names $(label).

        Code `$(other)` stays code, and so does a [link](https://example.org/$(label)).

        ```julia
        $(label) = 1
        ```
        """
    @test L.findings("probe.md", L.markdown_prose(markdown)) == ["probe.md:1: `$(label)` in a prose"]
    toml = "key = \"$(label) # inside a string\" # comment $(pull)\n"
    @test L.findings("probe.toml", L.toml_prose(toml)) == ["probe.toml:1: `$(pull)` in a comment"]
    # Ordinary words, units and domain names without a capital-number label pass.
    @test isempty(L.findings("probe.md", [(1, "prose",
        "Float32 inputs, CO2, 5 %, IEC 60287, LineCableModels and the Commons guards")]))
    @test L.scope() ⊇ ["docs/src/developers.md", "test/README.md", "test/support/runner.jl"]
    @test "test/quality/architecture_baseline.toml" ∉ L.scope()
end
