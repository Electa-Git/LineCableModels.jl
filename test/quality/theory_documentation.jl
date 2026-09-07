@testitem "Theory / publication / contents and hidden formula pages" tags=[:quality] begin
    root = pkgdir(LineCableModels)
    source = joinpath(root, "docs", "theory")
    include(joinpath(root, "docs", "theory_pages.jl"))

    contents = read(joinpath(source, "contents.md"), String)
    expected = [
        "matrix_formulation.md",
        "internal_impedance.md",
        "insulation_parameters.md",
        "earth_return_impedance.md",
        "earth_return_admittance.md",
        "modal_decomposition.md",
        "earth_properties.md"
    ]
    links = [String(m[1]) for m in eachmatch(r"(?m)^- \[[^\]]+\]\(([^)]+)\)$", contents)]
    @test startswith(contents, "# Contents\n")
    @test links == expected
    @test all(page -> isfile(joinpath(source, page)), expected)

    files = theory_files()
    @test "contents.md" in files
    @test "verification.md" ∉ files
    @test "index-list-for-PR-report.md" ∉ files
    @test "gaps-append-to-pr-report.md" ∉ files
    @test !any(path -> occursin("PR-theory-tracker", path), files)
    before = Dict(path => read(joinpath(source, path), String) for path in files)
    mktempdir() do temporary
        output = joinpath(temporary, "theory")
        pages = build_theory_pages!(source, output)
        @test length(pages) == 165
        @test allunique(pages)
        @test joinpath("theory", "contents.md") ∉ pages
        @test all(page -> joinpath("theory", page) in pages, expected)
        @test sort(theory_files(output)) == files
        for path in files
            @test read(joinpath(source, path), String) == before[path]
            endswith(path, ".md") || continue
            @test !occursin(r"(?m)^\s*\$\$\s*$", before[path])
            published = read(joinpath(output, path), String)
            edit_path = relpath(joinpath(source, path), dirname(joinpath(output, path)))
            metadata = "```@meta\nEditURL = $(repr(edit_path))\n```\n\n"
            normalized = replace(before[path],
                "#identification-and-source" => "#Identification-and-source")
            path == "matrix_formulation.md" &&
                (normalized = documenter_matrix_math(normalized))
            @test published == metadata * normalized
            for link in eachmatch(r"\]\(([^)\n]+)\)", published)
                target = first(split(String(link[1]), '#'))
                occursin(r"^[A-Za-z][A-Za-z0-9+.-]*:", target) && continue
                endswith(target, ".md") || endswith(target, ".tsv") ||
                    endswith(target, ".bib") || continue
                @test isfile(normpath(joinpath(dirname(joinpath(output, path)), target)))
            end
        end
    end

    matrix = read(joinpath(source, "matrix_formulation.md"), String)
    rendered_matrix = documenter_matrix_math(matrix)
    original_equations = [m[1] for m in eachmatch(r"(?ms)^```math\n(.*?)^```$", matrix)]
    rendered_equations = [m[1] for m in eachmatch(r"(?ms)^```math\n(.*?)^```$", rendered_matrix)]
    @test length(original_equations) == 25
    @test rendered_equations == original_equations
    @test !occursin('$', rendered_matrix)
    prose = replace(matrix, r"(?ms)^```math\n.*?^```$" => "")
    for expression in eachmatch(r"(?<!\\)\$([^$\n]+)(?<!\\)\$", prose)
        @test occursin("``" * expression[1] * "``", rendered_matrix)
    end
    @test occursin(raw"\begin{bmatrix}1\\1\end{bmatrix}", rendered_matrix)
    @test occursin(raw"\begin{bmatrix}1\\1\\1\end{bmatrix}", rendered_matrix)
    @test occursin(raw"\begin{bmatrix}\mathbf I\\I_p\end{bmatrix}", rendered_matrix)
    @test occursin(raw"\begin{bmatrix}\mathbf V\\V_p\end{bmatrix}", rendered_matrix)

    navigation = read(joinpath(root, "docs", "make.jl"), String)
    @test occursin("\"Home\" => \"index.md\",\n        hide(\"Theory\" => \"theory/contents.md\", theory_pages)", navigation)
    @test occursin("theory_pages = build_theory_pages!()", navigation)
    @test !occursin("\"literature-survey\"", navigation)
    stylesheet = read(joinpath(root, "docs", "src", "assets", "custom.css"), String)
    @test occursin("li:has(> a.tocitem[href\$=\"contents.html\"]) > ul", stylesheet)
    @test occursin("li:has(> a.tocitem[href\$=\"contents/\"]) > ul", stylesheet)
end
