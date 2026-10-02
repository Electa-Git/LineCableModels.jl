@testitem "Quality / Fatou / advisory source diagnostics" tags=[:quality] default_imports=false begin
    using Test, TOML, JSON3

    root = normpath(joinpath(@__DIR__, "..", ".."))
    config = joinpath(root, "fatou.toml")
    lint = TOML.parsefile(config)["lint"]
    selected = lint["select"]
    @test !isempty(selected) && allunique(selected)
    @test Set(keys(lint["severity"])) == Set(selected)
    @test all(==("warning"), values(lint["severity"]))
    target = TOML.parsefile(joinpath(root, "Project.toml"))["compat"]["julia"]

    function run_fatou(command)
        output, errors = IOBuffer(), IOBuffer()
        process = run(pipeline(ignorestatus(command); stdout=output, stderr=errors))
        return (; stdout=String(take!(output)), stderr=String(take!(errors)),
            exitcode=process.exitcode, termsignal=process.termsignal)
    end

    function check_fatou_version(executable)
        executable === nothing && error("Fatou 0.22.0 is required on PATH; install it before testing")
        result = run_fatou(`$executable --version`)
        result.exitcode == 0 && result.termsignal == 0 &&
            isempty(result.stderr) && strip(result.stdout) == "fatou 0.22.0" ||
            error("Expected Fatou 0.22.0: $(repr(result))")
        return nothing
    end

    function fatou_diagnostics(result, selected; io=devnull, root=pwd())
        result.termsignal == 0 || error("Fatou terminated by signal $(result.termsignal)")
        result.exitcode in (0, 1) || error("Unexpected Fatou exit $(result.exitcode)")
        # JSON mode puts findings on stdout. Configuration warnings and errors
        # remain on stderr, including unknown rule IDs and invalid Julia targets.
        any(line -> startswith(strip(line), "warning:") || startswith(strip(line), "error:"),
            split(result.stderr, '\n')) && error("Fatou setup failed: $(result.stderr)")
        diagnostics = JSON3.read(result.stdout)
        diagnostics isa JSON3.Array || error("Expected a Fatou diagnostic array")
        for diagnostic in diagnostics
            diagnostic isa JSON3.Object || error("Invalid Fatou diagnostic")
            all(key -> haskey(diagnostic, key), (:rule, :severity, :path, :range, :message)) ||
                error("Incomplete Fatou diagnostic")
            diagnostic.rule in selected || diagnostic.rule == "parse-error" ||
                error("Unexpected Fatou rule $(diagnostic.rule)")
            diagnostic.severity == (diagnostic.rule == "parse-error" ? "error" : "warning") ||
                error("Unexpected Fatou diagnostic severity")
            diagnostic.path isa String && isfile(diagnostic.path) ||
                error("Fatou diagnostic has no readable source path")
            diagnostic.range isa JSON3.Object &&
                all(key -> haskey(diagnostic.range, key), (:start, :end)) &&
                diagnostic.range.start isa Integer && diagnostic.range.end isa Integer &&
                0 <= diagnostic.range.start <= diagnostic.range.end <= filesize(diagnostic.path) ||
                error("Invalid Fatou byte range")
            diagnostic.message isa JSON3.Object && haskey(diagnostic.message, :body) &&
                diagnostic.message.body isa String || error("Invalid Fatou message")
        end
        result.exitcode == (isempty(diagnostics) ? 0 : 1) ||
            error("Fatou exit status does not agree with its diagnostic output")
        counts = Dict(rule => 0 for rule in selected)
        parse_failures = Set{String}()
        for diagnostic in diagnostics
            line, column = fatou_location(read(diagnostic.path, String), diagnostic.range.start)
            status = diagnostic.rule == "parse-error" ? "scan failure" : "advisory"
            println(io, "Fatou ", status, " [", diagnostic.rule, "] ", relpath(diagnostic.path, root),
                ":", line, ":", column, ": ", diagnostic.message.body)
            if diagnostic.rule == "parse-error"
                push!(parse_failures, diagnostic.path)
            else
                counts[diagnostic.rule] += 1
            end
        end
        for rule in sort(selected)
            println(io, "Fatou advisory count [", rule, "]: ", counts[rule])
        end
        println(io, "Fatou files with parse failures: ", length(parse_failures))
        isempty(parse_failures) || error("Fatou could not parse: $(join(sort!(collect(parse_failures)), ", "))")
        return diagnostics
    end

    # Fatou uses zero-based UTF-8 byte offsets. Report one-based character columns.
    function fatou_location(source, offset)
        0 <= offset <= ncodeunits(source) || error("Fatou offset is outside the source")
        index = offset + 1
        index == ncodeunits(source) + 1 || isvalid(source, index) ||
            error("Fatou offset splits a UTF-8 character")
        prefix = SubString(source, 1, prevind(source, index))
        line = count(==('\n'), prefix) + 1
        newline = findlast(==('\n'), prefix)
        column = length(newline === nothing ? prefix : SubString(prefix, nextind(prefix, newline))) + 1
        return line, column
    end

    executable = Sys.which("fatou")
    check_fatou_version(executable)
    println("Fatou 0.22.0; Julia compatibility target: ", target)

    @testset "Fatou integration controls" begin
        @test_throws ErrorException check_fatou_version(nothing)
        @test_throws ErrorException check_fatou_version(Base.julia_cmd())
        @test_throws Base.IOError run_fatou(`$(joinpath(tempdir(), "missing-fatou-$(getpid())")) --version`)
        clean = (; stdout="[]", stderr="", exitcode=0, termsignal=0)
        @test isempty(fatou_diagnostics(clean, selected))
        for bad in (merge(clean, (; stdout="")), merge(clean, (; stdout="{")),
                merge(clean, (; stdout="{}")), merge(clean, (; stdout="[{}]")),
                merge(clean, (; exitcode=1)), merge(clean, (; exitcode=7)),
                merge(clean, (; termsignal=15)),
                merge(clean, (; stderr="warning: unknown rule `typo`")),
                merge(clean, (; stderr="warning: invalid Julia version")))
            @test_throws Exception fatou_diagnostics(bad, selected)
        end
        terminated = run_fatou(`$(Base.julia_cmd()) --startup-file=no -e 'print("[]"); exit(7)'`)
        @test_throws ErrorException fatou_diagnostics(terminated, selected)
        @test fatou_location("αβ\nγδ", 0) == (1, 1)
        @test fatou_location("αβ\nγδ", 2) == (1, 2)
        @test fatou_location("αβ\nγδ", 5) == (2, 1)
        @test_throws ErrorException fatou_location("αβ", 1)

        mktempdir() do directory
            valid = joinpath(directory, "valid.jl")
            write(valid, "\"\"\"Value \$(1).\"\"\"\nidentity_α(α) = α\n")
            result = run_fatou(`$executable --config $config lint --julia-version $target --output json $valid`)
            @test isempty(fatou_diagnostics(result, selected))

            nested = mkpath(joinpath(directory, "nested"))
            mkdir(joinpath(directory, ".git"))
            write(joinpath(directory, ".gitignore"), "nested/\n")
            duplicate = joinpath(nested, "duplicate.jl")
            unused = joinpath(directory, "unused.jl")
            write(duplicate, "duplicate(x) = x\nduplicate(x) = x + 1\n")
            write(unused, "function unused(x)\n    α = x + 1\n    return x\nend\n")
            result = run_fatou(`$executable --config $config lint --julia-version $target --output json $valid $duplicate $unused`)
            diagnostics = fatou_diagnostics(result, selected)
            @test result.exitcode == 1 # Findings are reported, but this control passes.
            @test any(d -> d.rule == "duplicate-method" && d.path == duplicate, diagnostics)
            @test any(d -> d.rule == "unused-binding" && d.path == unused, diagnostics)
            walked = run_fatou(`$executable --config $config lint --julia-version $target --output json $directory`)
            @test all(d -> d.path != duplicate, fatou_diagnostics(walked, selected))

            versioned = joinpath(directory, "versioned.jl")
            write(versioned, "public identity_α\nidentity_α(α) = α\n")
            result = run_fatou(`$executable --config $config lint --julia-version 1.10 --output json $versioned`)
            @test any(d -> d.rule == "julia-version-compat", fatou_diagnostics(result, selected))

            invalid = joinpath(directory, "invalid.jl")
            write(invalid, "function broken(\n")
            result = run_fatou(`$executable --config $config lint --output json $invalid`)
            @test any(d -> d.rule == "parse-error", JSON3.read(result.stdout))
            @test_throws ErrorException fatou_diagnostics(result, selected)

            badconfig = joinpath(directory, "fatou.toml")
            write(badconfig, "[lint]\nselect = [\"misspelled-rule\"]\n")
            result = run_fatou(`$executable --config $badconfig lint --output json $valid`)
            @test_throws ErrorException fatou_diagnostics(result, selected)
            write(badconfig, "[lint\n")
            result = run_fatou(`$executable --config $badconfig lint --output json $valid`)
            @test_throws ErrorException fatou_diagnostics(result, selected)
            result = run_fatou(`$executable --config $config lint --output json $(joinpath(directory, "absent.jl"))`)
            @test_throws ErrorException fatou_diagnostics(result, selected)
        end
    end

    @testset "Fatou production scan" begin
        failures = 0
        try
            files = String[]
            for name in ("src", "ext")
                directory = joinpath(root, name)
                isdir(directory) || error("Missing production source root: $directory")
                for (path, _, names) in walkdir(directory), name in names
                    endswith(name, ".jl") && push!(files, joinpath(path, name))
                end
            end
            sort!(files)
            isempty(files) && error("Empty production source inventory")
            println("Fatou requested source files: ", length(files))
            result = run_fatou(`$executable --config $config lint --julia-version $target --output json $files`)
            # Retain native output too, including parse failures and suggested
            # fixes. Suggestions are data only. The scanner never invokes a fix option.
            println("Fatou native diagnostics: ", result.stdout)
            print(result.stderr)
            fatou_diagnostics(result, selected; io=stdout, root)
            @test result.exitcode in (0, 1)
        catch
            failures += 1
            rethrow()
        finally
            println("Fatou production scan failures: ", failures)
        end
    end
end
