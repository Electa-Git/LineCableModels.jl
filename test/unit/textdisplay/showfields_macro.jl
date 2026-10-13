# The expansion recorded before the line-number walker of `@showfields` became methods. The
# input is parsed under the file name "golden".
@testitem "TextDisplay / @showfields / expansion is unchanged" tags=[:unit, :commons] begin
    code = """
        # Generated line numbers name the call on line 2. The body keeps its line 3.
        LineCableModels.TextDisplay.@showfields Point "Point" point -> begin
            (x = point.x, y = point.y)
        end"""
    expansion = macroexpand(LineCableModels,
        last(Meta.parseall(code; filename = "golden").args))
    lines = LineNumberNode[]
    collect!(x) = nothing
    collect!(x::LineNumberNode) = push!(lines, x)
    collect!(x::Expr) = foreach(collect!, x.args)
    collect!(expansion)
    # The body appears once in each show method, with the caller's line.
    @test count(==(LineNumberNode(3, :golden)), lines) == 2
    @test all(in((LineNumberNode(2, :golden), LineNumberNode(3, :golden))), lines)
    # The recorded text is the printing of Julia 1.12. Hygiene counters are removed.
    if VERSION.major == 1 && VERSION.minor == 12
        @test replace(string(expansion), r"#\d+#" => "#") == """
        begin
            #= golden:2 =#
            LineCableModels.TextDisplay.name(::LineCableModels.TextDisplay.Type{<:Point}) = begin
                    #= golden:2 =#
                    "Point"
                end
            #= golden:2 =#
            function (LineCableModels.TextDisplay.Base).summary(var"#io"::LineCableModels.TextDisplay.IO, var"#value"::Point)
                #= golden:2 =#
                #= golden:2 =#
                LineCableModels.TextDisplay.print(var"#io", LineCableModels.TextDisplay.name(LineCableModels.TextDisplay.typeof(var"#value")))
            end
            #= golden:2 =#
            function (LineCableModels.TextDisplay.Base).show(var"#io"::LineCableModels.TextDisplay.IO, var"#value"::Point)
                #= golden:2 =#
                #= golden:2 =#
                var"#displayed" = let point = var"#value"
                        #= golden:2 =#
                        begin
                            #= golden:2 =#
                            #= golden:3 =#
                            (x = point.x, y = point.y)
                        end
                    end
                #= golden:2 =#
                LineCableModels.TextDisplay.fields(var"#io", LineCableModels.TextDisplay.name(LineCableModels.TextDisplay.typeof(var"#value")), var"#displayed"; multiline = false)
            end
            #= golden:2 =#
            function (LineCableModels.TextDisplay.Base).show(var"#io"::LineCableModels.TextDisplay.IO, ::Base.Multimedia.MIME{Symbol("text/plain")}, var"#value"::Point)
                #= golden:2 =#
                #= golden:2 =#
                var"#displayed" = let point = var"#value"
                        #= golden:2 =#
                        begin
                            #= golden:2 =#
                            #= golden:3 =#
                            (x = point.x, y = point.y)
                        end
                    end
                #= golden:2 =#
                LineCableModels.TextDisplay.fields(var"#io", LineCableModels.TextDisplay.name(LineCableModels.TextDisplay.typeof(var"#value")), var"#displayed"; multiline = true)
            end
        end"""
    end
end
