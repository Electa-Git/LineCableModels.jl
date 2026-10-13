# Each expected expression is the expansion recorded before the selector check of `@observe`
# became methods. Inputs are parsed under the file name "golden".
@testitem "Commons / @observe / expansions are unchanged" tags=[:unit, :commons] begin
    expand(code) = macroexpand(LineCableModels,
        last(Meta.parseall(code; filename = "golden").args))
    observe = GlobalRef(LineCableModels.Commons, :observe)
    @test expand("@observe Z[1, 2]") == :((Z, 1, 2))
    @test expand("@observe (Z, abs)[:, 1]") == :((Z, abs, :, 1))
    @test expand("@observe source Z[1, 2, 3]") == :($observe(source, Z, 1, 2, 3))
    @test expand("@observe source (statistics, Z, mean)[1, 2]") ==
          :($observe(source, statistics, Z, mean, 1, 2))
end
