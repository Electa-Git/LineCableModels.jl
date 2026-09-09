import Bonito

struct XRayPolicyProbe
    calls::Base.RefValue{Int}
end
function ComponentXRay.inspection(component::XRayPolicyProbe)
    component.calls[] += 1
    return ComponentXRay.ComponentInspection(:policy_probe)
end

@testset "X-ray policy gates inspection hooks" begin
    session = Bonito.Session(Bonito.NoConnection(); asset_server=Bonito.NoServer())
    calls = Ref(0)
    component = XRayPolicyProbe(calls)
    node = Bonito.DOM.div("Uninstrumented content")
    try
        @test ComponentXRay.instrument(session, node, component) === node
        @test calls[] == 0
        ComponentXRay.set_policy!(session, ComponentXRay.XRayPolicy(true))
        @test ComponentXRay.instrument(session, node, component) !== node
        @test calls[] == 1
        ComponentXRay.set_policy!(session, ComponentXRay.XRayPolicy())
        @test ComponentXRay.instrument(session, node, component) === node
        @test calls[] == 1
    finally
        close(session)
    end
end

@testset "X-ray CSS preview contract" begin
    X = LineCableModelsPlayground.ComponentXRay
    @test X.XRayPolicy(true).css_preview
    @test !X.XRayPolicy(; permitted=true, css_preview=false).css_preview
    @test_throws ArgumentError X.CssEditor(:script)
    @test_throws ArgumentError X.CssEditor(:number; minimum=2, maximum=1)
    @test_throws ArgumentError X.CssEditor(:number; step=0)
    @test_throws ArgumentError X.CssEditor(:number; maximum=Inf)
    descriptor = X.ComponentInspection(:example; css_scopes=[".example"],
        css_overrides=Dict("gap" => X.CssEditor(:length; units=["px"], maximum=30, step=2)))
    payload = X.inspection_payload(descriptor)
    @test payload["css_editors"]["gap"]["maximum"] == 30
    @test payload["css_editors"]["gap"]["units"] == ["px"]
    @test isempty(X.css_editors(:example))
    @test !haskey(X.CSS_EDITORS, "transform")
    @test !haskey(X.CSS_EDITORS, "--lc-focus")
    @test X.CSS_EDITORS["display"].kind == :choice
    @test X.CSS_EDITORS["color"].kind == :color
    script = X.XRAY_SCRIPT
    pointer = match(r"function pointerMove\(event\) \{(.*?)\n  \}"s, script).captures[1]
    @test !occursin("render(", pointer)
    @test !occursin("clearSelection", pointer)
    @test occursin("bindingCells.get(entry.name).textContent", script)
    @test occursin("preview.destroy()", script)
    @test occursin("CSS.supports", script)
    @test !occursin("eval(", script)
end
