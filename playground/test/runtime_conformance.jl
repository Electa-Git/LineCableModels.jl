# Test-only fixtures. Production dispatch is the inventory, not another registry.
module RuntimeConformance
using Test, InteractiveUtils, Bonito, UUIDs, JSON3
using LineCableModelsPlayground
const LCM = LineCableModelsPlayground
const Science = LCM.ScientificViews
const XRay = LCM.ComponentXRay

# These pre-v1 broker widgets retain their existing legacy regression coverage.
const LEGACY_SOURCES = ("src/widgets/JobControls.jl", "src/widgets/WorkerStatus.jl")
function owned_method(method)
    file = replace(relpath(string(method.file), LCM.PLAYGROUND_ROOT), '\\'=>'/')
    return !(file in LEGACY_SOURCES) &&
        (startswith(file,"src/widgets/") || startswith(file,"src/scientific/"))
end
function concrete_families(type)
    isabstracttype(type) || return [type]
    return reduce(vcat, concrete_families.(subtypes(type)); init=Any[])
end
function registered_families(callable)
    types = Any[]
    for method in methods(callable)
        owned_method(method) || continue
        signature = Base.unwrap_unionall(method.sig)
        append!(types, concrete_families(last(signature.parameters)))
    end
    return Set(types)
end

fixtures(::Type{WorkerSelector}, session, client) =
    (WorkerSelector(client,:parameters;profiles=("line-parameters",)),)
fixtures(::Type{PreparationStatus}, session, client) =
    (PreparationStatus(client,:parameters;parameters=Dict("private"=>"conformance-secret")),)
fixtures(::Type{ScientificJob}, session, client) =
    (ScientificJob(client,:parameters,"line.frequency_scan";parameters=Dict("private"=>"conformance-secret")),)
fixtures(::Type{WorkerDiagnostics}, session, client) = (WorkerDiagnostics(client),)
fixtures(::Type{WorkerControlPanel}, session, client) =
    (WorkerControlPanel(client,WorkerSelector(client,:parameters;profiles=("line-parameters",))),)
fixtures(::Type{JuliaTerminal}, session, client) = (JuliaTerminal(client,:terminal),)
fixtures(::Type{Science.CableGeometry}, session, client) = (Science.CableGeometry(),)
fixtures(::Type{Science.StudyRuntime}, session, client) = (Science.StudyRuntime(client),)
fixtures(::Type{Science.ScientificView}, session, client) =
    Tuple(Science.ScientificView(session,case,client) for case in
        (Science.StudyCases.LineParameters(),Science.StudyCases.CorridorImpedance()))

@testset "registration-derived runtime component conformance" begin
    inventory = registered_families(Bonito.jsrender)
    @test !isempty(inventory)
    @test inventory == registered_families(XRay.inspection)
    for type in sort!(collect(inventory);by=string)
        @testset "$type" begin
            # Adding a renderer/inspection without its normal fixture is a failure.
            @test hasmethod(fixtures,Tuple{Type{type},Session,RuntimeClient})
            for permitted in (false,true), cycle in 1:2
                session = Session()
                client = RuntimeClient(uuid4())
                try
                    XRay.set_policy!(session,XRay.XRayPolicy(;permitted))
                    for component in fixtures(type,session,client)
                        descriptor = XRay.inspection(component)
                        @test descriptor !== nothing
                        @test descriptor.source.line > 0
                        @test isfile(joinpath(LCM.PLAYGROUND_ROOT,descriptor.source.file))
                        @test !isempty(descriptor.css_scopes)
                        metadata = JSON3.write(XRay.inspection_payload(descriptor))
                        @test !occursin(string(client.run_id),metadata)
                        @test !occursin("conformance-secret",metadata)
                        # Serialize the actual normal renderer, including its own
                        # authored styles. No gallery markup or CSS substitutes.
                        html = sprint(show,MIME"text/html"(),Bonito.jsrender(session,component))
                        @test occursin("data-lcm-css-source",html)
                        styles = join((m.captures[1] for m in eachmatch(r"<style\b[^>]*>(.*?)</style>"s,html)),"\n")
                        @test !isempty(styles)
                        @test all(scope->occursin(scope,styles),descriptor.css_scopes)
                        @test occursin("__lcmXRay",html) == permitted
                    end
                finally
                    close(session)
                end
            end
        end
    end
end
end
