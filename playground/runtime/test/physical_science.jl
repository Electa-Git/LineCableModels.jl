# Opt-in real numerical images under the same managed agent and live authority.
# Passive scientific inputs/projections only; the coordinator imports no solver.
include(joinpath(@__DIR__,"..","..","src","scientific","StudyCases.jl"))

function physical_scientific_profiles()
    definitions = Dict{String,Any}[]
    for (id,variable,operations) in (("line-parameters","LCM_TEST_LINE_IMAGE",["line.frequency_scan"]),
            ("power-flow","LCM_TEST_POWER_IMAGE",["powerflow.prepare","impedance.evaluate"]))
        haskey(ENV,variable) || continue
        image = ENV[variable]
        push!(definitions,Dict("id"=>id,"kind"=>"scientific","isolation"=>"container",
            "environment"=>image,"fingerprint"=>last(split(image,"@sha256:")),"operations"=>operations,
            "budget"=>Dict("cpus"=>1.0,"memory_bytes"=>4*1024^3,"pids"=>256,
                "scratch_bytes"=>128*1024^2,"prepare_seconds"=>600,"job_seconds"=>120)))
    end
    return definitions
end

function physical_scientific_checks(base,client_options,headers,admit,profiles)
    request(method,path,user="alice",body=nothing) = HTTP.request(method,base*path,
        method=="POST" ? [headers(user);"X-LCM-Request"=>"1"] : headers(user),body;
        client_options...,retry=false,status_exception=false,request_timeout=10)
    snapshot(path) = begin
        response = request("GET",path)
        response.status == 200 || error("Scientific query failed with HTTP $(response.status)")
        JSON3.read(response.body,Dict{String,Any})
    end
    results = Dict{String,Any}()
    for profile in profiles
        id = profile["id"]
        case = id=="line-parameters" ? StudyCases.LineParameters() : StudyCases.CorridorImpedance()
        lease_id = admit(id)
        path = "/runtime/api/assignments/$lease_id/science"
        @test request("GET",path,"bob").status == 404
        @test timedwait(()->snapshot(path)["preparation"]=="cold",15;pollint=0.2) == :ok
        prepare = JSON3.write((action="prepare",parameters=StudyCases.preparation_inputs(case),request_id=string(uuid4())))
        @test request("POST",path,"alice",prepare).status == 202
        prior = nothing
        readiness = timedwait(600;pollint=0.5) do
            state = snapshot(path)
            observed = (state["phase"],state["preparation"],state["failure"])
            if observed != prior
                println("Physical science ",id,": ",observed,"; elapsed=",state["elapsed_seconds"])
                flush(stdout)
                prior = observed
            end
            state["failure"] === nothing || error("Physical scientific preparation failed: $(state["failure"])")
            @test request("GET","/health").status == 200
            state["preparation"] == "ready"
        end
        @test readiness==:ok
        readiness==:ok || error("Physical scientific readiness deadline exceeded")
        prepared = snapshot(path)
        @test prepared["executor_id"] !== nothing && prepared["executor_generation"] == 1
        @test request("POST",path,"alice",prepare).status == 202 # exact retry does not replace the executor
        inputs = StudyCases.inputs(case;minimum_frequency_hz=100,maximum_frequency_hz=150,frequency_points=2)
        jobs = "/runtime/api/assignments/$lease_id/jobs"
        submission = JSON3.write((operation=StudyCases.operation(case),parameters=inputs,request_id=string(uuid4())))
        @test request("POST",jobs,"bob",submission).status == 404
        response = request("POST",jobs,"alice",submission)
        @test response.status == 202
        response.status == 202 || error("Physical scientific job was not accepted")
        receipt = JSON3.read(response.body,Dict{String,Any})
        job = "/runtime/api/jobs/$(receipt["id"])"
        completion = timedwait(135;pollint=0.2) do
            state = snapshot(job)["state"]
            state in ("queued","submitted") || state=="succeeded" || error("Physical job ended $state")
            state == "succeeded"
        end
        @test completion==:ok
        completion==:ok || error("Physical scientific job deadline exceeded")
        @test request("GET",job*"/result","bob").status == 404
        outcome = snapshot(job*"/result")
        @test outcome["result"]["result"]["failure"] === nothing
        value = outcome["result"]["result"]["inline_result"]
        series = StudyCases.result_series(case,value)
        @test series.frequency ≈ [100.0,150.0]
        @test all(curve->all(isfinite,curve.values),series.curves)
        @test JSON3.read(request("POST",jobs,"alice",submission).body)["id"] == receipt["id"]
        results[id] = (;inputs,receipt,value)
    end
    @test !any(id.name in ("Bonito","LineCableModelsPlayground","LineCableModels","PowerImpedance") for id in keys(Base.loaded_modules))
    return results
end
