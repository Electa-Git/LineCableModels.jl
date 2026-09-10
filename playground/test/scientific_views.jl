@testset "shared scientific cases and session-local views" begin
    LCM = LineCableModelsPlayground
    S = LCM.ScientificViews
    C = S.StudyCases
    line, corridor = C.LineParameters(), C.CorridorImpedance()
    for case in (line, corridor)
        @test C.inputs(case) == C.inputs(case)
        @test C.inputs(case) !== C.inputs(case)
        @test C.preparation_inputs(case) !== C.preparation_inputs(case)
        @test_throws ArgumentError C.inputs(case; frequency_points=2.5)
        @test_throws ArgumentError C.inputs(case; frequency_points=true)
        @test_throws ArgumentError C.inputs(case; frequency_points=201)
        @test_throws ArgumentError C.inputs(case; minimum_frequency_hz=0)
        @test_throws ArgumentError C.inputs(case; minimum_frequency_hz=2500)
        @test_throws ArgumentError C.inputs(case; maximum_frequency_hz=NaN)
    end
    @test C.inputs(line)["frequencies_hz"][[1,end]] == [1.0,2500.0]
    @test length(C.inputs(line)["frequencies_hz"]) == 40
    @test_throws ArgumentError C.inputs(line; separation_m=0)
    @test_throws ArgumentError C.inputs(corridor; ugc_share=0)
    @test_throws ArgumentError C.inputs(corridor; ugc_share=1)
    @test_throws ArgumentError C.inputs(corridor; length_error_percent=51)
    @test C.inputs(corridor)["specification"] == C.preparation_inputs(corridor)["specification"]
    @test C.inputs(corridor)["prepared_resource_key"] == ""

    wire(real, imag) = Dict{String,Any}("real"=>real,"imag"=>imag)
    tensor = [[[wire(1.0,2.0),wire(2.0,3.0)] for _ in 1:2] for _ in 1:2]
    value = Dict{String,Any}("frequencies_hz"=>[1.0,100.0],
        "series_impedance_ohm_per_m"=>tensor, "shunt_admittance_s_per_m"=>deepcopy(tensor))
    for (quantity,unit,ys) in (("resistance","Ω/m",[1.0,2.0]), ("reactance","Ω/m",[2.0,3.0]),
            ("conductance","S/m",[1.0,2.0]), ("susceptance","S/m",[2.0,3.0]))
        series = C.result_series(line,value,quantity)
        @test series.unit == unit
        @test series.curves[1].values == ys
        @test length(S.plot_coordinates(series).paths) == 2
    end
    @test_throws ArgumentError C.result_series(line,value,"inductance")
    @test_throws ArgumentError C.result_series(line,merge(value,Dict("frequencies_hz"=>[1.0,1.0])))
    invalid = deepcopy(value)
    invalid["shunt_admittance_s_per_m"][2][2][2]["imag"] = Inf
    @test_throws ArgumentError C.result_series(line,invalid,"resistance")
    invalid = deepcopy(value)
    pop!(invalid["series_impedance_ohm_per_m"][1][1])
    @test_throws ArgumentError C.result_series(line,invalid)
    flow = Dict("curves"=>Dict(key=>Dict("frequency_hz"=>[1.0,100.0],"magnitude_db_ohm"=>[1.0,2.0])
        for key in ("base_minus_error","base","base_plus_error")))
    @test length(C.result_series(corridor,flow).curves) == 3
    @test C.result_series(corridor,flow).unit == "dB re 1 Ω"
    flow["curves"]["base"]["frequency_hz"] = [1.0,99.0]
    @test_throws ArgumentError C.result_series(corridor,flow)
    huge = (frequency=[1.0,100.0], curves=[(label="finite",values=[-floatmax(Float64),floatmax(Float64)])],unit="Ω/m")
    @test !occursin(r"Inf|NaN", repr(S.plot_coordinates(huge)))

    client = LCM.RuntimeClient(LCM.uuid4())
    session = LCM.Bonito.Session()
    try
        views = [S.ScientificView(session,case,client) for case in (line,corridor)]
        twin = S.ScientificView(session,line,client)
        for view in views
            for field in values(view.fields)
                c = field.control
                steps = (c.initial-c.minimum)/c.step
                @test isapprox(steps,round(steps); atol=1e-8)
            end
            @test view.job.operation == C.operation(view.case)
            @test view.job.parameters[] == C.inputs(view.case)
            @test view.job.result[] === nothing
            @test LCM.Bonito.jsrender(session,view) !== nothing
            inspection = LCM.ComponentXRay.inspection(view)
            @test isempty(inspection.bindings)
            @test !occursin(string(client.run_id),repr(inspection))
        end
        first(views).fields.separation_m.control.value[] = 0.75
        @test first(views).job.parameters[]["separation_m"] == 0.75
        @test twin.job.parameters[]["separation_m"] == 0.5
        @test first(views).fields.separation_m.control.id != twin.fields.separation_m.control.id
        first(views).valid[] = false
        @test first(views).job.parameters[] === nothing
        first(views).valid[] = true
        first(views).fields.frequency_points.control.value[] = 2.5
        @test first(views).job.parameters[] === nothing
        first(views).fields.frequency_points.control.value[] = 2
        @test length(first(views).job.parameters[]["frequencies_hz"]) == 2
        for component in (S.CableGeometry(),S.StudyRuntime(client))
            @test LCM.Bonito.jsrender(session,component) !== nothing
            @test !isempty(LCM.ComponentXRay.inspection(component).css_scopes)
        end
        # Exercise empty -> populated after jsrender has attached all observers.
        # A NamedTuple inferred with series::Nothing cannot accept real data.
        view = first(views)
        fence = LCM.AssignmentFence(string(LCM.uuid4()),string(client.run_id),"scientific-test",view.job.role,
            "worker-a",string(LCM.uuid4()),string(LCM.uuid4()),"line-parameters","1.0.0",repeat("a",64),1)
        request = LCM.new_job_request(view.job.operation,view.job.parameters[];session_id=fence.run_id)
        result = LCM.JobResult("1.0",request.job_id,request.operation,"1.0",request.input_hash,
            "fixture",fence.fingerprint,fence.worker_id,"miss",LCM.utc_timestamp(),LCM.utc_timestamp(),
            value,nothing,nothing,String[])
        envelope = LCM.AssignedResult("2.0",fence,result,LCM.PreparedExecution(string(LCM.uuid4()),1,repeat("b",64)))
        draft = string(LCM.uuid4())
        packet = LCM.JSON3.write((draft_id=draft,current=true,provenance=envelope,value))
        @test LCM.apply_job_projection!(view.job,packet,draft)
        @test view.job.result[].current
        @test S.display_state(view,view.job.result[],"resistance").series !== nothing
        LCM.invalidate_job_projection!(view.job)
        @test !view.job.result[].current
        @test S.display_state(view,view.job.result[],"resistance").series !== nothing
        @test typeof(S.display_state(view,nothing,"resistance")) ==
            typeof(S.display_state(view,view.job.result[],"resistance"))
        application = LCM.CableStudy.Application(client)
        state = LCM.WorkbenchUI.initialize(application,session)
        shell = LCM.WorkbenchUI.compose(application,state)
        @test shell.namespace == :cable_study
        @test shell.navigation.footer isa LCM.NavigationButton
        @test shell.navigation.footer.label == "Home" && shell.navigation.footer.href == "/"
        dashboard = LCM.NavigationButton("Dashboard"; href="/workbenches/", icon=LCM.WorkbenchUI.icon(:workbench))
        custom = LCM.CableStudy.Application(client; return_button=dashboard)
        @test LCM.WorkbenchUI.compose(custom,state).navigation.footer === dashboard
        @test state.views.parameters isa S.ScientificView{C.LineParameters}
        @test state.views.corridor isa S.ScientificView{C.CorridorImpedance}
        @test state.views.terminal isa LCM.JuliaTerminal
        @test state.views.parameters.job.parameters[] == twin.job.parameters[]
        @test length(shell.workspace.views) == 5
        LCM.WorkbenchUI.handle!(application,state,LCM.CableStudy.SelectView(:parameters))
        @test state.active[] == :parameters
        @test_throws ArgumentError LCM.WorkbenchUI.handle!(application,state,LCM.CableStudy.SelectView(:unknown))
        @test LCM.Bonito.jsrender(session,LCM.WorkbenchUI.Runtime(application,state,shell,LCM.ComponentXRay.XRayPolicy(true))) !== nothing
    finally
        close(session)
    end
    @test length(LCM.Showcase.routes(client)) == 5
    @test all(last(route) isa LCM.Bonito.App for route in LCM.Showcase.routes(client))
    @test !any(id.name in ("LineCableModels","PowerImpedance","LineCableModelsExecutionCore") for id in keys(Base.loaded_modules))
    css = S.STYLES
    @test !occursin(r"#[0-9a-fA-F]{3,8}\b",css)
    @test !occursin(r"--lc-[a-z-]+\s*:",css)
    root = LCM.PLAYGROUND_ROOT
    deck = read(joinpath(root,"presentations","showcase.qmd"),String)
    @test count("requires-run=\"true\"",deck) == 5
    @test Set(m.captures[1] for m in eachmatch(r"route=\"([^\"]+)\"",deck)) == Set(first.(LCM.Showcase.routes(client)))
end
