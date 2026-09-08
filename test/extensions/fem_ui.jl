@testitem "Gmsh FEM / UI maps and cancellation preserve session ownership" tags=[:extension,:fem_numerical] begin
    using Gmsh, JSON3
    if isempty(get(ENV,"DISPLAY",""))
        @test_skip "An accessible display is required for UI execution"
    elseif !haskey(ENV,"LINECABLEMODELS_FEM_UI_CASE")
        # FLTK/Gmsh retain native GUI state after finalize. Exercise each closure
        # scenario in a fresh process, including when other tests used the GUI.
        for action in ("before_mesh","before_solve","during_solve","complete")
            code="include(" * repr(joinpath(pkgdir(LineCableModels),"test","runtests.jl")) * ")"
            command=`$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $code fem_ui.jl`
            @test success(addenv(command,"LINECABLEMODELS_FEM_UI_CASE"=>action))
        end
    else
        E=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
        gmsh=Gmsh.gmsh
        copper=Material(kind=:conductor,rho=1.72e-8)
        dielectric=Material(kind=:insulator,rho=1e8,eps_r=2.3)
        design=build(CableDesign,"ui-worker",Stack(
            Group(:core,Region(:metal,Disk(0.005),copper)),
            Region(:insulation,Annulus(0.005,0.01),dielectric),
            Group(:sheath,Region(:shield,Annulus(0.01,0.011),copper))))
        system=build(LineCableSystem,design,(0.0,-0.1);connections=Dict(:core=>1,:sheath=>0))
        problem=LineParametersProblem(system;frequencies=[50.0,1000.0],
            earth_props=homogeneous(rho=100.0,eps_r=10.0))
        form=Formulation(:LineCableModelsFEM;options=(ideal_transposition=false,),
            fem_options=(ui=true,plot_field_maps=true,keep_run_directory=true,
                gmsh_verbosity=0,getdp_verbosity=0))
        function await_condition(predicate)
            deadline=time()+90
            while !predicate()
                time()<deadline || error("UI test timed out")
                sleep(0.005)
            end
        end
        for action in (Symbol(ENV["LINECABLEMODELS_FEM_UI_CASE"]),)
            run_path=Ref("")
            Gmsh.initialize()
            gmsh.model.add("caller-owned-ui-model")
            gmsh.onelab.set_string("caller-owned-parameter",["preserve me"])
            driver=@async try
                await_condition(()->!isempty(gmsh.onelab.get_string(E._onelab_name("ui/run_directory"))))
                run_path[]=only(gmsh.onelab.get_string(E._onelab_name("ui/run_directory")))
                await_condition(()->Bool(gmsh.fltk.is_available()))
                if action !== :before_mesh
                    gmsh.onelab.set_string(E._onelab_name("ui/action"),["generate_mesh"])
                    await_condition(()->gmsh.onelab.get_string(E._onelab_name("ui/mesh_state"))==["ready"])
                end
                if action in (:during_solve,:complete)
                    gmsh.onelab.set_string(E._onelab_name("ui/action"),["run_model"])
                    if action === :during_solve
                        await_condition(()->JSON3.read(read(joinpath(run_path[],"run.json"),String)).getdp_invocations>0)
                    else
                        await_condition(()->gmsh.onelab.get_string(E._onelab_name("ui/solve_state"))==["completed"])
                        @test only(gmsh.onelab.get_number(E._onelab_name("ui/completed_columns")))==4
                        @test only(gmsh.onelab.get_number(E._onelab_name("ui/completed_frequencies")))==2
                        @test length(gmsh.view.get_tags())==36
                    end
                end
                gmsh.fltk.finalize()
            catch
                gmsh.fltk.finalize()
                rethrow()
            end
            result=try
                compute(problem,form;options=(trace=true,))
            catch exception
                exception
            end
            try
                wait(driver)
                @test Bool(gmsh.is_initialized())
                @test gmsh.model.get_current()=="caller-owned-ui-model"
                @test gmsh.onelab.get_string("caller-owned-parameter")==["preserve me"]
                @test isempty(gmsh.view.get_tags())
                if action===:complete
                    @test result isa LineParameters
                    @test length(result.details.fem.run.map_paths)==36
                    @test all(isfile,result.details.fem.run.map_paths)
                else
                    @test result isa LineCableModelsFEMError
                    @test result.category== (action===:during_solve ? :cancelled : :not_executed)
                    state=JSON3.read(read(joinpath(run_path[],"run.json"),String))
                    @test state.state== (action===:during_solve ? "cancelled" : "not_executed")
                    @test E._assert_no_live_attempts(E.FEMRun(run_path[],E.cancelled,"test",:none,""))===nothing
                end
            finally
                Gmsh.finalize()
                isempty(run_path[]) || rm(run_path[];recursive=true,force=true)
            end
        end
    end
end
