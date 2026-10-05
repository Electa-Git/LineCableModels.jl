@testitem "FEM / native solver diagnostics and warning targets" tags=[:extension] begin
    using Gmsh
    E = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    controls(; kwargs...) = E.computation_options(E.LineCableModelsFEM,
        ComputationOptions(; kwargs...)).data
    c = controls()
    @test (c.mumps_error_analysis, c.mumps_refinement_max,
        c.mumps_backward_error_tolerance, c.mumps_forward_error_tolerance) ==
        (2, 2, 1e-12, 0.01)
    @test (c.linear_solver, c.gmres_iterations_max,
        c.gmres_relative_tolerance, c.gmres_absolute_tolerance) == (:mumps, 20, 1e-12, 0.0)
    for (key, values) in (
        :mumps_error_analysis => (-1, 3, true, :full),
        :mumps_refinement_max => (-2, 1.5, true, big(typemax(Cint)) + 1),
        :mumps_backward_error_tolerance => (0, -1, Inf, NaN, true),
        :mumps_forward_error_tolerance => (0, -1, Inf, NaN, true),
        :linear_solver => (:cg, "gmres", true),
        :gmres_iterations_max => (0, -1, 1.5, true, big(typemax(Cint)) + 1),
        :gmres_relative_tolerance => (0, -1, Inf, NaN, true),
        :gmres_absolute_tolerance => (-1, Inf, NaN, true))
        for value in values
            @test_throws ArgumentError controls(; key => value)
        end
    end
    @test controls(mumps_error_analysis=0, mumps_refinement_max=0).mumps_refinement_max == 0
    @test controls(linear_solver=:gmres, gmres_absolute_tolerance=1e-10).linear_solver === :gmres
    @test_throws ArgumentError controls(mumps_refinement_steps=2)

    # Native PETSc 3.14 MUMPS view: two RHSs reuse one factorization. Unused
    # condition and forward-error fields are printed as zero by backward-only analysis.
    native(mode, omega1, omega2, forward, steps) = """
    KSP Object: 1 MPI processes
      ICNTL(11) (error analysis): $mode
      RINFOG(6) (inf norm of residual): 3.0573e-23
      RINFOG(7),RINFOG(8) (backward error est): $omega1, $omega2
      RINFOG(9) (error estimate): $forward
      RINFOG(10),RINFOG(11)(condition numbers): 2.45993e+08, 2.08798e+20
      INFOG(15) (number of steps of iterative refinement after solution): $steps
    """
    ds = E._solver_diagnostics(native(2, 1e-13, 1e-15, 0, 0) *
        native(2, 2e-12, 1e-15, 0, 2), Val(:mumps))
    @test length(ds) == 2
    @test ds[1].backward_error ≈ 1.01e-13
    @test ds[1].scaled_residual == 3.0573e-23
    @test ds[1].forward_error === ds[1].cond1 === ds[1].cond2 === nothing
    @test ds[2].refinement_steps == 2
    @test_logs E._warn_solver_diagnostics(ds[1], c, Val(:mumps); frequency_hz=0.1, basis=1, log="native.log")
    @test_logs (:warn, r"backward-error target not met") E._warn_solver_diagnostics(
        ds[2], c, Val(:mumps); frequency_hz=0.1, basis=2, log="native.log")
    full = only(E._solver_diagnostics(native(1, 1e-14, 1e-16, 0.4, 1), Val(:mumps)))
    @test full.cond1 == 2.45993e8
    @test full.forward_error == 0.4
    @test_logs (:warn, r"scaled-solution sensitivity exceeds budget") E._warn_solver_diagnostics(
        full, controls(mumps_error_analysis=1), Val(:mumps); frequency_hz=0.1, basis=1, log="native.log")
    bad = only(E._solver_diagnostics(native(2, "nan", "inf", 0, 2), Val(:mumps)))
    @test !isfinite(bad.backward_error)
    @test_logs (:warn, r"backward-error target not met") E._warn_solver_diagnostics(
        bad, c, Val(:mumps); frequency_hz=0.1, basis=1, log="native.log")
    missing = only(E._solver_diagnostics("KSP Object: 1 MPI processes\nother solver\n", Val(:mumps)))
    @test_logs (:warn, r"diagnostics unavailable") E._warn_solver_diagnostics(
        missing, c, Val(:mumps); frequency_hz=0.1, basis=1, log="native.log")
    @test isempty(E._solver_diagnostics("no completed solve", Val(:mumps)))

    gmres(reason, n, residual, relative) = """
      0 KSP unpreconditioned resid norm 1 true resid norm 1 ||r(i)||/||b|| 1
      $n KSP unpreconditioned resid norm 1e-16 true resid norm $residual ||r(i)||/||b|| $relative
    Linear solve $(startswith(reason, "CONVERGED_") ? "converged" : "did not converge") due to $reason iterations $n
    KSP Object: 1 MPI processes
      type: gmres
      ICNTL(10) (max num of refinements): 0
    FEM algebraic residual: frequency=1 basis=1 residual=1e-8 rhs=1.414
    """
    gs = E._solver_diagnostics(gmres("CONVERGED_RTOL",2,1e-13,1e-13) *
        gmres("DIVERGED_ITS",1,1e-8,1e-8), Val(:gmres))
    @test getproperty.(gs, :iterations) == [2,1]
    @test gs[1].scaled_relative_residual == 1e-13
    @test gs[1].original_residual_norm == 1e-8
    @test_logs E._warn_solver_diagnostics(gs[1],c,Val(:gmres);frequency_hz=0.1,basis=1,log="native.log")
    @test_logs (:warn,r"GMRES convergence or recomputed residual target not met") E._warn_solver_diagnostics(
        gs[2],c,Val(:gmres);frequency_hz=0.1,basis=2,log="native.log")
    # Native convergence can coexist with a larger explicitly recomputed residual.
    gap = only(E._solver_diagnostics(gmres("CONVERGED_RTOL",2,1e-8,1e-8),Val(:gmres)))
    @test_logs (:warn,r"recomputed residual target not met") E._warn_solver_diagnostics(
        gap,c,Val(:gmres);frequency_hz=0.1,basis=1,log="native.log")
    @test_logs E._warn_solver_diagnostics(gap,controls(gmres_absolute_tolerance=1e-7),Val(:gmres);
        frequency_hz=0.1,basis=1,log="native.log")
    quiet = replace(gmres("CONVERGED_RTOL",2,1e-13,1e-13),
        r"(?m)^FEM algebraic residual:.*\n" => "")
    quiet_records = E._solver_diagnostics(quiet * quiet, Val(:gmres))
    @test length(quiet_records) == 2
    @test all(x -> x.original_residual_norm === nothing && x.iterations == 2, quiet_records)
    missing_gmres = only(E._solver_diagnostics(gmres("CONVERGED_RTOL",2,"unknown","unknown"),Val(:gmres)))
    @test_logs (:warn,r"GMRES diagnostics unavailable") E._warn_solver_diagnostics(
        missing_gmres,c,Val(:gmres);frequency_hz=0.1,basis=1,log="native.log")
    @test isempty(E._solver_diagnostics("no GMRES monitor output",Val(:gmres)))

    # Retain only the attempts supplying each adopted column, including a scan
    # recovered from different attempts. Diagnostic collection must not mutate
    # its cached per-attempt records or lose the second excitation.
    mktempdir() do root
        model = (problem=(frequencies=[0.1],), terminal_ids=["one", "two"])
        for (basis, directory) in ((1, "attempts/old"), (2, "attempts/new"))
            path = joinpath(root, directory)
            mkpath(path)
            write(joinpath(path, "getdp.log"), native(2, basis*1e-13, 0, 0, basis-1))
            E._write_json_atomic(joinpath(path, "attempt.json"),
                (frequency_hz=0.1, requested_bases=[basis]))
            checkpoint = E._column_paths(root, 1, basis, false).checkpoint
            E._write_json_atomic(checkpoint, (attempt=directory,))
        end
        retained = E._retained_solver_diagnostics(root, model, Val(:mumps))
        @test length(retained) == 2
        @test getproperty.(retained, :basis) == [1, 2]
        @test getproperty.(retained, :frequency_index) == [1, 1]
        @test retained[2].backward_error == 2e-13
        for dir in ("old", "new")
            write(joinpath(root,"attempts",dir,"getdp.log"),gmres("CONVERGED_RTOL",2,1e-13,1e-13))
        end
        retained = E._retained_solver_diagnostics(root,model,Val(:gmres))
        @test length(retained) == 2
        @test all(x -> x.linear_solver === :gmres && x.iterations == 2, retained)
    end
end
