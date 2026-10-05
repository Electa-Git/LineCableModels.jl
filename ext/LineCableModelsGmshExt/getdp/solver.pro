// Shared native controls for managed execution and detached ONELAB.
If(!Exists(LinearSolver)) LinearSolver = 0; EndIf // 0: direct MUMPS, 1: GMRES
If(!Exists(GmresIterationsMax)) GmresIterationsMax = 20; EndIf
If(!Exists(GmresRelativeTolerance)) GmresRelativeTolerance = 1.e-12; EndIf
If(!Exists(GmresAbsoluteTolerance)) GmresAbsoluteTolerance = 0.; EndIf
If(!Exists(MumpsOrdering)) MumpsOrdering = -1; EndIf
If(!Exists(PetscPrealloc)) PetscPrealloc = 0; EndIf
If(!Exists(MumpsErrorAnalysis)) MumpsErrorAnalysis = 2; EndIf
If(!Exists(MumpsRefinementMax)) MumpsRefinementMax = 2; EndIf
If(!Exists(MumpsBackwardErrorTolerance)) MumpsBackwardErrorTolerance = 1.e-12; EndIf
If(!Exists(MumpsForwardErrorTolerance)) MumpsForwardErrorTolerance = 0.01; EndIf
If(LinearSolver != 0 && LinearSolver != 1)
  Error("LinearSolver must be 0 (MUMPS) or 1 (GMRES)");
EndIf
If(GmresIterationsMax < 1 || GmresIterationsMax > 2147483647 || Floor[GmresIterationsMax] != GmresIterationsMax)
  Error("GmresIterationsMax must be a positive Cint-representable integer");
EndIf
If(!(GmresRelativeTolerance > 0 && GmresRelativeTolerance <= 1.7976931348623157e308) ||
   !(GmresAbsoluteTolerance >= 0 && GmresAbsoluteTolerance <= 1.7976931348623157e308))
  Error("GMRES tolerances must be finite; relative must be positive and absolute nonnegative");
EndIf
If(MumpsErrorAnalysis != 0 && MumpsErrorAnalysis != 1 && MumpsErrorAnalysis != 2)
  Error("MumpsErrorAnalysis must be 0 (off), 1 (full), or 2 (backward errors)");
EndIf
If(MumpsRefinementMax < 0 || MumpsRefinementMax > 2147483647 || Floor[MumpsRefinementMax] != MumpsRefinementMax)
  Error("MumpsRefinementMax must be a nonnegative Cint-representable integer");
EndIf
If(!(MumpsBackwardErrorTolerance > 0 && MumpsBackwardErrorTolerance <= 1.7976931348623157e308) ||
   !(MumpsForwardErrorTolerance > 0 && MumpsForwardErrorTolerance <= 1.7976931348623157e308))
  Error("MUMPS error tolerances must be finite and positive");
EndIf
// Restore A and b after each solve so original-coordinate residuals and
// subsequent right-hand sides retain their physical units.
FEMSolverOptions = "-ksp_diagonal_scale -ksp_diagonal_scale_fix -pc_type lu -pc_factor_mat_solver_type mumps -ksp_knoll false -ksp_initial_guess_nonzero false -ksp_error_if_not_converged false";
If(LinearSolver == 0)
  FEMSolverOptions = StrCat[FEMSolverOptions, " -ksp_type preonly",
    Sprintf[" -mat_mumps_icntl_11 %g -mat_mumps_icntl_10 %g -mat_mumps_cntl_2 %.17g",
      MumpsErrorAnalysis,MumpsRefinementMax,MumpsBackwardErrorTolerance]];
Else
  // One GMRES solve from zero. Fixed LU preconditioning: no inner refinement
  // or per-application sensitivity analysis. No preliminary direct solve.
  // Native monitoring explicitly recomputes the scaled residual, including
  // at termination. Its difference from the Krylov estimate remains visible.
  FEMSolverOptions = StrCat[FEMSolverOptions,
    " -ksp_type gmres -ksp_pc_side right -ksp_norm_type unpreconditioned -mat_mumps_icntl_10 0 -mat_mumps_icntl_11 0 -ksp_converged_reason -ksp_monitor_true_residual",
    Sprintf[" -ksp_max_it %g -ksp_rtol %.17g -ksp_atol %.17g",
      GmresIterationsMax,GmresRelativeTolerance,GmresAbsoluteTolerance]];
EndIf
If(MumpsOrdering >= 0)
  FEMSolverOptions = StrCat[FEMSolverOptions,Sprintf[" -mat_mumps_icntl_7 %g",MumpsOrdering]];
EndIf
If(PetscPrealloc > 0)
  FEMSolverOptions = StrCat[FEMSolverOptions,Sprintf[" -petsc_prealloc %g",PetscPrealloc]];
EndIf
If(LinearSolver == 0 && MumpsErrorAnalysis)
  // PETSc prints the actual estimates after every RHS, even with GetDP -v 0.
  // GetDP 3.5 does not expose RINFOG values to runtime .pro expressions.
  FEMSolverOptions = StrCat[FEMSolverOptions," -ksp_view"];
EndIf
