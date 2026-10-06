// Included by model.pro when the Physics selector is 1 (Helmholtz).
Jacobian { { Name Plain; Case { { Region All; Jacobian Vol; } } } }



// Maxwell formulation for exp(j omega t - Gamma z), with prescribed Gamma.
// A_z=a, A_t=Gamma*bt, phi=Gamma*v is exact for nonzero Gamma when all
// Gamma^2 terms are retained. At exactly zero, the normalized limit remains
// regular without division by Gamma. Gamma^2 means the complex square.
// In finite metal, ur stores w = u - Gamma^2 V, so Ez = -j omega a - w.
// Cancel the identical metal u/V basis contributions before assembly. keeping
// them separate loses the small driven current through large-term cancellation.
// The physical drive u = w + Gamma^2 V is restored in K. At Gamma=0, w=u.
// PerfectConductors=1 instead uses exact PEC contour current constraints.
//
// One imposed axial current supplies both the magnetic equation and the
// normalized transverse leakage. There is no independent electric drive.
// The transverse domain treats metal contours as equipotential electrodes.
// The magnetic block either retains the specified finite metal conductivity
// or excludes metal interiors when PerfectConductors=1.
//
// Gauge-invariant voltage/Gamma = integral((grad v+j omega bt).dl).
// Independent, horizontally shifted sampling lines run upward from the same
// reference to the receiver's own metal. Stored fields exclude metal interiors.
// No endpoint potentials or extra stretch factor enter the measurement.
// Raw P has units ohm m (inverse admittance). analytical Pe = j omega P.
// The longitudinal drive K=-(U+Gamma^2 V)/I (U aliases w) and series coefficient obey
// Z=K+Gamma^2*P. For PEC contours K=(j omega A-Gamma^2 V)/I.
//
// Reference: G. Ciuprina and R. V. Sabriego, "Electric circuit element boundary
// conditions for electromagneto-quasistatic and full wave models in A, phi
// potentials and their finite element implementation", J. Math. Ind. 14, 27
// (2024), https://doi.org/10.1186/s13362-024-00165-6, Sec. 4 and Appendix B.
// The paper gives the full-vector A/phi equations, terminal conditions and
// gauge framework. This 2D longitudinal substitution and the vertical voltage
// paths are our reduction, not the paper's 3D ECE model or its Darwin model.
// See the LineCableModelsFEM Julia docstring for equations, units and limits.

// PerfectConductors=1 excludes metal interiors and imposes exact PEC
// current-carrying contours. Default 0 retains the working finite-metal a/u block.
If(!Exists(PerfectConductors))
  PerfectConductors = 0;
EndIf
If(!Exists(RunDirectory))
  RunDirectory = "";
EndIf
If(!Exists(FrequencyCount))
  FrequencyCount = 0;
EndIf
If(!Exists(FrequencyIndex))
  FrequencyIndex = 1;
EndIf
If(!Exists(FrequencyHz))
  FrequencyHz = 1.;
EndIf
If(!Exists(GammaRe)) GammaRe = 0.; EndIf
If(!Exists(GammaIm)) GammaIm = 0.; EndIf
FiniteGamma = GammaRe != 0. || GammaIm != 0.;
If(!Exists(BasisTerminal))
  BasisTerminal = 1;
EndIf
If(Exists(BasisListPath))
  Include BasisListPath;
ElseIf(!Exists(RequestedBases))
  RequestedBases() = {BasisTerminal};
EndIf
If(!Exists(ReuseFactorization))
  ReuseFactorization = 1;
EndIf
If(!Exists(PlotFieldMaps))
  PlotFieldMaps = 0;
EndIf
If(!Exists(RawDirectory)) RawDirectory = StrCat[RunDirectory, "/raw"]; EndIf
RawJobDirectory = StrCat[RawDirectory, "/jobs"];
MapDirectory = StrCat[RunDirectory, "/maps"];

Group {
  Air = Region[{AIR_EM}];
  Earth = Region[{EARTH_EM}];
  AirPml = Region[{AIR_PML}];
  EarthPml = Region[{EARTH_PML}];
  Terminals = Region[{}];
  TerminalContours = Region[{}];
  MeasurementLines = Region[{}];

  For t In {1:NumTerminals}
    Terminal~{t} = Region[{(TERMINAL + t - 1)}];
    TerminalContour~{t} = Region[{(TERMINAL_CONTOUR + t - 1)}];
    MeasurementLine~{t} = Region[{(MEASUREMENT_LINE + t - 1)}];
    MeasurementLines += Region[{MeasurementLine~{t}}];
    Terminals += Region[{Terminal~{t}}];
    TerminalContours += Region[{TerminalContour~{t}}];
  EndFor
  Sur_Dirichlet_Mag = Region[{OUTBND_EM}];
  Sur_Dirichlet_Ele = Region[{OUTBND_EM}];
  // Complete the gauge tree on every conductor and outer boundary with prescribed bt circulation.
  // Omitting the electrode contours would constrain physical loop degrees of
  // freedom also to the potential gauge.
  GaugeBoundary = Region[{Sur_Dirichlet_Mag, TerminalContours}];
}

Jacobian { { Name VoltageLine; Case { { Region All; Jacobian Sur; } } } }

Include "materials.pro";
Group {
  // Conductors act as equipotential electrodes. Electric fields use the media outside metal.
  DomainMedia_Ele = Region[{Air, AirPml, Earth, EarthPml, PassiveMaterialRegions}];
  If(PerfectConductors)
    DomainFields = Region[{DomainMedia_Ele}];
    DomainFieldLoss = Region[{Earth, EarthPml, LossyMaterialRegions}];
    If(AirSigma != 0)
      DomainFieldLoss += Region[{Air, AirPml}];
    EndIf
  Else
    DomainFields = Region[{Domain_Mag}];
    DomainFieldLoss = Region[{DomainLoss}];
  EndIf
}

Function {
  material_tag[Air] = AIR_EM; material_tag[Earth] = EARTH_EM;
  material_tag[AirPml] = AIR_PML; material_tag[EarthPml] = EARTH_PML;
  nu[#{Air, AirPml}] = 1. / AirMu;
  sigma[#{Air, AirPml}] = AirSigma;
  epsilon[#{Air, AirPml}] = AirEpsilon;
  mu[#{Air, AirPml}] = AirMu;

  nu[#{Earth, EarthPml}] = 1. / EarthMu(FrequencyIndex - 1);
  sigma[#{Earth, EarthPml}] = EarthSigma(FrequencyIndex - 1);
  epsilon[#{Earth, EarthPml}] = EarthEpsilon(FrequencyIndex - 1);
  mu[#{Earth, EarthPml}] = EarthMu(FrequencyIndex - 1);

  omega[] = 2. * Pi * $FEMFrequencyHz;
  se[] = Complex[sigma[], omega[] * epsilon[]];
  gamma[] = Complex[GammaRe,GammaIm];
  gamma2[] = gamma[] * gamma[];
  // Axial field reconstruction after the exact metal-drive substitution.
  gamma2Ez[DomainMedia_Ele] = gamma2[];
  gamma2Ez[ConductorMaterialRegions] = 0.;
}

Include "pml.pro";

Constraint {
  { Name FEMMagneticVectorPotential;
    Case {
      { Region Sur_Dirichlet_Mag; Value 0.; }
    }
  }
  { Name FEMTerminalCurrent;
    Case {
      For t In {1:NumTerminals}
        { Region Terminal~{t}; Value $FEM_I~{t}; }
      EndFor
    }
  }
  { Name FEMTransverseBoundary; Type Assign;
    Case { { Region Sur_Dirichlet_Mag; Value 0.; }
           { Region TerminalContours; Value 0.; } }
  }
  { Name FEMTransverseGauge; Type Assign;
    Case { { Region DomainMedia_Ele; SubRegion GaugeBoundary; Value 0.; } }
  }
  { Name FEMScalarPotential; Type Assign;
    Case {
      // Both air-side and earth-side outer boundaries use vanishing scalar potential.
      { Region Sur_Dirichlet_Ele; Value 0.; }
    }
  }
}

FunctionSpace {
  If(PerfectConductors)
    { Name Hcurl_a_FEM_2D; Type Form1P;
      BasisFunction {
        { Name se; NameOfCoef ae; Function BF_PerpendicularEdge;
          Support Domain_Mag; Entity NodesOf[All, Not Terminals]; }
        { Name sf; NameOfCoef af; Function BF_GroupOfPerpendicularEdges;
          Support Domain_Mag; Entity GroupsOfNodesOf[Terminals]; }
      }
      GlobalQuantity {
        { Name A; Type AliasOf; NameOfCoef af; }
        { Name I; Type AssociatedWith; NameOfCoef af; }
      }
      Constraint {
        { NameOfCoef ae; EntityType NodesOf; NameOfConstraint FEMMagneticVectorPotential; }
        { NameOfCoef I; EntityType Auto; NameOfConstraint FEMTerminalCurrent; }
      }
    }
  Else
  { Name Hcurl_a_FEM_2D; Type Form1P;
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_PerpendicularEdge;
        Support Domain_Mag; Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef ae; EntityType NodesOf;
        NameOfConstraint FEMMagneticVectorPotential; }
    }
  }

  { Name Hregion_u_FEM_2D; Type Form1P;
    BasisFunction {
      { Name sr; NameOfCoef ur; Function BF_GroupOfPerpendicularEdges;
        Support ConductorMaterialRegions;
        Entity GroupsOfNodesOf[Terminals]; }
    }
    GlobalQuantity {
      { Name U; Type AliasOf; NameOfCoef ur; }
      { Name I; Type AssociatedWith; NameOfCoef ur; }
    }
    Constraint {
      { NameOfCoef I; EntityType Auto; NameOfConstraint FEMTerminalCurrent; }
    }
  }

  EndIf
  { Name Hcurl_bt_FEM_2D; Type Form1;
    BasisFunction {
      { Name se; NameOfCoef be; Function BF_Edge;
        Support DomainMedia_Ele; Entity EdgesOf[All]; }
    }
    Constraint {
      { NameOfCoef be; EntityType EdgesOf;
        NameOfConstraint FEMTransverseBoundary; }
      { NameOfCoef be; EntityType EdgesOfTreeIn; EntitySubType StartingOn;
        NameOfConstraint FEMTransverseGauge; }
    }
  }
  { Name Hgrad_v_FEM_2D; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Domain_Mag; Entity NodesOf[All, Not Terminals]; }
      { Name sf; NameOfCoef vf; Function BF_GroupOfNodes;
        Support Domain_Mag; Entity GroupsOfNodesOf[Terminals]; }
    }
    GlobalQuantity {
      { Name V; Type AliasOf; NameOfCoef vf; }

    }
    Constraint {

      { NameOfCoef vn; EntityType NodesOf;
        NameOfConstraint FEMScalarPotential; }
    }
  }
}

Formulation {
  { Name FEM_Helmholtz_2D; Type FemEquation;
    Quantity {
      { Name a; Type Local; NameOfSpace Hcurl_a_FEM_2D; }
      If(PerfectConductors)
        { Name A; Type Global; NameOfSpace Hcurl_a_FEM_2D [A]; }
        { Name I; Type Global; NameOfSpace Hcurl_a_FEM_2D [I]; }
      Else
      { Name ur; Type Local; NameOfSpace Hregion_u_FEM_2D; }
      { Name U; Type Global; NameOfSpace Hregion_u_FEM_2D [U]; }
      { Name I; Type Global; NameOfSpace Hregion_u_FEM_2D [I]; }
      EndIf
      { Name bt; Type Local; NameOfSpace Hcurl_bt_FEM_2D; }
      { Name v; Type Local; NameOfSpace Hgrad_v_FEM_2D; }
      { Name V; Type Global; NameOfSpace Hgrad_v_FEM_2D [V]; }

    }
    Equation {
      If(PerfectConductors)
        Galerkin { [nuPml[] * Dof{d a}, {d a}];
          In DomainMedia_Ele; Jacobian Vol; Integration I1; }
        Galerkin { [Complex[0,omega[]] * sePml[] * Dof{a}, {a}];
          In DomainMedia_Ele; Jacobian Vol; Integration I1; }
        GlobalTerm { [-Dof{I}, {A}]; In Terminals; }
      Else
      Galerkin {
        [nuPml[] * Dof{d a}, {d a}];
        In Domain_Mag; Jacobian Vol; Integration I1;
      }
      // exp(+j omega t): j omega sigma - omega^2 epsilon = j omega se.
      // Sigma vanishes outside DomainLoss. ur is supported only in metal.
      // Consolidate the identical trial and test pairs before harmonic assembly.
      Galerkin { [Complex[0,omega[]] * seZ[] * Dof{a}, {a}];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      Galerkin { [seZ[] * Dof{ur}, {a}];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0,omega[]] * seZ[] * Dof{a}, {ur}];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      Galerkin { [seZ[] * Dof{ur}, {ur}];
        In Domain_Mag; Jacobian Vol; Integration I1; }

        GlobalTerm { [Dof{I}, {U}]; In DomainCWithI; }
      EndIf
      If(FiniteGamma)
        // B_t = curl_t(a zhat) - Gamma^2 zhat x bt.
        Galerkin { [-gamma2[] * nuPml[] * (Vector[0,0,1] /\ Dof{bt}), {d a}];
          In DomainMedia_Ele; Jacobian Vol; Integration I1; }
        // In media E_z = -j omega a + Gamma^2 v. In finite metal the
        // Gamma^2 V contribution is already absorbed into the unknown ur.
        Galerkin { [-gamma2[] * seZ[] * Dof{v} * Vector[0,0,1], {a}];
          In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      EndIf
      // Normalized transverse Ampere uses `A_t = Gamma * bt` and `phi = Gamma * v`.
      // `C*(nu C bt) + (j omega se-Gamma^2 nu) bt + se grad(v) - nu grad(a) = 0`
      // in isotropic physical media. PML needs the explicit tensor rotations.
      // The tree gauge removes gradient degrees of freedom from bt. The
      // continuity rows below supply that part of Ampere in the nodal space.
      Galerkin { [nuPml[] * Dof{d bt}, {d bt}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0,omega[]] * sePml[] * Dof{bt}, {bt}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      Galerkin { [sePml[] * Dof{d v}, {bt}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      // d a = (d_y a, -d_x a, 0). rotating restores grad_t(a).
      Galerkin { [-(Vector[0,0,1] /\ (nuPml[] * Dof{d a})), {bt}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      If(FiniteGamma)
        Galerkin { [gamma2[] * (Vector[0,0,1] /\ (nuPml[] * (Vector[0,0,1] /\ Dof{bt}))), {bt}];
          In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      EndIf
      // Normalized continuity. The source is the axial current, not a second drive.
      Galerkin { [sePml[] * Dof{d v}, {d v}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0,omega[]] * sePml[] * Dof{bt}, {d v}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0,omega[]] * seZ[] * (Dof{a} * Vector[0,0,1]), {v}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      If(FiniteGamma)
        Galerkin { [-gamma2[] * seZ[] * Dof{v}, {v}];
          In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      EndIf
      GlobalTerm { [-Dof{I}, {V}]; In Terminals; }
    }
  }
}

Macro FEMSetBasisCurrent
For t In {1:NumTerminals}
  Evaluate[$FEM_I~{t} = Complex[
    UnitSource * ($FEMBasisTerminal == t), 0.]];
EndFor
Return

Macro FEMSolveBasis
Evaluate[$FEMConstraintStart = GetWallClockTime[]];
Call FEMSetBasisCurrent;
UpdateConstraint[Sys_FEM];
Evaluate[$FEMAssemblyStart = GetWallClockTime[]];
Test[$FEMFirstSolve || !ReuseFactorization]{
  Generate[Sys_FEM];
  Evaluate[$FEMSolveStart = GetWallClockTime[]];
  Solve[Sys_FEM];
}{
  GenerateRHSGroup[Sys_FEM, Terminals];
  Evaluate[$FEMSolveStart = GetWallClockTime[]];
  SolveAgain[Sys_FEM];
}
GetResidual[Sys_FEM, $FEMResidualNorm];
GetNormRightHandSide[Sys_FEM, $FEMRHSNorm];
Print[{$FEMFrequencyIndex, $FEMBasisTerminal, $FEMResidualNorm, $FEMRHSNorm},
  Format "FEM algebraic residual: frequency=%g basis=%g residual=%.17g rhs=%.17g"];
Evaluate[$FEMOutputStart = GetWallClockTime[]];
Return

Macro FEMFieldSystem
{ Name Sys_FEM; NameOfFormulation FEM_Helmholtz_2D; Type Complex; Frequency 1.; }
Return

Macro FEMScan
      SetGlobalSolverOptions[FEMSolverOptions];
      If(LinearSolver == 0 && MumpsErrorAnalysis)
        Print[{MumpsBackwardErrorTolerance,FEMMumpsForwardErrorBudget},
          Format "MUMPS targets: backward error %.6g; scaled-solution forward estimate %.6g (full analysis only; not terminal accuracy)"];
      EndIf
      If(LinearSolver == 1)
        Print[{GmresIterationsMax,GmresRelativeTolerance,GmresAbsoluteTolerance},
          Format "GMRES controls: maximum iterations %g; scaled residual targets relative %.17g absolute %.17g"];
      EndIf
      CreateDir[RawDirectory];
      CreateDir[RawJobDirectory];
      If(PlotFieldMaps) CreateDir[MapDirectory]; EndIf
      // Native observations accompany every frequency's retained columns.
      Print[{FrequencyIndex,FrequencyHz,FEMPmlTarget,FEMPmlFloorActive~{0},FEMPmlFloorActive~{1},FEMCutoff~{0},FEMCutoff~{1},FEMPmlSideNetAttenuation~{0},FEMPmlSideNetAttenuation~{1},FEMPmlTopNetAttenuation,FEMPmlBottomNetAttenuation,FEMPmlFloorUsed,FEMPmlSideLayers,FEMPmlTopLayers,FEMPmlBottomLayers,PmlSideEta,PmlTopEta,PmlBottomEta,FEMEarthSizingCeilingActive,FEMEarthLayerThickness,FEMEarthLayerClippedOrOmitted,FEMCapActive(0),FEMCapActive(1),FEMCapActive(2)},
        Format "%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g	%.17g",
        File StrCat[RawJobDirectory,Sprintf["/pml-f%04g.tsv",FrequencyIndex]]];
      If(FEMEarthSizingCeilingActive)
        Print[{FrequencyHz},Format "FEM warning: earth too resistive for FEM domain sizing; results not qualified (f=%.17g Hz)"];
      EndIf
      If(FEMCutoff~{0} || FEMCutoff~{1})
        Print[{FrequencyHz},Format "FEM warning: exact transverse cutoff; the transverse problem is singular (f=%.17g Hz)"];
      EndIf
      If(FEMPmlFloorUsed)
        Print[{FrequencyHz},Format "PML sizing floor active (informational; f=%.17g Hz)"];
      EndIf
      InitSolution[Sys_FEM];
      Evaluate[$FEMFrequencyIndex = FrequencyIndex];
      Evaluate[$FEMFrequencyHz = FrequencyHz];
      SetFrequency[Sys_FEM, FrequencyHz];
      SetTimeStep[1];
      // GetDP time is an output-step identity, not the physical frequency.
      // The scan index prevents failures due to a time range. SetFrequency owns physics.
      SetTime[FrequencyIndex];
      Evaluate[$FEMFirstSolve = 1];
      For basis_slot In {0:#RequestedBases()-1}
        basis = RequestedBases(basis_slot);
        Evaluate[$FEMBasisTerminal = basis];
        Call FEMSolveBasis;
        PostOperation[FEMAppendRaw~{basis}];
        If(PlotFieldMaps)
          PostOperation[FEMWriteMaps~{basis}];
        EndIf
        Evaluate[$FEMOutputEnd = GetWallClockTime[]];
        PostOperation[FEMCompleteColumn~{basis}];
        Evaluate[$FEMFirstSolve = 0];
      EndFor
Return

If(!Exists(ExportLineParameters))
Resolution {
  { Name LineCableModelsFEMScan;
    System { Call FEMFieldSystem; }
    Operation { Call FEMScan; }
  }
}
EndIf

PostProcessing {
  { Name FEMFields; NameOfFormulation FEM_Helmholtz_2D; NameOfSystem Sys_FEM;
    PostQuantity {
      { Name bt_mesh; Value { Term { [{bt}]; In DomainMedia_Ele; Jacobian Plain; } } }
      { Name voltage_gradient; Value { Term {
        [{d v}+Complex[0,omega[]]*{bt}]; In DomainMedia_Ele; Jacobian Plain;
      } } }
      For field_slot In {0:#RequestedBases()-1}
        field_basis = RequestedBases(field_slot);
        { Name ReVoltageLine~{field_basis}; Value { Integral {
          [Re[CompY[ComplexVectorField[XYZ[]]{100000+field_basis}]]];
          In MeasurementLines; Jacobian VoltageLine; Integration I2;
        } } }
        { Name ImVoltageLine~{field_basis}; Value { Integral {
          [Im[CompY[ComplexVectorField[XYZ[]]{100000+field_basis}]]];
          In MeasurementLines; Jacobian VoltageLine; Integration I2;
        } } }
      EndFor
      { Name v_local; Value { Term { [{v}]; In DomainMedia_Ele; Jacobian Plain; } } }
      { Name hz_scaled; Value { Term { [CompZ[nuPml[]*{d bt}]]; In DomainMedia_Ele; Jacobian Vol; } } }

      { Name material_region; Value { Term { [material_tag[]]; In DomainFields; Jacobian Vol; } } }
      { Name pml_mask; Value { Term { [pmlX[] > 0. || pmlY[] > 0.]; In DomainFields; Jacobian Vol; } } }
      { Name jt_mesh; Value { Term { [-sePml[] * ({d v}+Complex[0,omega[]]*{bt})]; In DomainMedia_Ele; Jacobian Vol; } } }
      { Name jt; Value { Term { [-se[] * pmlInvS[] * ({d v}+Complex[0,omega[]]*{bt})]; In DomainMedia_Ele; Jacobian Vol; } } }
      { Name az; Value {
        Term { [CompZ[{a}]]; In DomainFields; Jacobian Vol; }
      }}
      { Name b; Value {
        Term { [pmlBToPhysical[] * ({d a}-gamma2[]*(Vector[0,0,1] /\ {bt})+gamma[]*{d bt})]; In DomainFields; Jacobian Vol; }
      }}
      { Name bm; Value {
        Term { [Norm[pmlBToPhysical[] * ({d a}-gamma2[]*(Vector[0,0,1] /\ {bt})+gamma[]*{d bt})]]; In DomainFields; Jacobian Vol; }
      }}
      // e is E_t/Gamma [V]. ez is the driven axial field [V/m].
      { Name e; Value {
        Term { [pmlInvS[] * (-{d v}-Complex[0,omega[]]*{bt})]; In DomainMedia_Ele; Jacobian Vol; }
      }}
      If(PerfectConductors)
        { Name ez; Value { Term { [-CompZ[Dt[{a}]]+gamma2[]*{v}]; In DomainFields; Jacobian Vol; } } }
      Else
      { Name ez; Value {
        Term { [-CompZ[Dt[{a}] + {ur}]+gamma2Ez[]*{v}]; In DomainFields; Jacobian Vol; }
      }}
      EndIf
      { Name em; Value {
        Term { [Norm[pmlInvS[] * (-{d v}-Complex[0,omega[]]*{bt})]]; In DomainMedia_Ele; Jacobian Vol; }
      }}
      If(PerfectConductors)
        { Name jz; Value { Term { [se[]*(-CompZ[Dt[{a}]]+gamma2[]*{v})]; In DomainFields; Jacobian Vol; } } }
      Else
      { Name jz; Value {
        Term { [se[]*(-CompZ[Dt[{a}] + {ur}]+gamma2Ez[]*{v})];
          In DomainFields; Jacobian Vol; }
      }}
      EndIf
      { Name jm; Value {
        Term { [Norm[-se[] * pmlInvS[] * ({d v}+Complex[0,omega[]]*{bt})]]; In DomainMedia_Ele; Jacobian Vol; }
      }}
      If(PerfectConductors)
        { Name rhoj2; Value { Term { [0.5 * sigma[] * (SquNorm[CompZ[Dt[{a}]]-gamma2[]*{v}]+SquNorm[gamma[]*pmlInvS[]*({d v}+Complex[0,omega[]]*{bt})])]; In DomainFieldLoss; Jacobian Vol; } } }
      Else
      { Name rhoj2; Value {
        Term { [0.5 * sigma[] * (SquNorm[CompZ[Dt[{a}] + {ur}]-gamma2Ez[]*{v}]+SquNorm[gamma[]*pmlInvS[]*({d v}+Complex[0,omega[]]*{bt})])];
          In DomainFieldLoss; Jacobian Vol; }
      }}
      EndIf
      If(PerfectConductors)
        { Name ReZ; Value { Term { [Re[(Complex[0,omega[]]*{A}-gamma2[]*{V})/UnitSource]]; In Terminals; } } }
        { Name ImZ; Value { Term { [Im[(Complex[0,omega[]]*{A}-gamma2[]*{V})/UnitSource]]; In Terminals; } } }
      Else
      { Name ReZ; Value {
        Term { [-Re[({U}+gamma2[]*{V}) / UnitSource]]; In DomainCWithI; }
      }}
      { Name ImZ; Value {
        Term { [-Im[({U}+gamma2[]*{V}) / UnitSource]]; In DomainCWithI; }
      }}
      EndIf
      { Name ReP; Value {
        Term { [Re[{V} / UnitSource]]; In Terminals; }
      }}
      { Name ImP; Value {
        Term { [Im[{V} / UnitSource]]; In Terminals; }
      }}
    }
  }
}

// Expand output declarations only for the requested terminal columns.
For basis_slot In {0:#RequestedBases()-1}
  basis = RequestedBases(basis_slot);
  ColumnStem = Sprintf["getdp-f%04.0f-b%04.0f", FrequencyIndex, basis];
  RawZPath = StrCat[RawJobDirectory, "/", ColumnStem, "-Z.tsv"];
  RawScalarPath = StrCat[RawJobDirectory, "/", ColumnStem, "-Pscalar.tsv"];
  RawPPath = StrCat[RawJobDirectory, "/", ColumnStem, "-P.tsv"];
  TimingPath = StrCat[RawJobDirectory, "/", ColumnStem, "-timing.tsv"];
  CompletePath = StrCat[RawJobDirectory, "/", ColumnStem, ".done"];
If(PlotFieldMaps)
  FieldMapSuffix = Sprintf["_f%04.0f_b%04.0f.pos", FrequencyIndex, basis];
  FieldMapLabel = StrCat[Sprintf["; f=%.8g Hz; basis=", FrequencyHz],
    Str[TerminalNames(basis - 1)], Sprintf["; axial drive %.17g A", UnitSource],
    Sprintf["; Gamma=(%.17g,%.17g) 1/m", GammaRe, GammaIm],
    "; PML values are analytic continuation; phasor=real,imag"];
  PostOperation {
    { Name FEMWriteMaps~{basis}; NameOfPostProcessing FEMFields;
      LastTimeStepOnly 1;
      Operation {
        Print[bt_mesh, OnElementsOf DomainMedia_Ele,
          Name StrCat["Pullback of A_t/Gamma [T m2]", FieldMapLabel],
          File StrCat[MapDirectory, "/bt_mesh", FieldMapSuffix]];
        Print[v_local, OnElementsOf DomainMedia_Ele,
          Name StrCat["phi/Gamma [V m]; gauge dependent", FieldMapLabel],
          File StrCat[MapDirectory, "/v_local", FieldMapSuffix]];
        Print[hz_scaled, OnElementsOf DomainMedia_Ele,
          Name StrCat["Hz/Gamma [A]", FieldMapLabel],
          File StrCat[MapDirectory, "/hz_scaled", FieldMapSuffix]];
        Print[material_region, OnElementsOf DomainFields, Name "Material physical tag",
          File StrCat[MapDirectory, "/material_region", FieldMapSuffix]];
        Print[pml_mask, OnElementsOf DomainFields, Name "PML mask: 0 physical, 1 analytic continuation",
          File StrCat[MapDirectory, "/pml_mask", FieldMapSuffix]];
        Print[jt_mesh, OnElementsOf DomainMedia_Ele, Name StrCat["J_t/Gamma: mesh flux density [A/m]", FieldMapLabel],
          File StrCat[MapDirectory, "/jt_mesh", FieldMapSuffix]];
        Print[jt, OnElementsOf DomainMedia_Ele, Name StrCat["J_t/Gamma: physical components [A/m]", FieldMapLabel],
          File StrCat[MapDirectory, "/jt", FieldMapSuffix]];
        Print[az, OnElementsOf DomainFields, Name StrCat["Az [T m]", FieldMapLabel],
          File StrCat[MapDirectory, "/az", FieldMapSuffix]];
        Print[b, OnElementsOf DomainFields, Name StrCat["B [T]", FieldMapLabel],
          File StrCat[MapDirectory, "/b", FieldMapSuffix]];
        Print[bm, OnElementsOf DomainFields, Name StrCat["|B| [T]", FieldMapLabel],
          File StrCat[MapDirectory, "/bm", FieldMapSuffix]];
        Print[e, OnElementsOf DomainMedia_Ele, Name StrCat["E_t/Gamma [V]", FieldMapLabel],
          File StrCat[MapDirectory, "/e", FieldMapSuffix]];
        Print[ez, OnElementsOf DomainFields, Name StrCat["Ez [V/m]", FieldMapLabel],
          File StrCat[MapDirectory, "/ez", FieldMapSuffix]];
        Print[em, OnElementsOf DomainMedia_Ele, Name StrCat["|E_t/Gamma| [V]", FieldMapLabel],
          File StrCat[MapDirectory, "/em", FieldMapSuffix]];
        // Total current includes displacement in air and lossless dielectrics.
        Print[jz, OnElementsOf DomainFields, Name StrCat["Jz (conduction + displacement) [A/m2]", FieldMapLabel],
          File StrCat[MapDirectory, "/jz", FieldMapSuffix]];
        Print[jm, OnElementsOf DomainMedia_Ele, Name StrCat["|J_t/Gamma| [A/m]; axial drive 1 A", FieldMapLabel],
          File StrCat[MapDirectory, "/jm", FieldMapSuffix]];
        Print[rhoj2, OnElementsOf DomainFieldLoss, Name StrCat["S [W/m3]", FieldMapLabel],
          File StrCat[MapDirectory, "/rhoj2", FieldMapSuffix]];
      }
    }
  }
EndIf

PostOperation {
  { Name FEMAppendRaw~{basis}; NameOfPostProcessing FEMFields;
    LastTimeStepOnly 1;
    Format Table;
    NoMesh 1;
    Operation {
      DeleteFile[RawZPath];
      DeleteFile[RawPPath];
      DeleteFile[RawScalarPath];
      Print[voltage_gradient, OnElementsOf DomainMedia_Ele, Format Gmsh,
        StoreInField 100000+basis, File "", NoMesh 0];
      For response_terminal In {1:NumTerminals}
        // Scalar potential is a diagnostic only. It is absent from P.
        Print[ReP, OnRegion Terminal~{response_terminal}, Format Table,
          File "", StoreInVariable $FEMScalarRe];
        Print[ImP, OnRegion Terminal~{response_terminal}, Format Table,
          File "", StoreInVariable $FEMScalarIm];
        Print[{$FEMFrequencyIndex,$FEMFrequencyHz,response_terminal,
            $FEMBasisTerminal,$FEMScalarRe,$FEMScalarIm},
          Format "%g	%.17g	%g	%g	%.17g	%.17g",
          File RawScalarPath, AppendToExistingFile 1];
        Print[ReZ, OnRegion Terminal~{response_terminal}, Format Table,
          File "", StoreInVariable $FEMRawRe];
        Print[ImZ, OnRegion Terminal~{response_terminal}, Format Table,
          File "", StoreInVariable $FEMRawIm];
        Print[{$FEMLongitudinalDrive = Complex[$FEMRawRe, $FEMRawIm]},
          Format "%g", File ""];

        Print[ReVoltageLine~{basis}[MeasurementLine~{response_terminal}], OnGlobal,
          File "", StoreInVariable $FEMLineRe];
        Print[ImVoltageLine~{basis}[MeasurementLine~{response_terminal}], OnGlobal,
          File "", StoreInVariable $FEMLineIm];
        Print[{$FEMFrequencyIndex, $FEMFrequencyHz, response_terminal,
            $FEMBasisTerminal,
            $FEMLineRe/UnitSource,
            $FEMLineIm/UnitSource},
          Format "%g	%.17g	%g	%g	%.17g	%.17g",
          File RawPPath, AppendToExistingFile 1];
        Print[{$FEMP~{response_terminal}~{basis} = Complex[$FEMLineRe/UnitSource, $FEMLineIm/UnitSource]},
          Format "%g", File ""];
        Print[{$FEMZ~{response_terminal}~{basis} = $FEMLongitudinalDrive + gamma2[]*$FEMP~{response_terminal}~{basis}},
          Format "%g", File ""];
        Print[{$FEMFrequencyIndex, $FEMFrequencyHz, response_terminal,
            $FEMBasisTerminal, Re[$FEMZ~{response_terminal}~{basis}], Im[$FEMZ~{response_terminal}~{basis}]},
          Format "%g	%.17g	%g	%g	%.17g	%.17g",
          File RawZPath, AppendToExistingFile 1];
      EndFor
    }
  }

}

PostOperation {
  { Name FEMCompleteColumn~{basis}; NameOfPostProcessing FEMFields;
    LastTimeStepOnly 1;
    Format Table;
    Operation {
      Print[{$FEMFrequencyIndex, $FEMBasisTerminal,
          $FEMAssemblyStart - $FEMConstraintStart,
          $FEMSolveStart - $FEMAssemblyStart,
          $FEMOutputStart - $FEMSolveStart,
          $FEMOutputEnd - $FEMOutputStart, $FEMFirstSolve || !ReuseFactorization},
        Format "%g	%g	%.17g	%.17g	%.17g	%.17g	%g",
        File TimingPath];
      // Written last, after raw output and optional maps have closed.
      Print[{2, $FEMFrequencyIndex, $FEMFrequencyHz, $FEMBasisTerminal,
          NumTerminals, PlotFieldMaps},
        Format "%g	%g	%.17g	%g	%g	%g", File CompletePath];
    }
  }
}
EndFor
