// Included by model.pro when the ONELAB Physics selector is 1 (quasi-fw).
Jacobian { { Name Plain; Case { { Region All; Jacobian Vol; } } } }



// First-order Maxwell formulation for exp(j omega t - Gamma z).
// A_z=a, A_t=Gamma*bt, phi=Gamma*v. The equations are the normalized
// Gamma -> 0 limit; all O(Gamma) transverse fields are retained. This is
// not a finite-Gamma modal solver. The finite-metal A_z/u branch is retained;
// PerfectConductors=1 instead uses exact PEC contour current constraints.
//
// One imposed axial current supplies both the magnetic equation and the
// normalized transverse leakage. There is no independent electric drive.
// The transverse domain treats metal surfaces as equipotential electrodes;
// the magnetic block either retains the specified finite metal conductivity
// or excludes metal interiors when PerfectConductors=1.
//
// Voltage/Gamma = v_i-v_ref + j omega integral(bt.dl). PathDataPath supplies
// mesh-coordinate integration points and oriented line weights from the
// reference at earth infinity to each electrode. Jacobian Plain evaluates
// the pulled-back 1-form: its circulation equals the physical circulation,
// including in the infinite-element shell. The backend generates these paths
// from the mesh; dev/run_quasi_full.jl also exposes a manual PEC CLI experiment.
// Raw P has units ohm m (inverse admittance); analytical Pe = j omega P.
//
// Reference: G. Ciuprina and R. V. Sabriego, "Electric circuit element boundary
// conditions for electromagneto-quasistatic and full wave models in A, phi
// potentials and their finite element implementation", J. Math. Ind. 14, 27
// (2024), https://doi.org/10.1186/s13362-024-00165-6, Sec. 4 and Appendix B.
// The paper gives the full-vector A/phi equations, terminal conditions and
// gauge framework. This 2D longitudinal expansion and the vertical voltage
// paths are our reduction, not the paper's 3D ECE model or its Darwin model.
// See the LineCableModelsFEM Julia docstring for equations, units and limits.

// PerfectConductors=1 excludes metal interiors and imposes exact PEC
// current-carrying contours. Default 0 retains the working finite-metal a/u block.
If(!Exists(PerfectConductors))
  PerfectConductors = 0;
EndIf
If(!Exists(PathDataPath))
  Error("Pass voltage paths with -setstring PathDataPath");
EndIf
Include PathDataPath;
If(#PathStart() != NumTerminals || #PathCount() != NumTerminals)
  Error("Provide one oriented voltage path per terminal");
EndIf
If(#PathX() != #PathY() || #PathX() != #PathDX() || #PathX() != #PathDY())
  Error("Voltage-path coordinates and line weights must have equal lengths");
EndIf
For terminal In {0:NumTerminals-1}
  If(PathCount(terminal) < 1 || PathStart(terminal) < 0 ||
     PathStart(terminal)+PathCount(terminal) > #PathX())
    Error("Invalid voltage-path point range");
  EndIf
EndFor
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
If(!Exists(BasisTerminal))
  BasisTerminal = 1;
EndIf
If(Exists(BasisListPath))
  Include BasisListPath;
Else
  RequestedBases() = {BasisTerminal};
EndIf
If(!Exists(ReuseFactorization))
  ReuseFactorization = 1;
EndIf
If(!Exists(PlotFieldMaps))
  PlotFieldMaps = 0;
EndIf
If(!Exists(RawOutputStem))
  RawOutputStem = "scan";
EndIf
RawDirectory = StrCat[RunDirectory, "/raw"];
RawJobDirectory = StrCat[RawDirectory, "/jobs"];
MapDirectory = StrCat[RunDirectory, "/maps"];

Group {
  Air = Region[{AIR_EM}];
  Earth = Region[{EARTH_EM}];
  AirInf = Region[{AIR_INF}];
  EarthInf = Region[{EARTH_INF}];
  Terminals = Region[{}];
  TerminalContours = Region[{}];

  For t In {1:NumTerminals}
    Terminal~{t} = Region[{(TERMINAL + t - 1)}];
    TerminalContour~{t} = Region[{(TERMINAL_CONTOUR + t - 1)}];
    Terminals += Region[{Terminal~{t}}];
    TerminalContours += Region[{TerminalContour~{t}}];
  EndFor
  Sur_Dirichlet_Mag = Region[{OUTBND_EM}];
  Sur_Insulation_Ele = Region[{OUTBND_ELE_INS}];
  Sur_Dirichlet_Ele = Region[{OUTBND_ELE_REF}];
  // Complete the gauge tree on EVERY boundary with prescribed bt circulation.
  // Omitting the electrode contours would constrain physical loop degrees of
  // freedom in addition to the potential gauge.
  GaugeBoundary = Region[{Sur_Dirichlet_Mag, TerminalContours}];
}

Include "materials.pro";
Group {
  // Conductors are equipotential electrodes, not electric field media.
  DomainMedia_Ele = Region[{Air, AirInf, Earth, EarthInf, PassiveMaterialRegions}];
  If(PerfectConductors)
    DomainFields = Region[{DomainMedia_Ele}];
    DomainFieldLoss = Region[{Earth, EarthInf, LossyMaterialRegions}];
  Else
    DomainFields = Region[{Domain_Mag}];
    DomainFieldLoss = Region[{DomainLoss}];
  EndIf
}

Function {
  nu[#{Air, AirInf}] = 1. / AirMu;
  sigma[#{Air, AirInf}] = 0.;
  epsilon[#{Air, AirInf}] = AirEpsilon;
  mu[#{Air, AirInf}] = AirMu;

  nu[#{Earth, EarthInf}] = 1. / EarthMu(FrequencyIndex - 1);
  sigma[#{Earth, EarthInf}] = EarthSigma(FrequencyIndex - 1);
  epsilon[#{Earth, EarthInf}] = EarthEpsilon(FrequencyIndex - 1);
  mu[#{Earth, EarthInf}] = EarthMu(FrequencyIndex - 1);

  omega[] = 2. * Pi * $FEMFrequencyHz;
  se[] = Complex[sigma[], omega[] * epsilon[]];
}

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
      // The earth-side far boundary is the electric reference. The air-side
      // outer shell is intentionally unconstrained: the Galerkin formulation
      // supplies its natural zero-normal-current boundary condition.
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
  { Name FEM_QuasiFull_2D; Type FemEquation;
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
        Galerkin { [nu[] * Dof{d a}, {d a}];
          In DomainMedia_Ele; Jacobian Vol; Integration I1; }
        Galerkin { [Complex[0,omega[]] * se[] * Dof{a}, {a}];
          In DomainMedia_Ele; Jacobian Vol; Integration I1; }
        GlobalTerm { [-Dof{I}, {A}]; In Terminals; }
      Else
      Galerkin {
        [nu[] * Dof{d a}, {d a}];
        In Domain_Mag; Jacobian Vol; Integration I1;
      }
      Galerkin {
        DtDof [sigma[] * Dof{a}, {a}];
        In DomainLoss; Jacobian Vol; Integration I1;
      }
      Galerkin {
        [sigma[] * Dof{ur}, {a}];
        In DomainLoss; Jacobian Vol; Integration I1;
      }

      Galerkin {
        DtDof [sigma[] * Dof{a}, {ur}];
        In DomainLoss; Jacobian Vol; Integration I1;
      }
      Galerkin {
        [sigma[] * Dof{ur}, {ur}];
        In DomainLoss; Jacobian Vol; Integration I1;
      }

      Galerkin {
        DtDtDof [epsilon[] * Dof{a}, {a}];
        In Domain_Mag; Jacobian Vol; Integration I1;
      }
      Galerkin {
        DtDof [epsilon[] * Dof{ur}, {a}];
        In Domain_Mag; Jacobian Vol; Integration I1;
      }

      Galerkin {
        DtDtDof [epsilon[] * Dof{a}, {ur}];
        In Domain_Mag; Jacobian Vol; Integration I1;
      }
      Galerkin {
        DtDof [epsilon[] * Dof{ur}, {ur}];
        In Domain_Mag; Jacobian Vol; Integration I1;
      }

        GlobalTerm { [Dof{I}, {U}]; In DomainCWithI; }
      EndIf
      // First-order transverse Ampere, A_t = Gamma * bt, phi = Gamma * v:
      // C*(nu C bt) + j omega se bt + se grad(v) - nu grad(a) = 0.
      // The tree gauge removes gradient degrees of freedom from bt; the
      // continuity rows below provide that part of Ampere in the nodal space.
      Galerkin { [nu[] * Dof{d bt}, {d bt}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0,omega[]] * se[] * Dof{bt}, {bt}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      Galerkin { [se[] * Dof{d v}, {bt}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      // d a = (d_y a, -d_x a, 0); rotating restores grad_t(a).
      Galerkin { [-nu[] * (Vector[0,0,1] /\ Dof{d a}), {bt}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      // Normalized continuity; the source is the axial current, not a second drive.
      Galerkin { [se[] * Dof{d v}, {d v}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0,omega[]] * se[] * Dof{bt}, {d v}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0,omega[]] * se[] * (Dof{a} * Vector[0,0,1]), {v}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1; }
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
Evaluate[$FEMOutputStart = GetWallClockTime[]];
Return

Resolution {
  { Name LineCableModelsFEMScan;
    System {
      { Name Sys_FEM; NameOfFormulation FEM_QuasiFull_2D;
        Type Complex; Frequency 1.; }
    }
    Operation {
      CreateDir[RawDirectory];
      CreateDir[RawJobDirectory];
      CreateDir[MapDirectory];
      InitSolution[Sys_FEM];
      Evaluate[$FEMFrequencyIndex = FrequencyIndex];
      Evaluate[$FEMFrequencyHz = FrequencyHz];
      SetFrequency[Sys_FEM, FrequencyHz];
      SetTimeStep[1];
      // GetDP time is an output-step identity, not the physical frequency.
      // The scan index avoids time-range failures; SetFrequency owns physics.
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
    }
  }
}

PostProcessing {
  { Name FEMFields; NameOfFormulation FEM_QuasiFull_2D; NameOfSystem Sys_FEM;
    PostQuantity {
      { Name ReBX; Value { Term { [Re[CompX[{bt}]]]; In DomainMedia_Ele; Jacobian Plain; } } }
      { Name ImBX; Value { Term { [Im[CompX[{bt}]]]; In DomainMedia_Ele; Jacobian Plain; } } }
      { Name ReBY; Value { Term { [Re[CompY[{bt}]]]; In DomainMedia_Ele; Jacobian Plain; } } }
      { Name ImBY; Value { Term { [Im[CompY[{bt}]]]; In DomainMedia_Ele; Jacobian Plain; } } }
      { Name bt_mesh; Value { Term { [{bt}]; In DomainMedia_Ele; Jacobian Plain; } } }
      { Name v_local; Value { Term { [{v}]; In DomainMedia_Ele; Jacobian Plain; } } }
      { Name hz_scaled; Value { Term { [nu[]*CompZ[{d bt}]]; In DomainMedia_Ele; Jacobian Vol; } } }

      { Name az; Value {
        Term { [CompZ[{a}]]; In DomainFields; Jacobian Vol; }
      }}
      { Name b; Value {
        Term { [{d a}]; In DomainFields; Jacobian Vol; }
      }}
      { Name bm; Value {
        Term { [Norm[{d a}]]; In DomainFields; Jacobian Vol; }
      }}
      // e is E_t/Gamma [V]; ez is the driven axial field [V/m].
      { Name e; Value {
        Term { [-{d v}-Complex[0,omega[]]*{bt}]; In DomainMedia_Ele; Jacobian Vol; }
      }}
      If(PerfectConductors)
        { Name ez; Value { Term { [-CompZ[Dt[{a}]]]; In DomainFields; Jacobian Vol; } } }
      Else
      { Name ez; Value {
        Term { [-CompZ[Dt[{a}] + {ur}]]; In DomainFields; Jacobian Vol; }
      }}
      EndIf
      { Name em; Value {
        Term { [Norm[-{d v}-Complex[0,omega[]]*{bt}]]; In DomainMedia_Ele; Jacobian Vol; }
      }}
      If(PerfectConductors)
        { Name jz; Value { Term { [-CompZ[se[] * Dt[{a}]]]; In DomainFields; Jacobian Vol; } } }
      Else
      { Name jz; Value {
        Term { [-CompZ[se[] * (Dt[{a}] + {ur})]];
          In DomainFields; Jacobian Vol; }
      }}
      EndIf
      { Name jm; Value {
        Term { [Norm[-se[] * ({d v}+Complex[0,omega[]]*{bt})]]; In DomainMedia_Ele; Jacobian Vol; }
      }}
      If(PerfectConductors)
        { Name rhoj2; Value { Term { [0.5 * sigma[] * SquNorm[Dt[{a}]]]; In DomainFieldLoss; Jacobian Vol; } } }
      Else
      { Name rhoj2; Value {
        Term { [0.5 * sigma[] * SquNorm[Dt[{a}] + {ur}]];
          In DomainFieldLoss; Jacobian Vol; }
      }}
      EndIf
      If(PerfectConductors)
        { Name ReZ; Value { Term { [Re[Complex[0,omega[]]*{A}/UnitSource]]; In Terminals; } } }
        { Name ImZ; Value { Term { [Im[Complex[0,omega[]]*{A}/UnitSource]]; In Terminals; } } }
      Else
      { Name ReZ; Value {
        Term { [-Re[{U} / UnitSource]]; In DomainCWithI; }
      }}
      { Name ImZ; Value {
        Term { [-Im[{U} / UnitSource]]; In DomainCWithI; }
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
  FieldMapLabel = StrCat[Sprintf["; f=%.17g Hz; basis=", FrequencyHz],
    Str[TerminalNames(basis - 1)]];
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
        Print[az, OnElementsOf DomainFields, Name StrCat["Az [T m]", FieldMapLabel],
          File StrCat[MapDirectory, "/az", FieldMapSuffix]];
        Print[b, OnElementsOf DomainFields, Name StrCat["B [T]", FieldMapLabel],
          File StrCat[MapDirectory, "/b", FieldMapSuffix]];
        Print[bm, OnElementsOf DomainFields, Name StrCat["|B| [T]", FieldMapLabel],
          File StrCat[MapDirectory, "/bm", FieldMapSuffix]];
        Print[e, OnElementsOf DomainMedia_Ele, Name StrCat["E_t/Gamma [V]; axial drive 1 A", FieldMapLabel],
          File StrCat[MapDirectory, "/e", FieldMapSuffix]];
        Print[ez, OnElementsOf DomainFields, Name StrCat["Ez [V/m]; axial drive 1 A", FieldMapLabel],
          File StrCat[MapDirectory, "/ez", FieldMapSuffix]];
        Print[em, OnElementsOf DomainMedia_Ele, Name StrCat["|E_t/Gamma| [V]; axial drive 1 A", FieldMapLabel],
          File StrCat[MapDirectory, "/em", FieldMapSuffix]];
        Print[jz, OnElementsOf DomainFieldLoss, Name StrCat["Jz [A/m2]", FieldMapLabel],
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
      For response_terminal In {1:NumTerminals}
        Print[ReZ, OnRegion Terminal~{response_terminal}, Format Table,
          File "", StoreInVariable $FEMRawRe];
        Print[ImZ, OnRegion Terminal~{response_terminal}, Format Table,
          File "", StoreInVariable $FEMRawIm];
        Print[{$FEMFrequencyIndex, $FEMFrequencyHz, response_terminal,
            $FEMBasisTerminal, $FEMRawRe, $FEMRawIm},
          Format "%g	%.17g	%g	%g	%.17g	%.17g",
          File RawZPath, AppendToExistingFile 1];

        // Store the scalar trace separately, to make the vector contribution
        // inspectable. It is gauge-dependent and is NOT an inverse admittance.
        Print[ReP, OnRegion Terminal~{response_terminal}, Format Table,
          File "", StoreInVariable $FEMScalarRe];
        Print[ImP, OnRegion Terminal~{response_terminal}, Format Table,
          File "", StoreInVariable $FEMScalarIm];
        Print[{$FEMFrequencyIndex, $FEMFrequencyHz, response_terminal,
            $FEMBasisTerminal, $FEMScalarRe, $FEMScalarIm},
          Format "%g	%.17g	%g	%g	%.17g	%.17g",
          File RawScalarPath, AppendToExistingFile 1];
        Print[{$FEMLineRe = 0., $FEMLineIm = 0.}, Format "%g %g", File ""];
        For point In {PathStart(response_terminal-1):PathStart(response_terminal-1)+PathCount(response_terminal-1)-1}
          Print[ReBX, OnPoint {PathX(point),PathY(point),0}, File "", StoreInVariable $FEMBxRe];
          Print[ImBX, OnPoint {PathX(point),PathY(point),0}, File "", StoreInVariable $FEMBxIm];
          Print[ReBY, OnPoint {PathX(point),PathY(point),0}, File "", StoreInVariable $FEMByRe];
          Print[ImBY, OnPoint {PathX(point),PathY(point),0}, File "", StoreInVariable $FEMByIm];
          Print[{$FEMLineRe = $FEMLineRe + PathDX(point)*$FEMBxRe + PathDY(point)*$FEMByRe,
                 $FEMLineIm = $FEMLineIm + PathDX(point)*$FEMBxIm + PathDY(point)*$FEMByIm},
            Format "%.17g %.17g", File ""];
        EndFor
        Print[{$FEMFrequencyIndex, $FEMFrequencyHz, response_terminal,
            $FEMBasisTerminal,
            $FEMScalarRe - 2*Pi*$FEMFrequencyHz*$FEMLineIm/UnitSource,
            $FEMScalarIm + 2*Pi*$FEMFrequencyHz*$FEMLineRe/UnitSource},
          Format "%g	%.17g	%g	%g	%.17g	%.17g",
          File RawPPath, AppendToExistingFile 1];
      EndFor
    }
  }

}

PostOperation {
  { Name FEMCompleteColumn~{basis}; NameOfPostProcessing FEMFields;
    LastTimeStepOnly 1;
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
