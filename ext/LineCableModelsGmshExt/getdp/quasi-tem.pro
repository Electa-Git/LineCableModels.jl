// Quasi-TEM series impedance and scalar electrodynamic potential coefficients.
// The independent A_z/u_r and v blocks share one assembled system and one
// factorization per frequency. Z uses unit axial current [A]; P uses unit
// outward transverse terminal current [A/m]. No small Gamma is required.
// Each invocation owns one mesh/frequency and reuses its operator across
// terminal excitations. Julia owns frequency scheduling and checkpoints.

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
}

Include "materials.pro";
Group {
  // Conductors are equipotential electrodes, not electric field media.
  DomainMedia_Ele = Region[{Air, AirInf, Earth, EarthInf, PassiveMaterialRegions}];
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
  { Name FEMTransverseCurrent;
    Case {
      For t In {1:NumTerminals}
        { Region Terminal~{t}; Value $FEM_Q~{t}; }
      EndFor
    }
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

  { Name Hgrad_v_FEM_2D; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Domain_Mag; Entity NodesOf[All, Not Terminals]; }
      { Name sf; NameOfCoef vf; Function BF_GroupOfNodes;
        Support Domain_Mag; Entity GroupsOfNodesOf[Terminals]; }
    }
    GlobalQuantity {
      { Name V; Type AliasOf; NameOfCoef vf; }
      { Name Q; Type AssociatedWith; NameOfCoef vf; }
    }
    Constraint {
      { NameOfCoef Q; EntityType Auto; NameOfConstraint FEMTransverseCurrent; }
      { NameOfCoef vn; EntityType NodesOf;
        NameOfConstraint FEMScalarPotential; }
    }
  }
}

Formulation {
  { Name FEM_Z_P_2D; Type FemEquation;
    Quantity {
      { Name a; Type Local; NameOfSpace Hcurl_a_FEM_2D; }
      { Name ur; Type Local; NameOfSpace Hregion_u_FEM_2D; }
      { Name U; Type Global; NameOfSpace Hregion_u_FEM_2D [U]; }
      { Name I; Type Global; NameOfSpace Hregion_u_FEM_2D [I]; }
      { Name v; Type Local; NameOfSpace Hgrad_v_FEM_2D; }
      { Name V; Type Global; NameOfSpace Hgrad_v_FEM_2D [V]; }
      { Name Q; Type Global; NameOfSpace Hgrad_v_FEM_2D [Q]; }
    }
    Equation {
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

      // div(se grad(v)) + se k^2 v = 0, k^2 = -j omega mu se.
      // The electric block retains diffusion/displacement at Gamma = 0.
      Galerkin {
        [se[] * Dof{d v}, {d v}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1;
      }
      Galerkin {
        [Complex[0, omega[]] * mu[] * se[]^2 * Dof{v}, {v}];
        In DomainMedia_Ele; Jacobian Vol; Integration I1;
      }
      GlobalTerm {
        [Dof{Q}, {V}]; In Terminals;
      }
      GlobalTerm {
        [Dof{I}, {U}]; In DomainCWithI;
      }
    }
  }
}

Macro FEMSetBasisCurrent
For t In {1:NumTerminals}
  Evaluate[$FEM_I~{t} = Complex[
    UnitSource * ($FEMBasisTerminal == t), 0.]];
  Evaluate[$FEM_Q~{t} = Complex[
    -UnitTransverseSource * ($FEMBasisTerminal == t), 0.]];
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
      { Name Sys_FEM; NameOfFormulation FEM_Z_P_2D;
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
  { Name FEMFields; NameOfFormulation FEM_Z_P_2D; NameOfSystem Sys_FEM;
    PostQuantity {
      { Name az; Value {
        Term { [CompZ[{a}]]; In Domain_Mag; Jacobian Vol; }
      }}
      { Name b; Value {
        Term { [{d a}]; In Domain_Mag; Jacobian Vol; }
      }}
      { Name bm; Value {
        Term { [Norm[{d a}]]; In Domain_Mag; Jacobian Vol; }
      }}
      // These are two independent terminal excitations. Do not combine their
      // axial and transverse electric fields into a fictitious full-wave field.
      { Name e; Value {
        Term { [-{d v}]; In DomainMedia_Ele; Jacobian Vol; }
      }}
      { Name ez; Value {
        Term { [-CompZ[Dt[{a}] + {ur}]]; In Domain_Mag; Jacobian Vol; }
      }}
      { Name em; Value {
        Term { [Norm[-{d v}]]; In DomainMedia_Ele; Jacobian Vol; }
      }}
      { Name jz; Value {
        Term { [-CompZ[se[] * (Dt[{a}] + {ur})]];
          In Domain_Mag; Jacobian Vol; }
      }}
      { Name jm; Value {
        Term { [Norm[-se[] * {d v}]]; In DomainMedia_Ele; Jacobian Vol; }
      }}
      { Name rhoj2; Value {
        Term { [0.5 * sigma[] * SquNorm[Dt[{a}] + {ur}]];
          In DomainLoss; Jacobian Vol; }
      }}
      { Name ReZ; Value {
        Term { [-Re[{U} / UnitSource]]; In DomainCWithI; }
      }}
      { Name ImZ; Value {
        Term { [-Im[{U} / UnitSource]]; In DomainCWithI; }
      }}
      { Name ReP; Value {
        Term { [Re[{V} / UnitTransverseSource]]; In Terminals; }
      }}
      { Name ImP; Value {
        Term { [Im[{V} / UnitTransverseSource]]; In Terminals; }
      }}
    }
  }
}

// Expand output declarations only for the requested terminal columns.
For basis_slot In {0:#RequestedBases()-1}
  basis = RequestedBases(basis_slot);
  ColumnStem = Sprintf["getdp-f%04.0f-b%04.0f", FrequencyIndex, basis];
  RawZPath = StrCat[RawJobDirectory, "/", ColumnStem, "-Z.tsv"];
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
        Print[az, OnElementsOf Domain_Mag, Name StrCat["Az [T m]", FieldMapLabel],
          File StrCat[MapDirectory, "/az", FieldMapSuffix]];
        Print[b, OnElementsOf Domain_Mag, Name StrCat["B [T]", FieldMapLabel],
          File StrCat[MapDirectory, "/b", FieldMapSuffix]];
        Print[bm, OnElementsOf Domain_Mag, Name StrCat["|B| [T]", FieldMapLabel],
          File StrCat[MapDirectory, "/bm", FieldMapSuffix]];
        Print[e, OnElementsOf DomainMedia_Ele, Name StrCat["E_t [V/m]; transverse drive 1 A/m", FieldMapLabel],
          File StrCat[MapDirectory, "/e", FieldMapSuffix]];
        Print[ez, OnElementsOf Domain_Mag, Name StrCat["Ez [V/m]; axial drive 1 A", FieldMapLabel],
          File StrCat[MapDirectory, "/ez", FieldMapSuffix]];
        Print[em, OnElementsOf DomainMedia_Ele, Name StrCat["|E_t| [V/m]; transverse drive 1 A/m", FieldMapLabel],
          File StrCat[MapDirectory, "/em", FieldMapSuffix]];
        Print[jz, OnElementsOf DomainLoss, Name StrCat["Jz [A/m2]", FieldMapLabel],
          File StrCat[MapDirectory, "/jz", FieldMapSuffix]];
        Print[jm, OnElementsOf DomainMedia_Ele, Name StrCat["|J_t| [A/m2]; transverse drive 1 A/m", FieldMapLabel],
          File StrCat[MapDirectory, "/jm", FieldMapSuffix]];
        Print[rhoj2, OnElementsOf DomainLoss, Name StrCat["S [W/m3]", FieldMapLabel],
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

        Print[ReP, OnRegion Terminal~{response_terminal}, Format Table,
          File "", StoreInVariable $FEMRawRe];
        Print[ImP, OnRegion Terminal~{response_terminal}, Format Table,
          File "", StoreInVariable $FEMRawIm];
        Print[{$FEMFrequencyIndex, $FEMFrequencyHz, response_terminal,
            $FEMBasisTerminal, $FEMRawRe, $FEMRawIm},
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
