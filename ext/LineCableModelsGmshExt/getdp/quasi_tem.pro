// Coupled quasi-TEM A_z / u_r / phi formulation.
// Each LineCableModelsFEMScan invocation solves exactly one frequency and
// basis-terminal job. Julia owns the outer scan and process isolation.

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
If(!Exists(PlotFieldMaps))
  PlotFieldMaps = 0;
EndIf
If(!Exists(RawOutputStem))
  RawOutputStem = "scan";
EndIf
RawDirectory = StrCat[RunDirectory, "/raw"];
RawJobDirectory = StrCat[RawDirectory, "/jobs"];
RawZPath = StrCat[RawJobDirectory, "/", RawOutputStem, "-Z.tsv"];
RawPPath = StrCat[RawJobDirectory, "/", RawOutputStem, "-P.tsv"];
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

Include MaterialFunctionsPath;
Include FieldMapDispatchPath;

Function {
  nu[#{Air, AirInf}] = 1. / mu0;
  sigma_dc[#{Air, AirInf}] = 0.;
  epsilon[#{Air, AirInf}] = eps0;
  mu[#{Air, AirInf}] = mu0;
  tan_delta[#{Air, AirInf}] = 0.;

  nu[#{Earth, EarthInf}] = 1. / mu_earth;
  sigma_dc[#{Earth, EarthInf}] = sigma_earth;
  epsilon[#{Earth, EarthInf}] = eps_earth;
  mu[#{Earth, EarthInf}] = mu_earth;
  tan_delta[#{Earth, EarthInf}] = 0.;

  gamma_prop[] = Complex[GammaQuasiTEMRe, GammaQuasiTEMIm];
  inv_gamma[] = 1. / gamma_prop[];
  omega[] = 2. * Pi * $FEMFrequencyHz;
  sigma[] = sigma_dc[] + omega[] * epsilon[] * tan_delta[];
  se[] = Complex[sigma[], omega[] * epsilon[]];
}

Constraint {
  { Name FEMMagneticVectorPotential;
    Case {
      { Region Sur_Dirichlet_Mag; Value 0.; }
    }
  }
  { Name FEMVoltageReference;
    Case {
      { Region Earth; Value 0.; }
    }
  }
  { Name FEMTerminalCurrent;
    Case {
      For t In {1:NumTerminals}
        { Region Terminal~{t}; Value $FEM_I~{t}; }
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
      { NameOfCoef U; EntityType Auto; NameOfConstraint FEMVoltageReference; }
      { NameOfCoef I; EntityType Auto; NameOfConstraint FEMTerminalCurrent; }
    }
  }

  { Name Hgrad_phi_FEM_2D; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef phin; Function BF_Node;
        Support Domain_Mag; Entity NodesOf[All, Not Terminals]; }
      { Name sf; NameOfCoef phif; Function BF_GroupOfNodes;
        Support Domain_Mag; Entity GroupsOfNodesOf[Terminals]; }
    }
    GlobalQuantity {
      { Name Phi; Type AliasOf; NameOfCoef phif; }
    }
    Constraint {
      { NameOfCoef phin; EntityType NodesOf;
        NameOfConstraint FEMScalarPotential; }
    }
  }
}

Formulation {
  { Name FEM_a_phi_2D; Type FemEquation;
    Quantity {
      { Name a; Type Local; NameOfSpace Hcurl_a_FEM_2D; }
      { Name ur; Type Local; NameOfSpace Hregion_u_FEM_2D; }
      { Name U; Type Global; NameOfSpace Hregion_u_FEM_2D [U]; }
      { Name I; Type Global; NameOfSpace Hregion_u_FEM_2D [I]; }
      { Name phi; Type Local; NameOfSpace Hgrad_phi_FEM_2D; }
      { Name Phi; Type Global; NameOfSpace Hgrad_phi_FEM_2D [Phi]; }
    }
    Equation {
      Galerkin {
        [nu[] * Dof{d a}, {d a}];
        In Domain_Mag; Jacobian Vol; Integration I1;
      }
      Galerkin {
        DtDof [sigma[] * Dof{a}, {a}];
        In DomainC; Jacobian Vol; Integration I1;
      }
      Galerkin {
        [sigma[] * Dof{ur}, {a}];
        In DomainC; Jacobian Vol; Integration I1;
      }
      Galerkin {
        [-gamma_prop[] * sigma[] * (Vector[0,0,1] * Dof{phi}), {a}];
        In DomainC; Jacobian Vol; Integration I1;
      }
      Galerkin {
        DtDof [sigma[] * Dof{a}, {ur}];
        In DomainC; Jacobian Vol; Integration I1;
      }
      Galerkin {
        [sigma[] * Dof{ur}, {ur}];
        In DomainC; Jacobian Vol; Integration I1;
      }
      Galerkin {
        [-gamma_prop[] * sigma[] * (Vector[0,0,1] * Dof{phi}), {ur}];
        In DomainC; Jacobian Vol; Integration I1;
      }
      Galerkin {
        DtDof [gamma_prop[] * sigma[] * (Dof{a} * Vector[0,0,1]), {phi}];
        In DomainC; Jacobian Vol; Integration I1;
      }
      Galerkin {
        [gamma_prop[] * sigma[] * (Dof{ur} * Vector[0,0,1]), {phi}];
        In DomainC; Jacobian Vol; Integration I1;
      }
      Galerkin {
        [-gamma_prop[] * gamma_prop[] * sigma[] * Dof{phi}, {phi}];
        In DomainC; Jacobian Vol; Integration I1;
      }
      Galerkin {
        [sigma[] * Dof{d phi}, {d phi}];
        In DomainC; Jacobian Vol; Integration I1;
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
        DtDof [-gamma_prop[] * epsilon[] *
          (Vector[0,0,1] * Dof{phi}), {a}];
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
      Galerkin {
        DtDof [-gamma_prop[] * epsilon[] *
          (Vector[0,0,1] * Dof{phi}), {ur}];
        In Domain_Mag; Jacobian Vol; Integration I1;
      }
      Galerkin {
        DtDtDof [gamma_prop[] * epsilon[] *
          (Dof{a} * Vector[0,0,1]), {phi}];
        In Domain_Mag; Jacobian Vol; Integration I1;
      }
      Galerkin {
        DtDof [gamma_prop[] * epsilon[] *
          (Dof{ur} * Vector[0,0,1]), {phi}];
        In Domain_Mag; Jacobian Vol; Integration I1;
      }
      Galerkin {
        DtDof [-gamma_prop[] * gamma_prop[] * epsilon[] * Dof{phi}, {phi}];
        In Domain_Mag; Jacobian Vol; Integration I1;
      }
      Galerkin {
        DtDof [epsilon[] * Dof{d phi}, {d phi}];
        In Domain_Mag; Jacobian Vol; Integration I1;
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
EndFor
Return

Macro FEMSolveBasis
Call FEMSetBasisCurrent;
UpdateConstraint[Sys_FEM];
Generate[Sys_FEM];
Solve[Sys_FEM];
Return

Resolution {
  { Name LineCableModelsFEMScan;
    System {
      { Name Sys_FEM; NameOfFormulation FEM_a_phi_2D;
        Type Complex; Frequency 1.; }
    }
    Operation {
      CreateDir[RawDirectory];
      CreateDir[MapDirectory];
      InitSolution[Sys_FEM];
      Evaluate[$FEMFrequencyIndex = FrequencyIndex];
      Evaluate[$FEMFrequencyHz = FrequencyHz];
      SetFrequency[Sys_FEM, FrequencyHz];
      SetTimeStep[1];
      // GetDP time is an output-step identity, not the physical frequency.
      // The scan index avoids time-range failures; SetFrequency owns physics.
      SetTime[FrequencyIndex];
      Evaluate[$FEMBasisTerminal = BasisTerminal];
      Call FEMSolveBasis;
      PostOperation[FEMAppendRaw];
      If(PlotFieldMaps)
        Call FEMWriteMaps;
      EndIf
    }
  }
}

PostProcessing {
  { Name FEMFields; NameOfFormulation FEM_a_phi_2D; NameOfSystem Sys_FEM;
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
      { Name e; Value {
        Term { [-(Dt[{a}] + {ur} - gamma_prop[] *
          (Vector[0,0,1] * {phi}) + {d phi})];
          In Domain_Mag; Jacobian Vol; }
      }}
      { Name ez; Value {
        Term { [-CompZ[Dt[{a}] + {ur} - gamma_prop[] *
          (Vector[0,0,1] * {phi})]];
          In Domain_Mag; Jacobian Vol; }
      }}
      { Name em; Value {
        Term { [Norm[-(Dt[{a}] + {ur} - gamma_prop[] *
          (Vector[0,0,1] * {phi}) + {d phi})]];
          In Domain_Mag; Jacobian Vol; }
      }}
      { Name jz; Value {
        Term { [CompZ[-se[] * (Dt[{a}] + {ur} - gamma_prop[] *
          (Vector[0,0,1] * {phi}) + {d phi})]];
          In Domain_Mag; Jacobian Vol; }
      }}
      { Name jm; Value {
        Term { [Norm[-se[] * (Dt[{a}] + {ur} - gamma_prop[] *
          (Vector[0,0,1] * {phi}) + {d phi})]];
          In Domain_Mag; Jacobian Vol; }
      }}
      { Name rhoj2; Value {
        Term { [0.5 * sigma[] * SquNorm[Dt[{a}] + {ur} - gamma_prop[] *
          (Vector[0,0,1] * {phi}) + {d phi}]];
          In DomainC; Jacobian Vol; }
      }}
      { Name ReZ; Value {
        Term { [-Re[{U} / UnitSource]]; In DomainCWithI; }
      }}
      { Name ImZ; Value {
        Term { [-Im[{U} / UnitSource]]; In DomainCWithI; }
      }}
      { Name ReP; Value {
        Term { [Re[{Phi} * inv_gamma[] / UnitSource]]; In Terminals; }
      }}
      { Name ImP; Value {
        Term { [Im[{Phi} * inv_gamma[] / UnitSource]]; In Terminals; }
      }}
    }
  }
}

Include FieldMapOperationsPath;

PostOperation {
  { Name FEMAppendRaw; NameOfPostProcessing FEMFields;
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
