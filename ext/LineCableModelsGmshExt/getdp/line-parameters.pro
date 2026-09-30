// Native line-parameter algebra. Primitive entries are $FEMZ_i_j and $FEMP_i_j.
// Matrix index 0 denotes Z [ohm/m], index 1 denotes P [ohm m].
// The connection permutation, bundle basis and Schur complements match the
// line-parameter engine. All coefficients remain runtime complex values.

FEMPhases() = {}; FEMFirsts() = {}; FEMGround() = {};
For t In {1:NumTerminals}
  phase = Connections(t-1);
  If(phase < 0 || Floor[phase] != phase)
    Error("Connections must contain nonnegative integer phase IDs");
  EndIf
  If(phase == 0)
    FEMGround() += {t};
  Else
    found = 0;
    If(#FEMPhases())
      For k In {0:#FEMPhases()-1}
        If(FEMPhases(k) == phase) found = 1; EndIf
      EndFor
    EndIf
    If(!found) FEMPhases() += {phase}; FEMFirsts() += {t}; EndIf
  EndIf
EndFor
FEMPermutation() = FEMFirsts();
If(#FEMPhases())
  For k In {0:#FEMPhases()-1}
    For t In {1:NumTerminals}
      If(Connections(t-1) == FEMPhases(k) && t != FEMFirsts(k))
        FEMPermutation() += {t};
      EndIf
    EndFor
  EndFor
EndIf
FEMPermutation() += FEMGround();
FEMKeep() = {}; FEMEliminate() = {}; FEMPhaseMap() = {};
For i In {0:NumTerminals-1}
  phase = Connections(FEMPermutation(i)-1); FEMBundleFirst~{i} = -1;
  If(i > 0 && phase > 0)
    For j In {0:i-1}
      If(FEMBundleFirst~{i} < 0 && Connections(FEMPermutation(j)-1) == phase)
        FEMBundleFirst~{i} = j;
      EndIf
    EndFor
  EndIf
  If(ReduceBundle && FEMBundleFirst~{i} >= 0) phase = 0;
  ElseIf(ReduceBundle && !KronReduction && phase == 0) phase = -1;
  EndIf
  If((ReduceBundle || KronReduction) && phase == 0)
    FEMEliminate() += {i+1};
  Else
    FEMKeep() += {i+1}; FEMPhaseMap() += {phase};
  EndIf
EndFor
FEMRetained = #FEMKeep(); FEMRemoved = #FEMEliminate();
If(!FEMRetained) Error("The connection/reduction settings retain no terminal"); EndIf

// A carrier group only allocates global algebraic unknowns; it adds no PDE.
Group { FEMAlgebraCarrier = Region[TERMINAL]; }
FunctionSpace {
  For i In {1:FEMRetained}
    { Name HAdmittance~{i}; Type Form0;
      BasisFunction { { Name s; NameOfCoef a; Function BF_GroupOfNodes;
        Support FEMAlgebraCarrier; Entity GroupsOfNodesOf[FEMAlgebraCarrier]; } }
      GlobalQuantity { { Name value; Type AliasOf; NameOfCoef a; } }
    }
  EndFor
  If(FEMRemoved)
    For matrix In {0:1}
      For i In {1:FEMRemoved}
        { Name HSchur~{matrix}~{i}; Type Form0;
          BasisFunction { { Name s; NameOfCoef a; Function BF_GroupOfNodes;
            Support FEMAlgebraCarrier; Entity GroupsOfNodesOf[FEMAlgebraCarrier]; } }
          GlobalQuantity { { Name value; Type AliasOf; NameOfCoef a; } }
        }
      EndFor
    EndFor
  EndIf
}
Formulation {
  { Name FEMAdmittance; Type FemEquation;
    Quantity {
      For i In {1:FEMRetained}
        { Name y~{i}; Type Global; NameOfSpace HAdmittance~{i}[value]; }
      EndFor
    }
    Equation {
      For i In {1:FEMRetained}
        For j In {1:FEMRetained}
          GlobalTerm { [$FEMReduced~{1}~{i}~{j} * Dof{y~{j}}, {y~{i}}]; In FEMAlgebraCarrier; }
        EndFor
        GlobalTerm { [-($FEMAlgebraColumn == i), {y~{i}}]; In FEMAlgebraCarrier; }
      EndFor
    }
  }
  If(FEMRemoved)
    For matrix In {0:1}
      { Name FEMSchur~{matrix}; Type FemEquation;
        Quantity {
          For i In {1:FEMRemoved}
            { Name x~{i}; Type Global; NameOfSpace HSchur~{matrix}~{i}[value]; }
          EndFor
        }
        Equation {
          For i In {1:FEMRemoved}
            row = FEMEliminate(i-1);
            For j In {1:FEMRemoved}
              col = FEMEliminate(j-1);
              GlobalTerm { [$FEMWork~{matrix}~{row}~{col} * Dof{x~{j}}, {x~{i}}]; In FEMAlgebraCarrier; }
            EndFor
            GlobalTerm { [-GetVariable[matrix,row,$FEMAlgebraColumn]{$FEMWork}, {x~{i}}]; In FEMAlgebraCarrier; }
          EndFor
        }
      }
    EndFor
  EndIf
}

Macro FEMMatrixSystems
{ Name Sys_Admittance; NameOfFormulation FEMAdmittance; Type Complex; Frequency 1.; }
If(FEMRemoved)
  For matrix In {0:1}
    { Name Sys_Schur~{matrix}; NameOfFormulation FEMSchur~{matrix}; Type Complex; Frequency 1.; }
  EndFor
EndIf
Return

Macro FEMInitializeMatrices
// Complex initialization precedes every complex runtime expression.
InitSolution[Sys_Admittance];
If(FEMRemoved)
  For matrix In {0:1}
    InitSolution[Sys_Schur~{matrix}];
  EndFor
EndIf
Return

Macro FEMComputeMatrices
For i In {1:NumTerminals}
  row = FEMPermutation(i-1);
  For j In {1:NumTerminals}
    col = FEMPermutation(j-1);
    Evaluate[$FEMWork~{0}~{i}~{j} = $FEMZ~{row}~{col},
             $FEMWork~{1}~{i}~{j} = $FEMP~{row}~{col}];
  EndFor
EndFor
For matrix In {0:1}
  If(ReduceBundle)
    // First transform columns, then rows, using the same bundle basis.
    For j In {1:NumTerminals}
      first = FEMBundleFirst~{j-1}+1;
      If(first > 0)
        For i In {1:NumTerminals}
          Evaluate[$FEMWork~{matrix}~{i}~{j} = $FEMWork~{matrix}~{i}~{j} - $FEMWork~{matrix}~{i}~{first}];
        EndFor
      EndIf
    EndFor
    For i In {1:NumTerminals}
      first = FEMBundleFirst~{i-1}+1;
      If(first > 0)
        For j In {1:NumTerminals}
          Evaluate[$FEMWork~{matrix}~{i}~{j} = $FEMWork~{matrix}~{i}~{j} - $FEMWork~{matrix}~{first}~{j}];
        EndFor
      EndIf
    EndFor
  EndIf
  For j In {1:FEMRetained}
    col = FEMKeep(j-1);
    If(FEMRemoved)
      Evaluate[$FEMAlgebraColumn = col];
      Generate[Sys_Schur~{matrix}]; Solve[Sys_Schur~{matrix}];
      PostOperation[FEMStoreSchur~{matrix}];
    EndIf
    For i In {1:FEMRetained}
      row = FEMKeep(i-1);
      Evaluate[$FEMReduced~{matrix}~{i}~{j} = $FEMWork~{matrix}~{row}~{col}];
      If(FEMRemoved)
        For k In {1:FEMRemoved}
          removed = FEMEliminate(k-1);
          Evaluate[$FEMReduced~{matrix}~{i}~{j} = $FEMReduced~{matrix}~{i}~{j} -
            $FEMWork~{matrix}~{row}~{removed} * $FEMSchurValue~{k}];
        EndFor
      EndIf
    EndFor
  EndFor
  If(IdealTransposition)
    For offset In {0:FEMRetained-1}
      Evaluate[$FEMMean~{offset} = 0.];
      For i In {1:FEMRetained}
        j = (i-1+offset) % FEMRetained + 1;
        Evaluate[$FEMMean~{offset} = $FEMMean~{offset} + $FEMReduced~{matrix}~{i}~{j}/FEMRetained];
      EndFor
    EndFor
    For i In {1:FEMRetained}
      For j In {1:FEMRetained}
        offset = (j-i+FEMRetained) % FEMRetained;
        Evaluate[$FEMReduced~{matrix}~{i}~{j} = $FEMMean~{offset}];
      EndFor
    EndFor
  EndIf
EndFor
For j In {1:FEMRetained}
  Evaluate[$FEMAlgebraColumn = j];
  Generate[Sys_Admittance]; Solve[Sys_Admittance];
  PostOperation[FEMStoreAdmittance~{j}];
EndFor
Evaluate[$FEMInverseResidual = 0.];
For i In {1:FEMRetained}
  For j In {1:FEMRetained}
    Evaluate[$FEMProduct = -(i==j)];
    For k In {1:FEMRetained}
      Evaluate[$FEMProduct = $FEMProduct + $FEMReduced~{1}~{i}~{k} * $FEMY~{k}~{j}];
    EndFor
    // Nonfinite entries must never be published as completed matrices.
    Evaluate[$FEMFinite = Re[$FEMProduct]-Re[$FEMProduct] + Im[$FEMProduct]-Im[$FEMProduct]];
    Test[!($FEMFinite == 0.)] { Error["Nonfinite line-parameter matrix result"]; }
    Evaluate[$FEMInverseResidual = Max[$FEMInverseResidual, Norm[$FEMProduct]]];
  EndFor
EndFor
Print[{$FEMInverseResidual}, Format "Native P*Y-I maximum entry residual: %.17g"];
Return

PostProcessing {
  { Name FEMAdmittanceValues; NameOfFormulation FEMAdmittance; NameOfSystem Sys_Admittance;
    Quantity {
      For i In {1:FEMRetained}
        { Name value~{i}; Value { Term { [{y~{i}}]; In FEMAlgebraCarrier; } } }
      EndFor
    }
  }
  If(FEMRemoved)
    For matrix In {0:1}
      { Name FEMSchurValues~{matrix}; NameOfFormulation FEMSchur~{matrix}; NameOfSystem Sys_Schur~{matrix};
        Quantity {
          For i In {1:FEMRemoved}
            { Name value~{i}; Value { Term { [{x~{i}}]; In FEMAlgebraCarrier; } } }
          EndFor
        }
      }
    EndFor
  EndIf
}
PostOperation {
  For j In {1:FEMRetained}
    { Name FEMStoreAdmittance~{j}; NameOfPostProcessing FEMAdmittanceValues;
      LastTimeStepOnly 1; Operation {
        For i In {1:FEMRetained}
          Print[value~{i}, OnRegion FEMAlgebraCarrier, Format Table, File "", StoreInVariable $FEMY~{i}~{j}];
        EndFor
      }
    }
  EndFor
  If(FEMRemoved)
    For matrix In {0:1}
      { Name FEMStoreSchur~{matrix}; NameOfPostProcessing FEMSchurValues~{matrix};
        LastTimeStepOnly 1; Operation {
          For i In {1:FEMRemoved}
            Print[value~{i}, OnRegion FEMAlgebraCarrier, Format Table, File "", StoreInVariable $FEMSchurValue~{i}];
          EndFor
        }
      }
    EndFor
  EndIf
}
