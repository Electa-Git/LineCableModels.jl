// Native detached entry. The exported entry sets ModelDataPath and ProjectDirectory.
ExportLineParameters = 1;
If(!Exists(Physics)) Physics = 1; EndIf
PhysicsName = "helmholtz";
ResultRoot = StrCat[ProjectDirectory,"results/"];
RunDirectory = StrCat[ResultRoot,Sprintf["f%04g-",FrequencyIndex],PhysicsName,Sprintf["-b%04g",BasisTerminal]];
MatrixDirectory = StrCat[RunDirectory,"/matrices"];
CompleteFile = StrCat[RunDirectory,"/completed.txt"];
// Invalidate publication even when a subsequent input check fails.
DeleteFile[CompleteFile];
UndefineConstant["GetDP/{Output files"];
FEMPublished = DefineNumber[0, Name "Results/00Status", Choices{0="Not computed",1="Not completed",2="Completed",3="Diagnostic column",4="Mesh ready"}, ReadOnly 1];
If(FrequencyIndex < 1 || FrequencyIndex > FrequencyCount || Floor[FrequencyIndex] != FrequencyIndex)
  Error("Select an exported frequency case");
EndIf
If(BasisTerminal < 0 || BasisTerminal > NumTerminals || Floor[BasisTerminal] != BasisTerminal)
  Error("Basis must be zero (complete matrix) or a terminal index");
EndIf
If(#Connections() != NumTerminals || #TerminalNames() != NumTerminals)
  Error("The connection map and terminal names must match NumTerminals");
EndIf
If(UnitSource == 0)
  Error("Current normalization must be nonzero");
EndIf
Include "parameters.pro";
If(Physics != 1)
  Error("Physics must be 1 (Helmholtz)");
EndIf
If(BasisTerminal == 0) RequestedBases() = {1:NumTerminals};
Else RequestedBases() = {BasisTerminal}; EndIf
FEMPublishedText = DefineString[Str[RunDirectory], Name "Results/01Directory", ReadOnly 1];
FEMPublished = DefineNumber[FrequencyHz, Name "Results/02Frequency [Hz]", ReadOnly 1];
SetString["GetDP/1ResolutionChoices", "LineCableModelsFEM"];
SetString["GetDP/2PostOperationChoices", ""];
If(GetDPVerbosity < -1 || GetDPVerbosity > 5 || Floor[GetDPVerbosity] != GetDPVerbosity)
  Error("GetDPVerbosity must be -1 (stage default) or an integer from zero through five");
EndIf
FEMGetDPVerbosity = GetDPVerbosity >= 0 ? GetDPVerbosity : (RunAction == 0 ? 3 : 4);
SetString["GetDP/9ComputeCommand", StrChoice[RunAction == 0, Sprintf["-v %g",FEMGetDPVerbosity], Sprintf["-solve -v %g -nt %g",FEMGetDPVerbosity,GetDPThreads]]];
If(RunAction == 0 && !StrCmp[OnelabAction,"compute"])
  SetNumber["Results/00Status",4];
EndIf
Include "model.pro";
Include "line-parameters.pro";
FEMOutputScale = OutputTotal ? LineLength : 1.;
FEMMatrixNames() = Str["Z-primitive","P-primitive","Z","P","Y"];
FEMMatrixUnits() = Str["ohm/m","ohm m",StrChoice[OutputTotal,"ohm","ohm/m"],"ohm m",StrChoice[OutputTotal,"S","S/m"]];
// Stable output keys prevent old connection labels or larger reduced matrices
// from leaving stale visible values in the ONELAB tree.
For kind In {0:4}
  For i In {1:NumTerminals}
    For j In {1:NumTerminals}
      FEMVisible = kind < 2 ? (BasisTerminal == 0 || BasisTerminal == j) :
        (BasisTerminal == 0 && i <= FEMRetained && j <= FEMRetained);
      FEMOutputName = StrCat["Results/",Str[FEMMatrixNames(kind)],Sprintf["/entry (%g,%g)",i,j]];
      FEMPublished = DefineNumber[0, Name StrCat[FEMOutputName," real"], Visible FEMVisible, ReadOnly 1];
      FEMPublished = DefineNumber[0, Name StrCat[FEMOutputName," imag"], Visible FEMVisible, ReadOnly 1];
    EndFor
  EndFor
EndFor
FEMPublished = DefineNumber[0, Name "Results/03Inversion residual", Visible (BasisTerminal == 0), ReadOnly 1];

Resolution {
  { Name LineCableModelsFEM;
    System { Call FEMFieldSystem; Call FEMMatrixSystems; }
    Operation {
      If(RunAction)
        PostOperation[FEMClearViews];
        CreateDir[RunDirectory]; CreateDir[MatrixDirectory];
        // A map-disabled rerun must not leave maps from an older solution.
        For t In {1:NumTerminals}
          For q In {0:#FieldMapsHelmholtz()-1}
            DeleteFile[StrCat[MapDirectory,"/",Str[FieldMapsHelmholtz(q)],Sprintf["_f%04g_b%04g.pos",FrequencyIndex,t]]];
          EndFor
        EndFor
        Evaluate[SetNumberRunTime[1]{"Results/00Status"}];
        Call FEMScan;
        If(BasisTerminal == 0)
          Call FEMInitializeMatrices;
          Call FEMComputeMatrices;
        EndIf
        PostOperation[FEMWriteMatrices];
        If(BasisTerminal == 0)
          Evaluate[SetNumberRunTime[$FEMInverseResidual]{"Results/03Inversion residual"}];
        EndIf
        Print[{FrequencyIndex,FrequencyHz,Physics,BasisTerminal},
          Format "%g %.17g %g %g", File CompleteFile];
        Evaluate[SetNumberRunTime[2+(BasisTerminal != 0)]{"Results/00Status"}];
      EndIf
    }
  }
}

PostOperation {
  { Name FEMClearViews; NameOfPostProcessing FEMFields;
    Operation { SendMergeFileRequest[StrCat[ProjectDirectory,"views.geo"]]; }
  }
  { Name FEMWriteMatrices; NameOfPostProcessing FEMFields; LastTimeStepOnly 1;
    Format Table;
    Operation {
      For kind In {0:4}
        If(kind < 2 || BasisTerminal == 0)
          n = kind < 2 ? NumTerminals : FEMRetained;
          MatrixFile = StrCat[MatrixDirectory,"/",Str[FEMMatrixNames(kind)],".tsv"];
          Echo[StrCat["# ",Str[FEMMatrixUnits(kind)],"; exp(+j omega t); response rows, excitation columns; ",Sprintf["f=%.17g Hz",FrequencyHz]],File MatrixFile];
          Echo["row	column	response	excitation	real	imag", File MatrixFile, AppendToExistingFile 1];
          For i In {1:n}
            For j In {1:n}
              If(BasisTerminal == 0 || BasisTerminal == j)
                If(kind < 2)
                  row = i; col = j;
                Else
                  row = FEMPermutation(FEMKeep(i-1)-1);
                  col = FEMPermutation(FEMKeep(j-1)-1);
                EndIf
                // Escape literal percent signs in names before Print's format.
                FEMLabels() = Str[Str[TerminalNames(row-1)],Str[TerminalNames(col-1)]];
                // Labels are escaped once when the post-operation is parsed.
                For label In {0:1}
                  FEMText = Str[FEMLabels(label)]; FEMSafe = "";
                  For c In {0:StrLen[FEMText]-1}
                    letter = StrSub[FEMText,c,1];
                    FEMSafe = StrCat[FEMSafe,StrChoice[!StrCmp[letter,"%"],"%%",letter]];
                  EndFor
                  FEMEscaped~{label} = Str[FEMSafe];
                EndFor
                FEMFormat = StrCat["%g	%g	",FEMEscaped~{0},"	",FEMEscaped~{1},"	%.17g	%.17g"];
                If(kind == 0)
                  Print[{$FEMEntry = $FEMZ~{i}~{j}}, Format "%g", File ""];
                ElseIf(kind == 1)
                  Print[{$FEMEntry = $FEMP~{i}~{j}}, Format "%g", File ""];
                ElseIf(kind == 2)
                  Print[{$FEMEntry = FEMOutputScale*$FEMReduced~{0}~{i}~{j}}, Format "%g", File ""];
                ElseIf(kind == 3)
                  Print[{$FEMEntry = $FEMReduced~{1}~{i}~{j}}, Format "%g", File ""];
                Else
                  Print[{$FEMEntry = FEMOutputScale*$FEMY~{i}~{j}}, Format "%g", File ""];
                EndIf
                Print[{i,j,Re[$FEMEntry],Im[$FEMEntry]}, Format FEMFormat, File MatrixFile, AppendToExistingFile 1];
                FEMOutputName = StrCat["Results/",Str[FEMMatrixNames(kind)],Sprintf["/entry (%g,%g)",i,j]];
                Print[{SetNumberRunTime[Re[$FEMEntry]]{StrCat[FEMOutputName," real"]},
                       SetNumberRunTime[Im[$FEMEntry]]{StrCat[FEMOutputName," imag"]}},Format "%g %g",File ""];
              EndIf
            EndFor
          EndFor
        EndIf
      EndFor
    }
  }
}
