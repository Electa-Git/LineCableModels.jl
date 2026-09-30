// The production native algebra, exercised with independent complex matrices.
DefineConstant[NumTerminals=4,ReduceBundle=0,KronReduction=0,IdealTransposition=0,
               ConnectionCase=0,Singular=0];
TERMINAL = 100;
Connections() = {1:NumTerminals};
If(ConnectionCase == 1) Connections() = {2,1,2,0}; EndIf
Include AlgebraPath;
Resolution { { Name Matrices;
  System { Call FEMMatrixSystems; }
  Operation {
    DeleteFile["matrices.tsv"];
    Call FEMInitializeMatrices;
    For pass In {1:2}
      For i In {1:NumTerminals}
        For j In {1:NumTerminals}
          Evaluate[$FEMP~{i}~{j} = (1-Singular)*pass*Complex[(i==j)*5+.1*i+.2*j,.15*i-.07*j],
                   $FEMZ~{i}~{j} = pass*Complex[(i==j)*7+.2*i-.1*j,.12*i+.06*j]];
        EndFor
      EndFor
      Call FEMComputeMatrices;
      For i In {1:FEMRetained}
        For j In {1:FEMRetained}
          Print[{pass,i,j,Re[$FEMReduced~{0}~{i}~{j}],Im[$FEMReduced~{0}~{i}~{j}],
                 Re[$FEMReduced~{1}~{i}~{j}],Im[$FEMReduced~{1}~{i}~{j}],Re[$FEMY~{i}~{j}],Im[$FEMY~{i}~{j}]},
            Format "%g %g %g %.17g %.17g %.17g %.17g %.17g %.17g", File "matrices.tsv"];
        EndFor
      EndFor
    EndFor
  }
} }
