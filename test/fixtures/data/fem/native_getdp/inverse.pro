// A small algebraic system P y_j = e_j, solved by GetDP's own linear solver.
// Four complex, nonsymmetric rows demonstrate this is not a 3x3 Tensor inverse.
N = 4;
Group { Carrier = Region[100]; }
Function {
 For i In {1:N}
   For j In {1:N}
     p~{i}~{j}[] = Complex[(i==j)*5 + .1*i + .2*j, .15*i-.07*j];
   EndFor
 EndFor
}
FunctionSpace {
 For i In {1:N}
 {Name H~{i}; Type Form0;
   BasisFunction {
     {Name s; NameOfCoef a; Function BF_GroupOfNodes;
      Support Carrier; Entity GroupsOfNodesOf[Carrier];}
   }
   GlobalQuantity { {Name value; Type AliasOf; NameOfCoef a;} }
 }
 EndFor
}
Formulation { {Name Algebra; Type FemEquation;
 Quantity {
   For i In {1:N}
     {Name y~{i}; Type Global; NameOfSpace H~{i}[value];}
   EndFor
 }
 Equation {
   For i In {1:N}
     For j In {1:N}
       GlobalTerm { [p~{i}~{j}[] * Dof{y~{j}}, {y~{i}}]; In Carrier; }
     EndFor
     GlobalTerm { [-($Column==i), {y~{i}}]; In Carrier; }
   EndFor
 }
} }
Resolution { {Name Invert; System {
 {Name Matrix; NameOfFormulation Algebra; Type ComplexValue; Frequency 1;}
} Operation {
 DeleteFile["arbitrary.txt"];
 For column In {1:N}
   Evaluate[$Column=column];
   Generate[Matrix]; Solve[Matrix]; SaveSolution[Matrix];
   PostOperation[Column~{column}];
 EndFor
} } }
PostProcessing { {Name Values; NameOfFormulation Algebra; Quantity {
 For i In {1:N}
 {Name value~{i}; Value { Term { [{y~{i}}]; In Carrier; } } }
 EndFor
} } }
PostOperation {
 For column In {1:N}
 {Name Column~{column}; NameOfPostProcessing Values; LastTimeStepOnly 1;
   Operation {
     For i In {1:N}
       Print[value~{i}, OnRegion Carrier, Format Table, File "", StoreInVariable $Entry];
       Print[{i,column,Re[$Entry],Im[$Entry]}, Format "%g %g %.17g %.17g",
         File "arbitrary.txt", AppendToExistingFile 1];
     EndFor
   }
 }
 EndFor
}
