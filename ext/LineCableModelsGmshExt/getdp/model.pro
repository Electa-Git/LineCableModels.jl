// Helmholtz is the supported field model.
If(!Exists(Physics))
  Physics = 1;
EndIf
// Immutable inputs for one frequency and its requested terminal excitations.
If(!Exists(ModelDataPath))
  Error("Pass the immutable input path with -setstring ModelDataPath");
EndIf
Include ModelDataPath;
Include "parameters.pro";
If(#ReceiverInAir() != NumTerminals)
  Error("Provide one voltage reference convention per terminal");
EndIf
If(!Exists(UnitSource))
  UnitSource = 1.0;
EndIf
Include "jacobian.pro";
Include "integration.pro";
Include "solver.pro";
If(Physics == 1)
  Include "helmholtz.pro";
Else
  Error("Physics must be 1 (Helmholtz)");
EndIf
