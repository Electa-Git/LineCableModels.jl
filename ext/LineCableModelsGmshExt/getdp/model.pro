// Preserve the field-model selector; quasi-fw is the supported default.
If(!Exists(Physics))
  Physics = 1;
EndIf
// Immutable inputs for one frequency and its requested terminal excitations.
If(!Exists(ModelDataPath))
  Error("Pass the immutable input path with -setstring ModelDataPath");
EndIf
Include ModelDataPath;
If(#ReceiverInAir() != NumTerminals)
  Error("Provide one voltage reference convention per terminal");
EndIf
eps0 = 8.8541878128e-12;
mu0 = 1.2566370614359173e-6;
If(!Exists(UnitSource))
  UnitSource = 1.0;
EndIf
Include "jacobian.pro";
Include "integration.pro";
If(Physics == 1)
  Include "quasi-full.pro";
Else
  Error("Physics must be 1 (quasi-fw)");
EndIf
