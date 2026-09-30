// Finite Cartesian coordinate stretching for exp(+j omega t). Physical scalar
// materials remain defined by nu, mu, sigma, epsilon and se. The same side
// stretch is used in both media so the air/earth interface stays matched.
If(!Exists(PmlSlope)) PmlSlope = 1.; EndIf
Function {
  pmlX[] = (Fabs[X[]-Xcenter] > DomainHalfwidth) ?
    (Fabs[X[]-Xcenter]-DomainHalfwidth)/PmlSideThickness : 0.;
  pmlY[] = (Y[] > DomainHalfwidth) ?
    (Y[]-DomainHalfwidth)/PmlTopThickness :
    ((Y[] < -DomainHalfwidth) ? (-Y[]-DomainHalfwidth)/PmlBottomThickness : 0.);
  // Real stretching also damps evanescent modes; the negative imaginary
  // part damps outgoing propagating waves for the positive-time phasor.
  pmlSx[] = 1. + Complex[1., -PmlSlope]*PmlSideStrength*pmlX[]^3;
  pmlSy[] = 1. + Complex[1., -PmlSlope]*(Y[] >= 0. ? PmlTopStrength : PmlBottomStrength)*pmlY[]^3;
  pmlD[] = pmlSx[]*pmlSy[];
  // Native registers evaluate each stretch once within each tensor expression.
  pmlT[] = TensorDiag[(pmlSy[]#1)/(pmlSx[]#0), #0/#1, #0*#1];
  pmlInvT[] = TensorDiag[(pmlSx[]#0)/(pmlSy[]#1), #1/#0, 1./(#0*#1)];
  pmlInvS[] = TensorDiag[1./pmlSx[], 1./pmlSy[], 1.];
  pmlBToPhysical[] = TensorDiag[pmlSx[]/pmlD[], pmlSy[]/pmlD[], 1./pmlD[]];
  nuPml[All] = nu[]*pmlInvT[];
  sePml[All] = se[]*pmlT[];
  seZ[All] = se[]*pmlD[];
  // The coordinate stretch is exactly identity in every physical material.
  nuPml[#{Air, Earth, ConductorMaterialRegions, PassiveMaterialRegions}] = nu[]*TensorDiag[1.,1.,1.];
  sePml[#{Air, Earth, ConductorMaterialRegions, PassiveMaterialRegions}] = se[]*TensorDiag[1.,1.,1.];
  seZ[#{Air, Earth, ConductorMaterialRegions, PassiveMaterialRegions}] = se[];
}
