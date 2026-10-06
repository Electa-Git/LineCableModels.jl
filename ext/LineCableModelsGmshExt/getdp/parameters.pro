// Shared real constants for exp(j omega t - Gamma z). Both native parsers
// evaluate these expressions. Julia supplies physical inputs and controls.
// Native measurement sampling control. No conductor-interior size law.
If(!Exists(MeasurementLineSizeFactor)) MeasurementLineSizeFactor = .25; EndIf
If(MeasurementLineSizeFactor <= 0 || MeasurementLineSizeFactor > 1)
  Error("MeasurementLineSizeFactor must be in (0,1]");
EndIf
If(!Exists(FrequencyIndex)) FrequencyIndex = 1; EndIf
If(FrequencyIndex < 1 || FrequencyIndex > FrequencyCount || Floor[FrequencyIndex] != FrequencyIndex)
  Error("Select an exported frequency case");
EndIf
FrequencyHz = Frequencies(FrequencyIndex-1);
GammaRe = GammaReValues(FrequencyIndex-1);
GammaIm = GammaImValues(FrequencyIndex-1);

// Geometry arithmetic and numerical sizing are native responsibilities.
If(FrequencyHz <= 0 || DomainSizeFactor <= 0 || MeshSizeFactor <= 0 || ExteriorMeshSizeFactor < 1 || InterfaceRefinementFactor < 1)
  Error("Frequency and physical mesh controls must be positive; exterior/interface factors must be at least one");
EndIf
// Guard CLI and ONELAB edits as well as Julia-generated controls.
If(PhysicalVolumeQuadrature < 3 || (PhysicalVolumeQuadrature != 3 && PhysicalVolumeQuadrature != 4 && PhysicalVolumeQuadrature != 7 && PhysicalVolumeQuadrature != 12 && PhysicalVolumeQuadrature != 13))
  Error("PhysicalVolumeQuadrature must be 3, 4, 7, 12 or 13 (at least 3)");
EndIf
If(VolumeQuadrature != 4 && VolumeQuadrature != 7 && VolumeQuadrature != 12 && VolumeQuadrature != 13)
  Error("VolumeQuadrature must be 4, 7, 12 or 13");
EndIf
If(PmlQuadrature != 4 && PmlQuadrature != 9 && PmlQuadrature != 16)
  Error("PmlQuadrature must be 4, 9 or 16");
EndIf
If(PmlQuadrangles != 0 && PmlQuadrangles != 1) Error("PmlQuadrangles must be 0 or 1"); EndIf
FEMPositiveFactors() = {DomainSizeFactor,MeshSizeFactor,ExteriorMeshSizeFactor,InterfaceRefinementFactor,
  PmlSideThicknessFactor,PmlTopThicknessFactor,PmlBottomThicknessFactor,
  ConductorGeometryTolerance,ConductorSkinDepthElements,ConductorMeshGrowth,ConductorSkinDepths};
For FEMFactor In {0:#FEMPositiveFactors()-1}
  If(!(FEMPositiveFactors(FEMFactor) > 0 && FEMPositiveFactors(FEMFactor) <= 1.7976931348623157e308))
    Error("Native mesh factors must be finite and positive");
  EndIf
EndFor
FEMOmega = 2*Pi*FrequencyHz;
FEMLayoutRadius = 5;
For cable In {0:NumCables-1}
  FEMLayoutRadius = Max[FEMLayoutRadius,Sqrt[(CableX(cable)-Xcenter)^2+CableY(cable)^2]+CableRadius(cable)];
  For other In {cable+1:NumCables-1}
    FEMLayoutRadius = Max[FEMLayoutRadius,2*Sqrt[(CableX(cable)-CableX(other))^2+(CableY(cable)-CableY(other))^2]];
  EndFor
EndFor
For direction In {0:2}
  FEMLayers() = {PmlSideLayers,PmlTopLayers,PmlBottomLayers};
  FEMGrading() = {PmlSideGrading,PmlTopGrading,PmlBottomGrading};
  FEMThicknessFactor() = {PmlSideThicknessFactor,PmlTopThicknessFactor,PmlBottomThicknessFactor};
  If(FEMLayers(direction) < 1 || FEMLayers(direction) > 2147483646 || FEMLayers(direction) != Floor[FEMLayers(direction)] || !(FEMGrading(direction) >= 0 && FEMGrading(direction) <= 1.7976931348623157e308) || FEMThicknessFactor(direction) <= 0)
    Error("PML intervals must be positive integers; grading nonnegative and relative thickness positive");
  EndIf
EndFor
If(PmlReflection <= 0 || PmlReflection >= 1) Error("PML reflection must lie strictly between zero and one"); EndIf

FEMOmega = 2*Pi*FrequencyHz;
FEMMediumSigma() = {AirSigma, EarthSigma(FrequencyIndex-1)};
FEMMediumEpsilon() = {AirEpsilon, EarthEpsilon(FrequencyIndex-1)};
FEMMediumMu() = {AirMu, EarthMu(FrequencyIndex-1)};
// Crash guard: cap finite-PML strengths at ten times the Gamma=0 calibration.
// This sizing floor does not qualify near-cutoff conductance or change q^2.
PmlSizingFloor = 0.1;
For medium In {0:1}
  FEMBaseRe = -FEMOmega^2*FEMMediumMu(medium)*FEMMediumEpsilon(medium);
  FEMBaseIm = FEMOmega*FEMMediumMu(medium)*FEMMediumSigma(medium);
  FEMBaseAbs = Sqrt[FEMBaseRe^2+FEMBaseIm^2];
  FEMBaseB~{medium} = Sqrt[(FEMBaseAbs-FEMBaseRe)/2];
  FEMBaseA~{medium} = FEMBaseIm/(2*FEMBaseB~{medium});
  FEMRootDRe = FEMBaseRe-GammaRe^2+GammaIm^2;
  FEMRootDIm = FEMBaseIm-2*GammaRe*GammaIm;
  FEMRootDAbs = Sqrt[FEMRootDRe^2+FEMRootDIm^2];
  // Compute the larger root component directly, then obtain the smaller
  // component from 2*a*b = Im(q^2), avoiding cancellation near an axis.
  If(FEMRootDRe >= 0)
    FEMRootA~{medium} = Sqrt[(FEMRootDAbs+FEMRootDRe)/2];
    If(FEMRootA~{medium} == 0)
      FEMRootB~{medium} = 0;
    Else
      FEMRootB~{medium} = FEMRootDIm/(2*FEMRootA~{medium});
    EndIf
  Else
    FEMRootB~{medium} = Sqrt[(FEMRootDAbs-FEMRootDRe)/2];
    // This stable component evaluation initially gives Im(q) >= 0.
    FEMRootA~{medium} = FEMRootDIm/(2*FEMRootB~{medium});
  EndIf
  // Faster prescribed waves use Im(q) >= 0. Otherwise use Re(q) >= 0.
  // Re(q^2) alone cannot distinguish radiation from a slow lossy wave.
  If((GammaIm < FEMBaseB~{medium} && FEMRootB~{medium} < 0) || (GammaIm >= FEMBaseB~{medium} && FEMRootA~{medium} < 0))
    FEMRootA~{medium} = -FEMRootA~{medium};
    FEMRootB~{medium} = -FEMRootB~{medium};
  EndIf
  FEMPmlNeed~{medium} = 0;
  If(FEMRootB~{medium} > 0)
    FEMPmlNeed~{medium} = Min[1,Max[0,(FEMRootB~{medium}-FEMRootA~{medium})/FEMRootB~{medium}]];
  EndIf
  FEMPmlEtaBound~{medium} = 1;
  If(FEMRootB~{medium} < 0)
    FEMPmlEtaBound~{medium} = Min[1,2*FEMRootA~{medium}/(9*Fabs[FEMRootB~{medium}])];
  EndIf
  FEMCutoff~{medium} = FEMRootDAbs == 0;
  FEMBaseMagnitude~{medium} = Sqrt[FEMBaseA~{medium}^2+FEMBaseB~{medium}^2];
  FEMRootMagnitude~{medium} = Sqrt[FEMRootA~{medium}^2+FEMRootB~{medium}^2];
EndFor
// Resolve air waves only when cable fields can reach the interface.
FEMAirExcited = FEMRootA~{1} == 0;
For cable In {0:NumCables-1}
  If(CableY(cable)+CableRadius(cable) > 0)
    FEMAirExcited = 1;
  ElseIf(FEMRootA~{1} > 0)
    If(-CableY(cable)-CableRadius(cable) < 6/FEMRootA~{1}) FEMAirExcited = 1; EndIf
  EndIf
EndFor
// Earth base-root length: true decay length or one wavelength, whichever is shorter.
FEMEarthBaseLength = 2*Pi/FEMBaseMagnitude~{1};
If(FEMBaseA~{1} > 0)
  FEMEarthBaseLength = Min[FEMEarthBaseLength,1/FEMBaseA~{1}];
EndIf
// Safety ceiling for near-lossless quasi-static earth. Physical coefficients remain unchanged.
FEMEarthSizingCeiling = Sqrt[2*1e5/(FEMOmega*EarthMu(FrequencyIndex-1))];
FEMEarthSizingCeilingActive = FEMEarthBaseLength > FEMEarthSizingCeiling;
FEMEarthSizingLength = Min[FEMEarthBaseLength,FEMEarthSizingCeiling];
FEMResolutionRadius = Max[FEMLayoutRadius,FEMEarthSizingLength];
// Imaginary stretching is used only where a participating medium needs it.
// The side uses the same eta in air and earth to keep their interface matched.
PmlSideEta = Min[Max[FEMPmlNeed~{0},FEMPmlNeed~{1}],Min[FEMPmlEtaBound~{0},FEMPmlEtaBound~{1}]];
PmlTopEta = Min[FEMPmlNeed~{0},FEMPmlEtaBound~{0}];
PmlBottomEta = Min[FEMPmlNeed~{1},FEMPmlEtaBound~{1}];
FEMPmlEtas() = {PmlSideEta,PmlTopEta,PmlBottomEta};
// Gamma=0 keeps the earth base-root sizing length. A larger transverse rate
// shortens its numerical length.  A small or cutoff root cannot enlarge it.
FEMEarthDomainLength = FEMEarthSizingLength*(FEMBaseMagnitude~{1}/Max[FEMBaseMagnitude~{1},FEMRootMagnitude~{1}]);
DomainHalfwidth = Max[FEMLayoutRadius,DomainSizeFactor*FEMEarthDomainLength];
PmlSideThickness = PmlSideThicknessFactor*DomainHalfwidth;
PmlTopThickness = PmlTopThicknessFactor*DomainHalfwidth;
PmlBottomThickness = PmlBottomThicknessFactor*DomainHalfwidth;
For medium In {0:1}
  FEMCalibration~{medium} = (FEMBaseA~{medium}+FEMBaseB~{medium})/FEMBaseB~{medium};
  FEMPmlFloorActive~{medium} = 0;
EndFor
For direction In {0:2}
  For medium In {0:1}
    If(direction == 0 || (direction == 1 && medium == 0) || (direction == 2 && medium == 1))
      FEMPmlRate~{direction}~{medium} = (FEMRootA~{medium}+FEMPmlEtas(direction)*FEMRootB~{medium})/FEMCalibration~{medium};
      FEMPmlSizingRate~{direction}~{medium} = Max[FEMPmlRate~{direction}~{medium},PmlSizingFloor*FEMBaseB~{medium}];
    EndIf
  EndFor
EndFor
FEMPmlTarget = -Log[PmlReflection]/2;

// Quasi-static engineering extent cap. The specified rate floor caps added stretch;
// the original geometrical layer thickness remains part of its full length.
FEMCapFactor = 100;
FEMCapRate = FEMPmlTarget/(FEMCapFactor*DomainHalfwidth);
FEMCapQuasistatic = FEMBaseMagnitude~{0}*FEMCapFactor*DomainHalfwidth <= .1;
FEMCapActive() = {0,0,0};
FEMCapThicknesses() = {PmlSideThickness,PmlTopThickness,PmlBottomThickness};
For direction In {0:2}
  FEMCapOldRate = 1e300;
  For medium In {0:1}
    If(direction == 0 || (direction == 1 && medium == 0) || (direction == 2 && medium == 1))
      FEMCapOldRate = Min[FEMCapOldRate,FEMPmlSizingRate~{direction}~{medium}];
    EndIf
  EndFor
  FEMCapActive(direction) = FEMCapQuasistatic && FEMCapOldRate < FEMCapRate;
  If(FEMCapActive(direction))
    FEMPmlEtas(direction) = 0;
    For medium In {0:1}
      If(direction == 0 || (direction == 1 && medium == 0) || (direction == 2 && medium == 1))
        FEMPmlRate~{direction}~{medium} = FEMRootA~{medium}/FEMCalibration~{medium};
        FEMPmlSizingRate~{direction}~{medium} = Max[FEMCapRate,Max[FEMPmlRate~{direction}~{medium},PmlSizingFloor*FEMBaseB~{medium}]];
      EndIf
    EndFor
  EndIf
EndFor
PmlSideEta = FEMPmlEtas(0); PmlTopEta = FEMPmlEtas(1); PmlBottomEta = FEMPmlEtas(2);
For direction In {0:2}
  For medium In {0:1}
    If(direction == 0 || (direction == 1 && medium == 0) || (direction == 2 && medium == 1))
      FEMPmlFloorActive~{direction}~{medium} = FEMPmlRate~{direction}~{medium} < (1 - 1e-9)*PmlSizingFloor*FEMBaseB~{medium};
      FEMPmlFloorActive~{medium} = FEMPmlFloorActive~{medium} || FEMPmlFloorActive~{direction}~{medium};
    EndIf
  EndFor
EndFor


PmlSideStrength = 4*FEMPmlTarget/(PmlSideThickness*Min[FEMPmlSizingRate~{0}~{0},FEMPmlSizingRate~{0}~{1}]);
PmlTopStrength = 4*FEMPmlTarget/(PmlTopThickness*FEMPmlSizingRate~{1}~{0});
PmlBottomStrength = 4*FEMPmlTarget/(PmlBottomThickness*FEMPmlSizingRate~{2}~{1});
For medium In {0:1}
  FEMPmlSideAttenuation~{medium} = PmlSideThickness*PmlSideStrength*(FEMRootA~{medium}+PmlSideEta*FEMRootB~{medium})/4;
  // x_tilde_end = D + L*(1 + A/4 - j*eta*A/4), in outward coordinates.
  FEMPmlSideNetAttenuation~{medium} = (DomainHalfwidth+PmlSideThickness)*FEMRootA~{medium}+FEMPmlSideAttenuation~{medium};
EndFor
FEMPmlTopAttenuation = PmlTopThickness*PmlTopStrength*(FEMRootA~{0}+PmlTopEta*FEMRootB~{0})/4;
FEMPmlBottomAttenuation = PmlBottomThickness*PmlBottomStrength*(FEMRootA~{1}+PmlBottomEta*FEMRootB~{1})/4;
FEMPmlTopNetAttenuation = (DomainHalfwidth+PmlTopThickness)*FEMRootA~{0}+FEMPmlTopAttenuation;
FEMPmlBottomNetAttenuation = (DomainHalfwidth+PmlBottomThickness)*FEMRootA~{1}+FEMPmlBottomAttenuation;
FEMPmlFloorUsed = FEMPmlFloorActive~{0} || FEMPmlFloorActive~{1};

// Resolve the stretched layer until its field has decayed by the target T.
// Native calibrated medium weighting. prescribed interval floors are unchanged.
PmlPointsPerWavelength = 10;
PmlPPWClampMin = 1; PmlPPWClampMax = 3;
For medium In {0:1}
  FEMPmlPPWFactor~{medium} = PmlPPWClampMax;
  If(FEMRootA~{medium} > 0) FEMPmlPPWFactor~{medium} = Min[PmlPPWClampMax,Max[PmlPPWClampMin,Sqrt[.1*FEMRootMagnitude~{medium}/FEMRootA~{medium}]]]; EndIf
  FEMPmlPPW~{medium} = PmlPointsPerWavelength*FEMPmlPPWFactor~{medium}/MeshSizeFactor;
EndFor
For direction In {0:2}
  FEMPmlStrengths() = {PmlSideStrength,PmlTopStrength,PmlBottomStrength};
  FEMPmlThicknesses() = {PmlSideThickness,PmlTopThickness,PmlBottomThickness};
  FEMPmlXRe = FEMPmlThicknesses(direction)*(1+FEMPmlStrengths(direction)/4);
  FEMPmlXIm = -FEMPmlEtas(direction)*FEMPmlThicknesses(direction)*FEMPmlStrengths(direction)/4;
  FEMPmlXAbs = Sqrt[FEMPmlXRe^2+FEMPmlXIm^2];
  FEMPmlMaxPhase = 0;
  For medium In {0:1}
    If((medium != 0 || FEMAirExcited) && (direction == 0 || (direction == 1 && medium == 0) || (direction == 2 && medium == 1)))
      FEMPmlLayerExponent = FEMRootA~{medium}*FEMPmlXRe-FEMRootB~{medium}*FEMPmlXIm;
      FEMPmlPhase = FEMRootMagnitude~{medium}*FEMPmlXAbs;
      If(FEMPmlLayerExponent > 0)
        FEMPmlPhase = FEMPmlPhase*Min[1,FEMPmlTarget/FEMPmlLayerExponent];
      EndIf
      FEMPmlMaxPhase = Max[FEMPmlMaxPhase,FEMPmlPPW~{medium}*FEMPmlPhase];
    EndIf
  EndFor
  FEMPmlEffectiveLayers~{direction} = Max[FEMLayers(direction),Ceil[FEMPmlMaxPhase/(2*Pi)]];
EndFor
FEMPmlSideLayers = FEMPmlEffectiveLayers~{0};
FEMPmlTopLayers = FEMPmlEffectiveLayers~{1};
FEMPmlBottomLayers = FEMPmlEffectiveLayers~{2};

// Physical mesh targets: one native law for managed and detached execution.
MeshBulk = MeshSizeFactor*FEMResolutionRadius/20;
MeshGrowth = 1.2; MeshGrowthSlope = MeshSizeFactor*(MeshGrowth-1);
MeshExteriorStart = 2*FEMResolutionRadius;
MeshRemote = Min[ExteriorMeshSizeFactor*MeshBulk,MeshBulk+MeshGrowthSlope*Max[0,DomainHalfwidth-MeshExteriorStart]];
For medium In {0:1}
  FEMWaveLimit~{medium} = 1e300;
  If(FEMRootMagnitude~{medium} > 0) FEMWaveLimit~{medium} = MeshSizeFactor/(8*FEMRootMagnitude~{medium}); EndIf
  FEMWaveMesh~{medium} = Min[MeshBulk,FEMWaveLimit~{medium}];
  FEMDecayRadius~{medium} = 2*FEMResolutionRadius;
  If(FEMRootA~{medium} > 0) FEMDecayRadius~{medium} = Min[FEMDecayRadius~{medium},6/FEMRootA~{medium}]; EndIf
EndFor
MeshWaveAir = FEMWaveMesh~{0}; MeshWaveEarth = FEMWaveMesh~{1};
MeshDecayAir = FEMDecayRadius~{0}; MeshDecayEarth = FEMDecayRadius~{1};
// A wave field covering the physical box must also constrain its native
// transfinite outer boundary counts. This bound covers the current
// projected-conductor footprint, not a new decay or wavelength prescription.
FEMWaveBoxRadius = 1e300;
For cable In {0:NumCables-1}
  FEMCableWaveFootprint~{cable} = InterfaceRefinementFactor*(Fabs[CableY(cable)]+CableRadius(cable));
  FEMWaveBoxHorizontal = Max[0,Fabs[Xcenter-CableX(cable)]+DomainHalfwidth-FEMCableWaveFootprint~{cable}];
  FEMWaveBoxRadius = Min[FEMWaveBoxRadius,Sqrt[DomainHalfwidth^2+FEMWaveBoxHorizontal^2]];
EndFor
MeshRemoteAir = ExteriorMeshSizeFactor == 1 ? MeshBulk : MeshRemote;
If(FEMAirExcited && FEMWaveBoxRadius <= MeshDecayAir)
  MeshRemoteAir = Min[MeshRemoteAir,MeshWaveAir];
EndIf
MeshRemoteEarth = MeshRemote;
If(FEMWaveBoxRadius <= MeshDecayEarth)
  MeshRemoteEarth = Min[MeshRemoteEarth,MeshWaveEarth];
EndIf
If(!FEMAirExcited) MeshWaveAir = MeshRemoteAir; EndIf
MeshRemoteMax = Max[MeshBulk,Max[MeshRemoteAir,MeshRemoteEarth]];
MeshFine = 1e300;
For region In {0:NumPhysicalRegions-1}
  FEMRegionTarget = 1e300;
  // Several source sections may form one coalesced material region. Native
  // sizing retains the smallest source target without altering its topology.
  For member In {0:#RegionDimensionKinds~{region}()-1}
    FEMKind = RegionDimensionKinds~{region}(member);
    FEMD1 = RegionDimension1~{region}(member); FEMD2 = RegionDimension2~{region}(member);
    FEMRepeated = RegionRepeated~{region}(member);
    If(FEMKind == 1) FEMSize = FEMD1/(FEMRepeated ? 1 : 5);
    ElseIf(FEMKind == 2) FEMSize = (FEMD1-FEMD2)/2;
    ElseIf(FEMKind == 3) FEMSize = Min[FEMD1,FEMD2]/(FEMRepeated ? 2 : 5);
    ElseIf(FEMKind == 4) FEMSize = Min[FEMD1,FEMD2]/(FEMRepeated ? 1 : 5);
    Else FEMSize = FEMD1/FEMD2;
    EndIf
    FEMRegionTarget = Min[FEMRegionTarget,FEMSize];
  EndFor
  RegionSize~{region} = MeshSizeFactor*FEMRegionTarget;
  MeshFine = Min[MeshFine,RegionSize~{region}];
EndFor
MeshInterface = MeshBulk;
For cable In {0:NumCables-1}
  CableSize~{cable} = 1e300;
  For member In {0:#CableSizingRegions~{cable}()-1}
    region = CableSizingRegions~{cable}(member);
    CableSize~{cable} = Min[CableSize~{cable},RegionSize~{region}];
  EndFor
  MeshCableInterface~{cable} = Min[MeshBulk,CableSize~{cable}+MeshGrowthSlope*Max[0,Fabs[CableY(cable)]-CableRadius(cable)]];
  MeshInterface = Min[MeshInterface,MeshCableInterface~{cable}];
EndFor

// Limit the interface layer thickness to avoid buried cable exterior mesh neighbourhoods.
FEMEarthLayerThickness = MeshDecayEarth;
For cable In {0:NumCables-1}
  If(CableY(cable) < 0)
    FEMEarthLayerThickness = Min[FEMEarthLayerThickness,
      .5*(Fabs[CableY(cable)]-CableRadius(cable)-CableSize~{cable})];
  EndIf
EndFor
FEMEarthLayerActive = (MeshWaveEarth < MeshRemoteEarth) && (FEMEarthLayerThickness >= 2*MeshWaveEarth);
If(!FEMEarthLayerActive) FEMEarthLayerThickness = 0; EndIf
FEMEarthLayerClippedOrOmitted = (MeshWaveEarth < MeshRemoteEarth) && (FEMEarthLayerThickness < MeshDecayEarth);

// Derived ONELAB values are observations of the native calculation.
FEMBoundaryValue = DefineNumber[DomainHalfwidth, Name "Boundary/Derived/01Physical half-width [m]", ReadOnly 1];
SetNumber["Boundary/Derived/01Physical half-width [m]",DomainHalfwidth];
FEMBoundaryValue = DefineNumber[PmlSideThickness, Name "Boundary/Derived/02Side thickness [m]", ReadOnly 1];
SetNumber["Boundary/Derived/02Side thickness [m]",PmlSideThickness];
FEMBoundaryValue = DefineNumber[PmlTopThickness, Name "Boundary/Derived/03Top thickness [m]", ReadOnly 1];
SetNumber["Boundary/Derived/03Top thickness [m]",PmlTopThickness];
FEMBoundaryValue = DefineNumber[PmlBottomThickness, Name "Boundary/Derived/04Bottom thickness [m]", ReadOnly 1];
SetNumber["Boundary/Derived/04Bottom thickness [m]",PmlBottomThickness];
FEMBoundaryValue = DefineNumber[PmlSideEta, Name "Boundary/Derived/05Side eta", ReadOnly 1];
SetNumber["Boundary/Derived/05Side eta",PmlSideEta];
FEMBoundaryValue = DefineNumber[PmlTopEta, Name "Boundary/Derived/18Top eta", ReadOnly 1];
SetNumber["Boundary/Derived/18Top eta",PmlTopEta];
FEMBoundaryValue = DefineNumber[PmlBottomEta, Name "Boundary/Derived/19Bottom eta", ReadOnly 1];
SetNumber["Boundary/Derived/19Bottom eta",PmlBottomEta];
FEMBoundaryValue = DefineNumber[PmlSizingFloor, Name "Boundary/Derived/06Sizing-rate floor ratio", ReadOnly 1];
SetNumber["Boundary/Derived/06Sizing-rate floor ratio",PmlSizingFloor];
FEMBoundaryValue = DefineNumber[PmlSideStrength, Name "Boundary/Derived/07Side strength", ReadOnly 1];
SetNumber["Boundary/Derived/07Side strength",PmlSideStrength];
FEMBoundaryValue = DefineNumber[PmlTopStrength, Name "Boundary/Derived/08Top strength", ReadOnly 1];
SetNumber["Boundary/Derived/08Top strength",PmlTopStrength];
FEMBoundaryValue = DefineNumber[PmlBottomStrength, Name "Boundary/Derived/09Bottom strength", ReadOnly 1];
SetNumber["Boundary/Derived/09Bottom strength",PmlBottomStrength];
FEMBoundaryValue = DefineNumber[FEMPmlTopAttenuation, Name "Boundary/Derived/10Added top normal exponent", ReadOnly 1];
SetNumber["Boundary/Derived/10Added top normal exponent",FEMPmlTopAttenuation];
FEMBoundaryValue = DefineNumber[FEMPmlBottomAttenuation, Name "Boundary/Derived/11Added bottom normal exponent", ReadOnly 1];
SetNumber["Boundary/Derived/11Added bottom normal exponent",FEMPmlBottomAttenuation];
FEMBoundaryValue = DefineNumber[FEMPmlTopNetAttenuation, Name "Boundary/Derived/12Net top normal exponent", ReadOnly 1];
SetNumber["Boundary/Derived/12Net top normal exponent",FEMPmlTopNetAttenuation];
FEMBoundaryValue = DefineNumber[FEMPmlBottomNetAttenuation, Name "Boundary/Derived/13Net bottom normal exponent", ReadOnly 1];
SetNumber["Boundary/Derived/13Net bottom normal exponent",FEMPmlBottomNetAttenuation];
FEMBoundaryValue = DefineNumber[FEMPmlSideLayers, Name "Boundary/Derived/15Effective side intervals", ReadOnly 1];
SetNumber["Boundary/Derived/15Effective side intervals",FEMPmlSideLayers];
FEMBoundaryValue = DefineNumber[FEMPmlTopLayers, Name "Boundary/Derived/16Effective top intervals", ReadOnly 1];
SetNumber["Boundary/Derived/16Effective top intervals",FEMPmlTopLayers];
FEMBoundaryValue = DefineNumber[FEMPmlBottomLayers, Name "Boundary/Derived/17Effective bottom intervals", ReadOnly 1];
SetNumber["Boundary/Derived/17Effective bottom intervals",FEMPmlBottomLayers];
FEMBoundaryValue = DefineNumber[PmlPPWClampMin, Name "Boundary/Derived/18PPW clamp minimum", ReadOnly 1];
SetNumber["Boundary/Derived/18PPW clamp minimum",PmlPPWClampMin];
FEMBoundaryValue = DefineNumber[PmlPPWClampMax, Name "Boundary/Derived/19PPW clamp maximum", ReadOnly 1];
SetNumber["Boundary/Derived/19PPW clamp maximum",PmlPPWClampMax];
For medium In {0:1}
  FEMBoundaryValue = DefineNumber[FEMRootA~{medium}, Name Sprintf["Boundary/Derived/Medium %g/01Root real",medium], Label "Root real [1/m]", ReadOnly 1];
  SetNumber[Sprintf["Boundary/Derived/Medium %g/01Root real",medium],FEMRootA~{medium}];
  FEMBoundaryValue = DefineNumber[FEMRootB~{medium}, Name Sprintf["Boundary/Derived/Medium %g/02Root imag",medium], Label "Root imag [1/m]", ReadOnly 1];
  SetNumber[Sprintf["Boundary/Derived/Medium %g/02Root imag",medium],FEMRootB~{medium}];
  FEMBoundaryValue = DefineNumber[FEMCutoff~{medium}, Name Sprintf["Boundary/Derived/Medium %g/03Exact-cutoff flag",medium], ReadOnly 1];
  SetNumber[Sprintf["Boundary/Derived/Medium %g/03Exact-cutoff flag",medium],FEMCutoff~{medium}];
  FEMBoundaryValue = DefineNumber[FEMPmlFloorActive~{medium}, Name Sprintf["Boundary/Derived/Medium %g/04Sizing-floor flag",medium], ReadOnly 1];
  SetNumber[Sprintf["Boundary/Derived/Medium %g/04Sizing-floor flag",medium],FEMPmlFloorActive~{medium}];
  FEMBoundaryValue = DefineNumber[FEMPmlSideAttenuation~{medium}, Name Sprintf["Boundary/Derived/Medium %g/05Added side normal exponent",medium], ReadOnly 1];
  SetNumber[Sprintf["Boundary/Derived/Medium %g/05Added side normal exponent",medium],FEMPmlSideAttenuation~{medium}];
  FEMBoundaryValue = DefineNumber[FEMPmlSideNetAttenuation~{medium}, Name Sprintf["Boundary/Derived/Medium %g/06Net side normal exponent",medium], ReadOnly 1];
  SetNumber[Sprintf["Boundary/Derived/Medium %g/06Net side normal exponent",medium],FEMPmlSideNetAttenuation~{medium}];
EndFor

// A finite strength is not evidence of outgoing attenuation. Publish the
// evaluated values before rejecting a non-attenuating non-cutoff case.
// Exact q=0 is allowed to solve. Its zero exponent is below target and G is
// explicitly unqualified. floor activation alone is informational.
For medium In {0:1}
  If(FEMPmlSideNetAttenuation~{medium} <= 0 && !FEMCutoff~{medium} && !(FEMCapActive(0) && FEMPmlSideNetAttenuation~{medium} == 0))
    Error(Sprintf["PML failure: medium %g net side exponent %.17g is non-positive (f=%.17g, Gamma=%.17g+j*%.17g)",medium,FEMPmlSideNetAttenuation~{medium},FrequencyHz,GammaRe,GammaIm]);
  EndIf
EndFor
If(FEMPmlTopNetAttenuation <= 0 && !FEMCutoff~{0} && !(FEMCapActive(1) && FEMPmlTopNetAttenuation == 0))
  Error(Sprintf["PML failure: air net top exponent %.17g is non-positive (f=%.17g, Gamma=%.17g+j*%.17g)",FEMPmlTopNetAttenuation,FrequencyHz,GammaRe,GammaIm]);
EndIf
If(FEMPmlBottomNetAttenuation <= 0 && !FEMCutoff~{1} && !(FEMCapActive(2) && FEMPmlBottomNetAttenuation == 0))
  Error(Sprintf["PML failure: earth net bottom exponent %.17g is non-positive (f=%.17g, Gamma=%.17g+j*%.17g)",FEMPmlBottomNetAttenuation,FrequencyHz,GammaRe,GammaIm]);
EndIf

FEMEarthSizingValue = DefineNumber[FEMEarthBaseLength, Name "Boundary/Derived/21Earth base-root length [m]", ReadOnly 1];
SetNumber["Boundary/Derived/21Earth base-root length [m]",FEMEarthBaseLength];
FEMEarthSizingValue = DefineNumber[FEMEarthSizingLength, Name "Boundary/Derived/22Earth sizing length [m]", ReadOnly 1];
SetNumber["Boundary/Derived/22Earth sizing length [m]",FEMEarthSizingLength];
FEMEarthSizingValue = DefineNumber[FEMEarthSizingCeilingActive, Name "Boundary/Derived/23Earth sizing ceiling active", ReadOnly 1];
SetNumber["Boundary/Derived/23Earth sizing ceiling active",FEMEarthSizingCeilingActive];

FEMBoundaryValue = DefineNumber[FEMEarthLayerThickness, Name "Mesh/Derived/01Earth interface layer thickness [m]", ReadOnly 1];
SetNumber["Mesh/Derived/01Earth interface layer thickness [m]",FEMEarthLayerThickness];
FEMBoundaryValue = DefineNumber[FEMEarthLayerClippedOrOmitted, Name "Mesh/Derived/02Earth layer clipped or omitted", ReadOnly 1];
SetNumber["Mesh/Derived/02Earth layer clipped or omitted",FEMEarthLayerClippedOrOmitted];


FEMBoundaryValue = DefineNumber[FEMCapActive(0), Name "Boundary/Derived/24Side quasi-static extent cap active", ReadOnly 1];
SetNumber["Boundary/Derived/24Side quasi-static extent cap active",FEMCapActive(0)];

FEMBoundaryValue = DefineNumber[FEMCapActive(1), Name "Boundary/Derived/25Top quasi-static extent cap active", ReadOnly 1];
SetNumber["Boundary/Derived/25Top quasi-static extent cap active",FEMCapActive(1)];

FEMBoundaryValue = DefineNumber[FEMCapActive(2), Name "Boundary/Derived/26Bottom quasi-static extent cap active", ReadOnly 1];
SetNumber["Boundary/Derived/26Bottom quasi-static extent cap active",FEMCapActive(2)];

FEMBoundaryValue = DefineNumber[FEMAirExcited, Name "Mesh/Derived/03Air excited", ReadOnly 1];
SetNumber["Mesh/Derived/03Air excited",FEMAirExcited];
