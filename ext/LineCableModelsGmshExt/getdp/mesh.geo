// Shared native mesh prescription. Every field uses current physical inputs.
Mesh.MshFileVersion = 4.1; Mesh.Binary = 1; Mesh.SaveAll = 1;
Mesh.MeshSizeMin = Min(MeshFine,Min(MeshWaveAir,MeshWaveEarth));
Mesh.MeshSizeMax = MeshRemoteMax;
Mesh.MeshSizeFromPoints = 1; Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFactor = 1; Mesh.ElementOrder = 1;
// Delaunay resolves the large conductor/earth size gradients.
Mesh.Algorithm = 5;
For region In {0:NumPhysicalRegions-1}
  For member In {0:#FEMRegionPoints~{region}()-1}
    FEMPoint = FEMRegionPoints~{region}(member);
    If(!Exists(FEMPointSize~{FEMPoint})) FEMPointSize~{FEMPoint} = RegionSize~{region};
    Else FEMPointSize~{FEMPoint} = Min(FEMPointSize~{FEMPoint},RegionSize~{region}); EndIf
    Characteristic Length {FEMPoint} = FEMPointSize~{FEMPoint};
  EndFor
EndFor
For i In {1:#FEMInterfacePoints()-2}
  Characteristic Length {FEMInterfacePoints(i)} = FEMInterfaceSizes(i);
EndFor
FEMBackground() = {}; FEMMeasurementCommon() = {}; FEMNextField = 1;
For cable In {0:NumCables-1}
  FEMDistance = FEMNextField; FEMNextField += 1; Field[FEMDistance] = Distance;
  Field[FEMDistance].CurvesList = {FEMCableCurves~{cable}()}; Field[FEMDistance].Sampling = 100;
  FEMField = FEMNextField; FEMNextField += 1; Field[FEMField] = Threshold;
  Field[FEMField].InField = FEMDistance; Field[FEMField].SizeMin = CableSize~{cable};
  Field[FEMField].SizeMax = MeshRemoteMax; Field[FEMField].DistMin = 0;
  Field[FEMField].DistMax = Max(CableSize~{cable},(MeshRemoteMax-CableSize~{cable})/MeshGrowthSlope);
  FEMBackground() += {FEMField}; FEMMeasurementCommon() += {FEMField};
EndFor
FEMConstant = FEMNextField; FEMNextField += 1; Field[FEMConstant] = MathEval; Field[FEMConstant].F = Sprintf("%.17g",MeshBulk);
FEMField = FEMNextField; FEMNextField += 1; Field[FEMField] = Restrict;
Field[FEMField].InField = FEMConstant; Field[FEMField].IncludeBoundary = 1;
Field[FEMField].SurfacesList = {FEMMaterialSurfaces()}; FEMBackground() += {FEMField};
For medium In {0:1}
  FEMMeasurementSources~{medium}() = FEMMeasurementCommon();
  If(medium == 0)
    FEMSurfaces() = {FEMAirSurface,FEMAirPml()}; FEMRemote = MeshRemoteAir; FEMWave = MeshWaveAir; FEMDecay = MeshDecayAir;
  Else
    FEMSurfaces() = {FEMEarthSurface,FEMEarthPml()}; FEMRemote = MeshRemoteEarth; FEMWave = MeshWaveEarth; FEMDecay = MeshDecayEarth;
  EndIf
  FEMRadial = FEMNextField; FEMNextField += 1; Field[FEMRadial] = MathEval;
  Field[FEMRadial].F = Sprintf("Min(%.17g,%.17g+%.17g*Max(0,Sqrt((x-(%.17g))^2+y^2)-%.17g))",FEMRemote,MeshBulk,MeshGrowthSlope,Xcenter,MeshExteriorStart);
  FEMMeasurementSources~{medium}() += {FEMRadial};
  FEMField = FEMNextField; FEMNextField += 1; Field[FEMField] = Restrict;
  Field[FEMField].InField = FEMRadial; Field[FEMField].IncludeBoundary = 1;
  Field[FEMField].SurfacesList = {FEMSurfaces()}; FEMBackground() += {FEMField};
  If(FEMWave < FEMRemote)
    FEMDistance = FEMNextField; FEMNextField += 1; Field[FEMDistance] = Distance;
    Field[FEMDistance].CurvesList = {FEMCableCurves()}; Field[FEMDistance].Sampling = 200;
    FEMSources() = {FEMDistance};
    // The interface layer replaces the earth footprint wave-zone sources.
    If(medium == 0 || !FEMEarthLayerActive)
    For cable In {0:NumCables-1}
      FEMFootprint = FEMNextField; FEMNextField += 1; Field[FEMFootprint] = MathEval;
      Field[FEMFootprint].F = Sprintf("Sqrt(y^2+Max(Abs(x-(%.17g))-(%.17g),0)^2)",CableX(cable),FEMCableWaveFootprint~{cable});
      FEMSources() += {FEMFootprint};
    EndFor
    EndIf
    FEMDistance = FEMNextField; FEMNextField += 1; Field[FEMDistance] = Min; Field[FEMDistance].FieldsList = {FEMSources()};
    // Grow the lossy-earth exterior wave target geometrically with distance.
    If(medium == 1 && FEMRootA~{1} > 0)
    FEMWaveField = FEMNextField; FEMNextField += 1; Field[FEMWaveField] = MathEval;
    // MathEval is not reentrant: inline footprints instead of evaluating their fields.
    If(FEMEarthLayerActive)
      // This minimum contains only the cable Distance field.
      Field[FEMWaveField].F = Sprintf("Min(%.17g,%.17g*Exp(F%g*%.17g/2))",FEMRemote,FEMWave,FEMDistance,FEMRootA~{1});
    Else
      FEMWaveDistanceExpression = Sprintf("F%g",FEMSources(0));
      For cable In {0:NumCables-1}
        FEMWaveDistanceExpression = StrCat["Min(",FEMWaveDistanceExpression,
          Sprintf(",Sqrt(y^2+Max(Abs(x-(%.17g))-(%.17g),0)^2))",CableX(cable),FEMCableWaveFootprint~{cable})];
      EndFor
      Field[FEMWaveField].F = StrCat[Sprintf("Min(%.17g,%.17g*Exp((",FEMRemote,FEMWave),
        FEMWaveDistanceExpression,Sprintf(")*%.17g/2))",FEMRootA~{1})];
    EndIf
    Else
    FEMWaveField = FEMNextField; FEMNextField += 1; Field[FEMWaveField] = Threshold;
    Field[FEMWaveField].InField = FEMDistance; Field[FEMWaveField].SizeMin = FEMWave;
    Field[FEMWaveField].SizeMax = FEMRemote; Field[FEMWaveField].DistMin = FEMDecay;
    Field[FEMWaveField].DistMax = 2*FEMDecay;
    EndIf
    FEMMeasurementSources~{medium}() += {FEMWaveField};
    FEMField = FEMNextField; FEMNextField += 1; Field[FEMField] = Restrict; Field[FEMField].InField = FEMWaveField;
    Field[FEMField].IncludeBoundary = 1; Field[FEMField].SurfacesList = {FEMSurfaces()}; FEMBackground() += {FEMField};
  EndIf
EndFor
If(ConductorGeometryTolerance <= 0 || ConductorSkinDepthElements <= 0 || ConductorMeshGrowth < 1 || ConductorSkinDepths <= 0 || ConductorThicknessElements < 1 || ConductorThicknessElements != Floor(ConductorThicknessElements))
  Error("Conductor mesh controls must be positive; growth must be at least one");
EndIf
// One resolution scale also controls angular, first-normal and interior targets.
FEMCircleSegments = 12*2^Max(0,Ceil(Log(Sqrt(2*Pi^2/(3*ConductorGeometryTolerance))/12)/Log(2)));
For region In {0:NumPhysicalRegions-1}
  If(FEMConductor~{region})
    Mesh.MeshSizeMin = 0;
    material = FEMConductorMaterial~{region};
    FEMSigma = MaterialSigma~{material}(FrequencyIndex-1);
    FEMB = FEMOmega*MaterialEpsilon~{material}(FrequencyIndex-1);
    FEMH = Sqrt(FEMSigma^2+FEMB^2); FEMMu = MaterialMu(material-1);
    If(FEMB >= 0 && FEMH+FEMB > 0) FEMDecay = Abs(FEMSigma)*Sqrt(FEMOmega*FEMMu/(2*(FEMH+FEMB)));
    Else FEMDecay = Sqrt(FEMOmega*FEMMu*(FEMH-FEMB)/2); EndIf
    FEMPhase = Sqrt(FEMOmega*FEMMu*(FEMH+FEMB)/2);
    FEMDelta = 1e300; If(FEMDecay > 0) FEMDelta = 1/FEMDecay; EndIf
    FEMWidth = FEMConductorWidth~{region};
    FEMContourSize = RegionSize~{region};
    FEMDivisions = FEMConductorKind~{region} == 2 ? ConductorThicknessElements : 5;
    FEMFraction = FEMConductorKind~{region} == 1 ? .8 : (FEMConductorKind~{region} == 2 ? .45 : .9);
    FEMCap = MeshSizeFactor*FEMWidth/FEMDivisions; FEMBulk = Min(RegionSize~{region},FEMCap);
    If(FEMDecay*FEMWidth <= ConductorSkinDepths && FEMPhase > 0) FEMBulk = Min(FEMBulk,MeshSizeFactor*2*Pi/(12*FEMPhase)); EndIf
    FEMExtent = Min(ConductorSkinDepths*FEMDelta,FEMFraction*FEMWidth);
    FEMFirst = Min(MeshSizeFactor*FEMDelta/ConductorSkinDepthElements,FEMCap);
    If(FEMConductorSector~{region})
      FEMBulk = Min(FEMBulk,MeshSizeFactor*FEMWidth/30); FEMSegments = 2*FEMCircleSegments; FEMDepth = FEMSectorDepth~{region};
      If(ConductorMeshGrowth == 1) FEMLayers = Ceil(FEMDepth/FEMFirst);
      Else FEMLayers = Ceil(Log(1+(ConductorMeshGrowth-1)*FEMDepth/FEMFirst)/Log(ConductorMeshGrowth)); EndIf
      For arc In {0:#FEMSectorOuter~{region}()-1}
        FEMAngular = FEMSectorTurns~{region}(arc)*FEMSegments/(2*Pi);
        FEMTangential = FEMSectorLengths~{region}(arc)/((FEMWidth/5)*96/FEMSegments);
        Transfinite Curve {FEMSectorOuter~{region}(arc),FEMSectorInner~{region}(arc)} = Max(2,Max(Ceil(FEMAngular),Ceil(FEMTangential)))+1;
      EndFor
      For spoke In {0:#FEMSectorSpokes~{region}()-1}
        FEMCurve = FEMSectorSpokes~{region}(spoke);
        FEMRatio = FEMCurve > 0 ? ConductorMeshGrowth : 1/ConductorMeshGrowth;
        Transfinite Curve {Abs(FEMCurve)} = FEMLayers+1 Using Progression FEMRatio;
      EndFor
    Else
      If(ConductorMeshGrowth == 1)
        FEMLayers = Ceil(FEMExtent/FEMFirst); FEMFirst = FEMExtent/FEMLayers;
      Else
        FEMLayers = Ceil(Log(1+(ConductorMeshGrowth-1)*FEMExtent/FEMFirst)/Log(ConductorMeshGrowth));
        FEMFirst = FEMExtent*(ConductorMeshGrowth-1)/(ConductorMeshGrowth^FEMLayers-1);
      EndIf
      For arc In {0:#FEMRegionCurves~{region}()-1}
        // Allow 64 Float64 ulps when rounding angular interval counts.
        FEMAngularIntervals = FEMCircleSegments*FEMConductorArcFractions~{region}(arc);
        FEMAngularEpsilon = 0;
        If(FEMAngularIntervals > 0)
          FEMAngularEpsilon = 2^(Floor(Log(FEMAngularIntervals)/Log(2))-52);
        EndIf
        FEMAngular = Ceil(FEMAngularIntervals-64*FEMAngularEpsilon)+1;
        FEMLocal = Ceil(FEMConductorArcLengths~{region}(arc)/RegionSize~{region})+1;
        FEMArcNodes = Max(2,Max(FEMAngular,FEMLocal));
        Transfinite Curve {FEMRegionCurves~{region}(arc)} = FEMArcNodes;
        FEMContourSize = Min(FEMContourSize,FEMConductorArcLengths~{region}(arc)/(FEMArcNodes-1));
      EndFor
      FEMLayer = FEMNextField; FEMNextField += 1; Field[FEMLayer] = BoundaryLayer;
      Field[FEMLayer].CurvesList = {FEMRegionCurves~{region}()};
      Field[FEMLayer].Size = FEMFirst; Field[FEMLayer].Ratio = ConductorMeshGrowth;
      Field[FEMLayer].Thickness = FEMExtent*(1+1e-8); Field[FEMLayer].Quads = 0;
      FEMExcluded() = Surface{:};
      If(FEMDelta < FEMWidth)
        FEMExcluded() -= {FEMRegionSurfaces~{region}()};
        Field[FEMLayer].ExcludedSurfacesList = {FEMExcluded()}; BoundaryLayer Field = FEMLayer;
      Else Field[FEMLayer].ExcludedSurfacesList = {FEMExcluded()}; EndIf
    EndIf
    FEMConstant = FEMNextField; FEMNextField += 1; Field[FEMConstant] = MathEval; Field[FEMConstant].F = Sprintf("%.17g",FEMBulk);
    FEMField = FEMNextField; FEMNextField += 1; Field[FEMField] = Restrict; Field[FEMField].InField = FEMConstant;
    Field[FEMField].SurfacesList = {FEMConductorBulkSurfaces~{region}()}; FEMBackground() += {FEMField};
    // Reuse the native transfinite contour spacing and exterior growth law.
    // Boundary-size Extend fields are unavailable while floating 1D curves
    // are meshed, so this target is declared from the same arc counts above.
    FEMContourDistance = FEMNextField; FEMNextField += 1; Field[FEMContourDistance] = Distance;
    Field[FEMContourDistance].CurvesList = {FEMRegionCurves~{region}()}; Field[FEMContourDistance].Sampling = 200;
    FEMContourTarget = FEMNextField; FEMNextField += 1; Field[FEMContourTarget] = Threshold;
    Field[FEMContourTarget].InField = FEMContourDistance;
    Field[FEMContourTarget].SizeMin = FEMContourSize; Field[FEMContourTarget].SizeMax = MeshRemoteMax;
    Field[FEMContourTarget].DistMin = 0; Field[FEMContourTarget].DistMax = Max(FEMContourSize,(MeshRemoteMax-FEMContourSize)/MeshGrowthSlope);
    FEMMeasurementSources~{0}() += {FEMContourTarget}; FEMMeasurementSources~{1}() += {FEMContourTarget};
    // These passive constants are also retained for mesh-parity inspection.
    FEMConductorDelta~{region} = FEMDelta; FEMConductorFirst~{region} = FEMFirst;
    FEMConductorExtent~{region} = FEMExtent; FEMConductorBulk~{region} = FEMBulk;
  EndIf
EndFor
For cable In {0:NumCables-1}
  If(#FEMCablePassiveSurfaces~{cable}())
    FEMExtension = FEMNextField; FEMNextField += 1; Field[FEMExtension] = Extend;
    Field[FEMExtension].CurvesList = {FEMCablePassiveCurves~{cable}()};
    Field[FEMExtension].SizeMax = CableSize~{cable}; Field[FEMExtension].DistMax = CableSize~{cable}/MeshGrowthSlope;
    FEMField = FEMNextField; FEMNextField += 1; Field[FEMField] = Restrict; Field[FEMField].InField = FEMExtension;
    Field[FEMField].SurfacesList = {FEMCablePassiveSurfaces~{cable}()}; FEMBackground() += {FEMField};
  EndIf
EndFor
// Earth-normal skin resolution, independent of the horizontal footprint law.
If(FEMEarthLayerActive)
  FEMEarthLayer = FEMNextField; FEMNextField += 1; Field[FEMEarthLayer] = BoundaryLayer;
  Field[FEMEarthLayer].CurvesList = {FEMFiniteInterface()};
  Field[FEMEarthLayer].Size = MeshWaveEarth; Field[FEMEarthLayer].Ratio = 1.4;
  Field[FEMEarthLayer].Thickness = FEMEarthLayerThickness; Field[FEMEarthLayer].Quads = 0;
  FEMExcluded() = Surface{:}; FEMExcluded() -= {FEMEarthSurface};
  Field[FEMEarthLayer].ExcludedSurfacesList = {FEMExcluded()};
  BoundaryLayer Field = FEMEarthLayer;
  FEMLayerSize = FEMNextField; FEMNextField += 1; Field[FEMLayerSize] = MathEval;
  Field[FEMLayerSize].F = Sprintf("Min(%.17g,%.17g+%.17g*Max(0,Abs(y)-%.17g))",MeshRemoteEarth,MeshWaveEarth,MeshGrowthSlope,FEMEarthLayerThickness);
  FEMMeasurementSources~{1}() += {FEMLayerSize};
  // The shared interface also refines the adjacent air triangles. Floating
  // air lines must resolve that boundary neighbourhood without sizing surfaces.
  FEMInterfaceLineSize = FEMNextField; FEMNextField += 1; Field[FEMInterfaceLineSize] = MathEval;
  Field[FEMInterfaceLineSize].F = Sprintf("Min(%.17g,%.17g+%.17g*Abs(y))",MeshRemoteAir,MeshWaveEarth,MeshGrowthSlope);
  FEMMeasurementSources~{0}() += {FEMInterfaceLineSize};
EndIf
// Floating curves sample the minimum of all exterior targets; these fields never
// participate on surface entities. PML curves retain their transfinite law.
For medium In {0:1}
  If(#FEMMeasurementPhysical~{medium}() > 0)
    FEMLocal = FEMNextField; FEMNextField += 1; Field[FEMLocal] = Min;
    Field[FEMLocal].FieldsList = {FEMMeasurementSources~{medium}()};
    FEMHalf = FEMNextField; FEMNextField += 1; Field[FEMHalf] = Threshold;
    Field[FEMHalf].InField = FEMLocal; Field[FEMHalf].DistMin = 0;
    Field[FEMHalf].DistMax = MeshRemoteMax; Field[FEMHalf].SizeMin = 0;
    Field[FEMHalf].SizeMax = MeasurementLineSizeFactor*MeshRemoteMax; Field[FEMHalf].Sigmoid = 0;
    FEMField = FEMNextField; FEMNextField += 1; Field[FEMField] = Restrict;
    Field[FEMField].InField = FEMHalf; Field[FEMField].IncludeBoundary = 0;
    Field[FEMField].CurvesList = {FEMMeasurementPhysical~{medium}()};
    FEMBackground() += {FEMField};
  EndIf
EndFor
FEMField = FEMNextField; FEMNextField += 1; Field[FEMField] = Min; Field[FEMField].FieldsList = {FEMBackground()};
Background Field = FEMField;

// Native values are recorded for inspection; Julia never reconstructs them.
If(Exists(MeshMetadataPath))
  Printf("pml_ppw_clamp_bounds %.17g %.17g",PmlPPWClampMin,PmlPPWClampMax) >> MeshMetadataPath;
  Printf("pml_medium_points_per_wavelength %.17g %.17g",FEMPmlPPW~{0},FEMPmlPPW~{1}) >> MeshMetadataPath;
  Printf("physical_domain_halfwidth_m %.17g",DomainHalfwidth) > MeshMetadataPath;
  Printf("measurement_line_size_factor %.17g",MeasurementLineSizeFactor) >> MeshMetadataPath;
  Printf("pml_thickness_m %.17g %.17g %.17g",PmlSideThickness,PmlTopThickness,PmlBottomThickness) >> MeshMetadataPath;
  Printf("minimum_mesh_size_m %.17g",Mesh.MeshSizeMin) >> MeshMetadataPath;
  Printf("domain_mesh_size_m %.17g",MeshBulk) >> MeshMetadataPath;
  Printf("earth_interface_layer %.17g %g",FEMEarthLayerThickness,FEMEarthLayerClippedOrOmitted) >> MeshMetadataPath;
  Printf("air_excited %g",FEMAirExcited) >> MeshMetadataPath;
  Printf("interface_mesh_size_m %.17g",MeshInterface) >> MeshMetadataPath;
  Printf("exterior_mesh_sizes_m %.17g %.17g",MeshRemoteAir,MeshRemoteEarth) >> MeshMetadataPath;
  Printf("exterior_start_radius_m %.17g",MeshExteriorStart) >> MeshMetadataPath;
  Printf("wave_mesh_sizes_m %.17g %.17g",MeshWaveAir,MeshWaveEarth) >> MeshMetadataPath;
  Printf("wave_decay_radii_m %.17g %.17g",MeshDecayAir,MeshDecayEarth) >> MeshMetadataPath;
  Printf("adjacent_growth_factor %.17g",MeshGrowth) >> MeshMetadataPath;
  Printf("pml_eta %.17g %.17g %.17g",PmlSideEta,PmlTopEta,PmlBottomEta) >> MeshMetadataPath;
  Printf("pml_strength %.17g %.17g %.17g",PmlSideStrength,PmlTopStrength,PmlBottomStrength) >> MeshMetadataPath;
  Printf("transverse_roots %.17g %.17g %.17g %.17g",FEMRootA~{0},FEMRootB~{0},FEMRootA~{1},FEMRootB~{1}) >> MeshMetadataPath;
  Printf("cutoff_flags %g %g",FEMCutoff~{0},FEMCutoff~{1}) >> MeshMetadataPath;
  Printf("sizing_floor_flags %g %g",FEMPmlFloorActive~{0},FEMPmlFloorActive~{1}) >> MeshMetadataPath;
  Printf("pml_effective_layers %g %g %g",FEMPmlSideLayers,FEMPmlTopLayers,FEMPmlBottomLayers) >> MeshMetadataPath;
  Printf("pml_side_attenuation %.17g %.17g",FEMPmlSideAttenuation~{0},FEMPmlSideAttenuation~{1}) >> MeshMetadataPath;
  Printf("pml_top_attenuation %.17g",FEMPmlTopAttenuation) >> MeshMetadataPath;
  Printf("pml_bottom_attenuation %.17g",FEMPmlBottomAttenuation) >> MeshMetadataPath;
  Printf("pml_side_net_attenuation %.17g %.17g",FEMPmlSideNetAttenuation~{0},FEMPmlSideNetAttenuation~{1}) >> MeshMetadataPath;
  Printf("pml_top_net_attenuation %.17g",FEMPmlTopNetAttenuation) >> MeshMetadataPath;
  Printf("pml_bottom_net_attenuation %.17g",FEMPmlBottomNetAttenuation) >> MeshMetadataPath;
  For region In {0:NumPhysicalRegions-1}
    Printf("region_mesh_sizes_m %.17g",RegionSize~{region}) >> MeshMetadataPath;
  EndFor
  For cable In {0:NumCables-1}
    Printf("cable_outer_mesh_sizes_m %.17g",CableSize~{cable}) >> MeshMetadataPath;
    Printf("cable_interface_mesh_sizes_m %.17g",MeshCableInterface~{cable}) >> MeshMetadataPath;
  EndFor
  For region In {0:NumPhysicalRegions-1}
    If(FEMConductor~{region})
      Printf("conductor_targets %.17g %.17g %.17g %.17g",FEMConductorDelta~{region},FEMConductorFirst~{region},FEMConductorExtent~{region},FEMConductorBulk~{region}) >> MeshMetadataPath;
    Else Printf("conductor_targets null") >> MeshMetadataPath; EndIf
  EndFor
EndIf
