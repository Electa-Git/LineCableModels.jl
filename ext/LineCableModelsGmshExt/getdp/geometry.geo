// Cartesian exterior around the serialized physical CAD. Dimensions and
// resolution targets have already been evaluated in parameters.pro.
// GEO compares distances relative to the model's bounding-box diagonal. Keep
// the declared physical CAD tolerance when the exterior grows near DC.
FEMBoxWidth = 2*(DomainHalfwidth+PmlSideThickness);
FEMBoxHeight = 2*DomainHalfwidth+PmlTopThickness+PmlBottomThickness;
Geometry.Tolerance = CADGeometryTolerance/Max(1,Sqrt(FEMBoxWidth^2+FEMBoxHeight^2));
Geometry.ToleranceBoolean = 0;
Macro FEMEdgeMesh
  FEMEdgeFirst = Min(FEMEdgeFirst,FEMEdgeLast);
  If(FEMEdgeFirst == FEMEdgeLast)
    FEMEdgeCount = Max(2,Ceil(FEMEdgeLength/FEMEdgeLast)); FEMEdgeRatio = 1;
  Else
    FEMEdgeGap = (FEMEdgeLast-FEMEdgeFirst)/FEMEdgeFirst;
    FEMEdgeLog = FEMEdgeGap < 1e-6 ? FEMEdgeGap*(1-FEMEdgeGap/2+FEMEdgeGap^2/3) : Log(FEMEdgeLast/FEMEdgeFirst);
    FEMEdgeCount = Max(2,Ceil(FEMEdgeLength*FEMEdgeLog/(FEMEdgeLast-FEMEdgeFirst)));
    FEMEdgeRatio = Exp(FEMEdgeLog/(FEMEdgeCount-1));
  EndIf
  If(FEMEdgeReverse) FEMEdgeRatio = 1/FEMEdgeRatio; EndIf
  Transfinite Curve {FEMEdgeCurve} = FEMEdgeCount+1 Using Progression FEMEdgeRatio;
Return

FEMLeft = Xcenter-DomainHalfwidth; FEMRight = Xcenter+DomainHalfwidth;
FEMInnerX() = {FEMLeft,FEMRight};
FEMGridX() = {FEMLeft-PmlSideThickness,FEMInnerX(),FEMRight+PmlSideThickness};
FEMGridY() = {-DomainHalfwidth-PmlBottomThickness,-DomainHalfwidth,0,DomainHalfwidth,DomainHalfwidth+PmlTopThickness};
FEMNx = #FEMGridX(); FEMNy = #FEMGridY();
For ix In {0:FEMNx-1}
  For iy In {0:FEMNy-1}
    FEMGridPoint~{ix}~{iy} = newp;
    Point(FEMGridPoint~{ix}~{iy}) = {FEMGridX(ix),FEMGridY(iy),0,MeshBulk};
  EndFor
EndFor
For ix In {0:FEMNx-2}
  For iy In {0:FEMNy-1}
    // Interior horizontal lines are needed only at the bottom and top.
    If(ix == 0 || ix == FEMNx-2 || iy == 0 || iy == 1 || iy == 3 || iy == 4)
      FEMHorizontal~{ix}~{iy} = newl;
      Line(FEMHorizontal~{ix}~{iy}) = {FEMGridPoint~{ix}~{iy},FEMGridPoint~{ix+1}~{iy}};
      If(ix == 0 || ix == FEMNx-2)
        FEMRatio = Exp(PmlSideGrading/FEMPmlSideLayers);
        If(ix == 0) FEMRatio = 1/FEMRatio; EndIf
        Transfinite Curve {FEMHorizontal~{ix}~{iy}} = FEMPmlSideLayers+1 Using Progression FEMRatio;
      Else
        FEMRemote = iy < 2 ? MeshRemoteEarth : MeshRemoteAir;
        Transfinite Curve {FEMHorizontal~{ix}~{iy}} = Max(2,Ceil((FEMGridX(ix+1)-FEMGridX(ix))/FEMRemote))+1;
      EndIf
    EndIf
  EndFor
EndFor
For ix In {0:FEMNx-1}
  For iy In {0:FEMNy-2}
    If(iy == 0 || iy == 3 || ix == 0 || ix == 1 || ix == FEMNx-2 || ix == FEMNx-1)
      FEMVertical~{ix}~{iy} = newl;
      Line(FEMVertical~{ix}~{iy}) = {FEMGridPoint~{ix}~{iy},FEMGridPoint~{ix}~{iy+1}};
      If(iy == 0)
        Transfinite Curve {FEMVertical~{ix}~{iy}} = FEMPmlBottomLayers+1 Using Progression Exp(-PmlBottomGrading/FEMPmlBottomLayers);
      ElseIf(iy == 3)
        Transfinite Curve {FEMVertical~{ix}~{iy}} = FEMPmlTopLayers+1 Using Progression Exp(PmlTopGrading/FEMPmlTopLayers);
      Else
        FEMEdgeCurve = FEMVertical~{ix}~{iy}; FEMEdgeLength = DomainHalfwidth;
        FEMEdgeLast = iy == 1 ? MeshRemoteEarth : MeshRemoteAir;
        FEMWave = iy == 1 ? MeshWaveEarth : MeshWaveAir;
        FEMEdgeFirst = Min(FEMEdgeLast,Min(MeshBulk,FEMWave));
        FEMEdgeReverse = iy == 1;
        Call FEMEdgeMesh;
      EndIf
    EndIf
  EndFor
EndFor

// One seed source: cable centres, plus the layout's existing central seeds.
FEMInterfaceCandidates() = {CableX(),Xcenter,Xcenter-2,Xcenter+2};
FEMInterfaceCandidates() = Unique[FEMInterfaceCandidates()];
FEMInterfaceTolerance = Max(1e-9,1e-12*(FEMRight-FEMLeft));
FEMInterfaceX() = {FEMLeft}; FEMInterfaceSizes() = {Min(MeshBulk,MeshInterface)};
For i In {0:#FEMInterfaceCandidates()-1}
  FEMX = FEMInterfaceCandidates(i); FEMSize = MeshBulk;
  If(Abs(FEMX-Xcenter) <= FEMInterfaceTolerance) FEMSize = Min(FEMSize,MeshInterface); EndIf
  For cable In {0:NumCables-1}
    If(Abs(FEMX-CableX(cable)) <= FEMInterfaceTolerance) FEMSize = Min(FEMSize,MeshCableInterface~{cable}); EndIf
  EndFor
  FEMLast = #FEMInterfaceX()-1;
  FEMMergeDistance = Max(FEMInterfaceTolerance,.1*Min(FEMSize,FEMInterfaceSizes(FEMLast)));
  If(FEMX-FEMInterfaceX(FEMLast) > FEMMergeDistance && FEMRight-FEMX > Max(FEMInterfaceTolerance,.1*FEMSize))
    FEMInterfaceX() += {FEMX}; FEMInterfaceSizes() += {FEMSize};
  ElseIf(Abs(FEMX-FEMInterfaceX(FEMLast)) <= FEMMergeDistance)
    FEMInterfaceSizes(FEMLast) = Min(FEMInterfaceSizes(FEMLast),FEMSize);
  EndIf
EndFor
FEMInterfaceX() += {FEMRight}; FEMInterfaceSizes() += {Min(MeshBulk,MeshInterface)};
FEMInterfacePoints() = {FEMGridPoint~{1}~{2}};
For i In {1:#FEMInterfaceX()-2}
  FEMPoint = newp; Point(FEMPoint) = {FEMInterfaceX(i),0,0,MeshInterface};
  FEMInterfacePoints() += {FEMPoint};
EndFor
FEMInterfacePoints() += {FEMGridPoint~{FEMNx-2}~{2}};
FEMFiniteInterface() = {};
For i In {0:#FEMInterfacePoints()-2}
  If(FEMInterfaceX(i+1)-FEMInterfaceX(i) < FEMInterfaceTolerance) Error("Interface segment shorter than native tolerance"); EndIf
  FEMCurve = newl; Line(FEMCurve) = {FEMInterfacePoints(i),FEMInterfacePoints(i+1)};
  FEMFiniteInterface() += {FEMCurve};
EndFor
FEMAirOutline() = {FEMVertical~{FEMNx-2}~{2}};
FEMEarthOutline() = {-FEMVertical~{1}~{1}};
For ix In {FEMNx-3:1:-1}
  FEMAirOutline() += {-FEMHorizontal~{ix}~{3}};
EndFor
For ix In {1:FEMNx-3}
  FEMEarthOutline() += {FEMHorizontal~{ix}~{1}};
EndFor
FEMAirOutline() += {-FEMVertical~{1}~{2}};
FEMEarthOutline() += {FEMVertical~{FEMNx-2}~{1}};
FEMLoop = newll; Curve Loop(FEMLoop) = {FEMFiniteInterface(),FEMAirOutline()};
FEMAirSurface = news; Plane Surface(FEMAirSurface) = {FEMLoop,FEMAirHoles()};
FEMReverseInterface() = {};
For i In {#FEMFiniteInterface()-1:0:-1}
  FEMReverseInterface() += {-FEMFiniteInterface(i)};
EndFor
FEMLoop = newll; Curve Loop(FEMLoop) = {FEMReverseInterface(),FEMEarthOutline()};
FEMEarthSurface = news; Plane Surface(FEMEarthSurface) = {FEMLoop,FEMEarthHoles()};
FEMAirPml() = {}; FEMEarthPml() = {}; FEMOuterAir() = {}; FEMOuterEarth() = {};
FEMInterfaceCurves() = FEMFiniteInterface();
For ix In {0:FEMNx-2}
  For iy In {0:FEMNy-2}
    If(ix == 0 || ix == FEMNx-2 || iy == 0 || iy == 3)
      FEMLower = FEMHorizontal~{ix}~{iy}; FEMUpper = FEMHorizontal~{ix}~{iy+1};
      FEMRightCurve = FEMVertical~{ix+1}~{iy}; FEMLeftCurve = FEMVertical~{ix}~{iy};
      FEMLoop = newll; Curve Loop(FEMLoop) = {FEMLower,FEMRightCurve,-FEMUpper,-FEMLeftCurve};
      FEMSurface = news; Plane Surface(FEMSurface) = {FEMLoop};
      Transfinite Surface {FEMSurface} = {FEMGridPoint~{ix}~{iy},FEMGridPoint~{ix+1}~{iy},FEMGridPoint~{ix+1}~{iy+1},FEMGridPoint~{ix}~{iy+1}} AlternateLeft;
      If(PmlQuadrangles) Recombine Surface {FEMSurface}; EndIf
      FEMOuter() = {};
      If(ix == 0) FEMOuter() += {FEMLeftCurve}; EndIf
      If(ix == FEMNx-2) FEMOuter() += {FEMRightCurve}; EndIf
      If(iy == 0) FEMOuter() += {FEMLower}; EndIf
      If(iy == 3) FEMOuter() += {FEMUpper}; EndIf
      If(iy >= 2)
        FEMAirPml() += {FEMSurface}; FEMOuterAir() += {FEMOuter()};
      Else
        FEMEarthPml() += {FEMSurface}; FEMOuterEarth() += {FEMOuter()};
      EndIf
      If(iy == 1 && (ix == 0 || ix == FEMNx-2)) FEMInterfaceCurves() += {FEMUpper}; EndIf
    EndIf
  EndFor
EndFor

// Independent sampling curves have their own points and no surface embedding.
FEMMeasurementPhysical~{0}() = {}; FEMMeasurementPhysical~{1}() = {};
For region In {0:NumPhysicalRegions-1}
  FEMMeasurementBounds~{region}() = BoundingBox Curve {FEMRegionCurves~{region}()};
EndFor
For terminal In {0:NumTerminals-1}
  FEMShift = Min(1e-5,.01*FEMReceiverMetalDimension(terminal));
  FEMX = FEMReceiverX(terminal)+FEMShift;
  // The restricted measurement field supplies the local line targets.
  FEMEnd = newp; Point(FEMEnd) = {FEMX,FEMReceiverY(terminal),0,MeshRemoteMax};
  FEMStartY = 0;
  If(!ReceiverInAir(terminal)) FEMStartY = -DomainHalfwidth; EndIf
  FEMStart = newp; Point(FEMStart) = {FEMX,FEMStartY,0,MeshRemoteMax};
  // Seed every crossed cable neighbourhood. These are independent points on
  // the floating line. These points do not constrain the surface mesh.
  // Without seeds, 1D size integration can miss a narrow exterior field dip.
  FEMMeasureY() = {FEMStartY,FEMReceiverY(terminal)};
  For region In {0:NumPhysicalRegions-1}
    FEMBounds() = FEMMeasurementBounds~{region}();
    If(FEMX >= FEMBounds(0) && FEMX <= FEMBounds(3))
      For side In {0:1}
        FEMSeedY = FEMBounds(1+3*side);
        If(FEMSeedY > FEMStartY+CADGeometryTolerance && FEMSeedY < FEMReceiverY(terminal)-CADGeometryTolerance)
          FEMMeasureY() += {FEMSeedY};
        EndIf
      EndFor
    EndIf
  EndFor
  FEMMeasureY() = Unique(FEMMeasureY());
  FEMMedium = ReceiverInAir(terminal) ? 0 : 1;
  FEMMeasurement() = {}; FEMPrevious = FEMStart;
  For segment In {1:#FEMMeasureY()-1}
    If(segment == #FEMMeasureY()-1) FEMNext = FEMEnd;
    Else FEMNext = newp; Point(FEMNext) = {FEMX,FEMMeasureY(segment),0,MeshRemoteMax}; EndIf
    FEMPhysical = newl; Line(FEMPhysical) = {FEMPrevious,FEMNext};
    FEMMeasurementPhysical~{FEMMedium}() += {FEMPhysical};
    FEMMeasurement() += {FEMPhysical}; FEMPrevious = FEMNext;
  EndFor
  If(!ReceiverInAir(terminal))
    FEMOuter = newp; Point(FEMOuter) = {FEMX,-DomainHalfwidth-PmlBottomThickness,0,MeshRemoteEarth};
    FEMPmlLine = newl; Line(FEMPmlLine) = {FEMOuter,FEMStart};
    Transfinite Curve {FEMPmlLine} = Ceil(FEMPmlBottomLayers/MeasurementLineSizeFactor)+1 Using Progression Exp(-PmlBottomGrading/Ceil(FEMPmlBottomLayers/MeasurementLineSizeFactor));
    FEMMeasurement() += {FEMPmlLine};
  EndIf
  Physical Curve(Sprintf("LCM/measurement_line/%04g",terminal+1),MEASUREMENT_LINE+terminal) = {FEMMeasurement()};
EndFor
Physical Surface("LCM/domain/air",AIR_EM) = {FEMAirSurface};
Physical Surface("LCM/domain/earth",EARTH_EM) = {FEMEarthSurface};
Physical Surface("LCM/domain/air_pml",AIR_PML) = {FEMAirPml()};
Physical Surface("LCM/domain/earth_pml",EARTH_PML) = {FEMEarthPml()};
Physical Surface("LCM/domain/pml",DOMAIN_PML) = {FEMAirPml(),FEMEarthPml()};
Physical Curve("LCM/boundary/magnetic_dirichlet",OUTBND_EM) = {FEMOuterAir(),FEMOuterEarth()};
Physical Curve("LCM/boundary/electric_reference_air",OUTBND_ELE_AIR) = {FEMOuterAir()};
Physical Curve("LCM/boundary/electric_reference_earth",OUTBND_ELE_REF) = {FEMOuterEarth()};
Physical Curve("LCM/boundary/pml_inner",INNER_PML_BND) = {FEMAirOutline(),FEMEarthOutline()};
Physical Curve("LCM/interface/air_earth",INTERFACE_AIR_SOIL) = {FEMInterfaceCurves()};
Physical Surface("LCM/domain/field_maps",6003) = {FEMMaterialSurfaces(),FEMAirSurface,FEMEarthSurface,FEMAirPml(),FEMEarthPml()};
