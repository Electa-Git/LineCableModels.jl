// Prescribed rules only. Zero means inherit the existing triangle rule.
If(!Exists(VolumeQuadrature)) VolumeQuadrature = 12; EndIf
If(!Exists(PhysicalVolumeQuadrature)) PhysicalVolumeQuadrature = 0; EndIf
If(!Exists(PmlQuadrature)) PmlQuadrature = 9; EndIf
If(!Exists(PmlQuadrangles)) PmlQuadrangles = 0; EndIf
PhysicalTrianglePoints = PhysicalVolumeQuadrature;
If(PhysicalTrianglePoints == 0) PhysicalTrianglePoints = VolumeQuadrature; EndIf
Function {
  VolumeRule[All] = 0;
  VolumeRule[Region[{AIR_PML, EARTH_PML}]] = 1;
}
Integration {
  { Name I1; Criterion VolumeRule[];
    Case {
      { Type Gauss;
        Case {
          { GeoElement Point; NumberOfPoints 1; }
          { GeoElement Line; NumberOfPoints 4; }
          { GeoElement Triangle; NumberOfPoints PhysicalTrianglePoints; }
          { GeoElement Quadrangle; NumberOfPoints 4; }
          { GeoElement Triangle2; NumberOfPoints 7; }
        }
      }
      If(PmlQuadrangles)
        { Type GaussLegendre;
          Case { { GeoElement Quadrangle; NumberOfPoints PmlQuadrature; } }
        }
      Else
        { Type Gauss;
          Case { { GeoElement Triangle; NumberOfPoints VolumeQuadrature; } }
        }
      EndIf
    }
  }
  { Name I2;
    Case {
      { Type Gauss;
        Case { { GeoElement Line; NumberOfPoints 4; } }
      }
    }
  }
}
