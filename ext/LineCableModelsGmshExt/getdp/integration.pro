// Quadrature controls are supplied by the native model data.
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
          { GeoElement Triangle; NumberOfPoints PhysicalVolumeQuadrature; }
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
