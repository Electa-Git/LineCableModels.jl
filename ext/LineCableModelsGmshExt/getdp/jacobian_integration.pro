Jacobian {
  { Name Vol;
    Case {
      { Region DomainInf;
        Jacobian VolSphShell{Val_Rint, Val_Rext, Xcenter, Ycenter, Zcenter}; }
      { Region All; Jacobian Vol; }
    }
  }
  { Name Sur;
    Case {
      { Region All; Jacobian Sur; }
    }
  }
}

Integration {
  { Name I1;
    Case {
      { Type Gauss;
        Case {
          { GeoElement Point; NumberOfPoints 1; }
          { GeoElement Line; NumberOfPoints 4; }
          { GeoElement Triangle; NumberOfPoints 4; }
          { GeoElement Quadrangle; NumberOfPoints 4; }
          { GeoElement Triangle2; NumberOfPoints 7; }
        }
      }
    }
  }
}
