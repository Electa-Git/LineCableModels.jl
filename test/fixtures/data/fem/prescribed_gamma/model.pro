// Independent manufactured Maxwell problem, developed before production changes.
// Unit square, exp(j omega t - Gamma z), normalized transverse potentials.
// The exact potentials are a=sin(pi*x)sin(pi*y), v=nu*a/kappa,
// b=(-cos(pi*x)sin(pi*y), sin(pi*x)cos(pi*y)). Compare E and B, not gauge.
If(!Exists(GammaRe)) GammaRe = 0.; EndIf
If(!Exists(GammaIm)) GammaIm = 0.; EndIf
If(!Exists(StretchX)) StretchX = 0.; EndIf
If(!Exists(StretchY)) StretchY = 0.; EndIf
Group { Domain = Region[1]; Boundary = Region[2]; }
Function {
  gamma[] = Complex[GammaRe,GammaIm]; g2[] = gamma[]*gamma[];
  nu[] = 0.8; omega[] = 1.7; se[] = Complex[0.7,0.4*1.7];
  sx[] = Complex[1,StretchX]; sy[] = Complex[1,StretchY];
  detS[] = sx[]*sy[]; tx[] = sy[]/sx[]; ty[] = sx[]/sy[];
  nuTensor[] = nu[]*TensorDiag[1/tx[],1/ty[],1/detS[]];
  seTensor[] = se[]*TensorDiag[tx[],ty[],detS[]];
  aExact[] = Sin[Pi*X[]]*Sin[Pi*Y[]];
  gradA[] = Pi*Vector[Cos[Pi*X[]]*Sin[Pi*Y[]],Sin[Pi*X[]]*Cos[Pi*Y[]],0];
  bExact[] = Vector[-Cos[Pi*X[]]*Sin[Pi*Y[]]/tx[],Sin[Pi*X[]]*Cos[Pi*Y[]]/ty[],0];
  curlB[] = Pi*(1/tx[]+1/ty[])*Cos[Pi*X[]]*Cos[Pi*Y[]];
  lambda[] = Pi^2*nu[]*(1/sx[]^2+1/sy[]^2)+Complex[0,omega[]]*se[]-g2[]*nu[];
  fa[] = detS[]*lambda[]*aExact[];
  fb[] = lambda[]*TensorDiag[tx[],ty[],detS[]]*bExact[];
  etExact[] = -nu[]/se[]*gradA[]-Complex[0,omega[]]*bExact[];
  ezExact[] = (-Complex[0,omega[]]+g2[]*nu[]/se[])*aExact[];
  btExact[] = -(Vector[0,0,1] /\ (gradA[]+g2[]*bExact[]));
}
Jacobian { { Name Vol; Case { { Region All; Jacobian Vol; } } } }
Integration { { Name I1; Case { { Type Gauss; Case {
  { GeoElement Triangle; NumberOfPoints 12; }
} } } } }
Constraint {
  { Name Zero; Case { { Region Boundary; Value 0.; } } }
  { Name Gauge; Type Assign; Case {
    { Region Domain; SubRegion Boundary; Value 0.; }
  } }
}
FunctionSpace {
  { Name Az; Type Form1P;
    BasisFunction { { Name s; NameOfCoef a; Function BF_PerpendicularEdge;
      Support Domain; Entity NodesOf[All]; } }
    Constraint { { NameOfCoef a; EntityType NodesOf; NameOfConstraint Zero; } }
  }
  { Name Bt; Type Form1;
    BasisFunction { { Name s; NameOfCoef b; Function BF_Edge;
      Support Domain; Entity EdgesOf[All]; } }
    Constraint {
      { NameOfCoef b; EntityType EdgesOf; NameOfConstraint Zero; }
      { NameOfCoef b; EntityType EdgesOfTreeIn; EntitySubType StartingOn;
        NameOfConstraint Gauge; }
    }
  }
  { Name V; Type Form0;
    BasisFunction { { Name s; NameOfCoef v; Function BF_Node;
      Support Domain; Entity NodesOf[All]; } }
    Constraint { { NameOfCoef v; EntityType NodesOf; NameOfConstraint Zero; } }
  }
}
Formulation {
  { Name Maxwell; Type FemEquation;
    Quantity {
      { Name a; Type Local; NameOfSpace Az; }
      { Name bt; Type Local; NameOfSpace Bt; }
      { Name v; Type Local; NameOfSpace V; }
    }
    Equation {
      Galerkin { [nuTensor[]*Dof{d a}, {d a}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [-g2[]*nuTensor[]*(Vector[0,0,1] /\ Dof{bt}), {d a}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0,omega[]]*se[]*detS[]*Dof{a}, {a}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [-g2[]*se[]*detS[]*Dof{v}*Vector[0,0,1], {a}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [-fa[]*Vector[0,0,1], {a}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [nuTensor[]*Dof{d bt}, {d bt}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0,omega[]]*seTensor[]*Dof{bt}, {bt}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [g2[]*(Vector[0,0,1] /\ (nuTensor[]*(Vector[0,0,1] /\ Dof{bt}))), {bt}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [seTensor[]*Dof{d v}, {bt}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [-(Vector[0,0,1] /\ (nuTensor[]*Dof{d a})), {bt}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [-fb[], {bt}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [seTensor[]*Dof{d v}, {d v}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0,omega[]]*seTensor[]*Dof{bt}, {d v}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0,omega[]]*se[]*detS[]*(Dof{a}*Vector[0,0,1]), {v}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [-g2[]*se[]*detS[]*Dof{v}, {v}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [-fa[], {v}]; In Domain; Jacobian Vol; Integration I1; }
    }
  }
}
Resolution {
  { Name Check; System { { Name Sys; NameOfFormulation Maxwell; Type Complex; } }
    Operation { Generate[Sys]; Solve[Sys]; SaveSolution[Sys]; PostOperation[Errors]; }
  }
}
PostProcessing {
  { Name Fields; NameOfFormulation Maxwell; NameOfSystem Sys;
    PostQuantity {
      { Name EtError; Value { Integral {
        [SquNorm[-{d v}-Complex[0,omega[]]*{bt}-etExact[]]];
        In Domain; Jacobian Vol; Integration I1;
      } } }
      { Name EzError; Value { Integral {
        [SquNorm[-Complex[0,omega[]]*CompZ[{a}]+g2[]*{v}-ezExact[]]];
        In Domain; Jacobian Vol; Integration I1;
      } } }
      { Name BtError; Value { Integral {
        [SquNorm[{d a}-g2[]*(Vector[0,0,1] /\ {bt})-btExact[]]];
        In Domain; Jacobian Vol; Integration I1;
      } } }
      { Name BzError; Value { Integral {
        [SquNorm[CompZ[{d bt}]-curlB[]]];
        In Domain; Jacobian Vol; Integration I1;
      } } }
    }
  }
}
PostOperation {
  { Name Errors; NameOfPostProcessing Fields;
    Operation {
      Print[EtError[Domain], OnGlobal, Format Table, File "et.tsv"];
      Print[EzError[Domain], OnGlobal, Format Table, File "ez.tsv"];
      Print[BtError[Domain], OnGlobal, Format Table, File "bt.tsv"];
      Print[BzError[Domain], OnGlobal, Format Table, File "bz.tsv"];
    }
  }
}
