// File created with GetDP.jl: https://github.com/Electa-Git/GetDP.jl.

DefineConstant[
active_con = {1, Choices{1,9999}, Name "Input/Active conductor", Visible 1}];

Group{
  DomainInf = Region[ {300100201,  300200102} ];// Domain transformation to infinity
  Con_1 = Region[ {100101103} ];// cable_1_core_con_copper
  Con_2 = Region[ {100201103} ];// cable_2_core_con_copper
  Conductors = Region[ {100101103,  100201103} ];
  DomainC = Region[ {} ];// All conductor materials
  DomainCC = Region[ {} ];// All non-conductor materials
  DomainActive = Region[ {Con~{active_con}} ];// Sources
  DomainInactive = Region[ {Conductors - Con~{active_con}} ];// Conductors set to zero energization
  DomainC += Region[ {100101103} ];// cable_1_core_con_copper
  DomainC += Region[ {100201103} ];// cable_2_core_con_copper
  DomainC += Region[ {200200102} ];// layer_2_earth_con_material_rho=100.0_epsr=10.0_mu=1.0
  DomainC += Region[ {300200102} ];// infshell_2_earth_con_material_rho=100.0_epsr=10.0_mu=1.0
  DomainCC += Region[ {100101201} ];// cable_1_core_ins_air
  DomainCC += Region[ {100101204} ];// cable_1_core_ins_semicon1
  DomainCC += Region[ {100101205} ];// cable_1_core_ins_pe
  DomainCC += Region[ {100101206} ];// cable_1_core_ins_semicon2
  DomainCC += Region[ {100101207} ];// cable_1_core_ins_polyacrylate
  DomainCC += Region[ {100201201} ];// cable_2_core_ins_air
  DomainCC += Region[ {100201204} ];// cable_2_core_ins_semicon1
  DomainCC += Region[ {100201205} ];// cable_2_core_ins_pe
  DomainCC += Region[ {100201206} ];// cable_2_core_ins_semicon2
  DomainCC += Region[ {100201207} ];// cable_2_core_ins_polyacrylate
  DomainCC += Region[ {200100201} ];// layer_1_air_ins_air
  DomainCC += Region[ {300100201} ];// infshell_1_air_ins_air
  Domain_Mag = Region[ {DomainCC,  DomainC} ];
  Sur_Dirichlet_Mag = Region[ {1002201,  1001201} ];

}
Function{
  // Material properties for region 100101207: cable_1_core_ins_polyacrylate
  nu[Region[{100101207}]] = 795774.7154594767;
  sigma[Region[{100101207}]] = 0.00018867924528301886;
  epsilon[Region[{100101207}]] = 2.8599026635344e-10;
  // Material properties for region 100201206: cable_2_core_ins_semicon2
  nu[Region[{100201206}]] = 795774.7154594767;
  sigma[Region[{100201206}]] = 0.002;
  epsilon[Region[{100201206}]] = 8.854187812800001e-9;
  // Material properties for region 100101103: cable_1_core_con_copper
  nu[Region[{100101103}]] = 795779.4901364174;
  sigma[Region[{100101103}]] = 5.8001276028072625e7;
  epsilon[Region[{100101103}]] = 8.8541878128e-12;
  /*  Material properties for region 200200102:
  layer_2_earth_con_material_rho=100.0_epsr=10.0_mu=1.0 */
  nu[Region[{200200102}]] = 795774.7154594767;
  sigma[Region[{200200102}]] = 0.01;
  epsilon[Region[{200200102}]] = 8.8541878128e-11;
  // Material properties for region 100101205: cable_1_core_ins_pe
  nu[Region[{100101205}]] = 795774.7154594767;
  sigma[Region[{100101205}]] = 5.0761421319796954e-15;
  epsilon[Region[{100101205}]] = 2.036463196944e-11;
  // Material properties for region 100201204: cable_2_core_ins_semicon1
  nu[Region[{100201204}]] = 795774.7154594767;
  sigma[Region[{100201204}]] = 0.001;
  epsilon[Region[{100201204}]] = 8.854187812800001e-9;
  // Material properties for region 300100201: infshell_1_air_ins_air
  nu[Region[{300100201}]] = 795774.7154594767;
  sigma[Region[{300100201}]] = 0.0;
  epsilon[Region[{300100201}]] = 8.8541878128e-12;
  // Material properties for region 100101201: cable_1_core_ins_air
  nu[Region[{100101201}]] = 795774.7154594767;
  sigma[Region[{100101201}]] = 0.0;
  epsilon[Region[{100101201}]] = 8.8541878128e-12;
  // Material properties for region 100101204: cable_1_core_ins_semicon1
  nu[Region[{100101204}]] = 795774.7154594767;
  sigma[Region[{100101204}]] = 0.001;
  epsilon[Region[{100101204}]] = 8.854187812800001e-9;
  // Material properties for region 100201103: cable_2_core_con_copper
  nu[Region[{100201103}]] = 795779.4901364174;
  sigma[Region[{100201103}]] = 5.8001276028072625e7;
  epsilon[Region[{100201103}]] = 8.8541878128e-12;
  // Material properties for region 100201201: cable_2_core_ins_air
  nu[Region[{100201201}]] = 795774.7154594767;
  sigma[Region[{100201201}]] = 0.0;
  epsilon[Region[{100201201}]] = 8.8541878128e-12;
  // Material properties for region 100201205: cable_2_core_ins_pe
  nu[Region[{100201205}]] = 795774.7154594767;
  sigma[Region[{100201205}]] = 5.0761421319796954e-15;
  epsilon[Region[{100201205}]] = 2.036463196944e-11;
  // Material properties for region 200100201: layer_1_air_ins_air
  nu[Region[{200100201}]] = 795774.7154594767;
  sigma[Region[{200100201}]] = 0.0;
  epsilon[Region[{200100201}]] = 8.8541878128e-12;
  // Material properties for region 100101206: cable_1_core_ins_semicon2
  nu[Region[{100101206}]] = 795774.7154594767;
  sigma[Region[{100101206}]] = 0.002;
  epsilon[Region[{100101206}]] = 8.854187812800001e-9;
  /*  Material properties for region 300200102:
  infshell_2_earth_con_material_rho=100.0_epsr=10.0_mu=1.0 */
  nu[Region[{300200102}]] = 795774.7154594767;
  sigma[Region[{300200102}]] = 0.01;
  epsilon[Region[{300200102}]] = 8.8541878128e-11;
  // Material properties for region 100201207: cable_2_core_ins_polyacrylate
  nu[Region[{100201207}]] = 795774.7154594767;
  sigma[Region[{100201207}]] = 0.00018867924528301886;
  epsilon[Region[{100201207}]] = 2.8599026635344e-10;

}
Function{
  Freq = 0.001;
  UnitAmplitude = 1.0;

}
Constraint{
  { Name MagneticVectorPotential_2D; Type Assign;
    Case {
      { Region Sur_Dirichlet_Mag; Value 0.0; }
    }
  }
  { Name Voltage_2D; Type Assign;
    Case {
    }
  }
  { Name Current_2D; Type Assign;
    Case {
      { Region DomainInactive; Value 0.0; }
      { Region Con~{active_con}; Value UnitAmplitude; }
    }
  }
}
FunctionSpace {
  { Name Hcurl_a_Mag_2D; Type Form1P;
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_PerpendicularEdge; Support Domain_Mag; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef ae; EntityType NodesOf; NameOfConstraint MagneticVectorPotential_2D;}
    }
  }
  { Name Hregion_u_Mag_2D; Type Form1P;
    BasisFunction {
      { Name sr; NameOfCoef ur; Function BF_RegionZ; Support DomainC; Entity DomainC; }
    }
    GlobalQuantity {
      { Name U; Type AliasOf; NameOfCoef ur;}
      { Name I; Type AssociatedWith; NameOfCoef ur;}
    }
    Constraint {
      { NameOfCoef U; EntityType Region; NameOfConstraint Voltage_2D;}
      { NameOfCoef I; EntityType Region; NameOfConstraint Current_2D;}
    }
  }
}

Jacobian{
  { Name Vol; Case {
   {  Region DomainInf; Jacobian VolSphShell{5000.0, 6250.0, 0.0, 0.0, 0.0}; } 
   {  Region All; Jacobian Vol; } 
  } }
  { Name Sur; Case {
   {  Region All; Jacobian Sur; } 
  } }
}
Integration {
  { Name I1;
    Case {
      { Type Gauss;
      Case {
         {  GeoElement Point; NumberOfPoints 1; } 
         {  GeoElement Line; NumberOfPoints 4; } 
         {  GeoElement Triangle; NumberOfPoints 4; } 
         {  GeoElement Quadrangle; NumberOfPoints 4; } 
      }
      }
    }
  }
}
Formulation {
  { Name Darwin_a_2D; Type FemEquation;
      Quantity {
        { Name a; Type Local; NameOfSpace Hcurl_a_Mag_2D; }
        { Name ur; Type Local; NameOfSpace Hregion_u_Mag_2D; }
        { Name I; Type Global; NameOfSpace Hregion_u_Mag_2D [I]; }
        { Name U; Type Global; NameOfSpace Hregion_u_Mag_2D [U]; }
      }
      Equation {
        Galerkin { [ nu[] * Dof{d a} , {d a} ];
            In Domain_Mag;
            Jacobian Vol;
            Integration I1;
        }
        Galerkin { DtDof [ sigma[] * Dof{a} , {a} ];
            In DomainC;
            Jacobian Vol;
            Integration I1;
        }
        Galerkin { [ sigma[] * Dof{ur}, {a} ];
            In DomainC;
            Jacobian Vol;
            Integration I1;
        }
        Galerkin { DtDof [ sigma[] * Dof{a} , {ur} ];
            In DomainC;
            Jacobian Vol;
            Integration I1;
        }
        Galerkin { [ sigma[] * Dof{ur}, {ur}];
            In DomainC;
            Jacobian Vol;
            Integration I1;
        }
        //  Darwin approximation term
        Galerkin { DtDtDof [ epsilon[] * Dof{a} , {a}];
            In DomainC;
            Jacobian Vol;
            Integration I1;
        }
        Galerkin { DtDof[ epsilon[] * Dof{ur}, {a} ];
            In DomainC;
            Jacobian Vol;
            Integration I1;
        }
        Galerkin { DtDtDof [ epsilon[] * Dof{a} , {ur}];
            In DomainC;
            Jacobian Vol;
            Integration I1;
        }
        Galerkin { DtDof[ epsilon[] * Dof{ur}, {ur} ];
            In DomainC;
            Jacobian Vol;
            Integration I1;
        }
        GlobalTerm {  [ Dof{I} , {U} ] ;
            In Conductors;
        }
      }
  }
}
Resolution {
  { Name Darwin;
    System {
      { Name Sys_Mag; NameOfFormulation Darwin_a_2D; Type Complex; Frequency Freq; }
    }
    Operation {
      CreateDir["results/darwin"];
      InitSolution[Sys_Mag];
      Generate[Sys_Mag];
      Solve[Sys_Mag];
      SaveSolution[Sys_Mag];
      PostOperation[LineParams];
    }
  }
}

PostProcessing {
  { Name Darwin_a_2D; NameOfFormulation Darwin_a_2D;
      PostQuantity {
      { Name a; Value {
          Term { [ {a} ];
              In Domain_Mag;
              Jacobian Vol;
          }
      }}
      { Name az; Value {
          Term { [ CompZ[{a}] ];
              In Domain_Mag;
              Jacobian Vol;
          }
      }}
      { Name b; Value {
          Term { [ {d a} ];
              In Domain_Mag;
              Jacobian Vol;
          }
      }}
      { Name bm; Value {
          Term { [ Norm[{d a}] ];
              In Domain_Mag;
              Jacobian Vol;
          }
      }}
      { Name j; Value {
          Term { [ -sigma[]*(Dt[{a}]+{ur}) ];
              In DomainC;
              Jacobian Vol;
          }
      }}
      { Name jz; Value {
          Term { [ CompZ[-sigma[]*(Dt[{a}]+{ur})] ];
              In DomainC;
              Jacobian Vol;
          }
      }}
      { Name jm; Value {
          Term { [ Norm[-sigma[]*(Dt[{a}]+{ur})] ];
              In DomainC;
              Jacobian Vol;
          }
      }}
      { Name d; Value {
          Term { [ epsilon[] * Dt[Dt[{a}]+{ur}] ];
              In DomainC;
              Jacobian Vol;
          }
      }}
      { Name dz; Value {
          Term { [ CompZ[epsilon[] * Dt[Dt[{a}]+{ur}]] ];
              In DomainC;
              Jacobian Vol;
          }
      }}
      { Name dm; Value {
          Term { [ Norm[epsilon[] * Dt[Dt[{a}]+{ur}]] ];
              In DomainC;
              Jacobian Vol;
          }
      }}
      { Name rhoj2; Value {
          Term { [ 0.5*sigma[]*SquNorm[Dt[{a}]+{ur}] ];
              In DomainC;
              Jacobian Vol;
          }
      }}
      { Name U; Value {
          Term { [ {U} ];
              In DomainC;
          }
      }}
      { Name I; Value {
          Term { [ {I} ];
              In DomainC;
          }
      }}
      { Name Z; Value {
          Term { [ -{U} ];
              In DomainC;
          }
      }}
      }
  }
}

PostOperation {
  { Name Field_Maps; NameOfPostProcessing Darwin_a_2D;
      Operation {
          Print[ az, OnElementsOf Domain_Mag, Smoothing 1, Name "flux lines: Az [T m]", File StrCat[ "results/darwin/az_", Sprintf("%g",active_con), ".pos" ] ];
          Print[ b, OnElementsOf Domain_Mag, Smoothing 1, Name "B [T]", File StrCat[ "results/darwin/b_", Sprintf("%g",active_con), ".pos" ] ];
          Print[ bm, OnElementsOf Domain_Mag, Smoothing 1, Name "|B| [T]", File StrCat[ "results/darwin/bm_", Sprintf("%g",active_con), ".pos" ] ];
          Print[ jz, OnElementsOf Region[{DomainC}], Smoothing 1, Name "jz [A/m²] Conducting domain", File StrCat[ "results/darwin/jz_", Sprintf("%g",active_con), ".pos" ] ];
          Print[ rhoj2, OnElementsOf Region[{DomainC}], Smoothing 1, Name "Power density", File StrCat[ "results/darwin/rhoj2_", Sprintf("%g",active_con), ".pos" ] ];
          Print[ jm, OnElementsOf DomainC, Smoothing 1, Name "|j| [A/m²] Conducting domain", File StrCat[ "results/darwin/jm_", Sprintf("%g",active_con), ".pos" ] ];
          Print[ dm, OnElementsOf DomainC, Smoothing 1, Name "|D| [A/m²]", File StrCat[ "results/darwin/dm_", Sprintf("%g",active_con), ".pos" ] ];
      }
  }
  { Name LineParams; NameOfPostProcessing Darwin_a_2D;
      Operation {
          Print[ Z, OnRegion Conductors, Format Table, File "results/darwin/Z.dat", AppendToExistingFile (active_con > 1 ? 1 : 0) ];
      }
  }
}
