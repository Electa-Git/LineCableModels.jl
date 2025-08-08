// File created with GetDP.jl: https://github.com/Electa-Git/GetDP.jl.

DefineConstant[
active_con = {1, Choices{1,9999}, Name "Input/Active conductor", Visible 1}];

Group{
  DomainInf = Region[ {300100201,  300200102} ];// Domain transformation to infinity
  Con_1 = Region[ {100101103} ];// cable_1_core_con_copper
  Con_2 = Region[ {100102108} ];// cable_1_sheath_con_lead
  Con_3 = Region[ {100103110} ];// cable_1_armor_con_steel
  Con_4 = Region[ {100201103} ];// cable_2_core_con_copper
  Con_5 = Region[ {100202108} ];// cable_2_sheath_con_lead
  Con_6 = Region[ {100203110} ];// cable_2_armor_con_steel
  Conductors = Region[ {100101103,  100102108,  100103110,  100201103,  100202108,  100203110} ];
  DomainC = Region[ {} ];// All conductor materials
  DomainCC = Region[ {} ];// All non-conductor materials
  DomainActive = Region[ {Con~{active_con}} ];// Sources
  DomainInactive = Region[ {Conductors - Con~{active_con}} ];// Conductors set to zero energization
  DomainC += Region[ {100101103} ];// cable_1_core_con_copper
  DomainC += Region[ {100102108} ];// cable_1_sheath_con_lead
  DomainC += Region[ {100103110} ];// cable_1_armor_con_steel
  DomainC += Region[ {100201103} ];// cable_2_core_con_copper
  DomainC += Region[ {100202108} ];// cable_2_sheath_con_lead
  DomainC += Region[ {100203110} ];// cable_2_armor_con_steel
  DomainC += Region[ {200200102} ];// layer_2_earth_con_material_rho=100.0_epsr=10.0_mu=1.0
  DomainC += Region[ {300200102} ];// infshell_2_earth_con_material_rho=100.0_epsr=10.0_mu=1.0
  DomainCC += Region[ {100101201} ];// cable_1_core_ins_air
  DomainCC += Region[ {100101204} ];// cable_1_core_ins_semicon1
  DomainCC += Region[ {100101205} ];// cable_1_core_ins_pe
  DomainCC += Region[ {100101206} ];// cable_1_core_ins_semicon2
  DomainCC += Region[ {100101207} ];// cable_1_core_ins_polyacrylate
  DomainCC += Region[ {100102205} ];// cable_1_sheath_ins_pe
  DomainCC += Region[ {100102209} ];// cable_1_sheath_ins_pp
  DomainCC += Region[ {100103201} ];// cable_1_armor_ins_air
  DomainCC += Region[ {100103209} ];// cable_1_armor_ins_pp
  DomainCC += Region[ {100201201} ];// cable_2_core_ins_air
  DomainCC += Region[ {100201204} ];// cable_2_core_ins_semicon1
  DomainCC += Region[ {100201205} ];// cable_2_core_ins_pe
  DomainCC += Region[ {100201206} ];// cable_2_core_ins_semicon2
  DomainCC += Region[ {100201207} ];// cable_2_core_ins_polyacrylate
  DomainCC += Region[ {100202205} ];// cable_2_sheath_ins_pe
  DomainCC += Region[ {100202209} ];// cable_2_sheath_ins_pp
  DomainCC += Region[ {100203201} ];// cable_2_armor_ins_air
  DomainCC += Region[ {100203209} ];// cable_2_armor_ins_pp
  DomainCC += Region[ {200100201} ];// layer_1_air_ins_air
  DomainCC += Region[ {300100201} ];// infshell_1_air_ins_air
  Domain_Ele = Region[ {DomainCC,  DomainC} ];
  Sur_Dirichlet_Ele = Region[ {1002201,  1001201} ];

}
Function{
  // Material properties for region 100101207: cable_1_core_ins_polyacrylate
  nu[Region[{100101207}]] = 795774.7154594767;
  sigma[Region[{100101207}]] = 0.00018867924528301886;
  epsilon[Region[{100101207}]] = 2.8599026635344e-10;
  // Material properties for region 100201205: cable_2_core_ins_pe
  nu[Region[{100201205}]] = 795774.7154594767;
  sigma[Region[{100201205}]] = 5.0761421319796954e-15;
  epsilon[Region[{100201205}]] = 2.036463196944e-11;
  // Material properties for region 100201206: cable_2_core_ins_semicon2
  nu[Region[{100201206}]] = 795774.7154594767;
  sigma[Region[{100201206}]] = 0.002;
  epsilon[Region[{100201206}]] = 8.854187812800001e-9;
  // Material properties for region 100201207: cable_2_core_ins_polyacrylate
  nu[Region[{100201207}]] = 795774.7154594767;
  sigma[Region[{100201207}]] = 0.00018867924528301886;
  epsilon[Region[{100201207}]] = 2.8599026635344e-10;
  // Material properties for region 100202205: cable_2_sheath_ins_pe
  nu[Region[{100202205}]] = 795774.7154594767;
  sigma[Region[{100202205}]] = 5.0761421319796954e-15;
  epsilon[Region[{100202205}]] = 2.036463196944e-11;
  // Material properties for region 100203209: cable_2_armor_ins_pp
  nu[Region[{100203209}]] = 795774.7154594767;
  sigma[Region[{100203209}]] = 1.0e-15;
  epsilon[Region[{100203209}]] = 2.479172587584e-11;
  // Material properties for region 100202108: cable_2_sheath_con_lead
  nu[Region[{100202108}]] = 795788.2438596223;
  sigma[Region[{100202108}]] = 4.672897196261682e6;
  epsilon[Region[{100202108}]] = 8.8541878128e-12;
  // Material properties for region 100203201: cable_2_armor_ins_air
  nu[Region[{100203201}]] = 795774.7154594767;
  sigma[Region[{100203201}]] = 0.0;
  epsilon[Region[{100203201}]] = 8.8541878128e-12;
  // Material properties for region 100103110: cable_1_armor_con_steel
  nu[Region[{100103110}]] = 2652.5823848649225;
  sigma[Region[{100103110}]] = 7.246376811594203e6;
  epsilon[Region[{100103110}]] = 8.8541878128e-12;
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
  // Material properties for region 100202209: cable_2_sheath_ins_pp
  nu[Region[{100202209}]] = 795774.7154594767;
  sigma[Region[{100202209}]] = 1.0e-15;
  epsilon[Region[{100202209}]] = 2.479172587584e-11;
  // Material properties for region 100203110: cable_2_armor_con_steel
  nu[Region[{100203110}]] = 2652.5823848649225;
  sigma[Region[{100203110}]] = 7.246376811594203e6;
  epsilon[Region[{100203110}]] = 8.8541878128e-12;
  // Material properties for region 300100201: infshell_1_air_ins_air
  nu[Region[{300100201}]] = 795774.7154594767;
  sigma[Region[{300100201}]] = 0.0;
  epsilon[Region[{300100201}]] = 8.8541878128e-12;
  // Material properties for region 100101201: cable_1_core_ins_air
  nu[Region[{100101201}]] = 795774.7154594767;
  sigma[Region[{100101201}]] = 0.0;
  epsilon[Region[{100101201}]] = 8.8541878128e-12;
  // Material properties for region 100102209: cable_1_sheath_ins_pp
  nu[Region[{100102209}]] = 795774.7154594767;
  sigma[Region[{100102209}]] = 1.0e-15;
  epsilon[Region[{100102209}]] = 2.479172587584e-11;
  // Material properties for region 100102108: cable_1_sheath_con_lead
  nu[Region[{100102108}]] = 795788.2438596223;
  sigma[Region[{100102108}]] = 4.672897196261682e6;
  epsilon[Region[{100102108}]] = 8.8541878128e-12;
  // Material properties for region 100101204: cable_1_core_ins_semicon1
  nu[Region[{100101204}]] = 795774.7154594767;
  sigma[Region[{100101204}]] = 0.001;
  epsilon[Region[{100101204}]] = 8.854187812800001e-9;
  // Material properties for region 100103201: cable_1_armor_ins_air
  nu[Region[{100103201}]] = 795774.7154594767;
  sigma[Region[{100103201}]] = 0.0;
  epsilon[Region[{100103201}]] = 8.8541878128e-12;
  // Material properties for region 100103209: cable_1_armor_ins_pp
  nu[Region[{100103209}]] = 795774.7154594767;
  sigma[Region[{100103209}]] = 1.0e-15;
  epsilon[Region[{100103209}]] = 2.479172587584e-11;
  // Material properties for region 100201103: cable_2_core_con_copper
  nu[Region[{100201103}]] = 795779.4901364174;
  sigma[Region[{100201103}]] = 5.8001276028072625e7;
  epsilon[Region[{100201103}]] = 8.8541878128e-12;
  // Material properties for region 100102205: cable_1_sheath_ins_pe
  nu[Region[{100102205}]] = 795774.7154594767;
  sigma[Region[{100102205}]] = 5.0761421319796954e-15;
  epsilon[Region[{100102205}]] = 2.036463196944e-11;
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
  // Material properties for region 100201201: cable_2_core_ins_air
  nu[Region[{100201201}]] = 795774.7154594767;
  sigma[Region[{100201201}]] = 0.0;
  epsilon[Region[{100201201}]] = 8.8541878128e-12;

}
Function{
  Freq = 0.001;
  UnitAmplitude = 1.0;

}
Constraint{
  { Name ScalarPotential_2D; Type Assign;
    Case {
      { Region DomainInactive; Value 0.0; }
      { Region Con~{active_con}; Value UnitAmplitude; }
      { Region Sur_Dirichlet_Ele; Value 0.0; }
    }
  }
  { Name Charge_2D; Type Assign;
    Case {
    }
  }
}
FunctionSpace {
  { Name Hgrad_v_Ele; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node; Support Domain_Ele; Entity NodesOf[ All, Not Conductors ]; }
      { Name sf; NameOfCoef vf; Function BF_GroupOfNodes; Support Domain_Ele; Entity GroupsOfNodesOf[ Conductors ]; }
    }
    GlobalQuantity {
      { Name U; Type AliasOf; NameOfCoef vf;}
      { Name Q; Type AssociatedWith; NameOfCoef vf;}
    }
    Constraint {
      { NameOfCoef U; EntityType Region; NameOfConstraint ScalarPotential_2D;}
      { NameOfCoef Q; EntityType Region; NameOfConstraint Charge_2D;}
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint ScalarPotential_2D;}
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
  { Name Electrodynamics_v; Type FemEquation;
      Quantity {
        { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
        { Name U; Type Global; NameOfSpace Hgrad_v_Ele [U]; }
        { Name Q; Type Global; NameOfSpace Hgrad_v_Ele [Q]; }
      }
      Equation {
        Galerkin { [ sigma[] * Dof{d v} , {d v} ];
            In Domain_Ele;
            Jacobian Vol;
            Integration I1;
        }
        Galerkin { DtDof[ epsilon[] * Dof{d v} , {d v} ];
            In DomainCC;
            Jacobian Vol;
            Integration I1;
        }
        GlobalTerm {  [ Dof{Q} , {U} ] ;
            In Conductors;
        }
      }
  }
}
Resolution {
  { Name Electrodynamics;
    System {
      { Name Sys_Ele; NameOfFormulation Electrodynamics_v; Type Complex; Frequency Freq; }
    }
    Operation {
      CreateDir["results/electrodynamics"];
      Generate[Sys_Ele];
      Solve[Sys_Ele];
      SaveSolution[Sys_Ele];
      PostOperation[LineParams];
    }
  }
}

PostProcessing {
  { Name EleDyn_v; NameOfFormulation Electrodynamics_v;
      PostQuantity {
      { Name v; Value {
          Term { [ {v} ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name e; Value {
          Term { [ -{d v} ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name em; Value {
          Term { [ Norm[-{d v}] ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name d; Value {
          Term { [ -epsilon[] * {d v} ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name dm; Value {
          Term { [ Norm[-epsilon[] * {d v}] ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name j; Value {
          Term { [ -sigma[] * {d v} ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name jm; Value {
          Term { [ Norm[-sigma[] * {d v}] ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name jtot; Value {
          Term { Type Global; [ -sigma[] * {d v} ];
              In Domain_Ele;
              Jacobian Vol;
          }
          Term { Type Global; [ -epsilon[] * Dt[{d v}] ];
              In Domain_Ele;
              Jacobian Vol;
          }
      }}
      { Name U; Value {
          Term { [ {U} ];
              In Domain_Ele;
          }
      }}
      { Name Q; Value {
          Term { [ {Q} ];
              In Domain_Ele;
          }
      }}
      { Name Y; Value {
          Term { [ -{Q} ];
              In Domain_Ele;
          }
      }}
      }
  }
}

PostOperation {
  { Name Field_Maps; NameOfPostProcessing EleDyn_v;
      Operation {
          Print[ v, OnElementsOf Domain_Ele, File StrCat[ "results/electrodynamics/v_", Sprintf("%g",active_con), ".pos" ] ];
          Print[ em, OnElementsOf Domain_Ele, Name "|E| [V/m]", File StrCat[ "results/electrodynamics/em_", Sprintf("%g",active_con), ".pos" ] ];
          Print[ dm, OnElementsOf Domain_Ele, Name "|D| [A/m²]", File StrCat[ "results/electrodynamics/dm_", Sprintf("%g",active_con), ".pos" ] ];
          Print[ e, OnElementsOf Domain_Ele, Name "E [V/m]", File StrCat[ "results/electrodynamics/e_", Sprintf("%g",active_con), ".pos" ] ];
      }
  }
  { Name LineParams; NameOfPostProcessing EleDyn_v;
      Operation {
          Print[ Y, OnRegion Conductors, Format Table, File "results/electrodynamics/Y.dat", AppendToExistingFile (active_con > 1 ? 1 : 0) ];
      }
  }
}
