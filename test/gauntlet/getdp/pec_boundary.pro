// Exact PEC boundary limit of the existing quasi-TEM A_z / phi equations.
// Conductor interiors are excluded from Domain_Mag. Each terminal has constant
// traces a=A_z and psi=phi/Gamma; their boundary reactions prescribe I and I.
// At Gamma=0, the normalized exterior weak equations are
//   (nu grad a, grad v) + (j omega se a,v) = sum I_p v_p,
//   (se grad psi,grad w) + (j omega se a,w) = sum I_p w_p.
// Z=j omega a_p/I and raw P=psi_p/I, so Y=inv(raw P).
// No conductor conductivity, internal impedance, or analytical earth kernel
// occurs in this operator. The electric reference is the existing earth boundary.
// This is not the original thin_wire.pro Electric Helmholtz formulation.
// It must not be used as a reproduction of that historical FEM test; see
// docs/notes/thin-wire-fem-provenance.md for the original operator and replay.

FunctionSpace {
  { Name Hgrad_a_PEC; Type Form0;
    BasisFunction {
      { Name an; NameOfCoef an; Function BF_Node;
        Support Domain_Mag; Entity NodesOf[All, Not Terminals]; }
      { Name af; NameOfCoef af; Function BF_GroupOfNodes;
        Support Domain_Mag; Entity GroupsOfNodesOf[Terminals]; }
    }
    GlobalQuantity {
      { Name A; Type AliasOf; NameOfCoef af; }
      { Name I; Type AssociatedWith; NameOfCoef af; }
    }
    Constraint {
      { NameOfCoef an; EntityType NodesOf;
        NameOfConstraint FEMMagneticVectorPotential; }
      { NameOfCoef I; EntityType Auto; NameOfConstraint FEMTerminalCurrent; }
    }
  }
  { Name Hgrad_psi_PEC; Type Form0;
    BasisFunction {
      { Name pn; NameOfCoef pn; Function BF_Node;
        Support Domain_Mag; Entity NodesOf[All, Not Terminals]; }
      { Name pf; NameOfCoef pf; Function BF_GroupOfNodes;
        Support Domain_Mag; Entity GroupsOfNodesOf[Terminals]; }
    }
    GlobalQuantity {
      { Name Psi; Type AliasOf; NameOfCoef pf; }
      { Name J; Type AssociatedWith; NameOfCoef pf; }
    }
    Constraint {
      { NameOfCoef pn; EntityType NodesOf; NameOfConstraint FEMScalarPotential; }
      { NameOfCoef J; EntityType Auto; NameOfConstraint FEMTerminalCurrent; }
    }
  }
}

Formulation {
  { Name FEM_a_phi_2D; Type FemEquation;
    Quantity {
      { Name a; Type Local; NameOfSpace Hgrad_a_PEC; }
      { Name A; Type Global; NameOfSpace Hgrad_a_PEC [A]; }
      { Name I; Type Global; NameOfSpace Hgrad_a_PEC [I]; }
      { Name psi; Type Local; NameOfSpace Hgrad_psi_PEC; }
      { Name Psi; Type Global; NameOfSpace Hgrad_psi_PEC [Psi]; }
      { Name J; Type Global; NameOfSpace Hgrad_psi_PEC [J]; }
    }
    Equation {
      Galerkin { [nu[] * Dof{d a}, {d a}];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0.,omega[]] * se[] * Dof{a}, {a}];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      GlobalTerm { [-Dof{I}, {A}]; In Terminals; }

      Galerkin { [se[] * Dof{d psi}, {d psi}];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0.,omega[]] * se[] * Dof{a}, {psi}];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      GlobalTerm { [-Dof{J}, {Psi}]; In Terminals; }
    }
  }
}
