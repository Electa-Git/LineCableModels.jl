// Bind resolved numerical inputs to their primary mesh material regions.
Group {
  ConductorMaterialRegions = Region[{}];
  PassiveMaterialRegions = Region[{}];
  LossyMaterialRegions = Region[{}];
  For material In {1:NumMaterialRegions}
    // Parentheses force a scalar tag: without them GetDP expands the list.
    MaterialRegion~{material} = Region[{(MaterialRegionTags(material - 1))}];
    If(MaterialIsConductor(material - 1))
      ConductorMaterialRegions += Region[{MaterialRegion~{material}}];
    Else
      PassiveMaterialRegions += Region[{MaterialRegion~{material}}];
      If(MaterialHasLoss(material - 1))
        LossyMaterialRegions += Region[{MaterialRegion~{material}}];
      EndIf
    EndIf
  EndFor
  DomainCWithI = Region[{Terminals}];
  DomainC = Region[{ConductorMaterialRegions, Earth, EarthInf}];
  DomainCC = Region[{Air, AirInf, PassiveMaterialRegions}];
  DomainLoss = Region[{DomainC, LossyMaterialRegions}];
  Domain_Mag = Region[{DomainC, DomainCC}];
}

Function {
  For material In {1:NumMaterialRegions}
    nu[MaterialRegion~{material}] = 1. / MaterialMu(material - 1);
    mu[MaterialRegion~{material}] = MaterialMu(material - 1);
    // Sigma is effective conductivity: selected dielectric losses are already
    // included in Re(kappa). Do not add a second tan(delta) contribution.
    sigma[MaterialRegion~{material}] = MaterialSigma~{material}(FrequencyIndex - 1);
    epsilon[MaterialRegion~{material}] = MaterialEpsilon~{material}(FrequencyIndex - 1);
  EndFor
}
