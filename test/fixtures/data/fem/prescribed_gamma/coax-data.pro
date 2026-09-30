// Finite circular boundary with zero potentials; exact cylindrical reference.
NumTerminals=1; NumCables=1; NumMaterialRegions=1; FrequencyCount=1;
AIR_EM=1001; EARTH_EM=1002; AIR_PML=1003; EARTH_PML=1004;
OUTBND_EM=2001; TERMINAL=3001; TERMINAL_CONTOUR=4001;
VOLTAGE_PATH=7001; VOLTAGE_REFERENCE=8001;
ReceiverInAir()={0}; TerminalNames()=Str["inner"];
MaterialRegionTags()={9001}; MaterialIsConductor()={1}; MaterialHasLoss()={1};
MaterialMu()={1.25}; MaterialSigma_1()={1.}; MaterialEpsilon_1()={0.};
EarthMu()={1.25}; EarthSigma()={.7}; EarthEpsilon()={.4};
AirMu=1.25; AirEpsilon=.4;
Xcenter=0.; Ycenter=0.; DomainHalfwidth=2.;
PmlSideThickness=1.; PmlTopThickness=1.; PmlBottomThickness=1.;
PmlSideStrength=0.; PmlTopStrength=0.; PmlBottomStrength=0.;
PerfectConductors=1; FrequencyHz=1.7/(2*Pi); VolumeQuadrature=12;
