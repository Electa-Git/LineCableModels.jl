using LineCableModels

# The caller supplies the filename; the backend does not search for one.
config_path = isempty(ARGS) ? ENV["LINECABLEMODELS_PSCAD_CONFIG"] : only(ARGS)
station = PSCAD.RemoteConfig(config_path)

copper = Material(:conductor, 1.72e-8, 1.0)
dielectric = Material(:insulator, 1e14, 2.3)
design = build(CableDesign, "insulated-wire", terminal(:core,
    solid(copper, Disk(0.004)), insulation(dielectric; t=0.002)))
system = build(LineCableSystem, [design, design], [Pose2(0, -1), Pose2(1, -1)];
    connections=[Dict(:core=>1), Dict(:core=>2)])
problem = LineParametersProblem(system; earth_props=homogeneous(rho=100.0),
    frequencies=10.0 .^ range(-1, 6; length=101))
selected = Formulation(:pscad; earth_impedance=:wedepohl1973)
@time result = compute(problem, selected;
    options=(remote=station, verbosity=(default=0, PSCAD=1)))
Zresult = Z(result)
Yresult = Y(result)
