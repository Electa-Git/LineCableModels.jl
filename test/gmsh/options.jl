using Gmsh
using getdp_jll
using LineCableModels.Commons: ComputationOptions, defaults
using LineCableModels.Engine: GmshBackend

@testset "Gmsh computation options" begin
    options = ComputationOptions(GmshBackend; mesh_size_factor = 2.0)
    @test getfield(options, :values) ===
          merge(defaults(GmshBackend), (mesh_size_factor = 2.0,))
    @test_throws ArgumentError ComputationOptions(GmshBackend; mesh_policy = :reuse)

    native = (
        PmlSideThicknessFactor = 1.0,
        PmlTopThicknessFactor = 1.0,
        PmlBottomThicknessFactor = 1.0,
        PmlSideLayers = 16,
        PmlTopLayers = 16,
        PmlBottomLayers = 16,
        PmlSideGrading = log(20),
        PmlTopGrading = log(20),
        PmlBottomGrading = log(20),
        ExteriorMeshSizeFactor = 1.0,
        InterfaceRefinementFactor = 1.0,
        ConductorGeometryTolerance = 1e-3,
        ConductorSkinDepthElements = 6.0,
        ConductorMeshGrowth = sqrt(1.25),
        ConductorSkinDepths = 5.0,
        ConductorThicknessElements = 4,
        VolumeQuadrature = 12,
        PhysicalVolumeQuadrature = 3,
        PmlQuadrangles = 1,
        PmlQuadrature = 4,
        MumpsOrdering = 0,
        PetscPrealloc = 0,
        MumpsErrorAnalysis = 2,
        MumpsRefinementMax = 2,
        MumpsBackwardErrorTolerance = 1e-12,
        ThinLayerStrips = 1)
    @test defaults(GmshBackend).overrides == native

    options = ComputationOptions(GmshBackend; overrides = (PmlSideLayers = 48,))
    @test options.overrides === merge(native, (PmlSideLayers = 48,))
    @test ComputationOptions(GmshBackend;
        overrides = (PmlSideLayers = 48.0,)).overrides.PmlSideLayers === 48

    for (overrides, name) in (((Foo = 1,), "Foo"),
        ((PmlSideLayers = 0,), "PmlSideLayers"),
        ((ConductorMeshGrowth = 0.5,), "ConductorMeshGrowth"),
        ((ConductorSkinDepths = Inf,), "ConductorSkinDepths"),
        ((VolumeQuadrature = 5,), "VolumeQuadrature"))
        @test_throws ArgumentError ComputationOptions(GmshBackend; overrides)
        @test_throws name ComputationOptions(GmshBackend; overrides)
    end
end

# With the extension loaded, the generic computation builds the options and refuses a
# formulation that the Gmsh backend does not compute.
@testset "the Gmsh backend refuses the coaxial formulation" begin
    problem = Dict(REFERENCES)["two_bare_wires"]
    @test_throws ArgumentError("GmshBackend does not implement CoaxialFormulation") compute(
        GmshBackend(), problem, CoaxialFormulation())
end
