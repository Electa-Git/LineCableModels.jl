using Test
using LineCableModels
using DataFrames
using Measurements

# Access internal components for testing
const LCM = LineCableModels
const EP = LCM.EarthProps



@testset "FDEM Formulations" begin
    @testset "CPEarth" begin
        cp_formulation = EP.CPEarth()
        @test cp_formulation isa EP.AbstractFDEMFormulation
        @test EP._get_description(cp_formulation) == "CP model"
    end
end

@testset "_calc_earth_properties" begin
    frequencies = [50.0, 60.0, 1000.0]
    base_rho_g = 100.0
    base_epsr_g = 10.0
    base_mur_g = 1.0
    formulation = EP.CPEarth()

    @testset "Float64 inputs" begin
        rho, epsilon, mu = EP._calc_earth_properties(frequencies, base_rho_g, base_epsr_g, base_mur_g, formulation)

        @test length(rho) == length(frequencies)
        @test all(r -> r == base_rho_g, rho)

        @test length(epsilon) == length(frequencies)
        @test all(e -> isapprox(e, LCM.ε₀ * base_epsr_g), epsilon)

        @test length(mu) == length(frequencies)
        @test all(m -> isapprox(m, LCM.μ₀ * base_mur_g), mu)
    end

    @testset "Measurement inputs" begin
        rho_m = 100.0 ± 5.0
        epsr_m = 10.0 ± 0.5
        mur_m = 1.0 ± 0.01

        rho, epsilon, mu = EP._calc_earth_properties(frequencies, rho_m, epsr_m, mur_m, formulation)

        @test length(rho) == length(frequencies)
        @test all(r -> r == rho_m, rho)

        @test length(epsilon) == length(frequencies)
        @test all(e -> e.val ≈ (LCM.ε₀ * epsr_m).val, epsilon)
        @test all(e -> e.err ≈ (LCM.ε₀ * epsr_m).err, epsilon)

        @test length(mu) == length(frequencies)
        @test all(m -> m.val ≈ (LCM.μ₀ * mur_m).val, mu)
        @test all(m -> m.err ≈ (LCM.μ₀ * mur_m).err, mu)
    end
end

@testset "EarthLayer Constructor" begin
    frequencies = [50.0, 60.0]
    base_rho_g = 100.0
    base_epsr_g = 10.0
    base_mur_g = 1.0
    t = 5.0
    formulation = EP.CPEarth()

    layer = EP.EarthLayer(frequencies, base_rho_g, base_epsr_g, base_mur_g, t, formulation)

    @test layer.base_rho_g == base_rho_g
    @test layer.base_epsr_g == base_epsr_g
    @test layer.base_mur_g == base_mur_g
    @test layer.t == t
    @test length(layer.rho_g) == length(frequencies)
    @test all(layer.rho_g .== base_rho_g)
end

@testset "EarthModel Constructor" begin
    frequencies = [50.0, 60.0]
    rho_g = 100.0
    epsr_g = 10.0
    mur_g = 1.0

    @testset "Homogeneous Model" begin
        model = EarthModel(frequencies, rho_g, epsr_g, mur_g)
        @test length(model.layers) == 2 # Air + 1 earth layer
        @test model.vertical_layers == false
        @test model.FDformulation isa EP.CPEarth
        @test isinf(model.layers[1].t) # Air layer
        @test isinf(model.layers[2].t) # Homogeneous earth
        @test model.layers[2].base_rho_g == rho_g
    end

    @testset "Finite Thickness Layer" begin
        model = EarthModel(frequencies, rho_g, epsr_g, mur_g, t=20.0)
        @test length(model.layers) == 2
        @test model.layers[2].t == 20.0
    end

    @testset "Vertical Layers" begin
        model = EarthModel(frequencies, rho_g, epsr_g, mur_g, vertical_layers=true)
        @test model.vertical_layers == true
    end

    @testset "Input Validation" begin
        @test_throws AssertionError EarthModel([-50.0], rho_g, epsr_g, mur_g)
        @test_throws AssertionError EarthModel(frequencies, -100.0, epsr_g, mur_g)
        @test_throws AssertionError EarthModel(frequencies, rho_g, -10.0, mur_g)
        @test_throws AssertionError EarthModel(frequencies, rho_g, epsr_g, -1.0)
        @test_throws AssertionError EarthModel(frequencies, rho_g, epsr_g, mur_g, t=-5.0)
    end
end

@testset "add! for EarthModel" begin
    frequencies = [50.0, 60.0]

    @testset "Horizontal Layering" begin
        model = EarthModel(frequencies, 100.0, 10.0, 1.0, t=20.0)
        @test length(model.layers) == 2

        add!(model, frequencies, 200.0, 15.0, 1.0, t=50.0)
        @test length(model.layers) == 3
        @test model.layers[3].base_rho_g == 200.0
        @test model.layers[3].t == 50.0

        add!(model, frequencies, 500.0, 20.0, 1.0, t=Inf)
        @test length(model.layers) == 4
        @test isinf(model.layers[4].t)

        # Test invalid addition
        @test_throws ErrorException add!(model, frequencies, 1000.0, 25.0, 1.0, t=Inf)
    end

    @testset "Input Validation in add!" begin
        model = EarthModel(frequencies, 100.0, 10.0, 1.0, t=20.0)
        @test_throws AssertionError add!(model, [-50.0], 200.0, 15.0, 1.0)
        @test_throws AssertionError add!(model, frequencies, -200.0, 15.0, 1.0)
        @test_throws AssertionError add!(model, frequencies, 200.0, -15.0, 1.0)
        @test_throws AssertionError add!(model, frequencies, 200.0, 15.0, -1.0)
        @test_throws AssertionError add!(model, frequencies, 200.0, 15.0, 1.0, t=-5.0)
    end
end

@testset "Consecutive Infinite Layer Checks" begin
    frequencies = [50.0]
    @testset "Horizontal Model" begin
        model = EarthModel(frequencies, 100.0, 10.0, 1.0, t=Inf)
        # It's an error to add any layer after an infinite one in a horizontal model
        @test_throws ErrorException add!(model, frequencies, 200.0, 15.0, 1.0, t=10.0)
        @test_throws ErrorException add!(model, frequencies, 200.0, 15.0, 1.0, t=Inf)
    end

    @testset "Vertical Model" begin
        # Setup: model with two earth layers, the second being infinite
        model = EarthModel(frequencies, 100.0, 10.0, 1.0, t=Inf, vertical_layers=true)
        add!(model, frequencies, 150.0, 12.0, 1.0, t=20.0) # Add one more layer
        add!(model, frequencies, 150.0, 12.0, 1.0, t=Inf) # Add one more layer, now Inf

        # It's an error to add another infinite layer
        @test_throws ErrorException add!(model, frequencies, 200.0, 15.0, 1.0, t=Inf)

        # It should not be posssible to add a finite layer after an infinite one
        @test_throws ErrorException add!(model, frequencies, 300.0, 20.0, 1.0, t=5.0)
        @test length(model.layers) == 4
        @test model.layers[4].base_rho_g == 150.0
        @test isinf(model.layers[4].t) # The last layer should still be infinite
    end
end

@testset "show method for EarthModel" begin
    frequencies = [50.0]
    # Homogeneous model
    model_homo = EarthModel(frequencies, 100.0, 10.0, 1.0)
    s_homo = sprint(show, "text/plain", model_homo)
    @test contains(s_homo, "EarthModel with 1 horizontal earth layer (homogeneous)")
    @test contains(s_homo, "└─ Layer 2: [rho_g=100.0, epsr_g=10.0, mur_g=1.0, t=∞]")

    # Multilayer horizontal model
    model_multi_h = EarthModel(frequencies, 100.0, 10.0, 1.0, t=20.0)
    add!(model_multi_h, frequencies, 200.0, 15.0, 1.0, t=Inf)
    s_multi_h = sprint(show, "text/plain", model_multi_h)
    @test contains(s_multi_h, "EarthModel with 2 horizontal earth layers (multilayer)")
    @test contains(s_multi_h, "├─ Layer 2: [rho_g=100.0, epsr_g=10.0, mur_g=1.0, t=20.0]")
    @test contains(s_multi_h, "└─ Layer 3: [rho_g=200.0, epsr_g=15.0, mur_g=1.0, t=∞]")

    # Multilayer vertical model
    model_multi_v = EarthModel(frequencies, 100.0, 10.0, 1.0, t=Inf, vertical_layers=true)
    add!(model_multi_v, frequencies, 200.0, 15.0, 1.0, t=30.0)
    s_multi_v = sprint(show, "text/plain", model_multi_v)
    @test contains(s_multi_v, "EarthModel with 2 vertical earth layers (multilayer)")
    @test contains(s_multi_v, "├─ Layer 2: [rho_g=100.0, epsr_g=10.0, mur_g=1.0, t=∞]")
    @test contains(s_multi_v, "└─ Layer 3: [rho_g=200.0, epsr_g=15.0, mur_g=1.0, t=30.0]")
end

@testset "DataFrame for EarthModel" begin
    frequencies = [50.0]
    model = EarthModel(frequencies, 100.0, 10.0, 1.0, t=20.0)
    add!(model, frequencies, 200.0, 15.0, 1.0, t=Inf)

    df = DataFrame(model)
    @test df isa DataFrame
    @test names(df) == ["rho_g", "epsr_g", "mur_g", "thickness"]
    @test nrow(df) == 3

    # Air layer
    @test isinf(df.rho_g[1])
    @test df.epsr_g[1] == 1.0
    @test df.mur_g[1] == 1.0
    @test isinf(df.thickness[1])

    # First earth layer
    @test df.rho_g[2] == 100.0
    @test df.epsr_g[2] == 10.0
    @test df.thickness[2] == 20.0

    # Second earth layer
    @test df.rho_g[3] == 200.0
    @test df.epsr_g[3] == 15.0
    @test isinf(df.thickness[3])
end

@testset "show for EarthModel" begin
    frequencies = [50.0, 60.0]

    @testset "Homogeneous Horizontal" begin
        model = EarthModel(frequencies, 100.0, 10.0, 1.0)
        str_repr = sprint(show, "text/plain", model)
        @test occursin("EarthModel with 1 horizontal earth layer (homogeneous) and 2 frequency samples", str_repr)
        @test occursin("└─ Layer 2:", str_repr)
        @test occursin("t=∞", str_repr)
        @test occursin("Frequency-dependent model: CP model", str_repr)
    end

    @testset "Multilayer Horizontal" begin
        model = EarthModel(frequencies, 100.0, 10.0, 1.0, t=20.0)
        add!(model, frequencies, 200.0, 15.0, 1.0, t=Inf)
        str_repr = sprint(show, "text/plain", model)
        @test occursin("EarthModel with 2 horizontal earth layers (multilayer) and 2 frequency samples", str_repr)
        @test occursin("├─ Layer 2:", str_repr)
        @test occursin("└─ Layer 3:", str_repr)
        @test occursin("t=20", str_repr)
    end

    @testset "Multilayer Vertical" begin
        model = EarthModel(frequencies, 100.0, 10.0, 1.0, t=Inf, vertical_layers=true)
        add!(model, frequencies, 200.0, 15.0, 1.0, t=5.0)
        str_repr = sprint(show, "text/plain", model)
        @test occursin("EarthModel with 2 vertical earth layers (multilayer) and 2 frequency samples", str_repr)
        @test occursin("t=5", str_repr)
    end
end
