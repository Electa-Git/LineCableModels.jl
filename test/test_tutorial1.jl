@testitem "examples/tutorial1.jl tests" setup = [defaults] begin

    # Helpers
    function material_approx_equal(m::Material, rho, eps_r, mu_r, T0, alpha; atol=1e-12, rtol=1e-8)
        return isapprox(m.rho, rho; atol=atol, rtol=rtol) &&
               isapprox(m.eps_r, eps_r; atol=atol, rtol=rtol) &&
               isapprox(m.mu_r, mu_r; atol=atol, rtol=rtol) &&
               isapprox(m.T0, T0; atol=atol, rtol=rtol) &&
               isapprox(m.alpha, alpha; atol=atol, rtol=rtol)
    end

    @testset "initialize and inspect" begin
        materials = MaterialsLibrary()  # default initialization as in the tutorial
        @test materials !== nothing

        # DataFrame conversion
        df = DataFrame(materials)
        @test isa(df, DataFrame)
        @test nrow(df) >= 0  # should be defined (>=0); more specific tests below

        # Check expected columns if present
        expected = ["name", "rho", "eps_r", "mu_r", "T0", "alpha"]
        @test all(x -> x in string.(names(df)), expected)
    end

    @testset "add materials from tutorial" begin
        materials = MaterialsLibrary(add_defaults=false)  # start clean for deterministic tests

        # Define tutorial materials (subset representative of file)
        copper_corrected = Material(1.835e-8, 1.0, 0.999994, 20.0, 0.00393)
        aluminum_corrected = Material(3.03e-8, 1.0, 0.999994, 20.0, 0.00403)
        epr = Material(1e15, 3.0, 1.0, 20.0, 0.005)
        pvc = Material(1e15, 8.0, 1.0, 20.0, 0.1)

        add!(materials, "copper_corrected", copper_corrected)
        add!(materials, "aluminum_corrected", aluminum_corrected)
        add!(materials, "epr", epr)
        add!(materials, "pvc", pvc)

        # Verify that keys exist
        for name in ("copper_corrected", "aluminum_corrected", "epr", "pvc")
            @test haskey(materials, name)
        end

        # DataFrame contains the names
        df = DataFrame(materials)
        @test "copper_corrected" in df.name
        @test "epr" in df.name
    end

    @testset "remove duplicate" begin
        materials = MaterialsLibrary(add_defaults=false)
        epr = Material(1e15, 3.0, 1.0, 20.0, 0.005)
        add!(materials, "epr", epr)

        # Add duplicate and then remove it
        add!(materials, "epr_dupe", epr)
        @test haskey(materials, "epr_dupe")
        delete!(materials, "epr_dupe")
        @test !haskey(materials, "epr_dupe")
    end

    @testset "save and load round-trip (temp file)" begin
        materials = MaterialsLibrary(add_defaults=false)

        # Add a small set of materials
        copper_corrected = Material(1.835e-8, 1.0, 0.999994, 20.0, 0.00393)
        epr = Material(1e15, 3.0, 1.0, 20.0, 0.005)
        add!(materials, "copper_corrected", copper_corrected)
        add!(materials, "epr", epr)

        tmpfile = tempname() * ".json"
        try
            # Save to temporary file
            save(materials, file_name=tmpfile)
            @test isfile(tmpfile)

            # Load into a fresh library
            loaded = MaterialsLibrary(add_defaults=false)
            load!(loaded, file_name=tmpfile)

            # Keys present after load
            @test haskey(loaded, "copper_corrected")
            @test haskey(loaded, "epr")

            # Retrieve and compare properties
            copper_loaded = get(loaded, "copper_corrected")
            @test isa(copper_loaded, Material)
            @test material_approx_equal(copper_loaded, 1.835e-8, 1.0, 0.999994, 20.0, 0.00393)

            epr_loaded = get(loaded, "epr")
            @test isa(epr_loaded, Material)
            @test material_approx_equal(epr_loaded, 1e15, 3.0, 1.0, 20.0, 0.005)

        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end

    @testset "error handling" begin
        # Fresh empty library (no defaults) for deterministic error behavior
        empty_lib = MaterialsLibrary(add_defaults=false)

        # get on non-existent key should display alert
        @test get(empty_lib, "non_existent_material") === nothing


        # delete! on non-existent key should throw KeyError
        @test_throws KeyError delete!(empty_lib, "non_existent_material")

        # load! from a non-existent file should throw an I/O-related error (SystemError / IOError)
        bad_file_lib = MaterialsLibrary(add_defaults=false)
        @test_throws Exception load!(bad_file_lib, file_name="this_file_should_not_exist_hopefully_0123456789.json")
    end

    @testset "integration-like workflow (safe, uses temp files)" begin
        # Recreate the main tutorial workflow but using temporary save path
        materials = MaterialsLibrary(add_defaults=false)

        # Add the full tutorial list used in examples (representative)
        add!(materials, "copper_corrected", Material(1.835e-8, 1.0, 0.999994, 20.0, 0.00393))
        add!(materials, "aluminum_corrected", Material(3.03e-8, 1.0, 0.999994, 20.0, 0.00403))
        add!(materials, "lead", Material(21.4e-8, 1.0, 0.999983, 20.0, 0.00400))
        add!(materials, "steel", Material(13.8e-8, 1.0, 300.0, 20.0, 0.00450))
        add!(materials, "bronze", Material(3.5e-8, 1.0, 1.0, 20.0, 0.00300))
        add!(materials, "stainless_steel", Material(70.0e-8, 1.0, 500.0, 20.0, 0.0))
        add!(materials, "epr", Material(1e15, 3.0, 1.0, 20.0, 0.005))
        add!(materials, "pvc", Material(1e15, 8.0, 1.0, 20.0, 0.1))
        add!(materials, "laminated_paper", Material(1e15, 2.8, 1.0, 20.0, 0.0))
        add!(materials, "carbon_pe", Material(0.06, 1e3, 1.0, 20.0, 0.0))
        add!(materials, "conductive_paper", Material(18.5, 8.6, 1.0, 20.0, 0.0))

        # Duplicate add and delete
        add!(materials, "epr_dupe", get(materials, "epr"))
        @test haskey(materials, "epr_dupe")
        delete!(materials, "epr_dupe")
        @test !haskey(materials, "epr_dupe")

        tmpfile = tempname() * ".json"
        try
            save(materials, file_name=tmpfile)
            @test isfile(tmpfile)

            reloaded = MaterialsLibrary(add_defaults=false)
            load!(reloaded, file_name=tmpfile)

            # verify a representative sample of materials exists after reload
            for name in ("copper_corrected", "pvc", "stainless_steel")
                @test haskey(reloaded, name)
            end

            # Verify copper properties after reload
            copper = get(reloaded, "copper_corrected")
            @test material_approx_equal(copper, 1.835e-8, 1.0, 0.999994, 20.0, 0.00393)

        finally
            isfile(tmpfile) && rm(tmpfile)
        end
    end

end