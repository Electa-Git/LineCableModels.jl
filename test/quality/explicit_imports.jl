@testitem "Quality / explicit imports / package ownership" tags = [:quality] begin
    using ExplicitImports: test_explicit_imports, improper_qualified_accesses
    # These adapters participate in the numerical/UQ paths. Check the loaded
    # extensions too; a cold root-only scan cannot see their ownership errors.
    import Measurements, Distributions, Gmsh, Calculus, XLSX, CairoMakie
    for name in (:LineCableModelsMeasurementsExt, :LineCableModelsDistributionsExt,
            :LineCableModelsGmshExt, :LineCableModelsXLSXExt)
        @test Base.get_extension(LineCableModels, name) !== nothing
    end
    renderer = Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
    cairo = Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt)
    @test renderer !== nothing
    @test cairo !== nothing
    import Logging, JSON3

    # Gmsh's generated API deliberately uses qualified calls without export/
    # public annotations. Accept its documented bindings, not arbitrary private
    # names. The same applies to the documented logger and JSON protocols below.
    function documented_fem_access(owner::Module, name::Symbol)
        # Generated API/build metadata needed to identify the numerical library.
        owner === Gmsh.gmsh && name in (:GMSH_API_VERSION, :lib) && return true
        ancestor = owner
        while ancestor !== Gmsh.gmsh && parentmodule(ancestor) !== ancestor
            ancestor = parentmodule(ancestor)
        end
        allowed = (ancestor === Gmsh.gmsh && !startswith(string(name), "_")) ||
            (owner === Gmsh && name === :finalize) ||
            (owner === Logging && name in (:catch_exceptions, :handle_message,
                :min_enabled_level, :shouldlog)) ||
            (owner === JSON3 && name in (:read, :write, :pretty))
        allowed && isdefined(owner, name) || return false
        value = getfield(owner, name)
        return Base.Docs.hasdoc(owner, name) ||
            (value isa Function && Base.Docs.hasdoc(parentmodule(value), nameof(value)))
    end

    # Exact external contracts, classified in docs/src/developers.md. Public
    # annotations alone miss documented qualified APIs. The one layout removal
    # workaround is explicitly recommended by its upstream maintainer; it is
    # not a general permission to consume Makie or package-owned internals.
    function external_contract(consumer, owner, name)
        consumer === cairo && owner === CairoMakie && name === :activate! && return true
        consumer === renderer && owner === CairoMakie.Makie && name in (
            :automatic, :current_backend, :get_ticks, :get_tickvalues,
            :pseudolog10, :fast_string_boundingboxes) && return true
        consumer === renderer && owner === CairoMakie.Makie.GridLayoutBase &&
            name === :remove_from_gridlayout! && return true
        consumer === renderer && owner === Base && name in (:require, :IOError) && return true
        return false
    end

    # All other ExplicitImports checks remain unchanged, including ownership.
    test_explicit_imports(LineCableModels; all_qualified_accesses_are_public=false)
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    unexpected = String[]
    for (consumer, accesses) in improper_qualified_accesses(LineCableModels; skip=())
        for row in accesses
            row.public_access && continue
            row.self_qualified && continue # Covered by the existing separate check.
            row.accessing_from === Base && Base.ispublic(Core, row.name) && continue
            row.accessing_from === Logging && documented_fem_access(Logging, row.name) && continue
            consumer === extension && documented_fem_access(row.accessing_from, row.name) && continue
            external_contract(consumer, row.accessing_from, row.name) && continue
            push!(unexpected, "$(row.accessing_from).$(row.name) at $(row.location)")
        end
    end
    @test unexpected == String[]
    @test documented_fem_access(Gmsh.gmsh.model, :add_physical_group)
    @test !documented_fem_access(Gmsh.gmsh.model, :_unregistered_helper)
    @test !documented_fem_access(LineCableModels.Engine, :compute)
    @test !documented_fem_access(JSON3, :StructTypes)
    @test external_contract(renderer, CairoMakie.Makie, :get_ticks)
    @test !external_contract(LineCableModels.Engine, CairoMakie.Makie, :get_ticks)
    @test !external_contract(renderer, CairoMakie.Makie, :get_plot_visibilities)
    @test !external_contract(renderer, CairoMakie.Makie, :fast_string_boundingboxes_obs)
    @test !external_contract(renderer, LineCableModels.Engine, :compute)
end


@testitem "Quality / explicit imports / Gauntlet ownership" tags=[:quality] setup=[GauntletSupport] begin
    using ExplicitImports: test_explicit_imports, improper_qualified_accesses
    using .GauntletSupport: Gauntlet
    import Pkg
    import Logging
    using LinearAlgebra: BLAS

    file = joinpath(pkgdir(LineCableModels), "gauntlet", "Gauntlet.jl")
    test_explicit_imports(Gauntlet, file; all_qualified_accesses_are_public=false)
    # Native facilities with no public annotation: detect instrumented runs,
    # identify artifact hashes, restore recorded packages, and record BLAS settings.
    # IOError is the exception raised by native I/O, including interrupted pipes.
    # These are external interfaces; package-owned private access has no exception.
    native = (
        Base => (:JLOptions, :PkgId, :SHA1, :include, :loaded_modules, :require,
            :structdiff, :extension_parent_name, :IOError),
        Base.Filesystem => (:path_separator,),
        Pkg => (:dependencies,),
        BLAS => (:get_config, :get_num_threads),
        # Required AbstractLogger protocol methods have no public annotation.
        Logging => (:catch_exceptions, :handle_message, :min_enabled_level, :shouldlog),
    )
    unexpected = String[]
    for (consumer, accesses) in improper_qualified_accesses(Gauntlet, file; skip=())
        for row in accesses
            row.public_access && continue
            row.self_qualified && continue
            row.accessing_from === Base && Base.ispublic(Core, row.name) && continue
            any(owner === row.accessing_from && row.name in names for (owner, names) in native) && continue
            push!(unexpected, "$(row.accessing_from).$(row.name) at $(row.location)")
        end
    end
    @test unexpected == String[]
end
