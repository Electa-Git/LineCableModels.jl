@testitem "Quality / explicit imports / package ownership" tags = [:quality] begin
    using ExplicitImports: test_explicit_imports, improper_qualified_accesses
    # These adapters participate in the numerical/UQ paths. Check the loaded
    # extensions too; a cold root-only scan cannot see their ownership errors.
    import Measurements, Distributions, Gmsh
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

    # All other ExplicitImports checks remain unchanged, including ownership.
    test_explicit_imports(LineCableModels; all_qualified_accesses_are_public=false)
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    unexpected = String[]
    for (consumer, accesses) in improper_qualified_accesses(LineCableModels; skip=())
        for row in accesses
            row.public_access && continue
            row.self_qualified && continue # Covered by the existing separate check.
            row.accessing_from === Base && Base.ispublic(Core, row.name) && continue
            consumer === extension && documented_fem_access(row.accessing_from, row.name) && continue
            push!(unexpected, "$(row.accessing_from).$(row.name) at $(row.location)")
        end
    end
    @test isempty(unexpected)
    @test documented_fem_access(Gmsh.gmsh.model, :add_physical_group)
    @test !documented_fem_access(Gmsh.gmsh.model, :_unregistered_helper)
    @test !documented_fem_access(LineCableModels.Engine, :compute)
    @test !documented_fem_access(JSON3, :StructTypes)
end
