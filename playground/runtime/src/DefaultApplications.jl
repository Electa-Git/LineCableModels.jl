include(joinpath(@__DIR__, "..", "..", "src", "applications", "Catalogue.jl"))

"""
    CatalogueApplication

Bind a passive entry to the fixed owned UI driver. Only registered symbols
select implementations; neither browser data nor a request supplies a command.
"""
struct CatalogueApplication <: AbstractApplication
    "Immutable application declaration."
    definition::ApplicationDefinition
    "Approved driver selection."
    driver::Symbol
    "Whether owned-component diagnostics are enabled."
    xray::Bool
end
describe(app::CatalogueApplication) = app.definition
function ui_command(app::CatalogueApplication, ::HostContext)
    root = normpath(joinpath(@__DIR__, "..", ".."))
    driver = joinpath(root, "src", "applications", "ui_driver.jl")
    # Package installation/precompilation is an operator step. A run may use
    # existing caches or compile in its own process, never leave precompile
    # subprocesses behind when its supervisor is killed.
    return `$(Base.julia_cmd()) --startup-file=no --history-file=no --compiled-modules=existing --project=$root $driver $(String(app.driver)) $(app.xray)`
end

"""
    default_applications(; xray=false) -> ApplicationRegistry

Read the shared catalogue without loading a UI package. Entries whose live
implementation has not landed remain passive, with no implicit fallback.
"""
function default_applications(; xray::Bool=false)
    registry = ApplicationRegistry()
    for row in ApplicationCatalogue.entries()
        roles = Tuple(RuntimeRequirement(r.role, r.profiles; required=r.required) for r in row.requirements)
        definition = ApplicationDefinition(row.id, row.title, row.kind, row.entrypoint;
            version=row.version, visibility=row.visibility, entry_surface=row.entry_surface, requirements=roles)
        if row.ui in (:template, :starter, :specimen, :gallery, :showcase, :cable_study)
            register!(registry, CatalogueApplication(definition, row.ui, xray))
        else
            register!(registry, definition)
        end
    end
    return registry
end
