"""
    RuntimeRequirement(role, profiles; required=true)

Declare one runtime role and its compatible approved profile IDs. This record
does not reserve a worker or perform preparation.
"""
struct RuntimeRequirement
    "Application-local role identifier."
    role::String
    "Compatible registered profile identifiers, in preference order."
    profiles::Tuple{Vararg{String}}
    "Whether preparation is required for the application's live cases."
    required::Bool
    function RuntimeRequirement(role, profiles; required::Bool=true)
        choices = Tuple(checked_id.(profiles))
        !isempty(choices) && allunique(choices) ||
            throw(ArgumentError("a runtime role needs distinct compatible profiles"))
        return new(checked_id(role), choices, required)
    end
end

"""
    ApplicationDefinition(id, title, kind, entrypoint; version=v"1.0.0",
        visibility=:public, entry_surface=(kind == :presentation ? :published : :ui),
        requirements=())

Describe a registered presentation or workbench using passive metadata.
Registration never constructs its Bonito views or starts its workers.
Entry surface distinguishes a published document containing owned frames from
a page rendered directly by its UI host, independently of catalogue category.
"""
struct ApplicationDefinition
    "Stable application identifier."
    id::String
    "Human-readable catalogue title."
    title::String
    "Either presentation or workbench."
    kind::Symbol
    "Same-origin published entry path."
    entrypoint::String
    "Application contract version."
    version::VersionNumber
    "Public or developer catalogue visibility."
    visibility::Symbol
    "Published document or UI-host entry point."
    entry_surface::Symbol
    "Declared runtime roles; empty for a static application."
    requirements::Tuple{Vararg{RuntimeRequirement}}

    function ApplicationDefinition(id, title, kind, entrypoint;
            version::VersionNumber=v"1.0.0", visibility::Symbol=:public,
            entry_surface::Symbol=(kind == :presentation ? :published : :ui), requirements=())
        kind in (:presentation, :workbench) || throw(ArgumentError("invalid application kind"))
        visibility in (:public, :developer) || throw(ArgumentError("invalid visibility"))
        entry_surface in (:published, :ui) || throw(ArgumentError("invalid entry surface"))
        title isa AbstractString && 0 < length(title) <= 160 ||
            throw(ArgumentError("application title needs 1–160 characters"))
        entrypoint isa AbstractString && startswith(entrypoint, "/") &&
            !startswith(entrypoint, "//") && !occursin(r"[\\?#\s]", entrypoint) &&
            !any(==(".."), split(entrypoint, '/')) ||
            throw(ArgumentError("entrypoint must be a same-origin absolute path"))
        roles = Tuple(requirements)
        all(r -> r isa RuntimeRequirement, roles) && allunique(r.role for r in roles) ||
            throw(ArgumentError("runtime requirements must have unique roles"))
        return new(checked_id(id), String(title), kind, String(entrypoint),
            version, visibility, entry_surface, roles)
    end
end

"""
    HostContext

Describe one server-owned UI launch. It contains no user code or scientific
parameters. The run directory is created by the supervisor, not the browser.
"""
struct HostContext
    "Run identity used for authorization and route isolation."
    run_id::UUID
    "Same-origin route prefix for all assets and WebSockets."
    prefix::String
    "Private owned file where the child publishes readiness."
    ready_file::String
end

"""
    AbstractApplication

Extend describe and ui_command for a trusted application registration. The
coordinator owns admission and process cleanup independently of these hooks.
"""
abstract type AbstractApplication end

@required AbstractApplication begin
    describe(::AbstractApplication)
    ui_command(::AbstractApplication, ::HostContext)
end

"""
    LocalApplication(definition, project, script)

Use an approved on-disk Julia environment and script as an isolated UI host.
Neither path is taken from a browser request.
"""
struct LocalApplication <: AbstractApplication
    "Passive catalogue declaration."
    definition::ApplicationDefinition
    "Absolute approved Julia project directory."
    project::String
    "Absolute approved Julia startup script."
    script::String
    function LocalApplication(definition::ApplicationDefinition, project, script)
        project_path, script_path = abspath(project), abspath(script)
        isfile(joinpath(project_path, "Project.toml")) && isfile(script_path) ||
            throw(ArgumentError("application project and startup script must exist"))
        return new(definition, project_path, script_path)
    end
end
describe(app::LocalApplication) = app.definition
function ui_command(app::LocalApplication, context::HostContext)
    return `$(Base.julia_cmd()) --startup-file=no --history-file=no --compiled-modules=existing --project=$(app.project) $(app.script)`
end

"""
    ApplicationRegistry()

Keep trusted Julia application implementations by stable ID. Only passive
descriptions are exposed to the browser; implementations are never serialized.
"""
struct ApplicationRegistry
    "Approved definitions and their trusted implementation hooks."
    applications::Dict{String,AbstractApplication}
    "Immutable metadata captured once, without re-running declaration hooks."
    definitions::Dict{String,ApplicationDefinition}
end
ApplicationRegistry() = ApplicationRegistry(Dict{String,AbstractApplication}(),
    Dict{String,ApplicationDefinition}())

"""
    register!(registry, application) -> ApplicationDefinition

Verify required hooks, reject duplicate IDs and store a cheap declaration.
No UI command is invoked. Missing hooks or duplicate IDs raise ArgumentError.
"""
function register!(registry::ApplicationRegistry, app::AbstractApplication)
    RequiredInterfaces.check_interface_implemented(AbstractApplication, typeof(app)) === true ||
        throw(ArgumentError("application does not implement its required hooks"))
    description = describe(app)
    description isa ApplicationDefinition ||
        throw(ArgumentError("describe must return ApplicationDefinition"))
    haskey(registry.definitions, description.id) &&
        throw(ArgumentError("application ID is already registered"))
    registry.applications[description.id] = app
    registry.definitions[description.id] = description
    return description
end

"""
    register!(registry::ApplicationRegistry, definition::ApplicationDefinition)

Publish a passive application declaration with no installed UI implementation.
It remains readable in catalogues; launch is rejected before reserving capacity.
"""
function register!(registry::ApplicationRegistry, definition::ApplicationDefinition)
    haskey(registry.definitions, definition.id) && throw(ArgumentError("application ID is already registered"))
    registry.definitions[definition.id] = definition
    return definition
end
