module ComponentXRay

using Bonito
using JSON3
using Observables
using UUIDs

export ActionInspection,
    BindingInspection,
    ComponentInspection,
    PropertyInspection,
    SourceReference,
    XRayPolicy,
    inspectable,
    inspection,
    inspection_payload,
    install,
    instrument,
    set_policy!,
    source_reference,
    xray_policy

const PLAYGROUND_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const XRAY_SCRIPT_PATH = joinpath(@__DIR__, "component_xray.js")
const XRAY_STYLES_PATH = joinpath(@__DIR__, "component_xray.css")
include_dependency(XRAY_SCRIPT_PATH)
include_dependency(XRAY_STYLES_PATH)
const XRAY_SCRIPT = read(XRAY_SCRIPT_PATH, String)
const XRAY_STYLES = read(XRAY_STYLES_PATH, String)

"""Control whether component diagnostics are emitted and initially active."""
struct XRayPolicy
    permitted::Bool
    enabled::Bool

    function XRayPolicy(; permitted::Bool=false, enabled::Bool=false)
        enabled && !permitted && throw(ArgumentError(
            "X-ray diagnostics cannot be enabled when they are not permitted"
        ))
        return new(permitted, enabled)
    end
end

XRayPolicy(enabled::Bool) = XRayPolicy(; permitted=enabled, enabled)

struct SourceReference
    module_name::String
    file::String
    line::Int
end

function source_reference(
        module_name::Union{Module,AbstractString},
        file::AbstractString,
        line::Integer
    )
    normalized = normpath(string(file))
    relative = try
        relpath(normalized, PLAYGROUND_ROOT)
    catch
        normalized
    end
    startswith(relative, "..") && (relative = basename(normalized))
    return SourceReference(string(module_name), replace(relative, '\\' => '/'), Int(line))
end

function short_type(value)
    type = value isa Type ? value : typeof(value)
    name = string(nameof(type))
    owner = string(parentmodule(type))
    return owner in ("Base", "Core", "Main") ? name : "$owner.$name"
end

function summarize(value; limit::Int=96)
    rendered = if value === nothing
        "nothing"
    elseif value isa Symbol
        ":$(value)"
    elseif value isa AbstractString
        repr(string(value))
    elseif value isa Tuple
        "(" * join((summarize(item; limit=24) for item in value), ", ") * ")"
    elseif value isa AbstractVector
        "$(length(value))-element $(short_type(value))"
    elseif value isa AbstractDict
        "$(length(value))-entry $(short_type(value))"
    else
        sprint(show, value; context=:limit => true)
    end
    return length(rendered) <= limit ? rendered : first(rendered, limit - 1) * "…"
end

struct PropertyInspection
    name::String
    value::String
    julia_type::String
    origin::String
end

function PropertyInspection(
        name::Union{Symbol,AbstractString},
        value;
        origin::Union{Symbol,AbstractString}=:configured
    )
    return PropertyInspection(
        string(name),
        summarize(value),
        short_type(value),
        string(origin)
    )
end

struct BindingInspection
    name::String
    observable::Any
    julia_type::String
    notes::String
end

function BindingInspection(
        name::Union{Symbol,AbstractString},
        observable::Observables.AbstractObservable;
        notes::AbstractString=""
    )
    return BindingInspection(
        string(name),
        observable,
        short_type(observable),
        string(notes)
    )
end

struct ActionInspection
    event::String
    name::String
    callback_type::String
    owner::String
    disabled::Bool
end

function ActionInspection(
        event::Union{Symbol,AbstractString},
        name::Union{Symbol,AbstractString},
        callback;
        disabled::Bool=false
    )
    callback_type = callback === nothing ? "Nothing" : short_type(callback)
    owner = callback === nothing ? "—" : string(parentmodule(typeof(callback)))
    return ActionInspection(
        string(event),
        string(name),
        callback_type,
        owner,
        disabled
    )
end

struct ComponentInspection
    name::String
    julia_type::String
    source::SourceReference
    parameters::Vector{PropertyInspection}
    bindings::Vector{BindingInspection}
    actions::Vector{ActionInspection}
    css_scopes::Vector{String}
    notes::Vector{String}
end

function ComponentInspection(
        component;
        name::Union{Symbol,AbstractString}=nameof(typeof(component)),
        source::SourceReference=source_reference(
            parentmodule(typeof(component)),
            "unknown",
            0
        ),
        parameters=PropertyInspection[],
        bindings=BindingInspection[],
        actions=ActionInspection[],
        css_scopes=String[],
        notes=String[]
    )
    return ComponentInspection(
        string(name),
        short_type(component),
        source,
        collect(PropertyInspection, parameters),
        collect(BindingInspection, bindings),
        collect(ActionInspection, actions),
        string.(collect(css_scopes)),
        string.(collect(notes))
    )
end

"""Return diagnostic metadata for an owned component, or `nothing`."""
inspection(::Any) = nothing

inspectable(component) = !isnothing(inspection(component))

const POLICY_LOCK = ReentrantLock()
const SESSION_POLICIES = WeakKeyDict{Any,XRayPolicy}()

function set_policy!(session, policy::XRayPolicy)
    lock(POLICY_LOCK) do
        SESSION_POLICIES[session] = policy
    end
    return policy
end

function xray_policy(session)
    return lock(POLICY_LOCK) do
        get(SESSION_POLICIES, session, XRayPolicy())
    end
end

binding_value(binding::BindingInspection) = try
    summarize(binding.observable[])
catch error
    "unavailable ($(nameof(typeof(error))))"
end

function inspection_payload(
        inspection::ComponentInspection;
        id::AbstractString=string(uuid4())
    )
    return Dict{String,Any}(
        "id" => string(id),
        "name" => inspection.name,
        "julia_type" => inspection.julia_type,
        "source" => Dict(
            "module" => inspection.source.module_name,
            "file" => inspection.source.file,
            "line" => inspection.source.line,
        ),
        "parameters" => [
            Dict(
                "name" => parameter.name,
                "value" => parameter.value,
                "julia_type" => parameter.julia_type,
                "origin" => parameter.origin,
            ) for parameter in inspection.parameters
        ],
        "bindings" => [
            Dict(
                "name" => binding.name,
                "value" => binding_value(binding),
                "julia_type" => binding.julia_type,
                "notes" => binding.notes,
            ) for binding in inspection.bindings
        ],
        "actions" => [
            Dict(
                "event" => action.event,
                "name" => action.name,
                "callback_type" => action.callback_type,
                "owner" => action.owner,
                "disabled" => action.disabled,
            ) for action in inspection.actions
        ],
        "css_scopes" => inspection.css_scopes,
        "notes" => inspection.notes,
    )
end

function javascript_literal(value)
    literal = String(JSON3.write(value))
    literal = replace(literal, "</" => "<\\/")
    literal = replace(literal, '\u2028' => "\\u2028", '\u2029' => "\\u2029")
    return literal
end

function registration_script(payload)
    serialized = javascript_literal(payload)
    return """
    (() => {
        const script = document.currentScript;
        const node = script.parentElement.firstElementChild;
        const entry = {
            kind: 'register',
            node,
            metadata: $serialized
        };
        if (window.__lcmXRayDispatch) {
            window.__lcmXRayDispatch(entry);
        } else {
            (window.__lcmXRayPending ||= []).push(entry);
        }
    })();
    """
end

function binding_script(node, binding::BindingInspection)
    name = binding.name
    return js"""
    (() => {
        const target = $(node);
        const binding = $(binding.observable);
        binding.on(value => {
            const entry = {
                kind: 'binding',
                node: target,
                name: $(name),
                value
            };
            if (window.__lcmXRayDispatch) {
                window.__lcmXRayDispatch(entry);
            } else {
                (window.__lcmXRayPending ||= []).push(entry);
            }
        });
    })();
    """
end

function instrument(session, node, component)
    descriptor = inspection(component)
    isnothing(descriptor) && return node
    return instrument(session, node, descriptor)
end

function instrument(session, node, descriptor::ComponentInspection)
    policy = xray_policy(session)
    policy.permitted || return node

    id = string(uuid4())
    payload = inspection_payload(descriptor; id)
    scripts = Any[DOM.script(registration_script(payload))]
    append!(
        scripts,
        (DOM.script(binding_script(node, binding)) for binding in descriptor.bindings)
    )
    return DOM.div(
        node,
        scripts...;
        style="display: contents;",
        class="lc-xray-mount"
    )
end

function install(session, root)
    policy = xray_policy(session)
    policy.permitted || return nothing
    options = javascript_literal(Dict(
        "enabled" => policy.enabled,
        "shortcut" => "Ctrl+Shift+X",
    ))
    styles = javascript_literal(XRAY_STYLES)
    installer = """
    (() => {
        const page = document.currentScript.closest('.lc-wb-page');
        const root = page && page.querySelector('[data-lc-workbench]');
        if (!root) return;
        window.LineCableModelsComponentXRay.installXRay(
            root,
            $options,
            $styles
        );
    })();
    """
    return DOM.div(
        DOM.script(XRAY_SCRIPT),
        DOM.script(installer);
        style="display: contents;"
    )
end

end
