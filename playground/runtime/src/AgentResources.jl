using .Protocol: AssignmentFence

"""
    AbstractAgentResources

Implement installed_profiles, bind_agent!, recover_owned!, release_owned! and close for an
approved resource supervisor. These hooks handle owned process/container resources,
not scientific code evaluation. The agent root retains authority, scheduling and
cleanup acknowledgement; an application cannot replace those steps.
"""
abstract type AbstractAgentResources end

@required AbstractAgentResources begin
    installed_profiles(::AbstractAgentResources)
    bind_agent!(::AbstractAgentResources, ::AgentLeaseLedger)
    recover_owned!(::AbstractAgentResources)
    release_owned!(::AbstractAgentResources, ::AssignmentFence)
    Base.close(::AbstractAgentResources)
end

"""
    installed_profiles(resources) -> ProfileRegistry

Return the locally verified subset of approved profiles without preparing an
executor. An empty registry keeps the agent visible but unselectable.
"""
function installed_profiles end

"""
    bind_agent!(resources, ledger)

Bind the resource owner to the agent's exact lease authority before recovery.
This step must be passive and return nothing. A resource owner cannot later be
rebound to a different worker incarnation or a substitute lease ledger.
"""
function bind_agent! end

"""
    recover_owned!(resources)

Reconcile only this agent's receipt-identified resources before its first broker
announcement. Return nothing only on successful recovery; a failure prevents
startup. Never infer live authority from old process IDs or saved warm labels.
"""
function recover_owned! end

"""
    release_owned!(resources, fence) -> Bool

Stop the exact assignment's owned resources under their configured bounds.
Return true only after physical cleanup, false if unresolved. The agent runs one
cleanup task per occupied lease, so this work never blocks heartbeat handling.
Exceptions or false results keep capacity occupied and are retried at a bounded rate.
"""
function release_owned! end

function verified_agent_profiles(config::AgentConfig, resources::AbstractAgentResources)
    RequiredInterfaces.check_interface_implemented(AbstractAgentResources, typeof(resources)) === true ||
        throw(ArgumentError("agent resource supervisor does not implement its required hooks"))
    installed = installed_profiles(resources)
    installed isa ProfileRegistry || throw(ArgumentError("installed_profiles must return ProfileRegistry"))
    for (id, profile) in installed.definitions
        approved = get(config.profiles.definitions, id, nothing)
        approved !== nothing && all(name -> getfield(approved, name) == getfield(profile, name),
            fieldnames(ProfileDefinition)) || throw(ArgumentError("agent advertised an unapproved profile"))
    end
    return installed
end
