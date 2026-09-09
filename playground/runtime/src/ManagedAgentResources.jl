"""
    ManagedAgentResources(parent; terminal_limits=TerminalSessionLimits())

Compose scientific and private-terminal state owners over one managed physical
driver, journal and capacity limit. Both borrow the same live agent ledger.
Closing a partition does not close its sibling; closing this root always attempts
both partitions and the parent, retaining unresolved ownership for a retry.

Construction is passive and does not make unavailable profiles eligible. The
managed driver's installed-profile checks remain the admission authority.
"""
struct ManagedAgentResources <: AbstractAgentResources
    "Sole owner of physical handles, command runner and recovery journal."
    parent::ManagedResourceDriver
    "Scientific preparation, execution and cancellation state."
    science::ScientificResources{ManagedScientificView}
    "Private terminal byte streams and independent session deadlines."
    terminals::TerminalResources{ManagedTerminalView}
end
function ManagedAgentResources(parent::ManagedResourceDriver;terminal_limits=TerminalSessionLimits())
    ManagedAgentResources(parent,ScientificResources(ManagedScientificView(parent)),
        TerminalResources(ManagedTerminalView(parent);limits=terminal_limits))
end
Base.show(io::IO,::ManagedAgentResources)=print(io,"ManagedAgentResources(<owned>)")
installed_profiles(resources::ManagedAgentResources)=installed_profiles(resources.parent)

# The agent's scientific services consume only this partition. Adding a second
# process kind does not create another control loop or another lease ledger.
scientific_resources(::AbstractAgentResources)=nothing
scientific_resources(resources::ScientificResources)=resources
scientific_resources(resources::ManagedAgentResources)=resources.science
terminal_resources(::AbstractAgentResources)=nothing
terminal_resources(resources::TerminalResources)=resources
terminal_resources(resources::ManagedAgentResources)=resources.terminals

function bind_agent!(resources::ManagedAgentResources,ledger::AgentLeaseLedger)
    bind_agent!(resources.science,ledger)
    bind_agent!(resources.terminals,ledger)
    return nothing
end
function recover_owned!(resources::ManagedAgentResources)
    recover_owned!(resources.science)
    recover_owned!(resources.terminals)
    return nothing
end
function release_owned!(resources::ManagedAgentResources,fence::AssignmentFence)
    profile=get(installed_profiles(resources).definitions,fence.profile_id,nothing)
    profile !== nothing && (string(profile.version),profile.fingerprint)==(fence.profile_version,fence.fingerprint) ||
        throw(ArgumentError("managed cleanup profile differs"))
    return release_managed_kind!(Val(profile.kind),resources,fence)
end
release_managed_kind!(::Val{:scientific},resources,fence)=release_owned!(resources.science,fence)
release_managed_kind!(::Val{:terminal},resources,fence)=release_owned!(resources.terminals,fence)
function Base.close(resources::ManagedAgentResources)
    # Revoke both partitions before waiting for either physical teardown.
    lock(()->resources.science.closed=true,resources.science.lock)
    lock(()->resources.terminals.closed=true,resources.terminals.lock)
    completed=true
    for owner in (resources.science,resources.terminals,resources.parent)
        try close(owner) catch;completed=false end
    end
    completed || throw(ArgumentError("managed agent cleanup remains unresolved"))
    return nothing
end
