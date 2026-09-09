# Explicit no-executor fixture for the agent's resource-ownership contract.
mutable struct AgentResourceFixture <: AbstractAgentResources
    profiles::ProfileRegistry
    recovered::Int
    release_allowed::Base.RefValue{Bool}
    release_fails::Base.RefValue{Bool}
    recovery_fails::Bool
    closed::Bool
    released::Vector{RT.Protocol.AssignmentFence}
end
AgentResourceFixture(profiles) = AgentResourceFixture(profiles, 0, Ref(true), Ref(false),
    false, false, RT.Protocol.AssignmentFence[])
RT.installed_profiles(resource::AgentResourceFixture) = resource.profiles
RT.bind_agent!(::AgentResourceFixture, ::AgentLeaseLedger) = nothing
function RT.recover_owned!(resource::AgentResourceFixture)
    resource.recovery_fails && throw(ArgumentError("fixture recovery failed"))
    resource.recovered += 1
    return nothing
end
function RT.release_owned!(resource::AgentResourceFixture, fence::RT.Protocol.AssignmentFence)
    while !resource.release_allowed[] && !resource.closed
        sleep(0.01)
    end
    resource.release_fails[] && return false
    push!(resource.released, fence)
    return true
end
function Base.close(resource::AgentResourceFixture)
    resource.closed = true
    return nothing
end

struct IncompleteAgentResources <: AbstractAgentResources end
