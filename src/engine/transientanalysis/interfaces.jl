"""
    respond!(formula, workspace, parameters, options)

Fill the workspace buffers and return the time axis \\[s\\], voltages \\[V\\]
and currents \\[A\\].
"""
function respond! end

"""
$(TYPEDEF)

Own the input, plan and buffers for one transient computation.

$(TYPEDFIELDS)
"""
struct TransientWorkspace{I, P, B}
    "Physical input."
    input::I
    "Computation plan."
    plan::P
    "Preallocated arrays."
    buffers::B
end

initialize_buffers(::Formula, ::Type, input, plan, buffers) = buffers
