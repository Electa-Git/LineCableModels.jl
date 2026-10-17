"""
    fit(formula, workspace, parameters, options)

Fit the selected line model and return `(model, fitting_error, passivity, dc, diagnostics)`.
The `diagnostics` NamedTuple contains formula-specific values copied into the result details.
"""
function fit end

"""
$(TYPEDEF)

Own the input, plan and buffers for one line-model computation.

$(TYPEDFIELDS)
"""
struct LineModelWorkspace{I, P, B}
    "Physical input."
    input::I
    "Computation plan."
    plan::P
    "Preallocated arrays."
    buffers::B
end

initialize_buffers(::Formula, ::Type, input, plan, buffers) = buffers
