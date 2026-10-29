"""
$(TYPEDSIGNATURES)

Return the run directory of `problem` under `.linecablemodels/fem` in the current directory.
Replace each character outside `[A-Za-z0-9_-]` in the system identifier with `_`. Append the
local time and four random lower-case letters or digits. The function does not create the
directory.
"""
function directory(problem::LineParametersProblem)
    id = replace(problem.system.id, r"[^A-Za-z0-9_-]" => "_")
    stamp = Libc.strftime("%Y%m%d-%H%M%S", time())
    suffix = String(rand(['a':'z'; '0':'9'], 4))
    return joinpath(pwd(), ".linecablemodels", "fem", id, stamp * "-" * suffix)
end

"""
$(TYPEDSIGNATURES)

Scan the problem's frequencies with GetDP and apply the bundle and Kron reductions to its
primitive matrices. Compute the shunt admittance from the reduced potential coefficient:

```math
Y = P_{\\mathrm{r}}^{-1}.
```

Here ``P_{\\mathrm{r}}`` is the reduced potential coefficient \\[Ω·m\\]. Apply the ideal
transposition to both result matrices when selected.

# Arguments

- `Z`: destination for the series impedances in \\[Ω/m\\], by row, column, and frequency.
- `Y`: destination shunt admittances \\[S/m\\], in the same order.
- `trace`: destinations for primitive `Z` \\[Ω/m\\] and `P` \\[Ω·m\\] by terminal.
- `directory`: run directory with inputs and meshes.
- `model`: field model with the problem and terminal order.
- `formulation`: material laws and reduction controls.
- `buffers`: preallocated `factor`, `coupling`, `rhs`, `reduced`, and `identity` arrays.
- `options`: native solver controls and the selection of the trace.

# Returns

- `Z`, with computed entries. Also fill `Y` and the selected trace arrays.
"""
function solve!(Z::AbstractArray{Complex{T}, 3}, Y::AbstractArray{Complex{T}, 3}, trace,
        directory::AbstractString, model::FieldModel{T}, formulation::FEM, buffers,
        options::ComputationOptions{GmshBackend}) where {T <: Real}
    primitives = scan(directory, model, options)
    if options.trace
        copyto!(trace.Z, primitives.Z)
        copyto!(trace.P, primitives.P)
    end
    (; bundles, kept, eliminated) = reductions(model.problem.system, formulation.options)
    (; factor, coupling, rhs, reduced, identity) = buffers
    for k in eachindex(model.problem.frequencies)
        Zk = view(Z,:,:,k)
        Yk = view(Y,:,:,k)
        primitive_Z = view(primitives.Z,:,:,k)
        primitive_P = view(primitives.P,:,:,k)
        bundle!(primitive_Z, bundles)
        bundle!(primitive_P, bundles)
        eliminate!(Zk, primitive_Z, kept, eliminated, factor, coupling, rhs)
        eliminate!(reduced, primitive_P, kept, eliminated, factor, coupling, rhs)
        ldiv!(Yk, lu!(reduced), identity)
        if formulation.options.transposition
            transposition!(Zk)
            transposition!(Yk)
        end
    end
    return Z
end
