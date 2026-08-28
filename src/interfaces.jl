"Add one owned value to a mutable collection."
function add! end

"""
    build(Target, declarations...; kwargs...)

Construct one completed domain object from complete physical declarations.

`build(CableDesign, ...)` resolves a physical cable root and terminal state.
`build(LineCableSystem, ...)` places completed designs and resolves global
connections. When a direct argument is a `Grid` or `Gridspace`, the same
surface returns a `Gridspace{Target}` whose callable invokes scalar `build`
after selecting one point.

# Arguments

- `Target`: Completed domain type to construct.
- `declarations`: Complete physical declarations owned by `Target`.

# Returns

One completed `Target`, or `Gridspace{Target}` for explicit finite inputs.
"""
function build end

"Return the physical storage basis of a result."
function basis end

function R end
function L end
function C end
function resistance end
function inductance end
function capacitance end
