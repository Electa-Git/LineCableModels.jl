using Gauntlet
using TOML

length(ARGS) in (1, 2) || error(
    "usage: source_consistency.jl <raw-dataset> [evidence.toml]",
)
record = Gauntlet.ordinary_consistency(:pscad, ARGS[1])
payload = Dict(String(key) => value for (key, value) in pairs(record))
payload["schema_version"] = 1
if length(ARGS) == 2
    open(ARGS[2], "w") do io
        TOML.print(io, payload; sorted = true)
    end
else
    TOML.print(stdout, payload; sorted = true)
end
