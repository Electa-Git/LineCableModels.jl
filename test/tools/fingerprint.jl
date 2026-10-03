# The recording program of the equivalence check (`equivalence.jl`). It records the
# preservation corpus of the revision that contains it. Each revision records with its
# own version of this file.
# Run `julia --project=test test/tools/fingerprint.jl DESTINATION`.
#
# For each scenario of `preservation_corpus(2)` it writes one line per node of the
# inputs and the result, with the scenario, the path, the type and the value. Numbers are
# written as their bits, a Measurement as value and error. Each type name is written as
# `Owner.Name`, where `Owner` is the module that defines it directly, whatever Main
# imports. Base and Core types are bare. Anonymous functions are written as
# `var"#closure"`, because their numbers change when unrelated source changes, and no
# type alias is used. Random identities (UUIDs) are numbered by first appearance within
# a scenario, which keeps which nodes share one. Dictionaries are walked in sorted key
# order. An object met again records the path where it was first met. A last pass
# counts, per parametric scenario, the points it materializes and the blueprint lowerings
# of its designs.
using LineCableModels, Measurements
import LineCableModels.Engine as EN

Base.include(@__MODULE__, joinpath(@__DIR__, "..", "support", "scenarios.jl"))
using .CurrentScenarios: preservation_corpus

# Counts blueprint lowerings of designs whose nominal data is a counter, then runs the
# production lowering, as `test/integration/formulation_grid.jl` does.
struct LoweringCounter
    calls::Base.RefValue{Int}
end
function EN.flatten(engine::LineCableModelsCoaxial,
        design::CableDesign{T, R, G, NamedTuple{(:counter,), Tuple{LoweringCounter}}},
        ::Type{S}, methods::NamedTuple, solutions::Vector, design_index::Int) where {
        T <: Real, R <: AbstractCablePart, G <: LineCableModels.DataModel.CableGeometry, S <: Real}
    design.nominal_data.counter.calls[] += 1
    return invoke(EN.flatten,
        Tuple{LineCableModelsCoaxial, CableDesign, Type{S}, NamedTuple, Vector, Int},
        engine, design, S, methods, solutions, design_index)
end

const CLOSURE = "var\"#closure\""

bare(m::Module) = Base.moduleroot(m) in (Base, Core)

# A type parameter that is a value, such as a Symbol or a tuple.
parameter(x) = x isa Type || x isa TypeVar ? typename(x) :
    x isa Tuple ? string("(", join(parameter.(x), ", "), length(x) == 1 ? ",)" : ")") :
    repr(x)

function typename(@nospecialize(T))
    T isa TypeVar && return string(T.name)
    T isa Union && return string("Union{", join(typename.(Base.uniontypes(T)), ", "), "}")
    T isa UnionAll && return string(typename(T.body), " where ", T.var.name)
    T isa Core.TypeofVararg &&
        return string("Vararg{", isdefined(T, :T) ? typename(T.T) : "", "}")
    T isa DataType || return repr(T)
    name = string(T.name.name)
    owner = parentmodule(T)
    if startswith(name, "#")
        function_name = match(r"^#([^#]+)$", name)
        function_name === nothing && return CLOSURE
        return string("typeof(", bare(owner) ? "" : string(nameof(owner), "."),
            function_name[1], ")")
    end
    qualified = bare(owner) ? name : string(nameof(owner), ".", name)
    isempty(T.parameters) && return qualified
    return string(qualified, "{", join(parameter.(T.parameters), ", "), "}")
end

bits(x) = bytes2hex(reinterpret(UInt8, [x]))
flat(text) = replace(text, r"\s+" => " ")

# `state.seen` maps each mutable object met to its first path. `state.uuids` numbers the
# random identities (gridpoint source ids) by first appearance.
function walk!(io, path, x, state)
    line(value) = println(io, path, '\t', typename(typeof(x)), '\t', value)
    if x isa Base.UUID
        line(string("uuid #", get!(state.uuids, x, length(state.uuids) + 1)))
    elseif x isa Measurement
        line(string(bits(Measurements.value(x)), " ± ", bits(Measurements.uncertainty(x))))
    elseif x isa Complex
        line("")
        walk!(io, path * ".re", real(x), state)
        walk!(io, path * ".im", imag(x), state)
    elseif x isa Number && isbitstype(typeof(x))
        line(bits(x))
    elseif x isa Union{AbstractString, Symbol, Char, Nothing, Missing}
        line(repr(x))
    elseif x isa Type
        line(typename(x))
    elseif x isa Module
        line(string(nameof(x)))
    elseif x isa Ptr
        line("pointer")
    elseif ismutable(x) && haskey(state.seen, x)
        line("same as " * state.seen[x])
    else
        ismutable(x) && (state.seen[x] = path)
        if x isa AbstractDict
            line(string(length(x)))
            for key in sort!(collect(keys(x)); by = key -> sprint(show, key))
                walk!(io, string(path, "[", flat(sprint(show, key)), "]"), x[key], state)
            end
        elseif x isa AbstractArray
            line(string(size(x)))
            for index in eachindex(x)
                isassigned(x, index) && walk!(io, string(path, "[", index, "]"), x[index], state)
            end
        else
            line("")
            for name in fieldnames(typeof(x))
                isdefined(x, name) && walk!(io, string(path, ".", name), getfield(x, name), state)
            end
        end
    end
    return nothing
end

function counted(name)
    materialized = Ref(0)
    counter = LoweringCounter(Ref(0))
    counters = (materialized = () -> (materialized[] += 1), nominal_data = (; counter))
    compute(preservation_corpus(2; counters)[name]...)
    return (materialized[], counter.calls[])
end

function main(destination)
    open(destination, "w") do io
        corpus = preservation_corpus(2)
        for name in keys(corpus)
            arguments = corpus[name]
            out = IOBuffer()
            state = (seen = IdDict{Any, String}(), uuids = Dict{Base.UUID, Int}())
            try
                for (i, argument) in enumerate(arguments)
                    println(out, "inputs[", i, "].eltype\tType\t", typename(eltype(argument)))
                    walk!(out, "inputs[$i]", argument, state)
                end
                result = compute(arguments...)
                println(out, "result.eltype\tType\t", typename(eltype(result)))
                walk!(out, "result", result, state)
            catch error
                println(out, "result\tfailed\t", first(flat(sprint(showerror, error)), 300))
            end
            for text in eachline(IOBuffer(take!(out)))
                println(io, name, '\t', text)
            end
        end
        # The instrumented pass follows, so it cannot change the recorded values.
        for name in (:product, :zip, :linear_error, :monte_carlo)
            value = try
                join(counted(name), " materializations, ") * " lowerings"
            catch error
                "failed: " * first(flat(sprint(showerror, error)), 300)
            end
            println(io, name, "\tlaziness\tcounts\t", value)
        end
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(only(ARGS))
