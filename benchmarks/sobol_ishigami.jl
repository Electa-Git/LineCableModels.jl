using Distributions
using GlobalSensitivity
using LineCableModels
using QuasiMonteCarlo

const Grammar = LineCableModels.Grammar
const PB = LineCableModels.ParametricBuilder

struct IshigamiProblem <: Grammar.AbstractProblemDefinition
    x1::Float64
    x2::Float64
    x3::Float64
end

struct IshigamiFormulation <: Grammar.AbstractFormulation
    solves::Base.RefValue{Int}
end

struct IshigamiResult
    value::Float64
end

LineCableModels.validate(problem::IshigamiProblem) = problem

function Grammar.compute(
        problem::IshigamiProblem,
        formulation::IshigamiFormulation;
        options::NamedTuple = (;)
)
    isempty(options) || error("unexpected Ishigami computation options")
    formulation.solves[] += 1
    a = 7.0
    b = 0.1
    value = sin(problem.x1) + a * sin(problem.x2)^2 +
            b * problem.x3^4 * sin(problem.x1)
    return IshigamiResult(value)
end

Grammar.observe(result::IshigamiResult, ::typeof(R)) = result.value
Grammar.observables(::Type{IshigamiResult}) = (R,)

sigma = pi / sqrt(3.0)
space = PB.Gridspace{IshigamiProblem}(
    IshigamiProblem,
    (
        PB.Grid(0.0, PB.AbsoluteError(sigma)),
        PB.Grid(0.0, PB.AbsoluteError(sigma)),
        PB.Grid(0.0, PB.AbsoluteError(sigma))
    )
)
problem = ParametricProblem(space)
solves = Ref(0)
samples = 32_768
evaluations = 8samples
analysis = Sensitivity(
    IshigamiFormulation(solves),
    GlobalSensitivity.Sobol(order = [0, 1, 2]),
    (R,);
    samples,
    distribution = :uniform,
    input_labels = ("x1", "x2", "x3"),
    max_evaluations = evaluations,
    max_output_values = evaluations
)

elapsed = @elapsed result = compute(problem, analysis)
product = only(result)
estimated_first = only(product.first_order)
estimated_total = only(product.total_order)
estimated_second = only(product.second_order)

a = 7.0
b = 0.1
variance = 0.5 + a^2 / 8 + b * pi^4 / 5 + b^2 * pi^8 / 18
variance_first = 0.5 * (1 + b * pi^4 / 5)^2
variance_second = a^2 / 8
variance_interaction = 8b^2 * pi^8 / 225
reference_first = [
    variance_first / variance,
    variance_second / variance,
    0.0
]
reference_total = [
    (variance_first + variance_interaction) / variance,
    variance_second / variance,
    variance_interaction / variance
]
reference_pairs = [0.0, variance_interaction / variance, 0.0]
estimated_pairs = [
    estimated_second[1, 2],
    estimated_second[1, 3],
    estimated_second[2, 3]
]

first_error = maximum(abs.(estimated_first .- reference_first))
total_error = maximum(abs.(estimated_total .- reference_total))
second_error = maximum(abs.(estimated_pairs .- reference_pairs))

@assert solves[] == product.evaluations == evaluations
@assert product.inputs == ["x1", "x2", "x3"]
@assert product.active_indices == [1, 2, 3]
@assert first_error <= 5.0e-3
@assert total_error <= 5.0e-3
@assert second_error <= 5.0e-3

println("Ishigami Sobol benchmark")
println("samples=$samples evaluations=$evaluations elapsed_seconds=$elapsed")
println("first_order=$estimated_first")
println("total_order=$estimated_total")
println("second_order_pairs=$estimated_pairs")
println("max_errors=(first=$first_error, total=$total_error, second=$second_error)")
