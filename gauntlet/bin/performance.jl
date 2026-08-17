using Gauntlet
using TOML

corpus = Corpus(:smoke)
configuration = TOML.parsefile(joinpath(@__DIR__, "..", "config", "performance.toml"))
suite = Suite(
    :performance;
    corpus,
    ids = configuration["cases"],
    checks = (),
    performance = true
)
report = gauntlet(suite)
builtin_trials = Trial[gauntlet(
                           Gauntlet.performance_case(Val(Symbol(name)));
                           performance = true,
                           checks = ()
                       ) for name in get(configuration, "builtins", String[])]
report = Report(
    report.suite,
    vcat(builtin_trials, report.trials),
    report.started_at,
    report.finished_at
)
output = joinpath(@__DIR__, "..", "reports")
mkpath(output)
write_report(joinpath(output, "performance.json"), report)
display(report)
