isdefined(@__MODULE__,:PILOTS) || include("qualify.jl")

function spectrum!(variant)
    variant in ("harmonic","physical3","combined") || error("Select a qualified assembly prescription")
    cases=[(layout,f,gamma) for gamma in (0.,.99), layout in (:air,:soil,:mixed),
        f in (.1,1e3,1e6) if gamma==0 || f!=1e3]
    passed=true
    for case in cases
        problem,form=B.fixture(case...)
        base=prepare!(case,"baseline"); B.solve!(base)
        dir=prepare!(case,variant); B.solve!(dir)
        passed &= B.compare!(dirname(dir),problem,form;candidate=variant,tolerance=.02)
    end
    record(joinpath(ROOT,"spectrum-$variant.toml"),Dict("cases"=>length(cases),"passed"=>passed))
    say("COMPLETE spectrum ",variant," passed=",passed,"; fifteen prescribed cases")
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && spectrum!(only(ARGS))
