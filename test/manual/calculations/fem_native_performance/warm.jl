isdefined(@__MODULE__,:PILOTS) || include("qualify.jl")

function copied_run!(source,target)
    isfile(joinpath(target,"sources.toml")) && return target
    for f in [readlines(joinpath(source,".onelab-export-files"));
            ".onelab-export-files";"study.msh";"mesh.toml";"sources.toml"]
        dst=joinpath(target,f); mkpath(dirname(dst)); cp(joinpath(source,f),dst)
    end
    target
end

function warm!()
    passed=[v for v in VARIANTS if v ∉ ("binary","amd") && (v=="baseline" || all(PILOTS) do c
        TOML.parsefile(joinpath(ROOT,label(c),"comparison-$v.toml"))["passed"]
    end)]
    if all(v->v in passed,("harmonic","physical3"))
        combined=true
        for case in PILOTS
            dir=prepare!(case,"combined"); B.solve!(dir)
            problem,form=B.fixture(case...)
            combined &= B.compare!(dirname(dir),problem,form;candidate="combined",tolerance=.02)
        end
        if combined
            # Time the intended combined production change. Individual terms
            # already have separate pilot costs; do not repeat those campaigns.
            filter!(v->v ∉ ("harmonic","physical3"),passed)
            push!(passed,"combined")
        end
    end
    say("WARMED MATCHED RUNS: ",passed,"; all failed candidates excluded")
    for repeat in 1:2, case in PILOTS, variant in (isodd(repeat) ? passed : reverse(passed))
        dir=copied_run!(joinpath(ROOT,label(case),variant),
            joinpath(ROOT,"warm-$repeat",label(case),variant))
        B.solve!(dir;ordering=variant=="amd" ? 0 : nothing)
    end
    totals=Dict{String,Vector{Float64}}()
    worst=0.; flips=0
    open(joinpath(ROOT,"warmed-costs.csv"),"w") do io
        println(io,"repeat,case,variant,native_seconds,assembly_seconds,solve_seconds,peak_rss_kib,dofs")
        for variant in passed
            totals[variant]=Float64[]
            for repeat in 1:2
                total=0.
                for case in PILOTS
                    d=TOML.parsefile(joinpath(ROOT,"warm-$repeat",label(case),variant,"solve.toml"))
                    for quantity in ("Z","Y"), part in (real,imag)
                        a=part.(B.native_matrix(joinpath(ROOT,label(case),"baseline"),quantity))
                        b=part.(B.native_matrix(joinpath(ROOT,"warm-$repeat",label(case),variant),quantity))
                        for (x,y) in zip(a,b)
                            worst=max(worst,iszero(x) ? (iszero(y) ? 0. : Inf) : abs((y-x)/x))
                            flips+=sign(x)!=sign(y)
                        end
                    end
                    total+=d["seconds"]
                    println(io,join((repeat,label(case),variant,d["seconds"],d["assembly_seconds"],
                        d["solve_seconds"],d["peak_rss_kib"],d["dofs"]),','))
                end
                push!(totals[variant],total)
            end
        end
    end
    decisions=Dict{String,Any}()
    record(joinpath(ROOT,"warmed-preservation.toml"),Dict("maximum_component_relative_change"=>worst,
        "new_component_sign_changes"=>flips,"passed"=>worst<=.02 && flips==0))
    worst<=.02 && flips==0 || error("A warmed solve failed preservation")
    for variant in passed
        speed=1-median(totals[variant])/median(totals["baseline"])
        decisions[variant]=Dict("median_native_workload_seconds"=>median(totals[variant]),
            "relative_saving"=>speed,"passes_20_percent_target"=>speed>=.2)
        say("WARMED COST ",variant," ",decisions[variant])
    end
    record(joinpath(ROOT,"warmed-assessment.toml"),decisions)
    say("COMPLETE warmed native measurements; mesh conversion and Julia compilation reported separately")
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && warm!()
