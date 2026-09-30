# Controlled investigation, never called by the engine. Keep old FEM references.
include("finite_gamma_terminal_shift.jl")

function compare_shifted(dir, problem, form; reference="baseline-terminal-shift",
        candidate="localized-terminal-shift", label="terminal-shift")
    analytical = compute(problem, Formulation(
        earth_impedance=formula(:unified; options=(Γ=form.options.data.Γ,)),
        earth_admittance=formula(:unified; options=(Γ=form.options.data.Γ,)); options=REDUCTIONS))
    worst = 0.; flips = 0; old_change = 0.
    open(joinpath(dir,"$label-components.csv"),"w") do io
        println(io,"quantity,i,j,original_fem,reference,candidate,analytical,relative_change,sign_match")
        for (q,native,component) in (("R","Z",real),("X","Z",imag),("G","Y",real),("B","Y",imag))
            old = component.(native_matrix(joinpath(dir,"baseline"),native))
            a = component.(native_matrix(joinpath(dir,reference),native))
            b = component.(native_matrix(joinpath(dir,candidate),native))
            ref = component.(native=="Z" ? Z(analytical)[:,:,1] : Y(analytical)[:,:,1])
            for j in 1:2, i in 1:2
                relative = iszero(a[i,j]) ? (iszero(b[i,j]) ? 0. : Inf) : abs((b[i,j]-a[i,j])/a[i,j])
                old_change = max(old_change, abs((a[i,j]-old[i,j])/old[i,j]))
                match = sign(a[i,j])==sign(b[i,j]); flips += !match
                worst = max(worst,relative)
                println(io,join((q,i,j,old[i,j],a[i,j],b[i,j],ref[i,j],relative,match),','))
            end
        end
    end
    result = Dict("reference"=>reference,"candidate"=>candidate,
        "maximum_component_relative_change"=>worst,"new_sign_changes"=>flips,
        "reference_change_from_original"=>old_change,"relative_tolerance"=>.02,
        "passed"=>worst<=.02 && flips==0)
    record(joinpath(dir,"$label-comparison.toml"),result)
    say("SHIFTED COMPARISON ",basename(dir)," ",result)
    return result["passed"]
end

function qualify_terminal_shift()
    cases = [(:mixed,1e6,.99),(:mixed,.1,.99),(:air,.1,.99),(:air,1e6,.99),
        (:soil,.1,.99),(:soil,1e6,.99),(:mixed,.1,0.)]
    for (layout,f,fraction) in cases
        dir = joinpath(ROOT,"$layout-f$(f)-gamma$(fraction)")
        problem, form = fixture(layout,f,fraction)
        baseline = joinpath(dir,"baseline")
        target = joinpath(dir,"localized")
        if !isfile(joinpath(target,"sources.toml"))
            for file in [readlines(joinpath(baseline,".onelab-export-files"));
                    ".onelab-export-files";"sources.toml"]
                mkpath(dirname(joinpath(target,file)))
                cp(joinpath(baseline,file),joinpath(target,file);force=true)
            end
        end
        a = TOML.parsefile(joinpath(baseline,"mesh.toml"))
        @assert digest(joinpath(baseline,"study.msh"))==a["sha256"]
        b = mesh!(target,problem,true;media=:both)
        for key in ("pml_coordinates","path_coordinates","contour_coordinates")
            @assert a[key]==b[key]
        end
        terminal_shift_control(baseline)
        terminal_shift_control(target)
        compare_shifted(dir,problem,form) || error("Shifted-mesh qualification failed: $dir")
    end
    zero_case = joinpath(ROOT,"mixed-f0.1-gamma0.0")
    for candidate in ("baseline","localized"), quantity in ("Z","Y","P")
        @assert native_matrix(joinpath(zero_case,candidate),quantity) ==
            native_matrix(joinpath(zero_case,candidate*"-terminal-shift"),quantity)
    end
    for layout in (:air,:soil,:mixed), f in (.1,1e3,1e6)
        @assert TOML.parsefile(joinpath(ROOT,"$layout-f$(f)-gamma0.0",
            "two-percent-comparison.toml"))["passed"]
    end
    record(joinpath(ROOT,"both-media-corrected-selection.toml"),Dict(
        "passed"=>true,"candidate"=>"localized","relative_tolerance"=>.02,
        "zero_gamma_cases"=>9,"finite_gamma_cases"=>6,
        "finite_gamma_reference"=>"exact terminal-drive substitution on original frozen meshes",
        "zero_gamma_identity"=>true,"original_references_preserved"=>true))
    say("COMPLETE seven exact-substitution pairs; original references preserved")
end

if abspath(PROGRAM_FILE) == @__FILE__
    qualify_terminal_shift()
end
