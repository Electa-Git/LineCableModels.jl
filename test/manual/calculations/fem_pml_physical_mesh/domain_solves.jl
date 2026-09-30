# Four authorized serial native solves on the prepared, frozen domain meshes.
# All matrix formation/integration remains in GetDP; this only reads results.
using TOML, SHA, Printf, Dates
const ROOT = normpath(joinpath(@__DIR__,"../../../../.linecablemodels/fem/pml-physical-mesh"))
const PROBES = joinpath(ROOT,"domain-probes")
say(xs...) = (println(Dates.now()," ",xs...);flush(stdout))

function native_matrix(path)
    matrix = zeros(ComplexF64,2,2)
    for row in split.(readlines(path)[3:end],'\t')
        matrix[parse(Int,row[1]),parse(Int,row[2])] = complex(parse(Float64,row[5]),parse(Float64,row[6]))
    end
    matrix
end
function reference_matrix(path,quantity,frequency)
    matrix = fill(complex(NaN,NaN),2,2)
    for row in split.(readlines(path)[2:end],',')
        row[1]==quantity && parse(Float64,row[2])==frequency || continue
        matrix[parse(Int,row[3]),parse(Int,row[4])] = complex(parse(Float64,row[5]),parse(Float64,row[6]))
    end
    all(isfinite,matrix) || error("Missing saved $quantity reference at $frequency Hz")
    matrix
end

function solve(label)
    directory = joinpath(PROBES,label)
    mesh_info = TOML.parsefile(joinpath(directory,"mesh.toml"))
    mesh = joinpath(directory,"study.msh")
    mesh_info["mesh_sha256"] == bytes2hex(open(sha256,mesh)) || error("Prepared mesh changed: $label")
    index,frequency = mesh_info["frequency_index"],mesh_info["frequency_hz"]
    output = joinpath(directory,"results",@sprintf("f%04d-quasi-fw-b0000",index))
    marker = joinpath(directory,"solve.toml")
    if !isfile(joinpath(output,"completed.txt"))
        say("DOMAIN SOLVE BEGIN ",label," f=",frequency," Hz; two source columns")
        executable = joinpath(ROOT,"getdp-live")
        entry = joinpath(directory,"study.pro")
        command = `$executable $entry -msh $mesh -solve LineCableModelsFEM -setnumber FrequencyIndex $index -setnumber Physics 1 -setnumber PlotFieldMaps 0 -nt 1 -v 4`
        write(joinpath(directory,"command.txt"),string(command)*"\n")
        open(joinpath(directory,"solver.log"),"w") do io
            run(pipeline(`/usr/bin/time -f %e,%M -o $(joinpath(directory,"resources.csv")) $command`;
                stdout=io,stderr=io))
        end
    end
    parse.(Float64,split(read(joinpath(output,"completed.txt"),String))) ==
        [index,frequency,1,0] || error("Wrong native completion marker: $label")
    reference = joinpath(ROOT,"qualification/radius-0.085")
    signs = 0
    self_r_error = 0.
    g_relative_error = 0.
    g12 = NaN
    open(joinpath(directory,"components.csv"),"w") do io
        println(io,"quantity,frequency_hz,i,j,value,reference,signed_error,absolute_error,relative_error,retained_value,change_from_retained,sign_match")
        for (quantity,native) in (("R","Z"),("X","Z"),("G","Y"),("B","Y"))
            component = quantity in ("R","G") ? real : imag
            values = component.(native_matrix(joinpath(output,"matrices",native*".tsv")))
            expected = component.(reference_matrix(joinpath(reference,"reference.csv"),native,frequency))
            retained = component.(reference_matrix(joinpath(reference,"matrices.csv"),native,frequency))
            for j in 1:2,i in 1:2
                value,target = values[i,j],expected[i,j]
                relative = iszero(target) ? NaN : abs((value-target)/target)
                match = sign(value)==sign(target)
                quantity=="G" && !iszero(target) && !match && (signs+=1)
                quantity=="R" && i==j && (self_r_error=max(self_r_error,relative))
                quantity=="G" && (g_relative_error=max(g_relative_error,relative))
                quantity=="G" && i==1 && j==2 && (g12=value)
                println(io,join((quantity,frequency,i,j,value,target,value-target,abs(value-target),
                    relative,retained[i,j],value-retained[i,j],match),','))
            end
        end
    end
    transcript = read(joinpath(directory,"solver.log"),String)
    dofs = maximum(parse(Int,m.captures[1]) for m in eachmatch(r"System \d+/\d+: ([0-9]+) Dofs",transcript))
    resources = parse.(Float64,split(strip(read(joinpath(directory,"resources.csv"),String)),','))
    timings = [parse.(Float64,split(strip(read(joinpath(output,"raw/jobs",
        @sprintf("getdp-f%04d-b%04d-timing.tsv",index,basis)),String)))) for basis in 1:2]
    result = Dict("case"=>label,"frequency_hz"=>frequency,"G_sign_mismatches"=>signs,
        "G12"=>g12,"maximum_G_relative_error"=>g_relative_error,"maximum_self_R_relative_error"=>self_r_error,
        "native_seconds"=>resources[1],"peak_rss_kib"=>round(Int,resources[2]),"dofs"=>dofs,
        "constraint_seconds"=>sum(row[3] for row in timings),
        "assembly_seconds"=>sum(row[4] for row in timings),
        "solve_seconds"=>sum(row[5] for row in timings),
        "output_seconds"=>sum(row[6] for row in timings),
        "mesh_sha256"=>mesh_info["mesh_sha256"],"completed_source_columns"=>2)
    open(io -> TOML.print(io,result),marker*".tmp","w")
    mv(marker*".tmp",marker;force=true)
    say("DOMAIN SOLVE DONE ",label," G12=",g12," G sign mismatches=",signs,
        " self-R relative error=",self_r_error," native seconds=",resources[1]," DOFs=",dofs)
    result
end

# Default is the B pair only, allowing inspection before the C pair.
stage = isempty(ARGS) ? "b" : only(ARGS)
stage in ("b","c") || error("Select b (2delta / 24delta) or c (2delta / 2delta)")
extent = stage=="b" ? "l24" : "l2"
for name in ("low","mid")
    solve("d2-$extent-$name")
end
say("DOMAIN STAGE ",uppercase(stage)," COMPLETE; results preserved without clipping or sign correction")
