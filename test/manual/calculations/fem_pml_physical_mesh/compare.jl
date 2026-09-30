# Postprocessing only; never launches FEM or clips observations.
using TOML, Printf, Dates
function matrix_rows(path;frequency=nothing)
    result=Dict{Tuple{String,Int,Int},ComplexF64}()
    for line in eachline(path)
        startswith(line,"quantity,") && continue
        row=split(line,',')
        if length(row)==6
            frequency===nothing && error("frequency required for scan table")
            isapprox(parse(Float64,row[2]),frequency;rtol=1e-14) || continue
            result[(row[1],parse(Int,row[3]),parse(Int,row[4]))]=complex(parse(Float64,row[5]),parse(Float64,row[6]))
        else
            result[(row[1],parse(Int,row[2]),parse(Int,row[3]))]=complex(parse(Float64,row[4]),parse(Float64,row[5]))
        end
    end
    result
end
function compare_completed(root)
    baseline=joinpath(dirname(root),"pml-conductance-cost/qualification/radius-0.085")
    open(joinpath(root,"components.csv"),"w") do io
        println(io,"case,frequency_hz,quantity,i,j,value,reference,retained_144,signed_error,absolute_error,relative_error,change_from_144,sign_match")
        for directory in sort(readdir(root;join=true))
            marker=joinpath(directory,"complete.toml")
            isfile(marker) || continue
            completion=TOML.parsefile(marker);frequency=completion["frequency"]
            value=matrix_rows(joinpath(directory,"matrices.csv"))
            reference=matrix_rows(joinpath(baseline,"reference.csv");frequency)
            retained=matrix_rows(joinpath(baseline,"matrices.csv");frequency)
            for (quantity,primitive,component) in (("R","Z",real),("X","Z",imag),("G","Y",real),("B","Y",imag))
                errors=Float64[];wrong=0
                for j in 1:2,i in 1:2
                    key=(primitive,i,j)
                    v,r,b=component(value[key]),component(reference[key]),component(retained[key])
                    err=v-r;relative=abs(err)/abs(r);matched=sign(v)==sign(r)
                    push!(errors,relative);wrong+=!matched
                    println(io,join((basename(directory),frequency,quantity,i,j,v,r,b,err,abs(err),relative,v-b,matched),','))
                end
                println(Dates.now()," COMPARE ",basename(directory)," ",quantity,
                    " maximum_relative_error=",maximum(errors)," sign_disagreements=",wrong)
            end
        end
    end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    compare_completed(isempty(ARGS) ? joinpath(pwd(),".linecablemodels/fem/pml-physical-mesh") : abspath(only(ARGS)))
end
