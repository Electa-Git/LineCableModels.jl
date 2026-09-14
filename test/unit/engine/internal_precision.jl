@testitem "Engine / internal impedance / independent radial diffusion control" tags=[:unit] begin
    include(joinpath(pkgdir(LineCableModels),"test/support/radial_control.jl"))
    using .RadialControl
    using TOML
    diagnostics=Dict{String,Any}[]
    const II=LineCableModels.Engine.InternalImpedance
    selected=II.Formula(:default)
    for inner in (0.0,.003),f in (10.0,1000.0)
        references=map((128,256,512)) do bits
            setprecision(BigFloat,bits) do
                ri,ro,rho=BigFloat(inner),big".005",big"2e-8"
                RadialControl.surfaces(ri,ro,rho,4big(pi)*big"1e-7",2big(pi)*im*f)
            end
        end
        reference=last(references)
        for kind in (:inner,:outer,:transfer)
            expected=getproperty(reference,kind)
            uncertainty=reference.bound+abs(expected-getproperty(references[2],kind))
            @testset "$T $kind ri=$inner f=$f" for T in (Float32,Float64,BigFloat)
                inputs=(T(inner),T(.005),T(2e-8),one(T),2T(pi)*im*T(f))
                actual=selected(inputs...)(Val(kind))
                # The reference receives the exact represented binary arguments.
                # Permeability and pi remain independently evaluated constants.
                represented=map((256,512)) do bits
                    setprecision(BigFloat,bits) do
                        ri,ro,rho,mur,s=map(x->x isa Complex ? Complex{BigFloat}(x) : BigFloat(x),inputs)
                        RadialControl.surfaces(ri,ro,rho,4big(pi)*big"1e-7"*mur,s)
                    end
                end
                represented_expected=getproperty(last(represented),kind)
                represented_bound=last(represented).bound+
                    abs(represented_expected-getproperty(first(represented),kind))
                @test actual isa Complex{T}
                for component in (real,imag)
                    intended=component(expected)
                    target=component(represented_expected)
                    evaluation=BigFloat(component(actual))-target
                    representation=target-intended
                    observed=BigFloat(component(actual))-intended
                    push!(diagnostics,Dict("type"=>string(T),"kind"=>string(kind),
                        "inner"=>inner,"frequency"=>f,"component"=>string(component),
                        "observed"=>string(observed),"evaluation"=>string(evaluation),
                        "input_representation"=>string(representation),
                        "reference_bound"=>string(represented_bound),"target"=>string(1e-6abs(target))))
                    if iszero(target)
                        @test iszero(component(actual))
                    else
                        budget=1e-6*abs(target)
                        @test represented_bound<=budget/4
                        @test abs(component(actual)-target)+represented_bound<=budget
                    end
                end
            end
        end
    end
    directory=mktempdir(;prefix="lcm-radial-input-resolution-",cleanup=false)
    open(joinpath(directory,"decomposition.toml"),"w") do io
        TOML.print(io,Dict("cases"=>diagnostics))
    end
    println("Radial represented-input evidence: ",directory)
    # Derive the approach to DC from the same independently bounded radial
    # solution; no exact DC equality is imposed at a positive frequency.
    previous=Ref(Inf)
    for f in (1.0,.1,.01)
        reference=setprecision(BigFloat,256) do
            RadialControl.surfaces(big"0",big".005",big"2e-8",4big(pi)*big"1e-7",2big(pi)*im*f)
        end
        value=selected(0.0,.005,2e-8,1.0,2pi*im*f)(Val(:outer))
        rdc=2e-8/(pi*.005^2); ldc=4pi*1e-7/(8pi)
        remainder=abs(real(reference.outer)-rdc)+abs(imag(reference.outer)/(2pi*f)-ldc)
        @test remainder<previous[]
        @test abs(real(value)-rdc)<=abs(real(reference.outer)-rdc)+reference.bound+1e-6rdc
        @test abs(imag(value)/(2pi*f)-ldc)<=abs(imag(reference.outer)/(2pi*f)-ldc)+reference.bound/(2pi*f)+1e-6ldc
        previous[]=Float64(remainder)
    end
end
