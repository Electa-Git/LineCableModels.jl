@testitem "Engine / internal impedance / radial diffusion implementation" tags=[:unit] begin
    include(joinpath(pkgdir(LineCableModels),"test/support/radial_control.jl"))
    const II=LineCableModels.Engine.InternalImpedance
    selected=II.Formula(:default)
    for inner in (0.0,.003), f in (10.0,1000.0), T in (Float64,BigFloat)
        inputs=(T(inner),T(.005),T(2e-8),one(T),2T(pi)*im*T(f))
        # One reference per represented input serves all three coefficients.
        references=map((256,512)) do bits
            setprecision(BigFloat,bits) do
                ri,ro,rho,mur,s=map(x->x isa Complex ? Complex{BigFloat}(x) : BigFloat(x),inputs)
                RadialControl.surfaces(ri,ro,rho,4big(pi)*big"1e-7"*mur,s)
            end
        end
        evaluator=selected(inputs...)
        @testset "$T $kind ri=$inner f=$f" for kind in (:inner,:outer,:transfer)
            actual=evaluator(Val(kind))
            expected=getproperty(last(references),kind)
            bound=last(references).bound+abs(expected-getproperty(first(references),kind))
            @test actual isa Complex{T}
            for component in (real,imag)
                target=component(expected)
                if iszero(target)
                    @test iszero(component(actual))
                else
                    budget=1e-6*abs(target)
                    @test bound<=budget/4
                    @test abs(component(actual)-target)+bound<=budget
                end
            end
        end
    end
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


@testitem "Engine / internal impedance / Float32 input compatibility" tags=[:unit] begin
    const II=LineCableModels.Engine.InternalImpedance
    selected=II.Formula(:default)
    # Float32 is an input compatibility contract, with no extra accuracy target.
    for inner in (0f0,.003f0)
        evaluator=selected(inner,.005f0,2f-8,1f0,20Float32(pi)*im)
        for kind in (:inner,:outer,:transfer)
            value=evaluator(Val(kind))
            @test value isa ComplexF32
            @test isfinite(value)
        end
    end
end
