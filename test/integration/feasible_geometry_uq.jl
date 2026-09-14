@testitem "UQ / current joint radial scale and spacing / actual line aggregation and derivatives" tags=[:integration,:extension] begin
    using Measurements,Random,Statistics,TOML
    function problem(scale,spacing)
        scale>0 || throw(DomainError(scale,"radial dimensions must be positive"))
        conductor=Material(kind=:conductor,rho=2e-8,eps_r=1.0,mu_r=1.0,T0=20.0,alpha=.004)
        dielectric=Material(kind=:insulator,rho=1e8,eps_r=3.0)
        designs=[build(CableDesign,"joint-wire-$i",Stack(
            terminal(:core,Region(:metal,Disk(.005scale),conductor)),
            Region(:insulation,Shell(.005scale),dielectric))) for i in 1:2]
        system=build(LineCableSystem,designs,[(0.0,-1.0),(spacing,-1.3)];
            connections=[Dict(:core=>1),Dict(:core=>2)])
        return LineParametersProblem(system;frequencies=[50.0,1000.0],
            earth_props=homogeneous(rho=100.0,eps_r=10.0),temperature=20.0)
    end
    # Uniform inputs have half-width sqrt(3)*sigma. Even at both worst-case
    # bounds the circular surfaces are separated before clearance adjustment.
    @test hypot(.2-sqrt(3)*.002,.3)-2*.01*(1.1+sqrt(3)*.01)>0
    calls=Tuple{Float64,Float64}[]
    log_inputs=Ref(false)
    builder=(scale,spacing)->begin
        log_inputs[] && !(scale isa Measurement) && !(spacing isa Measurement) && push!(calls,(scale,spacing))
        problem(scale,spacing)
    end
    space=Gridspace{LineParametersProblem}(builder,
        (Grid((1.0,1.1),AbsoluteError(.01)),Grid(.2,AbsoluteError(.002))))
    inner=Formulation(options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    parametric=ParametricProblem(space)
    linear=compute(parametric,LinearError(inner))
    @test length(linear.values)==2
    calibration=get(ENV,"LINECABLEMODELS_VALIDATION_PHASE","final")=="calibration"
    N=calibration ? 16 : 128
    seed=calibration ? 103 : 2029
    method=MonteCarlo(inner;trials=N,seed,distribution=:uniform,
        return_samples=true,retain_details=true)
    log_inputs[]=true
    sampled=compute(parametric,method)
    log_inputs[]=false
    @test sampled.trial_counts==[N,N]
    @test all(isempty,sampled.details.failures)
    @test length(unique(sampled.point_seeds))==2
    # The builder log contains the inputs that reached actual scalar compute.
    # Reconstruct the exact retained draws through the current sampler, then
    # compare the resulting channel arrays and independently aggregate them.
    recorded=copy(calls)
    for (index,point) in enumerate(LineCableModels.points(space))
        rng=Xoshiro(sampled.point_seeds[index])
        expected=NamedTuple[]
        realized=Tuple{Float64,Float64}[]
        for trial in 1:N
            arguments=LineCableModels.realize_arguments(rng,point,:uniform)
            push!(realized,Tuple(arguments))
            value=compute(problem(arguments...),inner)
            push!(expected,(R=R(value),L=L(value),C=C(value),G=G(value)))
        end
        @test all(pair->count(==(pair),recorded)==1,realized)
        for channel in (:R,:L,:C,:G)
            values=cat((getproperty(value,channel) for value in expected)...;dims=4)
            @test getproperty(sampled.sample_values[index],channel)==values
            average=dropdims(sum(values;dims=4)./N;dims=4)
            @test getproperty(LineCableModels,channel)(sampled.values[index]) ≈ average rtol=1e-10 atol=0
        end
    end
    # Direct Measurements derivatives are compared to independently evaluated
    # central differences. The fixed refinement sequence cannot select a lucky h.
    directory=mktempdir(;prefix="lcm-uq-derivatives-provisional-",cleanup=false)
    println("UQ derivative evidence: ",directory);flush(stdout)
    diagnostics=Dict{String,Any}[]
    for nominal_scale in (1.0,1.1), coordinate in 1:2
        nominal=[nominal_scale,.2]
        # A unit-uncertainty spacing probe activates the clearance reserve and
        # changes the nominal geometry. Use the prescribed interior support.
        variable=measurement(nominal[coordinate],coordinate==1 ? .01 : .002)
        inputs=coordinate==1 ? (variable,nominal[2]) : (nominal[1],variable)
        differentiated_problem=problem(inputs...)
        @test [(LineCableModels.nominal(p.x),LineCableModels.nominal(p.y))
            for p in differentiated_problem.system.positions] == [(0.0,-1.0),(.2,-1.3)]
        differentiated=compute(differentiated_problem,inner)
        previous=nothing;resolved=0
        for relative in (.01,.005,.0025,.00125)
            step=relative*nominal[coordinate]
            plus=copy(nominal);minus=copy(nominal)
            plus[coordinate]+=step;minus[coordinate]-=step
            hi=compute(problem(plus...),inner);lo=compute(problem(minus...),inner)
            derivatives=map((R,L,C,G)) do quantity
                (quantity(hi).-quantity(lo))./(2step)
            end
            if previous!==nothing
                all_resolved=true
                for (quantity,current,coarse) in zip((R,L,C,G),derivatives,previous)
                    extrapolated=(4current.-coarse)./3
                    actual=map(x->x isa Measurement ? Measurements.derivative(x,variable) : 0.0,quantity(differentiated))
                    budget=.001abs.(extrapolated)
                    uncertainty=abs.(extrapolated.-current)
                    all_resolved &= all(uncertainty .<= budget./4)
                    all_resolved &= all(abs.(actual.-extrapolated).+uncertainty .<= budget)
                    for entry in eachindex(actual)
                        push!(diagnostics,Dict("scale"=>nominal_scale,"coordinate"=>coordinate,
                            "relative_step"=>relative,"quantity"=>string(quantity),"entry"=>entry,
                            "actual"=>actual[entry],"reference"=>extrapolated[entry],
                            "uncertainty"=>uncertainty[entry],"budget"=>budget[entry],
                            "reference_resolved"=>uncertainty[entry]<=budget[entry]/4,
                            "comparison_passed"=>abs(actual[entry]-extrapolated[entry])+uncertainty[entry]<=budget[entry]))
                    end
                end
                resolved=all_resolved ? resolved+1 : 0
            end
            previous=derivatives
        end
        open(joinpath(directory,"derivatives.toml"),"w") do io
            TOML.print(io,Dict("cases"=>diagnostics))
        end
        @test resolved>=2
    end
    @test_throws DomainError problem(-.1,.2)
end
