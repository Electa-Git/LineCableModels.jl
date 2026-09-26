@testitem "ModalAnalysis / finite observable components and vector retention" tags=[:unit] begin
    using LinearAlgebra
    using DataFrames
    using XLSX
    using Measurements
    import LineCableModels.Engine as E
    import LineCableModels.Grammar as G
    import LineCableModels.Units as U

    function finite_modal(::Type{T}; roots=Complex{T}[1+2im 3+4im 5+6im;
            2+3im 4+5im 6+7im] .* T(1e-3)) where {T<:AbstractFloat}
        f=T[10,20,30]
        maps=cat([Matrix{Complex{T}}(I,2,2) for _ in f]...;dims=3)
        z=zeros(eltype(roots),2,2,3)
        y=copy(z)
        for k in eachindex(f), mode in 1:2
            z[mode,mode,k]=roots[mode,k]
            y[mode,mode,k]=roots[mode,k]
        end
        modal=LineParameters(E.ModalDomain(ModalOperators(maps,maps),roots),
            SeriesImpedance(z),ShuntAdmittance(y),f,
            ComputationDetails((inputs=(system=(line_length=T(20),),),
                phase_coordinates=["core","sheath"])))
        return PropagationParameters(modal)
    end

    for T in (Float32,Float64)
        segment=finite_modal(T)
        @test alpha(segment)==real.(gamma(segment))
        @test beta(segment)==imag.(gamma(segment))
        @test velocity(segment)≈T(2)*T(π).*reshape(frequencies(segment),1,:)./beta(segment)
        @test eltype(velocity(segment))==T
        @test alpha(PropagationParameters(segment;line_length=T(5)))==alpha(segment)
        @test beta(PropagationParameters(segment;line_length=T(5)))==beta(segment)
        @test H(PropagationParameters(segment;line_length=T(5)))!=H(segment)
    end
    setprecision(256) do
        segment=finite_modal(BigFloat)
        @test eltype(velocity(segment))==BigFloat
        @test precision(velocity(segment)[1,1])>=256
        @test velocity(segment)[1,1]≈2*big(π)*frequencies(segment)[1]/beta(segment)[1,1]
    end
    @test U.quantity(alpha)==U.quantity(gamma,real)
    @test U.quantity(beta)==U.quantity(gamma,imag)
    @test U.quantity(gamma,real)!=U.quantity(gamma,imag)
    @test U.label(U.display_unit(U.quantity(alpha)))=="Np/km"
    @test U.label(U.display_unit(U.quantity(beta)))=="rad/km"
    @test U.scale_factor(U.native_unit(U.quantity(alpha)),U.display_unit(U.quantity(alpha)))==1000
    @test U.scale_factor(U.native_unit(U.quantity(beta)),U.display_unit(U.quantity(beta)))==1000
    @test U.label(U.display_unit(U.quantity(velocity)))=="m/s"
    @test U.label(U.quantity(Zc,abs))=="Characteristic impedance magnitude"
    @test U.label(U.quantity(Tv,angle))=="Modal-to-phase voltage transformation angle"

    segment=finite_modal(Float64)
    @test observe(segment,alpha)==alpha(segment)
    @test observe(segment,beta)==beta(segment)
    @test observe(segment,velocity)==velocity(segment)
    @test observe(segment,gamma,real)==alpha(segment)
    @test observe(segment,gamma,[2,1],[3,1])==gamma(segment)[[2,1],[3,1]]
    @test observe(segment,gamma,real,2,1:3)==real.(gamma(segment)[2,1:3])
    @test observe(segment.parameters,Tv,1,2,1:3)==Tv(segment)[1,2,1:3]
    phase_zc=Base.Fix2(Zc,(domain=PhaseDomain,))
    @test observe(segment,phase_zc,real,1,1,1:2)==
        real.(Zc(segment,PhaseDomain)[1,1,1:2])
    @test_throws ArgumentError observe(segment,
        Base.Fix2(Zc,(domain=:invalid,)),1,1,1)
    @test G.observation_requests(segment,(gamma,)).retained==
        ((gamma,real,Colon(),Colon()),(gamma,imag,Colon(),Colon()))
    @test_throws ArgumentError ObservedResult(segment,((gamma,real),))
    @test length(ObservedResult(segment,((gamma,real),);complete_pairs=true).quantities)==2
    @test G.observation_requests(segment,(alpha,beta)).retained==
        ((gamma,real,Colon(),Colon()),(gamma,imag,Colon(),Colon()))
    @test_throws ArgumentError ObservedResult(segment,(gamma,alpha,beta))
    observed=ObservedResult(segment)
    @test length(observed.quantities)==13
    @test all(product -> product.coordinates.kind===:vector,
        filter(q -> q.family in (:gamma,:velocity,:Zc,:Yc,:H),observed.quantities))
    @test G.observation_product(observed,alpha)===G.observation_product(observed,(gamma,real))
    @test G.observation_product(observed,beta)===G.observation_product(observed,(gamma,imag))
    @test observe(observed,alpha)==G.observation_product(observed,(gamma,real)).values
    @test G.observation_requests(observed,(alpha,)).displayed==((gamma,real),)
    @test_throws ArgumentError G.observation_requests(observed,(alpha,(gamma,real)))
    @test_throws ArgumentError report(observed;values=(alpha,(gamma,real)))
    @test G.observation_product(observed,alpha).coordinates.positions==[1,2]
    @test size(G.observation_product(observed,alpha).values)==(2,3)
    selection=@observe alpha[[2,1],[3,1]]
    subset=G.observation_product(observed,selection)
    @test subset.coordinates.positions==[2,1]
    @test subset.coordinates.samples==[3,1]
    @test subset.values≈alpha(segment)[[2,1],[3,1]].*1000
    @test G.observation_product(ObservedResult(segment,(selection,);complete_pairs=true),selection).values==subset.values
    scalar_request=@observe alpha[2,2]
    @test G.observation_product(observed,scalar_request).values≈alpha(segment)[2,2]*1000
    @test G.observation_product(observed,scalar_request).coordinates.positions==[2]
    range_request=@observe alpha[1:2,2:3]
    @test G.observation_product(observed,range_request).values≈alpha(segment)[1:2,2:3].*1000
    converted=ObservedResult(observed,(alpha,);length_unit=:base)
    @test G.observation_product(converted,alpha).values≈alpha(segment)
    @test G.observation_product(converted,alpha).available==G.observation_product(observed,alpha).available
    native_alpha=U.native_unit(U.quantity(alpha))
    for overrides in ((alpha=native_alpha,),Dict(alpha=>native_alpha))
        acquired=ObservedResult(segment,(alpha,);complete_pairs=true,quantity_units=overrides)
        reexpressed=ObservedResult(observed,(alpha,);quantity_units=overrides)
        for result in (acquired,reexpressed)
            product=G.observation_product(result,alpha)
            @test product.unit==native_alpha
            @test product.values≈alpha(segment)
            @test product.available==G.observation_product(observed,alpha).available
        end
    end
    @test G.observation_product(ObservedResult(observed,(alpha,);
        quantity_units=(gamma=:kilo,alpha=native_alpha)),alpha).unit==native_alpha
    @test G.observation_product(ObservedResult(observed,(alpha,);
        quantity_units=Dict(alpha=>native_alpha,
            (gamma,real)=>U.display_unit(U.quantity(alpha)))),alpha).unit!=native_alpha
    @test G.observation_product(ObservedResult(observed,(alpha,);
        quantity_units=Dict{Any,Any}(:alpha=>U.display_unit(U.quantity(alpha)),
            alpha=>native_alpha)),alpha).unit==native_alpha
    indexed_alias=@observe alpha[2,:]
    indexed_gamma_real=@observe (gamma,real)[2,:]
    for result in (
            ObservedResult(segment,(indexed_alias,);complete_pairs=true,
                quantity_units=Dict(indexed_alias=>native_alpha)),
            ObservedResult(observed,(indexed_alias,);
                quantity_units=Dict(indexed_alias=>native_alpha)))
        product=G.observation_product(result,indexed_alias)
        @test product.unit==native_alpha
        @test product.values≈alpha(segment)[2,:]
        @test product.available==G.observation_product(observed,indexed_alias).available
    end
    @test G.observation_product(ObservedResult(observed,(indexed_alias,);
        quantity_units=Dict{Any,Any}(indexed_alias=>native_alpha,
            (gamma,real)=>U.display_unit(U.quantity(alpha)))),indexed_alias).unit==native_alpha
    @test G.observation_product(ObservedResult(observed,(indexed_alias,);
        quantity_units=Dict{Any,Any}(indexed_alias=>native_alpha,
            indexed_gamma_real=>U.display_unit(U.quantity(alpha)))),indexed_alias).unit!=native_alpha
    native_beta=U.native_unit(U.quantity(beta))
    @test G.observation_product(ObservedResult(observed,(beta,);
        quantity_units=(beta=native_beta,)),beta).values≈beta(segment)
    shorter=ObservedResult(segment[[1,3]],(alpha,);complete_pairs=true)
    aligned=G.observation_product([ObservedResult(segment,(alpha,);complete_pairs=true),shorter],alpha)
    @test aligned[1].coordinates.samples==[1,2,3]
    @test aligned[2].coordinates.samples==[1,2]
    @test aligned[2].coordinates.frequencies==[10.0,30.0]
    @test length(G.observation_product(observed,(Tv,real)).coordinates.rows)==2
    @test G.observation_product(observed,(Tv,real)).coordinates.column_domain===:ModalDomain

    zero_roots=copy(gamma(segment))
    zero_roots[1,2]=complex(real(zero_roots[1,2]),0.0)
    zero=finite_modal(Float64;roots=zero_roots)
    unavailable=G.observation_product(ObservedResult(zero,(velocity,)),velocity)
    @test ismissing(unavailable.values[1,2])
    @test unavailable.missing_reason[1,2]===:nonfinite_value
    @test !unavailable.available[1,2]

    shared=measurement(0.002,0.0001)
    uncertain_roots=map(root -> complex(shared,shared+imag(root)),gamma(segment))
    uncertain=finite_modal(Float64;roots=uncertain_roots)
    uncertain_observed=ObservedResult(uncertain,(alpha,beta,velocity))
    alpha_value=G.observation_product(uncertain_observed,alpha).values[1,1]
    beta_value=G.observation_product(uncertain_observed,beta).values[1,1]
    @test iszero(uncertainty(alpha_value-beta_value))
    @test uncertainty(G.observation_product(uncertain_observed,velocity).values[1,1])>0
    indexed_beta=@observe beta[2,:]
    uncertain_converted=ObservedResult(uncertain_observed,(indexed_alias,indexed_beta);
        quantity_units=Dict(indexed_alias=>native_alpha,indexed_beta=>native_beta))
    @test G.observation_product(uncertain_converted,indexed_alias).unit==native_alpha
    @test G.observation_product(uncertain_converted,indexed_beta).unit==native_beta
    @test iszero(uncertainty(
        G.observation_product(uncertain_converted,indexed_alias).values[1,1]-
        G.observation_product(uncertain_converted,indexed_beta).values[1,1]))

    artifact=report(TableReportDefinition((alpha,Tv)),observed)
    alpha_table=artifact[alpha]
    @test alpha_table===artifact[1,alpha]
    @test alpha_table===artifact[(gamma,real)]
    @test size(alpha_table)==(3,3)
    @test names(alpha_table)==["frequency","Mode 1","Mode 2"]
    tv_table=artifact[(Tv,real)]
    @test size(tv_table)==(3,5)
    @test all(name -> occursin("Conductor",name) && occursin("Mode",name),names(tv_table)[2:end])
    @test occursin("Attenuation constant",sprint(show,MIME"text/plain"(),artifact;context=:limit=>false))
    @test occursin("Mode 1",sprint(show,MIME"text/html"(),artifact))
    gamma_artifact=report(TableReportDefinition((gamma,)),observed)
    @test gamma_artifact[alpha]===gamma_artifact[(gamma,real)]
    @test gamma_artifact[beta]===gamma_artifact[(gamma,imag)]
    mktempdir() do directory
        for extension in ("json","jls")
            path=joinpath(directory,"modal.$extension")
            save(observed,path)
            restored=import_data(:observed,path)
            @test G.observation_product(restored,alpha).values==G.observation_product(observed,alpha).values
            @test observe(restored,beta)==G.observation_product(restored,(gamma,imag)).values
            @test report(TableReportDefinition((gamma,)),restored)[alpha] isa DataFrame
            @test G.observation_product(restored,(Tv,real)).coordinates.column_labels==["1","2"]
            uncertain_path=joinpath(directory,"uncertain.$extension")
            save(uncertain_observed,uncertain_path)
            uncertain_restored=import_data(:observed,uncertain_path)
            a=G.observation_product(uncertain_restored,alpha).values[1,1]
            b=G.observation_product(uncertain_restored,beta).values[1,1]
            @test iszero(uncertainty(a-b))
            setprecision(256) do
                precise=ObservedResult(finite_modal(BigFloat),(alpha,velocity);complete_pairs=true)
                precise_path=joinpath(directory,"precise.$extension")
                save(precise,precise_path)
                precise_restored=import_data(:observed,precise_path)
                @test precision(G.observation_product(precise_restored,velocity).values[1,1])>=256
            end
        end
        output=export_data(:xlsx,ObservedResult(segment,(alpha,);complete_pairs=true);
            file_name=joinpath(directory,"modal.xlsx"))
        @test all(isfile,output)
    end
end
