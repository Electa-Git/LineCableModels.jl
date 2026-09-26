@testitem "ModalAnalysis / fixed nominal branch retains first-order scalar dependencies" tags=[:unit] begin
    using Measurements
    import Measurements: derivative

    x=measurement(2.0,0.1)
    phase=LineParameters(reshape([complex(x,zero(x))],1,1,1),
        reshape([3.0+0im],1,1,1),[50.0])
    modal=compute(ModalAnalysisProblem(phase),ModalAnalysisFormulation())
    root=real(gamma(modal)[1,1])
    @test nominal(root)≈sqrt(6.0)
    @test derivative(root,x)≈3/(2sqrt(6.0))
    characteristic=real(Zc(modal)[1,1])
    @test derivative(characteristic,x)≈1/(2sqrt(6.0))
    segment=PropagationParameters(modal;line_length=10.0)
    factor=real(H(segment)[1,1])
    @test derivative(factor,x)≈-10exp(-10sqrt(6.0))*3/(2sqrt(6.0))
    @test uncertainty(factor)>0
end

@testitem "ModalAnalysis / undefined first-order magnitude is unavailable" tags=[:unit] begin
    using Measurements
    import LineCableModels.Engine as E
    x=measurement(1.0,0.1)
    q1=complex(x,zero(x))
    q2=complex(2-x,zero(x))
    zero_q=zero(q1)
    coefficients=reshape([q1,zero_q,zero_q,q2],2,2,1)
    basis_matrix=ComplexF64[1 1;1 -1]./sqrt(2)
    maps=ModalOperators(reshape(copy(basis_matrix),2,2,1),
        reshape(copy(basis_matrix),2,2,1))
    modal=LineParameters(E.ModalDomain(maps,reshape([q1,q2],2,1)),
        SeriesImpedance(coefficients),ShuntAdmittance(copy(coefficients)),[50.0],
        ComputationDetails())
    segment=PropagationParameters(modal;line_length=1.0)
    voltage_H=Base.Fix2(H,(domain=PhaseDomain,field=:voltage))
    observed=ObservedResult(segment,((voltage_H,abs,1,2,:),);complete_pairs=true)
    product=LineCableModels.Grammar.observation_product(observed,(voltage_H,abs,1,2,:))
    @test product.values[1]===missing
    @test product.available[1]==false
    @test product.missing_reason[1]===:undefined_first_order_magnitude
end
