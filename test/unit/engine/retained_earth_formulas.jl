# Scientific real-axis controls survive withdrawal of the built-in author
# implementations. They exercise the generic callable integration service, not
# a substitute formula implementation or regenerated reference expectation.
@testitem "Engine / callable integration / independent Wise real-axis control" tags=[:unit] begin
    using QuadGK
    const E=LineCableModels.Engine
    controls=(; rtol = 1e-9, atol = 0.0, maxevals = 10^6)
    epsilon=8.8541878128e-12 .* [1.0, 10.0];
    mu=fill(4pi*1e-7, 2);
    rho=[Inf, 100.0]
    for f in (50.0, 10000.0), self in (true, false)

        heights=self ? (1.0, 1.0) : (1.0, 1.5)
        y=self ? 0.0 : 0.4;
        H=sum(heights)
        distance=self ? 0.005 : hypot(y, heights[1]-heights[2])
        references=[setprecision(BigFloat, bits) do
                        s=2big(pi)*im*f
                        e0=big"8.8541878128e-12";
                        u0=4big(pi)*big"1e-7"
                        g0=s^2*u0*e0;
                        g1=s*u0*(big".01"+s*10*e0)
                        ratio=g1/g0;
                        scale=sqrt(g1-g0)
                        cutoff=big(40)/H
                        kernel=lambda->exp(-H*lambda)*cos(y*lambda)/(ratio*lambda+sqrt(lambda^2+g1-g0))
                        integral,
                        estimate=quadgk(kernel, big"0", abs(scale/ratio),
                            abs(scale), inv(big(H)), cutoff;
                            rtol = big"1e-12", maxevals = 10^6)
                        # Re(ratio)=10 and Re(sqrt(lambda²+g1-g0))>0, hence
                        # |denominator|>=10lambda: this bounds the omitted real-axis tail.
                        tail=exp(-H*cutoff)/(10H*cutoff)
                        value=(log(hypot(big(y), big(H))/big(distance))+2integral)/(2big(pi)*e0)
                        (value = value, bound = (estimate+tail)/(big(pi)*e0))
                    end
                    for bits in (128, 256, 512)]
        expected=last(references)
        uncertainty=expected.bound+abs(expected.value-references[2].value)
        s=2pi*im*f
        g0=s^2*mu[1]*epsilon[1]
        g1=s*mu[2]*(inv(rho[2])+s*epsilon[2])
        ratio=g1/g0
        scale=sqrt(g1-g0)
        kernel=E.SpectralIntegral(lambda ->
            exp(-H*lambda)*cos(y*lambda)/(ratio*lambda+sqrt(lambda^2+g1-g0)))
        for numerical in (nothing, E.integration_workspace(Float64, ComplexF64))
            value, estimate=E.integrate(Val(:quad), kernel, controls, numerical;
                points=(abs(scale/ratio), abs(scale), inv(H)))
            actual=(log(hypot(y,H)/distance)+2value)/(2pi*epsilon[1])
            for component in (real, imag)
                budget=1e-6*abs(component(expected.value))
                @test uncertainty<=budget/4
                @test abs(component(actual-expected.value))+uncertainty<=budget
            end
        end
    end
end

@testitem "Engine / callable integration / independent Carson real-axis control" tags=[:unit] begin
    using QuadGK
    const E=LineCableModels.Engine
    controls=(; rtol = 1e-9, atol = 0.0, maxevals = 10^6)
    rho=[Inf, 100.0];
    epsilon=8.8541878128e-12 .* [1, 10];
    mu=fill(4pi*1e-7, 2)
    for f in (50.0, 10000.0), self in (true, false)

        s=2pi*im*f;
        heights=self ? (1.0, 1.0) : (1.0, 1.5)
        distance=self ? 0.005 : hypot(0.4, 0.5)
        y=self ? 0.0 : 0.4;
        h=sum(heights);
        image_distance=hypot(h, y)
        refs=[setprecision(BigFloat, bits) do
                  sb=2big(pi)*im*f;
                  mub=4big(pi)*big"1e-7"
                  k2=sb*mub/big"100";
                  hb=BigFloat(h);
                  yb=BigFloat(y)
                  cutoff=big"40"/hb
                  kernel(t)=exp(-hb*t)*cos(yb*t)/(t+sqrt(t^2+k2))
                  value,
                  estimate=quadgk(kernel, big"0", abs(sqrt(k2)), inv(hb), cutoff;
                      rtol = bits==128 ? big"1e-8" : bits==256 ? big"1e-10" : big"1e-12", maxevals = 10^6)
                  # Re sqrt(t²+i a)>0, so |denominator|>=t.
                  tail=exp(-hb*cutoff)/(hb*cutoff)
                  expected=sb*mub/(2big(pi))*(log(hypot(hb, yb)/BigFloat(distance))+2value)
                  (value = expected, bound = abs(sb*mub/big(pi))*(estimate+tail))
              end
              for bits in (128, 256, 512)]
        expected=last(refs).value
        uncertainty=last(refs).bound+abs(last(refs).value-refs[2].value)
        k2=s*mu[2]/rho[2]
        kernel=E.SpectralIntegral(lambda ->
            exp(-h*lambda)*cos(y*lambda)/(lambda+sqrt(lambda^2+k2)))
        for numerical in (nothing, E.integration_workspace(Float64, ComplexF64))
            value, estimate=E.integrate(Val(:quad), kernel, controls, numerical;
                points=(abs(sqrt(k2)), inv(h)))
            actual=s*mu[1]/(2pi)*(log(image_distance/distance)+2value)
            for component in (real, imag)
                budget=1e-6*abs(component(expected))
                @test uncertainty<=budget/4
                @test abs(component(actual)-component(expected))+uncertainty<=budget
            end
        end
    end
end
