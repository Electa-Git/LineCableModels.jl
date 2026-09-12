@testitem "Measurements / Bessel derivatives at small complex arguments and scaling branches" tags=[:extension] begin
    using Measurements, SpecialFunctions
    for magnitude in (1e-8,1e-6,1e-4)
        x=measurement(magnitude,magnitude/100)
        y=measurement(magnitude,magnitude/100)
        z=complex(x,y)
        expected=-besselk(1,complex(magnitude,magnitude))
        result=besselk(0,z)
        @test complex(Measurements.derivative(real(result),x),
            Measurements.derivative(imag(result),x)) ≈ expected rtol=1e-12
        @test complex(Measurements.derivative(real(result),y),
            Measurements.derivative(imag(result),y)) ≈ im*expected rtol=1e-12
        # Stable small-argument series also check the regular and scaled I/J
        # slopes, which Float64 output subtraction cannot resolve reliably.
        z0=complex(magnitude,magnitude)
        for (kernel,constant,slope,scale,scale_x,scale_y) in (
                (besseli,1+z0^2/4+z0^4/64,z0/2+z0^3/16+z0^5/384,1.,0.,0.),
                (besselj,1-z0^2/4+z0^4/64,-z0/2+z0^3/16-z0^5/384,1.,0.,0.),
                (besselix,1+z0^2/4+z0^4/64,z0/2+z0^3/16+z0^5/384,exp(-magnitude),-1.,0.),
                (besseljx,1-z0^2/4+z0^4/64,-z0/2+z0^3/16-z0^5/384,exp(-magnitude),0.,-1.))
            output=kernel(0,z)
            for (input,derivative) in ((x,scale*(slope+scale_x*constant)),
                    (y,scale*(im*slope+scale_y*constant)))
                actual=complex(Measurements.derivative(real(output),input),
                    Measurements.derivative(imag(output),input))
                @test actual ≈ derivative rtol=1e-12
            end
        end
    end
    # Rebuilt relative-step Cartesian differences are independent of the adapter's
    # recurrence implementation. Exercise both signs of each scaling coordinate.
    kernels=(besseli,besselk,besselj,bessely,besselh,
        besselix,besselkx,besseljx,besselyx,besselhx)
    # At tiny arguments, regular I/J values are nearly constant and subtracting
    # two Float64 outputs cannot resolve their slopes. Use the exact K0 test
    # above for that regime; rebuilt differences below are well-scaled controls.
    for magnitude in (1.,20.), sx in (-1.,1.), sy in (-1.,1.),
            order in (0.,0.5,1.,2.), kernel in kernels
        nominal_z=magnitude*complex(1.25sx,0.75sy)
        q=measurement(0.,1e-4*magnitude)
        direction=1+2im
        z=complex(real(nominal_z)+q,imag(nominal_z)+2q)
        actual=kernel(order,z)
        h=1e-5*magnitude
        expected=(kernel(order,nominal_z+h*direction)-
            kernel(order,nominal_z-h*direction))/(2h)
        slope=complex(Measurements.derivative(real(actual),q),
            Measurements.derivative(imag(actual),q))
        @test slope ≈ expected rtol=2e-6 atol=1e-11*max(1.,abs(expected))
        @test Measurements.value(actual) ≈ kernel(order,nominal_z) rtol=1e-12
        @test Measurements.cov(real(actual),imag(actual)) ≈
            real(slope)*imag(slope)*Measurements.uncertainty(q)^2 rtol=1e-12
    end
end
