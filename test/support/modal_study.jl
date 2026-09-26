@testmodule ModalStudyFixtures begin
    using LineCableModels

    # Explicit toy, nondispersive properties. These are not measured or certified
    # broadband cable data. The armor is an annular equivalent, not discrete wires.
    function material(kind,rho,eps_r,mu_r,alpha;tan_delta=0.0)
        return Material(kind,rho,eps_r,mu_r,20.0,alpha;
            rho_thermal=0.0,theta_max=90.0,tan_delta,sigma_solar=0.0)
    end

    function armored_design()
        copper=material(:conductor,1.724e-8,1.0,1.0,0.00393)
        steel=material(:conductor,1.5e-7,1.0,50.0,0.006)
        semicon=material(:semicon,1.0,100.0,1.0,0.0)
        insulation_material=material(:insulator,1e14,2.3,1.0,0.0)
        bedding_material=material(:insulator,1e12,3.0,1.0,0.0)
        jacket_material=material(:insulator,1e12,3.0,1.0,0.0)
        return @cable "toy-armored-single-core" begin
            @terminal :core begin
                core(copper;r=0.012)
                screen(semicon;t=0.0005)
                insulation(insulation_material;t=0.006)
                screen(semicon;t=0.0005)
            end
            @terminal :screen begin
                sheath(copper;t=0.0008)
            end
            bedding(bedding_material;t=0.002)
            @terminal :armor begin
                sheath(steel;t=0.0015)
            end
            jacket(jacket_material;t=0.003)
        end
    end

    function study_problem(layout::Symbol; frequencies=10.0 .^ range(0,7;length=281))
        design=armored_design()
        centers=if layout===:horizontal
            ((-0.15,-1.2),(0.0,-1.2),(0.15,-1.2))
        elseif layout===:trefoil
            ((-0.075,-1.2),(0.075,-1.2),(0.0,-1.2+sqrt(3)*0.15/2))
        else
            throw(ArgumentError("layout must be :trefoil or :horizontal"))
        end
        connections=[Dict(:core=>3i-2,:screen=>3i-1,:armor=>3i) for i in 1:3]
        system=build(LineCableSystem,[design,design,design],
            [Pose2(x,y) for (x,y) in centers];
            connections,system_id="toy-armored-$(layout)",line_length=600.0,
            environment=homogeneous(rho=100.0,eps_r=10.0,mu_r=1.0))
        problem=LineParametersProblem(system;temperature=20.0,
            earth_props=homogeneous(rho=100.0,eps_r=10.0,mu_r=1.0),
            frequencies=collect(frequencies))
        return (;design,system,problem,connections,centers)
    end

    function study_formulation()
        return Formulation(earth_impedance=:unified,earth_admittance=:unified,
            shunt_model=:coaxial,insulation_admittance=:lossy,
            semicon_admittance=:lossy;
            options=(reduce_bundle=false,kron_reduction=false,
                ideal_transposition=false))
    end
end
