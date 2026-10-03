using LineCableModels

conductor = Material(
    kind = :conductor,
    rho = eps(Float64),
    eps_r = 1.0,
    mu_r = 1.0,
    T0 = 20.0,
    alpha = 0.0
)

insulator = Material(
    kind = :insulator,
    rho = 1.97e14,
    eps_r = Grid((2.3,), 1.0),
    mu_r = 1.0,
    T0 = 20.0,
    alpha = 0.0
)

design = @cable "two_bare_wires" begin
    @terminal :core begin
        core(conductor; r=0.0425, tag=:core_metal)
    end
    insulation(insulator; t=1.0e-3, tag=:core_insulation)
end

system = @system "two_bare_wires" line_length=1.0 combine=:zip begin
    @at design (0.0, -1.0) core=1
    @at design (1.0, -1.0) core=2
end

earth = @earth begin
    layer(rho=Grid((10.0, 100.0)), eps_r=1.0, mu_r=1.0)
end

problems = LineParametersProblem(
    system,
    earth;
    temperature = 20.0,
    frequencies = collect(10.0 .^ range(0, stop = 6, length = 101))
)

display(problems)
