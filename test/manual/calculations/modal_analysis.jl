# Manual study of the nine underground coaxial cables in Case 3, Fig. 17 and
# Table IV. Run with LineCableModels and GLMakie in the active REPL project:
#     include("test/manual/calculations/modal_analysis.jl")
# Completed values and interactive figures remain available in the REPL.
# This file does not activate a project or write images to disk.

using Revise
using LineCableModels, GLMakie, LinearAlgebra

frequency_grid = 10.0 .^ range(-1, 8; length = 321)
source_length = 2_000.0
phase_options = (output_basis = :pul, verbosity = (default = 0,))
modal_formulation = ModalAnalysisFormulation(:chrysochos2014)

# Table IV radii are measured from the cable center. The table does not give
# dielectric bulk resistivity or conductor temperature coefficients; infinite
# dielectric resistivity and zero temperature coefficient avoid adding losses
# beyond its stated insulation loss factor and metal resistivity.
core_radius = 12.38e-3
sheath_inner_radius = 24.69e-3
sheath_outer_radius = 25.36e-3
insulation_outer_radius = 27.358e-3
metal = Material(:conductor, 1.7e-8)
core_insulator = Material(:insulator, Inf, 2.3; tan_delta = 0.0005)
sheath_insulator = Material(:insulator, Inf, 2.48; tan_delta = 0.0005)

design = @cable "case-3-coaxial" begin
    @terminal :core begin
        core(metal; r = core_radius)
        insulation(core_insulator; t = sheath_inner_radius - core_radius)
    end
    @terminal :sheath begin
        sheath(metal; t = sheath_outer_radius - sheath_inner_radius)
    end
    jacket(sheath_insulator; t = insulation_outer_radius - sheath_outer_radius)
end

# Fig. 17: 1.1335 m to the upper cable centers, 0.1905 m across columns,
# and 0.19035 m down rows. Each core and sheath has its own terminal.
centers = [(column * 0.1905, -1.1335 - row * 0.19035)
           for row in 0:2 for column in -1:1]
connections = [Dict(:core => 2i - 1, :sheath => 2i)
               for i in eachindex(centers)]
earth = homogeneous(rho = 100.0)
system = build(LineCableSystem, fill(design, length(centers)),
    [Pose2(x, y) for (x, y) in centers]; connections,
    system_id = "underground-case-3", line_length = source_length,
    environment = earth)
problem = LineParametersProblem(system; temperature = 20.0, earth_props = earth,
    frequencies = collect(frequency_grid))
formulation = Formulation(earth_impedance = :unified, earth_admittance = :unified,
    shunt_model = :coaxial, insulation_admittance = :lossy;
    options = (reduce_bundle = false, kron_reduction = false,
        ideal_transposition = false))

@assert length(centers) == 9
@assert length(system.terminal_order) == 18
@assert sort(system.connection_order) == collect(1:18)
display(design)
display(system)
display(problem)
geometry_plot = preview(system; earth_model = earth, backend = :gl,
    display_plot = true, open_export = false)

# Two public calculations, followed by the finite line from Table IV.
phase = @time compute(problem, formulation; options = phase_options)
modal = @time compute(ModalAnalysisProblem(phase), modal_formulation)
restored_phase = transform(PhaseDomain, modal)
line = PropagationParameters(modal)
display(phase)
display(modal)
display(line)

@assert size(Z(phase)) == size(Y(phase)) == (18, 18, length(frequency_grid))
@assert size(Tv(modal)) == size(Ti(modal)) == (18, 18, length(frequency_grid))
@assert size(gamma(modal)) == (18, length(frequency_grid))
@assert line_length(line) == source_length

diagnostics = details(modal).data.modal.diagnostics
# println("Modal diagnostics: ", diagnostics)
println("Inverse-coordinate residuals: Z=", norm(Z(restored_phase)-Z(phase)),
    ", Y=", norm(Y(restored_phase)-Y(phase)))

roots = gamma(line)
voltage_basis = Tv(line)
current_basis = Ti(line)
modal_Zc = Zc(line)
modal_Yc = Yc(line)
modal_H = H(line)
phase_Zc = Zc(line, PhaseDomain)
phase_Yc = Yc(line, PhaseDomain)
phase_voltage_H = H(line, PhaseDomain; field = :voltage)
phase_current_H = H(line, PhaseDomain; field = :current)
@assert size(phase_Zc) == size(phase_Yc) == size(phase_voltage_H) ==
        size(phase_current_H) == (18, 18, length(frequency_grid))

# These bound selectors are ordinary public observation requests.
observed = observables(line)
display(report(observed))

# Every retained mode is a curve for scalar modal quantities. For Tv and Ti,
# every phase-row × mode-column coefficient is a curve, including off-diagonals.
# The observation owner supplies physical labels, units and original indices.
# real/imag of Zc and Yc are R꜀/X꜀ and G꜀/B꜀, respectively.
plot_options = (overlay = :coordinates, backend = :gl,
    display_plot = true, open_export = false)
plots = (
    cartesian = LineCableModels.plot(observed;
        ydata = ((Zc, real), (Zc, imag), (Yc, real), (Yc, imag)), plot_options...),
    polar = LineCableModels.plot(observed;
        ydata = ((Zc, abs), (Zc, angle), (Yc, abs), (Yc, angle)), plot_options...),
    propagation = LineCableModels.plot(observed;
        ydata = (alpha, beta, velocity), plot_options...)
)

transform_matrices = observables(line,
    ((Tv, abs), (Tv, angle), (Ti, abs), (Ti, angle)))

display(report(transform_matrices))
LineCableModels.plot(transform_matrices; plot_options...)

# Numerical quality is for inspection, not an acceptance threshold.
voltage_condition = [cond(@view voltage_basis[:, :, k]) for k in eachindex(frequency_grid)]
current_condition = [cond(@view current_basis[:, :, k]) for k in eachindex(frequency_grid)]
eigen_residual = [begin
                      assessed = filter(x -> x !== nothing, column)
                      isempty(assessed) ? NaN : maximum(assessed)
                  end
                  for column in eachcol(diagnostics.eigen_residual)]
diagnostic_figure = Figure(size = (1200, 800))
z_axis = Axis(diagnostic_figure[1, 1]; title = "Z modal coupling", xscale = log10,
    xlabel = "Frequency (Hz)")
y_axis = Axis(diagnostic_figure[1, 2]; title = "Y modal coupling", xscale = log10,
    xlabel = "Frequency (Hz)")
residual_axis = Axis(diagnostic_figure[2, 1]; title = "Largest returned eigen-residual",
    xscale = log10, xlabel = "Frequency (Hz)")
condition_axis = Axis(diagnostic_figure[2, 2]; title = "Basis conditioning",
    xscale = log10, xlabel = "Frequency (Hz)")
lines!(z_axis, frequency_grid, diagnostics.z_coupling)
lines!(y_axis, frequency_grid, diagnostics.y_coupling)
lines!(residual_axis, frequency_grid, eigen_residual)
lines!(condition_axis, frequency_grid, voltage_condition; label = "Tv")
lines!(condition_axis, frequency_grid, current_condition; label = "Ti")
axislegend(condition_axis)
display(diagnostic_figure)

println("Case 3: 9 cables, 18 modes, ", length(frequency_grid),
    " samples, 2 km cable length")
println("Inspect phase, modal, line, observed, table_report, plots.cartesian, ",
    "plots.polar, plots.propagation, plots.transformations, and diagnostic_figure.")
nothing
