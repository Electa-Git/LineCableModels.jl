# Two bare copper wires, FEM against Unified earth return, in four placements:
# both in air, both in earth, and each mixed arrangement, over every wire radius
# and soil resistivity of run_two_bare_wires_fem.jl.
#
# Run from the repository root with:
#     julia -i --project=test test/manual/calculations/run_two_bare_wires_placements.jl
#
# Results are left in fem_results and unified_results, in the order of points.
# e_Z and e_Y are the scaled errors max|ΔZ|/max|Z_ii| and max|ΔY|/max|Y_ii| per
# point (rows) and frequency (columns).

using LineCableModels
using LinearAlgebra
using Printf
using Gmsh
using GLMakie

frequencies = 10.0 .^ range(-1.0, 8.0; length = 91)  # 0.1 Hz to 100 MHz, 10 per decade
copper = Material(MaterialsLibrary(add_defaults = true), :copper)

radius_values = (0.1e-2, 1.0e-2, 4.25e-2, 8.5e-2)
rho_values = (0.1, 1.0, 100.0, 1000.0)
# Wire heights [m], zipped by @system: one column per placement.
placements = ["both in air", "both in earth", "1 air, 2 earth", "1 earth, 2 air"]
y1 = Grid((1.0, -1.0, 1.0, -1.0))
y2 = Grid((1.0, -1.0, -1.0, 1.0))

wire = @cable "bare_wire" begin
    @terminal :core begin
        core(copper; r = Grid(radius_values), tag = :core_metal)
    end
end
# Both placements vary radius fastest, so the zip pairs equal radii.
system = @system "two_bare_wires" line_length=1.0 combine=:zip begin
    @at wire (0.0, y1) core=1
    @at wire (1.0, y2) core=2
end
earth = @earth begin
    layer(rho = Grid(rho_values), eps_r = 1.0, mu_r = 1.0)
end
problems = LineParametersProblem(system, earth; temperature = 20.0, frequencies)
# Radius varies fastest, then placement, then resistivity.
points = vec(["$p, r = $(100r) cm, ρ = $ρ Ω·m"
              for r in radius_values, p in placements, ρ in rho_values])

reductions = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false)
unified = Formulation(earth_impedance = :unified, earth_admittance = :unified,
    options = reductions)
fem = Formulation(:fem; options = (; physics = :helmholtz, reductions...))
fem_options = ComputationOptions(
    frequency_workers = 12,
    keep_run_directory = true,
    resume_run_directory = :latest,
    verbosity = (default = 1,),
    gmsh_verbosity = 0,
    getdp_verbosity = 0,
    mesh_size_factor = 0.5
)

unified_results = collect(compute(ParametricProblem(problems), Combinatorial(unified)))
fem_results = collect(compute(ParametricProblem(problems, fem_options), Combinatorial(fem)))
fem_run_directories = [details(r).data.fem.run.run_directory for r in fem_results]

scaled_error(a, e) = maximum(abs, a - e) / maximum(abs, diag(e))
e_Z, e_Y = map((Z, Y)) do quantity
    [scaled_error(observe(f, quantity)[:, :, k], observe(u, quantity)[:, :, k])
     for (f, u) in zip(fem_results, unified_results), k in eachindex(frequencies)]
end

println("\nFEM vs Unified, worst frequency per point")
@printf("%-46s %9s %12s %9s %12s\n", "point", "e_Z %", "at f [Hz]", "e_Y %", "at f [Hz]")
for (i, point) in enumerate(points)
    kZ = argmax(e_Z[i, :])
    kY = argmax(e_Y[i, :])
    @printf("%-46s %9.3f %12.4g %9.3f %12.4g\n", point,
        100e_Z[i, kZ], frequencies[kZ], 100e_Y[i, kY], frequencies[kY])
end

# FEM (solid) against Unified (dashed) for one point; change plot_point and rerun.
plot_point = 5
point_plot = Makie.plot(fem_results[plot_point], unified_results[plot_point];
    ydata = (Z, Y),
    length_unit = :base,
    series_labels = ["FEM", "Unified"],
    series_attributes = [(linestyle = :solid,), (linestyle = :dash,)],
    xscale = :log10,
    yscale = :log10,
    backend = :gl,
    display_plot = true,
    figure_title = "Two bare wires — " * points[plot_point]
);

nothing
