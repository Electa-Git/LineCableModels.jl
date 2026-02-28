using Revise
include("CableBuilder.jl")
using .CableBuilder
using BenchmarkTools

ms_cu = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393, 0.0, 0.0, 0.0, 0.0)

parts = (
	Conductor.Solid(:core, ms_cu; r = 0.02315),
	Conductor.Tubular(:sheath, ms_cu; t = 0.005),
)

# Flawless compilation
cds = CableDesignSpec(parts)

ast = first(cds)
println(ast.payload)


# 1. Define the Blueprint (The User DSL)
# Notice we mix deterministic scalars, a relative sweep, and an absolute sweep
copper_spec = MaterialSpec(
	rho = Grid(1.68e-8:0.01e-8:1.72e-8, 2.0), # 2% manufacturing tolerance
	eps_r = 1.0,                              # Auto-promoted to DeterministicGrid
	mu_r = 1.0,
	T0 = Grid(20.0:10.0:90.0, AbsoluteError(2.0)), # ± 2.0°C absolute sensor error
)

@show typeof(copper_spec)
@code_warntype rand(copper_spec)

@show @allocated rand(copper_spec)
@btime @allocated rand($copper_spec)


# ---------------------------------------------------------
# Benchmark 1: The Array Comprehension (User Snippet)
# ---------------------------------------------------------
println("--- Benchmarking 1000 Monte Carlo Realizations (Array Allocation) ---")
# We interpolate $copper_spec to avoid global scope penalty in the benchmark
@btime mc_materials = [rand($copper_spec) for _ in 1:1000]

# ---------------------------------------------------------
# Benchmark 2: The Bare-Metal Loop (Zero-Allocation Target)
# ---------------------------------------------------------
# This tests the pure speed of your ntuple/Grid/rand architecture 
# without the overhead of Julia allocating a 1000-element Vector
function pure_monte_carlo_loop(spec, N)
	# We just draw the sample and discard it to test the engine's raw speed
	for _ in 1:N
		m = rand(spec)
	end
	return nothing
end

println("\n--- Benchmarking 1000 Monte Carlo Realizations (Engine Only) ---")
@btime pure_monte_carlo_loop($copper_spec, 1000)


# ---------------------------------------------------------
# Test 1: The Deterministic Single Build
# ---------------------------------------------------------
println("--- Benchmarking Single Deterministic Build ---")
# We use $ to interpolate the variable into the macro so it doesn't benchmark global scope lookup
@btime first($cds)

# ---------------------------------------------------------
# Test 2: The Combinatorial Iterator (The Real Test)
# ---------------------------------------------------------
# Let's create 12 distinct cable designs on the fly using Grids
grid_parts = (
	Conductor.Solid(:core, ms_cu; r = Grid([0.02, 0.025, 0.03])), # 3 variations
	Conductor.Tubular(:sheath, ms_cu; t = Grid([0.004, 0.005, 0.006, 0.007])), # 4 variations
)
grid_cds = CableDesignSpec(grid_parts)

# A function barrier to test the loop exactly how your solver will use it
function exhaust_generator(spec)
	count = 0
	for design in spec
		# The design is materialized here. 
		# If the compiler is happy, this loop will allocate ZERO memory.
		count += 1
	end
	return count
end

println("\n--- Benchmarking 12-Design Combinatorial Sweep ---")
@btime exhaust_generator($grid_cds)
