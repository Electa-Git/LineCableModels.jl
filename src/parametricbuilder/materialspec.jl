# Use lib/material nominal; kw is either percent-only or (value,pct)
_pair_from_nominal(nom, x) =
	x === nothing                 ? (nom, nothing) :
	(x isa Tuple && length(x)==2) ? x :
	(nom, x)

# -------------------- material spec --------------------

"""
MaterialSpec: pass specs for fields (value spec + optional %unc)

Example:
  MaterialSpec(; rho=(2.826e-8, nothing),
				 eps_r=(1.0, nothing),
				 mu_r=(1.0, nothing),
				 T0=(20.0, nothing),
				 alpha=(4.0e-3, nothing))
"""
struct MaterialSpec
	rho::Any;
	eps_r::Any;
	mu_r::Any;
	T0::Any;
	alpha::Any;
	kappa::Any;
end
MaterialSpec(; rho, eps_r, mu_r, T0, alpha, kappa) = MaterialSpec(rho, eps_r, mu_r, T0, alpha, kappa)

# --- 1) Ad-hoc numeric: values (or (value,pct)) ---
Material(; rho, eps_r = 1.0, mu_r = 1.0, T0 = 20.0, alpha = 0.0, kappa = 1.0) =
	MaterialSpec(
		rho = _spec(rho),
		eps_r = _spec(eps_r),
		mu_r = _spec(mu_r),
		T0 = _spec(T0),
		alpha = _spec(alpha),
		kappa = _spec(kappa),
	)

# --- 2) From an existing Material: append %unc by default, or override with (value,pct) ---
function Material(
	m::Materials.Material;
	rho = nothing,
	eps_r = nothing,
	mu_r = nothing,
	T0 = nothing,
	alpha = nothing,
	kappa = nothing,
)
	MaterialSpec(
		rho   = _pair_from_nominal(m.rho, rho),
		eps_r = _pair_from_nominal(m.eps_r, eps_r),
		mu_r  = _pair_from_nominal(m.mu_r, mu_r),
		T0    = _pair_from_nominal(m.T0, T0),
		alpha = _pair_from_nominal(m.alpha, alpha),
		kappa = _pair_from_nominal(m.kappa, kappa),
	)
end

# --- 3) From a MaterialsLibrary + name ---
Material(lib::Materials.MaterialsLibrary, name::AbstractString; kwargs...) =
	Material(get(lib, name); kwargs...)
Material(lib::Materials.MaterialsLibrary, name::Symbol; kwargs...) =
	Material(lib, String(name); kwargs...)


function _make_range(ms::MaterialSpec)
	ρs = _make_range(ms.rho[1]; pct = ms.rho[2])
	εs = _make_range(ms.eps_r[1]; pct = ms.eps_r[2])
	μs = _make_range(ms.mu_r[1]; pct = ms.mu_r[2])
	Ts  = _make_range(ms.T0[1]; pct = ms.T0[2])
	αs = _make_range(ms.alpha[1]; pct = ms.alpha[2])
	κs = _make_range(ms.kappa[1]; pct = ms.kappa[2])
	[Materials.Material(ρ, ε, μ, T, α, κ) for (ρ, ε, μ, T, α, κ) in product(ρs, εs, μs, Ts, αs, κs)]
end