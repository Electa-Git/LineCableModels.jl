
struct ScaledBessel <: InternalImpedanceFormulation end
get_description(::ScaledBessel) = "Scaled Bessel (Schelkunoff)"


@inline function (f::ScaledBessel)(
	form::Symbol,
	r_in::T,
	r_ex::T,
	rho_c::T,
	mur_c::T,
	jω::Complex{T},
) where {T <: REALSCALAR}
	Base.@nospecialize form
	return form === :inner  ? f(Val(:inner), r_in, r_ex, rho_c, mur_c, jω)  :
		   form === :outer  ? f(Val(:outer), r_in, r_ex, rho_c, mur_c, jω)  :
		   form === :mutual ? f(Val(:mutual), r_in, r_ex, rho_c, mur_c, jω) :
		   throw(ArgumentError("Unknown ScaledBessel form: $form"))
end

@inline function (f::ScaledBessel)(
	::Val{:inner},
	r_in::T,
	r_ex::T,
	rho_c::T,
	mur_c::T,
	jω::Complex{T},
) where {T <: REALSCALAR}
	# Constants
	mu_c = T(μ₀) * mur_c
	sigma_c = _to_σ(rho_c)

	# Calculate the reciprocal of the skin depth
	m = sqrt(jω * mu_c * sigma_c)
	w_ex = m * r_ex

	if isapprox(r_in, 0.0, atol = eps(T))
		return zero(Complex{T}) # not physical, but consistent with :outer - algorithmic shortcut for solids/tubular blending
	else

		w_in = m * r_in

		sc_in = exp(abs(real(w_in)) - w_ex)
		sc_ex = exp(abs(real(w_ex)) - w_in)
		sc = sc_in / sc_ex

		# Bessel function terms with uncertainty handling
		N =
			(besselkx(0, w_in)) *
			(besselix(1, w_ex)) +
			sc *
			(besselix(0, w_in)) *
			(besselkx(1, w_ex))

		D =
			(besselkx(1, w_in)) *
			(besselix(1, w_ex)) -
			sc *
			(besselix(1, w_in)) *
			(besselkx(1, w_ex))

		return Complex{T}((jω * mu_c / 2π) * (1 / w_in) * (N / D))
	end

end

@inline function (f::ScaledBessel)(
	::Val{:outer},
	r_in::T,
	r_ex::T,
	rho_c::T,
	mur_c::T,
	jω::Complex{T},
) where {T <: REALSCALAR}
	# Constants
	mu_c = T(μ₀) * mur_c
	sigma_c = _to_σ(rho_c)

	# Calculate the reciprocal of the skin depth
	m = sqrt(jω * mu_c * sigma_c)
	w_ex = m * r_ex

	if isapprox(r_in, 0.0, atol = eps(T)) # solid conductor
		@debug "Using closed form for solid conductor"
		N = besselix(0, w_ex)
		D = besselix(1, w_ex)

	else
		w_in = m * r_in

		sc_in = exp(abs(real(w_in)) - w_ex)
		sc_ex = exp(abs(real(w_ex)) - w_in)
		sc = sc_in / sc_ex

		# Bessel function terms with uncertainty handling
		N =
			(besselix(0, w_ex)) *
			(besselkx(1, w_in)) +
			sc *
			(besselkx(0, w_ex)) *
			(besselix(1, w_in))

		D =
			(besselix(1, w_ex)) *
			(besselkx(1, w_in)) -
			sc *
			(besselkx(1, w_ex)) *
			(besselix(1, w_in))
	end

	return Complex{T}((jω * mu_c / 2π) * (1 / w_ex) * (N / D))
end

@inline function (f::ScaledBessel)(
	::Val{:mutual},
	r_in::T,
	r_ex::T,
	rho_c::T,
	mur_c::T,
	jω::Complex{T},
) where {T <: REALSCALAR}
	# Constants
	mu_c = T(μ₀) * mur_c
	sigma_c = _to_σ(rho_c)

	# Calculate the reciprocal of the skin depth
	m = sqrt(jω * mu_c * sigma_c)
	w_ex = m * r_ex

	if isapprox(r_in, 0.0, atol = eps(T))

		return zero(Complex{T}) # not physical, but consistent with :outer - algorithmic shortcut for solids/tubular blending 
	# return f(Val(:outer), r_in, r_ex, rho_c, mur_c, freq)

	else

		w_in = m * r_in

		sc_in = exp(abs(real(w_in)) - w_ex)
		sc_ex = exp(abs(real(w_ex)) - w_in)
		sc = sc_in / sc_ex

		# Bessel function terms with uncertainty handling
		N = 1.0 / sc_ex

		D =
			(besselix(1, w_ex)) *
			(besselkx(1, w_in)) -
			sc *
			(besselix(1, w_in)) *
			(besselkx(1, w_ex))

		return Complex{T}((1 / (2π * r_in * r_ex * sigma_c)) * (N / D))

	end

end
