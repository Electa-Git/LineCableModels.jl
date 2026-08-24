struct Levenberg{T <: Real} <: AbstractTransformFormulation
    tol::T
end

Levenberg(; tol::Real = 1e-8) = Levenberg(float(tol))

function description(
        ::Levenberg,
)
    "Levenberg–Marquardt (frequency-tracked eigen decomposition)"
end

"""
$(TYPEDSIGNATURES)

Transform phase-domain line parameters with frequency-tracked
Levenberg–Marquardt modal decomposition.

# Arguments

- `lp`: Phase-domain line parameters.

# Returns

- A tuple containing the frequency-indexed transformation tensor and
  modal-domain [`LineParameters`](@ref). The returned impedance and admittance
  matrices are diagonal at each frequency and retain the input storage basis.

# Notes

`f.tol` controls the accepted off-diagonal ratio. A larger residual emits a
warning and retains the calculated matrix. An accepted residual is set to an
exact diagonal matrix.
"""
function (f::Levenberg)(
        lp::LineParameters{
        Tc, U, PhaseDomain, Basis},
) where {Tc <: Complex, U <: Real, Basis}
    n, n2, nfreq = size(lp.Z.values)
    n == n2 || throw(DimensionMismatch("Z must be square"))
    size(lp.Y.values) == (n, n, nfreq) || throw(DimensionMismatch("Y must be n×n×nfreq"))

    # 1) Deterministic eigen/LM on nominal arrays
    Z_nom = nominal(lp.Z.values)
    Y_nom = nominal(lp.Y.values)
    f_nom = nominal(lp.f)
    Ti, _g_nom = levenberg_transform(n, Z_nom, Y_nom, f_nom; tol = f.tol)
    _rot_min_imag!(Ti)

    Zm = similar(lp.Z.values)
    Ym = similar(lp.Y.values)

    Tk = zeros(Tc, n, n)
    Zk = zeros(Tc, n, n)
    Yk = zeros(Tc, n, n)
    invT = zeros(Tc, n, n)

    @inbounds for k in 1:nfreq
        Tk .= @view Ti[:, :, k]
        invT .= inv(Tk)
        @views begin # Enforce reciprocity.
            copyto!(Zk, lp.Z.values[:, :, k])
            reciprocity!(Zk)
            copyto!(Yk, lp.Y.values[:, :, k])
            reciprocity!(Yk)
        end
        # Modal matrices retain uncertainty values.
        @views Zm[:, :, k] .= transpose(Tk) * Zk * Tk
        @views Ym[:, :, k] .= invT * Yk * transpose(invT)

        fname = String(nameof(typeof(f)))
        offdiagZ = offdiagonal_ratio(Zm[:, :, k])
        if offdiagZ > f.tol
            @warn "$fname: transformed Z exceeds the off-diagonal tolerance" ratio = offdiagZ tolerance = f.tol
        else
            @views Zm[:, :, k] .= Diagonal(diag(Zm[:, :, k]))
        end
        offdiagY = offdiagonal_ratio(Ym[:, :, k])
        if offdiagY > f.tol
            @warn "$fname: transformed Y exceeds the off-diagonal tolerance" ratio = offdiagY tolerance = f.tol
        else
            @views Ym[:, :, k] .= Diagonal(diag(Ym[:, :, k]))
        end
    end
    return Ti, LineParameters(ModalDomain, Zm, Ym, lp.f; basis = basis(lp))
end

#= ---------------------------------------------------------------------------
Internals
-----------------------------------------------------------------------------=#

# Propagate γ with uncertainty WITHOUT eigen():
# γ̂_k = sqrt.( diag( inv(T_k) * (Y_k*Z_k) * T_k ) )
function gamma(
        Ti::AbstractArray{Tc, 3},
        Z::AbstractArray{Tu, 3},
        Y::AbstractArray{Tu, 3}
) where {Tc <: Complex, Tu <: Complex}
    n, n2, nfreq = size(Ti)
    n == n2 || throw(DimensionMismatch("Ti must be n×n×nfreq"))
    size(Z) == size(Y) == (n, n, nfreq) || throw(DimensionMismatch("Z,Y must be n×n×nfreq"))

    # Element type follows uncertain inputs
    Tγ = promote_type(eltype(Z), eltype(Y))
    Gdiag = zeros(Tγ, n, n, nfreq)  # store as diagonal matrices for consistency

    Tk = zeros(Tc, n, n)
    invT = zeros(Tc, n, n)

    @inbounds for k in 1:nfreq
        Tk .= @view Ti[:, :, k]
        invT .= inv(Tk)

        S_k = @view(Y[:, :, k]) * @view(Z[:, :, k])         # Complex{Measurement} ok
        λdiag = diag(invT * S_k * Tk)
        γdiag = sqrt.(λdiag)
        @views Gdiag[:, :, k] .= Diagonal(γdiag)
    end
    return Gdiag
end

# Frequency-tracked Levenberg–Marquardt eigen solution
function levenberg_transform(
        n::Int,
        Z::AbstractArray{T, 3},
        Y::AbstractArray{T, 3},
        f::AbstractVector{U};
        tol::Real = one(first(f)) * 1.0e-8
) where {T <: Complex, U <: Real}
    unit = one(first(f))
    ε0 = unit * 88541878128 * (unit * 10)^(-22)
    μ0 = unit * 4 * (unit * π) * (unit * 10)^(-7)

    nfreq = size(Z, 3)
    Ti = zeros(T, n, n, nfreq)
    g = zeros(T, n, n, nfreq)  # Store diagonal values in the shared n×n×nfreq shape.

    Zk = zeros(T, n, n)
    Yk = zeros(T, n, n)

    # k = 1 → plain eigen-decomposition seed
    Zk .= @view Z[:, :, 1]
    Yk .= @view Y[:, :, 1]
    S = Yk * Zk
    E = eigen(S)  # S*v = λ*v
    Ti[:, :, 1] .= E.vectors
    g[:, :, 1] .= Diagonal(sqrt.(E.values)) # γ = sqrt(λ)

    # k ≥ 2 → LM tracking
    ord_sq = n^2
    for k in 2:nfreq
        Zk .= @view Z[:, :, k]
        Yk .= @view Y[:, :, k]

        S = Yk * Zk

        # Normalise as (S / norm_val) - I.
        ω = 2 * (one(f[k]) * π) * f[k]
        nrm = -(ω^2) * ε0 * μ0
        S̃ = (S ./ nrm) - I

        # Seed from previous step
        Tseed = @view Ti[:, :, k - 1]
        gseed = @view g[:, :, k - 1]
        λseed = (diag(gseed) .^ 2 ./ nrm) .- 1  # since S̃*T = T*Λ with Λ = λ̃ = (λ/nrm)-1

        # Build real-valued unknown vector: [Re(T); Im(T); Re(λ); Im(λ)]
        x0 = [vec(real(Tseed));
              vec(imag(Tseed));
              real(λseed);
              imag(λseed)]

        function _residual!(
                F::AbstractVector{<:R},
                x::AbstractVector{<:R}
        ) where {R <: Real}
            # Unpack
            Tr = reshape(@view(x[1:ord_sq]), n, n)
            Ti_ = reshape(@view(x[(ord_sq + 1):(2 * ord_sq)]), n, n)

            λr = @view x[(2 * ord_sq + 1):(2 * ord_sq + n)]
            λi = @view x[(2 * ord_sq + n + 1):(2 * ord_sq + 2n)]

            Λr = Diagonal(λr)
            Λi = Diagonal(λi)

            Sr = real(S̃)
            Si = imag(S̃)

            # Residual of S̃*T - T*Λ = 0, split into real/imag
            Rr = (Sr*Tr - Si*Ti_) - (Tr*Λr - Ti_*Λi)
            Ri = (Sr*Ti_ + Si*Tr) - (Tr*Λi + Ti_*Λr)

            F[1:ord_sq] .= vec(Rr)
            F[(ord_sq + 1):(2 * ord_sq)] .= vec(Ri)

            # Column normalisation constraints.
            # For each column j: ||t_r||^2 - ||t_i||^2 = 1  and  t_r ⋅ t_i = 0
            c1 = sum(abs2.(Tr), dims = 1) .- sum(abs2.(Ti_), dims = 1) .- 1
            c2 = sum(Tr .* Ti_, dims = 1)
            idx = 2*ord_sq
            @inbounds for j in 1:n
                F[idx + 2j - 1] = c1[j]
                F[idx + 2j] = c2[j]
            end
            return nothing
        end

        sol = nlsolve(
            _residual!,
            x0;
            method = :trust_region,
            autodiff = :forward,
            xtol = convert(U, tol),
            ftol = convert(U, tol)
        )

        if !converged(sol)
            @warn "LM solver did not converge at k=$k, using seed eigen-decomposition fallback"
            E = eigen(S)
            Ti[:, :, k] .= E.vectors
            g[:, :, k] .= Diagonal(sqrt.(E.values))
            continue
        end

        x = sol.zero
        Tr = reshape(@view(x[1:ord_sq]), n, n)
        Ti_ = reshape(@view(x[(ord_sq + 1):(2 * ord_sq)]), n, n)
        T̂ = Tr .+ im .* Ti_

        λr = @view x[(2 * ord_sq + 1):(2 * ord_sq + n)]
        λi = @view x[(2 * ord_sq + n + 1):(2 * ord_sq + 2n)]
        λ̃ = λr .+ im .* λi   # Normalised eigenvalues.

        # Undo normalisation: λ = (λ̃ + 1) * nrm, then γ = sqrt(λ).
        λ = (λ̃ .+ one(eltype(λ̃))) .* nrm
        γ = sqrt.(λ)

        Ti[:, :, k] .= T̂
        g[:, :, k] .= Diagonal(γ)
    end

    return Ti, g
end

# Rotate columns in place to minimise their imaginary parts at each frequency.
function _rot_min_imag!(Ti::AbstractArray{T, 3}) where {T <: Complex}
    n, n2, nfreq = size(Ti)
    n == n2 || throw(DimensionMismatch("Ti must be n×n×nfreq"))
    tmp = zeros(T, n, n)
    @inbounds for k in 1:nfreq
        tmp .= @view Ti[:, :, k]
        rot!(tmp)                      # column-wise rotation in-place
        Ti[:, :, k] .= tmp
    end
    return Ti
end

# Full modal + characteristic + phase back-projection
# Returns:
#   Zm, Ym      :: n×n×nfreq   (modal-domain series/shunt matrices)
#   Zc_mod,Yc_mod :: n×n×nfreq (diagonal: per-mode characteristic)
#   Zch, Ych    :: n×n×nfreq   (phase-domain characteristic back-projected)
function modal_quantities(
        Ti::AbstractArray{Tc, 3},
        Z::AbstractArray{Tu, 3},
        Y::AbstractArray{Tu, 3}
) where {Tc <: Complex, Tu <: Complex}
    n, n2, nfreq = size(Ti)
    n == n2 || throw(DimensionMismatch("Ti must be n×n×nfreq"))
    size(Z) == size(Y) == (n, n, nfreq) || throw(DimensionMismatch("Z,Y must be n×n×nfreq"))

    Tz = promote_type(eltype(Z), eltype(Y))  # keep uncertainties
    Zm = zeros(Tz, n, n, nfreq)
    Ym = zeros(Tz, n, n, nfreq)
    Zc_mod = zeros(Tz, n, n, nfreq)
    Yc_mod = zeros(Tz, n, n, nfreq)
    Zch = zeros(Tz, n, n, nfreq)
    Ych = zeros(Tz, n, n, nfreq)

    Tk = zeros(Tc, n, n)
    Zk = zeros(Tz, n, n)
    Yk = zeros(Tz, n, n)
    invT = zeros(Tc, n, n)

    @inbounds for k in 1:nfreq
        Tk .= @view Ti[:, :, k]
        invT .= inv(Tk)
        Zk .= @view Z[:, :, k]
        Yk .= @view Y[:, :, k]

        # Modal matrices (carry uncertainties)
        @views Zm[:, :, k] .= transpose(Tk) * Zk * Tk
        @views Ym[:, :, k] .= invT * Yk * transpose(invT)

        # Characteristic per-mode (diagonal) in modal domain
        zc = sqrt.(diag(@view Zm[:, :, k])) ./ sqrt.(diag(@view Ym[:, :, k]))
        @views Zc_mod[:, :, k] .= Diagonal(zc)
        @views Yc_mod[:, :, k] .= Diagonal(inv.(zc))

        # Phase-domain characteristic back-projection
        @views Zch[:, :, k] .= transpose(invT) * Zc_mod[:, :, k] * invT
        @views Ych[:, :, k] .= Tk * Yc_mod[:, :, k] * transpose(Tk)
    end
    return Zm, Ym, Zc_mod, Yc_mod, Zch, Ych
end

# Rotate columns to minimise their imaginary parts.
function rot!(S::AbstractMatrix{T}) where {T <: Complex}
    n, m = size(S)
    n == m || throw(DimensionMismatch("Input must be square"))
    @inbounds for j in 1:n
        col = @view S[:, j]

        # optimal angle
        num = -2 * sum(real.(col) .* imag.(col))             # real
        den = sum(real.(col) .^ 2 .- imag.(col) .^ 2)           # real
        R = typeof(real(zero(T)))
        ang = R(0.5) * atan(num, den)               # real

        s1 = cis(ang)
        s2 = cis(ang + R(pi / 2))

        A = col .* s1
        B = col .* s2

        # all-real quadratic metrics
        Ar = real.(A)
        Ai = imag.(A)
        Br = real.(B)
        Bi = imag.(B)

        aaa1 = sum(Ai .^ 2)
        bbb1 = sum(Ar .* Ai)
        ccc1 = sum(Ar .^ 2)
        err1 = aaa1 * cos(ang)^2 + bbb1 * sin(2*ang) + ccc1 * sin(ang)^2   # real

        aaa2 = sum(Bi .^ 2)
        bbb2 = sum(Br .* Bi)
        ccc2 = sum(Br .^ 2)
        err2 = aaa2 * cos(ang)^2 + bbb2 * sin(2*ang) + ccc2 * sin(ang)^2   # real

        col .*= (err1 < err2 ? s1 : s2)
    end
    return S
end

# Compute the imaginary part in place for the metric term.
@inline function imag!(x::AbstractVector{<:Complex})
    @inbounds for i in eachindex(x)
        x[i] = imag(x[i])
    end
    return x
end
