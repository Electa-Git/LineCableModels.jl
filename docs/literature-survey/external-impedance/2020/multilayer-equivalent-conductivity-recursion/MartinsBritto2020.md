# Martins-Britto–Lopes–Rondineau multilayer equivalent-conductivity recursion

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Thin overhead conductors; conductor internals/insulation are outside the recursion. |
| Calculated quantities | Frequency-dependent real equivalent conductivity that maps an ``N``-layer earth into a homogeneous Carson earth-return model |
| Earth structure | ``N`` flat horizontal layers; bottom layer semi-infinite. |
| Model and approximation | Each adjacent pair is replaced by the source's two-layer real equivalent conductivity, recursively from the bottom. The resulting ``\sigma_{eq}(f)`` is inserted into a homogeneous Carson kernel; it is not an exact multilayer reflection factor. |
| Main source | A. G. Martins-Britto, F. V. Lopes, and S. R. M. J. Rondineau (2020) |
| Citation key(s) | `:MartinsBritto2020` |
| Evidence status | Publication page image verified |

**Description.** Bottom-up pairwise homogenization of any horizontally layered earth into a frequency-dependent real conductivity for use in homogeneous-earth Carson formulas or EMTP line-constants tools.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Carson quasi-TEM overhead parent; no independent longitudinal term in the equivalent-conductivity recursion. | Stated — §III-A. |
| Air propagation constant ``γ_air`` | Homogeneous Carson parent uses free-space constants. | Stated — (2)–(7). |
| Earth propagation constant ``γ_earth`` | Represented indirectly through skin-penetration factors ``\sqrt{\pi f\mu_n\sigma_n}``. | Stated — (14)–(22). |
| Earth permittivity and displacement current | Recursion omits ``\epsilon_n`` although the benchmark multilayer/Carson parents retain it; validation assumes vacuum permittivity. | Stated — §§III–V. |
| Range of validity | Tested 1 Hz–2 MHz on 2–6 layer models; reported errors grow sharply for high conductivity contrast and depend on layer depth relative to skin depth. | Stated — §§V–VI. |
| Earth permeability ``μ_earth`` | Each recursion step retains layer ``\mu_n``; validation assumes ``\mu_n=\mu_0``. | Stated — (20)–(22), §V. |
| Arrangement | Applied to overhead self/mutual Carson geometry; equivalent material is geometry-independent within the stated approximation. | Stated — Fig. 2 and §IV. |
| Earth structure | ``N`` flat horizontal layers; bottom layer semi-infinite. | Stated — §III-A, Fig. 2. |
| Conductor and insulation geometry | Thin overhead conductors; conductor internals/insulation are outside the recursion. | Stated — Fig. 2. |
| Constitutive and field assumptions | Linear isotropic layers; real-valued conductivity homogenization driven by conductive skin penetration. | Stated — §§II, IV. |
| Conventions | Frequency ``f`` in Hz; layers numbered top-to-bottom; pair replacements proceed bottom-to-top. | Stated — §IV. |

**Expression.** Initialize the bottom pair by

```math
\sigma_{N-1,N}=\sigma_{N-1}\left[
\frac{(\sqrt{\sigma_{N-1}}+\sqrt{\sigma_N})-(\sqrt{\sigma_{N-1}}-\sqrt{\sigma_N})e^{-2h_{N-1}\sqrt{\pi f\mu_{N-1}\sigma_{N-1}}}}
{(\sqrt{\sigma_{N-1}}+\sqrt{\sigma_N})+(\sqrt{\sigma_{N-1}}-\sqrt{\sigma_N})e^{-2h_{N-1}\sqrt{\pi f\mu_{N-1}\sigma_{N-1}}}}
\right]^2.
\qquad\text{(20)}
```

Then recurse upward,

```math
\sigma_{m-1,m}=\sigma_{m-1}\left[
\frac{(\sqrt{\sigma_{m-1}}+\sqrt{\sigma_{m-1,m}})-(\sqrt{\sigma_{m-1}}-\sqrt{\sigma_{m-1,m}})e^{-2h_{m-1}\sqrt{\pi f\mu_{m-1}\sigma_{m-1}}}}
{(\sqrt{\sigma_{m-1}}+\sqrt{\sigma_{m-1,m}})+(\sqrt{\sigma_{m-1}}-\sqrt{\sigma_{m-1,m}})e^{-2h_{m-1}\sqrt{\pi f\mu_{m-1}\sigma_{m-1}}}}
\right]^2,
\qquad\text{(21)}
```

and obtain ``\sigma_{eq}`` by the same top-layer replacement printed in (22).

**Approximation.** Each adjacent pair is replaced by the source's two-layer real equivalent conductivity, recursively from the bottom. The resulting ``\sigma_{eq}(f)`` is inserted into a homogeneous Carson kernel; it is not an exact multilayer reflection factor.

**Limitations.** Accuracy is source-tested rather than guaranteed. Large conductivity contrast, deep layers and displacement-dominated regimes can invalidate the approximation. It does not yield a complex equivalent permittivity or permeability.

**Reference.** [MartinsBritto2020](@cite), equations (1)–(23), especially (20)–(22), PDF pp. 4–5.

**Transcription source.** Publication page images. Every square root, contrast sign, exponent, layer index, outer square and recursion direction was visually verified.

## Source transcription

Equation (22) is the final top-layer instance of (21), with ``\sigma_1,h_1,\mu_1`` and the accumulated lower equivalent. The source then uses it in the homogeneous Carson equations (2)–(4).

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\sigma_n`` | unchanged | physical conductivity of layer ``n`` | ``\mathrm{S/m}`` |
| ``\sigma_{m-1,m}`` | unchanged | accumulated equivalent conductivity | ``\mathrm{S/m}`` |
| ``\sigma_{eq}`` | unchanged | final homogeneous equivalent | ``\mathrm{S/m}`` |
| ``h_n`` | unchanged | finite layer thickness | m |

No notation was renamed.

## Evidence and approximation sources

The two-layer expression (14) is explicitly attributed to Tsiamitros et al.; the paper's contribution is its recursive extension (15)–(22) and validation across measured multilayer profiles.

## Limitations and discrepancies

- Equation (21)'s accumulated symbol is unusual (``\sigma_{m-1,m}`` appears on both sides at different recursion stages). The printed indexing is retained rather than rewritten as an algorithmic temporary variable.
