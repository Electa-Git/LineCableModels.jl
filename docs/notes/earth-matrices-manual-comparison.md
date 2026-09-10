# Manual comparison of the earth matrices

This calculation evaluates only earth matrices for two equal bare wires:
separation 1 m, depth 1 m, radius 0.0425 m, earth resistivity 0.1 Ω·m,
relative permittivity and permeability 1 in both half-spaces, ideal air,
Γ = 0, and voltage referenced to infinite earth depth.
Frequencies are 0.1, 1, 10, 100, 1000, 10000, 100000, and 1000000 Hz.

All 80 supplied complex reference values are reproduced within 4.49e-8 relative
complex difference. They correspond to **isolated primary-current normalization**,
using Cf = 1/[κg r K1(κg r)]. The full manuscript's physical-current closure
uses the entire L matrix and is included as a separately labelled fourth curve.
The prototype is not changed to make these two normalizations coincide.

Every tabulated quantity is the complete symmetric matrix

$$M=\begin{bmatrix}M_{11}&M_{12}\\M_{12}&M_{11}\end{bmatrix}.$$

Ze is in Ω/m and Pe in m/F. Ye = jω Pe⁻¹ is in S/m. The inverse is a full
2×2 matrix solve, never the reciprocal of each entry. No internal impedance,
insulation admittance, conductor material, cable assembly, or FEM enters this run.

## Reproduction and sources

```sh
julia --startup-file=no --compiled-modules=existing --project=. test/gauntlet/earth_matrices_manual.jl
julia --startup-file=no --compiled-modules=existing --project=. test/gauntlet/plot_earth_matrices_manual.jl
```

The first script evaluates the manuscript's kernels directly and independently
calls the currently registered underground Xue equations in
`src/engine/earthimpedance/formulas/default.jl` and
`src/engine/earthadmittance/formulas/default.jl`. These indexed earth functors
are called directly with EarthPair; the cable solver is not called.
The registered formula identifiers remain unchanged in this manual test.

The user-supplied numbers are retained independently in
`test/gauntlet/earth_matrices_expected.tsv`. Results and figures are written to
`.linecablemodels/qa/earth-matrices-manual/`: `matrices.json` includes all
intermediate K, H, L, A, D values, method comparisons and reference checks;
`matrices.csv` contains every matrix entry; `expected-comparison.csv` compares
all supplied numbers. Figures are `earth-Ze-Ye`, `earth-Pe`, and
`earth-relative-differences` (PNG/PDF, also SVG for the first two).

The supplied LaTeX is the framework source, specifically its Lambda-self,
Lambda-mutual, Z-pair-same, H-gg, old-C-self-basis, directional-C,
D-implementation and global-PY equations. Xue's source is
[General Formulation and Accurate Evaluation of Earth-Return Parameters for
Overhead/Underground Cables (2018)](https://publications.polymtl.ca/3190/1/2018_HaoyanXue.pdf).
The kernel reductions below are algebraic checks of the equations in this
checkout, and are restricted to equal permeability and Γ = 0.

## Definitions used in this manual calculation

Write s = jω, σ̂0 = sε0, σ̂g = 1/ρ + sεg, κm² = sμσ̂m,
am = sqrt(λ² + κm²), A = I0(κg r), D = κg r K1(κg r), and Cf = 1/D.
The square roots use the outgoing/decaying branch. Define, for horizontal
separation y and equal depths h,

$$Q(y)=\int_0^\infty\frac{e^{-2ha_g}\cos(\lambda y)}{a_0+a_g}\,d\lambda,$$

$$U(y)=\int_0^\infty\frac{\hat\sigma_g a_0\,e^{-2ha_g}\cos(\lambda y)}
{a_g(\hat\sigma_g a_0+\hat\sigma_0a_g)}\,d\lambda.$$

For either B = Q or B = U, define the two dimensionless matrices by

$$T^a_{11}(B)=K_0(\kappa_g r)+A[-K_0(2\kappa_g h)+2B(0)],$$

$$T^a_{12}(B)=A[K_0(\kappa_g d)-K_0(\kappa_g\sqrt{d^2+4h^2})+2B(d)],$$

$$T^p_{11}(B)=K_0(\kappa_g r)-K_0(\kappa_g\sqrt{r^2+4h^2})+2B(r),$$

$$T^p_{12}(B)=K_0(\kappa_g d)-K_0(\kappa_g\sqrt{d^2+4h^2})+2B(d).$$

Here d = 1 m is the centre separation, h = 1 m, superscript a means a
circumferential average, and p means Xue's point-field sampling. In particular,
point-field self sampling puts r in both the horizontal cosine and image
distance. Merely dropping I0 from the averaged self kernel does not reproduce
that convention.

With cZ = sμ/(2π) and cP = s/(2πσ̂g):

| Curve | Ze | Pe | Ye |
|---|---|---|---|
| Field average · Cf | cZ Tᵃ(Q)/D | cP Tᵃ(U)/D | s Pe⁻¹ |
| Point field · Cf | cZ Tᵖ(Q)/D | cP Tᵖ(U)/D | s Pe⁻¹ |
| Xue | cZ Tᵖ(Q) | cP Tᵖ(U) | s Pe⁻¹ |
| Full manuscript | K L⁻¹ | H L⁻¹ | s L H⁻¹ |

In the last row K = cZ Tᵃ(Q), H = cP Tᵃ(U), and

$$L=A^{-1}I-FK,\qquad F=\frac{2\pi\hat\sigma_g r I_1(\kappa_g r)}{\kappa_g A}.$$

Cf is an earth-side source-current factor. It is not a conductor internal
impedance. In the complete physical-current elimination, a common source-column
Cf cancels between the field map and current map; it does not replace L⁻¹.

## Where the differences originate

Xue's two impedance integrands combine into Q because
(λ² + κg²)/ag² = 1. For potential, put β = κ0²/κg² = σ̂0/σ̂g. Then

$$\frac{\lambda^2}{a_g^2(a_0+\beta a_g)}+
\frac{\kappa_g^2}{a_g^2(a_0+a_g)}
=\frac{a_0}{a_g(a_0+\beta a_g)}.$$

The right side is exactly the U integrand before the exponential/cosine weight.
Thus Xue's S12 + κg²S13 already contains the a0/ag factor in this specialization.
The discrepancy here is not a missing a0/ag multiplier in Xue.

For point-field sampling the following matrix identities are exact:

$$Z_e^p=C_f Z_e^X,\qquad P_e^p=C_f P_e^X,\qquad Y_e^p=D Y_e^X.$$

For the averaged mutual entries, Ze12ᵃ/Ze12ˣ = Pe12ᵃ/Pe12ˣ = A/D.
There is no analogous independent entrywise rule for Ye12; the entire Pe must
be inverted. At 1 MHz A = 0.99968200 + j0.03565268,
D = 0.94697550 − j0.11433829, and Cf = 1.04082018 + j0.12566914.
That complex factor explains the visible high-frequency self differences.
As frequency decreases, A and D approach 1 and these corrections shrink.
The remaining small self-potential difference includes the horizontal-radius
sampling versus actual circumferential sampling near the interface.

At 1 MHz the computed mutual admittances are:

| Normalization | Ye12 (S/m) | Relative complex difference from Xue |
|---|---|---:|
| Field average · Cf | −0.0209285155405 − j0.0164119523789 | 9.44041% |
| Point field · Cf | −0.0214933369280 − j0.0156506349104 | 12.6035% |
| Xue | −0.0204038969801 − j0.0189905457748 | — |
| Full manuscript | −0.0204038329058 − j0.0189910530095 | 0.00183420% |

Consequently the supplied Ze and self-Ye tables, and nearly identical mutual Ye
to Xue, do not all select the same current normalization. The supplied tables
select Cf, whereas the very close high-frequency mutual Ye is obtained with the
complete L closure. The latter also changes mutual Ze: at 1 MHz it is
0.000153570976500 + j0.00106298092351 Ω/m, versus
0.000231620017139 + j0.000998261393900 Ω/m for field-average/Cf.
These are distinguishable model choices, not numerical integration errors.

At 1 MHz the Cf field-average differences from Xue are 13.2133% for self Ze,
16.6673% for mutual Ze, 12.6035% for self Ye, and 9.44041% for mutual Ye,
using the magnitude of the complex difference divided by the Xue magnitude.
Although mutual Ye is close on a whole-range linear plot, its 1 MHz relative
complex discrepancy is not negligible under Cf normalization.

## Numerical checks

All 80 supplied complex values pass relative tolerance 5e-8, with separate
real/imaginary component checks at 2e-6. The largest difference (4.48e-8)
is Xue mutual Ye at 100 kHz; most entries agree much more closely.
The reduced kernels agree with the independently evaluated registered
S11/S12/S13 implementation to 5.39e-11 relative matrix error or better.
The full-L matrices reproduce the earlier direct manuscript prototype to
8.97e-16 relative error; no previous FEM results were used.

All eight frequencies pass both trapz and CIM. Quad uses rtol=1e-10;
trapz uses rtol=1e-6; CIM uses rtol=1e-6 and explicit dimensionless integral
atol=1e-6. No method silently falls back to another method.
Maximum relative differences against quad across all four conventions and
Ze, Pe, Ye are:

| Method | Complete matrix | Mutual entry |
|---|---:|---:|
| trapz | 4.75e-8 | 6.60e-8 |
| CIM | 1.78e-7 | 4.29e-6 |

The largest normalized Ye Pe − sI residual is 8.19e-16.
Current-map right solves and the Cf/A identities are checked separately.
This is a manual numerical and algebraic comparison, not a claim that either
normalization is validated against a conductor-surface boundary-value solution.

## Complete matrices at the eight frequencies

Each row below supplies both distinct entries; M22 = M11 and M21 = M12.
Values are computed results, not copied reference values. Pe includes its s
factor, and should not be confused with an inverse-admittance matrix in Ω·m.

### Ze: self (Ω/m)

| Hz | Field average · Cf | Point field · Cf | Xue | Full manuscript · L⁻¹ |
|---:|---|---|---|---|
| 0.1 | 9.90248763839e-08 + j1.21216893381e-06 | 9.90248776264e-08 + j1.21216893364e-06 | 9.90249604731e-08 + j1.21216892003e-06 | 9.90248349247e-08 + j1.21216894335e-06 |
| 1 | 9.97139401069e-07 + j1.06677704681e-05 | 9.97139499476e-07 + j1.06677704504e-05 | 9.97145907638e-07 + j1.06677692488e-05 | 9.97136537528e-07 + j1.06677712672e-05 |
| 10 | 1.01736063594e-05 + j9.19861976001e-05 | 1.0173613606e-05 + j9.19861958251e-05 | 1.01740897703e-05 + j9.19860910146e-05 | 1.01734261074e-05 + j9.19862629602e-05 |
| 100 | 0.000106940444326 + j0.000768373130341 | 0.000106940907233 + j0.000768372952923 | 0.000106974250849 + j0.000768363925644 | 0.000106930999962 + j0.000768378291106 |
| 1000 | 0.00115614676161 + j0.00605096302851 | 0.00115616693564 + j0.00605094696412 | 0.00115827841705 + j0.00605019230925 | 0.00115586075167 + j0.00605130645918 |
| 10000 | 0.011116357366 + j0.0429419845247 | 0.0111163726206 + j0.0429412697053 | 0.0112293424876 + j0.0428863908291 | 0.0111189271426 + j0.0429499376184 |
| 100000 | 0.0915576276723 + j0.284970041703 | 0.0915562904363 + j0.284972784807 | 0.0966204143304 + j0.281599002407 | 0.0915604682382 + j0.28494947132 |
| 1e+06 | 0.753096795349 + j1.57307087634 | 0.753096832116 + j1.57307086926 | 0.893026473375 + j1.40355176677 | 0.753096655932 + j1.57307095354 |

### Ze: mutual (Ω/m)

| Hz | Field average · Cf | Point field · Cf | Xue | Full manuscript · L⁻¹ |
|---:|---|---|---|---|
| 0.1 | 9.90237336015e-08 + j8.15291478414e-07 | 9.90237365083e-08 + j8.15291478061e-07 | 9.90237920486e-08 + j8.15291466681e-07 | 9.90236751889e-08 + j8.15291490173e-07 |
| 1 | 9.97039189555e-07 + j6.69900464024e-06 | 9.97039428401e-07 + j6.6990046047e-06 | 9.97043431743e-07 + j6.69900362537e-06 | 9.97034956493e-07 + j6.69900566167e-06 |
| 10 | 1.01650164133e-05 + j5.22994453426e-05 | 1.01650350601e-05 + j5.22994417183e-05 | 1.01653033334e-05 + j5.22993591794e-05 | 1.01647317957e-05 + j5.22995329316e-05 |
| 100 | 0.000106232271476 + j0.000371605899324 | 0.000106233596397 + j0.000371605520561 | 0.000106249416951 + j0.000371598746592 | 0.000106215651598 + j0.000371613278422 |
| 1000 | 0.00110273172454 + j0.00209562402273 | 0.00110280644075 + j0.00209558470399 | 0.00110349836351 + j0.0020950707315 | 0.00110205467226 + j0.00209618183357 |
| 10000 | 0.00804020091534 + j0.00484349995878 | 0.00804192704768 + j0.00484063284841 | 0.00805086370794 + j0.00481560329801 | 0.0080340309126 + j0.00486665134368 |
| 100000 | 0.00854186089967 − j0.008800966623 | 0.00851040077042 − j0.00883133763271 | 0.00829041366028 − j0.00894869591373 | 0.00872821625835 − j0.00870879258607 |
| 1e+06 | 0.000231620017139 + j0.0009982613939 | 0.000266967476589 + j0.000989057803893 | 0.000365898832667 + j0.000906088901382 | 0.0001535709765 + j0.00106298092351 |

### Pe: self (m/F)

| Hz | Field average · Cf | Point field · Cf | Xue | Full manuscript · L⁻¹ |
|---:|---|---|---|---|
| 0.1 | 0.0157074556445 + j0.144461677871 | 0.0157074556269 + j0.144459420568 | 0.0157074654782 + j0.144459418678 | 0.0157074500502 + j0.144461679284 |
| 1 | 0.157038226791 + j1.21436399898 | 0.157038225023 + j1.21434142595 | 0.157038952043 + j1.21434126278 | 0.157037855865 + j1.21436411508 |
| 10 | 1.56758749933 + j9.84162491267 | 1.56758732239 + j9.84139918278 | 1.56763799782 + j9.84138545943 | 1.56756560015 + j9.84163405063 |
| 100 | 15.4797766873 + j75.4467139821 | 15.4797590733 + j75.4444570312 | 15.4830051052 + j75.4433506991 | 15.4787487753 + j75.4473776106 |
| 1000 | 143.394671961 + j529.52697505 | 143.393000323 + j529.504605964 | 143.575407066 + j529.423383595 | 143.369935357 + j529.563981971 |
| 10000 | 1030.19149055 + j3386.62708068 | 1030.0861571 + j3386.4659851 | 1038.90943815 + j3381.71225366 | 1030.46728874 + j3387.26372981 |
| 100000 | 7243.84416954 + j22670.9589997 | 7243.77313073 + j22671.403588 | 7646.87938152 + j22403.7821234 | 7244.07333855 + j22669.337786 |
| 1e+06 | 59930.25106 + j125180.696797 | 59930.2549985 + j125180.695421 | 71065.429131 + j111690.728586 | 59930.2399664 + j125180.70294 |

### Pe: mutual (m/F)

| Hz | Field average · Cf | Point field · Cf | Xue | Full manuscript · L⁻¹ |
|---:|---|---|---|---|
| 0.1 | 0.0157072067653 + j0.111763478907 | 0.0157072071638 + j0.111763478851 | 0.0157072147655 + j0.111763477144 | 0.0157071997741 + j0.111763480503 |
| 1 | 0.157017857071 + j0.88738478288 | 0.15701788871 + j0.887384777282 | 0.157018417616 + j0.887384632441 | 0.157017373313 + j0.887384917299 |
| 10 | 1.56600228389 + j6.57210973933 | 1.5660046271 + j6.57210918099 | 1.56603817751 + j6.5720972969 | 1.56597178724 + j6.5721207109 |
| 100 | 15.3663676779 + j42.7789967262 | 15.3665202017 + j42.7789419386 | 15.3683238833 + j42.7780235494 | 15.3647492782 + j42.7798439975 |
| 1000 | 136.504899461 + j205.430290138 | 136.512223732 + j205.425423004 | 136.578460647 + j205.364823004 | 136.448254048 + j205.485113404 |
| 10000 | 739.252627569 + j334.298754266 | 739.371747792 + j334.035149686 | 739.885256126 + j331.796048414 | 738.845035917 + j336.164604573 |
| 100000 | 666.469457574 − j683.958542898 | 664.024539061 − j686.328234078 | 646.913563404 − j695.500330552 | 681.276499315 − j676.605813265 |
| 1e+06 | 18.4274632844 + j79.4238821368 | 21.2397859599 + j78.6916504415 | 29.1110252827 + j72.0905440386 | 12.2165453114 + j84.5741344949 |

### Ye: self (S/m)

| Hz | Field average · Cf | Point field · Cf | Xue | Full manuscript · L⁻¹ |
|---:|---|---|---|---|
| 0.1 | 10.8157943415 + j0.148174450827 | 10.8164679699 + j0.148176911321 | 10.8164680203 + j0.148177656361 | 10.8157943959 + j0.148174800244 |
| 1 | 11.0700332429 + j0.219075255399 | 11.0707114232 + j0.219078967774 | 11.0707119104 + j0.219085688207 | 11.0700337878 + j0.219077951356 |
| 10 | 11.4554340342 + j0.357216166879 | 11.45611942 + j0.357219007181 | 11.4561239644 + j0.357279219226 | 11.4554394667 + j0.3572351252 |
| 100 | 12.104617469 + j0.682506405029 | 12.105307176 + j0.682487061136 | 12.1053447986 + j0.683025634549 | 12.104669782 + j0.682616482021 |
| 1000 | 13.4220890785 + j1.62767340455 | 13.4227074822 + j1.62754038688 | 13.4228720182 + j1.63245933913 | 13.4225072307 + j1.628030342 |
| 10000 | 17.0856646384 + j4.24671778425 | 17.0859461628 + j4.24654238478 | 17.0835786396 + j4.29636370965 | 17.0865340758 + j4.24599503976 |
| 100000 | 25.1139585055 + j8.0613785804 | 25.113792807 + j8.06126495956 | 25.0849974323 + j8.59959815582 | 25.1140530075 + j8.06179703956 |
| 1e+06 | 40.8337744755 + j19.5491744055 | 40.8337745847 + j19.5491748169 | 40.0438886142 + j25.4787209246 | 40.8337746262 + j19.5491734068 |

### Ye: mutual (S/m)

| Hz | Field average · Cf | Point field · Cf | Xue | Full manuscript · L⁻¹ |
|---:|---|---|---|---|
| 0.1 | -8.39989797728 + j0.148028192399 | -8.40055095266 + j0.148030877289 | -8.4005510099 + j0.148030300131 | -8.39989803052 + j0.148027720975 |
| 1 | -8.14582199641 + j0.217878171373 | -8.14647014658 + j0.217883682034 | -8.14647073485 + j0.217878757981 | -8.14582252709 + j0.217874236327 |
| 10 | -7.76204469598 + j0.347898625151 | -7.76268292304 + j0.347914996852 | -7.76268909285 + j0.347874527478 | -7.76204999186 + j0.347867039028 |
| 100 | -7.12877289671 + j0.615735907047 | -7.12938048308 + j0.615807741401 | -7.12944760025 + j0.615496248286 | -7.12882628282 + j0.615497583379 |
| 1000 | -5.95591566968 + j1.21572829549 | -5.95636568001 + j1.21608824803 | -5.95713588629 + j1.21401374507 | -5.95648058738 + j1.21413363116 |
| 10000 | -3.31389189285 + j2.30229269369 | -3.31322089018 + j2.30371972311 | -3.32145135911 + j2.29578500719 | -3.32022793044 + j2.29552877665 |
| 100000 | 0.187887969262 + j1.04152480807 | 0.191606114169 + j1.04081776705 | 0.172028204452 + j1.05001589629 | 0.171688523172 + j1.0502290511 |
| 1e+06 | -0.0209285155405 − j0.0164119523789 | -0.021493336928 − j0.0156506349104 | -0.0204038969801 − j0.0189905457748 | -0.0204038329058 − j0.0189910530095 |
