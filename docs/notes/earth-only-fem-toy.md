# Proposed FEM toy for earth-only operators

This is an experiment specification, not a completed FEM result. The first
experiment checks the supplied line-source framework independently of its
Sommerfeld integrals. A second experiment uses actual bare PEC boundaries to
assess the retained conductor approximation. These are different comparisons.
Neither experiment includes conductor volume loss or insulation.

## First experiment: independent fields and measurement contours

Use a 2D cross-section of air above z = 0 and earth below. Place two source
locations at (y,z) = (0,-1) and (1,-1) m. Define circular measurement contours
of radius 0.0425 m. Earth conductivity is 10 S/m, and both media have vacuum
permittivity and permeability. Use the same outgoing condition as the analytical
problem, implemented with a verified infinite-domain treatment or PML.

The circles in this diagnostic are measurement contours in the background
medium, not PEC scattering boundaries. This matches the manuscript's retained
electric-line-source fields. No conductor material is assigned. It must be
labelled as a Green-field FEM test, not a simulation of exact PEC cylinders.

Excite source 1 with unit axial line current and source 2 with zero, then
interchange them. A regularized compact source can be used provided its support
is much smaller than the measurement radius and the result is extrapolated as
that support shrinks. Alternatively subtract the known homogeneous primary
singularity and solve its interface reaction by FEM. In the latter case only
the local homogeneous singularity is supplied analytically; no reflection,
transmission, Sommerfeld kernel, or final earth coefficient is supplied.

## Exact normalized Gamma = 0 field equations

A two-field limit avoids dividing numerical fields by a tiny Gamma. Write
s = j omega, sigma_hat = sigma + s epsilon, nu = 1/mu, and

- a = -Ex/s at Gamma = 0;
- b = lim(Hx/Gamma) as Gamma tends to zero;
- u = lim(Et/Gamma).

For the suppressed factor exp(-Gamma x), Maxwell's equations give the following
2D equations, with transverse coordinates (y,z):

$$-\nabla_t\cdot(\nu\nabla_t a)+s\hat\sigma a=J_x,$$

$$\mathbf H_t=(\nu\partial_z a,-\nu\partial_y a),$$

$$-\nabla_t\cdot\{\hat\sigma^{-1}(\nabla_t b+\mathbf H_t)\}
+s\mu b=0,$$

$$\mathbf u=\hat\sigma^{-1}
(\partial_z b+H_z,-\partial_y b-H_y).$$

The second equation is required even though Hx itself vanishes at Gamma = 0:
its leading coefficient contributes to the normalized transverse electric
field and thus to the earth-reference voltage.

Use continuous a and b at the interface and their natural conservative flux
conditions. For b, that flux contains Ht; imposing a homogeneous condition on
its derivative alone would omit the interface coupling. These equations follow
by expanding the two curl equations to the required order in Gamma, rather
than substituting an electric scalar-potential model for the transverse field.

The weak forms, for test functions v and w, have volume terms

$$\int \nu\nabla a\cdot\nabla v+s\hat\sigma av
=\int J_xv,$$

$$\int \hat\sigma^{-1}\nabla b\cdot\nabla w+s\mu bw
=-\int \hat\sigma^{-1}\mathbf H_t\cdot\nabla w,$$

with the chosen exterior radiation terms added consistently. Two nodal fields
are sufficient for this Gamma = 0 specialization. GetDP supports coupled
spaces and weak formulations; see its [official documentation](https://getdp.info/doc/texinfo/getdp.html#Coupled-spaces).

## Measure the maps before eliminating source amplitudes

For each unit source q, measure on each circle p:

$$K_{pq}=-\langle E_x\rangle_p=s\langle a\rangle_p,$$

$$L_{pq}=\oint_{C_p}\mathbf H_t\cdot d\boldsymbol\ell,$$

$$H_{pq}=-s\left\langle\int_{-\infty}^{z_p(\theta)}
u_z(y_p(\theta),z)\,dz\right\rangle_\theta.$$

In the last equation, `u_z` is the normalized electric-field component defined
above (the symbol is u, not the permeability reciprocal nu). Average complete
vertical paths to circumference endpoints, using the same deep-earth
reference as the manuscript. Avoid angular samples whose path passes through
a line source. Stop numerical path integration before entering a PML and
verify convergence as the reference depth increases.

The independently integrated flux must satisfy

$$\oint_{C_p}\hat\sigma\mathbf u\cdot\mathbf n\,dl=L_{pq}.$$

The final physical-current earth matrices are then

$$Z_e=KL^{-1},\qquad P_e=HL^{-1},\qquad Y_e=sLH^{-1}.$$

Store K, H and L separately. A failed final coefficient can then be traced to
an axial-field kernel, a voltage functional, or the measured current map.

To reproduce the user's isolated-current table, replace the air with identical
earth in a calibration solve and measure the isolated source's circle current
D = Ip / I_tilde. Then form Z = K/D and P = H/D. That calibration should recover
D = kappa_g r K1(kappa_g r); it need not assume that Bessel formula. Point-field
sampling at the horizontal radius for self terms, and at the receiver centre
for mutual terms, provides the corresponding point-field comparison.

For equal wires all matrices have equal diagonal entries and equal mutual
entries, but compute both source experiments before using this symmetry as a
check. Always invert the complete 2-by-2 P matrix to obtain Y.

## Second experiment: actual bare PEC surfaces

Use the same centres and radii, remove both conductor disks from the domain,
and impose ideal conductor boundaries through an explicitly driven exterior
parameter-extraction formulation. Include no finite-conductivity surrogate and
no insulation. Prescribe total conductor currents for series excitation; retain
floating receiver conductors with zero net current, not zero local surface
current. Use the same earth-reference voltage paths and measure complete
conduction-plus-displacement flux for shunt extraction.

This requires a consistent impressed longitudinal drive/port definition. A
source-free Maxwell problem with homogeneous PEC boundaries at an arbitrary
fixed Gamma is a modal problem; one cannot simply prescribe arbitrary currents
and read nonzero total tangential E on a PEC. The induced series voltage must
be obtained from the excitation reaction or a consistent external voltage
functional. Finite-Gamma calculations must pair their Z and Y with the same
Gamma and approach zero using a controlled normalization, or derive that
normalization directly in the driven weak form.

Actual PEC boundaries allow angular surface-current redistribution and
wire-to-wire proximity effects. The manuscript retains only an electric line
source and mean-field projection. Equality with exact PEC results is therefore
not the acceptance criterion for the first diagnostic. In the physical test,
resolve these additional angular degrees of freedom and quantify their effect.

The existing `test/gauntlet/getdp/pec_boundary.pro` is an exterior quasi-TEM
specialization. It removes conductor loss, but its scalar-potential output
cannot establish equality with the framework's full vertical electric-field
voltage without an additional derivation. Reusing its geometry is reasonable;
reusing that output as the voltage oracle is not yet justified.

## Minimal run and convergence requirements

Start at 1 MHz, where the Cf and full-L mutual-Y predictions differ visibly.
At this frequency the ground skin depth is approximately 0.159 m. Then add
100 kHz and 10 kHz; expand to all eight requested frequencies after the first
three pass. At 0.1 Hz the skin depth is approximately 503 m, so the original
high-frequency outer-domain size cannot be reused without verification.

Converge the local mesh, contour/path quadrature, source regularization, and
outer-domain/radiation treatment independently. Air-side radiation and
interface propagation mean that earth skin depth alone does not determine a
safe outer boundary. Repeat with expanded domains or changed PML parameters.

For the initial 1 MHz discrimination require estimated numerical uncertainty
below 0.1% in each non-negligible complex self and mutual entry. Use absolute
error around zero crossings. Resolving the much smaller 0.00183% full-L versus
Xue mutual-Y difference requires substantially tighter convergence and is not
the initial toy's acceptance threshold.

First verify the identical-half-space calibration, the flux equality above,
and Y P = s I. Then compare independent FEM K, H and L against the prototype.
Only after these pass interpret differences from a full PEC surface solution
as conductor-model differences.
