# Matrix-based formulation of line and cable coaxial parameters

## 1. Scope and physical organization

The parameters of a coaxial cable system are most naturally organized by separating the response of each radially layered cable from the response of the regions shared by several cables. The latter comprise the surrounding air or earth and, when present, a common metallic pipe. This separation is the central organizing principle of Ametani's general formulation [Ametani1980](@cite) and of Chapter 2 of *Cable System Transients: Theory, Modeling and Simulation* [Ametani2015b](@cite). It permits the local conductor models and the surrounding-field models to be assembled without expanding their individual analytical expressions.

The formulation developed here retains that block structure but does not restrict a coaxial unit to a core, sheath, and armor. A unit may contain any number $N$ of longitudinal metallic conductors: a solid central conductor followed by $N-1$ concentric annular conductors, separated by dielectric regions. The familiar core/sheath and core/sheath/armor constructions are the particular cases $N=2$ and $N=3$. A common enclosing pipe is an additional conductor, not a condition fixing $N$. The arbitrary-$N$ notation and compact assembly operators below are a reformulation of the block patterns in [Ametani1980, Ametani2015b](@cite), rather than equations reproduced verbatim from those sources.

Two meanings of “internal” must be distinguished. **Material-internal impedance** describes electromagnetic diffusion, loss, and magnetic energy inside a metal. **Cable-internal impedance**, as used in Ametani's decomposition, includes both these metallic contributions and the series magnetic-field contributions of the dielectric annuli within the chosen cable boundary. Consequently, a cable-internal block is not simply a resistance or metal-only impedance matrix. Conversely, an *outer-surface impedance* remains a material-internal quantity: “outer” identifies the metal boundary at which its response is represented, not the surrounding air or earth [Ametani1980](@cite) (Sec. 2.1); [Ametani2015b](@cite) (Sec. 2.1.1).

The discussion uses the conventional linear, longitudinally uniform cable-constants approximation. Each local cable is treated as cylindrically symmetric; arbitrary proximity-induced current redistribution is not resolved by its local block. The pipe-cavity model can describe the relative positions of several coaxial units, but this does not make their local conductor models fully proximity-aware. These qualifications become important when replacing discrete wires or noncircular conductors by equivalent coaxial layers [Ametani2015b](@cite) (Sec. 2.5.1); [Gustavsen2001](@cite) (Sec. III); [Gustavsen2005](@cite) (Sec. II).

## 2. Conductor coordinates and general block assembly

At angular frequency $\omega$, with time dependence $e^{j\omega t}$, the multiconductor transmission-line equations are

$$
\frac{\mathrm d\mathbf V}{\mathrm dx}=-\mathbf Z\mathbf I,
\qquad
\frac{\mathrm d\mathbf I}{\mathrm dx}=-\mathbf Y\mathbf V.
$$

Here $\mathbf Z$ and $\mathbf Y$ are per-unit-length matrices, in $\Omega/\mathrm m$ and $\mathrm S/\mathrm m$, respectively. Voltages are initially referred to a common external reference, and all conductor currents are positive in the same longitudinal direction. These are conductor coordinates, not modal coordinates or preselected core-to-sheath loop coordinates [Ametani1980](@cite) (Eqs. (1)–(4)); [Ametani2015b](@cite) (Eqs. (2.1)–(2.4)).

Consider $m$ coaxial units. Unit $j$ contains $N_j$ metallic conductors, ordered from the center outward, and the number of conductors excluding a possible common pipe is

$$
n=\sum_{j=1}^{m}N_j.
$$

Let $\mathbf B_j\in\mathbb C^{N_j\times N_j}$ denote the complete local cable-internal impedance block, including its metal and internal dielectric-region series contributions. Its entries need not be exposed to assemble the system. Define

$$
\mathbf B=\operatorname{blockdiag}(\mathbf B_1,\ldots,\mathbf B_m),
\qquad
\mathbf S=\operatorname{blockdiag}
\left(\mathbf 1_{N_1},\ldots,\mathbf 1_{N_m}\right)
\in\mathbb R^{n\times m},
$$

where $\mathbf 1_{N_j}$ is a column of ones. The rectangular matrix $\mathbf S$ relates individual conductor currents to the total longitudinal current of each coaxial unit:

$$
\mathbf J=\mathbf S^{\mathsf T}\mathbf I,
\qquad
J_j=\sum_{k=1}^{N_j}I_{j,k}.
$$

Within the cylindrical cable approximation, the field outside a unit is driven by this total current. A shared-region impedance matrix $\mathbf H\in\mathbb C^{m\times m}$ therefore contributes $\mathbf S\mathbf H\mathbf S^{\mathsf T}$ in conductor coordinates. Its block between units $j$ and $k$ is

$$
\left(\mathbf S\mathbf H\mathbf S^{\mathsf T}\right)_{jk}
=H_{jk}\mathbf 1_{N_j}\mathbf 1_{N_k}^{\mathsf T}.
$$

Thus, the repeated entries in Ametani's external and pipe-related blocks express a physical common contribution to all conductors inside each coaxial boundary. They are not a special property of three-conductor cables. The same construction supports different values of $N_j$ within one installation [Ametani1980](@cite) (Secs. 2.1–2.2); [Ametani2015b](@cite) (Eqs. (2.7), (2.14), and (2.29)–(2.33)).

The construction of a local block is immaterial to the shared-region assembly. Once a cable assembly is represented in conductor coordinates and associated with an enclosing circular boundary, the surrounding-region response depends on the geometry exposed to that region and on the total current carried by the conductors enclosed by the boundary. The detailed electromagnetic mechanisms retained inside the local block therefore remain separated from the analytical kernel used for the surrounding region.

Without a common pipe, the impedance assembly is simply

$$
\boxed{\mathbf Z=\mathbf B+\mathbf S\mathbf Z_{\mathrm{env}}\mathbf S^{\mathsf T}},
$$

where $\mathbf Z_{\mathrm{env}}$ contains the self and mutual impedances of the regions outside the individual cable boundaries, including the appropriate air/earth-return response. The local part is block diagonal; coupling between distinct cables enters through the shared-region matrix. This distinction does **not** imply that $\mathbf B_j$ is diagonal: its own conductors are coupled through their nested fields and metallic layers.

A round line conductor is the limiting local case $N_j=1$. Several such conductors can have an arbitrary line arrangement through their shared external matrix. What must be coaxial here is the local layered representation, not the entire arrangement of all line or cable axes.

Carson's overhead-line integral and Pollaczek's induction formulas supply the classical homogeneous-earth contributions to $\mathbf Z_{\mathrm{env}}$ [Carson1926, Pollaczek1926](@cite). Subsequent formulations include complex-depth approximations for overhead lines [Deri1981](@cite) and spectral expressions for conductors above, within, or between horizontal earth layers [Tsiamitros2008](@cite). Their selection depends on conductor placement, earth stratification, and the treatment of permeability, displacement current, and longitudinal propagation. A source that gives only the correction for imperfect earth must be combined with its corresponding geometric field term before it represents the complete external entry used here.

Mixed overhead–buried arrangements require the corresponding mutual entries in the same external matrix. Uribe provides an evaluation of Pollaczek's mixed integral within its original induction assumptions [Uribe2008](@cite), whereas Martins-Britto retains earth displacement current and a longitudinal propagation parameter in a generalized homogeneous-earth mixed formulation [MartinsBritto2024](@cite). These choices change the external coefficients, not the conductor-current summation represented by $\mathbf S$.

## 3. Physical content of a coaxial impedance block

### 3.1 Metallic surface and transfer contributions

For a solid core surrounded by annular metallic layers, four types of material-internal quantity supply the local impedance model. Their analytical evaluation is separate from the assembly problem.

| Primitive | Symbol used here | Physical role and location |
| --- | --- | --- |
| Outer-surface impedance of a solid disk | $z_{1}^{\mathrm o}$ | Represents the solid core's internal electromagnetic response at its outer cylindrical boundary; belongs to the local cable block. |
| Inner-surface impedance of an annulus | $z_{k}^{\mathrm i}$ | Represents the annular metal's response associated with its inner boundary; contributes to the region facing the conductor or layers inside it. |
| Outer-surface impedance of an annulus | $z_{k}^{\mathrm o}$ | Represents the same annular metal's response associated with its outer boundary; contributes to the region outside that layer. |
| Transfer impedance across an annulus | $z_{k}^{\mathrm t}$ | Couples the responses at the inner and outer boundaries through the finite metal wall; belongs to the same annular conductor model. |

The surface and transfer terms in $\mathbf B_j$ follow the cylindrical solutions developed by Schelkunoff, which express solid-conductor and annular impedances through modified Bessel functions [Schelkunoff1934](@cite). Ametani incorporates these quantities into the cable matrix [Ametani1980](@cite), with later treatments addressing bonded conducting and semiconducting layers [Ametani2004](@cite). The radial solutions assume concentric homogeneous layers and linear materials; the usual conductor-diffusion approximation neglects displacement current in the metal. Skin effect is retained, while angular current redistribution requires a proximity-effect treatment. Surface-admittance methods address skin and proximity effects in solid and hollow conductors [Patel2014](@cite). Their coupled terms supply a conductor-coordinate block for the represented assembly, after separation of any shared-region contribution already accounted for elsewhere.

For annular conductors, $k=2,\ldots,N$. The “disk” describes the core cross-section; the longitudinal object is a solid cylindrical conductor. The quantities above are the per-unit-length impedances used in cable formulations, not unnormalized surface impedances defined solely as an electric-to-magnetic-field ratio. In particular, their units are $\Omega/\mathrm m$ [Ametani1980](@cite) (Sec. 2.1); [Ametani2015b](@cite) (Eqs. (2.9)–(2.10)).

The two surface responses and the transfer response describe one coupled annular-metal problem. They should not be interpreted as three independent resistors added along a physical current path. Transfer impedance is also distinct from the mutual impedance between two separate cables: it describes coupling *through the wall of one conductor*. Its placement and sign depend on the common longitudinal-current convention and on the conversion from surface quantities to conductor coordinates. In the present discussion that conversion remains encapsulated in $\mathbf B_j$.

Adding another concentric metallic layer therefore introduces another annular set $(z^{\mathrm i},z^{\mathrm o},z^{\mathrm t})$, enlarges the conductor-coordinate block, and introduces the associated intervening dielectric region. It does not require a new system-level assembly principle. A sheath and a continuous armor are distinguished by geometry and material properties, not by different types of matrix object.

### 3.2 Nonmetal regions inside and outside the cable boundary

The magnetic field in the insulation between adjacent metallic conductors also produces a longitudinal series voltage gradient. Although this contribution is external to the metals, it is internal to the cable geometry. The same reasoning applies to a protective jacket between the outermost metal and the chosen outer cable radius. Ametani includes these annular field terms in the local cable-internal block together with the metallic surface and transfer terms [Ametani1980](@cite) (Sec. 2.1); [Ametani2015b](@cite) (Eqs. (2.8)–(2.13)).

The distinction can be summarized as follows:

| Physical region or mechanism | Allocation in the assembly |
| --- | --- |
| Solid core and annular cable metals | Material-internal contributions within $\mathbf B_j$. |
| Dielectric annuli and jacket up to the chosen cable boundary | Nonmetal series-field contributions within $\mathbf B_j$. |
| Air or earth outside individual cables, with no common pipe | $\mathbf S\mathbf Z_{\mathrm{env}}\mathbf S^{\mathsf T}$. |
| Cavity between the individual cable boundaries and a common pipe's inner surface | Pipe-interior shared block, described below. |
| Pipe metal | Inner-surface response in the pipe-interior block; outer-surface and through-wall transfer contributions in the pipe-wall assembly. |
| Pipe outer jacket and the region beyond it | Pipe-jacket contribution and external return contribution, respectively. |

The cable boundary is an accounting boundary and must be used consistently. If a jacket's annular field contribution is already included in $\mathbf B_j$, the external or cavity calculation must begin at the jacket's outer radius, not again at the underlying metal surface. Otherwise the same field region is counted twice.

The series contribution of a dielectric region must also be distinguished from its shunt contribution: magnetic energy in that region contributes to $\mathbf Z$, whereas electric-field storage and dielectric leakage contribute to $\mathbf Y$. Calling a region “insulation” does not place all its electromagnetic effects exclusively in the admittance matrix.

Ametani's logarithmic annular expressions supply the insulation series-field contributions in the conventional coaxial approximation [Ametani1980](@cite). The survey also includes Wait's electrically thin jacket expression, which retains dependence on the longitudinal propagation constant [Wait1978](@cite). Such a term can supply the same physical region only when its electrical-thickness and propagation assumptions are consistent with those of the surrounding-field calculation.

## 4. Particular coaxial cases: $N=2$ and $N=3$

For a core/sheath coaxial unit, the conductor ordering is $(c,s)$. Its unexpanded local block $\mathbf B_{cs}^{(2)}$ contains the solid-core outer-surface response; the sheath inner-surface, outer-surface, and transfer responses; and the field contributions of the core–sheath insulation and any outer jacket. For one such unit without a pipe,

$$
\mathbf Z_{cs}
=\mathbf B_{cs}^{(2)}
+z_{\mathrm{env}}\mathbf u_2\mathbf u_2^{\mathsf T},
\qquad
\mathbf u_2=\begin{bmatrix}1\\1\end{bmatrix},
\qquad
\mathbf B_{cs}^{(2)}\in\mathbb C^{2\times2}.
$$

For a core/sheath/armor unit, the ordering is $(c,s,a)$. The block $\mathbf B_{csa}^{(3)}$ contains the same core and sheath primitives, an additional annular set for the armor, the sheath–armor dielectric-region contribution, and any jacket outside the armor:

$$
\mathbf Z_{csa}
=\mathbf B_{csa}^{(3)}
+z_{\mathrm{env}}\mathbf u_3\mathbf u_3^{\mathsf T},
\qquad
\mathbf u_3=\begin{bmatrix}1\\1\\1\end{bmatrix},
\qquad
\mathbf B_{csa}^{(3)}\in\mathbb C^{3\times3}.
$$

These are the block-level counterparts of the two- and three-conductor cases in [Ametani2015b](@cite) (Eqs. (2.8)–(2.12)). They illustrate the change in dimension and physical ingredients without expanding the local matrix entries. For arbitrary $N$, the same statement is $\mathbf Z^{(N)}=\mathbf B^{(N)}+z_{\mathrm{env}}\mathbf 1_N\mathbf 1_N^{\mathsf T}$. The value of $z_{\mathrm{env}}$ must, of course, correspond to the actual outer geometry of each construction.

No return-current condition has been imposed in these primitive matrices. A grounded sheath, a bonded armor, or a selected core–sheath excitation is a subsequent circuit or coordinate condition; it is not a reason to omit that conductor's electromagnetic contribution during parameter assembly.

## 5. Common enclosing pipe

### 5.1 Four assembly contributions

For a finite metallic pipe retained as an explicit conductor, use the ordering

$$
\widetilde{\mathbf I}
=\begin{bmatrix}\mathbf I\\I_p\end{bmatrix},
\qquad
\widetilde{\mathbf V}
=\begin{bmatrix}\mathbf V\\V_p\end{bmatrix}.
$$

The full matrices now have order $n+1$. Ametani's series decomposition is

$$
\mathbf Z=\mathbf Z_i+\mathbf Z_p+\mathbf Z_c+\mathbf Z_0,
$$

where the four terms describe the individual coaxial units, the pipe-interior response, the pipe-wall/exterior-surface connection, and the external return region, respectively [Ametani1980](@cite) (Eq. (3) and Sec. 2.2); [Ametani2015b](@cite) (Eqs. (2.28)–(2.40)). The subscripts are region-based bookkeeping labels, not a pure metal/nonmetal split.

Let $\mathbf H_p\in\mathbb C^{m\times m}$ be the impedance matrix between the inner coaxial units with respect to the pipe inner surface. Its diagonal elements describe self responses and its off-diagonal elements describe mutual responses in the common enclosure. Then

$$
\mathbf Z_i=
\begin{bmatrix}
\mathbf B&\mathbf 0\\
\mathbf 0^{\mathsf T}&0
\end{bmatrix},
\qquad
\mathbf Z_p=
\begin{bmatrix}
\mathbf S\mathbf H_p\mathbf S^{\mathsf T}&\mathbf 0\\
\mathbf 0^{\mathsf T}&0
\end{bmatrix}.
$$

The kernel $\mathbf H_p$ includes both the cavity magnetic-field contribution and the pipe's inner-surface metallic response. It is therefore not an ideal-pipe geometric inductance alone. For a centered coaxial configuration, these are the annular-gap and inner-surface contributions familiar from radial layering. For multiple or eccentric inner units, the kernel also represents their positions relative to the pipe and to one another; its analytical expansion is not required for the block assembly [Ametani2015b](@cite) (Eqs. (2.33)–(2.36)).

This partition is a bookkeeping convention rather than a unique physical decomposition. The pipe inner-surface term may equivalently be retained in a coupled two-surface pipe-wall operator and omitted from $\mathbf H_p$, provided that the same contribution is not counted twice. The assembled $\mathbf Z$ is unchanged by such a redistribution of terms between compatible local operators.

The pipe-type literature supplies the self and mutual coefficients needed to populate $\mathbf H_p$ and hence $\mathbf Z_p$. Kane's filament–shield construction includes selected proximity terms [Kane1995](@cite); Yang treats finite pipe-wall thickness [Yang2001](@cite), and Fortin evaluates the finite-wall eddy-current response [Fortin2005](@cite). De Silva subsequently provides corrected finite-pipe coefficients [DeSilva2019](@cite). These formulations differ in their treatment of wall thickness, eccentric conductor placement, and angular field harmonics. In the partition used here, their cavity and inner-wall contributions enter $\mathbf H_p$. A published complete core–pipe loop impedance must first be separated from the individual-conductor contributions already assigned to $\mathbf B$.

To describe the remaining pipe contributions, define

$$
a_p=z_p^{\mathrm o}+z_{p,\mathrm{jacket}},
\qquad
t_p=z_p^{\mathrm t},
\qquad
\mathbf u=\mathbf 1_n,
\qquad
\mathbf v=\begin{bmatrix}\mathbf u\\1\end{bmatrix}.
$$

Here $z_p^{\mathrm o}$ is the pipe-metal outer-surface impedance, $z_p^{\mathrm t}$ is its through-wall transfer impedance, and $z_{p,\mathrm{jacket}}$ is the series field contribution of its outer insulating jacket, if present. The pipe inner-surface response is already assigned to $\mathbf H_p$ and must not be added again to $a_p$.

With the common current convention used above, the remaining macro blocks are

$$
\mathbf Z_c=
\begin{bmatrix}
(a_p-2t_p)\mathbf u\mathbf u^{\mathsf T}&(a_p-t_p)\mathbf u\\
(a_p-t_p)\mathbf u^{\mathsf T}&a_p
\end{bmatrix},
\qquad
\mathbf Z_0=z_e\mathbf v\mathbf v^{\mathsf T}.
$$

The scalar $z_e$ is the external self/return impedance of the enclosing pipe, evaluated outside its chosen jacket boundary. Its common contribution is driven by the total current of the enclosed assembly, $\mathbf u^{\mathsf T}\mathbf I+I_p$. The combinations $a_p-2t_p$, $a_p-t_p$, and $a_p$ are the pipe-level coefficients of [Ametani2015b](@cite) (Eq. (2.37)), written for arbitrary $n$; they are not an expansion of the local cable blocks or of the underlying surface-impedance formulas.

For a buried pipe, the external coefficient in $\mathbf Z_0$ can be obtained from a homogeneous-earth buried-conductor formulation, such as Pollaczek's integral treatment [Pollaczek1926](@cite), or from an approximation with a stated range of validity, such as Saad's closed-form expressions [Saad1996](@cite). Stratified soil instead calls for a compatible layered-earth expression [Tsiamitros2008](@cite). The source's self-term prescription is applied at the chosen outer boundary, subject to its conductor-radius and insulation assumptions. This calculation supplies the field outside the pipe and jacket; their already allocated contributions are not added again.

Writing $g_p=a_p+z_e$ gives the complete assembly in one expression:

$$
\boxed{
\mathbf Z=
\begin{bmatrix}
\mathbf B+\mathbf S\mathbf H_p\mathbf S^{\mathsf T}
+(g_p-2t_p)\mathbf u\mathbf u^{\mathsf T}
&(g_p-t_p)\mathbf u\\
(g_p-t_p)\mathbf u^{\mathsf T}&g_p
\end{bmatrix}.
}
$$

This form locates all relevant mechanisms: each cable's surface and transfer impedances remain in its own $\mathbf B_j$; the common cavity and pipe inner surface enter through $\mathbf H_p$; the pipe outer surface and transfer response enter through $g_p$ and $t_p$; and the exterior return contribution is included once, through $z_e$. A block between distinct inner units contains their cavity mutual response and the common pipe/exterior contribution, but no local $\mathbf B_j$ contribution.

### 5.2 Conductor-to-pipe mutual impedance versus return-loop impedance

The off-diagonal block between the inner conductors and the pipe is

$$
\mathbf Z_{\mathrm{inner},p}=(g_p-t_p)\mathbf u.
$$

Each entry is a mutual impedance in the externally referenced conductor matrix. It is not, by itself, the impedance measured with an inner conductor carrying the outgoing current and the pipe carrying the return current. For conductor $\alpha$, that loop impedance is

$$
z_{\alpha\text{–}p}^{\mathrm{loop}}
=Z_{\alpha\alpha}+Z_{pp}-Z_{\alpha p}-Z_{p\alpha}.
$$

The corresponding transformation for all pipe-return loops is obtained by defining

$$
\mathbf K=
\begin{bmatrix}
\mathbf I_n\\-\mathbf u^{\mathsf T}
\end{bmatrix},
\qquad
\widetilde{\mathbf I}=\mathbf K\mathbf i,
\qquad
\mathbf w=\mathbf K^{\mathsf T}\widetilde{\mathbf V}
=\mathbf V-\mathbf uV_p,
$$

where $\mathbf I_n$ is the identity matrix. This excitation imposes $I_p=-\mathbf u^{\mathsf T}\mathbf i$: the pipe carries the sum of the return currents. Applying the transformation to the assembled matrix gives

$$
\boxed{
\mathbf Z_{\mathrm{loop},p}
=\mathbf K^{\mathsf T}\mathbf Z\mathbf K
=\mathbf B+\mathbf S\mathbf H_p\mathbf S^{\mathsf T}.
}
$$

This identity follows directly from the preceding block assembly. The common pipe-exterior and transfer bookkeeping terms cancel for this balanced pipe-return excitation; the pipe's inner-surface response remains in $\mathbf H_p$. The cancellation does not make the pipe a perfect conductor or establish zero total current for every possible operating condition. It distinguishes an externally referenced primitive matrix from a particular family of differential return loops.

For a single enclosed coaxial unit, let $h_p$ be its scalar pipe-interior kernel. The two examples become

$$
\mathbf Z_{cs,p}^{\mathrm{loop}}
=\mathbf B_{cs}^{(2)}+h_p\mathbf u_2\mathbf u_2^{\mathsf T},
\qquad
\mathbf Z_{csa,p}^{\mathrm{loop}}
=\mathbf B_{csa}^{(3)}+h_p\mathbf u_3\mathbf u_3^{\mathsf T}.
$$

Their full externally referenced matrices have orders $3$ and $4$, with orderings $(c,s,p)$ and $(c,s,a,p)$, respectively. In general, one $N$-conductor coaxial unit plus a pipe has order $N+1$, while $m$ enclosed units have order $1+\sum_jN_j$. A single centered unit and a concentric pipe can equivalently be viewed as a radial stack with one additional annular metal; the separate pipe notation is particularly useful when several units share the enclosure.

### 5.3 Applicability of the pipe kernel

Keeping the pipe conductor explicit does not, on its own, remove approximations in the selected cavity kernel. Chapter 2 distinguishes finite-pipe assembly from the thick-wall approximation used in its pipe-interior formulas, and discusses the relationship between pipe-wall thickness and electromagnetic penetration depth [Ametani2015b](@cite) (Sec. 2.5.1). A thin-wall or low-frequency application must use an appropriate kernel and compatible surface/transfer quantities; the block notation alone cannot guarantee accuracy outside their validity range.

Similarly, neglecting proximity-induced changes in the individual cable impedances and in the pipe outer-surface impedance remains an assumption of this construction. The infinite-thickness, pipe-referenced model in [Ametani2015b](@cite) (Eq. (2.27)) is a particular shielding/reference limit, not a justification for deleting the pipe row and column merely because a physical pipe is grounded at its terminals.

The low-frequency analyses of Høidalen distinguish a finite pipe wall from a wall first assumed infinitely thick [Hoidalen2013](@cite); inner-surface and transfer terms drawn from incompatible limits do not describe one physical pipe. More recent work isolates a method-of-moments proximity correction by subtracting its zero-harmonic response before adding the increment to a compatible finite-pipe base [Hoidalen2025](@cite). This separation preserves the existing skin-effect contribution. Where proximity couples individual conductor responses, the correction must be assembled in the coordinates and coupled entries it actually describes.

## 6. Shunt admittance and potential-coefficient assembly

The shunt parameters follow the same geometric partition, but their physical primitives are different. In the electrostatic, lossless-dielectric formulation, the potential-coefficient matrix $\mathbf P$ relates conductor voltages to per-unit-length charges $\mathbf q$:

$$
\mathbf V=\mathbf P\mathbf q,
\qquad
\mathbf C=\mathbf P^{-1},
\qquad
\mathbf Y=j\omega\mathbf P^{-1}.
$$

Thus $\mathbf P$ has units $\mathrm m/\mathrm F$, and $\mathbf C$ has units $\mathrm F/\mathrm m$. Ametani assembles $\mathbf P$ before taking its inverse [Ametani1980](@cite) (Eq. (4)); [Ametani2015b](@cite) (Secs. 2.1.2 and 2.2.2).

The same membership matrix used for the series assembly aggregates the conductor charges associated with each enclosed boundary:

$$
\mathbf Q=\mathbf S^{\mathsf T}\mathbf q,
\qquad
Q_j=\sum_{k=1}^{N_j}q_{j,k}.
$$

Accordingly, a shared dielectric region described by a potential-coefficient matrix $\mathbf H_P$ contributes in conductor coordinates as

$$
\boxed{
\Delta\mathbf P
=\mathbf S\mathbf H_P\mathbf S^{\mathsf T}
}
$$

in direct analogy with the series contribution

$$
\boxed{
\Delta\mathbf Z
=\mathbf S\mathbf H_Z\mathbf S^{\mathsf T}.
}
$$

The common matrix form reflects the same geometric membership, while $\mathbf H_Z$ and $\mathbf H_P$ represent different electromagnetic primitives and must not be identified with one another.

Let $\mathbf D_j$ be the unexpanded local potential-coefficient block of coaxial unit $j$, and set $\mathbf D=\operatorname{blockdiag}(\mathbf D_1,\ldots,\mathbf D_m)$. Without a common pipe,

$$
\mathbf P=\mathbf D+\mathbf S\mathbf P_{\mathrm{env}}\mathbf S^{\mathsf T}.
$$

With an explicit pipe, let $\mathbf H_{P,p}$ describe the cavity potential coefficients relative to the pipe inner surface, and let $p_{\mathrm{ext}}$ collect the pipe-jacket and any exterior-space potential coefficients. The assembly is

$$
\mathbf P=
\begin{bmatrix}
\mathbf D+\mathbf S\mathbf H_{P,p}\mathbf S^{\mathsf T}&\mathbf 0\\
\mathbf 0^{\mathsf T}&0
\end{bmatrix}
+p_{\mathrm{ext}}\mathbf v\mathbf v^{\mathsf T}.
$$

This is the arbitrary-dimension counterpart of $\mathbf P_i+\mathbf P_p+\mathbf P_c+\mathbf P_0$ in [Ametani2015b](@cite) (Eqs. (2.41)–(2.51)). In particular, the scalar in the book's $\mathbf P_c$ is associated with the pipe's outer insulation [Ametani2015b](@cite) (Eq. (2.50b)); it is not a capacitance through the conducting pipe wall analogous to the series transfer impedance. Surface and transfer impedances belong to the longitudinal metal response, whereas the potential coefficients describe electric fields in the dielectric and external regions.

The inversion is a system-level operation: the inverse of the sum is generally not the sum of the inverses. Consequently, the shunt-admittance matrix cannot be obtained by independently inverting and adding these geometric blocks. The reference must also be properly specified; a singular matrix produced by an ideal reference constraint must first be expressed in independent voltage coordinates.

For directly buried cables, the conventional electrostatic treatment takes the surrounding soil as the reference equipotential, so that the separate exterior potential-coefficient term vanishes. The dielectric jacket contributions remain. This simplification of $\mathbf P$ does not eliminate the earth-return contribution to $\mathbf Z$ [Ametani2015b](@cite) (Eqs. (2.17)–(2.18) and (2.42)–(2.43)). Dielectric losses, when required, are represented through $\mathbf Y=\mathbf G+j\omega\mathbf C$ or an appropriate complex-permittivity formulation; the lossless potential-coefficient model should not be read as a complete description of lossy insulation [Gustavsen2005](@cite) (Secs. III.B and IX.A.4).

For the local dielectric regions, the surveyed admittance formulas extend from series-connected insulation and semiconducting-screen paths [Weeks1984](@cite) to complex-permittivity descriptions of conducting screens [Ametani2004](@cite) and radial networks with multiple semiconducting layers [Ghosh2022](@cite). They supply the local shunt response associated with $\mathbf D_j$, after conversion from radial branch or loop quantities to the conductor coordinates used here. A radial branch admittance is not itself an entry of the full conductor admittance matrix. Material conductivity, dielectric loss, and any assumed frequency dependence must be consistent across the retained layers.

A finite-conductivity earth also admits an external potential-coefficient treatment in place of the equipotential-soil approximation. Wise's overhead-line expression is an early example [Wise1948](@cite); subsequent formulations address buried cables in homogeneous earth [Papadopoulos2010b](@cite), cables within the upper soil layer of a two-layer earth [Papadopoulos2011](@cite), and complete-field and quasi-TEM buried-cable descriptions [Xue2018b](@cite). De Conti provides a closed-form approximation to the homogeneous-earth quasi-TEM potential coefficients [DeConti2023a](@cite). In the compatible complex-potential formulation, these expressions supply $\mathbf P_{\mathrm{env}}$ or the exterior part of $p_{\mathrm{ext}}$ before the full matrix inversion. Their soil and propagation assumptions determine which expression can be used.

For mixed overhead–buried conductors, Martins-Britto provides the mutual potential coefficient for homogeneous air and earth half-spaces [MartinsBritto2024](@cite). It enters the external potential matrix alongside the appropriate self and other mutual terms, with insulation contributions retained. It is not an isolated admittance to be added after inversion.

## 7. Explicit and equivalent representations of discrete conductors

Actual cables may contain stranded cores, wire screens, helically laid armor, and other structures that are not continuous coaxial metals. Such structures need not always be homogenized before the block formulation is applied. When the individual conductors are approximately circular and an analytical kernel is available for their positions within a common enclosing conductive boundary, they may be retained explicitly as separate conductor coordinates. Their self and mutual interactions in the shared region are then represented by the corresponding shared-region matrix, while each conductor retains its own material-internal response.

The limiting case of a single round wire is already included by $N_j=1$. Consequently, a collection of round wires inside a common circular pipe can be represented directly by treating the wires as distinct inner units and using a pipe-interior kernel for their relative positions and radii. The surrounding pipe remains an explicit conductor and its inner-surface, outer-surface, and transfer responses are handled by the same pipe assembly described in Section 5. This representation preserves the individual wire currents and voltages rather than replacing them a priori by one equivalent conductor.

Homogenization nevertheless remains useful when the individual wires do not need to be resolved. Gustavsen describes the conversion of manufacturer data into equivalent cable-constants representations [Gustavsen2001](@cite) (Secs. III and V), and the IEEE PES Task Force gives corresponding guidance for practical cable models [Gustavsen2005](@cite) (Secs. II–IV and IX). A stranded core may be replaced by a homogeneous solid conductor with effective resistivity selected to reproduce metallic fill or specified DC resistance. A wire screen may be replaced by an equivalent tubular conductor preserving its total metallic cross-sectional area and a suitable radial position. Effective insulation data may similarly be chosen to preserve capacitance when semiconducting layers are not represented explicitly [Gustavsen2001](@cite) (Sec. V); [Gustavsen2005](@cite) (Secs. III–IV and IX.A).

The explicit and homogenized descriptions therefore correspond to different levels of retained conductor detail, not to different system-level assembly principles. Their interchangeability is limited by the physical quantity and frequency range of interest. Matching metallic area, DC resistance, or capacitance does not establish equivalence of skin, proximity, helical, magnetic, or transfer-impedance effects. In particular, a wired magnetic armor requires suitable effective material assumptions rather than an equal-area tube alone. When such omitted behavior is significant, the discrete conductors or a more detailed field model must be retained [Gustavsen2001](@cite) (Sec. III.C); [Gustavsen2005](@cite) (Secs. II–III).

The survey distinguishes such equivalent-conductor formulas from field solutions for the actual cross-section. Ametani's area-and-perimeter equivalent-annulus approximation covers noncircular shapes, including sector conductors, but does not resolve their corner fields or proximity effects [Ametani1992](@cite). Reluctance-network formulations instead calculate a discretized coupled response for arbitrary cross-sections [Bormann2013](@cite), while a later two-dimensional finite-element treatment approximates helical effects in armored three-core cables through effective material properties [Gustavsen2023](@cite). These are alternatives at different levels of geometric approximation; a coupled cross-sectional impedance cannot generally replace one scalar coaxial primitive.

## 8. Synthesis

The essential matrix structure is independent of the names or number of metallic layers. Each local cable assembly contributes a conductor-coordinate block containing the electromagnetic mechanisms retained inside its chosen boundary. Shared regions couple the exposed units through matrices constructed from their boundary geometry and aggregate conductor currents or charges. A common pipe adds its cavity and inner-surface response, its outer-surface and through-wall transfer contributions, and one exterior return contribution.

The cases $N=2$ and $N=3$ differ in local dimension and in the number of annular conductor models, not in the assembly principle. The same distinction also accommodates collections of explicitly retained round conductors inside a common circular enclosure. Maintaining the separation between local material response, shared-region coupling, and voltage/current reference makes the matrix assembly applicable to arbitrary combinations of layered coaxial and pipe-enclosed circular assemblies, provided that an appropriate analytical kernel is available for each shared region. None of these extensions requires the local block entries or their underlying surface-impedance formulas to be expanded at system level.

## References

```@bibliography
Pages = [@__FILE__]
```
