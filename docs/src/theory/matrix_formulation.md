# Matrix formulation of line and cable parameters

LineCableModels calculates line and cable parameters by assembling the per-unit-length series-impedance matrix ``\mathbf Z`` and shunt-admittance matrix ``\mathbf Y``. These matrices contain the conductor self and mutual coefficients in the multiconductor transmission-line equations. They determine the propagation of voltages and currents along the line.

The matrix entries combine contributions from conducting materials, insulation, and the external medium. Each contribution is evaluated with a selected formulation. Any formulation expressible in this framework can be incorporated with compatible conductor definitions, reference conventions, and physical assumptions. The formulation sections linked below document the available equations and their domains of validity.

## Transmission-line equations

At angular frequency ``\omega``, with time dependence ``e^{j\omega t}``, the multiconductor transmission-line equations are

```math
\frac{\mathrm d\mathbf V}{\mathrm dx}=-\mathbf Z\mathbf I,
\qquad
\frac{\mathrm d\mathbf I}{\mathrm dx}=-\mathbf Y\mathbf V.
```

The coordinate ``x`` measures distance along the line. The matrices ``\mathbf Z`` and ``\mathbf Y`` have units ``\Omega/\mathrm m`` and ``\mathrm S/\mathrm m``, respectively. Their order equals the number of retained electrical conductors [Ametani1980, Ametani2015b](@cite).

In the conductor representation, each component of ``\mathbf V`` and ``\mathbf I`` corresponds to one retained conductor. All voltages are referred to a common external reference, and all currents are positive in the same longitudinal direction. Matrix indices enumerate conductors independently of their physical cross-sectional coordinates.

Sheaths, armor, and enclosing pipes remain electrical conductors when retained in this representation. Bonding conditions, return-current constraints, and modal transformations are applied separately from parameter assembly.

### Potential coefficients and shunt admittance

Shunt formulations often provide Maxwell potential coefficients rather than admittance entries. Under electrostatic, lossless-dielectric assumptions, the potential-coefficient matrix ``\mathbf P`` relates conductor voltages to per-unit-length charges ``\mathbf q``:

```math
\mathbf V=\mathbf P\mathbf q,
\qquad
\mathbf C=\mathbf P^{-1},
\qquad
\mathbf Y=j\omega\mathbf P^{-1}.
```

The units of ``\mathbf P`` and ``\mathbf C`` are ``\mathrm m/\mathrm F`` and ``\mathrm F/\mathrm m``, respectively. The admittance calculation requires the inverse of the complete potential-coefficient matrix [Ametani1980, Ametani2015b](@cite). Matrix inversion couples the coefficients, so ``Y_{jk}`` cannot be calculated as ``j\omega/P_{jk}``.

Lossy formulations may use complex-permittivity potential coefficients or provide the shunt response directly. Their normalization determines the conversion to ``\mathbf Y=\mathbf G+j\omega\mathbf C`` [Gustavsen2005](@cite). Under the charge-based, lossless definition above, ``\mathbf P=j\omega\mathbf Y^{-1}`` for invertible matrices and ``\omega\ne0``. The potential-coefficient matrix is therefore distinct from an inverse-admittance matrix.

## Physical contributions and formulation families

The formulation families correspond to the physical contributions required by the impedance and admittance calculations.

| Physical contribution | Matrix contribution | Formulation section |
| --- | --- | --- |
| Electromagnetic response inside conducting materials | Local series response, including metallic surface and through-wall transfer impedances | [Internal impedance](internal_impedance.md) |
| Fields in insulation, jackets, and semiconducting regions | Series-field contributions and local potential-coefficient or shunt contributions | [Insulation parameters](insulation_parameters.md) |
| Series fields outside the local cable assemblies | External self and mutual impedances, including earth return | [External impedance and earth return](earth_return_impedance.md) |
| Electric fields in the external medium | External potential coefficients used in the shunt-admittance calculation | [External admittance and earth return](earth_return_admittance.md) |
| Earth conductivity, permittivity, and their frequency dependence | Material inputs to the external-field formulations | [Earth properties](earth_properties.md) |

Each contribution has a selectable formulation. Formulations used together must have compatible geometry, material properties, voltage references, current conventions, and propagation assumptions.

A formulation may provide a scalar coefficient, a coupled block, or a complete matrix. Its output must correspond to the physical contribution and conductor representation required by the assembly. An earth-return correction requires the associated geometric term to form a complete external impedance. Radial branch admittances and loop impedances require conversion before use in a conductor matrix. Contributions already included in one term must not be repeated in another.

## Block matrix assembly

In the cylindrical construction considered here, each local unit consists of one round conductor or a concentric assembly of metallic conductors and dielectric regions. The shared regions comprise the air, earth, or pipe cavity outside the individual units. The following arbitrary-dimension notation reformulates the block patterns in [Ametani1980, Ametani2015b](@cite). The cylindrical assumptions apply to this decomposition, rather than to the transmission-line matrix formulation itself.

Consider ``m`` local units. Unit ``j`` contains ``N_j`` longitudinal metallic conductors, ordered from the center outward when the unit is coaxial. The conductor count, excluding a possible common pipe, is

```math
n=\sum_{j=1}^{m}N_j.
```

A coaxial unit may contain a solid central conductor followed by any number of concentric annular conductors, separated by dielectric regions. The values of ``N_j`` need not be identical across the installation.

Let ``\mathbf B_j`` and ``\mathbf D_j`` denote the local series-impedance and potential-coefficient blocks of unit ``j``, respectively. Both have order ``N_j``. Define

```math
\mathbf B=\operatorname{blockdiag}(\mathbf B_1,\ldots,\mathbf B_m),
\qquad
\mathbf D=\operatorname{blockdiag}(\mathbf D_1,\ldots,\mathbf D_m).
```

Each local block contains the self and mutual coefficients of the conductors within its unit. The block-diagonal structure excludes interactions between different units from ``\mathbf B`` and ``\mathbf D``.

### Current and charge aggregation

Define the conductor-membership matrix

```math
\mathbf S=\operatorname{blockdiag}
\left(\mathbf 1_{N_1},\ldots,\mathbf 1_{N_m}\right)
\in\mathbb R^{n\times m},
```

where ``\mathbf 1_{N_j}`` is a column of ones. Column ``j`` identifies the conductors belonging to unit ``j``. The total longitudinal current of each unit is

```math
\mathbf J=\mathbf S^{\mathsf T}\mathbf I,
\qquad
J_j=\sum_{k=1}^{N_j}I_{j,k}.
```

The corresponding total per-unit-length charge is

```math
\mathbf Q=\mathbf S^{\mathsf T}\mathbf q,
\qquad
Q_j=\sum_{k=1}^{N_j}q_{j,k}.
```

Under the cylindrical approximation, the shared-region formulations use ``J_j`` and ``Q_j`` as source amplitudes. Multiplication by ``\mathbf S`` assigns the resulting voltage gradient or potential contribution to each conductor in the corresponding unit.

For a shared-region impedance matrix ``\mathbf H_Z`` and potential-coefficient matrix ``\mathbf H_P``, the conductor-matrix contributions are

```math
\Delta\mathbf Z=\mathbf S\mathbf H_Z\mathbf S^{\mathsf T},
\qquad
\Delta\mathbf P=\mathbf S\mathbf H_P\mathbf S^{\mathsf T}.
```

The same membership matrix applies to both quantities. For example, the impedance block between units ``j`` and ``k`` is

```math
\left(\mathbf S\mathbf H_Z\mathbf S^{\mathsf T}\right)_{jk}
=(H_Z)_{jk}\mathbf 1_{N_j}\mathbf 1_{N_k}^{\mathsf T}.
```

The repeated coefficients result from representing each unit by its aggregate source amplitude [Ametani1980, Ametani2015b](@cite).

Angular current or charge distributions omitted by the cylindrical approximation require a more detailed coupled formulation. When that formulation couples several units, the coupled response must remain in the combined matrix rather than in independent local blocks.

### Assembly without a common pipe

Let ``\mathbf Z_{\mathrm{env}}`` and ``\mathbf P_{\mathrm{env}}`` describe the shared medium outside the individual cable boundaries. The complete matrices are

```math
\boxed{
\begin{aligned}
\mathbf Z&=\mathbf B+\mathbf S\mathbf Z_{\mathrm{env}}\mathbf S^{\mathsf T},\\
\mathbf P&=\mathbf D+\mathbf S\mathbf P_{\mathrm{env}}\mathbf S^{\mathsf T}.
\end{aligned}
}
```

The local blocks contain the contributions within each unit. The external matrices contain the self and mutual contributions of the surrounding medium.

Under the lossless potential-coefficient convention,

```math
\boxed{
\mathbf Y=j\omega
\left(\mathbf D+\mathbf S\mathbf P_{\mathrm{env}}\mathbf S^{\mathsf T}\right)^{-1}.
}
```

Inversion applies to the full sum. In general, independently inverting the local and external potential blocks and adding their admittances does not recover ``\mathbf Y``.

## Local impedance and potential-coefficient blocks

### Metallic surface and transfer impedances

The [internal-impedance formulations](internal_impedance.md) describe electromagnetic diffusion, loss, and magnetic energy inside conducting materials. For a solid core and concentric annular metals, the local series block contains the following quantities.

| Quantity | Symbol | Physical role |
| --- | --- | --- |
| Solid-core outer-surface impedance | ``z_1^{\mathrm o}`` | Internal response of the solid conductor at its outer cylindrical surface |
| Annular inner-surface impedance | ``z_k^{\mathrm i}`` | Response of the annular metal at its inner cylindrical surface |
| Annular outer-surface impedance | ``z_k^{\mathrm o}`` | Response of the annular metal at its outer cylindrical surface |
| Annular transfer impedance | ``z_k^{\mathrm t}`` | Coupling between the inner and outer surface responses through the metal wall |

For annular conductors, ``k=2,\ldots,N_j``. The solid core has a disk cross-section. All quantities in the table are per-unit-length impedances in ``\Omega/\mathrm m``, rather than unnormalized electric-to-magnetic-field ratios [Schelkunoff1934, Ametani1980, Ametani2015b](@cite).

Inner- and outer-surface impedances both describe the response inside the metal. The surface designation specifies the boundary at which that response is represented. It does not assign an outer-surface impedance to the surrounding air or earth.

The surface and transfer impedances describe one coupled annular-metal problem. They are not independent series elements. Transfer impedance represents coupling through one conductor wall, distinct from mutual impedance between separate cables. Its sign and position in ``\mathbf B_j`` follow the conversion from surface quantities to conductor voltages and currents.

The conventional radial conductor solution assumes concentric homogeneous metals and linear materials. The conductor-diffusion approximation neglects displacement current in the metal. It retains skin effect but does not resolve arbitrary proximity-induced angular redistribution [Schelkunoff1934, Ametani2015b](@cite).

A skin- and proximity-effect formulation can provide a coupled conductor response [Patel2014](@cite). That response must retain its coupling and exclude shared-region contributions already present elsewhere in the assembly.

### Insulation and jackets

The local series block ``\mathbf B_j`` contains both metallic impedances and magnetic-field contributions from dielectric gaps and jackets. Its nonmetal contributions extend to the chosen outer cable boundary [Ametani1980, Ametani2015b](@cite).

Electric-field storage and dielectric leakage in the same regions contribute to the shunt problem. Their potential-coefficient representation enters ``\mathbf D_j``. The [insulation-parameter formulations](insulation_parameters.md) describe the series and shunt contributions, including semiconducting regions and material losses.

Radial branch or loop quantities require conversion to the conductor representation before assembly into the local blocks. This requirement applies to individual annuli, combined insulation paths, and multilayer radial networks [Weeks1984, Ametani2004, Ghosh2022](@cite).

The local and external calculations must use the same separating surface. When ``\mathbf B_j`` includes the outer jacket, the external series calculation begins at the jacket's outer radius. Beginning at the outer metal surface would include the jacket field twice. The potential-coefficient calculation requires the same consistency for its dielectric regions.

### Two- and three-conductor coaxial units

Each additional concentric metallic layer introduces an annular-metal response and the intervening dielectric region, increasing the order of the local blocks.

For a core/sheath unit, the conductor ordering is ``(c,s)``. For a core/sheath/armor unit, it is ``(c,s,a)``. Their series matrices are

```math
\mathbf Z_{cs}
=\mathbf B_{cs}^{(2)}
+z_{\mathrm{env}}\mathbf 1_2\mathbf 1_2^{\mathsf T},
\qquad
\mathbf Z_{csa}
=\mathbf B_{csa}^{(3)}
+z_{\mathrm{env}}\mathbf 1_3\mathbf 1_3^{\mathsf T}.
```

Each external coefficient corresponds to the outer geometry of its construction. The local blocks contain the metallic responses, dielectric gaps, and outer jacket. Both expressions are particular cases of

```math
\mathbf Z^{(N)}
=\mathbf B^{(N)}
+z_{\mathrm{env}}\mathbf 1_N\mathbf 1_N^{\mathsf T}.
```

A round line conductor corresponds to ``N=1``. The relative positions of several such conductors enter the external matrix. Cylindrical symmetry is assumed locally and does not require the complete installation to be coaxial [Ametani1980, Ametani2015b](@cite).

## External fields and earth properties

The [external-impedance formulations](earth_return_impedance.md) describe the series response outside the local cable boundaries. The [external-admittance formulations](earth_return_admittance.md) describe the corresponding electric-field response, commonly through potential coefficients. Where the external medium includes earth, both calculations use the selected [earth-property model](earth_properties.md).

The conductor arrangement determines the required self and mutual coefficients. Overhead, buried, and mixed installations require formulations for the corresponding conductor pairs. In a mixed installation, the overhead–buried mutual coefficients belong to the same external matrices as the overhead–overhead and buried–buried coefficients.

The selected formulations specify the treatment of stratification, permeability, displacement current, and longitudinal propagation. Under the cylindrical construction above, these choices determine the external coefficients while ``\mathbf S`` specifies conductor membership.

An imperfect-earth impedance correction must be combined with its compatible geometric contribution to form a complete entry of ``\mathbf Z_{\mathrm{env}}``. External potential coefficients enter ``\mathbf P_{\mathrm{env}}`` before inversion of the complete potential matrix.

The voltage reference is part of the external shunt formulation. The conventional electrostatic treatment of directly buried cables assumes equipotential soil. Under this assumption, the separate exterior potential-coefficient term vanishes, while the insulation and jacket contributions remain [Ametani1980, Ametani2015b](@cite).

A finite-conductivity-earth formulation can retain the external potential coefficients instead of imposing an equipotential soil. Neither shunt treatment removes the earth-return impedance from ``\mathbf Z``.

## Common metallic pipe

A metallic pipe enclosing several units introduces a common cavity, a finite conducting wall, and an exterior region. An outer insulating jacket may separate the pipe metal from the surrounding medium.

Retain the pipe as an additional conductor, with ordering

```math
\widetilde{\mathbf I}
=\begin{bmatrix}\mathbf I\\I_p\end{bmatrix},
\qquad
\widetilde{\mathbf V}
=\begin{bmatrix}\mathbf V\\V_p\end{bmatrix}.
```

The complete matrices have order ``n+1=1+\sum_jN_j``. For a single enclosed unit with ``N`` conductors, the order is ``N+1``.

### Series-impedance assembly

The series decomposition is

```math
\mathbf Z=\mathbf Z_i+\mathbf Z_p+\mathbf Z_c+\mathbf Z_0.
```

The four terms describe the individual units, pipe-interior response, remaining pipe-wall and jacket contributions, and exterior return response [Ametani1980, Ametani2015b](@cite).

Let ``\mathbf H_p\in\mathbb C^{m\times m}`` describe the pipe-interior impedances relative to the pipe inner surface. Then

```math
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
```

The matrix ``\mathbf H_p`` contains the cavity self and mutual impedances together with the pipe's inner-surface metallic response. Its coefficients depend on the positions and radii of the enclosed units and on the selected pipe formulation.

A formulation using a different partition requires conversion to this decomposition before assembly. In particular, a complete core–pipe loop impedance may already contain a conductor contribution assigned to ``\mathbf B``. That contribution must be excluded from ``\mathbf H_p``.

For the remaining terms, define

```math
a_p=z_p^{\mathrm o}+z_{p,\mathrm{jacket}},
\qquad
t_p=z_p^{\mathrm t},
\qquad
\mathbf u=\mathbf 1_n,
\qquad
\mathbf v=\begin{bmatrix}\mathbf u\\1\end{bmatrix}.
```

The quantities ``z_p^{\mathrm o}`` and ``z_p^{\mathrm t}`` are the pipe-metal outer-surface and through-wall transfer impedances. The term ``z_{p,\mathrm{jacket}}`` describes the series field in the outer jacket, when present. The pipe inner-surface response is included in ``\mathbf H_p``.

With the common longitudinal-current convention,

```math
\mathbf Z_c=
\begin{bmatrix}
(a_p-2t_p)\mathbf u\mathbf u^{\mathsf T}&(a_p-t_p)\mathbf u\\
(a_p-t_p)\mathbf u^{\mathsf T}&a_p
\end{bmatrix},
\qquad
\mathbf Z_0=z_e\mathbf v\mathbf v^{\mathsf T}.
```

The scalar ``z_e`` is the pipe's external self/return impedance, evaluated outside its jacket surface. The corresponding source current is ``\mathbf u^{\mathsf T}\mathbf I+I_p``. The external-field calculation excludes the pipe-metal and jacket contributions already assigned to the other terms [Ametani2015b](@cite).

With ``g_p=a_p+z_e``, the complete matrix is

```math
\boxed{
\mathbf Z=
\begin{bmatrix}
\mathbf B+\mathbf S\mathbf H_p\mathbf S^{\mathsf T}
+(g_p-2t_p)\mathbf u\mathbf u^{\mathsf T}
&(g_p-t_p)\mathbf u\\
(g_p-t_p)\mathbf u^{\mathsf T}&g_p
\end{bmatrix}.
}
```

The upper-left block contains the enclosed-conductor coefficients. Its cross-unit entries consist of the cavity mutual impedance and the common pipe/exterior contribution. The last row and column contain the pipe self impedance and its mutual impedances with the enclosed conductors. No pipe-return current has been prescribed.

### Potential-coefficient assembly

Let ``\mathbf H_{P,p}`` describe the cavity potential coefficients relative to the pipe inner surface. Let ``p_{\mathrm{ext}}`` contain the pipe-jacket and exterior-space potential coefficients. With the same conductor ordering,

```math
\boxed{
\mathbf P=
\begin{bmatrix}
\mathbf D+\mathbf S\mathbf H_{P,p}\mathbf S^{\mathsf T}&\mathbf 0\\
\mathbf 0^{\mathsf T}&0
\end{bmatrix}
+p_{\mathrm{ext}}\mathbf v\mathbf v^{\mathsf T}.
}
```

The local dielectric, cavity, and exterior contributions are assembled before inversion [Ametani2015b](@cite). The pipe-jacket term represents the electric field in the outer insulation. It does not represent capacitance through the conducting pipe wall.

The voltage reference must leave independent electrical variables for inversion. When an ideal reference constraint makes the displayed matrix singular, the potential problem must first be expressed in independent voltage variables.

### Pipe-model assumptions

The cavity kernel and pipe surface/transfer impedances must use consistent finite-wall or thick-wall assumptions. Retaining the pipe as an explicit conductor does not alter the approximations in these quantities. Incompatible wall assumptions can invalidate their low-frequency combination [Ametani2015b, Hoidalen2013](@cite).

A cavity model that includes eccentric positions does not necessarily resolve proximity-induced current redistribution within the enclosed conductors. A proximity correction must contain only the response absent from its compatible base formulation, excluding the base skin-effect contribution [Hoidalen2025](@cite).

## Discrete and equivalent conductor representations

Stranded cores, wire screens, and armor can be represented by individual wires or by equivalent continuous conductors. The representation determines which conductor voltages, currents, and geometric details are retained.

Circular wires inside a common circular enclosure can be treated as separate units with ``N_j=1``. The cavity matrix contains their self and mutual interactions under the selected formulation. Each wire retains its own voltage and current, and the enclosure retains a separate conductor row and column.

Equivalent conductors preserve selected aggregate properties. A solid equivalent may reproduce a stranded core's metallic fill or specified direct-current resistance. An equivalent annulus may reproduce a wire screen's metallic area and radial position. Effective dielectric data may preserve capacitance when semiconducting layers are not represented separately [Gustavsen2001, Gustavsen2005](@cite).

Matching metallic area or direct-current resistance alone does not establish equivalence of skin, proximity, helical, magnetic, or transfer-impedance behavior. Replacing magnetic wire armor by an equal-area tube requires additional effective-material assumptions [Gustavsen2001, Gustavsen2005](@cite).

The [internal-impedance formulations](internal_impedance.md) include equivalent-conductor models and coupled cross-sectional field formulations. A coupled solution must retain the matrix entries that describe its interactions. It cannot be substituted directly for a scalar coaxial impedance.

## Circuit constraints and modal decomposition

The assembled matrices retain each conductor's electromagnetic contribution before the application of return-current or terminal constraints. Grounded sheaths and bonded armor therefore remain in the parameter calculation. A reduced circuit representation requires the corresponding electrical constraints.

### Pipe-return loop impedance

For the explicit-pipe assembly, the conductor-to-pipe mutual block is

```math
\mathbf Z_{\mathrm{inner},p}=(g_p-t_p)\mathbf u.
```

For outgoing current in conductor ``\alpha`` and return current in the pipe, the loop impedance is

```math
z_{\alpha\text{–}p}^{\mathrm{loop}}
=Z_{\alpha\alpha}+Z_{pp}-Z_{\alpha p}-Z_{p\alpha}.
```

For all pipe-return loops, define

```math
\mathbf K=
\begin{bmatrix}
\mathbf I_n\\-\mathbf u^{\mathsf T}
\end{bmatrix},
\qquad
\widetilde{\mathbf I}=\mathbf K\mathbf i,
\qquad
\mathbf w=\mathbf K^{\mathsf T}\widetilde{\mathbf V}
=\mathbf V-\mathbf uV_p,
```

where ``\mathbf I_n`` is the identity matrix. The loop currents ``\mathbf i`` impose ``I_p=-\mathbf u^{\mathsf T}\mathbf i``, and ``\mathbf w`` contains the conductor-to-pipe voltage differences. The transformed impedance matrix is

```math
\boxed{
\mathbf Z_{\mathrm{loop},p}
=\mathbf K^{\mathsf T}\mathbf Z\mathbf K
=\mathbf B+\mathbf S\mathbf H_p\mathbf S^{\mathsf T}.
}
```

For this balanced pipe-return excitation, the common exterior and remaining pipe-wall terms cancel. The pipe's inner-surface impedance remains in ``\mathbf H_p``. The result retains finite pipe conductivity and depends on the imposed return-current distribution. Terminal grounding alone does not impose that distribution.

For one enclosed coaxial unit, ``\mathbf H_p`` reduces to a scalar ``h_p``. The core/sheath and core/sheath/armor loop matrices are

```math
\mathbf Z_{cs,p}^{\mathrm{loop}}
=\mathbf B_{cs}^{(2)}+h_p\mathbf 1_2\mathbf 1_2^{\mathsf T},
\qquad
\mathbf Z_{csa,p}^{\mathrm{loop}}
=\mathbf B_{csa}^{(3)}+h_p\mathbf 1_3\mathbf 1_3^{\mathsf T}.
```

The corresponding externally referenced matrices have orders three and four. The loop matrices have lower order because the pipe current is constrained and the voltages are referred to the pipe.

### Propagation equations

For a longitudinally uniform section, differentiating the transmission-line equations gives

```math
\frac{\mathrm d^2\mathbf V}{\mathrm dx^2}
=\mathbf Z\mathbf Y\mathbf V,
\qquad
\frac{\mathrm d^2\mathbf I}{\mathrm dx^2}
=\mathbf Y\mathbf Z\mathbf I.
```

The products ``\mathbf Z\mathbf Y`` and ``\mathbf Y\mathbf Z`` define the voltage and current propagation problems, respectively. Internal, insulation, and external formulations affect propagation through their contributions to these matrices. [Modal decomposition](modal_decomposition.md) develops the corresponding modal representations and propagation quantities.
