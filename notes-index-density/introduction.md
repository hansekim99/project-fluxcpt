---
header-includes:
    - \usepackage{physics}
    - \usepackage{amsmath}
    - \usepackage{amssymb}
    - \usepackage{mathtools}
    - \usepackage{slashed}
    - \usepackage{tikz-feynman}
    - \usepackage{tikz}
    - \usetikzlibrary{knots,decorations.markings}
    - \usepackage{graphicx}
    - \usepackage{subcaption}
    - \usepackage{simpler-wick}
    - \allowdisplaybreaks
---

# Introduction

In this work, we calculate lower bounds for the supersymmetric vacua in flux compactifications of Type IIB string theories on CY threefolds in the KS database. For ${ h^{1,2} \leq 5 }$ (goal?) the index density is analytically computed (in the continuous flux approximation) using an efficient combinatorical method; in numerical integration over the Kahler moduli space, we define regions of control (wording?) that are explicitly calculated, and apply this analysis to chambers of extended Kahler cones (of ... ?). For ${ 6 \leq h^{1,2} \leq ? }$, we (? : codimension 1 singularities? asymptotic Hodge theory?)

\paragraph{Introduction and literature review}

Notes

- 2501.03984 p1 : "To date, progress in exploring the string landscape has largely relied on searches in hand-selected examples or on statistical arguments employing suitable approximations."

\paragraph{Summary}

Monte Carlo integration of the index density over the Kahler moduli space is carried out for low ${ h^{1,2} }$ where it is computationally viable. The index density is combinatorically...

## Type IIb Compactifications

In this section a brief introduction to compactifications of type IIb string theories, alongside a reference to the notation is given.  

Let ${ X_{3} }$ refer to a Calabi-Yau threefold, upon which type IIb string theory is compactified, and its mirror dual ${ \tilde{X}_{3} }$. Hodge numbers are used in reference to ${ X_{3} }$; i.e. ${ h^{1,1} = h^{1,1}(X_{3}) = h^{1,2}(\tilde{X}_{3}), h^{1,2} = h^{1,2}(X_{3}) = h^{1,1}(\tilde{X}_{3}) }$.

- To be precise, compactify on orientifold ${ X_{3}/\mathcal{I} }$; 2501.03984 p3

On ${ X_{3} }$, we introduce a symplectic basis of cycles ${ \{\Sigma _{I}, \Sigma ^{I}\} \subset H_{3}(X_{3}, \mathbb{Z} ) }$ and the Poincare dual forms ${ \{\alpha ^{I}, \beta _{I}\} }$. Then we define the periods by integrating the holomoprhic 3-form ${ \Omega  }$ over these cycles, 
$${ Z^{I} = \int _{\Sigma _{I}}\Omega = \int _{X_{3}}\Omega \wedge \alpha ^{I}, \mathcal{F}_{I} = \int _{\Sigma ^{I}}\Omega = \int _{X_{3}} \Omega \wedge \beta _{I}. }$$

The periods ${ Z^{I} }$, ${ I = 0,\cdots ,h^{1,2} }$ define local homogeneous coordinates in the complex structure moduli space of ${ X_{3} }$; then define projective coordinates ${ z^{i} = Z^{i}/Z^{0} }$, ${ i = 1, \cdots ,h^{1,2} }$.

Mirror symmetry maps the type IIb theory on ${ X_{3} }$ to a type IIa theory on ${ \tilde{X}_{3} }$, where the large complex structure region (of the complex structure moduli space) of ${ X_{3} }$ is mapped to the large volume region (of the Kahler moduli space) of ${ \tilde{X}_{3} }$. In the large volume region, the prepotential can be expanded as

$${ \dfrac{F(z)}{(Z^{0})^{2} } = -\dfrac{1}{3!} \kappa _{ijk}z^{i}z^{j}z^{k} + \dfrac{1}{2}a_{ij} z^{i}z^{j} + b_{i}z^{i} + \xi + F_{\text{inst}} , }$$

where ${ \kappa _{ijk} }$ are the intersection numbers of ${ \tilde{X}_{3} }$, the parameters depend on the (1,1)-forms of ${ \tilde{X}_{3} }$ and the second Chern class ${ c_{2}(\tilde{X}_{3}) }$ as

$${ \begin{aligned}
\kappa _{ijk} &= \int _{\tilde{X}_{3}} J_{i} \wedge J_{j}\wedge J_{k}, & a_{ij} &= \dfrac{1}{2}\int _{\tilde{X}_{3}}J_{i} \wedge J_{j} \wedge J_{k}  \\
b_{i} &= \dfrac{1}{4!} \int _{\tilde{X}_{3}} c_{2}(\tilde{X}_{3}) \wedge J_{i}, & \xi &= \dfrac{i}{2} \dfrac{\zeta (3)\chi (\tilde{X}_{3})}{(2\pi )^{3} },
\end{aligned}
 }$$

and the nonperturbative term ${ F_{\text{inst}} }$ arising from the worldsheet instanton effects is given as

$${ F_{\text{inst}}(z) = -\dfrac{1}{(2\pi i)^{3} } \sum _{q \in \mathcal{M}(\tilde{X}_{3})} n_{q}^{0} \text{Li}_{3}(e^{2\pi i \ev{q, z}}). }$$

The sum runs over the effective curves (charges) in the Mori cone of the dual threefold ${ \mathcal{M}(\tilde{X}_{3}) }$. The type IIa effective bosonic action is characterised by the Kahler metric and the gauge kinetic metric,

$${ \begin{aligned}
K(z) &= - \log i( \bar{z} ^{I} F_{I} - z^{I} \bar{F}_{I} ) \\
g_{i\overline{j}}(z) &= \partial _{z^{i}} \bar{\partial} _{z^{j}} K\\
\end{aligned}
}$$

where ${ F_{I} = \partial _{Z^{I}} F }$ and the holomorphic 3-form ${ \Omega  }$ is fixed by the choice ${ Z^{0} = 1 }$ (in which ${ F_{I} }$ corresponds to the dual cycles ${ \mathcal{F}_{I} }$ and satisfies ${ F_{I} = 2F - z^{i}\partial _{z^{i}}F }$). The Kahler metric is independent of the real components of the moduli, and can be rewritten in terms of ${ t^{i} = \text{Im} z^{i} }$ and ${ \partial _{i} = \partial _{t^{i}} }$; 

<!-- $${ \begin{aligned} -->
<!-- F &=  - \dfrac{n}{(2\pi i)^{3} } (Z^{0})^{2}  \text{Li}_{3}(e^{2\pi i q_{i}Z^{i} / Z^{0}}) = - \dfrac{n}{(2\pi i)^{3} }(Z^{0})^{2}  \text{Li}_{3}(w)\\ -->
<!-- \partial _{Z^{i}}F &= -\dfrac{n}{(2\pi i)^{3} } (Z^{0})^{2} \dfrac{\partial w}{\partial Z^{i}} \partial _{w } \text{Li}_{3} (w) = -\dfrac{n}{(2\pi i)^{3} } (Z^{0})^{2} \dfrac{2\pi i q_{i}}{Z^{0}} w \partial _{w} \text{Li}_{3}(w) = \dfrac{n}{(2\pi )^{2} } Z^{0} q_{i} \text{Li}_{2}(w)\\ -->
<!-- \partial _{Z^{0} } F &= -\dfrac{n}{(2\pi i)^{3} } 2 Z^{0} \text{Li}_{3}(w) - \dfrac{n}{(2\pi i)^{3} } (Z^{0})^{2} \dfrac{\partial w}{\partial Z^{0}} \partial _{w} \text{Li}_{3}(w) \\ -->
<!-- &= \dfrac{2}{Z^{0}}F - \dfrac{n}{(2\pi i)^{3} } (Z^{0})^{2} \dfrac{-2\pi i q_{i}Z^{i}}{(Z^{0})^{2} }w \partial _{w} \text{Li}_{3}(w)= \dfrac{2}{Z^{0}}F- \dfrac{n}{(2\pi )^{2} }  q_{i}Z^{i} \text{Li}_{2}(w) \\ -->
<!-- \overline{Z}^{I} F_{I}|_{Z^{0}=1} &= \dfrac{n}{(2\pi )^{2} } q_{i}\bar{z}^{i} \text{Li}_{2}(w) + 2F - \dfrac{n}{(2\pi )^{2} } q_{i}z^{i} \text{Li}_{2}(w) = 2 F - 2i\dfrac{n}{(2\pi )^{2} } q_{i}t^{i} \text{Li}_{2}(w)\\ -->
<!-- (\bar{Z}^{I}F_{I} - Z^{I}\bar{F}_{I})|_{Z^{0}=1} &= 4i \text{Im}F - \dfrac{4in}{(2\pi )^{2} } q_{i}t^{i} \text{Re} (\text{Li}_{2}(w)) = -\dfrac{4in}{(2\pi )^{3} }\text{Re}(\text{Li}_{3}(w)) - \dfrac{4in}{(2\pi )^{2} }q_{i}t^{i} \text{Re}(\text{Li}_{2}(w)) \\ -->
<!-- &= -4in \left[\dfrac{1}{(2\pi )^{3} }\text{Li}_{3}(e^{-2\pi q_{i}t^{i}})\right] -->
<!-- \end{aligned} -->
<!--  }$$ -->

$${ \begin{aligned}
K(t) &= -\log \dfrac{4}{3} \kappa _{ijk} t^{i}t^{j}t^{k} = - \log \kappa \\
g_{i\overline{j}}(t) &= \dfrac{\partial _{i} \kappa \partial _{j} \kappa - \kappa  \partial _{i}\partial _{j} \kappa }{4\kappa ^2 }\\
y_{ijk} &= \partial _{i}\partial _{j}\partial _{k} F = \kappa _{ijk} + \sum _{q \in \mathcal{M}(\tilde{X}_{3})} n^{0}_{q} q_{i}q_{j}q_{k} \text{Li}_{0} (e^{-2\pi q_{i}t^{i}})
\end{aligned}
 }$$

Note that all ${ a_{ij}, b_{i} }$ dependence drops out in ${ K(z) }$. The perturbative 
The metric is positive-definite only in the Kahler cone of the dual threefold, ${ \mathcal{K}(\tilde{X}_{3}) }$, which is dual to the Mori cone. This corresponds to the region where the epxansion of the prepotential is valid. Within this region, the Cholesky decomposition of the inverse metric gives the orthonormal frame ${ g^{i \bar{j}} = e^{i}_{I}\delta ^{I \bar{J}} e^{\bar{j}}_{\bar{J}} }$ and rescaled Yukawa couplings

$${ \mathcal{F}_{IJK} = -i (e^{i}_{I})(e^{j}_{K})(e^{k}_{K}) y_{ijk} e^{K} }$$

where indices ${ I,J,K = 1,\cdots ,h^{1,2} }$ are used distinctively from the period basis above.

