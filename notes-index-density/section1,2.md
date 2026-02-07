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
# Research Notes : Index Density

In analytically evalauting ${ \rho (z) }$, 
$${ \begin{aligned}
\rho _{\text{ind}}(z)&=_{(1)} \pi ^{-2(n+1)}\int \mathrm{d} ^2 X\ \mathrm{d} ^{2n} \phi\  e^{-|X|^{2} -|\phi |^{2} } |X|^{2} \det \begin{pmatrix}
{\delta^{I}}_{J}\overline{X} - \overline{\phi } ^{I} \phi _{J} / X & \mathcal{F}_{IJK} \overline{\phi }^{K}\\ \overline{\mathcal{F}}^{IJK} \phi _K & {\delta _{I}}^{J}X - \phi _{I} \overline{\phi }^{J} / \overline{X}
\end{pmatrix}\\
&=_{(2)} \pi ^{-2(n+1)}\int \mathrm{d} ^2 X\ \mathrm{d} ^{2n} \phi\  e^{-|X|^{2} -|\phi |^{2} } |X|^{2} \det ({\delta ^{I}}_{J} \overline{X}- \overline{\phi }^{I}\phi _{J} /X ) \det ({\delta _{I}}^{J}X - \phi _{I}\overline{\phi }^{J} / \overline{X}) \\
&\phantom{=_{(3)} \pi ^{-2(n+1)}\int \mathrm{d} ^2 \mathrm{d} ^{2n} \phi\  e^{-|X|^{2} -|\phi |^{2} }}\times \det ({\delta _{I}}^{P} - |X|^{-2}  (\mathcal{F}_{IJK} \overline{\phi }^{K})({\delta _{J}}^{L} - \phi_{J} \overline{\phi }^{L} / |X|^{2} )^{-1} (\overline{\mathcal{F}}^{LMN}\phi _{N})({\delta ^{M}}_{P} - \overline{\phi }^{M}\phi _{P}/|X|^{2} )^{-1} )\\
&=_{(3)} \pi ^{-2(n+1)}\int \mathrm{d} ^2 X\ \mathrm{d} ^{2n} \phi\  e^{-|X|^{2} -|\phi |^{2} } |X|^{2n-2} (|X|^{2}  - |\phi |^{2} )^{2} \\
&\phantom{=_{(3)} \pi ^{-2(n+1)}\int \mathrm{d} ^2 X }\times \det \left({\delta _{I}}^{P} - |X|^{-2} (\mathcal{F}_{IJK} \overline{\phi }^{K})\left({\delta ^{J}}_{L} + \dfrac{\overline{\phi } ^{J}\phi _{L}} {|X|^{2}  - |\phi |^{2} }\right) (\overline{\mathcal{F}}^{LMN}\phi _{N})\left({\delta _{M}}^{P} + \dfrac{\phi _{M}\overline{\phi }^{P}}{|X|^{2}  - |\phi |^{2}} \right)\right)\\
&=_{(4)} \pi ^{-2(n+1)}\int \mathrm{d} ^2 X\ \mathrm{d} ^{2n} \phi\  e^{-|X|^{2} -|\phi |^{2} } |X|^{2n-2} D^4\det \biggr({\delta _{I}}^{P} -   \mathcal{F}_{IJK} \bar{\mathcal{F}}^{JPN} \overline{\phi }^{K} \phi _{N}|X|^{-2} -  \mathcal{F}_{IJK} \bar{\mathcal{F}}^{JMN} \overline{\phi }^{K} \overline{\phi }^{P} \phi _{M} \phi _{N} |X|^{-2} / D^2  \\
&\phantom{=_{(3)} \pi ^{-2(n+1)}} - \mathcal{F}_{IJK} \bar{\mathcal{F}}^{LPN} \overline{\phi }^{J} \overline{\phi }^{K} \phi _{L} \phi _{N} |X|^{-2} / D^2 - \mathcal{F}_{IJK} \bar{\mathcal{F}}^{LMN} \overline{\phi }^{J}\overline{\phi }^{K}\overline{\phi }^{P} \phi _{L}\phi _{M}\phi _{N} |X|^{-2} / D^4\biggr)\\
&=_{(5)} \pi ^{-2(n+1)}\int \mathrm{d} ^2 X\ \mathrm{d} ^{2n} \phi\  e^{-|X|^{2} -|\phi |^{2} } |X|^{2n-2} D^4 \exp \tr \log (\mathbb{I}- |X|^{-2}N)\\
%&=_{(6)}(\pi/2) ^{-n-1} \int \mathrm{d} X \mathrm{d}^{n}\phi \cdot X \cdot \Pi \phi _{A}\cdot  e^{-X^2 -|\phi |^{2} } |X|^{-2n+2} D^4 \cdot (\text{contracted terms})
\end{aligned}
 }$$

- (1) : The variables ${ X }$ and ${ \phi _{I} }$ are dummy variables corresponding to a continuous limit of the flux quanta, where ${ Z }$ from Denef-Douglas is written as ${ \phi  }$ to prevent confusion with the complex moduli ${ z }$. Each are complex variables, integrated over the entire complex plane. 
    - In the orthogonal frame, raising and lowering indices is trivial. The convention here is that holomorphic variables in lower indices, and antiholo in upper. ${ |\phi |^{2} = \sum_{I}^{} \phi _{I} \overline{\phi }^{I}  }$.
- (2) : Evaluation of the determinant of the block matrix, from $${ \det \begin{pmatrix}
A&B\\C&D
\end{pmatrix}
= \det A \det (D - BA^{-1} C). }$$
- (3) : The determinant of a rank one perturbation to identity is given as $${ \det (\mathbb{I} - \overline{\phi }^{I}\phi _{J}/|X|^{2} ) = 1 - \phi _{I}\overline{\phi }^{I} / |X|^{2} = 1-|\phi |^{2} /|X|^{2} . }$$ In addition, the inverse of a rank one perturbation to identity is given by the Sherman-Morrison formula $${ (\mathbb{I} - \phi _{J}\overline{\phi }^{K} / |X|^{2} )^{-1}  = \mathbb{I} + (\overline{\phi }^{J}\phi _{K}) / (|X|^{2}  - |\phi |^{2} ). }$$
- (4) : ${ D^2 \coloneqq  |X|^{2}  - |\phi |^{2} }$.
- (5) : ${ \log \det M = \tr \log M }$. All non-identity terms are collected into the matrix ${ N }$.

<!--
- (6) : After transforming to polar coordinates (with real components ${ X,\phi _{I} }$), all angular integrals impose a contraction on the field indices as ${ \phi _{A} \overline{\phi }^{B} = e^{i(\alpha _{A}-\alpha _{B})} \text{Re}(\phi _{A}\overline{\phi }^{B}) }$ vanishes unless ${ A=B }$. Then, the angles are integrated out trivially.
-->

<!--
Two important simplifications come from analysing the form of the integral. 

- After factoring out the ${ |X|^{2} -\phi ^2  }$ terms in the denominators, the determinant is a polynomial involving ${ \mathcal{F}_{IJK}, \bar{\mathcal{F}}^{IJK}, \phi_{I} }$ and ${ X }$, with homog. degree in ${ \mathcal{F}_{IJK} }$ at ${ n }$. The entire integral is a Gaussian integral involving this polynomial, so will in itself be a ${ n }$-th degree homog. polynomial involving ${ \mathcal{F}_{IJK} }$ and ${ \bar{\mathcal{F}}^{IJK} }$. 
- The integral must be invariant with respect to ${ U(n) }$ transformations in the indices ${ IJK }$, so it can only be a function of the contracted terms ${ |\mathcal{F}|^{2} , \mathcal{F}_{IJK} \bar{\mathcal{F}}^{IJL} \mathcal{F}_{LMN} \bar{\mathcal{F}}^{KMN} , \cdots }$. 

Combining these two, we come to the conclusion that ${ \rho  }$ is a polynomial in terms of the ${ U(n) }$ invariant contracted terms. 
\textbf{Proposition. } The determinant can be exactly obtained by expanding ${ \exp \tr \log (\mathbb{I}-|X| ^{-2}N) }$ up to the ${ n }$-th degree, as the determinant is an ${ n }$-th deg. homog. polynomial and each ${ N }$ contributes exactly one degree of ${ \mathcal{F} }$. 
-->
$${ \begin{aligned}
N _{I} ^{ P} &= \mathcal{F}_{IJK} \bar{\mathcal{F}}^{JPN} \overline{\phi }^{K} \phi _{N}+  \mathcal{F}_{IJK} \bar{\mathcal{F}}^{JMN} \overline{\phi }^{K} \overline{\phi }^{P} \phi _{M} \phi _{N} / D^2  \\
& + \mathcal{F}_{IJK} \bar{\mathcal{F}}^{LPN} \overline{\phi }^{J} \overline{\phi }^{K} \phi _{L} \phi _{N} / D^2 + \mathcal{F}_{IJK} \bar{\mathcal{F}}^{LMN} \overline{\phi }^{J}\overline{\phi }^{K}\overline{\phi }^{P} \phi _{L}\phi _{M}\phi _{N}  / D^4
\end{aligned}
}$$

We note that in (4), each component of the matrix ${ N _{I}^{P} }$ inside the determinant is a degree 1 polynomial with respect to terms of the form ${ \mathcal{F}_{IJK}\bar{\mathcal{F}}^{LMN} }$. Therefore, the determinant itself will be a degree ${ n }$ polynomial wrt. ${ \mathcal{F}\bar{\mathcal{F}} }$, and the following proposition holds;
\paragraph{Proposition 1. } The integrand can be exactly obtained by expanding ${ \exp \tr \log (\mathbb{I}-|X|^{-2}N) }$ up to ${ n }$-th degree.

\newpage

$${ \begin{aligned}
N _{I} ^{ P} &= \mathcal{F}_{IJK} \bar{\mathcal{F}}^{JPN} \overline{\phi }^{K} \phi _{N}+  \mathcal{F}_{IJK} \bar{\mathcal{F}}^{JMN} \overline{\phi }^{K} \overline{\phi }^{P} \phi _{M} \phi _{N} / D^2  \\
& + \mathcal{F}_{IJK} \bar{\mathcal{F}}^{LPN} \overline{\phi }^{J} \overline{\phi }^{K} \phi _{L} \phi _{N} / D^2 + \mathcal{F}_{IJK} \bar{\mathcal{F}}^{LMN} \overline{\phi }^{J}\overline{\phi }^{K}\overline{\phi }^{P} \phi _{L}\phi _{M}\phi _{N}  / D^4
\end{aligned}
}$$

The objective is evaluating ${ \rho  }$ up to the ${ d }$-th degree in ${ N }$ for ${ d \leq n }$. This article is divided into three sections; the first evaluates ${ \rho  }$ for degrees ${ d=1 }$ and ${ d=2 }$ to demonstrate the general process, the second demonstrates a diagrammatic/combinatoric method, and the last provides proofs for the propositions in the article.

## ${ \boldsymbol{d=1} }$

\paragraph{Step 1} First, the expression ${ \exp \tr \log (\mathbb{I} - |X|^{-2} N) }$ is evaluated up to 1st order in ${ N }$. 
$${ \log (1- |X|^{-2}N) = - |X|^{-2}N - \dfrac{1}{2}|X|^{-4} N^2 -\cdots  }$$
$${ \begin{aligned}
\exp \tr \log (1-|X|^{-2}N) &= \exp \left(-|X|^{-2}\tr N - \dfrac{1}{2}|X|^{-4}\tr N^2 -\cdots \right)\\
&= 1 - |X|^{-2}\tr N + \cdots = 1 - |X|^{-2} N_I ^{I} + \cdots 
\end{aligned}
 }$$
The first order contribution to the integrand is

$${ \begin{aligned}
&e^{-|X|^{2}-|\phi |^{2}} |X|^{2n-2}& &\hspace{-.3cm}D^4 \cdot |X|^{-2} N _{I} ^{I} \\
=&e^{-|X|^{2}-|\phi |^{2}} |X|^{2n-2}& &\hspace{-.3cm}D^4 \cdot |X|^{-2} \biggr( \mathcal{F}_{IJK} \bar{\mathcal{F}}^{JIN} \overline{\phi }^{K} \phi _{N} + \mathcal{F}_{IJK} \bar{\mathcal{F}}^{JMN} \overline{\phi }^{K} \overline{\phi }^{I} \phi _{M} \phi _{N} / D^2  \\
& & & + \mathcal{F}_{IJK} \bar{\mathcal{F}}^{LIN} \overline{\phi }^{J} \overline{\phi }^{K} \phi _{L} \phi _{N} / D^2 + \mathcal{F}_{IJK} \bar{\mathcal{F}}^{LMN} \overline{\phi }^{J}\overline{\phi }^{K}\overline{\phi }^{I} \phi _{L}\phi _{M}\phi _{N}  / D^4 \biggr)\\
=&e^{-|X|^{2}-|\phi |^{2}} |X|^{2n-4}& &\hspace{-.3cm}\biggr( \mathcal{F}_{IJK} \bar{\mathcal{F}}^{JIN} \overline{\phi }^{K} \phi _{N} D^4 + \mathcal{F}_{IJK} \bar{\mathcal{F}}^{JMN} \overline{\phi }^{K} \overline{\phi }^{I} \phi _{M} \phi _{N} D^2  \\
& & & + \mathcal{F}_{IJK} \bar{\mathcal{F}}^{LIN} \overline{\phi }^{J} \overline{\phi }^{K} \phi _{L} \phi _{N} D^2 + \mathcal{F}_{IJK} \bar{\mathcal{F}}^{LMN} \overline{\phi }^{J}\overline{\phi }^{K}\overline{\phi }^{I} \phi _{L}\phi _{M}\phi _{N} \biggr)
\end{aligned}
}$$

\paragraph{Step 2} Next, integrate out all ${ \phi _{I}, \overline{\phi }^{J} }$. The Gaussian integral ${ \int \mathrm{d} ^{2n} \phi e^{-|\phi |^{2} } \phi _{I} \overline{\phi }^{J} }$ is only nonzero when ${ I=J }$, and equals ${ \delta _{I}^{J} }$. When there are multiple pairs of ${ \phi ,\overline{\phi } }$, the final result is given by Wick's theorem, and corresponds to all possible ways to contract over the given indices;
$${ \begin{aligned}
&\pi ^{-n}\int \mathrm{d} ^{2n} \phi e^{-|\phi |^{2} }\phi _{I}\overline{\phi }^{J} = \delta _{I}^{J}\\
&\pi ^{-n}\int \mathrm{d} ^{2n} \phi e^{-|\phi |^{2} }\phi _{I}\phi _{J} \overline{\phi }^{K}\overline{\phi }^{L} = \delta _{I}^{K}\delta _{J}^{L} + \delta _{I}^{L}\delta _{J}^{K}\\
&\pi ^{-n}\int \mathrm{d} ^{2n} \phi e^{-|\phi |^{2} }\phi _{I}\phi _{J}\phi _{K} \overline{\phi }^{L}\overline{\phi }^{M} \overline{\phi }^{N} = \delta _{I}^{L}\delta _{J}^{M} \delta _{K}^{N} + \delta _{I}^{L}\delta _{J}^{N} \delta _{K}^{M} + \delta _{I}^{M}\delta _{J}^{N} \delta _{K}^{L} + \delta _{I}^{M}\delta _{J}^{L} \delta _{K}^{N} + \delta _{I}^{N}\delta _{J}^{L} \delta _{K}^{M} + \delta _{I}^{N}\delta _{J}^{M} \delta _{K}^{L}\\
\end{aligned}
 }$$
There are also ${ D^2  }$ and ${ D^4 }$ terms involved, which are functions of ${ \phi ,\overline{\phi } }$; these terms are a scalar multiple of the contraction terms for the ${ D^0 }$ case above;
$${ \delta _{K}^{K}= \delta _{K}^{L}\delta _{L}^{K}= n,\ \delta _{I}^{K}\delta _{K}^{J} = \delta _{I}^{K}\delta _{K}^{L}\delta _{L}^{J}= \delta _{I}^{J} }$$
$${ \begin{aligned}
\pi ^{-n}\int \mathrm{d} ^{2n}\phi  e^{-|\phi |^{2} }\phi _{I}\overline{\phi }^{J} D ^{4} &= \int \mathrm{d} ^{2n}\phi  e^{-|\phi |^{2} } \phi _{I}\overline{\phi }^{J} (|X|^{4} - 2|X|^{2} \phi _{K}\overline{\phi }^{K} + \phi _{K}\phi _{L}\overline{\phi }^{K}\overline{\phi }^{L})\\
&= |X|^{4} \delta _{I}^{J} - 2|X|^{2} (\delta _{I} ^{J}\delta _{K}^{K} + \delta _{I}^{K}\delta _{K}^{J}) + (\delta _{I}^{J}\delta _{K}^{K}\delta _{L}^{L}+\delta _{I}^{J}\delta _{K}^{L}\delta _{L}^{K}+\delta _{I}^{K}\delta _{K}^{L}\delta _{L}^{J}+\delta _{I}^{K}\delta _{K}^{J}\delta _{L}^{L}+\delta _{I}^{L}\delta _{K}^{J}\delta _{L}^{K}+\delta _{I}^{L}\delta _{K}^{K}\delta _{L}^{J})\\
&= |X|^{4} \delta _{I}^{J} - 2|X|^{2} (\delta _{I} ^{J}\cdot n + \delta _{I}^{J}) + (\delta _{I}^{J}\cdot n^2 +\delta _{I}^{J}\cdot n+\delta _{I}^{J}+\delta _{I}^{J}\cdot n+\delta _{I}^{J}+\delta _{I}^{J}\cdot n)\\
&= \delta _{I}^{J}(|X|^{4} - 2(n+1)|X|^{2} + n^2 +3n+2)
\end{aligned}
 }$$
 $${ \begin{aligned}
\pi ^{-n} \int \mathrm{d} ^{2n}\phi e^{-|\phi |^{2} } \phi _{I}\phi _{J}\overline{\phi }^{K}\overline{\phi }^{L} D^2 &= \int \mathrm{d} ^{2n}\phi  e^{-|\phi |^{2} } \phi _{I}\phi _{J}\overline{\phi }^{K}\overline{\phi }^{L} (|X|^{2} -\phi _{M}\overline{\phi }^{M}) \\
 &= |X|^{2} (\delta _{I}^{K}\delta _{J}^{L} + \delta _{I}^{L}\delta _{J}^{K}) - (\delta _{I} ^{K}\delta _{J}^{L}\delta _{M}^{M} +\delta _{I} ^{K}\delta _{J}^{M}\delta _{M}^{L} +\delta _{I} ^{L}\delta _{J}^{M}\delta _{M}^{K} +\delta _{I} ^{L}\delta _{J}^{K}\delta _{M}^{M} +\delta _{I} ^{M}\delta _{J}^{K}\delta _{M}^{L} +\delta _{I} ^{M}\delta _{J}^{L}\delta _{M}^{K} )\\
  &= |X|^{2}( \delta _{I}^{K}\delta _{J}^{L} + \delta _{I}^{L}\delta _{J}^{K}) - (\delta _{I} ^{K}\delta _{J}^{L}\cdot n +\delta _{I} ^{K}\delta _{J}^{L} +\delta _{I} ^{L}\delta _{J}^{K} +\delta _{I} ^{L}\delta _{J}^{K}\cdot n +\delta _{I} ^{L}\delta _{J}^{K} +\delta _{I} ^{K}\delta _{J}^{L} )\\
  &=(\delta _{I}^{K}\delta _{J}^{L}+\delta _{I}^{L}\delta _{J}^{K})(|X|^{2} - (n+2))
 \end{aligned}
  }$$

Then, evaluating the total integral,

$${ \begin{aligned}
&\pi ^{-2(n+1)}e^{-|X|^{2} }|X|^{2n-4}\int \mathrm{d} ^{2n}\phi e^{-|\phi |^{2} } & &\hspace{-.3cm}\biggr(\mathcal{F}_{IJK}\bar{\mathcal{F}}^{JIN}\overline{\phi }^{K}\phi _{N}D ^{4} + \mathcal{F}_{IJK}\bar{\mathcal{F}}^{JMN}\overline{\phi }^{K}\overline{\phi }^{I}\phi _{M}\phi _{N}D^2 \\
& & &+ \mathcal{F}_{IJK}\bar{\mathcal{F}}^{LIN}\overline{\phi }^{J}\overline{\phi }^{K}\phi _{L}\phi _{N}D^2 + \mathcal{F}_{IJK}\bar{\mathcal{F}}^{LMN}\overline{\phi }^{J}\overline{\phi }^{K}\overline{\phi }^{I}\phi _{L}\phi _{M}\phi _{N} \biggr)\\
=& \pi ^{-(n+2)}e^{-|X|^{2} }|X|^{2n-4}[\mathcal{F}_{IJK} \bar{\mathcal{F}}^{JIN} \delta ^{K}_{N} & & \hspace{-.4cm}(|X|^{4}- 2(n+1)|X|^{2} + n^2 +3n+2)  \\
& & &\hspace{-2.5cm}+ \mathcal{F}_{IJK}\bar{\mathcal{F}}^{JMN}(\delta _{M}^{K}\delta _{N}^{I}+\delta _{M}^{I}\delta _{N}^{K}) (|X|^{2} -(n+2))\\
& & &\hspace{-2.5cm}+ \mathcal{F}_{IJK}\bar{\mathcal{F}}^{LIN}(\delta _{L}^{J}\delta _{N}^{K}+\delta _{L}^{K}\delta _{N}^{J}) (|X|^{2} -(n+2))\\
& & &\hspace{-2.5cm}+\mathcal{F}_{IJK}\bar{\mathcal{F}}^{LMN}(\delta _{I}^{L}\delta _{J}^{M} \delta _{K}^{N} + \delta _{I}^{L}\delta _{J}^{N} \delta _{K}^{M} + \delta _{I}^{M}\delta _{J}^{N} \delta _{K}^{L} + \delta _{I}^{M}\delta _{J}^{L} \delta _{K}^{N} + \delta _{I}^{N}\delta _{J}^{L} \delta _{K}^{M} + \delta _{I}^{N}\delta _{J}^{M} \delta _{K}^{L})]\\
\end{aligned}
 }$$

\paragraph{Step 3} Third, contract indices on ${ \mathcal{F}\bar{\mathcal{F}} }$ into scalar ${ U(n) }$ invariants. For ${ d=1 }$, there is only one invariant ${ \mathcal{F}_{IJK}\bar{\mathcal{F}}^{IJK} = |\mathcal{F}|^{2} }$. Then, integrate out ${ |X| }$.

$${ \begin{aligned}
&\pi ^{-(n+2)}\int \mathrm{d} ^2 X e^{-|X|^{2} }|X|^{2n-4}|\mathcal{F}|^{2}[ (|X|^{4}-2(n+1)|X|^{2} +n^2 +3n+2)+4 (|X|^{2} -(n+2))+6]\\
=&\pi ^{-(n+2)}|\mathcal{F}|^{2}\int \mathrm{d} ^2 X e^{-|X|^{2} }|X|^{2n-4}[ |X|^{4}+(-2n+2)|X|^{2} +n^2 -n ]\\
=&\pi ^{-(n+2)}|\mathcal{F}|^{2} (\ev{|X|^{2n}} + (-2n+2)\ev{|X|^{2n-2}} + (n^2 -n)\ev{|X|^{2n-4}}) \\
\end{aligned}
 }$$

Here, it is noted that for ${ n=1 }$, the final term ${ \ev{|X|^{2n-4}} =\pi  (-1)! }$ is a singular term, but its total coefficient is 0. This is reflected in the fact that for ${ n=1 }$, there is no singular term in the integrand for the initial Douglas-Denef equation. For higher orders, the proposition below (which is proved in a later section) is used;

\paragraph{Proposition 2. } Expanding the integrand of ${ \rho  }$ in powers of ${ |X| }$, the only singular term is proportional to ${ |X|^{-2}}$, and its coefficient is 0.

## ${ \boldsymbol{d=2} }$

\paragraph{Step 1} First, the expression ${ \exp \tr \log (\mathbb{I} - |X|^{-2} N) }$ is evaluated up to 2nd order in ${ N }$. 
$${ \log (\mathbb{I}- |X|^{-2}N) =- |X|^{-2}N - \dfrac{1}{2}|X|^{-4} N^2 -\cdots  }$$
$${ \begin{aligned}
\exp \tr \log (\mathbb{I}-|X|^{-2}N) &= \exp \left(-|X|^{-2}\tr N - \dfrac{1}{2}|X|^{-4}\tr N^2 -\cdots \right)\\
&= \cdots  - \dfrac{1}{2}|X|^{-4}\tr N^2 + \dfrac{1}{2}|X|^{-4}(\tr N)^{2} \cdots \\
&= \cdots - \dfrac{1}{2}|X|^{-4}N_I ^{J}N_J^I + \dfrac{1}{2}|X|^{-4} N_I^I N_J^J\cdots 
\end{aligned}
 }$$
Both second order contributions are of the form ${ N _{I_1} ^{P_1} N _{I_2} ^{P_2} }$ with "external indices" contracted via
$${ \begin{aligned}
&N _{I}^{J}N _{J}^{I} =N _{I_1}^{P_1} N _{I_2}^{P_2}\delta _{I_1}^{P_2}\delta _{I_2}^{P_1},\\
&N_I^IN_J^J = N _{I_1}^{P_1}N _{I_2}^{P_2} \delta _{I_1}^{P_1}\delta _{I_2}^{P_2}.
\end{aligned}
 }$$

The two contractions above are all the possible contractions of ${ N _{I_1}^{P_1} N _{I_2}^{P_2} }$, and the final contribution is a sum over them with a total prefactor of ${ \dfrac{\pm 1}{2!} }$. While it isn't clear yet, the proposition below (which is proved in a later section) is used for higher orders;

\paragraph{Proposition 3. } In the series expansion of ${ \exp \tr \log (\mathbb{I}-|X|^{-2}N) }$, the ${ d }$-th order term equals the sum over all possible contractions of external indices, with a prefactor ${ (-1)^{\text{disj}(\sigma )} / d! }$, where ${ \text{disj} (\sigma ) }$ refers to the number of disjoint cycles in the permutation ${ \sigma  : I_i \to  }$ (index ${ I_i }$ is contracted with).

\newpage

For a general prefactor ${ A }$, the contribution of ${ N _{I_1}^{P_1} N _{I_2}^{P_2} }$ is

$${ \begin{aligned}
&e^{-|X|^{2}-|\phi |^{2}} |X|^{2n-2}D^4 \cdot |X|^{-4} \left(A N _{I_1}^{P_1}N _{I_2}^{P_2} \right) \\
=&e^{-|X|^{2}-|\phi |^{2}} |X|^{2n-6}A\biggr(\mathcal{F}_{I_1J_1K_1}\bar{\mathcal{F}}^{J_1P_1N_1} \mathcal{F}_{I_2J_2K_2}\bar{\mathcal{F}}^{J_2P_2N_2}\overline{\phi }^{K_1}\phi _{N_1}\overline{\phi }^{K_2}\phi _{N_2} D ^{4}\\
& + 2\mathcal{F}_{I_1J_1K_1}\bar{\mathcal{F}}^{J_1P_1N_1}\mathcal{F}_{I_2J_2K_2}\bar{\mathcal{F}}^{J_2M_2N_2}\overline{\phi }^{K_1}\phi _{N_1}\overline{\phi }^{K_2}\overline{\phi }^{I}\phi _{M_2}\phi _{N_2} D ^{2} + 2\mathcal{F}_{I_1J_1K_1}\bar{\mathcal{F}}^{J_1P_1N_1}\mathcal{F}_{I_2J_2K_2}\bar{\mathcal{F}}^{L_2P_2N_2}\overline{\phi }^{K_1}\phi _{N_1}\overline{\phi }^{J_2}\overline{\phi }^{K_2}\phi _{L_2}\phi _{N_2} D ^{2}\\
& +\mathcal{F}_{I_1J_1K_1}\bar{\mathcal{F}}^{J_1M_1N_1}\mathcal{F}_{I_2J_2K_2}\bar{\mathcal{F}}^{J_2M_2N_2}\overline{\phi }^{K_1}\overline{\phi }^{P_1}\phi _{M_1}\phi _{N_1}\overline{\phi }^{K_2}\overline{\phi }^{P_2}\phi _{M_2}\phi _{N_2}\\
&+2\mathcal{F}_{I_1J_1K_1}\bar{\mathcal{F}}^{J_1M_1N_1}\mathcal{F}_{I_2J_2K_2}\bar{\mathcal{F}}^{L_2P_2N_2}\overline{\phi }^{K_1}\overline{\phi }^{P_1}\phi _{M_1}\phi _{N_1}\overline{\phi }^{J_2}\overline{\phi }^{K_2}\phi _{L_2}\phi _{N_2}+ \cdots \biggr)  \\
\end{aligned}
}$$

Each term in the expanded expression precisely corresponds to a combination with repetition chosen from the 4 terms in ${ N_I^P }$. Each term also has two lower indices ${ I_1,I_2 }$ and two upper indices ${ P_1,P_2 }$, which summing over all contractions (with appropriate prefactors) will give the total contribution to the integrand. Finally, each term has a number of "internal indices" attached to the variables ${ \phi ,\overline{\phi } }$, which evaluating the integral in Step 2 is equivalent to summing over all contractions. 

\paragraph{Step 2} The same method as in ${ d=1 }$ is used to integrate out the variables ${ \phi ,\overline{\phi } }$ via contractions. However, unlike ${ d=1 }$, there are terms that are singular in ${ D }$, corresponding to the terms in ${ \exp \tr \log (\mathbb{I}-|X|^{-2} N) }$ proportional to ${ D ^{-6} }$ or worse. Attempting to regulate these singularities is an ill-defined problem. Instead, the following proposition (which is proved in a later section) is used;

\paragraph{Proposition 4. } Expanding the integrand of ${ \rho  }$ as a formal series in ${ D^2 = |X|^{2} -|\phi |^{2}  }$, the singular coefficients vanish.

$${ N_I^P = N _{(I)}{}_{I}^{P} + N _{(II)}{}_{I}^{P} D ^{-2}+N _{(I I I)}{}_{I}^{P}D ^{-2}+N _{(IV)}{}_{I}^{P} D ^{-4}}$$
Writing as above, there are only 7 remaining terms 
$${ \underbrace{N _{(I)}N _{(IV)}, N _{(I I)}N _{(I I)}, N _{(I I)}N _{(I I I)}, N _{(I I I)}N _{(I I I)}}_{D ^{-4}}, \underbrace{N _{(I )} N _{(I I)}, N _{ (I)} N _{(I I I)}}_{D ^{-2}}, \underbrace{N _{(I)} N _{(I)}}_{D^0} }$$
and this is true for any ${ d \geq 2 }$ (higher order terms can only include more ${ N _{(I)} }$).  
For demonstration purposes, the first term in the expanded equation above for the contraction ${ N_I^JN_J^I }$ is computed;
$${ \begin{aligned}
&\pi ^{-2(n+1)}e^{-|X|^{2} }|X|^{2n-6}\cdot -\dfrac{1}{2} \int \mathrm{d} ^{2n} \phi e^{-|\phi |^{2} } \mathcal{F}_{I_1J_1K_1}\bar{\mathcal{F}}^{J_1P_1N_1}\mathcal{F}_{P_1J_2K_2}\bar{\mathcal{F}}^{J_2I_1N_2}\overline{\phi }^{K_1}\phi _{N_1}\overline{\phi }^{K_2}\phi _{N_2}D^4\\
=&\pi ^{-(n+2)}e^{-|X|^{2} }|X|^{2n-6}\cdot -\dfrac{1}{2} \mathcal{F}_{I_1J_1K_1}\bar{\mathcal{F}}^{J_1P_1N_1}\mathcal{F}_{P_1J_2K_2}\bar{\mathcal{F}}^{J_2I_1N_2}(|X|^{4}+(-2n-4)|X|^{2} +(n^2 +5n+8))(\delta ^{K_1}_{N_1}\delta ^{K_2}_{N_2}+\delta ^{K_1}_{N_2}\delta ^{K_2}_{N_1}),
\end{aligned}
 }$$

\paragraph{Step 3} For demonstration purposes, contract indices on ${ \mathcal{F}\bar{\mathcal{F}} }$ for the term above. For ${ d=2 }$, there are two possible invariants, ${ |\mathcal{F}|^{4} = \mathcal{F}_{(3,0)} = \mathcal{F}_{IJK}\mathcal{F}_{LMN}\bar{\mathcal{F}}^{IJK}\bar{\mathcal{F}}^{LMN} }$ and ${ \mathcal{F}_{(2,1)} = \mathcal{F}_{IJK}\mathcal{F}_{LMN}\bar{\mathcal{F}}^{LJK}\bar{\mathcal{F}}^{IMN} }$.
$${ \begin{aligned}
&\pi ^{-(n+2)}\int \mathrm{d} ^2 X e^{-|X|^{2} }|X|^{2n-6}\cdot -\dfrac{1}{2} \mathcal{F}_{I_1J_1K_1}\bar{\mathcal{F}}^{J_1P_1N_1}\mathcal{F}_{P_1J_2K_2}\bar{\mathcal{F}}^{J_2I_1N_2}(|X|^{4}+(-2n-4)|X|^{2} +(n^2 +5n+8))(\delta ^{K_1}_{N_1}\delta ^{K_2}_{N_2}+\delta ^{K_1}_{N_2}\delta ^{K_2}_{N_1})\\
&=\pi ^{-(n+2)} \cdot -\dfrac{1}{2}(\mathcal{F}_{(2,1)}+\mathcal{F}_{(2,1)})\int \mathrm{d} ^2 X e^{-|X|^{2} }(|X|^{2n-2} + (-2n-4)|X|^{2n-4} + (n^2 +5n+8)|X|^{2n-6})\\
&= -\pi ^{-(n+2)}\mathcal{F}_{(2,1)} (\ev{|X|^{2n-2}} + (-2n-4)\ev{|X|^{2n-4}}+ (n^2 +5n+8)\ev{|X|^{2n-6}})
\end{aligned}
 }$$
Due to Proposition 1, there is no need to consider ${ n=1 }$ for ${ d=2 }$. Due to Proposition 2, the possible singular term at ${ n=2 }$, ${ \ev{|X|^{2n-6}} }$ has vanishing coefficient when adding together all coefficients.

\newpage

