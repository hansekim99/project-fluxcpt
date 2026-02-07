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

## Proofs of Propositions
\paragraph{Proof of Proposition 4.}
Assume that when formally expanding the integrand as a series in ${ D }$, there are singular terms with nonzero coefficients. Then for any fixed value of ${ |X|^{2}  }$ and any path ${ |\phi |^{2} \to |X|^{2}  }$, the limit of the integrand along this path would be divergent. Fix such a path, ${ \phi = \overline{\phi }= a (|X|, 0 , \cdots ) }$ such that the limit is realised by ${ a \to 1 }$. Then,
$${ \begin{aligned}
\begin{pmatrix}
\delta ^{I}_{J} \overline{X} - \overline{\phi }^{I}\phi _{J}/X & \mathcal{F}_{IJK} \overline{\phi }^{K}\\ \bar{\mathcal{F}}^{IJK}\phi _{K} & \delta _{I}^{J} X - \phi _{I}\overline{\phi }^{J}/\overline{X}
\end{pmatrix} &= \begin{pmatrix}
\overline{X} - a^2 |X|^{2} /X & 0 & \cdots & 0 & \\
0 & \overline{X} & \cdots & 0 & & \mathcal{F}_{IJ1} a|X| \\
\vdots & \vdots & \ddots & 0 & \\
0 & 0 & 0 & \overline{X} & \\
& & & & X - a^2 |X|^{2} /\overline{X} & 0 & \cdots & 0 & \\
& \bar{\mathcal{F}}^{IJ1} a |X| & & & 0 & X & \cdots & 0 & \\
& & & & \vdots & \vdots & \ddots & 0 & \\
& & & & 0 & 0 & 0 & X & 
\end{pmatrix}\\
& \to \begin{pmatrix}
\ 0 & 0 & \cdots & 0 & \\
\ 0 & \overline{X} & \cdots & 0 & & \mathcal{F}_{IJ1} |X| \\
\ \vdots & \vdots & \ddots & 0 & \\
\ 0 & 0 & 0 & \overline{X} & \\
& & & & 0 & 0 & \cdots & 0 & \\
& \bar{\mathcal{F}}^{IJ1}  |X| & & & 0 & X & \cdots & 0 & \\
& & & & \vdots & \vdots & \ddots & 0 & \\
& & & & 0 & 0 & 0 & X & 
\end{pmatrix}
\end{aligned}
 }$$
However, the determinant of the matrix above cannot be divergent for any value of ${ |X| }$, as it is a finite degree polynomial. Therefore, there cannot be any singular terms in the formal series of the integrand.

\paragraph{Proof of Proposition 2. }
$${ \det \begin{pmatrix}
\delta ^{I}_{J} \overline{X} - \overline{\phi }^{I}\phi _{J}/X & \mathcal{F}_{IJK} \overline{\phi }^{K}\\ \bar{\mathcal{F}}^{IJK}\phi _{K} & \delta _{I}^{J} X - \phi _{I}\overline{\phi }^{J}/\overline{X}
\end{pmatrix}
 }$$
The matrix above is manifestly invariant under ${ U(n) }$ rotations 
$${ \phi _{I} \to U _{I}^{J}\phi _{J},\ \overline{\phi }^{I} \to {U^{\dagger}}^{I}_{J}\overline{\phi }^{J},\ \mathcal{F}_{IJK} \to U _{I}^{I^{\prime}}U _{J}^{J^{\prime}}U _{K}^{K^{\prime}} \mathcal{F}_{I^{\prime}J^{\prime}K^{\prime}}, }$$
and so is the exponent ${ e^{-|\phi |^{2} } }$. Consider the rotation such that ${ \phi _{I} \to (|\phi |,0,\cdots ), \overline{\phi }_{I} \to (|\phi |,0,\cdots ) }$ and define ${ \mathcal{F}_{IJK} \to \mathcal{F}^{\prime}_{IJK} }$. 
$${ \begin{aligned}
\begin{pmatrix}
\delta ^{I}_{J} \overline{X} - \overline{\phi }^{I}\phi _{J}/X & \mathcal{F}_{IJK} \overline{\phi }^{K}\\ \bar{\mathcal{F}}^{IJK}\phi _{K} & \delta _{I}^{J} X - \phi _{I}\overline{\phi }^{J}/\overline{X}
\end{pmatrix} &= \begin{pmatrix}
\overline{X} - |\phi |^{2} /X & 0 & \cdots & 0 & \\
0 & \overline{X} & \cdots & 0 & & \mathcal{F}^{\prime}_{IJ1} |\phi | \\
\vdots & \vdots & \ddots & 0 & \\
0 & 0 & 0 & \overline{X} & \\
& & & & X - |\phi |^{2}  /\overline{X} & 0 & \cdots & 0 & \\
& \bar{\mathcal{F}}^{\prime IJ1} |\phi | & & & 0 & X & \cdots & 0 & \\
& & & & \vdots & \vdots & \ddots & 0 & \\
& & & & 0 & 0 & 0 & X & 
\end{pmatrix}\\
\end{aligned}
 }$$
The determinant of the matrix above then can at most have a singularity of ${ 1 / X \overline{X} = |X|^{-2} }$, which is cancelled out by the ${ |X|^{2}  }$ term in the integrand. Therefore, there is no singular term in ${ |X| }$ in the integrand.

\paragraph{Proof of Proposition 3. }
To motivate the general construction, the example of ${ d=3 }$ is considered;
$${ \begin{aligned}
\exp \tr \log (\mathbb{I}-|X|^{-2}N) &= \exp \left(-|X|^{-2}\tr N -\dfrac{1}{2}|X|^{-4}\tr N^2 - \dfrac{1}{3}|X|^{-6}\tr N^3 -\cdots \right)\\
& = \cdots  - \dfrac{1}{3}|X|^{-6}\tr N^3 + \dfrac{1}{2!} \cdot  \dfrac{1}{2} \cdot 2 |X|^{-6} \tr N \tr N^2 - \dfrac{1}{3!} |X|^{-6} (\tr N)^{3} + \cdots\\  
& = \cdots  - \dfrac{1}{3}|X|^{-6}\tr N^3 + \dfrac{1}{2} |X|^{-6} \tr N \tr N^2 - \dfrac{1}{6} |X|^{-6} (\tr N)^{3} + \cdots\\
\end{aligned}
 }$$
$${ \begin{aligned}
N _{I_1}^{P_1} N _{I_2}^{P_2} N _{I_3}^{P_3} \delta _{I_1}^{P_1}\delta _{I_2}^{P_2}\delta _{I_3}^{P_3} &= N _{I}^{I} N _{J}^{J} N _{K}^{K} = (\tr N)^{3} ,& &\text{disj}(\sigma ) = 3\\
N _{I_1}^{P_1} N _{I_2}^{P_2} N _{I_3}^{P_3} \delta _{I_1}^{P_1}\delta _{I_2}^{P_3}\delta _{I_3}^{P_2} &= N _{I}^{I} N_J^K N_K^J = \tr N \tr N^2,& &\text{disj}(\sigma ) = 2\\
N _{I_1}^{P_1} N _{I_2}^{P_2} N _{I_3}^{P_3} \delta _{I_1}^{P_2}\delta _{I_2}^{P_3}\delta _{I_3}^{P_1} &= N _{I}^{J} N _{J}^{K} N _{K}^{I} = \tr N^3 ,& &\text{disj}(\sigma ) = 1 \\
N _{I_1}^{P_1} N _{I_2}^{P_2} N _{I_3}^{P_3} \delta _{I_1}^{P_2}\delta _{I_2}^{P_1}\delta _{I_3}^{P_3} &= N _{I}^{J} N _{J}^{I} N _{K}^{K} = \tr N^2 \tr N ,& &\text{disj}(\sigma ) = 2 \\
N _{I_1}^{P_1} N _{I_2}^{P_2} N _{I_3}^{P_3} \delta _{I_1}^{P_3}\delta _{I_2}^{P_1}\delta _{I_3}^{P_2} &= N _{I}^{K} N _{J}^{I} N _{K}^{J} = \tr N^{3} ,& &\text{disj}(\sigma ) = 1\\
N _{I_1}^{P_1} N _{I_2}^{P_2} N _{I_3}^{P_3} \delta _{I_1}^{P_3}\delta _{I_2}^{P_2}\delta _{I_3}^{P_1} &= N _{I}^{K} N _{J}^{J} N _{K}^{I} = \tr N \tr N^2 ,& &\text{disj}(\sigma ) = 2 \\
\end{aligned}
 }$$
The coefficient is accounted for combinatorically by the multiplicity of permutations generating the cycles corresponding to each trace term. Therefore, the final expression is a sum over all possible contractions with a prefactor ${ (-1)^{\text{disj}(\sigma )}/d! }$. The generalised proof is below;
$${ \begin{aligned}
&|X| ^{2n-2}D^4\exp \tr\log (\mathbb{I}-|X|^{-2}N) = |X|^{2n-2} D^4\left(1 + \left(\sum_{k}^{} -\dfrac{1}{k}|X|^{-2k}\tr N^k \right)+ \dfrac{1}{2!}\left(\sum_{k}^{}-\dfrac{1}{k} |X|^{-2k}\tr N^k \right)^{2} + \cdots \right)\\
&= _{(1)}D^4\sum_{d}^{} \sum_{\substack{m = 1 , \cdots ,d\\ k_1 + \cdots + k_m = d}}^{}   \dfrac{1}{m!} \cdot \dfrac{|X|^{2n-2-2k_1 \cdots -2k_m}(-1)^{m}}{k_1 k_2 \cdots k_m} \dfrac{m!}{\prod _{k}(\text{number of identical }k_i)!} \cdot \tr N ^{k_1} \tr N ^{k_2} \cdots \tr N ^{k_m}\\
&= _{(2)} D^4\sum_{d}^{} \sum_{\substack{m = 1 , \cdots ,d\\ k_1 + \cdots + k_m = d}} \dfrac{|X|^{2n-2d-2}(-1)^{m}}{d!} N(k_1,\cdots ,k_m) \cdot \tr N ^{k_1} \tr N ^{k_2} \cdots \tr N ^{k_m} \\
&= _{(3)} D^4\sum_{d}^{} \sum_{\substack{m = 1 , \cdots ,d\\ k_1 + \cdots + k_m = d}} \dfrac{|X|^{2n-2d-2}(-1)^{m}}{d!} \cdot \sum_{\substack{\sigma \in S_d\\ \text{disj. cycles }\\ k_1,\cdots ,k_d}} N _{I_1} ^{I _{\sigma (1)}} \cdots N _{I_d} ^{I _{\sigma (d)}}\\
&= _{(4)} D^4\sum_{d}^{} \sum_{\sigma \in S_d} \dfrac{|X| ^{2n-2d-2}(-1)^{\text{disj}(\sigma )}}{d!} N _{I_1} ^{I _{\sigma (1)}} \cdots N _{I_d} ^{I _{\sigma (d)}}   
\end{aligned}
 }$$

- (1) : Obtained from directly expanding the order ${ d }$ coefficient of the ${ m }$-th term in the ${ \exp \tr \log  }$ equation.
- (2) : ${ N(k_1,\cdots ,k_m) }$ is the number of permutations in ${ S_d }$ that have ${ m }$ disjoint cycles of length ${ k_1, \cdots ,k_m }$. 
$${ N(k_1,\cdots ,k_m) = \dfrac{d!}{k_1! k_2! \cdots k_m!} \dfrac{1}{\prod _{k}(\text{number of identical }k_i)!} \cdot (k_1-1)!(k_2-1)!\cdots (k_m-1)! }$$
- (3) : If ${ \sigma \in S_d }$ is a permutation with disjoint cycles of length ${ k_1,\cdots ,k_m }$, the contracted sum ${ N _{I_1} ^{I _{\sigma (1)}}\cdots N _{I_d} ^{I _{\sigma (d)}} }$ equals ${ \tr N ^{k_1}\cdots \tr N ^{k_m} }$.
- (4) : ${ \text{disj}(\sigma ) }$ refers to the number of disjoint cycles in ${ \sigma  }$.

\newpage

## Results

Based on the notation introduced in the previous section, we classify all invariants up to order 5.

$${ \begin{aligned}
|\mathcal{F}|^{4} =\mathcal{F}^{(2)}_{1}&=((3, 0), (0, 3)),\ &\mathcal{F}_{(2,1)} = \mathcal{F}^{(2)}_{2}&=((2, 1), (1, 2))\\
|\mathcal{F}|^{6} = \mathcal{F}^{(3)}_{1}&=((3, 0, 0), (0, 3, 0), (0, 0, 3)),\ &|\mathcal{F}|^{2} \mathcal{F}_{(2,1)} = \mathcal{F}^{(3)}_{2}&=((3, 0, 0), (0, 2, 1), (0, 1, 2))\\
\mathcal{F}^{(3)}_{3}&=((2, 1, 0), (0, 2, 1), (1, 0, 2)),&\ \mathcal{F}^{(3)}_{4}&=((2, 1, 0), (1, 1, 1), (0, 1, 2))\\
\mathcal{F}^{(3)}_{5}&=((1, 1, 1), (1, 1, 1), (1, 1, 1))\\
 \end{aligned}}$$

$${ \begin{aligned}
|\mathcal{F}|^{8}=\mathcal{F}^{(4)}_{1}&=((3, 0, 0, 0), (0, 3, 0, 0), (0, 0, 3, 0), (0, 0, 
0, 3))\\
|\mathcal{F}|^{4}\mathcal{F}_{(2,1)}=\mathcal{F}^{(4)}_{2}&=((3, 0, 0, 0), (0, 3, 0, 0), (0, 0, 2, 1), (0, 0, 
1, 2))\\
|\mathcal{F}|^{2} \mathcal{F}^{(3)}_{3}=\mathcal{F}^{(4)}_{3}&=((3, 0, 0, 0), (0, 2, 1, 0), (0, 0, 2, 1), (0, 1, 
0, 2))\\
\mathcal{F}_{(2,1)}^{2}= \mathcal{F}^{(4)}_{4}&=((2, 1, 0, 0), (1, 2, 0, 0), (0, 0, 2, 1), (0, 0, 
1, 2))\\
\mathcal{F}^{(4)}_{5}&=((2, 1, 0, 0), (0, 2, 1, 0), (0, 0, 2, 1), (1, 0, 
0, 2))\\
|\mathcal{F}|^{2} \mathcal{F}^{(3)}_{4}=\mathcal{F}^{(4)}_{6}&=((3, 0, 0, 0), (0, 2, 1, 0), (0, 1, 1, 1), (0, 0, 
1, 2))\\
\mathcal{F}^{(4)}_{7}&=((2, 1, 0, 0), (0, 2, 1, 0), (1, 0, 1, 1), (0, 0, 
1, 2))\\
\mathcal{F}^{(4)}_{8}&=((2, 0, 1, 0), (0, 2, 0, 1), (1, 0, 1, 1), (0, 1, 
1, 1))\\
|\mathcal{F}|^{2} \mathcal{F}^{(3)}_{5}=\mathcal{F}^{(4)}_{9}&=((3, 0, 0, 0), (0, 1, 1, 1), (0, 1, 1, 1), (0, 1, 
1, 1))\\
\mathcal{F}^{(4)}_{10}&=((2, 1, 0, 0), (0, 1, 2, 0), (1, 0, 1, 1), (0, 1,
 0, 2))\\
\mathcal{F}^{(4)}_{11}&=((2, 1, 0, 0), (0, 1, 1, 1), (1, 0, 1, 1), (0, 1,
 1, 1))\\
\end{aligned}
 }$$
$${ \begin{aligned}
\mathcal{F}^{(5)}_{1}&=((3, 0, 0, 0, 0), (0, 3, 0, 0, 0), (0, 0, 3, 0, 0)
, (0, 0, 0, 3, 0), (0, 0, 0, 0, 3))\\
\mathcal{F}^{(5)}_{2}&=((3, 0, 0, 0, 0), (0, 3, 0, 0, 0), (0, 0, 3, 0, 0)
, (0, 0, 0, 2, 1), (0, 0, 0, 1, 2))\\
\mathcal{F}^{(5)}_{3}&=((3, 0, 0, 0, 0), (0, 3, 0, 0, 0), (0, 0, 2, 1, 0)
, (0, 0, 0, 2, 1), (0, 0, 1, 0, 2))\\
\mathcal{F}^{(5)}_{4}&=((3, 0, 0, 0, 0), (0, 2, 1, 0, 0), (0, 1, 2, 0, 0)
, (0, 0, 0, 2, 1), (0, 0, 0, 1, 2))\\
\mathcal{F}^{(5)}_{5}&=((3, 0, 0, 0, 0), (0, 2, 1, 0, 0), (0, 0, 2, 1, 0)
, (0, 0, 0, 2, 1), (0, 1, 0, 0, 2))\\
\mathcal{F}^{(5)}_{6}&=((2, 1, 0, 0, 0), (1, 2, 0, 0, 0), (0, 0, 2, 1, 0)
, (0, 0, 0, 2, 1), (0, 0, 1, 0, 2))\\
\mathcal{F}^{(5)}_{7}&=((2, 1, 0, 0, 0), (0, 2, 1, 0, 0), (0, 0, 2, 1, 0)
, (0, 0, 0, 2, 1), (1, 0, 0, 0, 2))\\
\mathcal{F}^{(5)}_{8}&=((3, 0, 0, 0, 0), (0, 3, 0, 0, 0), (0, 0, 2, 1, 0)
, (0, 0, 1, 1, 1), (0, 0, 0, 1, 2))\\
\mathcal{F}^{(5)}_{9}&=((3, 0, 0, 0, 0), (0, 2, 1, 0, 0), (0, 0, 2, 1, 0)
, (0, 1, 0, 1, 1), (0, 0, 0, 1, 2))\\
\mathcal{F}^{(5)}_{10}&=((3, 0, 0, 0, 0), (0, 2, 0, 1, 0), (0, 0, 2, 0, 1
), (0, 1, 0, 1, 1), (0, 0, 1, 1, 1))\\
\mathcal{F}^{(5)}_{11}&=((2, 1, 0, 0, 0), (1, 2, 0, 0, 0), (0, 0, 2, 1, 0
), (0, 0, 1, 1, 1), (0, 0, 0, 1, 2))\\
\mathcal{F}^{(5)}_{12}&=((2, 1, 0, 0, 0), (0, 2, 1, 0, 0), (0, 0, 2, 1, 0
), (1, 0, 0, 1, 1), (0, 0, 0, 1, 2))\\
\mathcal{F}^{(5)}_{13}&=((2, 1, 0, 0, 0), (0, 2, 0, 1, 0), (0, 0, 2, 0, 1
), (1, 0, 0, 1, 1), (0, 0, 1, 1, 1))\\
\mathcal{F}^{(5)}_{14}&=((3, 0, 0, 0, 0), (0, 3, 0, 0, 0), (0, 0, 1, 1, 1
), (0, 0, 1, 1, 1), (0, 0, 1, 1, 1))\\
\end{aligned}
 }$$
 $${ \begin{aligned}
 \mathcal{F}^{(5)}_{15}&=((3, 0, 0, 0, 0), (0, 2, 1, 0, 0), (0, 0, 1, 2, 0
), (0, 1, 0, 1, 1), (0, 0, 1, 0, 2))\\
\mathcal{F}^{(5)}_{16}&=((3, 0, 0, 0, 0), (0, 2, 1, 0, 0), (0, 0, 1, 1, 1
), (0, 1, 0, 1, 1), (0, 0, 1, 1, 1))\\
\mathcal{F}^{(5)}_{17}&=((2, 1, 0, 0, 0), (1, 2, 0, 0, 0), (0, 0, 1, 1, 1
), (0, 0, 1, 1, 1), (0, 0, 1, 1, 1))\\
\mathcal{F}^{(5)}_{18}&=((2, 1, 0, 0, 0), (0, 2, 1, 0, 0), (1, 0, 1, 1, 0
), (0, 0, 0, 2, 1), (0, 0, 1, 0, 2))\\
\mathcal{F}^{(5)}_{19}&=((2, 1, 0, 0, 0), (0, 2, 1, 0, 0), (0, 0, 1, 2, 0
), (1, 0, 0, 1, 1), (0, 0, 1, 0, 2))\\
\mathcal{F}^{(5)}_{20}&=((2, 1, 0, 0, 0), (0, 2, 1, 0, 0), (0, 0, 1, 1, 1
), (1, 0, 0, 1, 1), (0, 0, 1, 1, 1))\\
\mathcal{F}^{(5)}_{21}&=((2, 0, 1, 0, 0), (0, 2, 0, 1, 0), (1, 0, 1, 1, 0
), (0, 1, 0, 1, 1), (0, 0, 1, 0, 2))\\
\mathcal{F}^{(5)}_{22}&=((2, 0, 1, 0, 0), (0, 2, 0, 1, 0), (1, 0, 1, 1, 0
), (0, 0, 0, 1, 2), (0, 1, 1, 0, 1))\\
\mathcal{F}^{(5)}_{23}&=((2, 0, 1, 0, 0), (0, 2, 0, 1, 0), (0, 0, 1, 1, 1
), (1, 0, 0, 1, 1), (0, 1, 1, 0, 1))\\
\mathcal{F}^{(5)}_{24}&=((2, 1, 0, 0, 0), (1, 1, 1, 0, 0), (0, 0, 1, 1, 1
), (0, 0, 1, 1, 1), (0, 1, 0, 1, 1))\\
\mathcal{F}^{(5)}_{25}&=((2, 1, 0, 0, 0), (0, 1, 2, 0, 0), (0, 0, 1, 1, 1
), (1, 0, 0, 1, 1), (0, 1, 0, 1, 1))\\
\mathcal{F}^{(5)}_{26}&=((2, 1, 0, 0, 0), (0, 1, 1, 0, 1), (1, 0, 1, 1, 0
), (0, 0, 0, 2, 1), (0, 1, 1, 0, 1))\\
\mathcal{F}^{(5)}_{27}&=((1, 1, 1, 0, 0), (1, 1, 0, 1, 0), (0, 0, 1, 1, 1
), (1, 0, 0, 1, 1), (0, 1, 1, 0, 1))\\
\mathcal{F}^{(5)}_{28}&=((1, 1, 1, 0, 0), (1, 1, 1, 0, 0), (0, 0, 1, 1, 1
), (0, 1, 0, 1, 1), (1, 0, 0, 1, 1))\\
 \end{aligned}
  }$$

The index density can be expressed as a linear combination of the invariants defined above;
$${ \begin{aligned}
n=1: \rho _{\text{ind}} = \pi ^{-2} &(\, 2-\mathcal{F}^{(1)}_{1}) \\
n=2: \rho _{\text{ind}} = \pi ^{-3}&(
4-\mathcal{F}^{(1)}_{1}+\tfrac{1}{2}\mathcal{F}^{(2)}_{1}-\tfrac{1}{2}\mathcal{F}^{(2)}_{2}
)\\
n=3:\rho_ {\text{ind}} = \pi ^{-4}&(
12-4\mathcal{F}^{(1)}_{1}+\mathcal{F}^{(2)}_{1}-\mathcal{F}^{(2)}_{2}-\tfrac{1}{6}\mathcal{F}^{(3)}_{1}+\tfrac{1}{2}\mathcal{F}^{(3)}_{2}-\mathcal{F}^{(3)}_{3}+\mathcal{F}^{(3)}_{4}-\tfrac{1}{3}\mathcal{F}^{(3)}_{5}
)\\
n=4:\rho _{\text{ind}} = \pi ^{-5}&(
48-12\mathcal{F}^{(1)}_{1}+2\mathcal{F}^{(2)}_{1}-2\mathcal{F}^{(2)}_{2}-\tfrac{1}{3}\mathcal{F}^{(3)}_{1}+\mathcal{F}^{(3)}_{2}-2\mathcal{F}^{(3)}_{3}+2\mathcal{F}^{(3)}_{4}-\tfrac{2}{3}\mathcal{F}^{(3)}_{5}\\
&+\tfrac{1}{24}\mathcal{F}^{(4)}_{1}-\tfrac{1}{4}\mathcal{F}^{(4)}_{2}+\mathcal{F}^{(4)}_{3}+\tfrac{1}{8}\mathcal{F}^{(4)}_{4}-\tfrac{1}{4}\mathcal{F}^{(4)}_{5}-\mathcal{F}^{(4)}_{6}-2\mathcal{F}^{(4)}_{7}+2\mathcal{F}^{(4)}_{8}+\tfrac{1}{3}\mathcal{F}^{(4)}_{9}+\mathcal{F}^{(4)}_{10}-\mathcal{F}^{(4)}_{11}\\
n=5:\rho _{\text{ind}}= \pi ^{-6}&(240-48\mathcal{F}^{(1)}_{1}+6\mathcal{F}^{(2)}_{1}-6\mathcal{F}^{(2)}_{2}-\tfrac{2}{3}\mathcal{F}^{(3)}_{1}+2\mathcal{F}^{(3)}_{2}-4\mathcal{F}^{(3)}_{3}+4\mathcal{F}^{(3)}_{4}-\tfrac{4}{3}\mathcal{F}^{(3)}_{5}\\
&+\tfrac{1}{12}\mathcal{F}^{(4)}_{1}-\tfrac{1}{2}\mathcal{F}^{(4)}_{2}+2\mathcal{F}^{(4)}_{3}+\tfrac{1}{4}\mathcal{F}^{(4)}_{4}-\tfrac{1}{2}\mathcal{F}^{(4)}_{5}-2\mathcal{F}^{(4)}_{6}-4\mathcal{F}^{(4)}_{7}+4\mathcal{F}^{(4)}_{8}+\tfrac{2}{3}\mathcal{F}^{(4)}_{9}+2\mathcal{F}^{(4)}_{10}-2\mathcal{F}^{(4)}_{11}\\
&-\tfrac{1}{120}\mathcal{F}^{(5)}_{1}+\tfrac{1}{12}\mathcal{F}^{(5)}_{2}-\tfrac{1}{2}\mathcal{F}^{(5)}_{3}-\tfrac{1}{8}\mathcal{F}^{(5)}_{4}+\tfrac{1}{4}\mathcal{F}^{(5)}_{5}+\tfrac{1}{2}\mathcal{F}^{(5)}_{6}-\tfrac{3}{5}\mathcal{F}^{(5)}_{7}+\tfrac{1}{2}\mathcal{F}^{(5)}_{8}+2\mathcal{F}^{(5)}_{9}-2\mathcal{F}^{(5)}_{10}\\
&-\tfrac{1}{2}\mathcal{F}^{(5)}_{11}+2\mathcal{F}^{(5)}_{12}+4\mathcal{F}^{(5)}_{13}-\tfrac{1}{6}\mathcal{F}^{(5)}_{14}-\mathcal{F}^{(5)}_{15}+\mathcal{F}^{(5)}_{16}+\tfrac{1}{6}\mathcal{F}^{(5)}_{17}-3\mathcal{F}^{(5)}_{18}+\mathcal{F}^{(5)}_{19}\\
&-3\mathcal{F}^{(5)}_{20}-4\mathcal{F}^{(5)}_{21}+2\mathcal{F}^{(5)}_{22}-2\mathcal{F}^{(5)}_{23}+2\mathcal{F}^{(5)}_{24}+\mathcal{F}^{(5)}_{25}+\mathcal{F}^{(5)}_{26}+\tfrac{2}{5}\mathcal{F}^{(5)}_{27}-\mathcal{F}^{(5)}_{28}
)
\end{aligned}
}$$

\newpage
