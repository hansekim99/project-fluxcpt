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

## Diagrammatic Method for Order ${ \boldsymbol{d} }$
The combinatoric nature of the process above can be utilised to introduce a diagrammatic representation. For each one of the 4 terms in ${ N_I^P }$, the lower indices of ${ \mathcal{F} }$ are represented as "outer" lower dots, upper indices of ${ \bar{\mathcal{F}} }$ as "outer" upper dots, lower indices of ${ \phi  }$ as "inner" lower dots, and upper indices of ${ \overline{\phi } }$ as "inner" upper dots.

Contraction of two indices is represented as a line connecting the two corresponding dots. Dots corresponding to "external indices" (which will have no lines connected to them in the "noncontracted" state) are denoted as a circle for emphasis. Importantly, as ${ \mathcal{F}_{IJK} }$ is totally symmetric, the "outer" upper dots can be permuted in any order;

$${
\mathcal{F}_{IJK}\bar{\mathcal{F}}^{JPN} \overline{\phi }^{K}\phi _{N} =\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (0,0.5) -- (1,-1);
    \strand (0,-0.5) -- (1,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (0,0.5) circle [radius=1pt] (0,-0.5) circle [radius=1pt];
\draw (0,1) circle [radius=2pt] (0,-1) circle [radius=2pt];
\end{tikzpicture}}}}
\hspace{1cm} \mathcal{F}_{IJK}\bar{\mathcal{F}}^{JMN} \overline{\phi }^{K}\overline{\phi }^{P}\phi _{M} \phi _{N}D ^{-2}=\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (-0.3,0.5) -- (1,-1);
    \strand (-0.3,-0.5) -- (0,1);
    \strand (0.3,-0.5) -- (1,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.3,0.5) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (-0.3,-0.5) circle [radius=1pt];
\draw (0,-1) circle [radius=2pt] (0.3,0.5) circle [radius=2pt];
\end{tikzpicture}
}}}}\hspace{.2cm}D ^{-2}$$
$${\mathcal{F}_{IJK}\bar{\mathcal{F}}^{LPN} \overline{\phi }^{J}\overline{\phi }^{K}\phi _{L} \phi _{N}D ^{-2}=\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (-0.3, 0.5);
    \strand (1,-1) -- (0.3,0.5);
    \strand (-0.3,-0.5) -- (-1,1);
    \strand (0.3,-0.5) -- (1,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.3,0.5) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (0.3,0.5) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (-0.3,-0.5) circle [radius=1pt];
\draw (0,-1) circle [radius=2pt] (0,1) circle [radius=2pt];
\end{tikzpicture}}}}\hspace{.2cm}D ^{-2}
\hspace{1cm} \mathcal{F}_{IJK}\bar{\mathcal{F}}^{LMN} \overline{\phi }^{J}\overline{\phi }^{K}\overline{\phi }^{P}\phi _{L}\phi _{M} \phi _{N} D ^{-4}=\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (-0.5,0.5);
    \strand (-1,1) -- (-0.5,-0.5);
    \strand (0.5,-0.5) -- (1,1);
    \strand (0.5,0.5) -- (1,-1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.5,0.5) circle [radius=1pt] (0.5,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt] (-0.5,-0.5) circle [radius=1pt] (0.5,0.5) circle [radius=1pt] (0,-0.5) circle [radius=1pt];
\draw (0,-1) circle [radius=2pt] (0,0.5) circle [radius=2pt];
\draw (0,-0.5) .. controls (0.5,-0.25) and (0.5,0.75) .. (0,1);
\end{tikzpicture}}}}\hspace{.2cm} D ^{-4}
}$$

### Contracting External Indices
Calculating the order ${ d }$ contribution to the integrand can be done by drawing all combinations of ${ d }$ diagrams from the 4 in ${ N }$, and drawing all possible contractions of lower and upper circles. Only diagrams with singular degree up to ${ D ^{-4} }$ need to be considered.  

\subparagraph{${ \boldsymbol{d=1} }$}
$${-|X|^{2n-4} D^4\hspace{.3cm}\mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (0,0.5) -- (1,-1);
    \strand (0,-0.5) -- (1,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (0,0.5) circle [radius=1pt] (0,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\draw (0,-1) .. controls (-0.5,-0.5) and (-0.5,0.5) .. (0,1);
\end{tikzpicture}}}}\hspace{.3cm},\ -|X|^{2n-4}D^2 \hspace{.1cm}\mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (-0.3,0.5) -- (1,-1);
    \strand (-0.3,-0.5) -- (0,1);
    \strand (0.3,-0.5) -- (1,1);
    \strand (0,-1) -- (0.3,0.5);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.3,0.5) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (-0.3,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0.3,0.5) circle [radius=1pt];
\end{tikzpicture}
}}}\hspace{.3cm} ,\ -|X|^{2n-4}D^2 \hspace{.1cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (-0.3, 0.5);
    \strand (1,-1) -- (0.3,0.5);
    \strand (-0.3,-0.5) -- (-1,1);
    \strand (0.3,-0.5) -- (1,1);
    \strand (0,-1) -- (0,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.3,0.5) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (0.3,0.5) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (-0.3,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.3cm},\ -|X|^{2n-4} \hspace{.1cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (-0.5,0.5);
    \strand (-1,1) -- (-0.5,-0.5);
    \strand (0.5,-0.5) -- (1,1);
    \strand (0.5,0.5) -- (1,-1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.5,0.5) 
circle [radius=1pt] (0.5,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt] (-0.5,-0.5) circle [radius=1pt] (0.5,0.5) circle [radius=1pt] (0,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0,0.5) circle [radius=1pt];
\draw (0,-1) .. controls (-0.5, -0.75) and (-0.5, 0.25) .. (0, 0.5);
\draw (0,-0.5) .. controls (0.5,-0.25) and (0.5,0.75) .. (0,1);
\end{tikzpicture}}}}\hspace{.2cm} }$$
\subparagraph{${ \boldsymbol{d=2} }$}
$${
\dfrac{1}{2}|X|^{2n-6}D^4
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\draw (1,-1) .. controls (0.5,-0.5) and (0.5,0.5) .. (1,1);
\draw (4,-1) .. controls (3.5,-0.5) and (3.5,0.5) .. (4,1);
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0,1);
    \strand (2,-1) -- (1,0.5);
    \strand (1,-0.5) -- (2,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (4,0.5);
    \strand (4,-0.5) -- (5,1);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (1,0.5) circle [radius=1pt]
    (1,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3,1) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (4,0.5) circle [radius=1pt]
    (4,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (4,1) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm},
-\dfrac{1}{2}|X|^{2n-6}D^4
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0,1);
    \strand (2,-1) -- (1,0.5);
    \strand (1,-0.5) -- (2,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (4,0.5);
    \strand (4,-0.5) -- (5,1);
    \strand (1,-1) -- (4,1);
    \strand (4,-1) -- (1,1);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (1,0.5) circle [radius=1pt]
    (1,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3,1) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (4,0.5) circle [radius=1pt]
    (4,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (4,1) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (1,1) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm}
}$$
$${ 
|X|^{2n-6}D^2
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\draw (4,-1) .. controls (3.5,-0.5) and (3.5,0.5) .. (4,1);
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0,1);
    \strand (2,-1) -- (0.7,0.5);
    \strand (0.7,-0.5) -- (1,1);
    \strand (1.3,-0.5) -- (2,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (4,0.5);
    \strand (4,-0.5) -- (5,1);
    \strand (1,-1) -- (1.3,0.5);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (0.7,0.5) circle [radius=1pt]
    (0.7,-0.5) circle [radius=1pt] (1,1) circle [radius=1pt]
    (1.3,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3,1) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (4,0.5) circle [radius=1pt]
    (4,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (1.3,0.5) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (4,1) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm},
-|X|^{2n-6}D^2
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0,1);
    \strand (2,-1) -- (0.7,0.5);
    \strand (0.7,-0.5) -- (1,1);
    \strand (1.3,-0.5) -- (2,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (4,0.5);
    \strand (4,-0.5) -- (5,1);
    \strand (1,-1) -- (4,1);
    \strand (4,-1) -- (1.3,0.5);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (0.7,0.5) circle [radius=1pt]
    (0.7,-0.5) circle [radius=1pt] (1,1) circle [radius=1pt]
    (1.3,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3,1) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (4,0.5) circle [radius=1pt]
    (4,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (4,1) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (1.3,0.5) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm}
}$$
$${ |X|^{2n-6}D^2
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\draw (4,-1) .. controls (3.5,-0.5) and (3.5,0.5) .. (4,1);
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0.7,0.5);
    \strand (2,-1) -- (1.3,0.5);
    \strand (0.7,-0.5) -- (0,1);
    \strand (1.3,-0.5) -- (2,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (4,0.5);
    \strand (4,-0.5) -- (5,1);
    \strand (1,-1) -- (1,1);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0.7,0.5) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (1.3,0.5) circle [radius=1pt]
    (0.7,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt]
    (1.3,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3,1) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (4,0.5) circle [radius=1pt]
    (4,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (4,1) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm}
,-|X|^{2n-6}D^2
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0.7,0.5);
    \strand (2,-1) -- (1.3,0.5);
    \strand (0.7,-0.5) -- (0,1);
    \strand (1.3,-0.5) -- (2,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (4,0.5);
    \strand (4,-0.5) -- (5,1);
    \strand (1,-1) -- (4,1);
    \strand (4,-1) -- (1,1);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0.7,0.5) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (1.3,0.5) circle [radius=1pt]
    (0.7,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt]
    (1.3,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3,1) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (4,0.5) circle [radius=1pt]
    (4,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (4,1) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (1,1) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm}
 }$$
 $${ |X|^{2n-6}
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\draw (1,-0.5) .. controls (1.5,-0.125) and (1.5,0.625) .. (1,1);
\draw (1,-1) .. controls (0.5,-0.625) and (0.5,0.125) .. (1,0.5);
\draw (4,-1) .. controls (3.5,-0.5) and (3.5,0.5) .. (4,1);
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0.5,0.5);
    \strand (2,-1) -- (1.5,0.5);
    \strand (0.5,-0.5) -- (0,1);
    \strand (1.5,-0.5) -- (2,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (4,0.5);
    \strand (4,-0.5) -- (5,1);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0.5,0.5) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (1.5,0.5) circle [radius=1pt]
    (0.5,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt]
    (1,-0.5) circle [radius=1pt] (1,1) circle [radius=1pt]
    (1.5,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3,1) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (4,0.5) circle [radius=1pt]
    (4,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (1,0.5) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (4,1) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm},
-|X|^{2n-6}
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\draw (1,-0.5) .. controls (1.5,-0.125) and (1.5,0.625) .. (1,1);
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0.5,0.5);
    \strand (2,-1) -- (1.5,0.5);
    \strand (0.5,-0.5) -- (0,1);
    \strand (1.5,-0.5) -- (2,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (4,0.5);
    \strand (4,-0.5) -- (5,1);
    \strand (1,-1) -- (4,1);
    \strand (4,-1) -- (1,0.5);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0.5,0.5) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (1.5,0.5) circle [radius=1pt]
    (0.5,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt]
    (1,-0.5) circle [radius=1pt] (1,1) circle [radius=1pt]
    (1.5,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3,1) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (4,0.5) circle [radius=1pt]
    (4,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (4,1) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (1,0.5) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm}
 }$$
$${ \dfrac{1}{2}|X|^{2n-6}
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0,1);
    \strand (2,-1) -- (0.7,0.5);
    \strand (0.7,-0.5) -- (1,1);
    \strand (1.3,-0.5) -- (2,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (3.7,0.5);
    \strand (3.7,-0.5) -- (4,1);
    \strand (4.3,-0.5) -- (5,1);
    \strand (1,-1) -- (1.3,0.5);
    \strand (4,-1) -- (4.3,0.5);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (0.7,0.5) circle [radius=1pt]
    (0.7,-0.5) circle [radius=1pt] (1,1) circle [radius=1pt]
    (1.3,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3,1) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (3.7,0.5) circle [radius=1pt]
    (3.7,-0.5) circle [radius=1pt] (4,1) circle [radius=1pt]
    (4.3,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (1.3,0.5) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (4.3,0.5) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm},
-\dfrac{1}{2}|X|^{2n-6}
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0,1);
    \strand (2,-1) -- (0.7,0.5);
    \strand (0.7,-0.5) -- (1,1);
    \strand (1.3,-0.5) -- (2,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (3.7,0.5);
    \strand (3.7,-0.5) -- (4,1);
    \strand (4.3,-0.5) -- (5,1);
    \strand (1,-1) -- (4.3,0.5);
    \strand (4,-1) -- (1.3,0.5);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (0.7,0.5) circle [radius=1pt]
    (0.7,-0.5) circle [radius=1pt] (1,1) circle [radius=1pt]
    (1.3,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3,1) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (3.7,0.5) circle [radius=1pt]
    (3.7,-0.5) circle [radius=1pt] (4,1) circle [radius=1pt]
    (4.3,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (4.3,0.5) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (1.3,0.5) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm}
 }$$
$${ |X|^{2n-6}
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0.7,0.5);
    \strand (2,-1) -- (1.3,0.5);
    \strand (0.7,-0.5) -- (0,1);
    \strand (1.3,-0.5) -- (2,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (3.7,0.5);
    \strand (3.7,-0.5) -- (4,1);
    \strand (4.3,-0.5) -- (5,1);
    \strand (1,-1) -- (1,1);
    \strand (4,-1) -- (4.3,0.5);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0.7,0.5) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (1.3,0.5) circle [radius=1pt]
    (0.7,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt]
    (1.3,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3,1) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (3.7,0.5) circle [radius=1pt]
    (3.7,-0.5) circle [radius=1pt] (4,1) circle [radius=1pt]
    (4.3,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (4.3,0.5) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm},
-|X|^{2n-6}
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0.7,0.5);
    \strand (2,-1) -- (1.3,0.5);
    \strand (0.7,-0.5) -- (0,1);
    \strand (1.3,-0.5) -- (2,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (3.7,0.5);
    \strand (3.7,-0.5) -- (4,1);
    \strand (4.3,-0.5) -- (5,1);
    \strand (1,-1) -- (4.3,0.5);
    \strand (4,-1) -- (1,1);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0.7,0.5) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (1.3,0.5) circle [radius=1pt]
    (0.7,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt]
    (1.3,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3,1) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (3.7,0.5) circle [radius=1pt]
    (3.7,-0.5) circle [radius=1pt] (4,1) circle [radius=1pt]
    (4.3,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (4.3,0.5) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (1,1) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm}
 }$$
 $${ \dfrac{1}{2}|X|^{2n-6}
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0.7,0.5);
    \strand (2,-1) -- (1.3,0.5);
    \strand (0.7,-0.5) -- (0,1);
    \strand (1.3,-0.5) -- (2,1);
    \strand (3,-1) -- (3.7,0.5);
    \strand (5,-1) -- (4.3,0.5);
    \strand (3.7,-0.5) -- (3,1);
    \strand (4.3,-0.5) -- (5,1);
    \strand (1,-1) -- (1,1);
    \strand (4,-1) -- (4,1);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0.7,0.5) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (1.3,0.5) circle [radius=1pt]
    (0.7,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt]
    (1.3,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3.7,0.5) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (4.3,0.5) circle [radius=1pt]
    (3.7,-0.5) circle [radius=1pt] (3,1) circle [radius=1pt]
    (4.3,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (4,1) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm},
-\dfrac{1}{2}|X|^{2n-6}
\hspace{.1cm} \tikz[baseline=-0.5ex]{\begin{tikzpicture}
\begin{knot}[consider self intersections,
 clip width=3,flip crossing=2]
    \strand (0,-1) -- (0.7,0.5);
    \strand (2,-1) -- (1.3,0.5);
    \strand (0.7,-0.5) -- (0,1);
    \strand (1.3,-0.5) -- (2,1);
    \strand (3,-1) -- (3.7,0.5);
    \strand (5,-1) -- (4.3,0.5);
    \strand (3.7,-0.5) -- (3,1);
    \strand (4.3,-0.5) -- (5,1);
    \strand (1,-1) -- (4,1);
    \strand (4,-1) -- (1,1);
\end{knot}
\filldraw [black]
    (0,-1) circle [radius=1pt] (0.7,0.5) circle [radius=1pt]
    (2,-1) circle [radius=1pt] (1.3,0.5) circle [radius=1pt]
    (0.7,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt]
    (1.3,-0.5) circle [radius=1pt] (2,1) circle [radius=1pt]
    (3,-1) circle [radius=1pt] (3.7,0.5) circle [radius=1pt]
    (5,-1) circle [radius=1pt] (4.3,0.5) circle [radius=1pt]
    (3.7,-0.5) circle [radius=1pt] (3,1) circle [radius=1pt]
    (4.3,-0.5) circle [radius=1pt] (5,1) circle [radius=1pt]
    (1,-1) circle [radius=1pt] (4,1) circle [radius=1pt]
    (4,-1) circle [radius=1pt] (1,1) circle [radius=1pt];
\end{tikzpicture}}\hspace{.3cm}
}$$

### Contracting Inner Indices 
The integrand is integrated over all variables ${ \phi  }$ by performing Wick contraction; for each diagram from the previous step, draw all possible contractions of lower internal dots with upper internal dots, with all combinatoric coefficients carrying over. Then erase all internal dots to leave only connections between external dots. If there is ${ D^2  }$ or ${ D^4 }$ in the term, additional contact terms must be considered;
 $${ \begin{aligned}
\ev{\overline{\phi }^{A_1}\cdots \overline{\phi }^{A_m}\phi _{B_1}\cdots \phi _{B_m}D^2 } &= (\text{all contractions})\ev{D^2 } + \sum (\text{all contractions not involving }A_i,B_j)\ev{\partial _{A_i}\partial ^{B_j}D^2 }\\
&= (\text{all contractions})\cdot (-n+|X|^{2} ) - \sum (\text{all contractions not involving }A_i,B_j)\cdot \delta _{A_i}^{B_j}\\
&= (\text{all contractions})\cdot \left(-n-\dfrac{m^2 \cdot (m-1)!}{m!}+|X|^{2}\right)\\
&= (\text{all contractions})\cdot (-n-m+|X|^{2} )\coloneqq (\text{all contractions}) \cdot D_{2,m,n}
 \end{aligned}
 }$$
$${ \begin{aligned}
\ev{\overline{\phi }^{A_1}\cdots \overline{\phi }^{A_m}\phi _{B_1}\cdots \phi _{B_m} D^4} &= (\text{all contractions})\ev{D^4} \\
&+ \sum (\text{all contractions not involving } A_i, B_j) \ev{\partial _{A_i}\partial ^{B_j}D^4}\\
&+ \sum (\text{all contractions not involving } A_i, A_k,B_j,B_l) \ev{\partial _{A_i}\partial _{A_k}\partial ^{B_j}\partial ^{B_l}D^4}\\
&= (\text{all contractions}) \cdot \left[2 \cdot n + n(n-1) - 2n |X|^{2} + |X|^{4}\right]\\
&+ \sum (\text{all contractions not involving } A_i, B_j)\cdot \delta ^{A_i}_{B_j} \cdot (2n+2-2  |X|^{2} )\\
&+ \sum (\text{all contractions not involving } A_i, A_k,B_j,B_l) \cdot (\delta ^{A_i}_{B_j}\delta ^{A_k}_{B_l}+\delta ^{A_i}_{B_l}\delta ^{A_k}_{B_j})\\
&= (\text{all contractions}) \cdot [(n^2 +n -2n|X|^{2} +|X|^{4})+(2n+2-2|X|^{2} )m+m(m-1)]\\
&= (\text{all contractions}) \cdot [(n^2 +n + m(2n+2) + m(m-1) )+ (-2n-2m)|X|^{2} +|X|^{4}]\\
&\coloneqq  (\text{all contractions})\cdot D_{4,m,n}
\end{aligned}
 }$$
The contact terms are only multiplicative in nature, and be exchanged with ${ D^2  }$ or ${ D^4 }$ respectively. Then, collect all terms wrt. powers in ${ |X|^{2}  }$, then integrate over ${ X }$ with ${ \ev{|X|^{2n}} = n! }$.
\newpage

