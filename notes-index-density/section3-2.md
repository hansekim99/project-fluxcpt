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

\paragraph{${ \boldsymbol{d=1} }$}
\begin{align*}
\mathrel{\vcenter{\hbox{\begin{tikzpicture}
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
\end{tikzpicture}}}}\hspace{.3cm}D^4 &\to &&\hspace{.3cm}\mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (0,0.5) -- (1,-1);
    \strand (0,-0.5) -- (1,1);
    \strand (0,-0.5) -- (0,0.5);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (0,0.5) circle [radius=1pt] (0,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\draw (0,-1) .. controls (-0.5,-0.5) and (-0.5,0.5) .. (0,1);
\end{tikzpicture}}}}\hspace{.3cm} D _{4,1,n} &&= \hspace{.3cm}\mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (1,1) -- (1,-1);
    \strand (0,-1) -- (0,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt]  (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.3cm} D _{4,1,n} \\
\hspace{.3cm}\mathrel{\vcenter{\hbox{\begin{tikzpicture}
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
}}}\hspace{.2cm}D^2 &\to (1)&&\hspace{.3cm}\mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (-0.3,0.5) -- (1,-1);
    \strand (-0.3,-0.5) -- (0,1);
    \strand (0.3,-0.5) -- (1,1);
    \strand (0,-1) -- (0.3,0.5);
    \strand (-0.3,-0.5) -- (-0.3,0.5);
    \strand (0.3,-0.5) -- (0.3,0.5);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.3,0.5) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (-0.3,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0.3,0.5) circle [radius=1pt];
\end{tikzpicture}
}}}\hspace{.2cm}D _{2,2,n}&&=\hspace{.3cm}\mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (1,1) -- (0,-1);
    \strand (1,-1) -- (0,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt]  (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.3cm} D_{2,2,n}\\
&\to(2) &&\hspace{.3cm}\mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (-0.3,0.5) -- (1,-1);
    \strand (-0.3,-0.5) -- (0,1);
    \strand (0.3,-0.5) -- (1,1);
    \strand (0,-1) -- (0.3,0.5);
    \strand (-0.3,-0.5) -- (0.3,0.5);
    \strand (0.3,-0.5) -- (-0.3,0.5);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.3,0.5) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (-0.3,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0.3,0.5) circle [radius=1pt];
\end{tikzpicture}
}}}\hspace{.2cm}D _{2,2,n}&&=\hspace{.3cm}\mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (1,1) -- (1,-1);
    \strand (0,-1) -- (0,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt]  (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.3cm} D _{2,2,n}\\
\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
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
\end{tikzpicture}}}}\hspace{.2cm}D ^{2}&\to (1)&&\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (-0.3, 0.5);
    \strand (1,-1) -- (0.3,0.5);
    \strand (-0.3,-0.5) -- (-1,1);
    \strand (0.3,-0.5) -- (1,1);
    \strand (0,-1) -- (0,1);
    \strand (-0.3,-0.5) -- (-0.3,0.5);
    \strand (0.3,0.5) -- (0.3,-0.5);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.3,0.5) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (0.3,0.5) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (-0.3,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.2cm}D _{2,2,n} &&=\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (1,1) -- (1,-1);
    \strand (0,-1) -- (0,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt]  (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.3cm} D _{2,2,n}\\
&\to (2)&&\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (-0.3, 0.5);
    \strand (1,-1) -- (0.3,0.5);
    \strand (-0.3,-0.5) -- (-1,1);
    \strand (0.3,-0.5) -- (1,1);
    \strand (0,-1) -- (0,1);
    \strand (-0.3,-0.5) -- (0.3,0.5);
    \strand (-0.3,0.5) -- (0.3,-0.5);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.3,0.5) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (0.3,0.5) circle [radius=1pt] (0.3,-0.5) circle [radius=1pt] (-0.3,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.2cm}D _{2,2,n} &&= \hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (1, -1);
    \strand (-1,-1) -- (1,1);
    \strand (0,-1) -- (0,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt]  (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.3cm} D _{2,2,n}\\
\mathrel{\vcenter{\hbox{\begin{tikzpicture}
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
\end{tikzpicture}}}}\hspace{.2cm}
\hspace{.2cm}&\to (1)&&\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (-0.5,0.5);
    \strand (-1,1) -- (-0.5,-0.5);
    \strand (0.5,-0.5) -- (1,1);
    \strand (0.5,0.5) -- (1,-1);
    \strand (-0.5,-0.5) -- (-0.5,0.5);
    \strand (0,-0.5) -- (0,0.5);
    \strand (0.5,-0.5) -- (0.5,0.5);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.5,0.5) 
circle [radius=1pt] (0.5,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt] (-0.5,-0.5) circle [radius=1pt] (0.5,0.5) circle [radius=1pt] (0,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0,0.5) circle [radius=1pt];
\draw (0,-1) .. controls (-0.5, -0.75) and (-0.5, 0.25) .. (0, 0.5);
\draw (0,-0.5) .. controls (0.5,-0.25) and (0.5,0.75) .. (0,1);
\end{tikzpicture}}}}\hspace{.2cm} &&=\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (1,1) -- (1,-1);
    \strand (0,-1) -- (0,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt]  (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.3cm} \\
&\to (1)&&\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (-0.5,0.5);
    \strand (-1,1) -- (-0.5,-0.5);
    \strand (0.5,-0.5) -- (1,1);
    \strand (0.5,0.5) -- (1,-1);
    \strand (-0.5,-0.5) -- (0,0.5);
    \strand (0,-0.5) -- (-0.5,0.5);
    \strand (0.5,-0.5) -- (0.5,0.5);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.5,0.5) 
circle [radius=1pt] (0.5,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt] (-0.5,-0.5) circle [radius=1pt] (0.5,0.5) circle [radius=1pt] (0,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0,0.5) circle [radius=1pt];
\draw (0,-1) .. controls (-0.5, -0.75) and (-0.5, 0.25) .. (0, 0.5);
\draw (0,-0.5) .. controls (0.5,-0.25) and (0.5,0.75) .. (0,1);
\end{tikzpicture}}}}\hspace{.2cm} &&=\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (0,1) -- (-1, -1);
    \strand (1,1) -- (1,-1);
    \strand (0,-1) -- (-1,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt]  (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.3cm} \\
&\to (1)&&\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (-0.5,0.5);
    \strand (-1,1) -- (-0.5,-0.5);
    \strand (0.5,-0.5) -- (1,1);
    \strand (0.5,0.5) -- (1,-1);
        \strand (-0.5,-0.5) -- (-0.5,0.5);
    \strand (0,-0.5) -- (0.5,0.5);
    \strand (0.5,-0.5) -- (0,0.5);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.5,0.5) 
circle [radius=1pt] (0.5,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt] (-0.5,-0.5) circle [radius=1pt] (0.5,0.5) circle [radius=1pt] (0,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0,0.5) circle [radius=1pt];
\draw (0,-1) .. controls (-0.5, -0.75) and (-0.5, 0.25) .. (0, 0.5);
\draw (0,-0.5) .. controls (0.5,-0.25) and (0.5,0.75) .. (0,1);
\end{tikzpicture}}}}\hspace{.2cm} &&=\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,1) -- (-1, -1);
    \strand (1,1) -- (0,-1);
    \strand (1,-1) -- (0,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt]  (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.3cm}\\
&\to (1)&&\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (-0.5,0.5);
    \strand (-1,1) -- (-0.5,-0.5);
    \strand (0.5,-0.5) -- (1,1);
    \strand (0.5,0.5) -- (1,-1);    
    \strand (-0.5,-0.5) -- (0.5,0.5);
    \strand (0,-0.5) -- (0,0.5);
    \strand (0.5,-0.5) -- (-0.5,0.5);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.5,0.5) 
circle [radius=1pt] (0.5,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt] (-0.5,-0.5) circle [radius=1pt] (0.5,0.5) circle [radius=1pt] (0,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0,0.5) circle [radius=1pt];
\draw (0,-1) .. controls (-0.5, -0.75) and (-0.5, 0.25) .. (0, 0.5);
\draw (0,-0.5) .. controls (0.5,-0.25) and (0.5,0.75) .. (0,1);
\end{tikzpicture}}}}\hspace{.2cm} &&=\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (1,1);
    \strand (0,-1) -- (0,1);
    \strand (1,-1) -- (-1,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt]  (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.3cm} \\
&\to (1)&&\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (-0.5,0.5);
    \strand (-1,1) -- (-0.5,-0.5);
    \strand (0.5,-0.5) -- (1,1);
    \strand (0.5,0.5) -- (1,-1);
        \strand (-0.5,-0.5) -- (0.5,0.5);
    \strand (0,-0.5) -- (-0.5,0.5);
    \strand (0.5,-0.5) -- (0,0.5);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.5,0.5) 
circle [radius=1pt] (0.5,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt] (-0.5,-0.5) circle [radius=1pt] (0.5,0.5) circle [radius=1pt] (0,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0,0.5) circle [radius=1pt];
\draw (0,-1) .. controls (-0.5, -0.75) and (-0.5, 0.25) .. (0, 0.5);
\draw (0,-0.5) .. controls (0.5,-0.25) and (0.5,0.75) .. (0,1);
\end{tikzpicture}}}}\hspace{.2cm} &&=\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (0,1);
    \strand (0,-1) -- (1,1);
    \strand (1,-1) -- (-1,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt]  (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.3cm}\\
&\to (1)&&\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (-0.5,0.5);
    \strand (-1,1) -- (-0.5,-0.5);
    \strand (0.5,-0.5) -- (1,1);
    \strand (0.5,0.5) -- (1,-1);
        \strand (-0.5,-0.5) -- (0,0.5);
    \strand (0,-0.5) -- (0.5,0.5);
    \strand (0.5,-0.5) -- (-0.5,0.5);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt] (-0.5,0.5) 
circle [radius=1pt] (0.5,-0.5) circle [radius=1pt] (0,1) circle [radius=1pt] (-0.5,-0.5) circle [radius=1pt] (0.5,0.5) circle [radius=1pt] (0,-0.5) circle [radius=1pt] (0,-1) circle [radius=1pt] (0,0.5) circle [radius=1pt];
\draw (0,-1) .. controls (-0.5, -0.75) and (-0.5, 0.25) .. (0, 0.5);
\draw (0,-0.5) .. controls (0.5,-0.25) and (0.5,0.75) .. (0,1);
\end{tikzpicture}}}}\hspace{.2cm} &&=\hspace{.3cm} \mathrel{\vcenter{\hbox{\begin{tikzpicture}
\begin{knot}[
    consider self intersections,
    clip width=3,
    flip crossing=2]
    \strand (-1,-1) -- (1,1);
    \strand (0,-1) -- (-1,1);
    \strand (1,-1) -- (0,1);
\end{knot}
\filldraw [black] (-1,1) circle [radius=1pt] (1,-1) circle [radius=1pt] (1,1) circle [radius=1pt] (-1,-1) circle [radius=1pt]  (0,-1) circle [radius=1pt] (0,1) circle [radius=1pt];
\end{tikzpicture}}}}\hspace{.3cm} 
\end{align*}


### Contracting Outer Indices : Classification of Diagrams into Invariants of ${ U(n) }$

The resulting diagrams are bipartite graphs, and represent invariants of ${ U(n) }$ that arise from contracting ${ d }$ copies of ${ \mathcal{F} }$ and ${ \bar{\mathcal{F}} }$. A large number of permutations create equivalent invariants; first note the term of the form 
$${ \mathcal{F}_{I_1J_1K_1}\cdots \mathcal{F}_{I_dJ_dK_d} \bar{\mathcal{F}}^{I^{\prime}_{1} J^{\prime}_{1}K^{\prime}_{1} } \cdots \bar{\mathcal{F}} ^{I^{\prime}_{d}J^{\prime}_{d}K^{\prime}_{d}} }$$
is identical under permutations between any ${ (I_i,J_i,K_i) }$ and ${ (I_j,J_j,K_j) }$, between any ${ (I^{\prime}_{i},J^{\prime}_{i},K^{\prime}_{i}) }$ and ${ (I^{\prime}_{j},J^{\prime}_{j},K^{\prime}_{j}) }$, or within any ${ (I_i,J_i,K_i) }$ or ${ (I^{\prime}_{i},J^{\prime}_{i},K^{\prime}_{i}) }$ (due to total symmtery of ${ \mathcal{F} }$). 
<!--The permutations within each block (within ${ (I_i,J_i,K_i) }$ and vice versa) can be factored out by representing the bipartite graph as a ${ d \times d }$ matrix, defining ${ a _{ij} }$ as the number of elements mapped to block ${ j }$ from block ${ i }$.-->

The permutations within each block (within ${ (I_i,J_i,K_i) }$ and vice versa) can be factored out by representing the bipartite graph as a weighted bipartite graph of ${ d  }$ upper and lower vertices each, where the edge connecting the lower vertex ${ i }$ and upper vertex ${ j }$ are the number of elements mapped from block ${ i }$ to block ${ j }$ in the initial graph.
$${\tikz[baseline=-0.5ex, x=6ex, y=6ex]{
  \path[use as bounding box] (-0.5,-1.6) rectangle (14.5,1.6);
  \begin{knot}[consider self intersections, clip width=2, flip crossing=2]
    \strand (0,-1) -- (11,1);
    \strand (2,-1) -- (5,1);
    \strand (3,-1) -- (3,1);
    \strand (5,-1) -- (2,1);
    \strand (6,-1) -- (6,1);
    \strand (8,-1) -- (1,1);
    \strand (9,-1) -- (9,1);
    \strand (11,-1) -- (0,1);
    \strand (1,-1) -- (10,1);
    \strand (4,-1) -- (7,1);
    \strand (7,-1) -- (4,1);
    \strand (10,-1) -- (8,1);
    \strand (12,-1) -- (12,1);
    \strand (13,-1) -- (13,1);
    \strand (14,-1) -- (14,1);
  \end{knot}
  \foreach \P in {(0,-1),(11,1),(2,-1),(5,1),(3,-1),(3,1),
                  (5,-1),(2,1),(6,-1),(6,1),(8,-1),(1,1),
                  (9,-1),(9,1),(11,-1),(0,1),(1,-1),(10,1),
                  (4,-1),(7,1),(7,-1),(4,1),(10,-1),(8,1),
                  (12,-1),(12,1),(13,-1),(13,1),(14,-1),(14,1)}
    \filldraw \P circle[radius=1pt];
}}$$
For example, the diagram above can be represents a ${ d=5 }$ diagram, which would be represented by the graph
$${ \tikz[baseline=-0.5ex, x = 12ex, y = 9ex]{
    \foreach \P in {(0,-1),(0,1),(1,-1),(1,1),(2,-1),(2,1),(3,-1),(3,1),(4,-1),(4,1)}
   \filldraw \P circle[radius=2pt];
   \draw (0,-1) -- (3,1) node [near start, fill = white] {2};
    \draw (1,-1) -- (0,1) node [near end, fill = white] {1};
    \draw (1,-1) -- (1,1) node [midway, fill = white] {1};
   \draw (1,-1) -- (2,1) node [near start, fill = white] {1};
   \draw (2,-1) -- (0,1) node [near end, fill = white] {1};
   \draw (2,-1) -- (1,1) node [near end, fill = white] {1};
   \draw (2,-1) -- (2,1) node [midway, fill = white] {1};
   \draw (3,-1) -- (0,1) node [near end, fill = white] {1};
   \draw (3,-1) -- (2,1) node [midway, fill = white] {1};
   \draw (3,-1) -- (3,1) node [midway, fill = white] {1};
   \draw (4,-1) -- (4,1) node [midway, fill = white] {3};
   \draw (0,-1) -- (1,1) node [near start, fill = white] {1};
} }$$

Then, the permutations between ${ \mathcal{F} }$ or ${ \bar{\mathcal{F}} }$ blocks correspond to graph isomorphisms. The problem of converting the initial bipartite graph into the weighted bipartite graph, and classifying each weighted bipartite graph by equivalence class with respect to graph isomorphisms is implemented using ``igraph``, and is at least a ${ \mathcal{O}(d!) }$ process taking up the majority of runtime.

It is still possible that different isomoprhism classes of the weighted bipartite graphs are equivalent scalar invariants, but a further reduction is not explored in this article. Numerically, it can be shown that the two invariants for ${ d=2 }$ and five invariants for ${ d=3 }$ differ given a sufficiently general setup (for example, ${ \mathcal{F}_{IJK} }$ given by ${ [1,1,1,6,9] }$). 

$Finally, notation to denote these graphs are introduced, as a matrix where each column denotes the connection weights for each lower node. Then the diagram above would be in the form
$${ ((0,1,0,2,0),(1,1,1,0,0),(1,1,1,0,0),(1,0,1,1,0),(0,0,0,0,3)). }$$

\newpage
