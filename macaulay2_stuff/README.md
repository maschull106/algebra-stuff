# Points in double nested Quot schemes and tangent spaces

This repository provides abstractions to represent points in a double nested Quot scheme of points and includes a method for computing the associated tangent space. The implementation is done in **Macaulay2**.


## Theory

Let $X$ be a quasi-projective variety and $E$ a coherent sheaf on $X$.
Below are key concepts needed to define the Double nested Quot scheme of points.

> ### Young diagram
> A *Young diagram* $\lambda$ is a finite collection of boxes arranged in descending left-aligned rows with nonincreasing length.
> It can be described by the lengths of its rows (which correspond to a partition).
> 
> **Example:**\
> For $\lambda = (5,4,4,2,1)$, the Young diagram is:
> <pre>
> <div style="line-height: 1.2;">
> ┌─┬─┬─┬─┬─┐
> ├─┼─┼─┼─┼─┘
> ├─┼─┼─┼─┤
> ├─┼─┼─┴─┘
> ├─┼─┘
> └─┘
> </div>
> </pre>



> ### Reverse plane partition
> Let $\lambda$ be a Young diagram.
> A *reverse plane partition* of shape $\lambda$ is a filling $\newcommand{\bn}{\boldsymbol{n}} \bn = (\bn_\square)_{\square\in\lambda}$ of the Young diagram with integers, such that the entries are nondecreasing along rows and columns.
>
> **Example:**\
> For $\lambda=(4,2,2)$, a reverse plane partition is:
> <pre>
> <div style="line-height: 1.2;">
> ┌───┬───┬───┬───┐
> | 0 | 2 | 3 | 5 |
> ├───┼───┼───┴───┘
> | 1 | 4 |
> ├───┼───┤
> | 2 | 4 |
> └───┴───┘
> </div>
> </pre>




> ### Double nested Quot scheme
> Let $\lambda$ be a Young diagram and $\boldsymbol{n}$ a reverse plane partition of shape $\lambda$.
> The *double nested Quot scheme of points* with underlying reverse plane partition $\boldsymbol{n}$ parameterises quotients $(E\rightarrow T_\square)_{\square\in\lambda}$ nested along the structure of $\lambda$, as follows
> $$
> \begin{matrix}
>    T_{0,0} & \leftarrow & T_{1,0} & \leftarrow & T_{2,0} & \leftarrow & T_{3,0} & \leftarrow & \dots\\
>    \uparrow && \uparrow && \uparrow && \uparrow\\
>    T_{0,1} & \leftarrow & T_{1,1} & \leftarrow & T_{2,1} & \leftarrow & T_{3,1} & \leftarrow & \dots\\
>    \uparrow && \uparrow && \uparrow && \uparrow\\
>   T_{0,2} & \leftarrow & T_{1,2} & \leftarrow & T_{2,2} & \leftarrow & \ddots\\
>    \uparrow && \uparrow && \uparrow\\
>   \vdots && \vdots && \vdots
> \end{matrix}
> $$
> where $T_\square$ is zero-dimensional of length $\boldsymbol{n}_\square$ for every $\square\in\lambda$.


The current code implementation focuses on the affine case $X = \mathrm{Spec}\, R$.

<!-- TODO: tangent space -->



## Code description


### How to create objects

A point in a double Quot scheme is represented as a directed graph, with a vertex for every 
quotient and an edge between quotients corresponding to adjacent boxes of $\lambda$,
in the direction opposite to the arrow.

The Macaulay2 class `DoubleNestedQuotSchemePoint` creates the abstraction for representing 
those objects.
There are several methods for instantiating the `DoubleNestedQuotSchemePoint` class.
The first one is illustrated below (we will show later 
other ones that simplify the process).\
We define a point nested as 
$$
\begin{matrix}
    T_0 & \leftarrow & T_1\\
    \uparrow && \uparrow\\
    T_2 & \leftarrow & T_3.
\end{matrix}
$$

The specific point we will define is
$$
\begin{matrix}
    && \begin{pmatrix}0\end{pmatrix}\\
    &0 & \leftarrow & \frac{k[x]}{x}\\
    \begin{pmatrix}0\end{pmatrix} & \uparrow && \uparrow & \begin{pmatrix}1&0\end{pmatrix}\\
    &\frac{k[x]}{x-1} & \leftarrow & \frac{k[x]}{x} \oplus \frac{k[x]}{x-1}\\
    && \begin{pmatrix}0&1\end{pmatrix}\\
    &&&&\nwarrow & \begin{pmatrix}1&0\\0&1\end{pmatrix}\\
    &&&&& k[x]^2
\end{matrix}
$$
where the matrices for every arrow describe the maps.

> **Remark.** Although conceptually we work over an algebraically closed field $k$, the computations are done over $\mathbb{Q}$ due to how Macaulay2 works.


```haskell
-- load the double_nested_quot_scheme.m2 file (which must be in the same directory)
load "double_nested_quot_scheme.m2"

-- creation of the ring R and an R-module F for defining Quot schemes
R = QQ[x]
E = R^2

-- two points P,Q in X = Spec R and the zero module (1_R is the element 1 in the ring R)
zeroModule = R^1/(1_R)
P = R^1/x
Q = R^1/(x-1)

-- the four modules that will constitute the nesting
T0 = zeroModule
T1 = P
T2 = Q
T3 = P++Q

-- the nested quotients T3 -> T2, T3 -> T2, T1 -> T0 and T2 -> T0
-- note: idMat(n, R) returns the identity n*n matrix over the ring R
-- note2: for a matrix M, the syntax M**R allows interpreting M as a matrix over R
f31 = map(T1, T3, matrix{{1,0}}**R)
f32 = map(T2, T3, matrix{{0,1}}**R)
f10 = map(T0, T1, idMat(1, R))
f20 = map(T0, T2, idMat(1, R))

-- the quotient maps from F to the four modules
q3 = map(T3, E, idMat(2, R))
q1 = map(T1, E, f31*q3)
q2 = map(T2, E, f32*q3)
q0 = map(T0, E, f10*q1)

-- the graph nodes constructing the point of the double nested Quot scheme
node3 = quotNode(q3)
node2 = quotNode(q2, Right=>nodeInfo(node3, f32))
node1 = quotNode(q1, Down=>nodeInfo(node3, f31))
node0 = quotNode(
    q0, 
    Right=>nodeInfo(node1, f10),
    Down=>nodeInfo(node2, f20)
)

-- the double nested Quot scheme point is constructed from the node corresponding to the unique minimal element of the Young diagram
quotPoint = doubleNestedQuotSchemePoint(node0)
Tspace = tangentSpace quotPoint
```


The above method for instantiating the `DoubleNestedQuotSchemePoint` class is fastidious.
We now describe a more compact way to define the same point.

We consider nested $R$-modules, quotients of another module $E$.
The modules and quotient maps are arranged in a matrix, call it $M$.
The entries $M_{i,j}$ are
- If $i,j \equiv 0 \pmod 2\:$ then $M_{i,j}$ is an $R$-module. Call the module $T$.
If $T$ corresponds to a maximal element of the underlying Young diagram, then instead of
just the module, a tuple must be given, containing first the module $T$ and second a matrix
describing the quotient $E \twoheadrightarrow T$.
- Otherwise, $M_{i,j}$ is a matrix describing the quotient map between the two adjacent modules in $M$.

The class `DoubleNestedQuotSchemePoint` can be instantiated by calling the `doubleNestedQuotSchemePoint` method with the module $E$ and the matrix $M$ described above.

As an example, the same point as before defined with this method.

```haskell
-- load the double_nested_quot_scheme.m2 file (which must be in the same directory)
load "double_nested_quot_scheme.m2"

-- creation of the ring R and an R-module F for defining Quot schemes
R = QQ[x]
E = R^2

-- two points P,Q in X = Spec R and the zero module (1_R is the element 1 in the ring R)
z = R^1/(1_R)
P = R^1/x
Q = R^1/(x-1)

quotPoint = doubleNestedQuotSchemePoint(
    E,
    {
        {z,     {{0}},      Q},
        {{{0}},             {{0,1}}},
        {P,     {{1,0}},    (P++Q, idMat(2))}
    }
)
```

> **Remark.** Note that the second row of the matrix (with two matrices describing maps) only contains two elements, whereas the firt and third row contain three elements. This is because only two of the tree elements in the latter rows correspond to modules, hence need maps defined between them.

> **Remark.** This whole construction is cosmetic. If it is formated with the appropriate spaces between entries of the big matrix, the different quotient modules and the maps between them can easily be read.

Another example of instantiation of `DoubleNestedQuotSchemePoint`, for the reverse plane partition
<pre>
<div style="line-height: 1.2;">
┌───┬───┬───┐
| 0 | 1 | 2 |
├───┼───┼───┤
| 1 | 2 | 3 |
├───┼───┼───┤
| 2 | 3 | 4 |
└───┴───┴───┘
</div>
</pre>
is given below.


```haskell
load "double_nested_quot_scheme.m2"

R = QQ[x]
z = R^1/(ideal(1)**R)
P = R^1/x
P2 = R^1/x^2

quotPoint = doubleNestedQuotSchemePoint(
    R^2,
    {
        {z,         {{0}},      P,         {{1,0}},         P++P},
        {{{0}},                 {{1}},                      idMat(2)},
        {P,         {{1}},      P2,        {{1,0}},         P2++P},
        {{{1,0}},               {{1,0}},                    idMat(2)},
        {P++P,      idMat(2),   P2++P,      idMat(2),       (P2++P2, idMat(2))}
    }
)
```

<!-- TODO: special cases of nested Quot schemes and nested Hilbert schemes -->