Code providing abstractions capable of representing a point in a double nested Quot scheme of points, together with a method for computing the associated tangent space.


## Theory

Let $X$ be a quasi-projective variety and $E$ a coherent sheaf on $X$.

> ### Young diagram
> A Young diagram $\lambda$ is a finite collection of boxes arranged in descending left-aligned rows with nonincreasing length.
> It can be described by the length of its rows (the associated partition).
> For instance, $\lambda = (5,4,4,2,1)$ is
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
> A reverse plane partition of shape $\lambda= (4,2,1)$.
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
> Parameterises quotients $(E\rightarrow T_\square)_{\square\in\lambda}$ nested along the structure of a Young diagram $\lambda$, as follows
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


In pratcite, the code only treats the affine case $X = \mathrm{Spec}\, R$.

TODO: tangent space



## Code description

Implementation in Macaulay2.

### How to create objects

A point in a double Quot scheme is represented as a graph, with a vertex for every quotient and an edge between quotients corresponding to adjacent boxes of $\lambda$.

The most basic method for creating such an object is illustrated below (we will show later more abstractions that simplify this process).\
We define a point nested as 
$$
\begin{matrix}
    T_0 & \leftarrow & T_1\\
    \uparrow && \uparrow\\
    T_2 & \leftarrow & T_3.
\end{matrix}
$$

```haskell
-- load the double_nested_quot_scheme.m2 file (which must be in the same directory)
load "double_nested_quot_scheme.m2"

-- creation of the ring R and an R-module F for defining Quot schemes
R = QQ[x]
F = R^2

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
q3 = map(T3, F, idMat(2, R))
q1 = map(T1, F, f31*q3)
q2 = map(T2, F, f32*q3)
q0 = map(T0, F, f10*q1)

-- the graph nodes constructing the point of the double nested Quot scheme
node3 = quotNode(q3)
node2 = quotNode(q2, Right=>nodeInfo(node3, f32))
node1 = quotNode(q1, Down=>nodeInfo(node3, f31))
node0 = quotNode(
    q0, 
    Right=>nodeInfo(node1, f10),
    Down=>nodeInfo(node2, f20)
)
quotPoint = doubleNestedQuotSchemePoint(node0)
Tspace = tangentSpace quotPoint
```



