-- Some verified examples where dimension of the tangent space is known from other sources or by hand computation 
-- to check code is still working well after each modifiction

needs "double_nested_quot_scheme.m2"

count = 0


-- Cheah's nested Hilbert schemes singularity examples

count = count+1
R = QQ[x,y]
I1 = ideal(x, y); I2 = ideal(x^2, y); I3 = ideal(x^2, x*y, y^2)
F = R^1
T1 = F/I1; T2 = F/I2; T3 = F/I3
quotPoint1 = nestedHilbSchemePoint2({F, T3, T2, T1})
T = tangentSpace quotPoint1
assert (dim T == 7)


count = count+1
quotPoint2 = nestedHilbSchemePoint2({F, T3, T1})
T = tangentSpace quotPoint2
assert (dim T == 8)


count = count+1
R = QQ[x,y,z]
I1 = ideal(x^2, y^2, z^2, x*y, x*z, y*z); I2 = ideal(x^2, y^2, z^2)
F = R^1
T1 = F/I1; T2 = F/I2
quotPoint3 = nestedHilbSchemePoint2({F, T2, T1})
T = tangentSpace quotPoint3
assert (dim T == 30)




-- Singularity example for nested Quot schemes by S. Monavari and A. T. Ricolfi

count = count+1
R = QQ[x, y]
F = R^2
P = R^1/ideal(x, y)
T2 = P ++ P
T1 = P
p2 = map(T2, F, matrix{{1,0},{0,1}}**R)
f = map(T1, T2, matrix{{1, 1}}**R)
quotPoint4 = nestedQuotSchemePoint2({p2, f})
T = tangentSpace quotPoint4
assert (dim T == 7)






-- Double nested Hilbert scheme on the 2*2 square diagram singular and non-singular examples

R = QQ[x]
one = R^1/(ideal(1)**R); P = R^1/ideal(x); Q = R^1/ideal(x-1); P2 = R^1/ideal(x^2); PQ = R^1/ideal(x*(x-1));


count = count+1
T = tangentSpace doubleNestedQuotSchemePoint(
    R^1,
    {
        {one,    {{1}},    Q},
        {{{1}},           {{0,1}}},
        {P,    {{1,0}},    (P++Q, {{1},{1}})}
    }
)
assert (dim T == 2)

-- same as last one but different representation
count = count+1
T = tangentSpace doubleNestedQuotSchemePoint(
    R^1,
    {
        {one,    {{1}},    Q},
        {{{1}},           {{1}}},
        {P,    {{1}},    (PQ, {{1}})}
    }
)
assert (dim T == 2)

count = count+1
T = tangentSpace doubleNestedQuotSchemePoint(
    R^1,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1}}},
        {P,    {{1}},    (PQ, {{1}})}
    }
)
assert (dim T == 2)

count = count+1
T = tangentSpace doubleNestedQuotSchemePoint(
    R^1,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1}}},
        {P,    {{1}},    (P2, {{1}})}
    }
)
assert (dim T == 3)



-- Double nested Quot scheme on the 2*2 square diagram singular and non-singular examples (TODO: verify them more thoroughly by hand)
-- ┌───┬───┐
-- | 0 | 1 |
-- ├───┼───┤
-- | 1 | 2 |
-- └───┴───┘

count = count+1
T = tangentSpace doubleNestedQuotSchemePoint(
    R^2,
    {
        {one,    {{1}},    Q},
        {{{1}},           {{0,1}}},
        {P,    {{1,0}},    (P++Q, idMatrix(2))}
    }
)
assert (dim T == 4)

count = count+1
T = tangentSpace doubleNestedQuotSchemePoint(
    R^2,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1,0}}},
        {P,    {{1,0}},    (P++Q, idMatrix(2))}
    }
)
assert (dim T == 4)

count = count+1
T = tangentSpace doubleNestedQuotSchemePoint(
    R^2,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1,0}}},
        {P,    {{1,0}},    (P++P, idMatrix(2))}
    }
)
assert (dim T == 5)

count = count+1
T = tangentSpace doubleNestedQuotSchemePoint(
    R^2,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1,0}}},
        {P,    {{1,1}},    (P++P, idMatrix(2))}
    }
)
assert (dim T == 4)

count = count+1
T = tangentSpace doubleNestedQuotSchemePoint(
    R^2,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1}}},
        {P,    {{1}},    (P2, {{1,1}})}
    }
)
assert (dim T == 5)



-- Double nested Quot scheme on 
-- ┌───┬───┐
-- | 0 | 2 |
-- ├───┼───┤
-- | 2 | 4 |
-- └───┴───┘

count = count+1
T = tangentSpace doubleNestedQuotSchemePoint(
    R^2,
    {
        {one,         {{0,0}},            P++P},
        {{{0,0}},                       {{1,0},{0,1}}},
        {P++P,      {{1,0},{0,1}},      (P2++P2, idMatrix(2))}
    }
)
assert (dim T == 12)


count = count+1
T = tangentSpace doubleNestedQuotSchemePoint(
    R^2,
    {
        {one,         {{0}},            P2},
        {{{0}},                       {{1,0}}},
        {P2,      {{1,0}},      (P2++P2, idMatrix(2))}
    }
)
assert (dim T == 10)


count = count+1
B1 = matrix{{x,1},{0,x}}
B2 = matrix{{x,-1},{0,x}}
B3 = matrix{{x,-1},{0,x}}
B4 = matrix{{x,1},{0,x}}
rootMat = idMatrix(2, R)
n11 = quotKernelNode()
n01 = quotKernelNode(Right=>kernelNodeInfo(n11, B2))
n10 = quotKernelNode(Down=>kernelNodeInfo(n11, B4))
n00 = quotKernelNode(Right=>kernelNodeInfo(n10, B3), Down=>kernelNodeInfo(n01, B1))
T = tangentSpace doubleNestedQuotSchemePoint(rootMat, n00)
assert (dim T == 9)


count = count+1
p = x+1
B1 = matrix{{x^2,p},{0,1}}
B2 = matrix{{1,-p},{0,x^2}}
B3 = matrix{{1,0},{0,x^2}}
B4 = matrix{{x^2,0},{0,1}}
rootMat = idMatrix(2, R)
n11 = quotKernelNode()
n01 = quotKernelNode(Right=>kernelNodeInfo(n11, B2))
n10 = quotKernelNode(Down=>kernelNodeInfo(n11, B4))
n00 = quotKernelNode(Right=>kernelNodeInfo(n10, B3), Down=>kernelNodeInfo(n01, B1))
T = tangentSpace doubleNestedQuotSchemePoint(rootMat, n00)
assert (dim T == 8)




-- Double nested Quot scheme on 
-- ┌───┬───┬───┐
-- | 0 | 1 | 2 |
-- ├───┼───┼───┤
-- | 1 | 2 | 3 |
-- ├───┼───┼───┤
-- | 2 | 3 | 4 |
-- └───┴───┴───┘

count = count+1
T = tangentSpace doubleNestedQuotSchemePoint(
    R^2,
    {
        {one,       {{1}},      P,         {{1,0}},         P++P},
        {{{1}},                 {{1}},                      idMat(2)},
        {P,         {{1}},      P2,        {{1,0}},         P2++P},
        {{{1,0}},               {{1,0}},                    idMat(2)},
        {P++P,      idMat(2),   P2++P,      idMat(2),       (P2++P2, idMat(2))}
    }
)
-- assert (dim T == 11)
