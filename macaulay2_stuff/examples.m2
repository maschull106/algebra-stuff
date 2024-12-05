needs "double_nested_quot_scheme.m2"


R = QQ[x]
F = R^2
zeroo = R^1/(ideal(1)**R)
P = R^1/ideal(x)
S = R^1/ideal(x-1)
P2 = R^1/ideal(x^2)
P3 = R^1/x^3
P4 = R^1/x^4
PS = R^1/ideal(x*(x-1))


T0 = tangentSpace doubleNestedQuotSchemePoint(
    F,
    {
        {zeroo,    {{1}},    P},
        {{{1}},           {{1,0}}},
        {P,    {{1,0}},    (P++P, idMatrix(2))}
    }
)

T1 = tangentSpace doubleNestedQuotSchemePoint(
    F,
    {
        {zeroo,         {{0,0}},            P++P},
        {{{0,0}},                       {{1,0},{0,1}}},
        {P++P,      {{1,0},{0,1}},      (P2++P2, idMatrix(2))}
    }
)

T2 = tangentSpace doubleNestedQuotSchemePoint(
    F,
    {
        {P,         {{1,0}},            P++P},
        {{{1}},                       {{1,0},{0,1}}},
        {P2,      {{1,0}},      (P2++P, idMatrix(2))}
    }
)

T3 = tangentSpace doubleNestedQuotSchemePoint(
    F,
    {
        {zeroo,    {{1}},    P},
        {{{1}},           {{-1,1}}},
        {P3,    {{-1,x^2+x}},    (P4++P2, {{x,-1},{1,0}})}
    }
)


B1 = matrix{{x,1},{0,x}}
B2 = matrix{{x,0},{0,x}}
B3 = matrix{{x,0},{0,x}}
B4 = matrix{{x,2},{0,x}}
rootMat = idMatrix(2, R)
n11 = quotKernelNode()
n01 = quotKernelNode(Right=>kernelNodeInfo(n11, B2))
n10 = quotKernelNode(Down=>kernelNodeInfo(n11, B4))
n00 = quotKernelNode(Right=>kernelNodeInfo(n10, B3), Down=>kernelNodeInfo(n01, B1))
quotPoint = doubleNestedQuotSchemePoint(rootMat, n00)
Tspace = tangentSpace quotPoint





Q1 = zeroo
Q2 = P
Q3 = P3
Q4 = P4++P2

phi21 = map(Q1, Q2, matrix{{1}}**R)
phi31 = map(Q1, Q3, matrix{{1}}**R)
phi42 = map(Q2, Q4, matrix{{-1,1}}**R)
phi43 = map(Q3, Q4, matrix{{-1, x^2+x}}**R)
assert (phi31*phi43 == phi21*phi42)
-- q4 = map(Q4, F, matrix{{x,-1},{1,0}}**R)
q4 = map(Q4, F, matrix{{1,0},{0,1}}**R)
q2 = phi42*q4
q3 = phi43*q4
q1 = phi21*q2

K1 = ker q1
K2 = ker q2
K3 = ker q3
K4 = ker q4

C1 = generators K1
C2 = generators K2
C3 = generators K3
C4 = generators K4


(D1,P1,QC1) = smithNormalForm C1
(D2,P2,QC2) = smithNormalForm C2
(D3,P3,QC3) = smithNormalForm C3
(D4,P4,QC4) = smithNormalForm C4

A1 = C1*QC1
A2 = C2*QC2
A3 = C3*QC3
A4 = C4*QC4
