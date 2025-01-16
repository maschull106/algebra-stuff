needs "double_nested_quot_scheme.m2"

R = QQ[x];

-- construction sheaf for the Quot scheme
F = R^2;

-- another sheaf Q
I = ideal(x^2);
Iprime = ideal(x-1);
P = R^1/I;
Pprime = R^1/Iprime;
Q = P ++ Pprime;
T2 = Q
T1 = Q

-- let's define a point [F --> Q] of the Quot scheme
M = matrix{{1,0},{1,0}}
M = M**R -- force M to be over R
p2 = map(T2, F, M)
N = matrix{{0, 1},{1,0}}
N = N**R
f = map(T1, T2, N)




F = R^1;
Q = R^1/ideal(x-1);
P = R^1/ideal(x);
T2 = P ++ Q;
T1 = P;
p2 = map(T2, F, matrix{{1},{1}}**R);
f = map(T1, T2, matrix{{1, 0}}**R);




R=QQ[x,y,z]; I=ideal(x^2, x*y^2, x*y*z, x*z^2, y^2*z^2, y*z^3, z^4, y^3-x*z)
F = R^1
O = R^1/I
f = map(O, F, matrix{{1}}**R)

-- quotPoint = nestedQuotSchemePoint({f})
-- T = tangentSpace quotPoint




-- Cheah's nested Hilbert schemes singularity examples

R = QQ[x,y]
I1 = ideal(x, y); I2 = ideal(x^2, y); I3 = ideal(x^2, x*y, y^2)
F = R^1
T1 = F/I1; T2 = F/I2; T3 = F/I3
quotPoint1 = nestedHilbSchemePoint2({F, T3, T2, T1})
Tspace1 = tangentSpace quotPoint1


quotPoint2 = nestedHilbSchemePoint2({F, T3, T1})
Tspace2 = tangentSpace quotPoint2


R = QQ[x,y, z]
I1 = ideal(x^2, y^2, z^2, x*y, x*z, y*z); I2 = ideal(x^2, y^2, z^2)
F = R^1
T1 = F/I1; T2 = F/I2
quotPoint3 = nestedHilbSchemePoint2({F, T2, T1})
Tspace3 = tangentSpace quotPoint3




-- Singularity examples for nested Quot schemes

R = QQ[x, y]
F = R^2
P = R^1/ideal(x, y)
T2 = P ++ P
T1 = P
p2 = map(T2, F, matrix{{1,0},{0,1}}**R)
f = map(T1, T2, matrix{{1, 1}}**R)
quotPoint4 = nestedQuotSchemePoint2({p2, f})
Tspace4 = tangentSpace quotPoint4








-- Q0 <--- Q1
--  ^       ^
--  |       |
-- Q2 <--- Q3


R = QQ[x]
F = R^1
I0 = ideal(1)**R; I1 = ideal(x); I2 = ideal(x); I3 = ideal(x^2)
T0 = F/I0; T1 = F/I1; T2 = F/I2; T3 = F/I3
idMat = idMatrix(1)**R;
q0 = map(T0, F, idMat); q1 = map(T1, F, idMat); q2 = map(T2, F, idMat); q3 = map(T3, F, idMat)
f20 = map(T0, T2, idMat)
f10 = map(T0, T1, idMat)
f32 = map(T2, T3, idMat)
f31 = map(T1, T3, idMat)
node3 = quotNode(q3)
node2 = quotNode(q2, Right=>nodeInfo(node3, f32))
node1 = quotNode(q1, Down=>nodeInfo(node3, f31))
node0 = quotNode(
    q0, 
    Right=>nodeInfo(node1, f10),
    Down=>nodeInfo(node2, f20)
)
quotPoint5 = doubleNestedQuotSchemePoint(node0)
Tspace5 = tangentSpace quotPoint5







R = QQ[x]
F = R^2
I0 = ideal(1)**R; I1 = ideal(x); I2 = ideal(x); I3 = ideal(x^2)
T0 = R^1/I0; T1 = R^1/I1; T2 = R^1/I2; T3 = T1++T1 --T3 = R^1/I3
idMat = idMatrix(1)**R;
pMat = matrix{{1, 1}}**R
q0 = map(T0, F, pMat); q1 = map(T1, F, pMat); q2 = map(T2, F, pMat); q3 = map(T3, F, matrix{{1,0},{0,1}}**R)
f20 = map(T0, T2, idMat)
f10 = map(T0, T1, idMat)
f32 = map(T2, T3, pMat)
f31 = map(T1, T3, pMat)
node3 = quotNode(q3)
node2 = quotNode(q2, Right=>nodeInfo(node3, f32))
node1 = quotNode(q1, Down=>nodeInfo(node3, f31))
node0 = quotNode(
    q0, 
    Right=>nodeInfo(node1, f10),
    Down=>nodeInfo(node2, f20)
)
quotPoint6 = doubleNestedQuotSchemePoint(node0)
Tspace6 = tangentSpace quotPoint6


-- {
--     {T0,    f10,    T1,     f41,    T4},
--     {f20,           f31},
--     {T2,    f32,    T3},
--     {f25},
--     {(T5, f)}
-- }



R = QQ[x]
F = R^2
one = R^1/(ideal(1)**R); P = R^1/ideal(x); Q = R^1/ideal(x-1); P2 = R^1/ideal(x^2);
PQ = R^1/ideal(x*(x-1));
q3 = map(P++Q, F, matrix{{1,0},{0,1}}**R)
-- q3 = map(PQ, F, matrix{{1,1}}**R)



quotPoint7 = doubleNestedQuotSchemePoint(
    F,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1,0}}},
        {P,    {{1,0}},    (P++P, idMatrix(2))}
    }
)
Tspace7 = tangentSpace quotPoint7



quotPoint8 = doubleNestedQuotSchemePoint(
    F,
    {
        {   one,    {{1}},      P           },
        {   {{1}},              {{1}}       },
        {   P,      {{1}},      (P2, {{1, 1}}) }
    }
)
Tspace8 = tangentSpace quotPoint8




p10 = nestedHilbSchemePoint2({R^1, P, one})
T10 = tangentSpace p10



quotPoint9 = doubleNestedQuotSchemePoint(
    R^1,
    {
        {   one,    {{1}},      (P, {{1}} )          }
    }
)
Tspace9 = tangentSpace quotPoint9
