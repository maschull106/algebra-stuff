needs "quot_scheme.m2"

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




-- nesting = quotNesting(T2, T1, p2, f)
-- quotPoint = nestedQuotSchemePoint(F, {nesting})
-- T = tangentSpace quotPoint





-- p1 = f*p2
-- K1 = ker p1
-- K2 = ker p2
-- -- p = quotSchemePoint(f);

-- -- T = tangentSpace p;


-- phi = f;
-- psi = inclusionMap(K1, K2);

-- H1 = Hom(K1, T1);
-- H2 = Hom(K2, T2);
-- H3 = Hom(K2, T1);

-- basis1 = basis H1;
-- basis2 = basis H2;
-- basis3 = basis H3;

-- L = constraints(K1, K2, T1, T2, psi, phi, basis1, basis2, basis3);

-- -- phi = map(T1, K2, matrix{{x+1, 2-x}})








-- Q0 <--- Q1
--  ^
--  |
-- Q2


R = QQ[x]
F = R^1
I0 = ideal(x); I1 = ideal(x^2); I2 = ideal(x*(x-1))
T0 = F/I0; T1 = F/I1; T2 = F/I2
idMat = idMatrix(1)**R;
q0 = map(T0, F, idMat); q1 = map(T1, F, idMat); q2 = map(T2, F, idMat)
f20 = map(T0, T2, idMat)
f10 = map(T0, T1, idMat)
node2 = graphNode(q2)
node1 = graphNode(q1)
node0 = graphNode(
    q0, 
    Right=>nodeInfo(node1, f10),
    Down=>nodeInfo(node2, f20)
)
p = doubleNestedQuotSchemePoint(node0)
T = tangentSpace p
