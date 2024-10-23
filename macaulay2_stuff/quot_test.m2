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


nesting = quotNesting(T2, T1, p2, f)
quotPoint = simpleNestedQuotSchemePoint(F, nesting)
T = tangentSpace quotPoint





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
