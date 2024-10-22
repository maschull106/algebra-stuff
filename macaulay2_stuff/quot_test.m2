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
M = matrix{{1,1},{0,0}}
M = M**R -- force M to be over R
p2 = map(T2, F, M)
N = matrix{{0, 1},{1,0}}
N = N**R
f = map(T1, T2, N)
p1 = f*p2
K1 = ker p1
K2 = ker p2
-- p = quotSchemePoint(f);

-- T = tangentSpace p;


phi = f;
psi = inclusionMap(K1, K2);

H1 = Hom(K1, T1);
H2 = Hom(K2, T2);
H3 = Hom(K2, T1);

base1 = basis H1;
base2 = basis H2;
base3 = basis H3;

-- phi = map(T1, K2, matrix{{x+1, 2-x}})
