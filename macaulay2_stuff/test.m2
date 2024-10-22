needs "hilbert_scheme.m2"

R = QQ[x, y];
I = ideal(x^2, y^2);

p = hilbertSchemePoint I;
