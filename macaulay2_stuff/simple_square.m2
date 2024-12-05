R = QQ[x]
F = R^2
one = R^1/(ideal(1)**R)
P = R^1/ideal(x)
Q = R^1/ideal(x-1)
P2 = R^1/ideal(x^2)
PQ = R^1/ideal(x*(x-1))




tangentDim = p -> dim(tangentSpace p)


F1 = (p,q) -> (
    return doubleNestedQuotSchemePoint(
        F,
        {
            {one,    {{1}},    P},
            {{{1}},           {{1}}},
            {P,    {{1}},    (P2, {{p,q}})}
        }
    );
)


F2 = A -> (
    return doubleNestedQuotSchemePoint(
        F,
        {
            {one,    {{1}},    P},
            {{{1}},           {{1,0}}},
            {P,    {{1,0}},    (P++P, matrix A)}
        }
    );
)


F3 = A -> (
    return doubleNestedQuotSchemePoint(
        F,
        {
            {one,    {{1}},    P},
            {{{1}},           {{1,0}}},
            {P,    {{0,1}},    (P++P, matrix A)}
        }
    );
)
