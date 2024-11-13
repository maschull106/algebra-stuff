needs "double_nested_quot_scheme.m2"


R = QQ[x]
F = R^2
one = R^1/(ideal(1)**R)
P = R^1/ideal(x)
Q = R^1/ideal(x-1)
P2 = R^1/ideal(x^2)
PQ = R^1/ideal(x*(x-1))


T1 = tangentSpace doubleNestedQuotSchemePoint(
    F,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1,0}}},
        {P,    {{1,0}},    (P++P, idMatrix(2))}
    }
)


T2 = tangentSpace doubleNestedQuotSchemePoint(
    F,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1,0}}},
        {P,    {{0,1}},    (P++P, idMatrix(2))}
    }
)


T3 = tangentSpace doubleNestedQuotSchemePoint(
    F,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1}}},
        {P,    {{1}},    (P2, {{1,1}})}
    }
)


T4 = tangentSpace doubleNestedQuotSchemePoint(
    F,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1}}},
        {P,    {{1}},    (P2, {{1,0}})}
    }
)


T5 = tangentSpace doubleNestedQuotSchemePoint(
    F,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1}}},
        {P,    {{1}},    (PQ, {{1,0}})}
    }
)


T6 = tangentSpace doubleNestedQuotSchemePoint(
    F,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1}}},
        {P,    {{1}},    (PQ, {{x, x-1}})}
    }
)


T7 = tangentSpace doubleNestedQuotSchemePoint(
    F,
    {
        {one,    {{1}},    P},
        {{{1}},           {{1}}},
        {P,    {{1}},    (P2, {{1, 0}})}
    }
)

print (dim T7)