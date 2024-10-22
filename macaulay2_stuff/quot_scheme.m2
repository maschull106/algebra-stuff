needs "hilbert_scheme.m2"

QuotSchemePoint = new Type of HashTable
NestedQuotSchemePoint = new Type of HashTable
QuotSchemeTangentSpace = new Type of TangentSpace
NestedQuotSchemeTangentSpace = new Type of TangentSpace


quotSchemePoint = method(TypicalValue => QuotSchemePoint)

quotSchemePoint Matrix := f -> (
    if dim(target f) != 0 then error "the quotient must be zero dimensional";

    return new QuotSchemePoint from {
        base => ring source f,
        constructionSheaf => source f,
        quotientSheaf => target f,
        quotientMap => f,
        kernelSheaf => kernel f,
        len => degree target f
    };
)


tangentSpace QuotSchemePoint := p -> (
    return new QuotSchemeTangentSpace from {
        originPoint => p,
        space => Hom(p.kernelSheaf, p.quotientSheaf)
    }
)


nestedQuotSchemePoint = method(TypicalValue => NestedQuotSchemePoint)

nestedQuotSchemePoint Matrix := f -> (
    if dim(target f) != 0 then error "the quotient must be zero dimensional";

    return new NestedQuotSchemePoint from {
        base => ring source f,
        constructionSheaf => source f,
        quotientSheaf => target f,
        quotientMap => f,
        kernelSheaf => kernel f,
        len => degree target f
    };
)


zeroList = n -> for i from 1 to n list 0;
zeroMutableList = n -> new MutableList from zeroList(n);

isZero = v -> v == 0_(module v);

mSize = m -> length flatten entries m;

-- toVector = v -> vector flatten entries v; -- (transpose matrix {flatten entries v})_0;
toVector = v -> vector flatten entries transpose matrix entries v;

-- reduceToBasis (Matrix, Matrix) := opts -> (v, Basis) -> (
--     v = (transpose matrix {flatten entries phi})_0;
--     return reduceToBasis(v, Basis, coords=>opts.coords);
-- )

leadCoefRatio = (p1, p2) -> leadCoefficient(p1) / leadCoefficient(p2)

reduceOnIndex = (v, Basis, i, coords) -> (
    while v_i != 0 do (
        for j from 0 to numcols(Basis)-1 do 
            if v_i == 0 then return v
            else if Basis_(i,j) == 0 then continue
            else if leadMonomial(v_i) == leadMonomial(Basis_(i,j)) then (
                c = leadCoefRatio(v_i, Basis_(i,j)); coords#j = coords#j + c; v = v - c*Basis_j;
            )
    );
    return v;
)

reduceToBasis = method(TypicalValue => List)
reduceToBasis (Matrix, Matrix) := (v, Basis) -> (
    originalV = v;
    v = toVector(v);
    coords = zeroMutableList(numcols Basis);
    makeReturn = x -> new List from coords;
    if isZero(v) then return makeReturn();
    
    for i from 0 to mSize(v)-1 do 
        if isZero(v) then return makeReturn()
        else v = reduceOnIndex(v, Basis, i, coords);

    if not isZero(v) then (print originalV; print v; error "Reducing to basis wrong!!!";);
    return makeReturn();
)




inclusionMap = (Target, Source) -> map(Target, Source, gens Source)


asMorphism = (v, Target, Source) -> (
    m = numgens Source;
    n = numgens Target;
    M = for i from 0 to n-1 list for j from 0 to m-1 list v_(j*m+i);
    return map(Target, Source, M);
)
