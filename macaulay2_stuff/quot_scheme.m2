needs "hilbert_scheme.m2"
needs "protect.m2"


QuotSchemePoint = new Type of HashTable
NestedQuotSchemePoint = new Type of HashTable
QuotSchemeTangentSpace = new Type of TangentSpace
NestedQuotSchemeTangentSpace = new Type of TangentSpace

QuotNesting = new Type of HashTable
SimpleNestedQuotSchemePoint = new Type of HashTable
SimpleNestedQuotSchemeTangentSpace = new Type of TangentSpace


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

-- isZero = v -> v == 0_(module v);    -- doesn't always work correctly for some reason
isZero = v -> v == (v-v);

mSize = m -> length flatten entries m;

basisLength = B -> numcols B

-- toVector = v -> vector flatten entries v; -- (transpose matrix {flatten entries v})_0;
toVector = v -> vector flatten entries transpose matrix entries v;

-- reduceToBasis (Matrix, Matrix) := opts -> (v, Basis) -> (
--     v = (transpose matrix {flatten entries phi})_0;
--     return reduceToBasis(v, Basis, coords=>opts.coords);
-- )

leadCoefRatio = (p1, p2) -> leadCoefficient(p1) / leadCoefficient(p2)

reduceOnIndex = (v, Basis, i, coords) -> (
    while v_i != 0 do (
        for j from 0 to basisLength(Basis)-1 do 
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
    coords = zeroMutableList(basisLength Basis);
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


zeroMutableMatrix = (n, m) -> new MutableMatrix from matrix(for i from 1 to n list zeroList(m))




quotNesting = method(TypicalValue => QuotNesting)

quotNesting (Module, Module, Matrix, Matrix) := (T2, T1, morphFT2, morphT2T1) -> (
    if dim(T2) != 0 or dim(T1) != 0 then error "the quotients must be zero dimensional";

    quot2 = morphFT2;
    quot1 = morphT2T1 * morphFT2;
    K2 = ker quot2;
    K1 = ker quot1;

    -- psi: K2 -> K1
    -- phi: T2 -> T1
    phi = morphT2T1;
    psi = inclusionMap(K1, K2);

    return new QuotNesting from {
        Module2 => T2,
        Module1 => T1,
        Quotient2 => quot2,
        Quotient1 => quot1,
        Phi => phi,
        Psi => psi,
        Kernel2 => K2,
        Kernel1 => K1,
        Length => {degree T1, degree T2},
        HomBasis1 => basis Hom(K1, T1),
        HomBasis2 => basis Hom(K2, T2),
        HomCrossBasis => basis Hom(K2, T1)
    };
)


nestedQuotTangentConstraints = method(TypicalValue => Matrix)

nestedQuotTangentConstraints QuotNesting := nest -> (
    m1 = basisLength nest.HomBasis1;
    m2 = basisLength nest.HomBasis2;
    m = m1 + m2;
    n = basisLength nest.HomCrossBasis;

    -- C = zeroMutableMatrix(n, m);

    psiCompose = v -> (
        fK1T1 = asMorphism(v, nest.Module1, nest.Kernel1);
        fK2T1 = fK1T1 * nest.Psi;
        coords = reduceToBasis(fK2T1, nest.HomCrossBasis);
        return coords;
    );
    phiCompose = v -> (
        fK2T2 = asMorphism(v, nest.Module2, nest.Kernel2);
        fK2T1 = - nest.Phi * fK2T2;
        coords = reduceToBasis(fK2T1, nest.HomCrossBasis);
        return coords;
    );

    C1 = for i from 0 to m1-1 list psiCompose(nest.HomBasis1_i);

    C2 = for i from 0 to m2-1 list phiCompose(nest.HomBasis2_i);

    -- C1 = matrix C1;
    -- C1 = transpose C1;
    -- C2 = matrix C2;
    -- C2 = transpose C2;
    -- return C1 | C2;
    return transpose(matrix(C1)) | transpose(matrix(C2));
)


simpleNestedQuotSchemePoint = method(TypicalValue => SimpleNestedQuotSchemePoint)

simpleNestedQuotSchemePoint (Module, QuotNesting) := (F, nesting) -> (
    return new SimpleNestedQuotSchemePoint from {
        ConstructionSheaf => F,
        Nesting => nesting
    };
)


tangentSpace SimpleNestedQuotSchemePoint := p -> (
    C = nestedQuotTangentConstraints(p.Nesting);

    return new SimpleNestedQuotSchemeTangentSpace from {
        OriginPoint => p,
        Space => ker C
    };
)


dim SimpleNestedQuotSchemeTangentSpace := T -> rank T.Space


