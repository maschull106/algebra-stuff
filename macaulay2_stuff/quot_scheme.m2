needs "hilbert_scheme.m2"
needs "protect.m2"


QuotSchemePoint = new Type of HashTable
NestedQuotSchemePoint = new Type of HashTable
QuotSchemeTangentSpace = new Type of TangentSpace
NestedQuotSchemeTangentSpace = new Type of TangentSpace

QuotNesting = new Type of HashTable
SimpleNestedQuotSchemePoint = new Type of NestedQuotSchemePoint
SimpleNestedQuotSchemeTangentSpace = new Type of NestedQuotSchemeTangentSpace


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

toVector = v -> vector flatten entries v; -- (transpose matrix {flatten entries v})_0;
-- toVector = v -> vector flatten entries transpose matrix entries v;

-- reduceToBasis (Matrix, Matrix) := opts -> (v, Basis) -> (
--     v = (transpose matrix {flatten entries phi})_0;
--     return reduceToBasis(v, Basis, coords=>opts.coords);
-- )

leadCoefRatio = (p1, p2) -> leadCoefficient(p1) / leadCoefficient(p2)

reduceOnIndex = (v, Basis, i, coords) -> (
    while v_i != 0 do (
        for j from 0 to basisLength(Basis)-1 do 
            if v_i == 0 then return v
            else if Basis_j_i == 0 then continue
            else if leadMonomial(v_i) == leadMonomial(Basis_j_i) then (
                c = leadCoefRatio(v_i, Basis_j_i); coords#j = coords#j + c; v = v - c*toVector(Basis_j);
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



-- TODO: investigate this function more
-- inclusionMap = (Target, Source) -> try map(Target, Source, gens Source) else inducedMap(Target, Source)
inclusionMap = (Target, Source) -> try inducedMap(Target, Source) else (print "warning, classical induced map didn't work!!"; map(Target, Source, gens Source))


asMorphism = (v, Target, Source) -> (
    m = numgens Source;
    n = numgens Target;
    M = for i from 0 to n-1 list for j from 0 to m-1 list v_(j*n+i);
    return map(Target, Source, M);
)


zeroMatrix = (n, m) -> if n == 0 then if m == 0 then matrix{{}}**Q else transpose zeroMatrix(m, n) else matrix(for i from 1 to n list zeroList(m))**QQ
zeroMutableMatrix = (n, m) -> (new MutableMatrix from zeroMatrix(n, m))

idMatrix = n -> matrix(for i from 1 to n list for j from 1 to n list if i==j then 1 else 0)



quotNesting = method(TypicalValue => QuotNesting)

quotNesting (Module, Module, Matrix, Matrix) := (T2, T1, morphFT2, morphT2T1) -> (
    if dim(T2) > 0 or dim(T1) > 0 then error "the quotients must be zero dimensional";

    quot2 = morphFT2;
    quot1 = morphT2T1 * morphFT2;
    K2 = ker quot2;
    K1 = ker quot1;

    -- psi: K2 -> K1
    -- phi: T2 -> T1
    phi = morphT2T1;
    psi = inclusionMap(K1, K2);

    return new QuotNesting from {
        ModuleSource => T2,
        ModuleTarget => T1,
        QuotientSource => quot2,
        QuotientTarget => quot1,
        Phi => phi,
        Psi => psi,
        KernelSource => K2,
        KernelTarget => K1,
        Length => {degree T1, degree T2},
        HomBasisTarget => basis Hom(K1, T1),
        HomBasisSource => basis Hom(K2, T2),
        HomCrossBasis => basis Hom(K2, T1)
    };
)


ConstraintsCouple = new Type of HashTable

nestedQuotTangentConstraints = method(TypicalValue => ConstraintsCouple)

nestedQuotTangentConstraints QuotNesting := nest -> (
    m1 = basisLength nest.HomBasisTarget;
    m2 = basisLength nest.HomBasisSource;
    m = m1 + m2;
    n = basisLength nest.HomCrossBasis;

    -- C = zeroMutableMatrix(n, m);

    -- H = Hom(nest.KernelSource, nest.ModuleTarget);
    psiCompose = v -> (
        fK1T1 = asMorphism(v, nest.ModuleTarget, nest.KernelTarget);
        fK2T1 = fK1T1 * nest.Psi;
        coords = reduceToBasis(fK2T1, nest.HomCrossBasis);
        return coords;
    );
    phiCompose = v -> (
        fK2T2 = asMorphism(v, nest.ModuleSource, nest.KernelSource);
        fK2T1 = - nest.Phi * fK2T2;
        coords = reduceToBasis(fK2T1, nest.HomCrossBasis);
        return coords;
    );

    C1 = for i from 0 to m1-1 list psiCompose(nest.HomBasisTarget_i);

    C2 = for i from 0 to m2-1 list phiCompose(nest.HomBasisSource_i);

    if length C1 == 0 then C1 = {{}};
    if length C2 == 0 then C2 = {{}};
    C1 = transpose(matrix(C1));
    C2 = transpose(matrix(C2));

    return new ConstraintsCouple from {ConstraintOnTarget => C1, ConstraintOnSource => C2};
    -- return transpose(matrix(C1)) | transpose(matrix(C2));
)


simpleNestedQuotSchemePoint = method(TypicalValue => SimpleNestedQuotSchemePoint)

simpleNestedQuotSchemePoint (Module, QuotNesting) := (F, nesting) -> (
    return new SimpleNestedQuotSchemePoint from {
        ConstructionSheaf => F,
        Nesting => nesting
    };
)


tangentSpace SimpleNestedQuotSchemePoint := p -> (
    constraints = nestedQuotTangentConstraints(p.Nesting);
    C = constraints.ConstraintOnTarget | constraints.ConstraintOnSource;

    return new SimpleNestedQuotSchemeTangentSpace from {
        OriginPoint => p,
        Space => ker C
    };
)


dim NestedQuotSchemeTangentSpace := T -> rank T.Space


nestedHilbSchemePoint = method(TypicalValue => NestedQuotSchemePoint)

nestedHilbSchemePoint List := modules -> (
    R = ring modules_0;
    r = rank modules_0;
    idMat = idMatrix(r)**R;
    morphisms = for i from 0 to length(modules)-2 list (
        map(modules_(i+1), modules_i, idMat)
    );
    return nestedQuotSchemePoint(morphisms);
)


nestedQuotSchemePoint List := morphisms -> (
    -- TODO: surjectivity check 
    F = source morphisms_0;
    morphFT2 = morphisms_0;
    nestings = for i from 1 to length(morphisms)-1 list (
        T2 = source morphisms_i;
        T1 = target morphisms_i;
        -- if i > 1 then (print target morphFT2; print source morphisms_(i-1); morphFT2 = morphisms_(i-1) * morphFT2;);
        if i > 1 then morphFT2 = morphisms_(i-1) * morphFT2;
        quotNesting(T2, T1, morphFT2, morphisms_i)
    );
    return nestedQuotSchemePoint(F, nestings);
)

nestedQuotSchemePoint (Module, List) := (F, nestings) -> (
    return new NestedQuotSchemePoint from {
        ConstructionSheaf => F,
        Nestings => nestings
    };
)


tangentSpace NestedQuotSchemePoint := p -> (
    widthSoFar = 0;
    for i from 0 to length(p.Nestings)-1 do (
        constr = nestedQuotTangentConstraints(p.Nestings_i);
        C = constr.ConstraintOnSource | constr.ConstraintOnTarget;
        w2 = numcols constr.ConstraintOnSource;
        w1 = numcols constr.ConstraintOnTarget;
        h = numrows C;
        if widthSoFar == 0 then (totalC = C; widthSoFar = w2; continue;);
        currentH = numrows totalC;
        totalC = (totalC | zeroMatrix(currentH, w1)) || (zeroMatrix(h, widthSoFar) | C);
        widthSoFar = widthSoFar + w2;
    );

    return new NestedQuotSchemeTangentSpace from {
        OriginPoint => p,
        Space => ker totalC
    };
)
