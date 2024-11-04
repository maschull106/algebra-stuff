needs "protect.m2"
needs "quot_scheme.m2"

DoubleNestedQuotSchemePoint = new Type of HashTable
DoubleNestedQuotSchemeTangentSpace = new Type of NestedQuotSchemeTangentSpace

NestingCell = new Type of HashTable
CellInfo = new Type of HashTable

GraphNode = new Type of HashTable
NodeInfo = new Type of HashTable


doubleNestedQuotSchemePoint = method(TypicalValue => DoubleNestedQuotSchemePoint)


cellInfo = method(TypicalValue => CellInfo)

cellInfo (NestingCell, QuotNesting) := (cell, nest) -> new CellInfo from {Cell => cell, Nesting => nest}

nestingCell = method(TypicalValue => NestingCell, Options => {CellInfoRight=>null, CellInfoDown=>null})

nestingCell (Module, Matrix) := opts -> (Q, p) -> (
    -- Q <<--- QRight
    -- ^
    -- |
    -- QDown

    -- p : F -->> Q

    -- TODO: sanity checks
    return new NestingCell from {
        TargetModule => Q,
        QuotientMap => p,
        HomBasis => basis Hom(ker p, Q),
        CellRight => if opts.CellInfoRight === null then null else opts.CellInfoRight.Cell,
        NestingRight => if opts.CellInfoRight === null then null else opts.CellInfoRight.Nesting,
        CellDown => if opts.CellInfoDown === null then null else opts.CellInfoDown.Cell,
        NestingDown => if opts.CellInfoDown === null then null else opts.CellInfoDown.Nesting
    }
)



nodeInfo = method(TypicalValue => NodeInfo)

nodeInfo (GraphNode, Matrix) := (n, f) -> new NodeInfo from {Node=>n, MapFromNode=>f}

graphNode = method(TypicalValue => GraphNode, Options => {Right=>null, Down=>null})


-- graphNode = {Right=>null, Down=>null, Quot=>null} >> opts -> () -> (
graphNode Matrix := opts -> quotMap -> (
    return new GraphNode from {
        NodeRight => if opts.Right === null then null else opts.Right.Node,
        MapFromRight => if opts.Right === null then null else opts.Right.MapFromNode,
        NodeDown => if opts.Down === null then null else opts.Down.Node,
        MapFromDown => if opts.Down === null then null else opts.Down.MapFromNode,
        QuotMap => quotMap
    }
)


constructNestingCell = method(TypicalValue => NestingCell)

constructNestingCell GraphNode := node -> (
    if node === null then return null;

    -- TODO: get rid of repeating code
    if not instance(cellsRight, List) then cellsRight = {};
    cellRight = if node.NodeRight === null then null else constructNestingCell(node.NodeRight);
    cellsRight = append(cellsRight, cellRight);
    cellDown = if node.NodeDown === null then null else constructNestingCell(node.NodeDown);
    cellRight = cellsRight_-1;
    cellsRight = drop(cellsRight, -1);

    if not(node.MapFromRight === null) then (
        morph = node.MapFromRight;
        T2 = source morph;
        T1 = target morph;
        morphFT2 = node.NodeRight.QuotMap;
        nestingRight = quotNesting(T2, T1, morphFT2, morph);
    ) else nestingRight = null;
    infoRight = if cellRight === null then null else cellInfo(cellRight, nestingRight);

    if not(node.MapFromDown === null) then (
        morph = node.MapFromDown;
        T2 = source morph;
        T1 = target morph;
        morphFT2 = node.NodeDown.QuotMap;
        nestingDown = quotNesting(T2, T1, morphFT2, morph);
    ) else nestingDown = null;
    infoDown = if cellDown === null then null else cellInfo(cellDown, nestingDown);

    cell = nestingCell(target node.QuotMap, node.QuotMap, CellInfoRight=>infoRight, CellInfoDown=>infoDown);
    return cell;
)

doubleNestedQuotSchemePoint GraphNode := node -> (
    F = source node.QuotMap;
    nestBase = constructNestingCell(node);
    return doubleNestedQuotSchemePoint(F, nestBase);
)


-- TODO: better name for this methods
constructNestingArrays = method(TypicalValue => List, Options => {OnRows => true})

constructNestingArrays NestingCell := opts -> nestCell -> (
    nextArray = nestBase;
    arrays = while not(nextArray === null) list (
        nextCell = nextArray;
        nextArray = nextArray#(if opts.OnRows then CellDown else CellRight);
        while not(nextCell === null) list (
            nesting = nextCell#(if opts.OnRows then NestingRight else NestingDown);
            nextCell = nextCell#(if opts.OnRows then CellRight else CellDown);
            if nesting === null then break;
            nesting
        )
    );
    return arrays;
)

doubleNestedQuotSchemePoint (Module, NestingCell) := (F, nestBase) -> (
    rows = constructNestingArrays(nestBase, OnRows=>true);
    cols = constructNestingArrays(nestBase, OnRows=>false);
    l = homBasesLengths(nestBase);

    return new DoubleNestedQuotSchemePoint from {
        ConstructionSheaf => F,
        NestingBase => nestBase,
        NestedRows => rows,
        NestedCols => cols,
        BasisLengthsCumulative => l_0,
        CellFlatIndices => l_1
    }
)


iterateArrays = method(Options => {OnRows=>true, PreOp=>()->null, LoopOp=>(cell, row, col)->null, PostOp=>()->null})
iterateRows = method(Options => {PreOp=>()->null, LoopOp=>(cell, row, col)->null, PostOp=>()->null})
iterateCols = method(Options => {PreOp=>()->null, LoopOp=>(cell, row, col)->null, PostOp=>()->null})

iterateRows NestingCell := opts -> nestBase -> return iterateArrays(nestBase, OnRows=>true, PreOp=>opts.PreOp, LoopOp=>opts.LoopOp, PostOp=>opts.PostOp)
iterateCols NestingCell := opts -> nestBase -> return iterateArrays(nestBase, OnRows=>false, PreOp=>opts.PreOp, LoopOp=>opts.LoopOp, PostOp=>opts.PostOp)

iterateArrays NestingCell := opts -> nestBase -> (
    opts.PreOp();

    nextRow = nestBase;
    row = 0;
    while not(nextRow === null) do (
        nextCell = nextRow;
        col = 0;
        while not(nextCell === null) do (
            if opts.OnRows then opts.LoopOp(nextCell, row, col) else opts.LoopOp(nextCell, col, row);
            col = col+1;
            nextCell = nextCell#(if opts.OnRows then CellRight else CellDown);
        );
        row = row+1;
        nextRow = nextRow#(if opts.OnRows then CellDown else CellRight);
    );
    return opts.PostOp();
)


homBasesLengths = method(TypicalValue => List)

homBasesLengths NestingCell := nestBase -> (
    lengths = {};
    lengthsCumulative = {};
    flatIndices = new MutableHashTable from {};

    iterateRows(
        nestBase, 
        LoopOp => (cell, row, col) -> (
                flatIndices#(row, col) = length(lengths);
                lengthsCumulative = append(lengthsCumulative, sum(lengths));
                lengths = append(lengths, basisLength(nextCell.HomBasis));
            )
        );
    lengthsCumulative = append(lengthsCumulative, sum(lengths));
    return {lengthsCumulative, flatIndices};
)


doubleNestingConstraint = method(TypicalValue => Matrix)

doubleNestingConstraint (DoubleNestedQuotSchemePoint, QuotNesting, ZZ, ZZ) := (p, nesting, sourceIndex, targetIndex) -> (
    totalM = p.BasisLengthsCumulative_-1;
    constr = nestedQuotTangentConstraints(nesting);
    n = numrows constr.ConstraintOnSource;
    C = zeroMutableMatrix(n, totalM);
    matrixWriteSlice(C, 0, p.BasisLengthsCumulative_sourceIndex, constr.ConstraintOnSource);
    matrixWriteSlice(C, 0, p.BasisLengthsCumulative_targetIndex, constr.ConstraintOnTarget);
    return new Matrix from C;
)


tangentSpace DoubleNestedQuotSchemePoint := p -> (
    totalM = p.BasisLengthsCumulative_-1;
    totalC = zeroMatrix(0, totalM);

    nextRow = p.NestingBase;
    row = 0;
    while not(nextRow === null) do (
        nextCell = nextRow;
        col = 0;
        while not(nextCell === null) do (
            -- nesting to the right
            if not(nextCell.CellRight === null) then (
                targetIndex = p.CellFlatIndices#(row, col);
                sourceIndex = p.CellFlatIndices#(row, col+1);
                C = doubleNestingConstraint(p, nextCell.NestingRight, sourceIndex, targetIndex);
                totalC = totalC || C;
            );
            if not(nextCell.CellDown === null) then (
                targetIndex = p.CellFlatIndices#(row, col);
                sourceIndex = p.CellFlatIndices#(row+1, col);
                C = doubleNestingConstraint(p, nextCell.NestingDown, sourceIndex, targetIndex);
                totalC = totalC || C;
            );

            col = col+1;
            nextCell = nextCell.CellRight;
        );
        row = row+1;
        nextRow = nextRow.CellDown;
    );

    return new DoubleNestedQuotSchemeTangentSpace from {
        OriginPoint => p,
        Space => ker totalC
    };
)



matrixWriteSlice = (M, rowStart, colStart, Mslice) -> (
    for i from 0 to numrows(Mslice)-1 do for j from 0 to numcols(Mslice)-1 do M_(rowStart+i,colStart+j) = Mslice_(i,j);
)



nestedQuotSchemePoint2 = method(TypicalValue => DoubleNestedQuotSchemePoint)
nestedQuotSchemePoint2 List := morphisms -> (
    -- TODO: surjectivity check
    p = morphisms_0;
    F = source p;
    node = graphNode(p);
    nestings = for i from 1 to length(morphisms)-1 list (
        p = morphisms_i * p;
        node = graphNode(p, Right=>nodeInfo(node, morphisms_i));
    );
    return doubleNestedQuotSchemePoint(node);
)


nestedHilbSchemePoint2 = method(TypicalValue => DoubleNestedQuotSchemePoint)

nestedHilbSchemePoint2 List := modules -> (
    R = ring modules_0;
    r = rank modules_0;
    idMat = idMatrix(r)**R;
    morphisms = for i from 0 to length(modules)-2 list (
        map(modules_(i+1), modules_i, idMat)
    );
    return nestedQuotSchemePoint2(morphisms);
)






makeNode = method(TypicalValue => GraphNode)
makeNode (ZZ, ZZ, List, MutableHashTable) := (row, col, data, memory) -> (
    memoryFetch = (r, c) -> (if not member((r, c), keys memory) then memory#(r, c) = makeNode(r, c, data, memory) else print "fetching from memory"; return memory#(r, c););
    getModule = (r, c) -> (if r%2 != 0 or c%2 != 0 then error("invalid indices"); cell = data_r_c; if instance(cell, Module) then return cell; if instance(cell_0, Module) then return cell_0; error("couldn't get module"));

    -- TODO: get rid of repeating code
    if length data_row > col+1 then (
        n = memoryFetch(row, col+2);
        trgt = data_row_col;
        src = getModule(row, col+2);
        R = ring trgt;
        f = map(trgt, src, data_row_(col+1)**R);
        nodeRight = nodeInfo(n, f);
        qRight = f * n.QuotMap;
        -- temporary save of right info
        memory#(row, col) = graphNode(qRight, Right=>nodeRight);
    ) else (nodeRight = null; qRight = null;);
    if length data_col > row+1 then (
        n = memoryFetch(row+2, col);
        trgt = data_row_col;
        src = getModule(row+2, col);
        R = ring trgt;
        mat = data_(row+1)_(col//2);
        f = map(trgt, src, mat**R);
        down = nodeInfo(n, f);
        qDown = f * n.QuotMap;
    ) else (down = null; qDown = null;);

    hasRightBranch = member((row, col), keys memory);
    if not hasRightBranch and down === null then (
        cell = data_row_col;
        Q = cell_0; q = cell_1;
        return graphNode(q);
    );

    q = qDown;
    if hasRightBranch then (
        n = memory#(row, col); 
        qRight = n.QuotMap;
        q = qRight;
        nodeRight = nodeInfo(n.NodeRight, n.MapFromRight);
        if qLeft =!= null and qRight != qLeft then error("noncommutative");
    ) else (nodeRight = null);

    return graphNode(q, Right=>nodeRight, Down=>down);
)

doubleNestedQuotSchemePoint (Module, List) := (F, data) -> (
    node = makeNode(0, 0, data, new MutableHashTable from {});
    return doubleNestedQuotSchemePoint(node);
)
