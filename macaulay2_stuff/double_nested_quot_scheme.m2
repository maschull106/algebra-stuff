needs "protect.m2"
needs "quot_scheme.m2"

DoubleNestedQuotSchemePoint = new Type of HashTable
DoubleNestedQuotSchemeTangentSpace = new Type of NestedQuotSchemeTangentSpace

NestingCell = new Type of HashTable
CellInfo = new Type of HashTable

GraphNode = new Type of HashTable
NodeInfo = new Type of HashTable

QuotKernelNode = new Type of MutableHashTable
KernelNodeInfo = new Type of HashTable


doubleNestedQuotSchemePoint = method(TypicalValue => DoubleNestedQuotSchemePoint)


cellInfo = method(TypicalValue => CellInfo)

cellInfo (NestingCell, QuotNesting) := (cell, nest) -> new CellInfo from {Cell => cell, Nesting => nest}

nestingCell = method(TypicalValue => NestingCell, Options => {CellInfoRight=>null, CellInfoDown=>null})

nestingCell (Module, Matrix) := opts -> (quotModule, quotMap) -> (
    -- Q <<--- QRight
    -- ^
    -- |
    -- QDown

    -- p : F -->> Q

    -- TODO: sanity checks
    return new NestingCell from {
        TargetModule => quotModule,
        QuotientMap => quotMap,
        HomBasis => basis Hom(ker quotMap, quotModule),
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






kernelNodeInfo = method(TypicalValue => KernelNodeInfo)

kernelNodeInfo (QuotKernelNode, Matrix) := (n, B) -> new KernelNodeInfo from {Node=>n, KernelMat=>B}

-- quotKernelNode = method(TypicalValue => QuotKernelNode, Options => {Right=>null, Down=>null})

quotKernelNode = {Right=>null, Down=>null} >> opts -> () -> (
    return new QuotKernelNode from {
        NodeRight => if opts.Right === null then null else opts.Right.Node,
        KernelMatRight => if opts.Right === null then null else opts.Right.KernelMat,
        NodeDown => if opts.Down === null then null else opts.Down.Node,
        KernelMatDown => if opts.Down === null then null else opts.Down.KernelMat,
        GraphNodeConversion => null
    }
)


isLeafNode = method(TypicalValue => Boolean)
isLeafNode QuotKernelNode := n -> (n.NodeRight === null and n.NodeDown === null)


kernelNodeToGraphNode = method(TypicalValue => GraphNode)

kernelNodeToGraphNode (Matrix, QuotKernelNode) := (rootKernelMat, n) -> (
    if n.GraphNodeConversion =!= null then return n.GraphNodeConversion;
    
    K := image rootKernelMat;
    F := target rootKernelMat;
    Q := F/K;
    q := inducedMap(Q, F);

    if isLeafNode(n) then (
        conversion := graphNode(q);
        n.GraphNodeConversion = conversion;
        return conversion;
    );

    infoRight := null;
    infoDown := null;
    if n.NodeRight =!= null then (
        nodeRight := kernelNodeToGraphNode(rootKernelMat*n.KernelMatRight, n.NodeRight);
        f := inducedMap(Q, target nodeRight.QuotMap);
        infoRight = nodeInfo(nodeRight, f);
    );
    if n.NodeDown =!= null then (
        nodeDown := kernelNodeToGraphNode(rootKernelMat*n.KernelMatDown, n.NodeDown);
        f := inducedMap(Q, target nodeDown.QuotMap);
        infoDown = nodeInfo(nodeDown, f);
    );

    conversion := graphNode(q, Right=>infoRight, Down=>infoDown);
    n.GraphNodeConversion = conversion;
    return conversion;
)


doubleNestedQuotSchemePoint (Matrix, QuotKernelNode) := (rootKernelMat, n) -> (
    gNode = kernelNodeToGraphNode(rootKernelMat, n);
    return doubleNestedQuotSchemePoint(gNode);
)





constructNestingCell = method(TypicalValue => NestingCell)

constructNestingCell GraphNode := node -> (
    if node === null then return null;

    cellRight := if node.NodeRight === null then null else constructNestingCell(node.NodeRight);
    cellDown := if node.NodeDown === null then null else constructNestingCell(node.NodeDown);

    nestingRight = null;
    if not(node.MapFromRight === null) then (
        morph := node.MapFromRight;
        T2 := source morph;
        T1 := target morph;
        morphFT2 := node.NodeRight.QuotMap;
        nestingRight = quotNesting(T2, T1, morphFT2, morph);
    );
    infoRight := if cellRight === null then null else cellInfo(cellRight, nestingRight);

    nestingDown = null;
    if not(node.MapFromDown === null) then (
        morph := node.MapFromDown;
        T2 := source morph;
        T1 := target morph;
        morphFT2 := node.NodeDown.QuotMap;
        nestingDown = quotNesting(T2, T1, morphFT2, morph);
    );
    infoDown := if cellDown === null then null else cellInfo(cellDown, nestingDown);

    cell := nestingCell(target node.QuotMap, node.QuotMap, CellInfoRight=>infoRight, CellInfoDown=>infoDown);
    return cell;
)

doubleNestedQuotSchemePoint GraphNode := node -> (
    F := source node.QuotMap;
    nestBase = constructNestingCell(node);
    return doubleNestedQuotSchemePoint(F, nestBase);
)


-- TODO: better name for this methods
constructNestingArrays = method(TypicalValue => List, Options => {OnRows => true})

constructNestingArrays NestingCell := opts -> nestCell -> (
    nextArray := nestBase;
    arrays := while not(nextArray === null) list (
        nextCell := nextArray;
        nextArray = nextArray#(if opts.OnRows then CellDown else CellRight);
        while not(nextCell === null) list (
            nesting := nextCell#(if opts.OnRows then NestingRight else NestingDown);
            nextCell = nextCell#(if opts.OnRows then CellRight else CellDown);
            if nesting === null then break;
            nesting
        )
    );
    return arrays;
)

doubleNestedQuotSchemePoint (Module, NestingCell) := (F, nestBase) -> (
    rows := constructNestingArrays(nestBase, OnRows=>true);
    cols := constructNestingArrays(nestBase, OnRows=>false);
    l := homBasesLengths(nestBase);

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

    nextRow := nestBase;
    row := 0;
    while not(nextRow === null) do (
        nextCell := nextRow;
        col := 0;
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
    lengths := {};
    lengthsCumulative := {};
    flatIndices := new MutableHashTable from {};

    iterateRows(
        nestBase, 
        LoopOp => (cell, row, col) -> (
                flatIndices#(row, col) = length(lengths);
                lengthsCumulative = append(lengthsCumulative, sum(lengths));
                lengths = append(lengths, basisLength(cell.HomBasis));
            )
        );
    lengthsCumulative = append(lengthsCumulative, sum(lengths));
    return {lengthsCumulative, flatIndices};
)


doubleNestingConstraint = method(TypicalValue => Matrix)

doubleNestingConstraint (DoubleNestedQuotSchemePoint, QuotNesting, ZZ, ZZ) := (p, nesting, sourceIndex, targetIndex) -> (
    totalM := p.BasisLengthsCumulative_-1;
    constr := nestedQuotTangentConstraints(nesting);
    n := numrows constr.ConstraintOnSource;
    C := zeroMutableMatrix(n, totalM);
    matrixWriteSlice(C, 0, p.BasisLengthsCumulative_sourceIndex, constr.ConstraintOnSource);
    matrixWriteSlice(C, 0, p.BasisLengthsCumulative_targetIndex, constr.ConstraintOnTarget);
    return new Matrix from C;
)


tangentSpace DoubleNestedQuotSchemePoint := p -> (
    totalM := p.BasisLengthsCumulative_-1;
    totalC := zeroMatrix(0, totalM);

    nextRow := p.NestingBase;
    row := 0;
    while not(nextRow === null) do (
        nextCell := nextRow;
        col := 0;
        while not(nextCell === null) do (
            -- nesting to the right
            if not(nextCell.CellRight === null) then (
                targetIndex := p.CellFlatIndices#(row, col);
                sourceIndex := p.CellFlatIndices#(row, col+1);
                C := doubleNestingConstraint(p, nextCell.NestingRight, sourceIndex, targetIndex);
                totalC = totalC || C;
            );
            if not(nextCell.CellDown === null) then (
                targetIndex := p.CellFlatIndices#(row, col);
                sourceIndex := p.CellFlatIndices#(row+1, col);
                C := doubleNestingConstraint(p, nextCell.NestingDown, sourceIndex, targetIndex);
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
    p := morphisms_0;
    F := source p;
    node := graphNode(p);
    nestings = for i from 1 to length(morphisms)-1 list (
        p = morphisms_i * p;
        node = graphNode(p, Right=>nodeInfo(node, morphisms_i));
    );
    return doubleNestedQuotSchemePoint(node);
)


nestedHilbSchemePoint2 = method(TypicalValue => DoubleNestedQuotSchemePoint)

nestedHilbSchemePoint2 List := modules -> (
    R := ring modules_0;
    r := rank modules_0;
    idR := idMatrix(r)**R;
    morphisms := for i from 0 to length(modules)-2 list (
        map(modules_(i+1), modules_i, idR)
    );
    return nestedQuotSchemePoint2(morphisms);
)






makeNode = method(TypicalValue => GraphNode, Options => {Memory=>new MutableHashTable from {}})
makeNode (Module, ZZ, ZZ, List) := opts -> (F, row, col, data) -> (
    memoryFetch := (r, c) -> (if not member((r, c), keys opts.Memory) then opts.Memory#(r, c) = makeNode(F, r, c, data, Memory=>opts.Memory); return opts.Memory#(r, c););
    getModule := (r, c) -> (if r%2 != 0 or c%2 != 0 then error("invalid indices"); cell = data_r_c; if instance(cell, Module) then return cell; if instance(cell_0, Module) then return cell_0; error("couldn't get module"));
    nextModule := (r, c, onRow) -> (if onRow then (r, c+2) else (r+2, c));

    createAdjacentNode = (r, c, onRow) -> (
        n := memoryFetch(nextModule(r, c, onRow));
        trgt := data_r_c;
        src := getModule(nextModule(r, c, onRow));
        R := ring trgt;
        mat := if onRow then data_r_(c+1) else data_(r+1)_(c//2);
        if instance(mat, List) then mat = matrix mat;
        f := map(trgt, src, mat**R);
        if not isWellDefined f then (print f; print (source f); print (target f); error "invalid map";);
        nodeAdj := nodeInfo(n, f);
        qAdj := f * n.QuotMap;
        return (nodeAdj, qAdj);
    );
    
    -- right branch
    nodeRight := null; qRight := null;
    if length data_row > col+1 then (
        l := createAdjacentNode(row, col, true);
        nodeRight = l_0; qRight = l_1;
    );
    -- down branch
    nodeDown := null; qDown := null;
    if length data > row+1 and length data_(row+2) > col then (
        l := createAdjacentNode(row, col, false);
        nodeDown = l_0; qDown = l_1;
    );

    isLeaf := nodeRight === null and nodeDown === null;
    if isLeaf then (
        cell := data_row_col;
        QModule := cell_0; q := cell_1;
        if instance(q, List) then q = matrix q;
        R := ring QModule;
        q = map(QModule, F, q**R);
        return graphNode(q);
    );

    if qRight =!= null and qDown =!= null and qRight != qDown then error("noncommutative");
    q := if qRight =!= null then qRight else qDown;

    return graphNode(q, Right=>nodeRight, Down=>nodeDown);
)

doubleNestedQuotSchemePoint (Module, List) := (F, data) -> (
    node := makeNode(F, 0, 0, data, Memory=>new MutableHashTable from {});
    return doubleNestedQuotSchemePoint(node);
)





YoungDiagram = new Type of List

getShape = method(TypicalValue=>YoungDiagram)
getBoxModule = method(TypicalValue=>Module)
getBoxQuotientMap = method(TypicalValue=>Matrix)

getShape DoubleNestedQuotSchemePoint := p -> (
    return for row in p.NestedRows list length(row)+1;
)

getBoxModule (DoubleNestedQuotSchemePoint, ZZ, ZZ) := (p, i, j) -> (
    -- TODO: make sure it can handle corner cases like simple nesting
    row := p.NestedRows_j;
    if i == 0 then (
        return (row_0).ModuleTarget;
    );
    return (row_(i-1)).ModuleSource;
)

getBoxQuotientMap (DoubleNestedQuotSchemePoint, ZZ, ZZ) := (p, i, j) -> (
    -- TODO: make sure it can handle corner cases like simple nesting
    row := p.NestedRows_j;
    if i == 0 then (
        return (row_0).QuotientTarget;
    );
    return (row_(i-1)).QuotientSource;
)


diagonal = method(TypicalValue=>Module)
diagonal Module := M -> (
    F := ambient M;
    A := relations M;
    (D, P, Q) := smithNormalForm A;
    return F/image(D);
)



-- toString DoubleNestedQuotSchemePoint := p -> (
--     for i from 1 to length(p.NestedRows) list 0;
-- )

-- net DoubleNestedQuotSchemePoint := p -> (
--     for i from 1 to length(p.NestedRows) list for j from 1 to lengthnet(getBoxModule());
-- )