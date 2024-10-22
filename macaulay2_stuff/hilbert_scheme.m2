HilbertSchemePoint = new Type of HashTable
TangentSpace = new Type of HashTable
HilbertSchemeTangentSpace = new Type of TangentSpace


hilbertSchemePoint = method(TypicalValue => HilbertSchemePoint)

hilbertSchemePoint Ideal := I -> (
    if dim I != 0 then error "the ideal must correspond to a zero dimensional subscheme";

    R = ring I;
    return new HilbertSchemePoint from {
        base => R,
        structureIdeal => I,
        coordinateRing => R/I,
        len => degree I
    };
)


tangentSpace = method(TypicalValue => TangentSpace)

tangentSpace HilbertSchemePoint := p -> (
    return new HilbertSchemeTangentSpace from {
        originPoint => p,
        space => Hom(p.structureIdeal, p.coordinateRing)
    }
)


dim TangentSpace := T -> ( return degree T.space; )



superscripts = new HashTable from {
    "0" => "⁰",
    "1" => "¹",
    "2" => "²",
    "3" => "³",
    "4" => "⁴",
    "5" => "⁵",
    "6" => "⁶",
    "7" => "⁷",
    "8" => "⁸",
    "9" => "⁹"
};

toSuperscript = method(TypicalValue => String)
toSuperscript ZZ := x -> (
    digitsList = for digit in toString(x) list superscripts#digit;
    return concatenate(digitsList);
)

toString HilbertSchemePoint := p -> ( return "Point of Hilb" | toSuperscript(p.len) | "(A" | toSuperscript(dim p.base) | ") defined by ideal \n" | toString(p.structureIdeal); )

net HilbertSchemePoint := p -> (
    "Point of Hilb" | toSuperscript(p.len) | "(A" | toSuperscript(dim p.base) | ") defined by " | toString(p.structureIdeal)
)