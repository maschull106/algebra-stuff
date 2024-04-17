from __future__ import annotations
from .monomial import *
import sympy
import copy


class GroebnerPolynomial:
    def __init__(self, monomials: List[MonomialWithCoef], symbols: List[Symbol], order: MonomialOrder):
        self.order = order
        self.monomials = self._monom_sort(copy.deepcopy(monomials))
        self.symbols = symbols
    
    @classmethod
    def make(cls, poly: Union[GroebnerPolynomial, sympy.Expr, Scalar], order: MonomialOrder = degrevlex, symbols: List[Symbol] = None) -> GroebnerPolynomial:
        if isinstance(poly, cls):
            return poly
        elif isinstance(poly, MonomialWithCoef):
            return cls([poly], symbols if symbols is not None else poly.symbols, order=order)
        elif isinstance(poly, Monomial):
            return cls([MonomialWithCoef(1, poly)], symbols if symbols is not None else poly.symbols, order=order)
        elif isinstance(poly, sympy.Expr):
            poly = poly.expand()
        elif isinstance(poly, int):
            poly = sympy.Integer(poly)
        elif isinstance(poly, float):
            poly = sympy.Float(poly)
        elif isinstance(poly, complex):
            poly = poly.real + poly.imag*sympy.I
        poly = poly.expand()
        terms, p_symbols = poly.as_terms()
        if symbols is None:
            symbols = p_symbols
        else:
            if not set(p_symbols).issubset(set(symbols)):
                raise ValueError("encountered unknown symbols")
            symbol2index = {s: i for i, s in enumerate(symbols)}
            for i in range(len(terms)):
                expr, (coef, degrees, _) = terms[i]
                new_degrees = [0]*len(symbols)
                for s, d in zip(p_symbols, degrees):
                    new_degrees[symbol2index[s]] = d
                terms[i] = (expr, (coef, tuple(new_degrees), _))
            
        monomials = [
            MonomialWithCoef(
                coef=Scalar.make(t[1][0]),
                monomial=Monomial(symbols, list(t[1][1]))
            )
            for t in terms
        ]
        return cls(monomials, symbols, order)
    
    def _leading(self) -> MonomialWithCoef:
        if self.is_zero():
            return MonomialWithCoef(0, Monomial(self.symbols, []))
        return self.monomials[0]
    
    def get_coef(self, i: int) -> Scalar:
        return self.monomials[i].coef
    
    def get_monomial(self, i: int) -> Monomial:
        return self.monomials[i].monomial
    
    @property
    def lc(self) -> Scalar:
        return self._leading().coef
    
    @property
    def lm(self) -> Monomial:
        return self._leading().monomial
    
    def degree(self) -> int:
        return self.lm.total_degree()
    
    def _monom_sort(self, monomials: List[MonomialWithCoef]) -> List[MonomialWithCoef]:
        # works both in place and not in place
        monomials.sort(key=lambda m: self.order.eval(m.monomial), reverse=True)
        return monomials

    def _monom_collapse(self, monomials: List[MonomialWithCoef]) -> List[MonomialWithCoef]:
        # does not work in place
        monomials = copy.deepcopy(monomials)
        i = 0
        while i < len(monomials):
            j = i+1
            while j < len(monomials) and self.order.eq(monomials[i].monomial, monomials[j].monomial):
                j += 1
            if j > i+1:
                monomials[i].coef = sum(monomials[k].coef for k in range(i, j))
                if monomials[i].coef == 0:
                    i -= 1
                monomials = monomials[:i+1] + monomials[j:]
            i += 1
        return monomials
    
    def _monom_arrange(self, monomials: List[MonomialWithCoef]) -> List[MonomialWithCoef]:
        self._monom_sort(monomials)
        return self._monom_collapse(monomials)
    
    def _make_and_check(self, other) -> GroebnerPolynomial:
        other = GroebnerPolynomial.make(other, order=self.order, symbols=self.symbols)
        if not self.symbols == other.symbols:
            raise ValueError("symbols between polynomials do not match")
        return other
    
    def __eq__(self, other: GroebnerPolynomial) -> bool:
        return self.monomials == other.monomials
    
    def __add__(self, other) -> GroebnerPolynomial:
        other = self._make_and_check(other)
        monomials = self.monomials + other.monomials
        monomials = self._monom_arrange(monomials)
        return GroebnerPolynomial(monomials, self.symbols, self.order)
    
    def __radd__(self, other):
        return self + other
    
    def __neg__(self):
        monomials = [MonomialWithCoef(-m.coef, m.monomial) for m in self.monomials]
        return GroebnerPolynomial(monomials, self.symbols, self.order)
    
    def __sub__(self, other):
        return self.__add__(-other)
    
    def __mul__(self, other):
        other = self._make_and_check(other)
        monomials = [m1*m2 for m1 in self.monomials for m2 in other.monomials]
        monomials = self._monom_arrange(monomials)
        return GroebnerPolynomial(monomials, self.symbols, self.order)
    
    def __rmul__(self, other):
        return self * other
    
    def __truediv__(self, other: Scalar):
        if isinstance(other, (int, float, complex)):
            return GroebnerPolynomial([m / other for m in self.monomials], self.symbols, self.order)
        raise TypeError
    
    def __len__(self):
        return len(self.monomials)
    
    def is_zero(self) -> bool:
        if len(self.monomials) > 0 and all(m.coef == 0 for m in self.monomials):
            if len(self.monomials) > 1:
                print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
            return True
        return len(self.monomials) == 0
    
    def to_monic(self):
        c = self.lc
        return self / c
        
    def __repr__(self):
        # TODO: display - instead of +- when coefficient is negative
        if self.is_zero():
            return "0"
        sep_plus = " + " if Params.verbose else "+"
        sep_minus = " - " if Params.verbose else "-"
        s_poly = repr(self.monomials[0])
        for i in range(1, len(self.monomials)):
            s = repr(self.monomials[i])
            if s[0] == "-":
                s_poly += sep_minus
                s = s[1:]
            else:
                s_poly += sep_plus
            s_poly += s
        return s_poly
        #return sep_plus.join(map(repr, self.monomials))
    
    def __hash__(self):
        monom_hash = tuple(hash(m.monomial) for m in self.monomials)
        coef_hash = tuple(hash(m.coef) for m in self.monomials)
        return hash((monom_hash, coef_hash))