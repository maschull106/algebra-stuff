from __future__ import annotations
from sympy import Symbol
import sympy
from typing import TYPE_CHECKING, List, Union
from fractions import Fraction
import algebra_stuff.groebner_polynomial as gb  # real import
if TYPE_CHECKING:   # fake import, only for annotations
    from .groebner_polynomial import GroebnerPolynomial
#from .common import *
from .common import Params, ExecTimes, init_globals, focus_poly_ring, focus_base_ring, infer_poly_ring, infer_base_ring, get_global_scope, set_global_scope, revert_global_scope, filter_zero, list_add, null_space, exact_null_space, decide_dtype

# class GroebnerPolynomial:
#     def __new__(self, *args, **kwargs):
#         raise RuntimeError("cannot instantiate not initialized class GroebnerPolynomial")


class MonomialOrder:
    def eval(self, monomial: Monomial):
        pass
    
    def eq(self, m1: Monomial, m2: Monomial) -> bool:
        return self.eval(m1) == self.eval(m2)
    
    def gt(self, m1: Monomial, m2: Monomial) -> bool:
        return self.eval(m1) > self.eval(m2)
    
    def lt(self, m1: Monomial, m2: Monomial) -> bool:
        return self.eval(m1) < self.eval(m2)


class LexicographicOrder(MonomialOrder):
    def eval(self, monomial: Monomial):
        return monomial.degrees


class GradedReverseLexicographicOrder(MonomialOrder):
    def eval(self, monomial: Monomial):
        return (monomial.total_degree(), [-monomial.degrees[i] for i in range(monomial.var_count()-1, -1, -1)])


lex = LexicographicOrder()
degrevlex = GradedReverseLexicographicOrder()


def base_decomp(n: int, b: int = 10, descending_order: bool = True):
    decomp = []
    while n != 0:
        decomp.append(n%b)
        n //= b
    if descending_order: decomp.reverse()
    return decomp


class Scalar:
    FLOAT = 0
    COMPLEX = 1
    FRACTION = 2
    MODE = FLOAT
    TYPES = (int, float, complex, Fraction, sympy.Float, sympy.Integer, sympy.core.numbers.One, sympy.core.numbers.Zero)

    @staticmethod
    def make(t):
        #t = t[0]
        if Scalar.MODE == Scalar.FLOAT:
            return t
        if Scalar.MODE == Scalar.FRACTION:
            if isinstance(t, sympy.Float):
                t = float(t)
            return Fraction(t)
        if Scalar.MODE == Scalar.COMPLEX:
            print("Warning: COMPLEX scalar mode not 100% functional")
            return t
            return t[0] + 1j*t[1]
        raise ValueError
    
    
class Monomial:
    superscript_numbers = {0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹'}
    
    def __init__(self, symbols: List[Symbol], degrees: List[int]):
        # make sure degrees has the same size as symbols
        degrees = degrees[:len(symbols)]
        degrees.extend([0]*(len(symbols)-len(degrees)))
        
        self.symbols = symbols
        self.degrees = degrees
    
    def total_degree(self) -> int:
        return sum(self.degrees)
    
    def var_count(self) -> int:
        return len(self.symbols)
    
    def __eq__(self, other: Monomial):
        return self.symbols == other.symbols and self.degrees == other.degrees
    
    def __mul__(self, other: Union[Monomial, MonomialWithCoef, GroebnerPolynomial, Scalar]):
        if isinstance(other, Monomial):
            degrees = [d1+d2 for d1, d2 in zip(self.degrees, other.degrees)]
            return Monomial(self.symbols, degrees)
        if isinstance(other, MonomialWithCoef):
            return MonomialWithCoef(coef=1, monomial=self) * other
        if isinstance(other, gb.GroebnerPolynomial):
            return other * self
        if isinstance(other, Scalar.TYPES):
            return MonomialWithCoef(other, self)
        print("Encountered unexpected type in multiplication:", other, type(other))
        raise TypeError
    
    def __rmul__(self, other):
        return self * other
    
    def __truediv__(self, other: Monomial):
        degrees = []
        for d1, d2 in zip(self.degrees, other.degrees):
            if d1 < d2:
                raise ValueError("the first monomial is not a multiple of the second")
            degrees.append(d1-d2)
        return Monomial(self.symbols, degrees)
    
    def gcd(self, other: Monomial):
        degrees = [min(d1, d2) for d1, d2 in zip(self.degrees, other.degrees)]
        return Monomial(self.symbols, degrees)
    
    def lcm(self, other: Monomial):
        degrees = [max(d1, d2) for d1, d2 in zip(self.degrees, other.degrees)]
        return Monomial(self.symbols, degrees)
    
    def is_multiple(self, other: Monomial):
        return all(d1 >= d2 for d1, d2 in zip(self.degrees, other.degrees))
    
    def is_constant(self) -> bool:
        return all(d == 0 for d in self.degrees)
    
    def __repr__(self):
        if self.is_constant():
            return "1"
        exp_format = lambda e: "" if e == 1 else "".join(map(Monomial.superscript_numbers.get, base_decomp(e)))
        return "".join(repr(x) + exp_format(e) for x, e in zip(self.symbols, self.degrees) if e > 0)
    
    def __hash__(self):
        def enum():
            for symb, d in zip(self.symbols, self.degrees):
                for _ in range(d):
                    yield symb
        
        t = tuple(hash(symb.name) for symb in enum())
        return hash(t)
    
    def macaulay2_repr(self):
        if self.is_constant():
            return "1"
        return "*".join(repr(x) + "^" + repr(e) for x, e in zip(self.symbols, self.degrees) if e > 0)


class MonomialWithCoef:
    def __init__(self, coef: Scalar, monomial: Monomial):
        self.coef = coef
        self.monomial = monomial
    
    def __eq__(self, other: MonomialWithCoef):
        return self.coef == other.coef and self.monomial == other.monomial
    
    def __mul__(self, other: Union[MonomialWithCoef, Monomial, GroebnerPolynomial]):
        if isinstance(other, MonomialWithCoef):
            return MonomialWithCoef(self.coef*other.coef, self.monomial*other.monomial)
        if isinstance(other, Monomial):
            return other*self
        if isinstance(other, gb.GroebnerPolynomial):
            return other*self
        raise TypeError
    
    def __rmul__(self, other):
        return self * other
    
    def __truediv__(self, other: Union[MonomialWithCoef, Monomial, Scalar]):
        if isinstance(other, MonomialWithCoef):
            return MonomialWithCoef(self.coef/other.coef, self.monomial/other.monomial)
        if isinstance(other, Monomial):
            return MonomialWithCoef(self.coef, self.monomial/other)
        if isinstance(other, Scalar.TYPES):
            return MonomialWithCoef(self.coef/other, self.monomial)
        raise TypeError
    
    def coef_repr(self):
        if self.coef == 1 and not self.monomial.is_constant():
            s_coef = ""
        elif self.coef == -1 and not self.monomial.is_constant():
            s_coef = "-"
        else:
            if isinstance(self.coef, complex) and self.coef.imag == 0.0:
                if (self.coef.real.is_integer()):
                    s_coef = repr(int(self.coef.real))
                else:
                    s_coef = repr(self.coef.real)
            elif isinstance(self.coef, float) and self.coef.is_integer():
                s_coef = repr(int(self.coef))
            elif isinstance(self.coef, Fraction) and self.coef.as_integer_ratio()[1] == 1:
                s_coef = repr(int(self.coef))
            else:
                s_coef = repr(self.coef)
        return s_coef
    
    def __repr__(self):
        if self.monomial.is_constant():
            s_mon = ""
        else:
            s_mon = repr(self.monomial)
        s_coef = self.coef_repr()
        return s_coef + s_mon

    def macaulay2_repr(self):
        if self.monomial.is_constant():
            return self.coef_repr()
        s_mon = self.monomial.macaulay2_repr()
        #s_coef = self.coef_repr()
        #s_coef = str(int(self.coef))    # only works for integer coefs
        n, d = Fraction(self.coef).as_integer_ratio()
        s_coef = str(n) + "/" + str(d)
        return s_coef + "*" + s_mon
        