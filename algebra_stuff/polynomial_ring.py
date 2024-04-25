from __future__ import annotations
from .groebner_polynomial import *
from typing import TYPE_CHECKING
import algebra_stuff.module as module   # real import
if TYPE_CHECKING:   # fake import, only for annotations
    from .module import QuotientRing, RingQuotientModule, IdealQuotientModule


class PolyRing:
    def __init__(self, base=sympy.CC, n: int = 1, order: MonomialOrder = degrevlex, symbols: List[sympy.Symbol]=None, make_symbols_global_vars: bool = True):
        if n < 0:
            raise ValueError
        self.base = base
        self.n = n
        if symbols is None or len(symbols) < n:
            if n <= 4:
                s_symbols = ",".join("xyzw"[:n])
            else:
                s_symbols = ",".join([f"x{i}" for i in range(1, n+1)])
            symbols = sympy.symbols(s_symbols)
        self.symbols = symbols[: n]
        self.order = order
        focus_poly_ring(self)
        if make_symbols_global_vars:
            self._make_symbols_global_vars()
    
    def _make_symbols_global_vars(self):
        # convenient but very dirty and potentially dangerous
        # TODO: doesn't work outside of the module's work
        for symbol in self.symbols:
            print(f"Setting global symbol '{symbol.name}'")
            get_global_scope()[symbol.name] = symbol
    
    def sort_list(self, l: List[GroebnerPolynomial]) -> List[GroebnerPolynomial]:
        return sorted(l, key=lambda f: self.order.eval(f.lm), reverse=True)
        
    def ideal(self, *gens: GroebnerPolynomial) -> PolyRingIdeal:
        return PolyRingIdeal(self, gens)
    
    def max_ideal(self):
        # give the ideal (x1, ..., xn)
        gens = self.symbols
        return PolyRingIdeal(self, gens)
    
    def __repr__(self):
        if Params.focus_on_display:
            focus_poly_ring(self)
        sep = ", " if Params.verbose else ","
        return f"{repr(self.base)}[{sep.join(map(repr, self.symbols))}]"
    
    def __floordiv__(self, ideal: PolyRingIdeal) -> QuotientRing:
        # TODO: implement some sanity checks
        if not isinstance(ideal, PolyRingIdeal):
            raise TypeError
        return module.QuotientRing(ideal)
    
    def __truediv__(self, ideal: PolyRingIdeal) -> RingQuotientModule:
        # TODO: implement some sanity checks
        if not isinstance(ideal, PolyRingIdeal):
            raise TypeError
        base_ring = infer_base_ring()
        return module.RingQuotientModule(base_ring=base_ring, ideal=ideal)


class PolyRingIdeal:
    def __init__(self, base: PolyRing, gens: List[GroebnerPolynomial]):
        self.base = base
        self.symbols = base.symbols
        self.order = base.order
        self.gens = [GroebnerPolynomial.make(poly, order=self.order, symbols=self.symbols) for poly in gens]
        self.groebner_basis: List[GroebnerPolynomial] = None
        self._compute_groebner_basis()
    
    def _compute_groebner_basis(self):
        self.groebner_basis = self.buchberger_algo(self.gens)
        self.groebner_basis = self.basis_reduce(self.groebner_basis)
        self.groebner_basis = self.sort_list(self.groebner_basis)
    
    def sort_list(self, l: List[GroebnerPolynomial]) -> List[GroebnerPolynomial]:
        return sorted(l, key=lambda f: self.order.eval(f.lm), reverse=True)
    
    def contains(self, f: GroebnerPolynomial) -> bool:
        f = GroebnerPolynomial.make(f, order=self.order, symbols=self.symbols)
        f_reduced = self.reduce(f, self.groebner_basis)
        return f_reduced.is_zero()
    
    def contains_ideal(self, ideal: PolyRingIdeal) -> bool:
        return all(self.contains(f) for f in ideal.groebner_basis)
    
    def has_max_radical(self) -> bool:
        vars_represented = [False]*len(self.symbols)
        for f in self.groebner_basis:
            for i in range(len(self.symbols)):
                if f.lm.degrees[i] > 0:
                    if f.lm.degrees[i] == f.lm.total_degree():
                        vars_represented[i] = True
                    break
        return all(vars_represented)
    
    def degree(self) -> List[GroebnerPolynomial]:
        # TODO: find better name for this function + combine with degree_for_base function
        if not self.has_max_radical():
            return float("inf")
        def is_multiple(mon_degrees):
            monomial = Monomial(self.symbols, mon_degrees)
            for f in self.groebner_basis:
                if monomial.is_multiple(f.lm):
                    return True
            quotient_gens.append(monomial)
            return False
        
        quotient_gens = []
        self._xplore_monomials(is_multiple)
        quotient_gens = list(map(GroebnerPolynomial.make, quotient_gens))
        quotient_gens = self.sort_list(quotient_gens)
        return quotient_gens
    
    def _xplore_monomials(self, stop_operation, current_degrees=None, current=0, new=True):
        if current_degrees is None:
            current_degrees = [0]*len(self.symbols)
        if new and stop_operation(current_degrees):
            return
        if current == len(current_degrees):
            return
        current_degrees[current] += 1
        self._xplore_monomials(stop_operation, current_degrees, current, new=True)
        current_degrees[current] -= 1
        self._xplore_monomials(stop_operation, current_degrees, current+1, new=False)
    
    def degree_for_base(self, base_ideal: PolyRingIdeal) -> List[GroebnerPolynomial]:
        if not self.has_max_radical():
            return float("inf")
        
        def is_mutliple(poly: GroebnerPolynomial):
            for f in self.groebner_basis:
                if poly.lm.is_multiple(f.lm):
                    return True
            quotient_gens.append(poly)
            return False
        def stop_cond(mon_degrees):
            monomial = Monomial(self.symbols, mon_degrees)
            stop = [is_mutliple(g*monomial) for g in base_ideal.groebner_basis]
            return all(stop)
        
        quotient_gens = []
        self._xplore_monomials(stop_cond)
        quotient_gens = list(map(GroebnerPolynomial.make, quotient_gens))
        quotient_gens = self.basis_soft_minimize(quotient_gens)
        quotient_gens = self.sort_list(quotient_gens)
        return quotient_gens
    
    def reduced(self, f: GroebnerPolynomial):
        return self.reduce(f, self.groebner_basis)
    
    def to_quotient_basis(self, f: GroebnerPolynomial):
        f = self.reduced(f)
        basis = self.degree()
        coefs = []
        for g in basis:
            c = 0
            if f.lm.is_multiple(g.lm):
                c = f.lc/g.lc
                f = f - c*g
            coefs.append(c)
        if not f.is_zero():
            print("WTFFFF", f)
        return coefs
    
    def quotient_square_representation_matrix(self, g):
        pass
    
    def __repr__(self):
        sep = ", " if Params.verbose else ","
        s = f"<{sep.join(map(repr, self.gens))}>"
        if len(s) > Params.long:
            return "[long ideal]"
        return s
    
    @staticmethod
    def reduction_step(f: GroebnerPolynomial, g: GroebnerPolynomial) -> GroebnerPolynomial:
        return f - f.lc/g.lc*f.lm/g.lm*g
    
    @staticmethod
    def general_reduction_step(f: GroebnerPolynomial, g: GroebnerPolynomial, i: int = 0):
        c = f.get_coef(i)
        m = f.get_monomial(i)
        return f - c/g.lc*m/g.lm*g
    
    @staticmethod
    def reduce(f: GroebnerPolynomial, G: List[GroebnerPolynomial]) -> GroebnerPolynomial:
        reducible = True
        while reducible:
            reducible = False
            for g in G:
                if f.lm.is_multiple(g.lm):
                    reducible = True
                    f = PolyRingIdeal.reduction_step(f, g)
                    if f.is_zero():
                        return f
                    break
        return f
    
    @staticmethod
    def total_reduce(f: GroebnerPolynomial, G: List[GroebnerPolynomial]) -> GroebnerPolynomial:
        i = 0
        while len(f) > i:
            reducible = True
            while reducible:
                if len(f) <= i:
                    break
                reducible = False
                for g in G:
                    if f.get_monomial(i).is_multiple(g.lm):
                        reducible = True
                        f = PolyRingIdeal.general_reduction_step(f, g, i)
                        break
            i += 1
        return f

    @staticmethod
    def buchberger_step(f: GroebnerPolynomial, g: GroebnerPolynomial, F: List[GroebnerPolynomial], new: List[GroebnerPolynomial]):
        gcd = f.lm.gcd(g.lm)
        S = 1/f.lc * g.lm / gcd * f - 1/g.lc * f.lm / gcd * g
        S = PolyRingIdeal.reduce(S, F)
        if not S.is_zero():
            new.append(S)
            
    @staticmethod
    def buchberger_algo(F: List[GroebnerPolynomial]):
        done = False
        G = F
        old = []
        new = F
        while not done:
            G = old + new
            next_new = []
            for f in old:
                for g in new:
                    PolyRingIdeal.buchberger_step(f, g, G, next_new)
            for i in range(len(new)-1):
                for j in range(i+1, len(new)):
                    f, g = new[i], new[j]
                    PolyRingIdeal.buchberger_step(f, g, G, next_new)
            done = len(next_new) == 0
            old = G
            new = next_new
        return G
    
    @staticmethod
    def basis_minimize(G: List[GroebnerPolynomial]):
        for i, f in enumerate(G):
            for j, g in enumerate(G):
                if i != j and f.lm.is_multiple(g.lm):
                    G.remove(f)
                    return PolyRingIdeal.basis_minimize(G)
        return G
    
    @staticmethod
    def basis_soft_minimize(G: List[GroebnerPolynomial]):
        for i, f in enumerate(G):
            for j, g in enumerate(G):
                if i != j and f.lm == g.lm:
                    G.remove(f)
                    return PolyRingIdeal.basis_soft_minimize(G)
        return G
            
    @staticmethod
    def basis_reduce(F: List[GroebnerPolynomial], minimize=True) -> List[GroebnerPolynomial]:
        if minimize:
            PolyRingIdeal.basis_minimize(F)
        G: List[GroebnerPolynomial] = []
        for i, f in enumerate(F):
            _F = F[:i] + F[i+1:]
            f = PolyRingIdeal.total_reduce(f, _F)
            G.append(f)
        for i in range(len(G)):
            G[i] = G[i].to_monic()
        return G
    
    def __eq__(self, other: PolyRingIdeal):
        if not isinstance(other, PolyRingIdeal):
            raise TypeError
        return self.groebner_basis == other.groebner_basis

    def __pow__(self, d: int):
        def prod(set):
            p = 1
            for g in set:
                p *= g
            return p
        gens = [prod(s) for s in choices(self.gens, d)]
        return PolyRingIdeal(self.base, gens)
    
    def __truediv__(self, other: PolyRingIdeal) -> IdealQuotientModule:
        # TODO: implement some sanity checks (some alr implemented in IdealQuotientModule, maybe move them here)
        if not isinstance(other, PolyRingIdeal):
            raise TypeError
        base_ring = infer_base_ring()
        return module.IdealQuotientModule(base_ring, top_ideal=self, bot_ideal=other)
    

def choices(S, size, start_index=0, current_choices=None):
    if current_choices is None:
        current_choices = []
    if size == 0:
        yield current_choices.copy()
        return
    if start_index >= len(S):
        return
    for i in range(start_index, len(S)):
        current_choices.append(S[i])
        yield from choices(S, size-1, start_index=i, current_choices=current_choices)
        current_choices.pop()


def ideal(*gens: GroebnerPolynomial) -> PolyRingIdeal:
    poly_ring = infer_poly_ring()
    return PolyRingIdeal(poly_ring, gens)
