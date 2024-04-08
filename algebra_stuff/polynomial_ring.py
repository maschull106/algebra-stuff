from __future__ import annotations
from .groebner_polynomial import *


class PolyRing:
    def __init__(self, base=sympy.CC, n: int = 1, order: MonomialOrder = degrevlex):
        if n < 0:
            raise ValueError
        self.base = base
        self.n = n
        if n <= 4:
            s_symbols = ",".join("xyzw"[:n])
        else:
            s_symbols = ",".join([f"x{i}" for i in range(1, n+1)])
        self.symbols = sympy.symbols(s_symbols)
        self.order = order
        
    def ideal(self, *gens: GroebnerPolynomial) -> PolyRingIdeal:
        return PolyRingIdeal(self.symbols, gens, self.order)
    
    def __repr__(self):
        return f"{repr(self.base)}[{", ".join(map(repr, self.symbols))}]"


class PolyRingIdeal:
    def __init__(self, symbols: List[Symbol], gens: List[GroebnerPolynomial], order: MonomialOrder = degrevlex):
        self.symbols = symbols
        self.order = order
        self.gens = [GroebnerPolynomial.make(poly, order=order, symbols=symbols) for poly in gens]
        self.groebner_basis: List[GroebnerPolynomial] = None
        self._compute_groebner_basis()
    
    def _compute_groebner_basis(self):
        self.groebner_basis = self.buchberger_algo(self.gens)
        self.groebner_basis = self.basis_reduce(self.groebner_basis)
    
    def contains(self, f: GroebnerPolynomial) -> bool:
        f = GroebnerPolynomial.make(f, order=self.order, symbols=self.symbols)
        f_reduced = self.reduce(f, self.groebner_basis)
        return f_reduced.is_zero()
    
    def has_max_radical(self) -> bool:
        vars_represented = [False]*len(self.symbols)
        for f in self.groebner_basis:
            for i in range(len(self.symbols)):
                if f.lm.degrees[i] > 0:
                    if f.lm.degrees[i] == f.lm.total_degree:
                        vars_represented[i] = True
                    break
        return all(vars_represented)
    
    def degree(self):
        def is_contained(mon_degrees):
            monomial = Monomial(self.symbols, mon_degrees)
            if self.contains(monomial):
                return False
            quotient_gens.append(monomial)
            return True
        
        quotient_gens = []
        self._xplore_monomials(is_contained)
        return quotient_gens
    
    def _xplore_monomials(self, operation, current_degrees=None, current=0):
        # TODO: is bugged: need to go through all the changes at current at once, otherwise some monomials will be considered multiple times
        if current_degrees is None:
            current_degrees = [0]*len(self.symbols)
        if not operation(current_degrees):
            return
        if current == len(current_degrees):
            return
        current_degrees[current] += 1
        self._xplore_monomials(operation, current_degrees, current)
        current_degrees[current] -= 1
        self._xplore_monomials(operation, current_degrees, current+1)
    
    def __repr__(self):
        return f"<{", ".join(map(repr, self.gens))}>"
    
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
    def basis_reduce(F: List[GroebnerPolynomial]) -> List[GroebnerPolynomial]:
        PolyRingIdeal.basis_minimize(F)
        G = []
        for i, f in enumerate(F):
            _F = F[:i] + F[i+1:]
            f = PolyRingIdeal.total_reduce(f, _F)
            G.append(f)
        for i in range(len(G)):
            G[i] = G[i].to_monic()
        return G
