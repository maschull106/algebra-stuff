from __future__ import annotations
from .groebner_polynomial import *
import numpy as np


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
    
    def sort_list(self, l: List[GroebnerPolynomial]) -> List[GroebnerPolynomial]:
        return sorted(l, key=lambda f: self.order.eval(f.lm), reverse=True)
        
    def ideal(self, *gens: GroebnerPolynomial) -> PolyRingIdeal:
        return PolyRingIdeal(self.symbols, gens, self.order)
    
    def __repr__(self):
        return f"{repr(self.base)}[{', '.join(map(repr, self.symbols))}]"


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
        self.groebner_basis = self.sort_list(self.groebner_basis)
    
    def sort_list(self, l: List[GroebnerPolynomial]) -> List[GroebnerPolynomial]:
        return sorted(l, key=lambda f: self.order.eval(f.lm), reverse=True)
    
    def contains(self, f: GroebnerPolynomial) -> bool:
        f = GroebnerPolynomial.make(f, order=self.order, symbols=self.symbols)
        f_reduced = self.reduce(f, self.groebner_basis)
        return f_reduced.is_zero()
    
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
    
    def degree_over_square(self) -> List[GroebnerPolynomial]:
        if not self.has_max_radical():
            return float("inf")
        
        def is_mutliple(poly: GroebnerPolynomial):
            for f in square.groebner_basis:
                if poly.lm.is_multiple(f.lm):
                    return True
            quotient_gens.append(poly)
            return False
        def stop_cond(mon_degrees):
            monomial = Monomial(self.symbols, mon_degrees)
            stop = [is_mutliple(g*monomial) for g in self.groebner_basis]
            return all(stop)
        
        square = self**2
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
        return f"<{', '.join(map(repr, self.gens))}>"
    
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

    def __pow__(self, d: int):
        def prod(set):
            p = 1
            for g in set:
                p *= g
            return p
        gens = [prod(s) for s in choices(self.gens, d)]
        return PolyRingIdeal(self.symbols, gens, self.order)
    

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


class Quotient:
    def __init__(self, ideal: PolyRingIdeal):
        self.ideal = ideal
        self.basis = self._get_basis()
        
    def _get_basis(self) -> List[GroebnerPolynomial]:
        return self.ideal.degree()


class ModuleFromIdeal:
    # module whose structure is defined from an ideal in some undefined way (see child classes for concrete examples)
    def __init__(self, base_ring: Quotient, ideal: PolyRingIdeal):
        self.base_ring = base_ring
        self.structure_ideal = ideal
        self.basis = self._get_basis()
        
    def _get_basis(self) -> List[GroebnerPolynomial]:
        pass
    
    def contains(self, f: GroebnerPolynomial) -> bool:
        return False
    
    def to_basis(self, f: GroebnerPolynomial) -> GroebnerPolynomial:
        f = GroebnerPolynomial.make(f, order=self.structure_ideal.order, symbols=self.structure_ideal.symbols)
        orig_f = f
        f = self.structure_ideal.reduced(f)
        coefs = []
        for g in self.basis:
            c = 0
            if f.lm.is_multiple(g.lm):
                c = f.lc/g.lc
                f = f - c*g
            coefs.append(c)
        if not f.is_zero():
            print("WTFFFF", orig_f, f)
        return coefs
    
    def get_matrix_representation(self, g: GroebnerPolynomial):
        # for g an element of the base ring, compute the matrix representation of the transformation induced by g
        matrix = [self.to_basis(g*f) for f in self.basis]
        matrix = np.array(matrix).T
        return matrix
    
    def get_matrices_representation(self):
        return [self.get_matrix_representation(g) for g in self.base_ring.basis]
    
    def construct_endo_matrix(self):
        n = len(self.basis)
        k = len(self.base_ring.basis)
        F = self.get_matrices_representation()
        M = np.zeros((n*n*k, n*n))
        for f_ind in range(k):
            for i in range(n):
                for j in range(n):
                    for l in range(n):
                        ind = f_ind*n**2 + i*n + j
                        M[ind, i*n+l] += F[f_ind][l, j]
                        M[ind, l*n+j] -= F[f_ind][i, l]
                        print(i*n+l)
                        print(l*n+j)
                        print()
        return M      


class RingQuotient(ModuleFromIdeal):
    # of the form R/I
    def __init__(self, base_ring: Quotient, ideal: PolyRingIdeal):
        super().__init__(base_ring, ideal)
    
    def _get_basis(self) -> List[GroebnerPolynomial]:
        return self.structure_ideal.degree()
    
    def contains(self, f: GroebnerPolynomial) -> bool:
        return True
    

class IdealQuotient(ModuleFromIdeal):
    # of the for I/I^2
    def __init__(self, base_ring: Quotient, ideal: PolyRingIdeal):
        self.ideal = ideal
        super().__init__(base_ring, ideal**2)
    
    def _get_basis(self) -> List[GroebnerPolynomial]:
        return self.ideal.degree_over_square()
    
    def contains(self, f: GroebnerPolynomial) -> bool:
        return self.ideal.contains(f)


def hom(M: ModuleFromIdeal, N: ModuleFromIdeal) -> np.ndarray:
    ring = M.base_ring
    m = len(M.basis)
    n = len(N.basis)
    k = len(ring.basis)
    FM = M.get_matrices_representation()
    FN = N.get_matrices_representation()
    A = np.zeros((n*m*k, n*m))
    for f_ind in range(k):
        for i in range(n):
            for j in range(m):
                ind = f_ind*n*m + i*m + j
                for l in range(m):
                    A[ind, i*m+l] += FM[f_ind][l, j]
                for l in range(n):
                    A[ind, l*m+j] -= FN[f_ind][i, l]
    return A


def hom_rank(M: ModuleFromIdeal, N: ModuleFromIdeal, tol: float = None) -> int:
    m = len(M.basis)
    n = len(N.basis)
    h = hom(M, N)
    return m*n - np.linalg.matrix_rank(h, tol=tol)
