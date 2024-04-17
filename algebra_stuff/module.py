from __future__ import annotations
from .polynomial_ring import *
import numpy as np
import scipy


class FiniteDimRing:
    def __init__(self):
        self.basis = self._get_basis()
        focus_base_ring(self)
    
    def _get_basis(self) -> List[GroebnerPolynomial]:
        pass
    
    def __eq__(self, other: FiniteDimRing):
        return False
    

class QuotientRing(FiniteDimRing):
    def __init__(self, ideal: PolyRingIdeal):
        self.ideal = ideal
        super().__init__()
    
    def _get_basis(self) -> List[GroebnerPolynomial]:
        return self.ideal.degree()
    
    def __eq__(self, other: FiniteDimRing):
        return isinstance(other, QuotientRing) and self.ideal == other.ideal


class Module:
    def __init__(self, base_ring: FiniteDimRing):
        # TODO: implement checking that it is indeed a module over the given ring
        self.base_ring = base_ring
        self.basis = self._get_basis()
    
    def _get_basis(self) -> List[GroebnerPolynomial]:
        pass
    
    def set_base(self, base_ring: FiniteDimRing):
        # TODO: implement checking that it is indeed a module over the given ring
        self.base_ring = base_ring
    
    def contains(self, f: GroebnerPolynomial) -> bool:
        return False
    
    def to_basis(self, f: GroebnerPolynomial) -> List[Scalar]:
        pass
    
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
        return M


class ModuleFromIdeal(Module):
    # module whose structure is defined from an ideal in some undefined way (see child classes for concrete examples)
    def __init__(self, base_ring: FiniteDimRing, ideal: PolyRingIdeal):
        self.structure_ideal = ideal
        super().__init__(base_ring)
    
    def to_basis(self, f: GroebnerPolynomial) -> List[Scalar]:
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


class RingQuotientModule(ModuleFromIdeal):
    # of the form R/I
    def __init__(self, base_ring: FiniteDimRing, ideal: PolyRingIdeal):
        super().__init__(base_ring, ideal)
    
    def _get_basis(self) -> List[GroebnerPolynomial]:
        return self.structure_ideal.degree()
    
    def contains(self, f: GroebnerPolynomial) -> bool:
        return True


class IdealQuotientModule(ModuleFromIdeal):
    def __init__(self, base_ring: FiniteDimRing, top_ideal: PolyRingIdeal, bot_ideal: PolyRingIdeal):
        if not top_ideal.contains_ideal(bot_ideal):
            raise ValueError
        self.top_ideal = top_ideal
        super().__init__(base_ring, bot_ideal)
    
    def _get_basis(self) -> List[GroebnerPolynomial]:
        return self.structure_ideal.degree_for_base(base_ideal=self.top_ideal)
    
    def contain(self, f: GroebnerPolynomial) -> bool:
        return self.top_ideal.contains(f)


def hom_calculation_complexity(M: ModuleFromIdeal, N: ModuleFromIdeal):
    ring = M.base_ring
    m = len(M.basis)
    n = len(N.basis)
    k = len(ring.basis)
    print("k =", k)
    print("m =", m)
    print("n =", n)
    print("total matrix entries:", k*m*n*m*n)
    

def hom_constraints(M: ModuleFromIdeal, N: ModuleFromIdeal) -> np.ndarray:
    if not M.base_ring == N.base_ring:
        raise TypeError("modules are not defined over the same ring")
    
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
    h = hom_constraints(M, N)
    return m*n - np.linalg.matrix_rank(h, tol=tol)


def null_space(A: np.ndarray):
    # TODO: make sure this is equivalent to just scipy.linalg.null_space(A)
    # just scipy.linalg.null_space(A) works but seems to be way slower for A with nb of rows (a lot) bigger than number of columns
    P, L, U = scipy.linalg.lu(A)
    basis = scipy.linalg.null_space(U)
    return basis


def hom(M: ModuleFromIdeal, N: ModuleFromIdeal) -> np.ndarray:
    C = hom_constraints(M, N)
    basis = null_space(C)
    return basis
