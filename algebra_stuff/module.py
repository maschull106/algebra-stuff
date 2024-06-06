from __future__ import annotations
from .polynomial_ring import *
import numpy as np
from typing import Callable


Morphism = Callable[[GroebnerPolynomial], GroebnerPolynomial]


class FiniteDimRing:
    def __init__(self, order: MonomialOrder = degrevlex, symbols: List[Symbol] = None):
        self.order = order
        self.symbols = symbols
        self.basis = self._get_basis()
        focus_base_ring(self)
    
    def _get_basis(self) -> List[GroebnerPolynomial]:
        pass
    
    @property
    def dim(self):
        """dimension as a vector space"""
        return len(self.basis)
    
    def __eq__(self, other: FiniteDimRing):
        return False
    

class QuotientRing(FiniteDimRing):
    def __init__(self, ideal: PolyRingIdeal):
        self.ideal = ideal
        super().__init__(ideal.order, ideal.symbols)
    
    def _get_basis(self) -> List[GroebnerPolynomial]:
        return self.ideal.degree()
    
    def __eq__(self, other: FiniteDimRing):
        return isinstance(other, QuotientRing) and self.ideal == other.ideal
    
    def __repr__(self):
        if Params.focus_on_display: # TODO: move this in the base class
            focus_base_ring(self)
        s_ring = repr(self.ideal.base)
        s_ideal = repr(self.ideal)
        length = max(len(s_ring), len(s_ideal))
        pad = lambda s: " "*((length-len(s))//2) + s + " "*((length-len(s)+1)//2)
        return pad(s_ring) + "\n" + "―"*length + "\n" + pad(s_ideal)


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
    
    @property
    def dim(self):
        """dimension as a vector space"""
        return len(self.basis)
    
    def contains(self, f: GroebnerPolynomial) -> bool:
        return False
    
    def random_element(self) -> GroebnerPolynomial:
        coefs = [np.random.randint(-5, 5) for _ in range(self.dim)]
        return self.from_basis(coefs)
    
    def to_basis(self, f: GroebnerPolynomial) -> List[Scalar]:
        pass
    
    def from_basis(self, vect: List[Scalar]) -> GroebnerPolynomial:
        zero = GroebnerPolynomial.make(0, order=self.base_ring.order, symbols=self.base_ring.symbols)
        return sum((c*g for c, g in zip(vect, self.basis)), start=zero)
    
    def get_matrix_representation(self, phi: Morphism, dtype="float64") -> np.ndarray:
        # for phi an endomorphism on the module, compute the matrix representation
        matrix = [self.to_basis(phi(f)) for f in self.basis]
        matrix = np.array(matrix, dtype=dtype).T
        return matrix
    
    def get_matrices_representation(self, dtype="float64") -> List[np.ndarray]:
        # for every g in the basis of the base ring, compute the matrix representation of the transformation induced by g
        return [self.get_matrix_representation(lambda f: g*f, dtype=dtype) for g in self.base_ring.basis]
    
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
            if f.lm == g.lm:
            #if f.lm.is_multiple(g.lm):
                c = f.lc/g.lc
                f = f - c*g
            coefs.append(c)
        if not f.is_zero():
            print("ERROR: NOT ZERO IN BASIS", orig_f, "  ||  ", f)
            additional_coefs = self.to_basis(f)
            coefs = list_add(coefs, additional_coefs)
        return coefs


class RingQuotientModule(ModuleFromIdeal):
    # of the form R/I
    def __init__(self, base_ring: FiniteDimRing, ideal: PolyRingIdeal):
        super().__init__(base_ring, ideal)
    
    def _get_basis(self) -> List[GroebnerPolynomial]:
        return self.structure_ideal.degree()
    
    def contains(self, f: GroebnerPolynomial) -> bool:
        return True
    
    def __repr__(self):
        s_ring = repr(self.structure_ideal.base)
        s_ideal = repr(self.structure_ideal)
        length = max(len(s_ring), len(s_ideal))
        pad = lambda s: " "*((length-len(s))//2) + s + " "*((length-len(s)+1)//2)
        s = pad(s_ring) + "\n" + "―"*length + "\n" + pad(s_ideal)
        s_base = repr(self.base_ring)
        return concat_multiline_strings(s, s_base)


class IdealQuotientModule(ModuleFromIdeal):
    def __init__(self, base_ring: FiniteDimRing, top_ideal: PolyRingIdeal, bot_ideal: PolyRingIdeal):
        if not top_ideal.contains_ideal(bot_ideal):
            raise ValueError
        self.top_ideal = top_ideal
        super().__init__(base_ring, bot_ideal)
    
    def _get_basis(self) -> List[GroebnerPolynomial]:
        return self.structure_ideal.degree_for_base(base_ideal=self.top_ideal)
    
    def contains(self, f: GroebnerPolynomial) -> bool:
        return self.top_ideal.contains(f)
    
    def __repr__(self):
        s_top_ideal = repr(self.top_ideal)
        s_bot_ideal = repr(self.structure_ideal)
        length = max(len(s_top_ideal), len(s_bot_ideal))
        pad = lambda s: " "*((length-len(s))//2) + s + " "*((length-len(s)+1)//2)
        s = pad(s_top_ideal) + "\n" + "―"*length + "\n" + pad(s_bot_ideal)
        s_base = repr(self.base_ring)
        return concat_multiline_strings(s, s_base)


def concat_multiline_strings(s1: str, s2: str, sep="module of", sep_index=1):
    sep = " "*4 + sep + " "*4
    _sep = " "*len(sep)
    l1 = s1.split("\n")
    l2 = s2.split("\n")
    length1 = len(l1[0]) if len(l1) > 0 else 0
    final_lines = []
    for i, (a, b) in enumerate(zip(l1, l2)):
        current_sep = sep if i == sep_index else _sep
        final_lines.append(a+current_sep+b)
    for i in range(len(l1), len(l2)):
        current_sep = sep if i == sep_index else _sep
        final_lines.append(" "*length1 + current_sep + l2[i])
    return "\n".join(final_lines)
