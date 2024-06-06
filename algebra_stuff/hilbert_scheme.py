from .hom import *


class HilbertScheme:
    """
    Representing Hilb(Spec R)
    """
    
    def __init__(self, R: PolyRing = None):
        if R is None:
            R = infer_poly_ring()
        self.R = R
    
    def tangent_space(self, I: PolyRingIdeal, use_scipy: bool = True):
        # TODO: sanity check that I is an ideal of R
        return HilbertSchemeTangentSpace(self, I, use_scipy=use_scipy)


class HilbertSchemeTangentSpace:
    def __init__(self, base: HilbertScheme, I: PolyRingIdeal, use_scipy: bool = True):
        """
        I ideal of a polynomial ring R, giving a zero dimensional subscheme Z of length n.
        Represents the tangent of the Hilbert scheme of points n at the point [Z]
        """
        self.base = base
        self.R = base.R
        self.I = I
        self.S = self.R//I  # all modules will be over this ring
        self.J = I/I**2
        self.O = self.R/I
        self.use_scipy = use_scipy
        self.hom_space = HomSpace(self.J, self.O, base=self.S, use_scipy=use_scipy)
    
    def _constraints(self) -> np.ndarray:
        return self.hom_space._get_constraints()
    
    def dim(self) -> int:
        return self.hom_space.dim()
        #return hom_rank(self.J, self.O, char=self.char)
    
    def basis(self) -> np.ndarray:
        return self.hom_space.basis()
        #return hom(self.J, self.O, char=self.char)
