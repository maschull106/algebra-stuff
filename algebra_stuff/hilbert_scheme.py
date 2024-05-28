from .hom import *


class HilbertScheme:
    """
    Representing Hilb(Spec R)
    """
    
    def __init__(self, R: PolyRing):
        self.R = R
    
    def tangent_space(self, I: PolyRingIdeal):
        # TODO: sanity check that I is an ideal of R
        return HilbertSchemeTangentSpace(self, I)


class HilbertSchemeTangentSpace:
    def __init__(self, base: HilbertScheme, I: PolyRingIdeal, char: int = 0):
        """
        I ideal of a polynomial ring R, giving a zero dimensional subscheme Z of length n.
        Represents the tangent of the Hilbert scheme of points n at the point [Z]
        """
        self.base = base
        self.char = char
        self.R = base.R
        self.I = I
        self.S = self.R//I  # all modules will be over this ring
        self.J = I/I**2
        self.O = self.R/I
    
    # TODO: precompute constraints
    
    def constraints(self) -> np.ndarray:
        return hom_constraints(self.J, self.O)
    
    def dim(self) -> int:
        return hom_rank(self.J, self.O, char=self.char)
    
    def basis(self) -> np.ndarray:
        return hom(self.J, self.O, char=self.char)
