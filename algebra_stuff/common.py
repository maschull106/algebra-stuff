from __future__ import annotations
from typing import TYPE_CHECKING, List
#import algebra_stuff as alg     # real import
if TYPE_CHECKING:   # fake import, only for annotations
    from algebra_stuff import *


class Globals:
    INSTANCE = None
    
    def __new__(cls):
        if cls.INSTANCE is None:
            cls.INSTANCE = super().__new__(cls)
            return cls.INSTANCE
        raise TypeError("'Globals' class can only be instantiated once")
    
    def initialize(self):
        self._scope = globals()  # scope where declaration of global variables will happen
        self._fields = []    # TODO: implement fields
        self._poly_rings: List[PolyRing] = []
        self._base_rings: List[FiniteDimRing] = []   # for module construction (only rings that are finite dimensional as vector spaces for now)
    
    def focus_poly_ring(self, ring: PolyRing):
        if ring in self._poly_rings:
            self._poly_rings.remove(ring)
        self._poly_rings.append(ring)
    
    def focus_base_ring(self, ring: FiniteDimRing):
        if ring in self._base_rings:
            self._base_rings.remove(ring)
        self._base_rings.append(ring)
    
    def infer_base_ring(self):
        # TODO: add context for infering, such as latest ring for which module can be constructed or something
        if len(self._base_rings) == 0:
            raise LookupError("no known rings")
        return self._base_rings[-1]
    
    def infer_poly_ring(self):
        # TODO: add context for infering
        if len(self._poly_rings) == 0:
            raise LookupError("no known polynomial rings")
        return self._poly_rings[-1]
    
    def get_scope(self) -> dict:
        return self._scope
    
    def set_scope(self, scope: dict):
        # set the scope where declaration of global variables will happen
        self._scope = scope


def init_globals():
    GLOBALS.initialize()

def focus_poly_ring(ring: PolyRing):
    GLOBALS.focus_poly_ring(ring)

def focus_base_ring(ring: FiniteDimRing):
    GLOBALS.focus_base_ring(ring)

def infer_poly_ring() -> PolyRing:
    return GLOBALS.infer_poly_ring()

def infer_base_ring() -> FiniteDimRing:
    return GLOBALS.infer_base_ring()

def get_global_scope() -> dict:
    return GLOBALS.get_scope()

def set_global_scope(scope: dict):
    GLOBALS.set_scope(scope)


GLOBALS = Globals()