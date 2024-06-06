from __future__ import annotations
from typing import TYPE_CHECKING, List, Dict, Callable
from dataclasses import dataclass
import time
import random
import numpy as np
import scipy
import sympy
from fractions import Fraction
#import algebra_stuff as alg     # real import
if TYPE_CHECKING:   # fake import, only for annotations
    from algebra_stuff import *


class Params:
    verbose = False
    long = 120
    focus_on_display = False


@dataclass
class TimeStep:
    t: float
    descr: str


class AA:
    def __init__(self, time_steps: List[TimeStep]):
        self.time_steps = time_steps
    
    def __repr__(self):
        return "\n".join(f"{step.descr}: {step.t} s" for step in self.time_steps)
    
    
class ExecTimes:
    _current_steps: Dict[str, List[TimeStep]] = {}
    _current_context: str = ""
    
    @classmethod
    def time_step(cls, description=""):
        # description should describe the task that will be executed until the next time step
        cls._current_steps.setdefault(cls._current_context, []).append(TimeStep(time.time(), description))

    @classmethod
    def track_time(cls, func: Callable):
        def time_tracked_f(*args, **kwargs):
            name = func.__name__
            context = name + str(random.randint(10**10, 10**11))
            past_context = cls._current_context
            cls._current_context = context
            res = func(*args, **kwargs)
            cls.time_step(description="end of func")
            time_steps = AA([
            TimeStep(cls._current_steps[context][i+1].t - cls._current_steps[context][i].t, cls._current_steps[context][i].descr)
                for i in range(len(cls._current_steps[context])-1)
            ])
            setattr(cls, name, time_steps)
            del cls._current_steps[context]
            cls._current_context = past_context
            return res
        
        return time_tracked_f


class Globals:
    INSTANCE = None
    
    def __new__(cls):
        if cls.INSTANCE is None:
            cls.INSTANCE = super().__new__(cls)
            return cls.INSTANCE
        raise TypeError("'Globals' class can only be instantiated once")
    
    def initialize(self):
        self._scope = globals()  # scope where declaration of global variables will happen
        self._scope_history = []    # for reverting to previous scope
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
        self._scope_history.append(self._scope)
        self._scope = scope
    
    def revert_global_scope(self):
        if len(self._scope_history) == 0:
            return
        self._scope = self._scope_history[-1]
        self._scope_history.pop()


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

def revert_global_scope():
    GLOBALS.revert_global_scope()


GLOBALS = Globals()


def filter_zero(C: np.ndarray):
    """removes rows that only contain zeros from the array"""
    if C.size == 0:
        return C
    try:
        return C[np.any(C, axis=1)]
    except:
        return C[np.any(C != 0, axis=1)]


def list_add(l1: list, l2: list):
    return [a+b for a, b in zip(l1, l2)]


def gaussian_elimination(A: np.ndarray, inplace=True):
    """put matrix A in column echelon form"""
    if not inplace:
        A = np.copy(A)
    n, m = A.shape
    limit = min(m, n)
    ech_col = 0
    for i in range(limit):
        if A[i, ech_col] == 0:
            nonzeros = np.where(A[i] != 0)[0]
            if nonzeros.size == 0:
                continue
            j = nonzeros[0]
            A[:, [ech_col, j]] = A[:, [j, ech_col]]
        c = A[i, ech_col]
        for j in range(m):
            if j != ech_col:
                d = A[i, j]
                A[: j] -= d/c * A[:, ech_col]
    
    return A


def null_space(A: np.ndarray) -> np.ndarray:
    # TODO: make sure this is equivalent to just scipy.linalg.null_space(A)
    # just scipy.linalg.null_space(A) works but seems to be way slower for A with nb of rows (a lot) bigger than number of columns
    P, L, U = scipy.linalg.lu(A)
    basis = scipy.linalg.null_space(U)
    return basis


def exact_null_space(A: np.ndarray) -> np.ndarray:
    n, m = A.shape
    A = np.concatenate([A, np.eye(m, dtype=decide_dtype(use_scipy=False))], axis=0)
    A = A + Fraction()  # force fraction if possible
    M = sympy.Matrix(A)
    M = ((M.T).echelon_form()).T
    M = np.array(M)
    basis = M[:,np.where(np.all(M[:n]==0, axis=0))[0]][n:]
    return basis


def decide_dtype(use_scipy: bool = True):
    return "float64" if use_scipy else "object"
