from dataclasses import dataclass
from typing import Tuple
from .module import *
import scipy


class HomSpace:
    def __init__(self, M: ModuleFromIdeal, N: ModuleFromIdeal):
        # TODO: make some sanity checks
        self.M = M
        self.N = N
        self.ring = M.base_ring
    
    def get_matrix_representation(self, phi: Morphism):
        # for a module morphism phi: M -> N, compute the matrix representation
        matrix = [self.N.to_basis(phi(f)) for f in self.M.basis]
        matrix = np.array(matrix).T
        return matrix


def hom_complexity(M: ModuleFromIdeal, N: ModuleFromIdeal):
    ring = M.base_ring
    m = len(M.basis)
    n = len(N.basis)
    k = len(ring.basis)
    print("k =", k)
    print("m =", m)
    print("n =", n)
    print("total matrix entries:", k*m*n*m*n)
    

@ExecTimes.track_time
def hom_constraints(M: ModuleFromIdeal, N: ModuleFromIdeal) -> np.ndarray:
    # give a matrix C such that for H a matrix k^m -> k^n,
    # H represents a morphism of modules M -> N iff CH_ = 0,
    # where H_ is a vectorized form of H (Hij = H_(i*m+j))
    if not M.base_ring == N.base_ring:
        raise TypeError("modules are not defined over the same ring")
    
    ExecTimes.time_step("get matrices of M")
    ring = M.base_ring
    m = len(M.basis)
    n = len(N.basis)
    k = len(ring.basis)
    FM = M.get_matrices_representation()
    ExecTimes.time_step("get matrices of N")
    FN = N.get_matrices_representation()
    ExecTimes.time_step(f"init zero matrix of dim {n*m*k} x {n*m}")
    A = np.zeros((n*m*k, n*m))
    ExecTimes.time_step("fill the matrix with constraints")
    for f_ind in range(k):
        for i in range(n):
            for j in range(m):
                ind = f_ind*n*m + i*m + j
                A[ind, i*m: i*m+m] += FM[f_ind][:, j]
                A[ind, j: n*m+j: m] -= FN[f_ind][i]
                # for l in range(m):
                #     A[ind, i*m+l] += FM[f_ind][l, j]
                # for l in range(n):
                #     A[ind, l*m+j] -= FN[f_ind][i, l]
    A = A[np.any(A, axis=1)]    # filter out some useless constraints
    return A


@ExecTimes.track_time
def hom_rank(M: ModuleFromIdeal, N: ModuleFromIdeal, tol: float = None, char: int = 0) -> int:
    m = len(M.basis)
    n = len(N.basis)
    ExecTimes.time_step("get constraints")
    C = hom_constraints(M, N)
    if char > 0:
        C %= char
    ExecTimes.time_step("calculate rank")
    return m*n - (np.linalg.matrix_rank(C, tol=tol) if C.size > 0 else 0)


def null_space(A: np.ndarray):
    # TODO: make sure this is equivalent to just scipy.linalg.null_space(A)
    # just scipy.linalg.null_space(A) works but seems to be way slower for A with nb of rows (a lot) bigger than number of columns
    P, L, U = scipy.linalg.lu(A)
    basis = scipy.linalg.null_space(U)
    return basis


def hom(M: ModuleFromIdeal, N: ModuleFromIdeal, char: int = 0) -> np.ndarray:
    C = hom_constraints(M, N, char=char)
    basis = null_space(C)
    return basis


@dataclass
class HilbertModules:
    S: FiniteDimRing
    J: IdealQuotientModule
    O: RingQuotientModule
    
    def components(self):
        return self.S, self.J, self.O
    
    def dims(self):
        return self.S.dim, self.J.dim, self.O.dim


class HilbertTangentSpace:
    def __init__(self, I: PolyRingIdeal, char: int = 0):
        """
        I ideal of a polynomial ring R, giving a zero dimensional subscheme Z of length n.
        Represents the tangent of the Hilbert scheme of points n at the point [Z]
        """
        self.char = char
        self.R = I.base
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
    
    
@dataclass
class NestedModules:
    S: FiniteDimRing
    J1: IdealQuotientModule
    J2: IdealQuotientModule
    O1: RingQuotientModule
    O2: RingQuotientModule
    
    def components(self):
        return self.S, self.J1, self.J2, self.O1, self.O2
    
    def dims(self):
        return self.S.dim, self.J1.dim, self.J2.dim, self.O1.dim, self.O2.dim
    
    
def get_nested_modules(I1: PolyRingIdeal, I2: PolyRingIdeal, S: FiniteDimRing = None) -> NestedModules:
    # I2 subset of I1
    R = I1.base
    if S is None:
        S = R//I2   # all modules will be over this ring
    else:
        focus_base_ring(S)
    J1 = I1/I1**2   # maybe should divide by I2**2 instead to be able to define morphism J2 -> J1 (no actually)
    O1 = R/I1
    J2 = I2/I2**2
    O2 = R/I2
    return NestedModules(S=S, J1=J1, J2=J2, O1=O1, O2=O2)


def nested_hom_constraints(I1: PolyRingIdeal=None, I2: PolyRingIdeal=None, nested_modules: NestedModules=None):
    # I2 subset of I1
    if nested_modules is None:
        if I1 is None or I2 is None:
            raise ValueError("not enough information to compute hom")
        nested_modules = get_nested_modules(I1, I2)
    S, J1, J2, O1, O2 = nested_modules.components()
    k, m1, m2, n1, n2 = nested_modules.dims()
    C1 = hom_constraints(J1, O1)
    C2 = hom_constraints(J2, O2)
    # TODO:
    # matrix representation of phi: I2 -> I1 and psi: O2 -> O1
    # matrix C of n1 * m2 constraints on n1*m1 + n2*m2 variables corresponding to H∘phi - psi∘K = 0
    # final constraint matrix
    # (C1 0)
    # (0 C1)
    # (C   )
    phi = HomSpace(J2, J1).get_matrix_representation(lambda f: f)
    psi = HomSpace(O2, O1).get_matrix_representation(lambda f: f)
    C = np.zeros((n1*m2, n1*m1 + n2*m2))
    for i in range(n1):
        for j in range(m2):
            ind = i*m2 + j
            C[ind, i*m1: i*m1+m1] += phi[:, j]
            C[ind, n1*m1 + j: n1*m1 + n2*m2 + j: m2] -= psi[i]
            # for l in range(m1):
            #     C[ind, i*m1+l] += phi[l, j]
            # for l in range(n1):
            #     C[ind, n1*m1 + l*m2+j] -= psi[i, l]
    C = np.concatenate(
        [
            np.concatenate([C1, np.zeros((C1.shape[0], n2*m2))], axis=1),
            np.concatenate([np.zeros((C2.shape[0], n1*m1)), C2], axis=1),
            C
        ],
        axis=0
    )
    C = C[np.any(C, axis=1)]    # filter out some useless constraints
    return C


def nested_hom_rank(I1: PolyRingIdeal, I2: PolyRingIdeal, tol: float = None) -> int:
    nested_modules = get_nested_modules(I1, I2)
    k, m1, m2, n1, n2 = nested_modules.dims()
    C = nested_hom_constraints(nested_modules=nested_modules)
    return n1*m1 + n2*m2 - np.linalg.matrix_rank(C, tol=tol)


class YoungDiagramIdeals:
    def __init__(self, diagram: List[int], nested_ideals: List[List[PolyRingIdeal]], R: PolyRing = None):
        self.diagram = diagram
        self.nested_ideals = nested_ideals
        if R is None:
            R = infer_poly_ring()
        self.R = R
        self.sanity_check()
    
    def sanity_check(self):
        if any(self.diagram[i+1] > self.diagram[i] for i in range(len(self.diagram)-1)):
            raise ValueError
        if len(self.nested_ideals) != len(self.diagram):
            raise ValueError
        if any(len(row) != row_len for row, row_len in zip(self.nested_ideals, self.diagram)):
            return ValueError
        if any(ideal.base != self.R for row in self.nested_ideals for ideal in row):
            raise ValueError
        
        inclusion_constraint = lambda seq: all(seq[i].contains_ideal(seq[i+1]) for i in range(len(seq)-1))
        if not all(map(inclusion_constraint, self.rows())):
            raise ValueError
        if not all(map(inclusion_constraint, self.columns())):
            raise ValueError

    def __getitem__(self, i: int, j: int):
        """row i, column j"""
        return self.nested_ideals[i][j]
    
    def __iter__(self):
        for row in self.nested_ideals:
            for ideal in row:
                yield ideal
    
    def size(self) -> int:
        return sum(map(len, self.nested_ideals))
    
    def row(self, i: int) -> List[PolyRingIdeal]:
        return self.nested_ideals[i]
    
    def column(self, i: int) -> List[PolyRingIdeal]:
        column = []
        row_ind = 0
        while row_ind < len(self.diagram) and self.diagram[row_ind] >= i+1:
            column.append(self.nested_ideals[i])
        return column
    
    def rows(self) -> List[List[PolyRingIdeal]]:
        return self.nested_ideals
    
    def columns(self) -> List[List[PolyRingIdeal]]:
        return [self.column(i) for i in range(self.diagram[0])]
    
    def index_mapping(self, i: int, j: int) -> int:
        """convert index (i, j) to the corresponding index if all the rows were concatenated"""
        return sum(len(self.nested_ideals[k]) for k in range(i)) + j
    
    def row_repr(self, first: bool = False, last: bool = False) -> str:
        return ""
    
    def __repr__(self):
        return ""


class DoubleNestedHilbertScheme:
    def __init__(self, diagram: List[int], R: PolyRing = None):
        """
        diagram: a representation of the Young diagram that defines the structure of the space
        For example, the list [3, 3, 2, 1, 1] corresponds to the diagram below.
        ┌─┬─┬─┐
        ├─┼─┼─┤
        ├─┼─┼─┘
        ├─┼─┘
        ├─┤
        └─┘
        """
        if R is None:
            R = infer_poly_ring()
        self.R = R
        self.diagram = diagram
    
    def tangent_space(self, nested_ideals: List[List[PolyRingIdeal]]):
        return DoubleNestedHilbertSchemeTangentSpace(self, nested_ideals)


class DoubleNestedHilbertSchemeTangentSpace:
    def __init__(self, base: DoubleNestedHilbertScheme, nested_ideals: List[List[PolyRingIdeal]]):
        """
        The shape of the nested_ideals list must coincide with the space's diagram
        """
        self.base = base
        self.diagram_ideals = YoungDiagramIdeals(self.base.diagram, nested_ideals)
        self.constraint_sizes: List[int] = []
        self.constraints = self._compute_constraints()
    
    @staticmethod
    def _nested_hom_constraints(I1: PolyRingIdeal = None, I2: PolyRingIdeal = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Required: I₂ subset of I₁
        Constraints for the morphisms I₁ -> O₁ and I₂ -> O₂ to respect the inclusion Z₁ -> Z₂,
        that is commutativity of the following diagram
        I₁ → I₂
        ↓    ↓
        O₁ → O₂
        
        Returns the constraint in two separate matrices, to be glued later
        """
        nested_modules = get_nested_modules(I1, I2)
        S, J1, J2, O1, O2 = nested_modules.components()
        k, m1, m2, n1, n2 = nested_modules.dims()
        phi = HomSpace(J2, J1).get_matrix_representation(lambda f: f)
        psi = HomSpace(O2, O1).get_matrix_representation(lambda f: f)
        C1 = np.zeros((n1*m2, n1*m1))
        C2 = np.zeros((n1*m2, n2*m2))
        for i in range(n1):
            for j in range(m2):
                ind = i*m2 + j
                C1[ind, i*m1: i*m1+m1] += phi[:, j]
                C2[ind, j: n2*m2 + j: m2] -= psi[i]
        return C1, C2

    def _zeros_before(self, ind):
        return sum(self.constraint_sizes[:ind])
    
    def _zeros_after(self, ind):
        return sum(self.constraint_sizes[ind+1:])
    
    def _zeros_between(self, ind1, ind2):
        return sum(self.constraint_sizes[ind1+1:ind2])
    
    def _morphism_constraints(self) -> np.ndarray:
        """
        Constraints for every map I -> O to be a morphism of R-modules.
        Also updates the constraint_sizes attribute
        """
        morphism_constraints = []
        for I in self.diagram_ideals:
            C = HilbertTangentSpace(I).constraints()
            morphism_constraints.append(C)
            size = C.shape[1]
            self.constraint_sizes.append(size)
        
        # padding with zeros
        for ind, C in enumerate(morphism_constraints):
            constraint_count = C.shape[0]
            padding = lambda dim: np.zeros((constraint_count, dim))
            morphism_constraints[ind] = np.concatenate([padding(self._zeros_before(ind)), C, padding(self._zeros_after(ind))], axis=1)
        C = np.concatenate(morphism_constraints, axis=0)
        return C
    
    def _kernel_constraint(self, pos1: Tuple[int, int], pos2: Tuple[int, int]) -> np.ndarray:
        """
        Constraint for the morphisms corresponding I₁ -> O₁ and I₂ -> O₂ to respect the inclusion Z₁ -> Z₂
        Basically a wrapper around _nested_hom_constraints glueing the constraints appropriately
        """
        get_ideal_ind = lambda i, j: (self.diagram_ideals.index_mapping(i, j), self.diagram_ideals[i, j])
        ind1, I1 = get_ideal_ind(*pos1)
        ind2, I2 = get_ideal_ind(*pos2)
        C1, C2 = self._nested_hom_constraints(I1, I2)
        constraint_count = C1.shape[0]
        
        # padding with zeros
        padding = lambda dim: np.zeros((constraint_count, dim))
        C = np.concatenate(
            [
                padding(self._zeros_before(ind1)),
                C1,
                padding(self._zeros_between(ind1, ind2)),
                C2,
                padding(self._zeros_after(ind2))
            ],
            axis=1
        )
        filter_zero = lambda C: C[np.any(C, axis=1)]    # filter out some useless constraints
        C = filter_zero(C)
        return C

    def _compute_constraints(self) -> np.ndarray:
        morphism_constraints = self._morphism_constraints()

        row_constraints = []
        for i, row in enumerate(self.diagram_ideals.rows()):
            for j in range(len(row)-1):
                C = self._kernel_constraint((i, j), (i, j+1))
                row_constraints.append(C)
        
        column_constraints = []
        for j, column in enumerate(self.diagram_ideals.columns()):
            for i in range(len(column)-1):
                C = self._kernel_constraint((i, j), (i+1, j))
                column_constraints.append(C)
        
        C = np.concatenate(
            [
                morphism_constraints,
                row_constraints,
                column_constraints
            ],
            axis = 0
        )
        return C
    
    def dim(self) -> int:
        max_rank = sum(self.constraint_sizes)
        return max_rank - (np.linalg.matrix_rank(self.constraints) if self.constraints.size > 0 else 0)
    
    def basis(self) -> np.ndarray:
        return null_space(self.constraints)
