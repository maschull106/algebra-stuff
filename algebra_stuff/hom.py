from __future__ import annotations
from dataclasses import dataclass
from .module import *
import scipy
import scipy.linalg


class HomSpace:
    def __init__(self, M: Module, N: Module, base: FiniteDimRing = None, precompute_constraints: bool = True, use_scipy: bool = True):
        # TODO: make some sanity checks (or try to convert modules to be over the given base)
        self.M = M
        self.N = N
        self.use_scipy = use_scipy #or Scalar.MODE != Scalar.FRACTION
        self.dtype = decide_dtype(self.use_scipy)
        if base is None:
            base = M.base_ring
        self.ring = base
        if self.M.base_ring != self.ring:
            self.M.set_base(self.ring)
        if self.N.base_ring != self.ring:
            self.N.set_base(self.ring)
        
        self.constraints_computed = False
        if precompute_constraints:
            self.constraints = self._compute_constraints()
        else:
            self.constraints = np.zeros((1, 1), dtype=self.dtype)
    
    def domain(self) -> Module:
        return self.M
    
    def codomain(self) -> Module:
        return self.N
    
    @ExecTimes.track_time
    def _compute_constraints(self) -> np.ndarray:
        # give a matrix C such that for H a matrix k^m -> k^n,
        # H represents a morphism of modules M -> N iff CH_ = 0,
        # where H_ is a vectorized form of H (Hij = H_(i*m+j))
        if not self.M.base_ring == self.N.base_ring:
            raise TypeError("modules are not defined over the same ring")
        
        ExecTimes.time_step("get matrices of M")
        m = self.M.dim
        n = self.N.dim
        k = self.ring.dim
        FM = self.M.get_matrices_representation(dtype=self.dtype)
        ExecTimes.time_step("get matrices of N")
        FN = self.N.get_matrices_representation(dtype=self.dtype)
        ExecTimes.time_step(f"init zero matrix of dim {n*m*k} x {n*m}")
        C = np.zeros((n*m*k, n*m), dtype=self.dtype)
        CCC = np.copy(C)
        ExecTimes.time_step("fill the matrix with constraints")
        for f_ind in range(k):
            for i in range(n):
                # TODO: try to get rid of one more for loop via numpy
                for j in range(m):
                    ind = f_ind*n*m + i*m + j
                    C[ind, i*m: i*m+m] += FM[f_ind][:, j]
                    C[ind, j: n*m+j: m] -= FN[f_ind][i]
        #C = C[np.any(C, axis=1)]    # filter out some useless constraints
        C = filter_zero(C)
        self.constraints_computed = True
        return C
    
    def _get_constraints(self) -> np.ndarray:
        if self.constraints_computed:
            return self.constraints
        self.constraints = self._compute_constraints()
        return self.constraints
    
    def dim(self) -> int:
        m = self.M.dim
        n = self.N.dim
        C = self._get_constraints()
        C = C if self.use_scipy else C.astype("float64")
        rank = (np.linalg.matrix_rank(C) if C.size > 0 else 0)
        return m*n - rank
    
    def basis(self) -> np.ndarray:
        C = self._get_constraints()
        if self.use_scipy:
            basis = null_space(C)
        else:
            basis = exact_null_space(C)
            if Scalar.MODE != Scalar.FRACTION:
                basis = basis.astype("float64")
        return basis
    
    def basis_as_matrices(self) -> List[np.ndarray]:
        basis = self.basis()
        m = self.M.dim
        n = self.N.dim
        basis_list = []
        for vect in basis.T:
            vect = vect.reshape((n, m))
            basis_list.append(vect)
        return basis_list
    
    def basis_as_morphisms(self) -> List[MorphismFromMatrix]:
        basis = self.basis_as_matrices()
        return [MorphismFromMatrix(self.M, self.N, matrix) for matrix in basis]
    
    def apply_morphism(self, phi: np.ndarray, f: GroebnerPolynomial):
        """apply a morphism phi given in matrix form to the element f"""
        vect = self.M.to_basis(f)
        vect = np.array(vect)
        res = phi @ vect
        return self.N.from_basis(res)
    
    def get_matrix_representation(self, phi: Morphism):
        """for a module morphism phi: M -> N, compute the matrix representation"""
        matrix = [self.N.to_basis(phi(f)) for f in self.M.basis]
        matrix = np.array(matrix, dtype=self.dtype).T
        return matrix


class MorphismFromMatrix:
    def __init__(self, M: Module, N: Module, matrix: np.ndarray):
        self.M = M
        self.N = N
        self.matrix = matrix
    
    def domain(self) -> Module:
        return self.M
    
    def codomain(self) -> Module:
        return self.N
    
    def __call__(self, f: GroebnerPolynomial):
        if not self.M.contains(f):
            raise ValueError
        vect = self.M.to_basis(f)
        vect = np.array(vect)
        res = self.matrix @ vect
        return self.N.from_basis(res)


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


def hom(M: ModuleFromIdeal, N: ModuleFromIdeal, char: int = 0) -> np.ndarray:
    C = hom_constraints(M, N, char=char)
    basis = null_space(C)
    return basis
    
    
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
    phi = HomSpace(J2, J1, precompute_constraints=False).get_matrix_representation(lambda f: f)
    psi = HomSpace(O2, O1, precompute_constraints=False).get_matrix_representation(lambda f: f)
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

