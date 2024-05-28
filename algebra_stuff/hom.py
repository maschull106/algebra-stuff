from dataclasses import dataclass
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

