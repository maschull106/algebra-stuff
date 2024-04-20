from .module import *


def hom_complexity(M: ModuleFromIdeal, N: ModuleFromIdeal):
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
    ExecTimes.track_time(hom_constraints)
    return A


def hom_rank(M: ModuleFromIdeal, N: ModuleFromIdeal, tol: float = None) -> int:
    m = len(M.basis)
    n = len(N.basis)
    C = hom_constraints(M, N)
    return m*n - np.linalg.matrix_rank(C, tol=tol)


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
