import numpy as np
import time
from scipy import optimize
from scipy import linalg

# A = np.array([[1.j, 1, -2], [-3, 7, -3], [3, -5, 5]])
A = np.load("mat/Z_PMCHWT.npy")
Ar = A.conj().T

b = np.random.rand(A.shape[1]) + 1j * np.random.rand(A.shape[1])
Ar_dot_b = np.dot(Ar, b)
norm_Ar_dot_b = np.linalg.norm(Ar_dot_b)


def f(vec):
    vec_real, vec_imag = np.split(vec, 2)
    vec_complex = vec_real + 1j * vec_imag
    resi = np.dot(A, vec_complex) - b
    return (np.vdot(resi, resi)).real * 0.5


def gradf(vec):
    vec_real, vec_imag = np.split(vec, 2)
    Ar_dot_A_dot_vec_real = np.dot(Ar, np.dot(A, vec_real))
    Ar_dot_A_dot_vec_imag = np.dot(Ar, np.dot(A, vec_imag))
    gradR_f = np.real(Ar_dot_A_dot_vec_real) - np.imag(Ar_dot_A_dot_vec_imag) \
              - np.real(Ar_dot_b)
    gradI_f = np.real(Ar_dot_A_dot_vec_imag) + np.imag(Ar_dot_A_dot_vec_real) \
              - np.imag(Ar_dot_b)
    return np.r_[gradR_f, gradI_f]


def gradf2(vec):
    vec_real, vec_imag = np.split(vec, 2)
    vec_complex = vec_real + 1j * vec_imag
    resi2 = np.dot(Ar, np.dot(A, vec_complex)) - Ar_dot_b
    return np.r_[resi2.real, resi2.imag]


x0 = np.zeros(A.shape[1] * 2)
res = optimize.minimize(f, x0, method="CG", jac=gradf2,
                        options={"gtol": norm_Ar_dot_b * 1e-3, "norm": 2})
x_real, x_imag = np.split(res.x, 2)
x_optimized = x_real + 1j * x_imag
# x = optimize.fmin_cg(f, x0, fprime=gradf, maxiter=10)
# x_real, x_imag = np.split(x, 2)
# x_optimized = x_real + 1j*x_imag

x_solved = linalg.solve(A, b)
x0 = np.r_[x_solved.real, x_solved.imag]
print("check_grad =", optimize.check_grad(f, gradf, x0))

# a = np.load("mat/Z_PMCHWT.npy")
# b = np.ones(a.shape[0], dtype=complex)
# x = linalg.solve(a, b)
#
# bnorm2 = np.real(np.vdot(b, b))
# f = lambda vec: np.real(np.dot(np.dot(a, vec) - b, (np.dot(a, vec) - b).conj())) / np.real(bnorm2)
# def f(vec):
#     vec_complex = vec[:vec.size/2] + 1j*vec[vec.size/2:]
#     resi = np.dot(a, vec_complex) - b
#     return np.real(np.vdot(resi, resi)) / bnorm2
# a_r = a.conj().T
# gradf = lambda vec: np.dot(a_r, np.dot(a, vec) - b)
# res = optimize.minimize(f, np.zeros(a.shape[1]*2))
# x_optimized = res.x[:res.x.size/2] + 1j*res.x[res.x.size/2:]

print("Complete!")
