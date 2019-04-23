import numpy as np


class ScalarComplexFunction:
    def __init__(self, f):
        self.f = f

    def __call__(self, x):
        x_real, x_imag = np.split(x, 2)
        cx = x_real + 1j * x_imag
        return self.f(cx)


class VectorComplexFunction:
    def __init__(self, f):
        self.f = f

    def __call__(self, x):
        x_real, x_imag = np.split(x, 2)
        cx = x_real + 1j * x_imag
        y = self.f(cx)
        return np.r_[y.real, y.imag]


class Residual:
    def __init__(self, a, b):
        self.a = a
        self.b = b
        self.aT_dot_b = a.rmatvec(b)
        self.x0 = np.zeros(self.a.shape[1], dtype=complex)
        self.a_dot_x = np.zeros(self.a.shape[0], dtype=complex)

    def f(self, vec):
        if np.linalg.norm(vec - self.x0) != 0.:
            self.a_dot_x = self.a(vec)
            self.x0 = vec
        residual = self.a_dot_x - self.b
        return np.vdot(residual, residual).real

    def fp(self, vec):
        if np.linalg.norm(vec - self.x0) != 0.:
            self.a_dot_x = self.a(vec)
            self.x0 = vec
        return 2 * (self.a.rmatvec(self.a_dot_x) - self.aT_dot_b)

    # def hessp(self, p):
    #     return 2 * self.a.rmatvec(self.a(p))


class Quadratic:
    def __init__(self, c, alpha):
        self.c = c
        self.alpha = alpha
        self.x0 = np.zeros(self.c.shape[1], dtype=complex)
        self.quad = 0.

    def calc_quad(self, x):
        return np.dot(x.real, self.c(x.real)).real + np.dot(x.imag, self.c(x.imag)).real

    def f(self, x):
        if np.linalg.norm(self.x0 - x) != 0.:
            self.quad = self.calc_quad(x)
            self.x0 = x
        return (self.quad - self.alpha) ** 2

    def fp(self, x):
        if np.linalg.norm(self.x0 - x) != 0.:
            self.quad = self.calc_quad(x)
            self.x0 = x
        return 4 * (self.quad - self.alpha) * self.c(x)

    # def hessp(self, x, p):
    #     cx = self.c(x)
    #     return 8. * np.dot(np.outer(cx, cx), p) \
    #            + 4. * (np.vdot(x, self.c(x)) - self.alpha) * self.c(p)


class PenaltyFunction:
    def __init__(self, a, b, c, alpha, mu, mu0=1., mu1=1.):
        self.residual = Residual(a, b)
        self.quadratic = Quadratic(c, alpha)
        self.mu = mu
        self.mu0 = mu0
        self.mu1 = mu1

    def f(self, x):
        f0 = self.mu0 * self.residual.f(x)
        f1 = self.mu1 * self.quadratic.f(x)
        print(f0, f1)
        return f0 + f1

    def fp(self, x):
        return self.mu0 * self.residual.fp(x) + self.mu1 * self.quadratic.fp(x)

    # def hessp(self, x, p):
    #     return self.residual.hessp(p) + self.quadratic.hessp(x, p) / self.mu


class QuadraticError:
    def __init__(self, c, alpha, tol):
        self.c = c
        self.alpha_lower = (1 - tol) * alpha
        self.alpha_upper = (1 + tol) * alpha
        self.x0 = np.zeros(self.c.shape[1], dtype=complex)
        self.quad = 0.

    def calc_quad(self, x):
        return np.dot(x.real, self.c(x.real)).real + np.dot(x.imag, self.c(x.imag)).real

    def f(self, x):
        if np.linalg.norm(self.x0 - x) != 0.:
            self.quad = self.calc_quad(x)
            self.x0 = x
        quad_half = -self.quad * .5
        if quad_half < self.alpha_lower:
            y = (quad_half - self.alpha_lower) ** 2
        elif self.alpha_upper < quad_half:
            y = (quad_half - self.alpha_upper) ** 2
        else:
            y = 0.
        return y

    def fp(self, x):
        if np.linalg.norm(self.x0 - x) != 0.:
            self.quad = self.calc_quad(x)
            self.x0 = x
        quad_half = -self.quad * .5
        if quad_half < self.alpha_lower:
            yp = (self.quad + self.alpha_lower * 2.) * self.c(x)
        elif self.alpha_upper < quad_half:
            yp = (self.quad + self.alpha_upper * 2.) * self.c(x)
        else:
            yp = np.zeros_like(self.x0)
        return yp


class PenaltyFunctionError:
    def __init__(self, a, b, c, alpha, mu, mu0=1., mu1=1., tol=0.):
        self.residual = Residual(a, b)
        self.quadratic = QuadraticError(c, alpha, tol)
        self.mu = mu
        self.mu0 = mu0
        self.mu1 = mu1

    def f(self, x):
        f0 = self.mu0 * self.residual.f(x)
        f1 = self.mu1 * self.quadratic.f(x)
        print(f0, f1)
        return f0 + f1

    def fp(self, x):
        fp0 = self.mu0 * self.residual.fp(x)
        fp1 = self.mu1 * self.quadratic.fp(x)
        # print(fp0, fp1)
        return fp0 + fp1

    # def hessp(self, x, p):
    #     return self.residual.hessp(p) + self.quadratic.hessp(x, p) / self.mu
