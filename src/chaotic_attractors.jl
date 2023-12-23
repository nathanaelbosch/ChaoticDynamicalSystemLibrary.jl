ATTRACTOR_DATA = JSON.parsefile(joinpath(@__DIR__, "chaotic_attractors.json"))

function Lorenz end
originalcode(::typeof(Lorenz)) = """
class Lorenz(DynSys):
    @staticjit
    def _rhs(x, y, z, t, beta, rho, sigma):
        xdot = sigma * y - sigma * x
        ydot = rho * x - x * z - y
        zdot = x * y - beta * z
        return xdot, ydot, zdot
    @staticjit
    def _jac(x, y, z, t, beta, rho, sigma):
        row1 = [-sigma, sigma, 0]
        row2 = [rho - z, -1, -x]
        row3 = [y, x, -beta]
        return [row1, row2, row3]
"""
@doc make_docstring(Lorenz) Lorenz
function Lorenz()
    function rhs(du, u, p, t)
        @unpack beta, rho, sigma = p
        du[1] = sigma * (u[2] - u[1])
        du[2] = u[1] * (rho - u[3]) - u[2]
        du[3] = u[1] * u[2] - beta * u[3]
    end
    function jac(J, u, p, t)
        @unpack beta, rho, sigma = p
        J[1, 1] = -sigma
        J[1, 2] = sigma
        J[1, 3] = 0
        J[2, 1] = rho - u[3]
        J[2, 2] = -1
        J[2, 3] = -u[1]
        J[3, 1] = u[2]
        J[3, 2] = u[1]
        J[3, 3] = -beta
    end
    u0 = Float64.(ATTRACTOR_DATA["Lorenz"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Lorenz"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs, jac = jac)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function LorenzBounded end
function originalcode(::typeof(LorenzBounded))
    """
class LorenzBounded(DynSys):
    @staticjit
    def _rhs(x, y, z, t, beta, r, rho, sigma):
        xdot = sigma * y - sigma * x - sigma/r**2 * y * x ** 2 - sigma/r**2 * y ** 3 - sigma/r**2 * y * z ** 2 + sigma/r**2 * x ** 3 + sigma/r**2 * x * y ** 2 + sigma/r**2 * x * z ** 2
        ydot = rho * x - x * z - y - rho/r**2 * x ** 3 - rho/r**2 * x * y ** 2 - rho/r**2 * x * z ** 2 + 1/r**2 * z * x ** 3 + 1/r**2 * x * z * y ** 2 + 1/r**2 * x * z ** 3 + 1/r**2 * y * x ** 2 + 1/r**2 * y ** 3 + 1/r**2 * y * z ** 2
        zdot = x * y - beta * z - 1/r**2 * y * x ** 3 - 1/r**2 * x * y ** 3 - 1/r**2 * x * y * z ** 2 + beta/r**2 * z * x ** 2 + beta/r**2 * z * y ** 2 + beta/r**2 * z ** 3
        return xdot, ydot, zdot
"""
end
@doc make_docstring(LorenzBounded) LorenzBounded
function LorenzBounded()
    function rhs(du, u, p, t)
        @unpack beta, r, rho, sigma = p
        x, y, z = u
        du[1] = sigma * y - sigma * x - sigma / r^2 * y * x^2 - sigma / r^2 * y^3 -
                sigma / r^2 * y * z^2 + sigma / r^2 * x^3 + sigma / r^2 * x * y^2 +
                sigma / r^2 * x * z^2
        du[2] = rho * x - x * z - y - rho / r^2 * x^3 - rho / r^2 * x * y^2 -
                rho / r^2 * x * z^2 + 1 / r^2 * z * x^3 + 1 / r^2 * x * z * y^2 +
                1 / r^2 * x * z^3 + 1 / r^2 * y * x^2 + 1 / r^2 * y^3 + 1 / r^2 * y * z^2
        du[3] = x * y - beta * z - 1 / r^2 * y * x^3 - 1 / r^2 * x * y^3 -
                1 / r^2 * x * y * z^2 + beta / r^2 * z * x^2 + beta / r^2 * z * y^2 +
                beta / r^2 * z^3
    end
    u0 = Float64.(ATTRACTOR_DATA["LorenzBounded"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["LorenzBounded"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function LorenzCoupled end
originalcode(::typeof(LorenzCoupled)) = """
class LorenzCoupled(DynSys):
    @staticjit
    def _rhs(x1, y1, z1, x2, y2, z2, t, beta, eps, rho, rho1, rho2, sigma):
        x1dot = sigma * y1 - sigma * x1
        y1dot = rho1 * x1 - x1 * z1 - y1
        z1dot = x1 * y1 - beta * z1
        x2dot = sigma * y2 - sigma * x2 + eps * x1 - eps * x2
        y2dot = rho2 * x2 - x2 * z2 - y2
        z2dot = x2 * y2 - beta * z2
        return x1dot, y1dot, z1dot, x2dot, y2dot, z2dot
"""
@doc make_docstring(LorenzCoupled) LorenzCoupled
function LorenzCoupled()
    function rhs(du, u, p, t)
        @unpack beta, eps, rho, rho1, rho2, sigma = p
        du[1] = sigma * (u[2] - u[1])
        du[2] = rho1 * u[1] - u[1] * u[3] - u[2]
        du[3] = u[1] * u[2] - beta * u[3]
        du[4] = sigma * (u[5] - u[4]) + eps * u[1] - eps * u[4]
        du[5] = rho2 * u[4] - u[4] * u[6] - u[5]
        du[6] = u[4] * u[5] - beta * u[6]
    end
    u0 = Float64.(ATTRACTOR_DATA["LorenzCoupled"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["LorenzCoupled"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Lorenz96 end
originalcode(::typeof(Lorenz96)) = """
class Lorenz96(DynSys):
    def rhs(self, X, t):
        Xdot = np.zeros_like(X)
        Xdot[0] = (X[1] - X[-2]) * X[-1] - X[0] + self.f
        Xdot[1] = (X[2] - X[-1]) * X[0] - X[1] + self.f
        Xdot[-1] = (X[0] - X[-3]) * X[-2] - X[-1] + self.f
        Xdot[2:-1] = (X[3:] - X[:-3]) * X[1:-2] - X[2:-1] + self.f
        return Xdot
"""
@doc make_docstring(Lorenz96) Lorenz96
function Lorenz96()
    function rhs(du, u, p, t)
        f = p[1]
        du[1] = (u[2] - u[end - 1]) * u[end] - u[1] + f
        du[2] = (u[3] - u[end]) * u[1] - u[2] + f
        du[end] = (u[1] - u[end - 2]) * u[end - 1] - u[end] + f
        @simd ivdep for i in 3:(length(u) - 1)
            du[i] = (u[i + 1] - u[i - 2]) * u[i - 1] - u[i] + f
        end
    end
    u0 = Float64.(ATTRACTOR_DATA["Lorenz96"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Lorenz96"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Lorenz84 end
originalcode(::typeof(Lorenz84)) = """
class Lorenz84(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, f, g):
        xdot = -a * x - y ** 2 - z ** 2 + a * f
        ydot = -y + x * y - b * x * z + g
        zdot = -z + b * x * y + x * z
        return xdot, ydot, zdot
"""
@doc make_docstring(Lorenz84) Lorenz84
function Lorenz84()
    function rhs(du, u, p, t)
        @unpack a, b, f, g = p
        du[1] = -a * u[1] - u[2]^2 - u[3]^2 + a * f
        du[2] = -u[2] + u[1] * u[2] - b * u[1] * u[3] + g
        du[3] = -u[3] + b * u[1] * u[2] + u[1] * u[3]
    end
    u0 = Float64.(ATTRACTOR_DATA["Lorenz84"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Lorenz84"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Rossler end
originalcode(::typeof(Rossler)) = """
class Rossler(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c):
        xdot = -y - z
        ydot = x + a * y
        zdot = b + z * x - c * z
        return xdot, ydot, zdot
"""
@doc make_docstring(Rossler) Rossler
function Rossler()
    function rhs(du, u, p, t)
        @unpack a, b, c = p
        du[1] = -u[2] - u[3]
        du[2] = u[1] + a * u[2]
        du[3] = b + u[3] * u[1] - c * u[3]
    end
    u0 = Float64.(ATTRACTOR_DATA["Rossler"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Rossler"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Thomas end
originalcode(::typeof(Thomas)) = """
class Thomas(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b):
        xdot = -a * x + b * np.sin(y)
        ydot = -a * y + b * np.sin(z)
        zdot = -a * z + b * np.sin(x)
        return xdot, ydot, zdot
"""
@doc make_docstring(Thomas) Thomas
function Thomas()
    function rhs(du, u, p, t)
        @unpack a, b = p
        du[1] = -a * u[1] + b * sin(u[2])
        du[2] = -a * u[2] + b * sin(u[3])
        du[3] = -a * u[3] + b * sin(u[1])
    end
    u0 = Float64.(ATTRACTOR_DATA["Thomas"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Thomas"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function DoublePendulum end
originalcode(::typeof(DoublePendulum)) = """
class DoublePendulum(DynSys):
    @staticjit
    def _rhs(th1, th2, p1, p2, t, d, m):
        g = 9.82
        pre = 6 / (m * d ** 2)
        denom = 16 - 9 * np.cos(th1 - th2) ** 2
        th1_dot = pre * (2 * p1 - 3 * np.cos(th1 - th2) * p2) / denom
        th2_dot = pre * (8 * p2 - 3 * np.cos(th1 - th2) * p1) / denom
        p1_dot = (
            -0.5
            * (m * d ** 2)
            * (th1_dot * th2_dot * np.sin(th1 - th2) + 3 * (g / d) * np.sin(th1))
        )
        p2_dot = (
            -0.5
            * (m * d ** 2)
            * (-th1_dot * th2_dot * np.sin(th1 - th2) + 3 * (g / d) * np.sin(th2))
        )
        return th1_dot, th2_dot, p1_dot, p2_dot
"""
@doc make_docstring(DoublePendulum) DoublePendulum
function DoublePendulum()
    function rhs(du, u, p, t)
        @unpack d, m = p
        g = 9.82
        pre = 6 / (m * d^2)
        denom = 16 - 9 * cos(u[1] - u[2])^2
        du[1] = pre * (2 * u[3] - 3 * cos(u[1] - u[2]) * u[4]) / denom
        du[2] = pre * (8 * u[4] - 3 * cos(u[1] - u[2]) * u[3]) / denom
        du[3] = -0.5 * (m * d^2) *
                (du[1] * du[2] * sin(u[1] - u[2]) + 3 * (g / d) * sin(u[1]))
        du[4] = -0.5 * (m * d^2) *
                (-du[1] * du[2] * sin(u[1] - u[2]) + 3 * (g / d) * sin(u[2]))
    end
    u0 = Float64.(ATTRACTOR_DATA["DoublePendulum"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["DoublePendulum"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SwingingAtwood end
originalcode(::typeof(SwingingAtwood)) = """
class SwingingAtwood(DynSys):
    @staticjit
    def _rhs(r, th, pr, pth, t, m1, m2):
        g = 9.82
        rdot = pr / (m1 + m2)
        thdot = pth / (m1 * r ** 2)
        prdot = pth ** 2 / (m1 * r ** 3) - m2 * g + m1 * g * np.cos(th)
        pthdot = -m1 * g * r * np.sin(th)
        return rdot, thdot, prdot, pthdot

    @staticjit
    def _postprocessing(r, th, pr, pth):
        return r, np.sin(th), pr, pth
"""
@doc make_docstring(SwingingAtwood) SwingingAtwood
function SwingingAtwood()
    function rhs(du, u, p, t)
        @unpack m1, m2 = p
        g = 9.82
        du[1] = u[3] / (m1 + m2)
        du[2] = u[4] / (m1 * u[1]^2)
        du[3] = u[4]^2 / (m1 * u[1]^3) - m2 * g + m1 * g * cos(u[2])
        du[4] = -m1 * g * u[1] * sin(u[2])
    end
    u0 = Float64.(ATTRACTOR_DATA["SwingingAtwood"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SwingingAtwood"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function GuckenheimerHolmes end
originalcode(::typeof(GuckenheimerHolmes)) = """
class GuckenheimerHolmes(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c, d, e, f):
        xdot = a * x - b * y + c * z * x + d * z * x ** 2 + d * z * y ** 2
        ydot = a * y + b * x + c * z * y
        zdot = e - z ** 2 - f * x ** 2 - f * y ** 2 - a * z ** 3
        return xdot, ydot, zdot
"""
@doc make_docstring(GuckenheimerHolmes) GuckenheimerHolmes
function GuckenheimerHolmes()
    function rhs(du, u, p, t)
        @unpack a, b, c, d, e, f = p
        du[1] = a * u[1] - b * u[2] + c * u[3] * u[1] + d * u[3] * u[1]^2 +
                d * u[3] * u[2]^2
        du[2] = a * u[2] + b * u[1] + c * u[3] * u[2]
        du[3] = e - u[3]^2 - f * u[1]^2 - f * u[2]^2 - a * u[3]^3
    end
    u0 = Float64.(ATTRACTOR_DATA["GuckenheimerHolmes"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["GuckenheimerHolmes"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HenonHeiles end
originalcode(::typeof(HenonHeiles)) = """
class HenonHeiles(DynSys):
    @staticjit
    def _rhs(x, y, px, py, t, lam):
        xdot = px
        ydot = py
        pxdot = -x - 2 * lam * x * y
        pydot = -y - lam * x ** 2 + lam * y ** 2
        return xdot, ydot, pxdot, pydot
"""
@doc make_docstring(HenonHeiles) HenonHeiles
function HenonHeiles()
    function rhs(du, u, p, t)
        lam = p[1]
        du[1] = u[3]
        du[2] = u[4]
        du[3] = -u[1] - 2 * lam * u[1] * u[2]
        du[4] = -u[2] - lam * u[1]^2 + lam * u[2]^2
    end
    u0 = Float64.(ATTRACTOR_DATA["HenonHeiles"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HenonHeiles"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Halvorsen end
originalcode(::typeof(Halvorsen)) = """
class Halvorsen(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b):
        xdot = -a * x - b * y - b * z - y ** 2
        ydot = -a * y - b * z - b * x - z ** 2
        zdot = -a * z - b * x - b * y - x ** 2
        return xdot, ydot, zdot
"""
@doc make_docstring(Halvorsen) Halvorsen
function Halvorsen()
    function rhs(du, u, p, t)
        @unpack a, b = p
        du[1] = -a * u[1] - b * u[2] - b * u[3] - u[2]^2
        du[2] = -a * u[2] - b * u[3] - b * u[1] - u[3]^2
        du[3] = -a * u[3] - b * u[1] - b * u[2] - u[1]^2
    end
    u0 = Float64.(ATTRACTOR_DATA["Halvorsen"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Halvorsen"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Chua end
originalcode(::typeof(Chua)) = """
class Chua(DynSys):
    @staticjit
    def _rhs(x, y, z, t, alpha, beta, m0, m1):
        ramp_x = m1 * x + 0.5 * (m0 - m1) * (np.abs(x + 1) - np.abs(x - 1))
        xdot = alpha * (y - x - ramp_x)
        ydot = x - y + z
        zdot = -beta * y
        return xdot, ydot, zdot
"""
@doc make_docstring(Chua) Chua
function Chua()
    function rhs(du, u, p, t)
        @unpack alpha, beta, m0, m1 = p
        ramp_x = m1 * u[1] + 0.5 * (m0 - m1) * (abs(u[1] + 1) - abs(u[1] - 1))
        du[1] = alpha * (u[2] - u[1] - ramp_x)
        du[2] = u[1] - u[2] + u[3]
        du[3] = -beta * u[2]
    end
    u0 = Float64.(ATTRACTOR_DATA["Chua"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Chua"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function MultiChua end
originalcode(::typeof(MultiChua)) = """
class MultiChua(DynSys):
    def diode(self, x):
        m, c = self.m, self.c
        total = m[-1] * x
        for i in range(1, 6):
            total += 0.5 * (m[i - 1] - m[i]) * (np.abs(x + c[i]) - np.abs(x - c[i]))
        return total

    def rhs(self, X, t):
        x, y, z = X
        xdot = self.a * (y - self.diode(x))
        ydot = x - y + z
        zdot = -self.b * y
        return (xdot, ydot, zdot)
"""
@doc make_docstring(MultiChua) MultiChua
function MultiChua()
    function diode(x, m, c)
        total = m[end] * x
        @simd ivdep for i in 1:5
            total += 0.5 * (m[i] - m[i + 1]) * (abs(x + c[i]) - abs(x - c[i]))
        end
        return total
    end
    function rhs(du, u, p, t)
        @unpack a, b, m, c = p
        du[1] = a * (u[2] - diode(u[1], m, c))
        du[2] = u[1] - u[2] + u[3]
        du[3] = -b * u[2]
    end
    u0 = Float64.(ATTRACTOR_DATA["MultiChua"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["MultiChua"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Duffing end
originalcode(::typeof(Duffing)) = """
class Duffing(DynSys):
    @staticjit
    def _rhs(x, y, z, t, alpha, beta, delta, gamma, omega):
        xdot = y
        ydot = -delta * y - beta * x - alpha * x ** 3 + gamma * np.cos(z)
        zdot = omega
        return xdot, ydot, zdot

    @staticjit
    def _postprocessing(x, y, z):
        return x, y, np.cos(z)
"""
@doc make_docstring(Duffing) Duffing
function Duffing()
    function rhs(du, u, p, t)
        @unpack alpha, beta, delta, gamma, omega = p
        du[1] = u[2]
        du[2] = -delta * u[2] - beta * u[1] - alpha * u[1]^3 + gamma * cos(u[3])
        du[3] = omega
    end
    u0 = Float64.(ATTRACTOR_DATA["Duffing"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Duffing"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

# class MackeyGlass(DynSysDelay):
#     @staticjit
#     def _rhs(x, xt, t, beta, gamma, n, tau):
#         xdot = beta * (xt / (1 + xt ** n)) - gamma * x
#         return xdot

# class IkedaDelay(DynSysDelay):
#     @staticjit
#     def _rhs(x, xt, t, c, mu, tau, x0):
#         xdot = mu * np.sin(xt - x0) - c * x
#         return xdot

# class SprottDelay(IkedaDelay):
#     pass

# class VossDelay(DynSysDelay):
#     @staticjit
#     def _rhs(x, xt, t, alpha, tau):
#         f = -10.44 * xt ** 3 - 13.95 * xt ** 2 - 3.63 * xt + 0.85
#         xdot = -alpha * x + f
#         return xdot

# class ScrollDelay(DynSysDelay):
#     @staticjit
#     def _rhs(x, xt, t, alpha, beta, tau):
#         f = np.tanh(10 * xt)
#         xdot = -alpha * xt + beta * f
#         return xdot

# class PiecewiseCircuit(DynSysDelay):
#     @staticjit
#     def _rhs(x, xt, t, alpha, beta, c, tau):
#         f = -((xt / c) ** 3) + 3 * xt / c
#         xdot = -alpha * xt + beta * f
#         return xdot

# # ## this was not chaotic
# # class ENSODelay(DynSysDelay):
# #     @staticjit
# #     def _rhs(x, xt, t, alpha, beta, tau):
# #         xdot = x - x**3 - alpha * xt + beta
# #         return xdot

function DoubleGyre end
originalcode(::typeof(DoubleGyre)) = """
class DoubleGyre(DynSys):
    @staticjit
    def _rhs(x, y, z, t, alpha, eps, omega):
        a = eps * np.sin(z)
        b = 1 - 2 * eps * np.sin(z)
        f = a * x ** 2 + b * x
        dx = -alpha * np.pi * np.sin(np.pi * f) * np.cos(np.pi * y)
        dy = alpha * np.pi * np.cos(np.pi * f) * np.sin(np.pi * y) * (2 * a * x + b)
        dz = omega
        return dx, dy, dz

    @staticjit
    def _postprocessing(x, y, z):
        return x, y, np.sin(z)
"""
@doc make_docstring(DoubleGyre) DoubleGyre
function DoubleGyre()
    function rhs(du, u, p, t)
        @unpack alpha, eps, omega = p
        a = eps * sin(u[3])
        b = 1 - 2 * eps * sin(u[3])
        f = a * u[1]^2 + b * u[1]
        du[1] = -alpha * pi * sin(pi * f) * cos(pi * u[2])
        du[2] = alpha * pi * cos(pi * f) * sin(pi * u[2]) * (2 * a * u[1] + b)
        du[3] = omega
    end
    u0 = Float64.(ATTRACTOR_DATA["DoubleGyre"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["DoubleGyre"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function BlinkingRotlet end
originalcode(::typeof(BlinkingRotlet)) = """
class BlinkingRotlet(DynSys):
    @staticjit
    def _rotlet(r, theta, a, b, bc):
        \"\"\"A rotlet velocity field\"\"\"
        kappa = a ** 2 + (b ** 2 * r ** 2) / a ** 2 - 2 * b * r * np.cos(theta)
        gamma = (1 - r ** 2 / a ** 2) * (a ** 2 - (b ** 2 * r ** 2) / a ** 2)
        iota = (b ** 2 * r) / a ** 2 - b * np.cos(theta)
        zeta = b ** 2 + r ** 2 - 2 * b * r * np.cos(theta)
        nu = a ** 2 + b ** 2 - (2 * b ** 2 * r ** 2) / a ** 2
        vr = b * np.sin(theta) * (-bc * (gamma / kappa ** 2) - 1 / kappa + 1 / zeta)
        vth = (
            bc * (gamma * iota) / kappa ** 2
            + bc * r * nu / (a ** 2 * kappa)
            + iota / kappa
            - (r - b * np.cos(theta)) / zeta
        )
        return vr, vth

    @staticjit
    def _protocol(t, tau, stiffness=20):
        return 0.5 + 0.5 * np.tanh(tau * stiffness * np.sin(2 * np.pi * t / tau))

    def rhs(self, X, t):
        r, theta, tt = X
        weight = self._protocol(tt, self.tau)
        dr1, dth1 = self._rotlet(r, theta, self.a, self.b, self.bc)
        dr2, dth2 = self._rotlet(r, theta, self.a, -self.b, self.bc)
        dr = weight * dr1 + (1 - weight) * dr2
        dth = (weight * dth1 + (1 - weight) * dth2) / r
        dtt = 1
        return self.sigma * dr, self.sigma * dth, dtt

    def _postprocessing(self, r, th, tt):
        return r * np.cos(th), r * np.sin(th), np.sin(2 * np.pi * tt / self.tau)
"""
@doc make_docstring(BlinkingRotlet) BlinkingRotlet
function BlinkingRotlet()
    function rotlet(r, theta, a, b, bc)
        kappa = a^2 + (b^2 * r^2) / a^2 - 2 * b * r * cos(theta)
        gamma = (1 - r^2 / a^2) * (a^2 - (b^2 * r^2) / a^2)
        iota = (b^2 * r) / a^2 - b * cos(theta)
        zeta = b^2 + r^2 - 2 * b * r * cos(theta)
        nu = a^2 + b^2 - (2 * b^2 * r^2) / a^2
        vr = b * sin(theta) * (-bc * (gamma / kappa^2) - 1 / kappa + 1 / zeta)
        vth = bc * (gamma * iota) / kappa^2 + bc * r * nu / (a^2 * kappa) + iota / kappa -
              (r - b * cos(theta)) / zeta
        return vr, vth
    end
    function protocol(t, tau, stiffness = 20)
        return 0.5 + 0.5 * tanh(tau * stiffness * sin(2 * pi * t / tau))
    end
    function rhs(du, u, p, t)
        @unpack a, b, bc, tau, sigma = p
        r, theta, tt = u
        weight = protocol(tt, tau)
        dr1, dth1 = rotlet(r, theta, a, b, bc)
        dr2, dth2 = rotlet(r, theta, a, -b, bc)
        dr = weight * dr1 + (1 - weight) * dr2
        dth = (weight * dth1 + (1 - weight) * dth2) / r
        dtt = 1
        du[1] = sigma * dr
        du[2] = sigma * dth
        du[3] = dtt
    end
    u0 = Float64.(ATTRACTOR_DATA["BlinkingRotlet"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["BlinkingRotlet"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function LidDrivenCavityFlow end
function originalcode(::typeof(LidDrivenCavityFlow))
    """
class LidDrivenCavityFlow(DynSys):
    @staticjit
    def _lid(x, y, a, b, tau, u1, u2):
        \"\"\"The velocity field when the left domain drives\"\"\"
        prefactor1 = 2 * u1 * np.sin(np.pi * x / a) / (2 * b * np.pi + a * np.sinh(2 * np.pi * b / a))
        prefactor2 = 2 * u2 * np.sin(2 * np.pi * x / a) / (4 * b * np.pi + a * np.sinh(4 * np.pi * b / a))
        vx1 = -b * np.pi * np.sinh(np.pi * b / a) * np.sinh(np.pi * y / a) + np.cosh(np.pi * b / a) * (np.pi * y * np.cosh(np.pi * y /a) + a * np.sinh(np.pi * y / a))
        vx2 = -2 * b * np.pi * np.sinh(2 * np.pi * b / a) * np.sinh(2 * np.pi * y / a) + np.cosh(2 * np.pi * b / a) * (2 * np.pi * y * np.cosh(2 * np.pi * y / a) + a * np.sinh(2 * np.pi * y / a))
        vx = prefactor1 * vx1 + prefactor2 * vx2

        prefactor1 = 2 * np.pi * u1 * np.cos(np.pi * x / a) / (2 * b * np.pi + a * np.sinh(2 * np.pi * b / a))
        prefactor2 = 4 * np.pi * u2 * np.cos(2 * np.pi * x / a) / (4 * b * np.pi + a * np.sinh(4 * np.pi * b / a))
        vy1 = b * np.sinh(np.pi * b / a) * np.cosh(np.pi * y / a) - np.cosh(np.pi * b / a) * y * np.sinh(np.pi * y / a)
        vy2 = b * np.sinh(2 * np.pi * b / a) * np.cosh(2 * np.pi * y / a) - np.cosh(2 * np.pi * b / a) * y * np.sinh(2 * np.pi * y / a)
        vy = prefactor1 * vy1 + prefactor2 * vy2

        # vy1 = b * np.sinh(np.pi * b / a) * np.cosh(np.pi * y / a) - np.cosh(np.pi * b / a) * y * np.sinh(np.pi * y / a)
        # vy2 = b * np.sinh(2 * np.pi * b / a) * np.cosh(2 * np.pi * y / a) - np.cosh(2 * np.pi * b / a) * y * np.sinh(2 * np.pi * y / a)
        # vy = np.pi * prefactor1 * vy1 + 2 * np.pi * prefactor2 * vy2

        return vx, vy

    # @staticjit
    # def _right(x, y, a, b, tau, u1, u2):
    #     \"\"\"The velocity field when the right domain drives\"\"\"
    #     prefactor1 = 2 * u1 * np.sin(np.pi * x / a) / (2 * b * np.pi + a * np.sinh(2 * np.pi * b / a))
    #     prefactor2 = 2 * u2 * np.sin(2 * np.pi * x / a) / (4 * b * np.pi + a * np.sinh(4 * np.pi * b / a))
    #     vx1 = -b * np.pi * np.sinh(np.pi * b / a) * np.sinh(np.pi * y / a) - np.cosh(np.pi * b / a) * (np.pi * y * np.cosh(np.pi * y /a) + a * np.sinh(np.pi * y /a))
    #     vx2 = -4 * b * np.pi * np.sinh(2 * np.pi * b / a) * np.sinh(2 * np.pi * y / a) - np.cosh(2 * np.pi * b / a) * (2 * np.pi * y * np.cosh(2 * np.pi * y /a) + a * np.sinh(2 * np.pi * y /a))
    #     vx = prefactor1 * vx1 - prefactor2 * vx2

    #     prefactor1 = 2 * np.pi * u1 * np.cos(np.pi * x / a) / (2 * b * np.pi + a * np.sinh(2 * np.pi * b / a))
    #     prefactor2 = 4 * np.pi * u2 * np.cos(2 * np.pi * x / a) / (4 * b * np.pi + a * np.sinh(4 * np.pi * b / a))
    #     vy1 = -b * np.sinh(np.pi * b / a) * np.cosh(np.pi * y / a) + np.cosh(np.pi * b / a) * y * np.sinh(np.pi * y / a)
    #     vy2 = -2 * b * np.sinh(2 * np.pi * b / a) * np.cosh(2 * np.pi * y / a) + np.cosh(2 * np.pi * b / a) * 2 * y * np.sinh(2 * np.pi * y / a)
    #     vy = prefactor1 * vy1 + prefactor2 * vy2

    #     return vx, vy

    @staticjit
    def _protocol(t, tau, stiffness=20):
        return 0.5 + 0.5 * np.tanh(tau * stiffness * np.sin(2 * np.pi * t / tau))

    def rhs(self, X, t):
        x, y, tt = X
        weight = self._protocol(tt, self.tau)
        dx1, dy1 = self._lid(x, y, self.a, self.b, self.tau, self.u1, self.u2)
        dx2, dy2 = self._lid(x, y, self.a, self.b, self.tau, -self.u1, self.u2)
        dx = weight * dx1 + (1 - weight) * dx2
        dy = weight * dy1 + (1 - weight) * dy2
        dtt = 1
        return dx, dy, dtt

    def _postprocessing(self, x, y, tt):
        return x, y, np.sin(2 * np.pi * tt / self.tau)
"""
end
@doc make_docstring(LidDrivenCavityFlow) LidDrivenCavityFlow
function LidDrivenCavityFlow()
    function lid(x, y, a, b, tau, u1, u2)
        prefactor1 = 2 * u1 * sin(pi * x / a) / (2 * b * pi + a * sinh(2 * pi * b / a))
        prefactor2 = 2 * u2 * sin(2 * pi * x / a) / (4 * b * pi + a * sinh(4 * pi * b / a))
        vx1 = -b * pi * sinh(pi * b / a) * sinh(pi * y / a) +
              cosh(pi * b / a) * (pi * y * cosh(pi * y / a) + a * sinh(pi * y / a))
        vx2 = -2 * b * pi * sinh(2 * pi * b / a) * sinh(2 * pi * y / a) +
              cosh(2 * pi * b / a) *
              (2 * pi * y * cosh(2 * pi * y / a) + a * sinh(2 * pi * y / a))
        vx = prefactor1 * vx1 + prefactor2 * vx2

        prefactor1 = 2 * pi * u1 * cos(pi * x / a) / (2 * b * pi + a * sinh(2 * pi * b / a))
        prefactor2 = 4 * pi * u2 * cos(2 * pi * x / a) /
                     (4 * b * pi + a * sinh(4 * pi * b / a))
        vy1 = b * sinh(pi * b / a) * cosh(pi * y / a) -
              cosh(pi * b / a) * y * sinh(pi * y / a)
        vy2 = b * sinh(2 * pi * b / a) * cosh(2 * pi * y / a) -
              cosh(2 * pi * b / a) * y * sinh(2 * pi * y / a)
        vy = prefactor1 * vy1 + prefactor2 * vy2

        return vx, vy
    end
    function protocol(t, tau, stiffness = 20)
        return 0.5 + 0.5 * tanh(tau * stiffness * sin(2 * pi * t / tau))
    end
    function rhs(du, u, p, t)
        @unpack a, b, tau, u1, u2 = p
        x, y, tt = u
        weight = protocol(tt, tau)
        dx1, dy1 = lid(x, y, a, b, tau, u1, u2)
        dx2, dy2 = lid(x, y, a, b, tau, -u1, u2)
        dx = weight * dx1 + (1 - weight) * dx2
        dy = weight * dy1 + (1 - weight) * dy2
        dtt = 1
        du[1] = dx
        du[2] = dy
        du[3] = dtt
    end
    u0 = Float64.(ATTRACTOR_DATA["LidDrivenCavityFlow"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["LidDrivenCavityFlow"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

# class BlinkingVortex(BlinkingRotlet):
#     pass

# class InteriorSquirmer(DynSys):

#     @staticjit
#     def _rhs_static(r, th, t, a, g, n):

#         nvals = np.arange(1, n + 1)
#         sinvals, cosvals = np.sin(th * nvals), np.cos(th * nvals)
#         rnvals = r ** nvals

#         vrn = g * cosvals + a * sinvals
#         vrn *= (nvals * rnvals * (r ** 2 - 1)) / r

#         vth = 2 * r + (r ** 2 - 1) * nvals / r
#         vth *= a * cosvals - g * sinvals
#         vth *= rnvals

#         return np.sum(vrn), np.sum(vth) / r

#     @staticjit
#     def _jac_static(r, th, t, a, g, n):

#         nvals = np.arange(1, n + 1)
#         sinvals, cosvals = np.sin(th * nvals), np.cos(th * nvals)
#         rnvals = r ** nvals
#         trigsum = a * sinvals + g * cosvals
#         trigskew = a * cosvals - g * sinvals

#         j11 = np.copy(trigsum)
#         j11 *= nvals * rnvals * (2 * r ** 2 + (r ** 2 - 1) * (nvals - 1))
#         j11 = (1 / r ** 2) * np.sum(j11)

#         j12 = np.copy(trigskew)
#         j12 *= -(nvals ** 2) * rnvals * (1 - r ** 2) / r
#         j12 = np.sum(j12)

#         j21 = 2 * rnvals * (2 * nvals + 1) * (-np.copy(trigskew))
#         j21 += (n * (1 - r ** 2) * rnvals * (nvals - 1) / r ** 2) * np.copy(
#             g * sinvals + a * cosvals
#         )
#         j21 = -np.sum(j21)

#         j22 = np.copy(trigsum)
#         j22 *= -nvals * rnvals * (2 * r + (r ** 2 - 1) * nvals / r)
#         j22 = np.sum(j22)
#         # (1 / r**2) *

#         ## Correct for polar coordinates
#         vth = np.copy(trigskew)
#         vth *= 2 * r + (r ** 2 - 1) * nvals / r
#         vth *= rnvals
#         vth = np.sum(vth) / r
#         j21 = j21 / r - vth / r
#         j22 /= r

#         return np.array([[j11, j12], [j21, j22]])

#     @staticjit
#     def _protocol(t, tau, stiffness=20):
#         return 0.5 + 0.5 * np.tanh(tau * stiffness * np.sin(2 * np.pi * t / tau))

#     def _postprocessing(self, r, th, tt):
#         return r * np.cos(th), r * np.sin(th), np.sin(2 * np.pi * tt / self.tau)

#     def jac(self, X, t):
#         r, th = X[0], X[1]
#         phase = self._protocol(t, self.tau)
#         return self._jac_static(r, th, t, self.a * phase, self.g * (1 - phase), self.n)

#     def rhs(self, X, t):
#         r, th, tt = X
#         phase = self._protocol(tt, self.tau)
#         dtt = 1
#         dr, dth = self._rhs_static(r, th, t, self.a * phase, self.g * (1 - phase), self.n)
#         return dr, dth, dtt

function OscillatingFlow end
originalcode(::typeof(OscillatingFlow)) = """
class OscillatingFlow(DynSys):
    @staticjit
    def _rhs(x, y, z, t, b, k, omega, u):
        f = x + b * np.sin(z)
        dx = u * np.cos(k * y) * np.sin(k * f)
        dy = -u * np.sin(k * y) * np.cos(k * f)
        dz = omega
        return dx, dy, dz

    def _postprocessing(self, x, y, z):
        return np.cos(self.k * x), y, np.sin(z)
"""
@doc make_docstring(OscillatingFlow) OscillatingFlow
function OscillatingFlow()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack b, k, omega, u = p
        f = x + b * sin(z)
        du[1] = u * cos(k * y) * sin(k * f)
        du[2] = -u * sin(k * y) * cos(k * f)
        du[3] = omega
    end
    u0 = Float64.(ATTRACTOR_DATA["OscillatingFlow"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["OscillatingFlow"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function BickleyJet end
originalcode(::typeof(BickleyJet)) = """
class BickleyJet(DynSys):
    @staticjit
    def _rhs(y, x, z, t, ell, eps, k, omega, sigma, u):
        sechy = 1 / np.cosh(y / ell)
        inds = np.arange(3)
        un = k[inds] * (x - z * sigma[inds])
        dx = u * sechy ** 2 * (-1 - 2 * np.dot(np.cos(un), eps) * np.tanh(y / ell))
        dy = ell * u * sechy ** 2 * np.dot(eps * k, np.sin(un))
        dz = omega
        return dy, dx, dz

    def _postprocessing(self, x, y, z):
        km = np.min(self.k)
        sm = np.min(self.sigma)
        return x, np.sin(km * y), np.sin(self.omega * z * km * sm)
"""
@doc make_docstring(BickleyJet) BickleyJet
function BickleyJet()
    function rhs(du, u, p, t)
        y, x, z = u
        @unpack ell, eps, k, omega, sigma, u = p
        sechy = 1 / cosh(y / ell)
        un = @. k * (x - z * sigma)
        du[1] = u * sechy^2 * (-1 - 2 * dot(cos.(un), eps) * tanh(y / ell))
        du[2] = ell * u * sechy^2 * dot(eps .* k, sin.(un))
        du[3] = omega
    end
    u0 = Float64.(ATTRACTOR_DATA["BickleyJet"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["BickleyJet"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function ArnoldBeltramiChildress end
originalcode(::typeof(ArnoldBeltramiChildress)) = """
class ArnoldBeltramiChildress(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c):
        dx = a * np.sin(z) + c * np.cos(y)
        dy = b * np.sin(x) + a * np.cos(z)
        dz = c * np.sin(y) + b * np.cos(x)
        return dx, dy, dz

    @staticjit
    def _postprocessing(x, y, z):
        return np.sin(x), np.cos(y), np.sin(z)
"""
@doc make_docstring(ArnoldBeltramiChildress) ArnoldBeltramiChildress
function ArnoldBeltramiChildress()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c = p
        du[1] = a * sin(z) + c * cos(y)
        du[2] = b * sin(x) + a * cos(z)
        du[3] = c * sin(y) + b * cos(x)
    end
    u0 = Float64.(ATTRACTOR_DATA["ArnoldBeltramiChildress"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["ArnoldBeltramiChildress"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function JerkCircuit end
originalcode(::typeof(JerkCircuit)) = """
class JerkCircuit(DynSys):
    @staticjit
    def _rhs(x, y, z, t, eps, y0):
        xdot = y
        ydot = z
        zdot = -z - x - eps * (np.exp(y / y0) - 1)
        return xdot, ydot, zdot
"""
@doc make_docstring(JerkCircuit) JerkCircuit
function JerkCircuit()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack eps, y0 = p
        du[1] = y
        du[2] = z
        du[3] = -z - x - eps * (exp(y / y0) - 1)
    end
    u0 = Float64.(ATTRACTOR_DATA["JerkCircuit"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["JerkCircuit"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function ForcedBrusselator end
originalcode(::typeof(ForcedBrusselator)) = """
class ForcedBrusselator(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, f, w):
        xdot = a + x ** 2 * y - (b + 1) * x + f * np.cos(z)
        ydot = b * x - x ** 2 * y
        zdot = w
        return xdot, ydot, zdot

    @staticjit
    def _postprocessing(x, y, z):
        return x, y, np.sin(z)
"""
@doc make_docstring(ForcedBrusselator) ForcedBrusselator
function ForcedBrusselator()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, f, w = p
        du[1] = a + x^2 * y - (b + 1) * x + f * cos(z)
        du[2] = b * x - x^2 * y
        du[3] = w
    end
    u0 = Float64.(ATTRACTOR_DATA["ForcedBrusselator"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["ForcedBrusselator"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function WindmiReduced end
originalcode(::typeof(WindmiReduced)) = """
class WindmiReduced(DynSys):
    @staticjit
    def _rhs(i, v, p, t, a1, b1, b2, b3, d1, vsw):
        idot = a1 * (vsw - v)
        vdot = b1 * i - b2 * p ** 1 / 2 - b3 * v
        pdot = (
            vsw ** 2 - p ** (5 / 4) * vsw ** (1 / 2) * (1 + np.tanh(d1 * (i - 1))) / 2
        )
        return idot, vdot, pdot
"""
@doc make_docstring(WindmiReduced) WindmiReduced
function WindmiReduced()
    function rhs(du, u, p, t)
        @unpack a1, b1, b2, b3, d1, vsw = p
        i, v, p = u
        du[1] = a1 * (vsw - v)
        du[2] = b1 * i - b2 * p^(1 / 2) - b3 * v
        du[3] = vsw^2 - p^(5 / 4) * vsw^(1 / 2) * (1 + tanh(d1 * (i - 1))) / 2
    end
    u0 = Float64.(ATTRACTOR_DATA["WindmiReduced"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["WindmiReduced"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function MooreSpiegel end
originalcode(::typeof(MooreSpiegel)) = """
class MooreSpiegel(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, eps):
        xdot = y
        ydot = a * z
        zdot = -z + eps * y - y * x ** 2 - b * x
        return xdot, ydot, zdot
"""
@doc make_docstring(MooreSpiegel) MooreSpiegel
function MooreSpiegel()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, eps = p
        du[1] = y
        du[2] = a * z
        du[3] = -z + eps * y - y * x^2 - b * x
    end
    u0 = Float64.(ATTRACTOR_DATA["MooreSpiegel"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["MooreSpiegel"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function CoevolvingPredatorPrey end
originalcode(::typeof(CoevolvingPredatorPrey)) = """
class CoevolvingPredatorPrey(DynSys):
    @staticjit
    def _rhs(x, y, alpha, t, a1, a2, a3, b1, b2, d1, d2, delta, k1, k2, k4, vv):
        xdot = x * (
            -((a3 * y) / (1 + b2 * x))
            + (a1 * alpha * (1 - k1 * x * (-alpha + alpha * delta))) / (1 + b1 * alpha)
            - d1
            * (
                1
                - k2 * (-(alpha ** 2) + (alpha * delta) ** 2)
                + k4 * (-(alpha ** 4) + (alpha * delta) ** 4)
            )
        )
        ydot = (-d2 + (a2 * x) / (1 + b2 * x)) * y
        alphadot = vv * (
            -((a1 * k1 * x * alpha * delta) / (1 + b1 * alpha))
            - d1 * (-2 * k2 * alpha * delta ** 2 + 4 * k4 * alpha ** 3 * delta ** 4)
        )
        return xdot, ydot, alphadot
"""
@doc make_docstring(CoevolvingPredatorPrey) CoevolvingPredatorPrey
function CoevolvingPredatorPrey()
    function rhs(du, u, p, t)
        x, y, alpha = u
        @unpack a1, a2, a3, b1, b2, d1, d2, delta, k1, k2, k4, vv = p
        du[1] = x * (-((a3 * y) / (1 + b2 * x))
                 +
                 (a1 * alpha * (1 - k1 * x * (-alpha + alpha * delta))) / (1 + b1 * alpha)
                 -
                 d1
                 *
                 (1
                  -
                  k2 * (-alpha^2 + (alpha * delta)^2)
                  +
                  k4 * (-alpha^4 + (alpha * delta)^4)))
        du[2] = (-d2 + (a2 * x) / (1 + b2 * x)) * y
        du[3] = vv * (-((a1 * k1 * x * alpha * delta) / (1 + b1 * alpha))
                 -
                 d1 * (-2 * k2 * alpha * delta^2 + 4 * k4 * alpha^3 * delta^4))
    end
    u0 = Float64.(ATTRACTOR_DATA["CoevolvingPredatorPrey"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["CoevolvingPredatorPrey"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function KawczynskiStrizhak end
originalcode(::typeof(KawczynskiStrizhak)) = """
class KawczynskiStrizhak(DynSys):
    @staticjit
    def _rhs(x, y, z, t, beta, gamma, kappa, mu):
        xdot = gamma * y - gamma * x ** 3 + 3 * mu * gamma * x
        ydot = -2 * mu * x - y - z + beta
        zdot = kappa * x - kappa * z
        return xdot, ydot, zdot
"""
@doc make_docstring(KawczynskiStrizhak) KawczynskiStrizhak
function KawczynskiStrizhak()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack beta, gamma, kappa, mu = p
        du[1] = gamma * y - gamma * x^3 + 3 * mu * gamma * x
        du[2] = -2 * mu * x - y - z + beta
        du[3] = kappa * x - kappa * z
    end
    u0 = Float64.(ATTRACTOR_DATA["KawczynskiStrizhak"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["KawczynskiStrizhak"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function BelousovZhabotinsky end
originalcode(::typeof(BelousovZhabotinsky)) = """
class BelousovZhabotinsky(DynSys):
    @staticjit
    def _rhs(
        x,
        z,
        v,
        t,
        c1,
        c10,
        c11,
        c12,
        c13,
        c2,
        c3,
        c4,
        c5,
        c6,
        c7,
        c8,
        c9,
        ci,
        kf,
        t0,
        y0,
        yb1,
        yb2,
        yb3,
        z0,
    ):
        ybar = (1 / y0) * yb1 * z * v / (yb2 * x + yb3 + kf)
        if x < 0.0:
            x = 0
        rf = (ci - z0 * z) * np.sqrt(x)
        xdot = c1 * x * ybar + c2 * ybar + c3 * x ** 2 + c4 * rf + c5 * x * z - kf * x
        zdot = (c6 / z0) * rf + c7 * x * z + c8 * z * v + c9 * z - kf * z
        vdot = c10 * x * ybar + c11 * ybar + c12 * x ** 2 + c13 * z * v - kf * v
        return xdot * t0, zdot * t0, vdot * t0
"""
@doc make_docstring(BelousovZhabotinsky) BelousovZhabotinsky
function BelousovZhabotinsky()
    function rhs(du, u, p, t)
        x, z, v = u
        @unpack c1, c10, c11, c12, c13, c2, c3, c4, c5, c6, c7, c8, c9, ci, kf, t0, y0, yb1, yb2, yb3, z0 = p
        ybar = (1 / y0) * yb1 * z * v / (yb2 * x + yb3 + kf)
        if x < 0.0
            x = 0
        end
        rf = (ci - z0 * z) * sqrt(x)
        du[1] = c1 * x * ybar + c2 * ybar + c3 * x^2 + c4 * rf + c5 * x * z - kf * x
        du[2] = (c6 / z0) * rf + c7 * x * z + c8 * z * v + c9 * z - kf * z
        du[3] = c10 * x * ybar + c11 * ybar + c12 * x^2 + c13 * z * v - kf * v
    end
    u0 = Float64.(ATTRACTOR_DATA["BelousovZhabotinsky"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["BelousovZhabotinsky"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function IsothermalChemical end
originalcode(::typeof(IsothermalChemical)) = """
class IsothermalChemical(DynSys):
    @staticmethod
    def _rhs(alpha, beta, gamma, t, delta, kappa, mu, sigma):
        alphadot = mu * (kappa + gamma) - alpha * beta ** 2 - alpha
        betadot = (alpha * beta ** 2 + alpha - beta) / sigma
        gammadot = (beta - gamma) / delta
        return alphadot, betadot, gammadot
"""
@doc make_docstring(IsothermalChemical) IsothermalChemical
function IsothermalChemical()
    function rhs(du, u, p, t)
        alpha, beta, gamma = u
        @unpack delta, kappa, mu, sigma = p
        du[1] = mu * (kappa + gamma) - alpha * beta^2 - alpha
        du[2] = (alpha * beta^2 + alpha - beta) / sigma
        du[3] = (beta - gamma) / delta
    end
    u0 = Float64.(ATTRACTOR_DATA["IsothermalChemical"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["IsothermalChemical"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function VallisElNino end
originalcode(::typeof(VallisElNino)) = """
class VallisElNino(DynSys):
    @staticmethod
    def _rhs(x, y, z, t, b, c, p):
        xdot = b * y - c * x - c * p
        ydot = -y + x * z
        zdot = -z - x * y + 1
        return xdot, ydot, zdot
"""
@doc make_docstring(VallisElNino) VallisElNino
function VallisElNino()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack b, c, p = p
        du[1] = b * y - c * x - c * p
        du[2] = -y + x * z
        du[3] = -z - x * y + 1
    end
    u0 = Float64.(ATTRACTOR_DATA["VallisElNino"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["VallisElNino"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function RabinovichFabrikant end
originalcode(::typeof(RabinovichFabrikant)) = """
class RabinovichFabrikant(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, g):
        xdot = y * z - y + y * x ** 2 + g * x
        ydot = 3 * x * z + x - x ** 3 + g * y
        zdot = -2 * a * z  - 2 * x * y * z
        return (xdot, ydot, zdot)
"""
@doc make_docstring(RabinovichFabrikant) RabinovichFabrikant
function RabinovichFabrikant()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, g = p
        du[1] = y * z - y + y * x^2 + g * x
        du[2] = 3 * x * z + x - x^3 + g * y
        du[3] = -2 * a * z - 2 * x * y * z
    end
    u0 = Float64.(ATTRACTOR_DATA["RabinovichFabrikant"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["RabinovichFabrikant"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function NoseHoover end
originalcode(::typeof(NoseHoover)) = """
class NoseHoover(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a):
        xdot = y
        ydot = -x + y * z
        zdot = a - y ** 2
        return xdot, ydot, zdot
"""
@doc make_docstring(NoseHoover) NoseHoover
function NoseHoover()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a = p
        du[1] = y
        du[2] = -x + y * z
        du[3] = a - y^2
    end
    u0 = Float64.(ATTRACTOR_DATA["NoseHoover"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["NoseHoover"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Dadras end
originalcode(::typeof(Dadras)) = """
class Dadras(DynSys):
    @staticjit
    def _rhs(x, y, z, t, c, e, o, p, r):
        xdot = y - p * x + o * y * z
        ydot = r * y - x * z + z
        zdot = c * x * y - e * z
        return xdot, ydot, zdot
"""
@doc make_docstring(Dadras) Dadras
function Dadras()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack c, e, o, p, r = p
        du[1] = y - p * x + o * y * z
        du[2] = r * y - x * z + z
        du[3] = c * x * y - e * z
    end
    u0 = Float64.(ATTRACTOR_DATA["Dadras"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Dadras"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function RikitakeDynamo end
originalcode(::typeof(RikitakeDynamo)) = """
class RikitakeDynamo(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, mu):
        xdot = -mu * x + y * z
        ydot = -mu * y - a * x + x * z
        zdot = 1 - x * y
        return xdot, ydot, zdot
"""
@doc make_docstring(RikitakeDynamo) RikitakeDynamo
function RikitakeDynamo()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, mu = p
        du[1] = -mu * x + y * z
        du[2] = -mu * y - a * x + x * z
        du[3] = 1 - x * y
    end
    u0 = Float64.(ATTRACTOR_DATA["RikitakeDynamo"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["RikitakeDynamo"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function NuclearQuadrupole end
function originalcode(::typeof(NuclearQuadrupole))
    """
class NuclearQuadrupole(DynSys):
    @staticjit
    def _rhs(q1, q2, p1, p2, t, a, b, d):
        q1dot = a * p1
        q2dot = a * p2
        p1dot = - a * q1 + 3 / np.sqrt(2) * b * q1 ** 2 - 3 / np.sqrt(2) * b * q2 ** 2 - d * q1 ** 3 - d * q1 * q2 ** 2
        p2dot = -a * q2 - 3 * np.sqrt(2) * b * q1 * q2 - d * q2 * q1 ** 2 - d * q2 ** 3
        return q1dot, q2dot, p1dot, p2dot
"""
end
@doc make_docstring(NuclearQuadrupole) NuclearQuadrupole
function NuclearQuadrupole()
    function rhs(du, u, p, t)
        q1, q2, p1, p2 = u
        @unpack a, b, d = p
        du[1] = a * p1
        du[2] = a * p2
        du[3] = -a * q1 + 3 / sqrt(2) * b * q1^2 - 3 / sqrt(2) * b * q2^2 - d * q1^3 -
                d * q1 * q2^2
        du[4] = -a * q2 - 3 * sqrt(2) * b * q1 * q2 - d * q2 * q1^2 - d * q2^3
    end
    u0 = Float64.(ATTRACTOR_DATA["NuclearQuadrupole"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["NuclearQuadrupole"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function PehlivanWei end
originalcode(::typeof(PehlivanWei)) = """
class PehlivanWei(DynSys):
    @staticjit
    def _rhs(x, y, z, t):
        xdot = y - y * z
        ydot = y + y * z - 2 * x
        zdot = 2 - x * y - y ** 2
        return xdot, ydot, zdot
"""
@doc make_docstring(PehlivanWei) PehlivanWei
function PehlivanWei()
    function rhs(du, u, p, t)
        x, y, z = u
        du[1] = y - y * z
        du[2] = y + y * z - 2 * x
        du[3] = 2 - x * y - y^2
    end
    u0 = Float64.(ATTRACTOR_DATA["PehlivanWei"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["PehlivanWei"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottTorus end
originalcode(::typeof(SprottTorus)) = """
class SprottTorus(DynSys):
    @staticjit
    def _rhs(x, y, z, t):
        xdot = y + 2 * x * y + x * z
        ydot = 1 - 2 * x ** 2 + y * z
        zdot = x - x ** 2 - y ** 2
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottTorus) SprottTorus
function SprottTorus()
    function rhs(du, u, p, t)
        x, y, z = u
        du[1] = y + 2 * x * y + x * z
        du[2] = 1 - 2 * x^2 + y * z
        du[3] = x - x^2 - y^2
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottTorus"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottTorus"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottJerk end
originalcode(::typeof(SprottJerk)) = """
class SprottJerk(DynSys):
    @staticjit
    def _rhs(x, y, z, t, mu):
        xdot = y
        ydot = z
        zdot = -x + y ** 2 - mu * z
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottJerk) SprottJerk
function SprottJerk()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack mu = p
        du[1] = y
        du[2] = z
        du[3] = -x + y^2 - mu * z
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottJerk"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottJerk"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

# ## Not chaotic
# # class JerkCircuit(DynSys):
# #     def rhs(self, X, t):
# #         x, y, z = X
# #         xdot = y
# #         ydot = z
# #         zdot = -z - x  - self.eps*(np.exp(y/self.y0) - 1)
# #         return (xdot, ydot, zdot)

function SprottA end
originalcode(::typeof(SprottA)) = """
class SprottA(DynSys):
    @staticjit
    def _rhs(x, y, z, t):
        xdot = y
        ydot = -x + y * z
        zdot = 1 - y ** 2
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottA) SprottA
function SprottA()
    function rhs(du, u, p, t)
        x, y, z = u
        du[1] = y
        du[2] = -x + y * z
        du[3] = 1 - y^2
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottA"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottA"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottB end
originalcode(::typeof(SprottB)) = """
class SprottB(DynSys):
    @staticjit
    def _rhs(x, y, z, t):
        xdot = y * z
        ydot = x - y
        zdot = 1 - x * y
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottB) SprottB
function SprottB()
    function rhs(du, u, p, t)
        x, y, z = u
        du[1] = y * z
        du[2] = x - y
        du[3] = 1 - x * y
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottB"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottB"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottC end
originalcode(::typeof(SprottC)) = """
class SprottC(DynSys):
    @staticjit
    def _rhs(x, y, z, t):
        xdot = y * z
        ydot = x - y
        zdot = 1 - x ** 2
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottC) SprottC
function SprottC()
    function rhs(du, u, p, t)
        x, y, z = u
        du[1] = y * z
        du[2] = x - y
        du[3] = 1 - x^2
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottC"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottC"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottD end
originalcode(::typeof(SprottD)) = """
class SprottD(DynSys):
    @staticjit
    def _rhs(x, y, z, t):
        xdot = -y
        ydot = x + z
        zdot = x * z + 3 * y ** 2
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottD) SprottD
function SprottD()
    function rhs(du, u, p, t)
        x, y, z = u
        du[1] = -y
        du[2] = x + z
        du[3] = x * z + 3 * y^2
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottD"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottD"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottE end
originalcode(::typeof(SprottE)) = """
class SprottE(DynSys):
    @staticjit
    def _rhs(x, y, z, t):
        xdot = y * z
        ydot = x ** 2 - y
        zdot = 1 - 4 * x
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottE) SprottE
function SprottE()
    function rhs(du, u, p, t)
        x, y, z = u
        du[1] = y * z
        du[2] = x^2 - y
        du[3] = 1 - 4 * x
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottE"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottE"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottF end
originalcode(::typeof(SprottF)) = """
class SprottF(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a):
        xdot = y + z
        ydot = -x + a * y
        zdot = x ** 2 - z
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottF) SprottF
function SprottF()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a = p
        du[1] = y + z
        du[2] = -x + a * y
        du[3] = x^2 - z
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottF"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottF"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottG end
originalcode(::typeof(SprottG)) = """
class SprottG(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a):
        xdot = a * x + z
        ydot = x * z - y
        zdot = -x + y
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottG) SprottG
function SprottG()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a = p
        du[1] = a * x + z
        du[2] = x * z - y
        du[3] = -x + y
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottG"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottG"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottH end
originalcode(::typeof(SprottH)) = """
class SprottH(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a):
        xdot = -y + z ** 2
        ydot = x + a * y
        zdot = x - z
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottH) SprottH
function SprottH()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a = p
        du[1] = -y + z^2
        du[2] = x + a * y
        du[3] = x - z
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottH"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottH"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottI end
originalcode(::typeof(SprottI)) = """
class SprottI(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a):
        xdot = -a * y
        ydot = x + z
        zdot = x + y ** 2 - z
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottI) SprottI
function SprottI()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a = p
        du[1] = -a * y
        du[2] = x + z
        du[3] = x + y^2 - z
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottI"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottI"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottJ end
originalcode(::typeof(SprottJ)) = """
class SprottJ(DynSys):
    @staticjit
    def _rhs(x, y, z, t):
        xdot = 2 * z
        ydot = -2 * y + z
        zdot = -x + y + y ** 2
        return (xdot, ydot, zdot)
"""
@doc make_docstring(SprottJ) SprottJ
function SprottJ()
    function rhs(du, u, p, t)
        x, y, z = u
        du[1] = 2 * z
        du[2] = -2 * y + z
        du[3] = -x + y + y^2
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottJ"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottJ"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottK end
originalcode(::typeof(SprottK)) = """
class SprottK(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a):
        xdot = x * y - z
        ydot = x - y
        zdot = x + a * z
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottK) SprottK
function SprottK()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a = p
        du[1] = x * y - z
        du[2] = x - y
        du[3] = x + a * z
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottK"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottK"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottL end
originalcode(::typeof(SprottL)) = """
class SprottL(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b):
        xdot = y + b * z
        ydot = a * x ** 2 - y
        zdot = 1 - x
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottL) SprottL
function SprottL()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b = p
        du[1] = y + b * z
        du[2] = a * x^2 - y
        du[3] = 1 - x
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottL"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottL"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottM end
originalcode(::typeof(SprottM)) = """
class SprottM(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a):
        xdot = -z
        ydot = -x ** 2 - y
        zdot = a + a * x + y
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottM) SprottM
function SprottM()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a = p
        du[1] = -z
        du[2] = -x^2 - y
        du[3] = a + a * x + y
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottM"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottM"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottN end
originalcode(::typeof(SprottN)) = """
class SprottN(DynSys):
    @staticjit
    def _rhs(x, y, z, t):
        xdot = -2 * y
        ydot = x + z ** 2
        zdot = 1 + y - 2 * z
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottN) SprottN
function SprottN()
    function rhs(du, u, p, t)
        x, y, z = u
        du[1] = -2 * y
        du[2] = x + z^2
        du[3] = 1 + y - 2 * z
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottN"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottN"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottO end
originalcode(::typeof(SprottO)) = """
class SprottO(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a):
        xdot = y
        ydot = x - z
        zdot = x + x * z + a * y
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottO) SprottO
function SprottO()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a = p
        du[1] = y
        du[2] = x - z
        du[3] = x + x * z + a * y
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottO"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottO"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottP end
originalcode(::typeof(SprottP)) = """
class SprottP(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a):
        xdot = a * y + z
        ydot = -x + y ** 2
        zdot = x + y
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottP) SprottP
function SprottP()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a = p
        du[1] = a * y + z
        du[2] = -x + y^2
        du[3] = x + y
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottP"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottP"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottQ end
originalcode(::typeof(SprottQ)) = """
class SprottQ(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b):
        xdot = -z
        ydot = x - y
        zdot = a * x + y ** 2 + b * z
        return (xdot, ydot, zdot)
"""
@doc make_docstring(SprottQ) SprottQ
function SprottQ()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b = p
        du[1] = -z
        du[2] = x - y
        du[3] = a * x + y^2 + b * z
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottQ"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottQ"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottR end
originalcode(::typeof(SprottR)) = """
class SprottR(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b):
        xdot = a - y
        ydot = b + z
        zdot = x * y - z
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottR) SprottR
function SprottR()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b = p
        du[1] = a - y
        du[2] = b + z
        du[3] = x * y - z
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottR"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottR"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottS end
originalcode(::typeof(SprottS)) = """
class SprottS(DynSys):
    @staticjit
    def _rhs(x, y, z, t):
        xdot = -x - 4 * y
        ydot = x + z ** 2
        zdot = 1 + x
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottS) SprottS
function SprottS()
    function rhs(du, u, p, t)
        x, y, z = u
        du[1] = -x - 4 * y
        du[2] = x + z^2
        du[3] = 1 + x
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottS"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottS"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SprottMore end
originalcode(::typeof(SprottMore)) = """
class SprottMore(DynSys):
    @staticjit
    def _rhs(x, y, z, t):
        xdot = y
        ydot = -x - np.sign(z) * y
        zdot = y ** 2 - np.exp(-(x ** 2))
        return xdot, ydot, zdot
"""
@doc make_docstring(SprottMore) SprottMore
function SprottMore()
    function rhs(du, u, p, t)
        x, y, z = u
        du[1] = y
        du[2] = -x - sign(z) * y
        du[3] = y^2 - exp(-(x^2))
    end
    u0 = Float64.(ATTRACTOR_DATA["SprottMore"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SprottMore"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Arneodo end
originalcode(::typeof(Arneodo)) = """
class Arneodo(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c, d):
        xdot = y
        ydot = z
        zdot = -a * x - b * y - c * z + d * x ** 3
        return xdot, ydot, zdot
"""
@doc make_docstring(Arneodo) Arneodo
function Arneodo()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c, d = p
        du[1] = y
        du[2] = z
        du[3] = -a * x - b * y - c * z + d * x^3
    end
    u0 = Float64.(ATTRACTOR_DATA["Arneodo"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Arneodo"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

# class Coullet(Arneodo):
#     pass

function Rucklidge end
originalcode(::typeof(Rucklidge)) = """
class Rucklidge(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b):
        xdot = -a * x + b * y - y * z
        ydot = x
        zdot = -z + y ** 2
        return xdot, ydot, zdot
"""
@doc make_docstring(Rucklidge) Rucklidge
function Rucklidge()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b = p
        du[1] = -a * x + b * y - y * z
        du[2] = x
        du[3] = -z + y^2
    end
    u0 = Float64.(ATTRACTOR_DATA["Rucklidge"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Rucklidge"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Sakarya end
originalcode(::typeof(Sakarya)) = """
class Sakarya(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c, h, p, q, r, s):
        xdot = a * x + h * y + s * y * z
        ydot = -b * y - p * x + q * x * z
        zdot = c * z - r * x * y
        return xdot, ydot, zdot
"""
@doc make_docstring(Sakarya) Sakarya
function Sakarya()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c, h, p, q, r, s = p
        du[1] = a * x + h * y + s * y * z
        du[2] = -b * y - p * x + q * x * z
        du[3] = c * z - r * x * y
    end
    u0 = Float64.(ATTRACTOR_DATA["Sakarya"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Sakarya"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

# class LiuChen(Sakarya):
#     pass

function RayleighBenard end
originalcode(::typeof(RayleighBenard)) = """
class RayleighBenard(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, r):
        xdot = a * y - a * x
        ydot = r * y - x * z
        zdot = x * y - b * z
        return xdot, ydot, zdot
"""
@doc make_docstring(RayleighBenard) RayleighBenard
function RayleighBenard()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, r = p
        du[1] = a * y - a * x
        du[2] = r * y - x * z
        du[3] = x * y - b * z
    end
    u0 = Float64.(ATTRACTOR_DATA["RayleighBenard"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["RayleighBenard"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Finance end
originalcode(::typeof(Finance)) = """
class Finance(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c):
        xdot = (1 / b - a) * x + z + x * y
        ydot = -b * y - x ** 2
        zdot = -x - c * z
        return xdot, ydot, zdot
"""
@doc make_docstring(Finance) Finance
function Finance()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c = p
        du[1] = (1 / b - a) * x + z + x * y
        du[2] = -b * y - x^2
        du[3] = -x - c * z
    end
    u0 = Float64.(ATTRACTOR_DATA["Finance"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Finance"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Bouali2 end
originalcode(::typeof(Bouali2)) = """
class Bouali2(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, bb, c, g, m, y0):
        xdot = a * y0 * x - a * x * y - b * z
        ydot = -g * y + g * y * x ** 2
        zdot = -1.5 * m * x + m * bb * x * z - c * z
        return xdot, ydot, zdot
"""
@doc make_docstring(Bouali2) Bouali2
function Bouali2()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, bb, c, g, m, y0 = p
        du[1] = a * y0 * x - a * x * y - b * z
        du[2] = -g * y + g * y * x^2
        du[3] = -1.5 * m * x + m * bb * x * z - c * z
    end
    u0 = Float64.(ATTRACTOR_DATA["Bouali2"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Bouali2"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

# class Bouali(Bouali2):
#     pass

function LuChenCheng end
originalcode(::typeof(LuChenCheng)) = """
class LuChenCheng(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c):
        xdot = -(a * b) / (a + b) * x - y * z + c
        ydot = a * y + x * z
        zdot = b * z + x * y
        return xdot, ydot, zdot
"""
@doc make_docstring(LuChenCheng) LuChenCheng
function LuChenCheng()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c = p
        du[1] = -(a * b) / (a + b) * x - y * z + c
        du[2] = a * y + x * z
        du[3] = b * z + x * y
    end
    u0 = Float64.(ATTRACTOR_DATA["LuChenCheng"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["LuChenCheng"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function LuChen end
originalcode(::typeof(LuChen)) = """
class LuChen(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c):
        xdot = a * y - a * x
        ydot = -x * z + c * y
        zdot = x * y - b * z
        return xdot, ydot, zdot
"""
@doc make_docstring(LuChen) LuChen
function LuChen()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c = p
        du[1] = a * y - a * x
        du[2] = -x * z + c * y
        du[3] = x * y - b * z
    end
    u0 = Float64.(ATTRACTOR_DATA["LuChen"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["LuChen"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function QiChen end
originalcode(::typeof(QiChen)) = """
class QiChen(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c):
        xdot = a * y - a * x + y * z
        ydot = c * x + y - x * z
        zdot = x * y - b * z
        return xdot, ydot, zdot
"""
@doc make_docstring(QiChen) QiChen
function QiChen()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c = p
        du[1] = a * y - a * x + y * z
        du[2] = c * x + y - x * z
        du[3] = x * y - b * z
    end
    u0 = Float64.(ATTRACTOR_DATA["QiChen"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["QiChen"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function ZhouChen end
originalcode(::typeof(ZhouChen)) = """
class ZhouChen(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c, d, e):
        xdot = a * x + b * y + y * z
        ydot = c * y - x * z + d * y * z
        zdot = e * z - x * y
        return xdot, ydot, zdot
"""
@doc make_docstring(ZhouChen) ZhouChen
function ZhouChen()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c, d, e = p
        du[1] = a * x + b * y + y * z
        du[2] = c * y - x * z + d * y * z
        du[3] = e * z - x * y
    end
    u0 = Float64.(ATTRACTOR_DATA["ZhouChen"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["ZhouChen"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function BurkeShaw end
originalcode(::typeof(BurkeShaw)) = """
class BurkeShaw(DynSys):
    @staticjit
    def _rhs(x, y, z, t, e, n):
        xdot = -n * x - n * y
        ydot = y - n * x * z
        zdot = n * x * y + e
        return xdot, ydot, zdot
"""
@doc make_docstring(BurkeShaw) BurkeShaw
function BurkeShaw()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack e, n = p
        du[1] = -n * x - n * y
        du[2] = y - n * x * z
        du[3] = n * x * y + e
    end
    u0 = Float64.(ATTRACTOR_DATA["BurkeShaw"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["BurkeShaw"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Chen end
originalcode(::typeof(Chen)) = """
class Chen(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c):
        xdot = a * y - a * x
        ydot = (c - a) * x - x * z + c * y
        zdot = x * y - b * z
        return xdot, ydot, zdot
"""
@doc make_docstring(Chen) Chen
function Chen()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c = p
        du[1] = a * y - a * x
        du[2] = (c - a) * x - x * z + c * y
        du[3] = x * y - b * z
    end
    u0 = Float64.(ATTRACTOR_DATA["Chen"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Chen"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function ChenLee end
originalcode(::typeof(ChenLee)) = """
class ChenLee(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c):
        xdot = a * x - y * z
        ydot = b * y + x * z
        zdot = c * z + 0.3333333333333333333333333 * x * y
        return xdot, ydot, zdot
"""
@doc make_docstring(ChenLee) ChenLee
function ChenLee()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c = p
        du[1] = a * x - y * z
        du[2] = b * y + x * z
        du[3] = c * z + 0.3333333333333333333333333 * x * y
    end
    u0 = Float64.(ATTRACTOR_DATA["ChenLee"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["ChenLee"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function WangSun end
originalcode(::typeof(WangSun)) = """
class WangSun(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, d, e, f, q):
        xdot = a * x + q * y * z
        ydot = b * x + d * y - x * z
        zdot = e * z + f * x * y
        return xdot, ydot, zdot
"""
@doc make_docstring(WangSun) WangSun
function WangSun()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, d, e, f, q = p
        du[1] = a * x + q * y * z
        du[2] = b * x + d * y - x * z
        du[3] = e * z + f * x * y
    end
    u0 = Float64.(ATTRACTOR_DATA["WangSun"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["WangSun"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function YuWang end
originalcode(::typeof(YuWang)) = """
class YuWang(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c, d):
        xdot = a * (y - x)
        ydot = b * x - c * x * z
        zdot = np.exp(x * y) - d * z
        return xdot, ydot, zdot
"""
@doc make_docstring(YuWang) YuWang
function YuWang()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c, d = p
        du[1] = a * (y - x)
        du[2] = b * x - c * x * z
        du[3] = exp(x * y) - d * z
    end
    u0 = Float64.(ATTRACTOR_DATA["YuWang"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["YuWang"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function YuWang2 end
originalcode(::typeof(YuWang2)) = """
class YuWang2(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c, d):
        xdot = a * (y - x)
        ydot = b * x - c * x * z
        zdot = np.cosh(x * y) - d * z
        return xdot, ydot, zdot
"""
@doc make_docstring(YuWang2) YuWang2
function YuWang2()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c, d = p
        du[1] = a * (y - x)
        du[2] = b * x - c * x * z
        du[3] = cosh(x * y) - d * z
    end
    u0 = Float64.(ATTRACTOR_DATA["YuWang2"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["YuWang2"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SanUmSrisuchinwong end
originalcode(::typeof(SanUmSrisuchinwong)) = """
class SanUmSrisuchinwong(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a):
        xdot = y - x
        ydot = -z * np.tanh(x)
        zdot = -a + x * y + np.abs(y)
        return xdot, ydot, zdot
"""
@doc make_docstring(SanUmSrisuchinwong) SanUmSrisuchinwong
function SanUmSrisuchinwong()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a = p
        du[1] = y - x
        du[2] = -z * tanh(x)
        du[3] = -a + x * y + abs(y)
    end
    u0 = Float64.(ATTRACTOR_DATA["SanUmSrisuchinwong"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SanUmSrisuchinwong"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function DequanLi end
originalcode(::typeof(DequanLi)) = """
class DequanLi(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, c, d, eps, f, k):
        xdot = a * y - a * x + d * x * z
        ydot = k * x + f * y - x * z
        zdot = c * z + x * y - eps * x ** 2
        return xdot, ydot, zdot
"""
@doc make_docstring(DequanLi) DequanLi
function DequanLi()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, c, d, eps, f, k = p
        du[1] = a * y - a * x + d * x * z
        du[2] = k * x + f * y - x * z
        du[3] = c * z + x * y - eps * x^2
    end
    u0 = Float64.(ATTRACTOR_DATA["DequanLi"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["DequanLi"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

# class PanXuZhou(DequanLi):
#     pass

# class Tsucs2(DequanLi):
#     pass

function ArnoldWeb end
originalcode(::typeof(ArnoldWeb)) = """
class ArnoldWeb(DynSys):
    @staticjit
    def _rhs(p1, p2, x1, x2, z, t, mu, w):
        denom = 4 + np.cos(z) + np.cos(x1) + np.cos(x2)
        p1dot = -mu * np.sin(x1) / denom ** 2
        p2dot = -mu * np.sin(x2) / denom ** 2
        x1dot = p1
        x2dot = p2
        zdot = w
        return p1dot, p2dot, x1dot, x2dot, zdot

    @staticjit
    def _postprocessing(p1, p2, x1, x2, z):
        return p1, p2, np.sin(x1), np.sin(x2), np.cos(z)
"""
@doc make_docstring(ArnoldWeb) ArnoldWeb
function ArnoldWeb()
    function rhs(du, u, p, t)
        p1, p2, x1, x2, z = u
        @unpack mu, w = p
        denom = 4 + cos(z) + cos(x1) + cos(x2)
        du[1] = -mu * sin(x1) / denom^2
        du[2] = -mu * sin(x2) / denom^2
        du[3] = p1
        du[4] = p2
        du[5] = w
    end
    u0 = Float64.(ATTRACTOR_DATA["ArnoldWeb"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["ArnoldWeb"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function NewtonLiepnik end
originalcode(::typeof(NewtonLiepnik)) = """
class NewtonLiepnik(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b):
        xdot = -a * x + y + 10 * y * z
        ydot = -x - 0.4 * y + 5 * x * z
        zdot = b * z - 5 * x * y
        return xdot, ydot, zdot
"""
@doc make_docstring(NewtonLiepnik) NewtonLiepnik
function NewtonLiepnik()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b = p
        du[1] = -a * x + y + 10 * y * z
        du[2] = -x - 0.4 * y + 5 * x * z
        du[3] = b * z - 5 * x * y
    end
    u0 = Float64.(ATTRACTOR_DATA["NewtonLiepnik"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["NewtonLiepnik"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HyperRossler end
originalcode(::typeof(HyperRossler)) = """
class HyperRossler(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a, b, c, d):
        xdot = -y - z
        ydot = x + a * y + w
        zdot = b + x * z
        wdot = -c * z + d * w
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(HyperRossler) HyperRossler
function HyperRossler()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d = p
        du[1] = -y - z
        du[2] = x + a * y + w
        du[3] = b + x * z
        du[4] = -c * z + d * w
    end
    u0 = Float64.(ATTRACTOR_DATA["HyperRossler"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HyperRossler"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HyperLorenz end
originalcode(::typeof(HyperLorenz)) = """
class HyperLorenz(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a, b, c, d):
        xdot = a * y - a * x + w
        ydot = -x * z + c * x - y
        zdot = -b * z + x * y
        wdot = d * w - x * z
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(HyperLorenz) HyperLorenz
function HyperLorenz()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d = p
        du[1] = a * y - a * x + w
        du[2] = -x * z + c * x - y
        du[3] = -b * z + x * y
        du[4] = d * w - x * z
    end
    u0 = Float64.(ATTRACTOR_DATA["HyperLorenz"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HyperLorenz"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HyperCai end
originalcode(::typeof(HyperCai)) = """
class HyperCai(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a, b, c, d, e):
        xdot = a * y - a * x
        ydot = b * x + c * y - x * z + w
        zdot = -d * z + y ** 2
        wdot = -e * x
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(HyperCai) HyperCai
function HyperCai()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d, e = p
        du[1] = a * y - a * x
        du[2] = b * x + c * y - x * z + w
        du[3] = -d * z + y^2
        du[4] = -e * x
    end
    u0 = Float64.(ATTRACTOR_DATA["HyperCai"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HyperCai"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HyperBao end
originalcode(::typeof(HyperBao)) = """
class HyperBao(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a, b, c, d, e):
        xdot = a * y - a * x + w
        ydot = c * y - x * z
        zdot = x * y - b * z
        wdot = e * x + d * y * z
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(HyperBao) HyperBao
function HyperBao()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d, e = p
        du[1] = a * y - a * x + w
        du[2] = c * y - x * z
        du[3] = x * y - b * z
        du[4] = e * x + d * y * z
    end
    u0 = Float64.(ATTRACTOR_DATA["HyperBao"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HyperBao"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HyperJha end
originalcode(::typeof(HyperJha)) = """
class HyperJha(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a, b, c, d):
        xdot = a * y - a * x + w
        ydot = -x * z + b * x - y
        zdot = x * y - c * z
        wdot = -x * z + d * w
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(HyperJha) HyperJha
function HyperJha()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d = p
        du[1] = a * y - a * x + w
        du[2] = -x * z + b * x - y
        du[3] = x * y - c * z
        du[4] = -x * z + d * w
    end
    u0 = Float64.(ATTRACTOR_DATA["HyperJha"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HyperJha"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HyperQi end
originalcode(::typeof(HyperQi)) = """
class HyperQi(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a, b, c, d, e, f):
        xdot = a * y - a * x + y * z
        ydot = b * x + b * y - x * z
        zdot = -c * z - e * w + x * y
        wdot = -d * w + f * z + x * y
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(HyperQi) HyperQi
function HyperQi()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d, e, f = p
        du[1] = a * y - a * x + y * z
        du[2] = b * x + b * y - x * z
        du[3] = -c * z - e * w + x * y
        du[4] = -d * w + f * z + x * y
    end
    u0 = Float64.(ATTRACTOR_DATA["HyperQi"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HyperQi"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Qi end
originalcode(::typeof(Qi)) = """
class Qi(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a, b, c, d):
        xdot = a * y - a * x + y * z * w
        ydot = b * x + b * y - x * z * w
        zdot = -c * z + x * y * w
        wdot = -d * w + x * y * z
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(Qi) Qi
function Qi()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d = p
        du[1] = a * y - a * x + y * z * w
        du[2] = b * x + b * y - x * z * w
        du[3] = -c * z + x * y * w
        du[4] = -d * w + x * y * z
    end
    u0 = Float64.(ATTRACTOR_DATA["Qi"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Qi"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function LorenzStenflo end
originalcode(::typeof(LorenzStenflo)) = """
class LorenzStenflo(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a, b, c, d):
        xdot = a * y - a * x + d * w
        ydot = c * x - x * z - y
        zdot = x * y - b * z
        wdot = -x - a * w
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(LorenzStenflo) LorenzStenflo
function LorenzStenflo()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d = p
        du[1] = a * y - a * x + d * w
        du[2] = c * x - x * z - y
        du[3] = x * y - b * z
        du[4] = -x - a * w
    end
    u0 = Float64.(ATTRACTOR_DATA["LorenzStenflo"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["LorenzStenflo"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HyperYangChen end
originalcode(::typeof(HyperYangChen)) = """
class HyperYangChen(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a=30, b=3, c=35, d=8):
        xdot = a * y - a * x
        ydot = c * x - x * z + w
        zdot = -b * z + x * y
        wdot = -d * x
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(HyperYangChen) HyperYangChen
function HyperYangChen()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d = p
        du[1] = a * y - a * x
        du[2] = c * x - x * z + w
        du[3] = -b * z + x * y
        du[4] = -d * x
    end
    u0 = Float64.(ATTRACTOR_DATA["HyperYangChen"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HyperYangChen"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HyperYan end
originalcode(::typeof(HyperYan)) = """
class HyperYan(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a=37, b=3, c=26, d=38):
        xdot = a * y - a * x
        ydot = (c - a) * x - x * z + c * y
        zdot = -b * z + x * y - y * z + x * z - w
        wdot = -d * w + y * z - x * z
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(HyperYan) HyperYan
function HyperYan()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d = p
        du[1] = a * y - a * x
        du[2] = (c - a) * x - x * z + c * y
        du[3] = -b * z + x * y - y * z + x * z - w
        du[4] = -d * w + y * z - x * z
    end
    u0 = Float64.(ATTRACTOR_DATA["HyperYan"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HyperYan"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HyperXu end
originalcode(::typeof(HyperXu)) = """
class HyperXu(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a=10, b=40, c=2.5, d=2, e=16):
        xdot = a * y - a * x + w
        ydot = b * x + e * x * z
        zdot = -c * z - x * y
        wdot = x * z - d * y
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(HyperXu) HyperXu
function HyperXu()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d, e = p
        du[1] = a * y - a * x + w
        du[2] = b * x + e * x * z
        du[3] = -c * z - x * y
        du[4] = x * z - d * y
    end
    u0 = Float64.(ATTRACTOR_DATA["HyperXu"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HyperXu"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HyperWang end
originalcode(::typeof(HyperWang)) = """
class HyperWang(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a=10, b=40, c=2.5, d=10.6, e=4):
        xdot = a * y - a * x
        ydot = -x * z + b * x + w
        zdot = -c * z + e * x ** 2
        wdot = -d * x
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(HyperWang) HyperWang
function HyperWang()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d, e = p
        du[1] = a * y - a * x
        du[2] = -x * z + b * x + w
        du[3] = -c * z + e * x^2
        du[4] = -d * x
    end
    u0 = Float64.(ATTRACTOR_DATA["HyperWang"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HyperWang"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HyperPang end
originalcode(::typeof(HyperPang)) = """
class HyperPang(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a=36, b=3, c=20, d=2):
        xdot = a * y - a * x
        ydot = -x * z + c * y + w
        zdot = x * y - b * z
        wdot = -d * x - d * y
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(HyperPang) HyperPang
function HyperPang()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d = p
        du[1] = a * y - a * x
        du[2] = -x * z + c * y + w
        du[3] = x * y - b * z
        du[4] = -d * x - d * y
    end
    u0 = Float64.(ATTRACTOR_DATA["HyperPang"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HyperPang"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HyperLu end
originalcode(::typeof(HyperLu)) = """
class HyperLu(DynSys):
    @staticjit
    def _rhs(x, y, z, w, t, a=36, b=3, c=20, d=1.3):
        xdot = a * y - a * x + w
        ydot = -x * z + c * y
        zdot = x * y - b * z
        wdot = d * w + x * z
        return xdot, ydot, zdot, wdot
"""
@doc make_docstring(HyperLu) HyperLu
function HyperLu()
    function rhs(du, u, p, t)
        x, y, z, w = u
        @unpack a, b, c, d = p
        du[1] = a * y - a * x + w
        du[2] = -x * z + c * y
        du[3] = x * y - b * z
        du[4] = d * w + x * z
    end
    u0 = Float64.(ATTRACTOR_DATA["HyperLu"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HyperLu"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function SaltonSea end
originalcode(::typeof(SaltonSea)) = """
class SaltonSea(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, d, k, lam, m, mu, r, th):
        xdot = r * x * (1 - (x + y) / k) - lam * x * y
        ydot = lam * x * y - m * y * z / (y + a) - mu * y
        zdot = th * y * z / (y + a) - d * z
        return xdot, ydot, zdot
"""
@doc make_docstring(SaltonSea) SaltonSea
function SaltonSea()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, d, k, lam, m, mu, r, th = p
        du[1] = r * x * (1 - (x + y) / k) - lam * x * y
        du[2] = lam * x * y - m * y * z / (y + a) - mu * y
        du[3] = th * y * z / (y + a) - d * z
    end
    u0 = Float64.(ATTRACTOR_DATA["SaltonSea"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["SaltonSea"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function ExcitableCell end
originalcode(::typeof(ExcitableCell)) = """
class ExcitableCell(DynSys):
    def rhs(self, X, t):
        v, n, c = X

        alpham = 0.1 * (25 + v) / (1 - np.exp(-0.1 * v - 2.5))
        betam = 4 * np.exp(-(v + 50) / 18)
        minf = alpham / (alpham + betam)

        alphah = 0.07 * np.exp(-0.05 * v - 2.5)
        betah = 1 / (1 + np.exp(-0.1 * v - 2))
        hinf = alphah / (alphah + betah)

        alphan = 0.01 * (20 + v) / (1 - np.exp(-0.1 * v - 2))
        betan = 0.125 * np.exp(-(v + 30) / 80)
        ninf = alphan / (alphan + betan)
        tau = 1 / (230 * (alphan + betan))

        ca = c / (1 + c)

        vdot = (
            self.gi * minf ** 3 * hinf * (self.vi - v)
            + self.gkv * n ** 4 * (self.vk - v)
            + self.gkc * ca * (self.vk - v)
            + self.gl * (self.vl - v)
        )
        ndot = (ninf - n) / tau
        cdot = self.rho * (minf ** 3 * hinf * (self.vc - v) - self.kc * c)
        return vdot, ndot, cdot
"""
@doc make_docstring(ExcitableCell) ExcitableCell
function ExcitableCell()
    function rhs(du, u, p, t)
        v, n, c = u
        @unpack gi, gkv, gkc, gl, rho, kc, vi, vk, vc, vl = p
        alpham = 0.1 * (25 + v) / (1 - exp(-0.1 * v - 2.5))
        betam = 4 * exp(-(v + 50) / 18)
        minf = alpham / (alpham + betam)

        alphah = 0.07 * exp(-0.05 * v - 2.5)
        betah = 1 / (1 + exp(-0.1 * v - 2))
        hinf = alphah / (alphah + betah)

        alphan = 0.01 * (20 + v) / (1 - exp(-0.1 * v - 2))
        betan = 0.125 * exp(-(v + 30) / 80)
        ninf = alphan / (alphan + betan)
        tau = 1 / (230 * (alphan + betan))

        ca = c / (1 + c)

        du[1] = gi * minf^3 * hinf * (vi - v) + gkv * n^4 * (vk - v) + gkc * ca * (vk - v) +
                gl * (vl - v)
        du[2] = (ninf - n) / tau
        du[3] = rho * (minf^3 * hinf * (vc - v) - kc * c)
    end
    u0 = Float64.(ATTRACTOR_DATA["ExcitableCell"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["ExcitableCell"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function CaTwoPlus end
originalcode(::typeof(CaTwoPlus)) = """
class CaTwoPlus(DynSys):
    def rhs(self, X, t):
        z, y, a = X
        Vin = self.V0 + self.V1 * self.beta
        V2 = self.Vm2 * (z ** 2) / (self.K2 ** 2 + z ** 2)
        V3 = (
            (self.Vm3 * (z ** self.m) / (self.Kz ** self.m + z ** self.m))
            * (y ** 2 / (self.Ky ** 2 + y ** 2))
            * (a ** 4 / (self.Ka ** 4 + a ** 4))
        )
        V5 = (
            self.Vm5
            * (a ** self.p / (self.K5 ** self.p + a ** self.p))
            * (z ** self.n / (self.Kd ** self.n + z ** self.n))
        )
        zdot = Vin - V2 + V3 + self.kf * y - self.k * z
        ydot = V2 - V3 - self.kf * y
        adot = self.beta * self.V4 - V5 - self.eps * a
        return (zdot, ydot, adot)
"""
@doc make_docstring(CaTwoPlus) CaTwoPlus
function CaTwoPlus()
    function rhs(du, u, p, t)
        z, y, a = u
        @unpack V0, V1, Vm2, Vm3, Vm5, V4, K2, Kz, Ky, Ka, K5, Kd, kf, k, beta, m, n, p, eps = p
        Vin = V0 + V1 * beta
        V2 = Vm2 * (z^2) / (K2^2 + z^2)
        V3 = (Vm3 * (z^m) / (Kz^m + z^m)) * (y^2 / (Ky^2 + y^2)) * (a^4 / (Ka^4 + a^4))
        V5 = Vm5 * (a^p / (K5^p + a^p)) * (z^n / (Kd^n + z^n))
        du[1] = Vin - V2 + V3 + kf * y - k * z
        du[2] = V2 - V3 - kf * y
        du[3] = beta * V4 - V5 - eps * a
    end
    u0 = Float64.(ATTRACTOR_DATA["CaTwoPlus"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["CaTwoPlus"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function CellCycle end
originalcode(::typeof(CellCycle)) = """
class CellCycle(DynSys):
    def rhs(self, X, t):
        c1, m1, x1, c2, m2, x2 = X
        Vm1, Um1 = 2 * [self.Vm1]
        vi1, vi2 = 2 * [self.vi]
        H1, H2, H3, H4 = 4 * [self.K]
        K1, K2, K3, K4 = 4 * [self.K]
        V2, U2 = 2 * [self.V2]
        Vm3, Um3 = 2 * [self.Vm3]
        V4, U4 = 2 * [self.V4]
        Kc1, Kc2 = 2 * [self.Kc]
        vd1, vd2 = 2 * [self.vd]
        Kd1, Kd2 = 2 * [self.Kd1]
        kd1, kd2 = 2 * [self.kd1]
        Kim1, Kim2 = 2 * [self.Kim]
        V1 = Vm1 * c1 / (Kc1 + c1)
        U1 = Um1 * c2 / (Kc2 + c2)
        V3 = m1 * Vm3
        U3 = m2 * Um3
        c1dot = vi1 * Kim1 / (Kim1 + m2) - vd1 * x1 * c1 / (Kd1 + c1) - kd1 * c1
        c2dot = vi2 * Kim2 / (Kim2 + m1) - vd2 * x2 * c2 / (Kd2 + c2) - kd2 * c2
        m1dot = V1 * (1 - m1) / (K1 + (1 - m1)) - V2 * m1 / (K2 + m1)
        m2dot = U1 * (1 - m2) / (H1 + (1 - m2)) - U2 * m2 / (H2 + m2)
        x1dot = V3 * (1 - x1) / (K3 + (1 - x1)) - V4 * x1 / (K4 + x1)
        x2dot = U3 * (1 - x2) / (H3 + (1 - x2)) - U4 * x2 / (H4 + x2)
        return c1dot, m1dot, x1dot, c2dot, m2dot, x2dot
"""
@doc make_docstring(CellCycle) CellCycle
function CellCycle()
    function rhs(du, u, p, t)
        c1, m1, x1, c2, m2, x2 = u
        @unpack V2, Vm1, V4, vi, Kim, vd, Kc, K, Kd1, kd1, Vm3 = p
        Vm1, Um1 = Vm1, Vm1
        vi1, vi2 = vi, vi
        H1, H2, H3, H4 = K, K, K, K
        K1, K2, K3, K4 = K, K, K, K
        V2, U2 = V2, V2
        Vm3, Um3 = Vm3, Vm3
        V4, U4 = V4, V4
        Kc1, Kc2 = Kc, Kc
        vd1, vd2 = vd, vd
        Kd1, Kd2 = Kd1, Kd1
        kd1, kd2 = kd1, kd1
        Kim1, Kim2 = Kim, Kim
        V1 = Vm1 * c1 / (Kc1 + c1)
        U1 = Um1 * c2 / (Kc2 + c2)
        V3 = m1 * Vm3
        U3 = m2 * Um3
        du[1] = vi1 * Kim1 / (Kim1 + m2) - vd1 * x1 * c1 / (Kd1 + c1) - kd1 * c1
        du[4] = vi2 * Kim2 / (Kim2 + m1) - vd2 * x2 * c2 / (Kd2 + c2) - kd2 * c2
        du[2] = V1 * (1 - m1) / (K1 + (1 - m1)) - V2 * m1 / (K2 + m1)
        du[5] = U1 * (1 - m2) / (H1 + (1 - m2)) - U2 * m2 / (H2 + m2)
        du[3] = V3 * (1 - x1) / (K3 + (1 - x1)) - V4 * x1 / (K4 + x1)
        du[6] = U3 * (1 - x2) / (H3 + (1 - x2)) - U4 * x2 / (H4 + x2)
    end
    u0 = Float64.(ATTRACTOR_DATA["CellCycle"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["CellCycle"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function CircadianRhythm end
originalcode(::typeof(CircadianRhythm)) = """
class CircadianRhythm(DynSys):
    @staticjit
    def _rhs(
        m,
        fc,
        fs,
        fn,
        th,
        t,
        Ki,
        k,
        k1,
        k2,
        kd,
        kdn,
        km,
        ks,
        n,
        vd,
        vdn,
        vm,
        vmax,
        vmin,
        v,
    ):
        vs = 2.5 * ((0.5 + 0.5 * np.cos(th)) + vmin) * (vmax - vmin)
        mdot = vs * (Ki ** n) / (Ki ** n + fn ** n) - vm * m / (km + m)
        fcdot = ks * m - k1 * fc + k2 * fn - k * fc
        fsdot = k * fc - vd * fs / (kd + fs)
        fndot = k1 * fc - k2 * fn - vdn * fn / (kdn + fn)
        thdot = 2 * np.pi / 24
        return mdot, fcdot, fsdot, fndot, thdot

    @staticjit
    def _postprocessing(m, fc, fs, fn, th):
        return m, fc, fs, fn, np.cos(th)
"""
@doc make_docstring(CircadianRhythm) CircadianRhythm
function CircadianRhythm()
    function rhs(du, u, p, t)
        m, fc, fs, fn, th = u
        @unpack Ki, k, k1, k2, kd, kdn, km, ks, n, vd, vdn, vm, vmax, vmin = p
        vs = 2.5 * ((0.5 + 0.5 * cos(th)) + vmin) * (vmax - vmin)
        du[1] = vs * (Ki^n) / (Ki^n + fn^n) - vm * m / (km + m)
        du[2] = ks * m - k1 * fc + k2 * fn - k * fc
        du[3] = k * fc - vd * fs / (kd + fs)
        du[4] = k1 * fc - k2 * fn - vdn * fn / (kdn + fn)
        du[5] = 2 * pi / 24
    end
    u0 = Float64.(ATTRACTOR_DATA["CircadianRhythm"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["CircadianRhythm"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function FluidTrampoline end
originalcode(::typeof(FluidTrampoline)) = """
class FluidTrampoline(DynSys):
    @staticmethod
    def _rhs(x, y, th, t, gamma, psi, w):
        xdot = y
        ydot = -1 - np.heaviside(-x, 0) * (x + psi * y * np.abs(y)) + gamma * np.cos(th)
        thdot = w
        return (xdot, ydot, thdot)

    @staticjit
    def _postprocessing(x, y, th):
        return x, y, np.cos(th)
"""
@doc make_docstring(FluidTrampoline) FluidTrampoline
function FluidTrampoline()
    function rhs(du, u, p, t)
        x, y, th = u
        @unpack gamma, psi, w = p
        du[1] = y
        du[2] = -1 - (x < 0) * (x + psi * y * abs(y)) + gamma * cos(th)
        du[3] = w
    end
    u0 = Float64.(ATTRACTOR_DATA["FluidTrampoline"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["FluidTrampoline"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Aizawa end
function originalcode(::typeof(Aizawa))
    """
class Aizawa(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c, d, e, f):
        xdot = x * z - b * x - d * y
        ydot = d * x + y * z - b * y
        zdot = c + a * z - 0.333333333333333333 * z ** 3 - x ** 2 - y ** 2 - e * z * x ** 2 - e * z * y ** 2 + f * z * x ** 3
        return xdot, ydot, zdot
"""
end
@doc make_docstring(Aizawa) Aizawa
function Aizawa()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c, d, e, f = p
        du[1] = x * z - b * x - d * y
        du[2] = d * x + y * z - b * y
        du[3] = c + a * z - 1 / 3 * z^3 - x^2 - y^2 - e * z * x^2 - e * z * y^2 +
                f * z * x^3
    end
    u0 = Float64.(ATTRACTOR_DATA["Aizawa"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Aizawa"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function AnishchenkoAstakhov end
originalcode(::typeof(AnishchenkoAstakhov)) = """
class AnishchenkoAstakhov(DynSys):
    def rhs(self, X, t):
        x, y, z = X
        mu, eta = self.mu, self.eta
        xdot = mu * x + y - x * z
        ydot = -x
        zdot = -eta * z + eta * np.heaviside(x, 0) * x ** 2
        return (xdot, ydot, zdot)
"""
@doc make_docstring(AnishchenkoAstakhov) AnishchenkoAstakhov
function AnishchenkoAstakhov()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack mu, eta = p
        du[1] = mu * x + y - x * z
        du[2] = -x
        du[3] = -eta * z + eta * (x > 0) * x^2
    end
    u0 = Float64.(ATTRACTOR_DATA["AnishchenkoAstakhov"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["AnishchenkoAstakhov"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function ShimizuMorioka end
originalcode(::typeof(ShimizuMorioka)) = """
class ShimizuMorioka(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b):
        xdot = y
        ydot = x - a * y - x * z
        zdot = -b * z + x ** 2
        return xdot, ydot, zdot
"""
@doc make_docstring(ShimizuMorioka) ShimizuMorioka
function ShimizuMorioka()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b = p
        du[1] = y
        du[2] = x - a * y - x * z
        du[3] = -b * z + x^2
    end
    u0 = Float64.(ATTRACTOR_DATA["ShimizuMorioka"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["ShimizuMorioka"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function GenesioTesi end
originalcode(::typeof(GenesioTesi)) = """
class GenesioTesi(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c):
        xdot = y
        ydot = z
        zdot = -c * x - b * y - a * z + x ** 2
        return xdot, ydot, zdot
"""
@doc make_docstring(GenesioTesi) GenesioTesi
function GenesioTesi()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c = p
        du[1] = y
        du[2] = z
        du[3] = -c * x - b * y - a * z + x^2
    end
    u0 = Float64.(ATTRACTOR_DATA["GenesioTesi"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["GenesioTesi"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function AtmosphericRegime end
originalcode(::typeof(AtmosphericRegime)) = """
class AtmosphericRegime(DynSys):
    @staticjit
    def _rhs(
        x, y, z, t, alpha, beta, mu1, mu2, omega, sigma
    ):
        xdot = mu1 * x + sigma * x * y
        ydot = mu2 * y + omega * z + alpha * y * z + beta * z ** 2 - sigma * x ** 2
        zdot = mu2 * z - omega * y - alpha * y ** 2 - beta * y * z
        return xdot, ydot, zdot
"""
@doc make_docstring(AtmosphericRegime) AtmosphericRegime
function AtmosphericRegime()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack alpha, beta, mu1, mu2, omega, sigma = p
        du[1] = mu1 * x + sigma * x * y
        du[2] = mu2 * y + omega * z + alpha * y * z + beta * z^2 - sigma * x^2
        du[3] = mu2 * z - omega * y - alpha * y^2 - beta * y * z
    end
    u0 = Float64.(ATTRACTOR_DATA["AtmosphericRegime"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["AtmosphericRegime"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Hadley end
originalcode(::typeof(Hadley)) = """
class Hadley(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, f, g):
        xdot = -y ** 2 - z ** 2 - a * x + a * f
        ydot = x * y - b * x * z - y + g
        zdot = b * x * y + x * z - z
        return xdot, ydot, zdot
"""
@doc make_docstring(Hadley) Hadley
function Hadley()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, f, g = p
        du[1] = -y^2 - z^2 - a * x + a * f
        du[2] = x * y - b * x * z - y + g
        du[3] = b * x * y + x * z - z
    end
    u0 = Float64.(ATTRACTOR_DATA["Hadley"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Hadley"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function ForcedVanDerPol end
originalcode(::typeof(ForcedVanDerPol)) = """
class ForcedVanDerPol(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, mu, w):
        ydot = mu * (1 - x ** 2) * y - x + a * np.sin(z)
        xdot = y
        zdot = w
        return xdot, ydot, zdot

    @staticjit
    def _postprocessing(x, y, z):
        return x, y, np.sin(z)
"""
@doc make_docstring(ForcedVanDerPol) ForcedVanDerPol
function ForcedVanDerPol()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, mu, w = p
        du[2] = mu * (1 - x^2) * y - x + a * sin(z)
        du[1] = y
        du[3] = w
    end
    u0 = Float64.(ATTRACTOR_DATA["ForcedVanDerPol"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["ForcedVanDerPol"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function ForcedFitzHughNagumo end
originalcode(::typeof(ForcedFitzHughNagumo)) = """
class ForcedFitzHughNagumo(DynSys):
    @staticjit
    def _rhs(v, w, z, t, a, b, curr, f, gamma, omega):
        vdot = v - v ** 3 / 3 - w + curr + f * np.sin(z)
        wdot = gamma * (v + a - b * w)
        zdot = omega
        return vdot, wdot, zdot

    @staticjit
    def _postprocessing(x, y, z):
        return x, y, np.sin(z)
"""
@doc make_docstring(ForcedFitzHughNagumo) ForcedFitzHughNagumo
function ForcedFitzHughNagumo()
    function rhs(du, u, p, t)
        v, w, z = u
        @unpack a, b, curr, f, gamma, omega = p
        du[1] = v - v^3 / 3 - w + curr + f * sin(z)
        du[2] = gamma * (v + a - b * w)
        du[3] = omega
    end
    u0 = Float64.(ATTRACTOR_DATA["ForcedFitzHughNagumo"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["ForcedFitzHughNagumo"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HindmarshRose end
originalcode(::typeof(HindmarshRose)) = """
class HindmarshRose(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c, d, s, tx, tz):
        xdot = -x + 1 / tx * y - a / tx * x ** 3 + b / tx * x ** 2 + 1 / tx * z
        ydot = -a * x ** 3 - (d - b) * x ** 2 + z
        zdot = -s / tz * x - 1 / tz * z + c / tz
        return xdot, ydot, zdot
"""
@doc make_docstring(HindmarshRose) HindmarshRose
function HindmarshRose()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c, d, s, tx, tz = p
        du[1] = -x + 1 / tx * y - a / tx * x^3 + b / tx * x^2 + 1 / tx * z
        du[2] = -a * x^3 - (d - b) * x^2 + z
        du[3] = -s / tz * x - 1 / tz * z + c / tz
    end
    u0 = Float64.(ATTRACTOR_DATA["HindmarshRose"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HindmarshRose"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Colpitts end
originalcode(::typeof(Colpitts)) = """
class Colpitts(DynSys):
    def rhs(self, X, t):
        x, y, z = X
        u = z - (self.e - 1)
        fz = -u * (1 - np.heaviside(u, 0))
        xdot = y - self.a * fz
        ydot = self.c - x - self.b * y - z
        zdot = y - self.d * z
        return (xdot, ydot, zdot)
"""
@doc make_docstring(Colpitts) Colpitts
function Colpitts()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c, d, e = p
        u = z - (e - 1)
        fz = -u * (1 - (u > 0))
        du[1] = y - a * fz
        du[2] = c - x - b * y - z
        du[3] = y - d * z
    end
    u0 = Float64.(ATTRACTOR_DATA["Colpitts"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Colpitts"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Laser end
originalcode(::typeof(Laser)) = """
class Laser(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, b, c, d, h, k):
        xdot = a * y - a * x + b * y * z ** 2
        ydot = c * x + d * x * z ** 2
        zdot = h * z + k * x ** 2
        return xdot, ydot, zdot
"""
@doc make_docstring(Laser) Laser
function Laser()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c, d, h, k = p
        du[1] = a * y - a * x + b * y * z^2
        du[2] = c * x + d * x * z^2
        du[3] = h * z + k * x^2
    end
    u0 = Float64.(ATTRACTOR_DATA["Laser"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Laser"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Blasius end
originalcode(::typeof(Blasius)) = """
class Blasius(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, alpha1, alpha2, b, c, k1, k2, zs):
        xdot = a * x - alpha1 * x * y / (1 + k1 * x)
        ydot = -b * y + alpha1 * x * y / (1 + k1 * x) - alpha2 * y * z / (1 + k2 * y)
        zdot = -c * (z - zs) + alpha2 * y * z / (1 + k2 * y)
        return xdot, ydot, zdot
"""
@doc make_docstring(Blasius) Blasius
function Blasius()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, alpha1, alpha2, b, c, k1, k2, zs = p
        du[1] = a * x - alpha1 * x * y / (1 + k1 * x)
        du[2] = -b * y + alpha1 * x * y / (1 + k1 * x) - alpha2 * y * z / (1 + k2 * y)
        du[3] = -c * (z - zs) + alpha2 * y * z / (1 + k2 * y)
    end
    u0 = Float64.(ATTRACTOR_DATA["Blasius"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Blasius"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function TurchinHanski end
originalcode(::typeof(TurchinHanski)) = """
class TurchinHanski(DynSys):
    @staticjit
    def _rhs(n, p, z, t, a, d, e, g, h, r, s):
        ndot = (
            r * (1 - e * np.sin(z)) * n
            - r * (n ** 2)
            - g * (n ** 2) / (n ** 2 + h ** 2)
            - a * n * p / (n + d)
        )
        pdot = s * (1 - e * np.sin(z)) * p - s * (p ** 2) / n
        zdot = 2 * np.pi
        return ndot, pdot, zdot

    @staticjit
    def _postprocessing(x, y, z):
        return x, y, np.sin(z)
"""
@doc make_docstring(TurchinHanski) TurchinHanski
function TurchinHanski()
    function rhs(du, u, p, t)
        @unpack a, d, e, g, h, r, s = p
        n, p, z = u
        du[1] = r * (1 - e * sin(z)) * n - r * n^2 - g * n^2 / (n^2 + h^2) -
                a * n * p / (n + d)
        du[2] = s * (1 - e * sin(z)) * p - s * p^2 / n
        du[3] = 2 * pi
    end
    u0 = Float64.(ATTRACTOR_DATA["TurchinHanski"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["TurchinHanski"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function StickSlipOscillator end
originalcode(::typeof(StickSlipOscillator)) = """
class StickSlipOscillator(DynSys):
    def _t(self, v):
        return self.t0 * np.sign(v) - self.alpha * v + self.beta * v ** 3

    @staticjit
    def _rhs(x, v, th, t, a, alpha, b, beta, eps, gamma, t0, vs, w):
        tq = t0 * np.sign(v - vs) - alpha * v + beta * (v - vs) ** 3
        xdot = v
        vdot = eps * (gamma * np.cos(th) - tq) + a * x - b * x ** 3
        thdot = w
        return xdot, vdot, thdot

    @staticjit
    def _postprocessing(x, v, th):
        return x, v, np.cos(th)
"""
@doc make_docstring(StickSlipOscillator) StickSlipOscillator
function StickSlipOscillator()
    function rhs(du, u, p, t)
        x, v, th = u
        @unpack a, alpha, b, beta, eps, gamma, t0, vs, w = p
        tq = t0 * sign(v - vs) - alpha * v + beta * (v - vs)^3
        du[1] = v
        du[2] = eps * (gamma * cos(th) - tq) + a * x - b * x^3
        du[3] = w
    end
    u0 = Float64.(ATTRACTOR_DATA["StickSlipOscillator"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["StickSlipOscillator"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function HastingsPowell end
originalcode(::typeof(HastingsPowell)) = """
class HastingsPowell(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a1, a2, b1, b2, d1, d2):
        xdot = x * (1 - x) - y * a1 * x / (1 + b1 * x)
        ydot = y * a1 * x / (1 + b1 * x) - z * a2 * y / (1 + b2 * y) - d1 * y
        zdot = z * a2 * y / (1 + b2 * y) - d2 * z
        return xdot, ydot, zdot
"""
@doc make_docstring(HastingsPowell) HastingsPowell
function HastingsPowell()
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a1, a2, b1, b2, d1, d2 = p
        du[1] = x * (1 - x) - y * a1 * x / (1 + b1 * x)
        du[2] = y * a1 * x / (1 + b1 * x) - z * a2 * y / (1 + b2 * y) - d1 * y
        du[3] = z * a2 * y / (1 + b2 * y) - d2 * z
    end
    u0 = Float64.(ATTRACTOR_DATA["HastingsPowell"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["HastingsPowell"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function CellularNeuralNetwork end
originalcode(::typeof(CellularNeuralNetwork)) = """
class CellularNeuralNetwork(DynSys):
    @staticjit
    def f(x):
        return 0.5 * (np.abs(x + 1) - np.abs(x - 1))

    def rhs(self, X, t):
        x, y, z = X
        xdot = -x + self.d * self.f(x) - self.b * self.f(y) - self.b * self.f(z)
        ydot = -y - self.b * self.f(x) + self.c * self.f(y) - self.a * self.f(z)
        zdot = -z - self.b * self.f(x) + self.a * self.f(y) + self.f(z)
        return (xdot, ydot, zdot)
"""
@doc make_docstring(CellularNeuralNetwork) CellularNeuralNetwork
function CellularNeuralNetwork()
    function _f(x)
        return 0.5 * (abs(x + 1) - abs(x - 1))
    end
    function rhs(du, u, p, t)
        x, y, z = u
        @unpack a, b, c, d = p
        du[1] = -x + d * _f(x) - b * _f(y) - b * _f(z)
        du[2] = -y - b * _f(x) + c * _f(y) - a * _f(z)
        du[3] = -z - b * _f(x) + a * _f(y) + _f(z)
    end
    u0 = Float64.(ATTRACTOR_DATA["CellularNeuralNetwork"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["CellularNeuralNetwork"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function BeerRNN end
originalcode(::typeof(BeerRNN)) = """
class BeerRNN(DynSys):
    @staticjit
    def _sig(x):
        return 1.0 / (1.0 + np.exp(-x))

    def rhs(self, X, t):
        Xdot = (-X + np.matmul(self.w, self._sig(X + self.theta))) / self.tau
        return Xdot
"""
@doc make_docstring(BeerRNN) BeerRNN
function BeerRNN()
    function sig(x)
        return 1.0 / (1.0 + exp(-x))
    end
    function rhs(du, u, p, t)
        @unpack w, theta, tau = p
        du .= (-u .+ w * sig.(u .+ theta)) ./ tau
    end
    u0 = Float64.(ATTRACTOR_DATA["BeerRNN"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["BeerRNN"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function Torus end
originalcode(::typeof(Torus)) = """
class Torus(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a, n, r):
        xdot = (-a * n * np.sin(n * t)) * np.cos(t) - (r + a * np.cos(n * t)) * np.sin(
            t
        )
        ydot = (-a * n * np.sin(n * t)) * np.sin(t) + (r + a * np.cos(n * t)) * np.cos(
            t
        )
        zdot = a * n * np.cos(n * t)
        return xdot, ydot, zdot
"""
@doc make_docstring(Torus) Torus
function Torus()
    function rhs(du, u, p, t)
        @unpack a, n, r = p
        du[1] = (-a * n * sin(n * t)) * cos(t) - (r + a * cos(n * t)) * sin(t)
        du[2] = (-a * n * sin(n * t)) * sin(t) + (r + a * cos(n * t)) * cos(t)
        du[3] = a * n * cos(n * t)
    end
    u0 = Float64.(ATTRACTOR_DATA["Torus"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Torus"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

# class CaTwoPlusQuasiperiodic(CaTwoPlus):
#     pass

function Hopfield end
originalcode(::typeof(Hopfield)) = """
class Hopfield(DynSys):
    def f(self, x):
        return (1 + np.tanh(x)) / 2

    def rhs(self, X, t):
        Xdot = -X / self.tau + self.f(self.eps * np.matmul(self.k, X)) - self.beta
        return Xdot
"""
@doc make_docstring(Hopfield) Hopfield
function Hopfield()
    function _f(x)
        return (1 + tanh(x)) / 2
    end
    function rhs(du, u, p, t)
        @unpack k, tau, eps, beta = p
        du .= -u ./ tau .+ _f.(eps .* k * u) .- beta
    end
    u0 = Float64.(ATTRACTOR_DATA["Hopfield"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["Hopfield"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function MacArthur end
originalcode(::typeof(MacArthur)) = """
class MacArthur(DynSys):
    def growth_rate(self, rr):
        u0 = rr / (self.k.T + rr)
        u = self.r * u0.T
        return np.min(u.T, axis=1)

    def rhs(self, X, t):
        nn, rr = X[:5], X[5:]
        mu = self.growth_rate(rr)
        nndot = nn * (mu - self.m)
        rrdot = self.d * (self.s - rr) - np.matmul(self.c, (mu * nn))
        return np.hstack([nndot, rrdot])
"""
@doc make_docstring(MacArthur) MacArthur
function MacArthur()
    function growth_rate(rr, k, r)
        u0 = rr ./ (k' .+ rr)
        u = r .* u0'
        return minimum(u, dims=2)
    end
    function rhs(du, u, p, t)
        @views nn, rr = u[1:5], u[6:10]
        @unpack m, r, k, d, s, c = p
        mu = growth_rate(rr, k, r)
        du[1:5] .= nn .* (mu .- m)
        du[6:10] .= d .* (s .- rr) .- c * (mu .* nn)
    end
    u0 = Float64.(ATTRACTOR_DATA["MacArthur"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["MacArthur"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function ItikBanksTumor end
originalcode(::typeof(ItikBanksTumor)) = """
class ItikBanksTumor(DynSys):
    @staticjit
    def _rhs(x, y, z, t, a12, a13, a21, a31, d3, k3, r2, r3):
        xdot = x * (1 - x) - a12 * x * y - a13 * x * z
        ydot = r2 * y * (1 - y) - a21 * x * y
        zdot = r3 * x * z / (x + k3) - a31 * x * z - d3 * z
        return xdot, ydot, zdot
"""
@doc make_docstring(ItikBanksTumor) ItikBanksTumor
function ItikBanksTumor()
    function rhs(du, u, p, t)
        @unpack a12, a13, a21, a31, d3, k3, r2, r3 = p
        x, y, z = u
        du[1] = x * (1 - x) - a12 * x * y - a13 * x * z
        du[2] = r2 * y * (1 - y) - a21 * x * y
        du[3] = r3 * x * z / (x + k3) - a31 * x * z - d3 * z
    end
    u0 = Float64.(ATTRACTOR_DATA["ItikBanksTumor"]["initial_conditions"])
    p = format_parameters(ATTRACTOR_DATA["ItikBanksTumor"]["parameters"])
    tspan = (0.0, 1.0)
    f = ODEFunction(rhs)
    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

# ## Doesn't match described dynamics
# # class CosmologyFriedmann(DynSys):
# #     @staticjit
# #     def _rhs(x, y, z, t, a, b, c, d, p):
# #         xdot = y
# #         ydot = -a * y**2 / x - b * x - c * x**3 + d * p * x
# #         zdot = 3 * (y / x) * (p + z)
# #         return xdot, ydot, zdot

# ## Doesn't match described dynamics
# # class MixMasterUniverse(DynSys):
# #     def rhs(self, X, t):
# #         a, b, g, adot_, bdot_, gdot_ = X
# #         adot = adot_
# #         bdot = bdot_
# #         gdot = gdot_
# #         addot = (np.exp(2*b) - np.exp(2*g))**2 - np.exp(4*a)
# #         bddot = (np.exp(2*g) - np.exp(2*a))**2 - np.exp(4*b)
# #         gddot = (np.exp(2*a) - np.exp(2*b))**2 - np.exp(4*g)
# #         return (adot, bdot, gdot, addot, bddot, gddot)

# ## Doesn't match described dynamics
# # class Universe(DynSys):
# #     def rhs(self, X, t):
# #         #Xdot = X * np.matmul(self.a, 1 - X)
# #         Xdot = self.r * X * (1 - np.matmul(self.a, X))
# #         return Xdot

# # class SeasonalSEIR:
# #     """
# #     This is extremely unstable for some reason
# #     Olsen, Schaffer. Science 1990

# #     eq = SeasonalSEIR()
# #     tpts = np.linspace(0, 1000, 100000)
# #     ic = (.579, .02, .001, 1e-6)
# #     sol = integrate_dyn(eq, ic, tpts)

# #     plt.plot(sol[0], sol[2], '.k', markersize=.1)

# #     """
# #     def __init__(self):
# #         pass

# #     def __call__(self, X, t):
# #         (s, e, i, th) = X

# #         ## measles
# #         m = 0.02
# #         a = 35.84
# #         g = 100
# #         b0 = 1800
# #         b1 = 0.28

# # #         ## chickenpox
# # #         m = 0.02
# # #         a = 36
# # #         g = 34.3
# # #         b0 = 537
# # #         b1 = 0.3

# #         b = b0*(1 + b1*np.cos(th))
# #         sdot = m*(1 - s) - b*s*i
# #         edot = b*s*i - (m + a)*e
# #         idot = a*e - (m + g)*i
# #         thdot = 2*np.pi
# #         return (sdot, edot, idot, thdot)

# # class SeasonalSEIR:
# #     """
# #     Seasonally forced SEIR model
# #     Zhang, Q., Liu, C., & Zhang, X. (2012). Analysis and Control of an SEIR Epidemic System with Nonlinear Transmission Rate. Lecture Notes in Control and Information Sciences, 203–225. doi:10.1007/978-1-4471-2303-3_14
# #     """
# #     def __init__(self, b=0.02, beta0=1800, alpha=35.84, gamma=100.0, beta1=0.28):
# #         self.b, self.beta0, self.alpha, self.gamma, self.beta1 = b, beta0, alpha, gamma, beta1
# #     def __call__(self, X, t):
# #         """
# #         The dynamical equation for the system
# #         - X : tuple corresponding to the three coordinates
# #         - t : float (the current time)
# #         alpha is a
# #         b is mu
# #         """
# #         s, e, i, th = X
# #         beta = self.beta0*(1 + self.beta1*np.cos(th)) # seasonal forcing
# #         sdot = self.b - self.b*s - beta*s*i
# #         edot = beta*s*i - (self.alpha + self.b)*e
# #         idot = self.alpha*e - (self.gamma + self.b)*i
# #         thdot = 2*np.pi
# #         return (sdot, edot, idot, thdot)

# # class SeasonalSEIR:
# #     """
# #     Seasonally forced SEIR model
# #     Zhang, Q., Liu, C., & Zhang, X. (2012). Analysis and Control of an SEIR Epidemic System with Nonlinear Transmission Rate. Lecture Notes in Control and Information Sciences, 203–225. doi:10.1007/978-1-4471-2303-3_14
# #     """
# #     def __init__(self, b=0.02, beta0=1800, alpha=35.84, gamma=100.0, beta1=0.28):
# #         self.b, self.beta0, self.alpha, self.gamma, self.beta1 = b, beta0, alpha, gamma, beta1
# #     def __call__(self, X, t):
# #         """
# #         The dynamical equation for the system
# #         - X : tuple corresponding to the three coordinates
# #         - t : float (the current time)
# #         alpha is a
# #         b is mu
# #         """
# #         s, e, i, th = X
# #         beta = self.beta0*(1 + self.beta1*np.cos(th)) # seasonal forcing
# #         sdot = self.b - self.b*s - beta*s*i
# #         edot = beta*s*i - (self.alpha + self.b)*e
# #         idot = self.alpha*e - (self.gamma + self.b)*i
# #         thdot = 2*np.pi
# #         return (sdot, edot, idot, thdot)

# # class SeasonalSEIR:
# #     """
# #     Seasonally forced SEIR model
# #     """
# #     def __init__(self, mu=0.02, b0=1800, a=35.84, g=100.0, d=0.28):
# #         self.mu, self.b0, self.a, self.g, self.d = mu, b0, a, g, d
# #     def __call__(self, X, t):
# #         """
# #         The dynamical equation for the system
# #         - X : tuple corresponding to the three coordinates
# #         - t : float (the current time)
# #         alpha is a
# #         b is mu
# #         """
# #         s, e, i, th = X
# #         b = self.b0*(1 + self.d*np.cos(th)) # seasonal forcing
# #         sdot = self.mu - self.mu*s - b*s*i
# #         edot = b*s*i - (self.mu + self.a)*e
# #         idot = self.a*e - (self.mu + self.g)*i
# # #         edot = 0
# # #         idot = 0
# #         thdot = 2*np.pi

# #         return (sdot, edot, idot, thdot)

# # class SeasonalSEIR:
# #     """
# #     Seasonally forced SEIR model
# # Bifurcation analysis of periodic SEIR and SIR epidemic models
# #     """
# #     def __init__(self, mu=0.02, b0=1884.95*5, a=35.842, g=100.0, d=0.255):
# #         self.mu, self.b0, self.a, self.g, self.d = mu, b0, a, g, d
# #     def __call__(self, X, t):
# #         """
# #         The dynamical equation for the system
# #         - X : tuple corresponding to the three coordinates
# #         - t : float (the current time)
# #         """
# #         s, e, i, th = X
# #         b = self.b0*(1 + self.d*np.cos(th)) # seasonal forcing
# #         sdot = self.mu - self.mu*s - b*s*i
# #         edot = b*s*i - (self.mu + self.a)*e
# #         idot = self.a*e - (self.mu + self.g)*i
# #         thdot = 2*np.pi
# #         return (sdot, edot, idot, thdot)

# # class SeasonalSIR:
# #     """
# #     Seasonally forced SEIR model
# #     """
# #     def __init__(self, mu=0.02, b0=1884.95, a=35.842, g=100, d=0.255):
# #         self.mu, self.b0, self.a, self.g, self.d = mu, b0, a, g, d
# #     def __call__(self, X, t):
# #         """
# #         The dynamical equation for the system
# #         - X : tuple corresponding to the three coordinates
# #         - t : float (the current time)
# #         """
# #         s, i, th = X
# #         b = self.b0*(1 + self.d*np.sin(th)) # seasonal forcing
# #         sdot = self.mu - self.mu*s - b*s*i
# #         idot =  b*s*i - (self.mu + self.g)*i
# #         thdot = 2*np.pi
# #         return (sdot, idot, thdot)

# # class Robinson:
# #     """
# #     C Robinson 1989 Nonlinearity
# #     Was unable to find published parameters for which this has a stable attractor,
# #     it may only be transient
# #     """
# #     def __init__(self, a=0.71, b=1.8587, v=1.0, gamma=0.7061, delta=0.1):
# #         self.a, self.b, self.v, self.gamma, self.delta = a, b, v, gamma, delta
# #     def __call__(self, X, t):
# #         """
# #         The dynamical equation for the system
# #         - X : tuple corresponding to the three coordinates
# #         - t : float (the current time)
# #         """
# #         x, y, z = X
# #         xdot = y
# #         ydot = x - 2*x**3 - self.a*y + (self.b*x**2)*y - self.v*y*z
# #         zdot = -self.gamma*z + self.delta*x**2
# #         return (xdot, ydot, zdot)

# # Sato: Cardiac model analogous to HH
# # https://www.sciencedirect.com/science/article/pii/S1007570416300016
# # hundreds of equations

# # class HodgkinHuxley:
# #     def __init__(self, i=7.92197):
# #         self.i = i

# #     def phi(self, x):
# #         return x/(np.exp(x) - 1)

# #     def __call__(self, X, t):
# #         """
# #         The dynamical equation for the system
# #         - X : tuple corresponding to the three coordinates
# #         - t : float (the current time)
# #         """
# #         v, m, n, h = X
# #         i = self.i
# #         vdot = i - (120*m**3*h*(v + 115) + 36*n**4*(v - 12) + 0.3*(v + 10.599))
# #         mdot = (1 - m)*self.phi((v + 25)/10) - m*(4*np.exp(v/18))
# #         ndot = (1 - n)*0.1*self.phi((v + 10)/10) - n*0.125*np.exp(v/80)
# #         hdot = (1 - h)*0.07*np.exp(v/20) - h/(1 + np.exp((v + 30)/10))
# #         return (vdot, mdot, ndot, hdot)

# ##############################
# ##
# ## Quasiperiodic systems
# ##
# ##############################

# # class SymmetricKuramoto(DynSys):
# #     def coupling(self, x, n=4):
# #         k = 1 + np.arange(n)
# #         return np.sum(self.a[:, None, None] * np.cos(k[:, None, None]*x[None, ...] - self.eta[:, None, None]), axis=0)
# #     def rhs(self, X, t):
# #         phase_diff = X[:, None] - X[None, :]
# #         Xdot = self.w + np.mean(self.coupling(phase_diff), axis=0) ## need to apply coupling element-wise
# #         return Xdot

# #     "SymmetricKuramoto": {
# #         "initial_conditions": [
# #             0.1,
# #             0.01,
# #             0.01,
# #             0.01
# #         ],
# #         "dt": 0.01,
# #         "parameters": {
# #             "w" : 0.0,
# #             "a": [-2, -2, -1,  -0.88],
# #             "eta": [0.1104, -0.1104, 0.669, 0.669]
# #         },
# #         "citation": "Bick, Christian, et al. Chaos in symmetric phase oscillator networks. Physical review letters 107.24 (2011): 244101.",
# #         "period": 7737.7
# #     }

# # class InterfacialFlight:
# #     """
# #     """
# #     def __init__(self):
# #         pass
# #     def __call__(self, X, t):
# #         x, y, z = X
# #         rclaw = 57 #um
# #         m = 2.2
# #         cl = 1.55
# #         ly = 73

# #         l0 = 137*1e6 # convert from uN/mg into um/s^2
# #         f = 116
# #         r = 0.15
# #         phi0 = np.pi/2

# #         cdleg = 3
# #         rhow = 1e-9 # water density in mg/um^3

# #         sigma = 72800 # water surface tension mg/s^2

# #         kinv = 2709 # um
# #         hclaw = rclaw*np.log(2*kinv/rclaw)

# #         sech = lambda pp : 1 / np.cosh(pp)

# #         phi_arg = 2*np.pi*f*t + phi0
# #         sin_term = np.sin(phi_arg)
# #         hsin_term = np.heaviside(sin_term, 0)
# #         zdot = (l0/m)*np.cos(phi_arg)*(hsin_term + r*(1 - hsin_term))
# #         zdot += -(8*np.pi*sigma*rclaw/m)*sech((x - hclaw)/rclaw)*np.sign(x)
# #         zdot += -(2*rhow*cdleg*np.pi*rclaw**2/m)*y*np.sign(y)

# #         xdot = 1
# #         ydot = x
# #         return (xdot, ydot, zdot)

# ## Not chaotic
# # class Robinson(DynSys):
# #     @staticjit
# #     def _rhs(x, y, z, t, a, b, c, d, v):
# #         xdot = y
# #         ydot = x - 2 * x**3 - a * y + b * x**2 * y - v * y * z
# #         zdot = -c * z + d * x**2
# #         return (xdot, ydot, zdot)
