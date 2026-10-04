"""Independent symbolic checks for the Maple focus/period algorithm.

Requires SymPy. This verifies the mathematics with exact closed-form
solutions; it does not execute or translate the Maple implementation.
"""

import sympy as sp

x, y, r, rho, theta = sp.symbols("x y r rho theta", real=True)
a, b = sp.symbols("a b", real=True)
c, s = sp.cos(theta), sp.sin(theta)


def equal(actual, expected):
    difference = sp.trigsimp(sp.cancel(sp.expand(actual - expected)))
    assert difference == 0, (actual, expected, difference)


def polar(field_x, field_y):
    replacement = {x: r * c, y: r * s}
    field_x = field_x.subs(replacement, simultaneous=True)
    field_y = field_y.subs(replacement, simultaneous=True)
    radial = sp.trigsimp(sp.expand(c * field_x + s * field_y))
    angular = sp.trigsimp(sp.expand((c * field_y - s * field_x) / r))
    return radial, angular


def powers_of_radius(radius_coefficients, order):
    powers = [[sp.S.Zero] * (order + 1) for _ in range(order + 1)]
    powers[0][0] = sp.S.One
    for exponent in range(1, order + 1):
        for degree in range(exponent, order + 1):
            powers[exponent][degree] = sp.expand(sum(
                powers[exponent - 1][degree - k] * radius_coefficients[k]
                for k in range(1, degree + 1)
            ))
    return powers


def recurrence(radial, angular, radius_coefficients, order=5, period_order=4):
    """Check reciprocal and radial recurrences against supplied exact u_n."""
    g_degree = sp.degree(radial, r)
    h_degree = sp.degree(angular, r)
    G = {j: sp.expand(radial).coeff(r, j) for j in range(2, g_degree + 1)}
    H = {j: sp.expand(angular).coeff(r, j) for j in range(h_degree + 1)}
    equal(H[0], -1)
    A = [1 / H[0]]
    for n in range(1, max(order - 2, period_order) + 1):
        A.append(sp.expand(-sum(H[j] * A[n - j]
                               for j in range(1, min(n, h_degree) + 1)) / H[0]))
    reciprocal = sum(A[n] * r**n for n in range(len(A)))
    product = sp.expand(angular * reciprocal)
    equal(product.coeff(r, 0), 1)
    for n in range(1, len(A)):
        equal(product.coeff(r, n), 0)
    powers = powers_of_radius(radius_coefficients, order)
    R = {j: sp.expand(sum(G[m] * A[j - m]
                         for m in range(2, min(j, g_degree) + 1)))
         for j in range(2, order + 1)}
    for n in range(2, order + 1):
        rhs = sum(R[j] * powers[j][n] for j in range(2, n + 1))
        equal(sp.diff(radius_coefficients[n], theta), rhs)
        equal(sp.sympify(radius_coefficients[n]).subs(theta, 0), 0)
    # P = - integral A(theta, r(theta,rho)) dtheta.
    time_coefficients = [sp.expand(-sum(A[j] * powers[j][n]
                                       for j in range(n + 1)))
                         for n in range(period_order + 1)]
    return A, R, time_coefficients


def trigonometric_integral(expression):
    """Exact integral on [0, 2 Pi] for a polynomial in sin(theta), cos(theta)."""
    cs, sn = sp.symbols("cs sn")
    polynomial = sp.Poly(sp.expand(expression).subs({c: cs, s: sn}), cs, sn)
    result = sp.S.Zero
    for (i, j), coefficient in polynomial.terms():
        if i % 2 == 0 and j % 2 == 0:
            moment = (2 * sp.pi * sp.factorial(i) * sp.factorial(j)
                      / (4**((i + j) // 2) * sp.factorial(i // 2)
                         * sp.factorial(j // 2) * sp.factorial((i + j) // 2)))
            result += coefficient * moment
    return sp.factor(result)


def check_focus():
    q = x*x + y*y
    field_x = y + a*x*q + b*y*q
    field_y = -x + a*y*q - b*x*q
    radial, angular = polar(field_x, field_y)
    equal(radial, a*r**3)
    equal(angular, -1 - b*r**2)
    radius = [0, 1, 0, -a*theta, 0,
              sp.Rational(3, 2)*a*a*theta*theta + a*b*theta]
    _, _, time = recurrence(radial, angular, radius)
    focus = [sp.simplify(sp.sympify(radius[n]).subs(theta, 2*sp.pi)/(2*sp.pi))
             for n in range(2, 6)]
    for actual, value in zip(focus, [0, -a, 0, a*b + 3*sp.pi*a*a]):
        equal(actual, value)
    period = [sp.integrate(time[n], (theta, 0, 2*sp.pi)) for n in range(1, 5)]
    expected = [0, -2*sp.pi*b, 0, 2*sp.pi*b*b + 4*sp.pi*sp.pi*a*b]
    for actual, value in zip(period, expected):
        equal(actual, value)
    # At a non-center, the positive-time clockwise passage has a different
    # coefficient. The prescribed P is a formal backward-passage quantity.
    clockwise_p4 = sp.integrate(-time[4], (theta, 0, -2*sp.pi))
    equal(clockwise_p4, 2*sp.pi*b*b - 4*sp.pi*sp.pi*a*b)
    print("FOCUS: G3=a; H0=-1; H2=-b")
    print("FOCUS: g2..g5 =", focus)
    print("FOCUS: formal p1..p4 =", expected)
    print("FOCUS: clockwise positive-time p4 =", clockwise_p4)


def check_exact_center():
    field_x = (1 + b*x)*(y + a*x*y)
    field_y = (1 + b*x)*(-x + a*y*y)
    radial, angular = polar(field_x, field_y)
    equal(radial, a*r*r*s + a*b*r**3*c*s)
    equal(angular, -1 - b*r*c)
    exact_radius = rho/(1 + a*rho*(1 - c))
    equal(sp.diff(exact_radius, theta),
          (radial/angular).subs(r, exact_radius))
    equal(exact_radius.subs(theta, 0), rho)
    equal(exact_radius.subs(theta, 2*sp.pi), rho)
    radius = [0] + [(-a)**(n - 1)*(1 - c)**(n - 1) for n in range(1, 6)]
    _, R, time = recurrence(radial, angular, radius)
    equal(R[2], -a*s)
    for n in range(3, 6):
        equal(R[n], 0)
    period = [trigonometric_integral(time[n]) for n in range(1, 5)]
    expected = [0, sp.pi*b*(b - a), -2*sp.pi*a*b*(b - a),
                3*sp.pi*b*(b - a)*(b*b - 2*a*b + 5*a*a)/4]
    for actual, value in zip(period, expected):
        equal(actual, value)
    exact_period = 2*sp.pi*(-a/(b - a) + b*(1 + a*rho)
                    / ((b - a)*sp.sqrt((1 + a*rho)**2 - (b - a)**2*rho*rho)))
    series = sp.series(exact_period, rho, 0, 5).removeO().expand()
    equal(series.coeff(rho, 0), 2*sp.pi)
    for n in range(1, 5):
        equal(series.coeff(rho, n), expected[n - 1])
    numeric_parameters = [sp.simplify(sp.sympify(value).subs({a: 1, b: 2}))
                          for value in expected]
    assert numeric_parameters == [0, 2*sp.pi, -4*sp.pi, sp.Rational(15, 2)*sp.pi]
    print("CENTER: G2=a*sin(theta); G3=a*b*cos(theta)*sin(theta)")
    print("CENTER: H0=-1; H1=-b*cos(theta); H2=0")
    print("CENTER: g_n=0 for all n (exact return)")
    print("CENTER: p1..p4 =", expected)
    print("CENTER a=1,b=2: p1..p4 =", numeric_parameters)


def check_zhang_p2():
    a1, a2, a3, a4, b1, b2, b3 = sp.symbols("a1 a2 a3 a4 b1 b2 b3")
    # Overall time reversal of equation (1) in Zhang, Hou and Zeng (2000).
    # This preserves each center orbit and its positive period, while matching
    # the requested clockwise linear part and keeping x=r*cos, y=r*sin.
    field_x = y - a1*x*x - a2*y*y - a3*x*x*y - a4*y**3
    field_y = -x - b1*x*y - b2*x**3 - b3*x*y*y
    radial, angular = polar(field_x, field_y)
    G2 = -(a1*c**3 + (a2 + b1)*s*s*c)
    G3 = -((a3 + b2)*c**3*s + (a4 + b3)*c*s**3)
    H1 = (a1 - b1)*c*c*s + a2*s**3
    H2 = (a3 - b3)*c*c*s*s - b2*c**4 + a4*s**4
    equal(radial, G2*r*r + G3*r**3)
    equal(angular, -1 + H1*r + H2*r*r)
    u2 = a1*s + (a2 + b1 - a1)*s**3/3
    equal(sp.diff(u2, theta), -G2)
    p2 = trigonometric_integral(H1*u2 + H1*H1 + H2)
    reference = sp.pi/12*(-5*a1*b1 + 4*a1*a1 - 3*b3 + 3*a3 + b1*b1
                          - a2*b1 + 10*a1*a2 - 9*b2 + 9*a4 + 10*a2*a2)
    equal(p2, reference)
    equal(p2.subs({a2: 0, a3: 0, a4: 0, b1: 0, b2: 0, b3: 0}), sp.pi*a1*a1/3)
    print("ZHANG (2000), page 774: full symbolic p2 agrees")
    print("ZHANG a1-only: p2=Pi*a1^2/3")


if __name__ == "__main__":
    check_focus()
    check_exact_center()
    check_zhang_p2()
    print("All independent symbolic checks passed.")
