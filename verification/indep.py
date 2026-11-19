"""Independent re-derivation of the coefficient systems (does NOT import search.py).
Chain rule on symbols p = psi, y = psi'; reduction of y^2 by polynomial remainder."""
import sympy as sp
p, y = sp.symbols('p y')

def total_derivative(f, Q, k):
    # d/dw f(p, y) with psi' = y, and y' = Q'(p)/2 when k = 2; for k = 1, y is eliminated (y = Q(p)).
    if k == 1:
        return sp.expand(sp.diff(f, p)*Q)
    return sp.expand(sp.diff(f, p)*y + sp.diff(f, y)*sp.diff(Q, p)/2)

def reduce_y(expr, Q, k):
    if k == 1: return sp.expand(expr)
    return sp.expand(sp.rem(sp.expand(expr), y**2 - Q, y))

def system(residual_fn, A, B, Q, k):
    U = [A + (B*y if k == 2 else 0)]
    for _ in range(3):
        U.append(reduce_y(total_derivative(U[-1], Q, k), Q, k))
    R = reduce_y(residual_fn(U), Q, k)
    eqs = []
    if k == 1:
        eqs = sp.Poly(R, p).coeffs()
    else:
        P = sp.Poly(R, y)
        for j in range(2):
            eqs += sp.Poly(P.coeff_monomial(y**j), p).coeffs()
    return [e for e in eqs if e != 0]
