"""Item 1: the affine change z -> al0*z + be0, Y -> al0*Y (k=2) is a ring isomorphism
A_Q -> A_Qt with Qt(z) = Q(al0 z + be0)/al0^k, and it commutes with the derivations. Checked symbolically
for generic Q of degree m = 1..5 and generic elements A(z) + B(z)Y."""
import sympy as sp
z, Y, al0, be0 = sp.symbols('z Y al0 be0')
def red(e, Q, k):
    e = sp.expand(e)
    return e if k == 1 else sp.expand(sp.rem(e, Y**2 - Q, Y))
def D(e, Q, k):
    if k == 1: return sp.expand(sp.diff(e, z)*Q)
    return red(sp.diff(e, z)*Y + sp.diff(e, Y)*sp.diff(Q, z)/2, Q, k)
ok = True
for k in (1, 2):
    for m in range(1, 6):
        q = sp.symbols(f'q0:{m+1}'); Q = sum(q[i]*z**i for i in range(m+1))
        Qt = sp.expand(Q.subs(z, al0*z+be0)/al0**k)
        sigma = lambda e: sp.expand(e.subs({z: al0*z+be0, Y: al0*Y}, simultaneous=True))
        a = sp.symbols('a0:4'); b = sp.symbols('b0:3')
        E = sum(a[i]*z**i for i in range(4)) + (sum(b[i]*z**i for i in range(3))*Y if k == 2 else 0)
        lhs = red(sigma(D(E, Q, k)), Qt, k)          # sigma(D E)
        rhs = D(sigma(E), Qt, k)                     # Dt(sigma E)
        good = sp.simplify(lhs - rhs) == 0
        ok &= good
        print(f'k={k} m={m}: sigma o D == Dt o sigma: {good}')
print('ALL OK' if ok else 'FAILURE')
