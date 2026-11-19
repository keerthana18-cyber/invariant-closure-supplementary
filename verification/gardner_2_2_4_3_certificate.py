"""Hand-checkable certificate for Gardner, (k,m)=(2,2), (d_A,d_B)=(4,3), where Singular's
lift did not return cofactors. The two top-degree equations already contradict the
exact-degree conditions a4, b3, q2, alpha != 0."""
import sympy as sp, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep import system, p
c, al, s = sp.symbols('c al s'); q2, q0 = sp.symbols('q2 q0')
a = sp.symbols('a0:5'); b = sp.symbols('b0:4'); a4, b3 = a[4], b[3]
f = lambda U: U[3]+6*U[0]*U[1]+6*al*U[0]**2*U[1]-c*U[1]
eqs = system(f, sum(a[i]*p**i for i in range(5)), sum(b[i]*p**i for i in range(4)), q2*p**2+q0, 2)
E0 = 24*al*b3*q2*(3*a4**2+b3**2*q2); E1 = 24*a4*al*(a4**2+3*b3**2*q2)
assert any(sp.expand(e-E0) == 0 for e in eqs) and any(sp.expand(e-E1) == 0 for e in eqs)
R = s*a4*b3*q2*al-1
c1 = s**3*a4**2*q2**2*al**2*3*b3/192; c0 = -s**3*a4**3*q2*al**2/192
h, r = sp.div(sp.expand(c1*E1+c0*E0-1), R, s); assert r == 0
print('c1*E1 + c0*E0 - h*R =', sp.expand(c1*E1 + c0*E0 - h*R))
