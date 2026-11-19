"""Numerical check (mpmath, 30 digits) of the Jacobi elliptic closed forms in the paper:
the KdV cnoidal profile and the Gardner sn^2 profile satisfy their auxiliary equations
(U')^2 = Q(U); by the first-integral proposition they then satisfy the reduced equations."""
import mpmath as mp
mp.mp.dps = 30
# KdV: (U')^2 = -2 (U - e1)(U - e2)(U - e3)
e1, e2, e3 = mp.mpf(-1), mp.mpf('0.3'), mp.mpf(2)
m = (e3-e2)/(e3-e1); a = mp.sqrt((e3-e1)/2)
U = lambda w: e2 + (e3-e2)*mp.ellipfun('cn', a*w, m=m)**2
for w in (mp.mpf('0.3'), mp.mpf('1.1')):
    print('KdV cnoidal residual  %.1e' % abs(mp.diff(U, w)**2 + 2*(U(w)-e1)*(U(w)-e2)*(U(w)-e3)))
# Gardner (alpha = -1, kappa = 1): (psi')^2 = q4 (psi-e1)(psi-e2)(psi-e3)(psi-e4)
e = [mp.mpf(x) for x in ('-1.5', '-0.2', '1.2', '2.5')]; q4 = 1
g = mp.sqrt(q4*(e[2]-e[0])*(e[3]-e[1]))/2
m = (e[2]-e[1])*(e[3]-e[0])/((e[2]-e[0])*(e[3]-e[1]))
def psi(w):
    s = mp.ellipfun('sn', g*w, m=m)**2
    return (s*(e[2]-e[1])*e[0]-(e[2]-e[0])*e[1])/(s*(e[2]-e[1])-(e[2]-e[0]))
Q = lambda z: q4*(z-e[0])*(z-e[1])*(z-e[2])*(z-e[3])
for w in (mp.mpf('0.3'), mp.mpf('1.1')):
    print('Gardner elliptic residual  %.1e' % abs(mp.diff(psi, w)**2 - Q(psi(w))))
