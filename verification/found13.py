"""Item 4: independent check of the 13 FOUND closures.
(a) exact residual with the independent derivation of indep.py (not search.py);
(b) numerical check: integrate psi'' = Q'(psi)/2 (k=2) or psi' = Q(psi) (k=1) with mpmath,
    form U = A(psi) + B(psi) psi' and evaluate the ORIGINAL reduced ODE by numerical differentiation."""
import json, random, sys, sympy as sp, mpmath as mpm
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
RESULTS = ROOT / 'results' / 'results612.jsonl'
OUTPUT = HERE / 'found13_check.json'
from indep import system, p
mpm.mp.dps = 50
z = sp.Symbol('z')
c, Lam, al = sp.symbols('c Lambda alpha')
EQ = {'KdV': lambda U: U[3]+6*U[0]*U[1]-c*U[1], 'mKdV': lambda U: U[3]-6*U[0]**2*U[1]-c*U[1],
      'Burgers': lambda U: U[2]-2*U[0]*U[1]+c*U[1], 'Boussinesq': lambda U: U[2]-Lam*U[0]+U[0]**2,
      'Gardner': lambda U: U[3]+6*U[0]*U[1]+6*al*U[0]**2*U[1]-c*U[1]}
q = sp.symbols('q0:6'); g2, g3 = sp.symbols('g2 g3')
QN = {(1,2): z**2+q[0], (2,3): 4*z**3-g2*z-g3, (2,4): z**4+q[2]*z**2+q[1]*z+q[0]}
random.seed(1)
rows = []
for l in RESULTS.open():
    r = json.loads(l)
    if r['status'] != 'FOUND': continue
    e, k, m, dA, dB = r['case']
    loc = {str(s): s for s in [c, Lam, al, g2, g3, *q, *sp.symbols('a0:6'), *sp.symbols('b0:6')]}
    so = {sp.Symbol(kk) if kk not in loc else loc[kk]: sp.sympify(v, locals=loc) for kk, v in r['detail'].items()}
    a = sp.symbols(f'a0:{dA+1}'); b = sp.symbols(f'b0:{dB+1}') if dB >= 0 else ()
    A = sp.Integer(0) + sum(a[i]*z**i for i in range(dA+1)); B = sp.Integer(0) + sum(b[i]*z**i for i in range(len(b)))
    Q = QN[(k, m)]
    expr = [A.subs(so), B.subs(so), Q.subs(so), c.subs(so), Lam.subs(so), al.subs(so)]
    free = sorted(set().union(*[x.free_symbols for x in expr]) - {z}, key=str)
    for attempt in range(50):
        val = {f: sp.Rational(random.randint(-9, 9), random.randint(1, 4)) for f in free}
        Av, Bv, Qv, cv, Lv, alv = [sp.expand(x.subs(val)) for x in expr]
        lead_ok = sp.Poly(Av, z).degree() == dA if dA >= 1 else True
        lead_ok = lead_ok and (dB < 0 or sp.Poly(Bv, z).degree() == dB) and Qv != 0 and sp.Poly(Qv, z).degree() == m
        if e == 'Gardner': lead_ok = lead_ok and alv != 0
        nonsq = True
        if k == 2:
            fac = sp.factor_list(Qv)[1]; nonsq = any(mult % 2 for _, mult in fac)
        if lead_ok and nonsq: break
    sub = {c: cv, Lam: Lv, al: alv}
    f = lambda U: EQ[e](U).subs(sub)
    # (a) exact, independent derivation
    exact = system(f, Av.subs(z, p), Bv.subs(z, p), Qv.subs(z, p), k)
    exact_ok = all(sp.simplify(x) == 0 for x in exact)
    # (b) numerical integration and original ODE
    Qf = sp.lambdify(z, Qv, 'mpmath'); dQf = sp.lambdify(z, sp.diff(Qv, z), 'mpmath')
    Af = sp.lambdify(z, Av, 'mpmath'); Bf = sp.lambdify(z, Bv, 'mpmath')
    psi0 = mpm.mpf('0.37')
    while abs(Qf(psi0)) < mpm.mpf('1e-3'): psi0 += mpm.mpf('0.11')
    if k == 1:
        sol = mpm.odefun(lambda w, Y: [Qf(Y[0])], 0, [psi0])
        U = lambda w: Af(sol(w)[0])
    else:
        sol = mpm.odefun(lambda w, Y: [Y[1], dQf(Y[0])/2], 0, [psi0, mpm.sqrt(Qf(psi0))])
        U = lambda w: Af(sol(w)[0]) + Bf(sol(w)[0])*sol(w)[1]
    res = []
    for w0 in [mpm.mpf('0.05'), mpm.mpf('0.15')]:
        d = [mpm.diff(U, w0, n, h=mpm.mpf('1e-10')) for n in range(4)]
        val_ = complex(sp.N(EQ[e](d).subs({c: cv, Lam: Lv, al: alv}), 30)) if False else None
        Ud = [sp.Float(str(mpm.re(x)), 35) + sp.I*sp.Float(str(mpm.im(x)), 35) for x in d]
        R = sp.N(EQ[e](Ud).subs(sub), 30)
        res.append(abs(complex(R))/ (1+max(abs(complex(x)) for x in Ud))**3)
    sc = max(abs(complex(sp.N(x))) for x in Ud) + 1
    rows.append((r['case'], {str(k_): str(v) for k_, v in val.items()}, str(sp.factor(Qv)), exact_ok, max(res), nonsq))
    print(r['case'], 'Q =', sp.factor(Qv), '| exact residual zero:', exact_ok, '| nonsquare:', nonsq, '| max numerical residual: %.1e' % max(res), flush=True)
json.dump(rows, OUTPUT.open('w'), indent=1)
print('written:', OUTPUT)
