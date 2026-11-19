"""Singular cofactor certificates for the supplementary audit (Fisher and KdV-Burgers, 80 cases),
using the independent derivation of indep.py. Same exact-degree and nonzero conditions as scripts/box.py
and scripts/boxI.py: c != 0 (Fisher, Family II) and nu != 0 (KdV-Burgers, both families)."""
import sympy as sp, sys, json, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from indep import system, p
from singular import run
c, nu, s = sp.symbols('c nu s')
q = sp.symbols('q0:5'); g2, g3 = sp.symbols('g2 g3')
EQ = {'Fisher': (lambda U: U[2]+c*U[1]+U[0]-U[0]**2, [c]),
      'KdV-Burgers': (lambda U: U[3]-nu*U[2]+U[0]*U[1]-c*U[1], [c, nu])}
QN = {(1, 2): (p**2+q[0], [q[0]]), (2, 3): (4*p**3-g2*p-g3, [g2, g3]), (2, 4): (p**4+q[2]*p**2+q[1]*p+q[0], [q[2], q[1], q[0]])}
cases = [(e, 1, 2, dA, -1) for e in EQ for dA in range(1, 5)] + \
        [(e, 2, m, dA, dB) for e in EQ for m in (3, 4) for dA in range(-1, 4) for dB in range(-1, 3) if not (dB == -1 and dA <= 0)]
def cert_ok(G, cof): return cof is not None and sp.expand(sum(a*b for a, b in zip(cof, G)) - 1) == 0
out = (HERE / 'verify_supp.jsonl').open('w')
for e, k, m, dA, dB in cases:
    t = time.time(); f, pars = EQ[e]; Q, qv = QN[(k, m)]
    a = sp.symbols(f'a0:{dA+1}') if dA >= 0 else (); b = sp.symbols(f'b0:{dB+1}') if dB >= 0 else ()
    A = sp.Integer(0) + sum(a[i]*p**i for i in range(len(a))); B = sp.Integer(0) + sum(b[i]*p**i for i in range(len(b)))
    eqs = system(f, A, B, Q, k)
    extra = (c if (e == 'Fisher' and k == 2) else 1) * (nu if e == 'KdV-Burgers' else 1)
    nz = sp.prod(([a[-1]] if a else []) + ([b[-1]] if b else [])) * extra
    G = eqs + [s*nz - 1]
    inc, cof, Gs = run(G, list(qv)+list(a)+list(b)+pars+[s], True)
    st = ('CERT_VERIFIED' if cert_ok(Gs, cof) else 'CERT_FAILED_CHECK') if inc else 'EXISTS'
    out.write(json.dumps({'case': [e, k, m, dA, dB], 'status': st, 'seconds': round(time.time()-t, 2)})+'\n'); out.flush()
    print((e, k, m, dA, dB), st, round(time.time()-t, 1), flush=True)
