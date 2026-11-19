"""Class B reproduction of the 612-case exact-degree audit.
No Stage 1 information is used: no pole order, no a0, no b, no (s1).
Affine normalisation of Q (psi -> alpha*psi + beta) is applied; it maps closures with given
(k, m, d_A, d_B) bijectively onto closures with the same data and preserves non-squareness."""
import sympy as sp, sys, os, json, time, multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from search import red, D
z, Y = sp.symbols('z Y')

def normQ(k, m):
    if k == 1 and m == 1: q1 = sp.Symbol('q1'); return q1*z, [q1], q1, []
    if k == 1 and m == 2: q0 = sp.Symbol('q0'); return z**2+q0, [q0], 1, []
    if k == 2 and m == 1: return z, [], 1, []
    if k == 2 and m == 2:
        q2, q0 = sp.symbols('q2 q0'); return q2*z**2+q0, [q2, q0], q2, [q0]
    if k == 2 and m == 3:
        g2, g3 = sp.symbols('g2 g3'); return 4*z**3-g2*z-g3, [g2, g3], 1, []
    if k == 2 and m == 4:
        q2, q1, q0 = sp.symbols('q2 q1 q0')
        return z**4+q2*z**2+q1*z+q0, [q2, q1, q0], 1, [q1, 4*q0-q2**2]
    if k == 2 and m == 5:
        q3, q2, q1, q0 = sp.symbols('q3 q2 q1 q0')
        return z**5+q3*z**3+q2*z**2+q1*z+q0, [q3, q2, q1, q0], 1, []
    raise ValueError

c, Lam, al = sp.symbols('c Lambda alpha')
EQ = {
 'KdV':        (lambda U: U[3]+6*U[0]*U[1]-c*U[1], [c], 1),
 'mKdV':       (lambda U: U[3]-6*U[0]**2*U[1]-c*U[1], [c], 1),
 'Burgers':    (lambda U: U[2]-2*U[0]*U[1]+c*U[1], [c], 1),
 'Boussinesq': (lambda U: U[2]-Lam*U[0]+U[0]**2, [Lam], 1),
 'Gardner':    (lambda U: U[3]+6*U[0]*U[1]+6*al*U[0]**2*U[1]-c*U[1], [c, al], al),
}

def build(eq, k, m, dA, dB):
    Q, qv, qlead, sqgens = normQ(k, m)
    a = sp.symbols(f'a0:{dA+1}') if dA >= 0 else ()
    b = sp.symbols(f'b0:{dB+1}') if dB >= 0 else ()
    U0 = sum(a[i]*z**i for i in range(len(a))) + (sum(b[i]*z**i for i in range(len(b)))*Y if b else 0)
    Us = [U0]
    for _ in range(3): Us.append(D(Us[-1], Q, k))
    f, pars, extra = EQ[eq]
    R = red(sp.expand(f(Us)), Q, k)
    eqs = []
    if k == 1: eqs = sp.Poly(R, z).coeffs()
    else:
        for j in range(2): eqs += sp.Poly(sp.Poly(R, Y).coeff_monomial(Y**j), z).coeffs()
    lead = ([a[-1]] if dA >= 0 else []) + ([b[-1]] if b else [])
    unk = list(qv)+list(a)+list(b)+list(pars)
    return eqs, unk, sp.prod(lead)*qlead*extra, sqgens, Q, a, b, f

def consistent(eqs, unk, nz):
    s = sp.Symbol('sigma_')
    G = sp.groebner(list(eqs)+[s*nz-1], *unk, s, order='grevlex')
    return list(G.exprs) != [1]

def extract(eq, k, m, dA, dB, eqs, unk, nz, sqgens, Q, a, b, f):
    sols = sp.solve(eqs, unk, dict=True)
    for so in sols:
        if sp.simplify(nz.subs(so)) == 0: continue
        Qs = sp.expand(Q.subs(so))
        if k == 2 and sqgens and all(sp.simplify(g.subs(so)) == 0 for g in sqgens): continue
        U0 = sum(ai*z**i for i, ai in enumerate(a)) + (sum(bi*z**i for i, bi in enumerate(b))*Y if b else 0)
        U0 = U0.subs(so); Us = [U0]
        for _ in range(3): Us.append(D(Us[-1], Qs, k))
        if sp.simplify(red(sp.expand(f(Us).subs(so)), Qs, k)) == 0:
            return {str(kk): str(v) for kk, v in so.items()}
    return None

def case(args, out):
    eq, k, m, dA, dB = args
    t = time.time()
    eqs, unk, nz, sqgens, Q, a, b, f = build(eq, k, m, dA, dB)
    if not consistent(eqs, unk, nz):
        out.put(('NO_SOLUTION_CERTIFIED', 'relaxed system inconsistent', time.time()-t)); return
    if k == 2 and sqgens:
        if not any(consistent(eqs, unk, nz*g) for g in sqgens):
            out.put(('NO_SOLUTION_CERTIFIED', 'all solutions on square locus', time.time()-t)); return
    out.put(('EXISTS_PENDING_EXTRACTION', '', time.time()-t))
    sol = extract(eq, k, m, dA, dB, eqs, unk, nz, sqgens, Q, a, b, f)
    out.put(('FOUND', sol, time.time()-t) if sol else ('UNRESOLVED', 'no validated candidate extracted', time.time()-t))

def cases():
    L = []
    E5 = list(EQ)
    for e in E5:
        L += [(e, 1, 1, dA, -1) for dA in range(1, 5)]
        L += [(e, 1, 2, dA, -1) for dA in range(1, 5)]
        for m in (1, 2):
            L += [(e, 2, m, dA, dB) for dA in range(-1, 5) for dB in range(-1, 4) if not (dB == -1 and dA <= 0)]
    for e, (m, mA, mB) in [(e, x) for e in ('KdV', 'Boussinesq') for x in [(3, 4, 3), (4, 3, 2)]] + \
                          [(e, x) for e in ('mKdV', 'Gardner') for x in [(3, 5, 4), (4, 3, 2)]] + \
                          [('Burgers', (m, 4, 3)) for m in (3, 4, 5)]:
        L += [(e, 2, m, dA, dB) for dA in range(-1, mA+1) for dB in range(-1, mB+1) if not (dB == -1 and dA <= 0)]
    return L

if __name__ == '__main__':
    TL = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    fresh = len(sys.argv) > 2 and sys.argv[2] == 'fresh'
    resf = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'results612.jsonl')
    done = set()
    if fresh and os.path.exists(resf): os.remove(resf)
    if os.path.exists(resf):
        for l in open(resf): done.add(tuple(json.loads(l)['case']))
    todo = [cs for cs in cases() if tuple(cs) not in done]
    print(len(cases()), 'cases;', len(todo), 'to do', flush=True)
    running = []
    while todo or running:
        while todo and len(running) < 2:
            cs = todo.pop(0); q = mp.Queue(); p = mp.Process(target=case, args=(cs, q)); p.start()
            running.append((cs, p, q, time.time(), []))
        time.sleep(0.5)
        for r in running[:]:
            cs, p, q, t0, msgs = r
            while not q.empty(): msgs.append(q.get())
            final = [x for x in msgs if x[0] != 'EXISTS_PENDING_EXTRACTION']
            if final or not p.is_alive() or time.time()-t0 > TL:
                if not final:
                    p.terminate()
                    st = ('UNRESOLVED', 'non-square existence shown; extraction timed out' if msgs else 'timeout before decision', time.time()-t0)
                else: st = final[0]
                p.join()
                with open(resf, 'a') as fh:
                    fh.write(json.dumps({'case': list(cs), 'status': st[0], 'detail': st[1], 'seconds': round(st[2], 1)})+'\n')
                print(cs, st[0], round(st[2], 1), flush=True)
                running.remove(r)
