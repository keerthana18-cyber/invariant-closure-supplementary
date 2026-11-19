import sympy as sp, sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep import system, p
from singular import run
mode = sys.argv[1]            # 'norm' or 'gen'; add 'fresh' as second argument to ignore earlier results
from pathlib import Path
HERE = Path(__file__).resolve().parent
c, Lam, al, s = sp.symbols('c Lam al s')
EQ = {'KdV': (lambda U: U[3]+6*U[0]*U[1]-c*U[1], [c], 1),
      'mKdV': (lambda U: U[3]-6*U[0]**2*U[1]-c*U[1], [c], 1),
      'Burgers': (lambda U: U[2]-2*U[0]*U[1]+c*U[1], [c], 1),
      'Boussinesq': (lambda U: U[2]-Lam*U[0]+U[0]**2, [Lam], 1),
      'Gardner': (lambda U: U[3]+6*U[0]*U[1]+6*al*U[0]**2*U[1]-c*U[1], [c, al], al)}

def Qnorm(k, m):   # exactly as stated in the manuscript
    q = sp.symbols('q0:6'); g2, g3 = sp.symbols('g2 g3')
    T = {(1,1): (q[1]*p, [q[1]], q[1], []),
         (1,2): (p**2+q[0], [q[0]], 1, []),
         (2,1): (p, [], 1, []),
         (2,2): (q[2]*p**2+q[0], [q[2], q[0]], q[2], [q[0]]),
         (2,3): (4*p**3-g2*p-g3, [g2, g3], 1, []),
         (2,4): (p**4+q[2]*p**2+q[1]*p+q[0], [q[2], q[1], q[0]], 1, [q[1], 4*q[0]-q[2]**2]),
         (2,5): (p**5+q[3]*p**3+q[2]*p**2+q[1]*p+q[0], list(q[:4]), 1, [])}
    return T[(k, m)]

SQ = {}
def Qgen(k, m):
    q = sp.symbols(f'q0:{m+1}')
    Q = sum(q[i]*p**i for i in range(m+1))
    gens = []
    if k == 2 and m % 2 == 0:
        if m not in SQ:   # ideal of the square locus by elimination in Singular
            import subprocess
            n = m//2; t = sp.symbols(f't0:{n+1}')
            S2 = sp.Poly(sp.expand(sum(t[i]*p**i for i in range(n+1))**2), p)
            rel = [q[i]-S2.coeff_monomial(p**i) for i in range(m+1)]
            code = f"ring r=0,({','.join(map(str,t))},{','.join(map(str,q))}),dp;ideal I={','.join(str(e).replace('**','^') for e in rel)};ideal E=eliminate(I,{'*'.join(map(str,t))});print(E);quit;"
            out = subprocess.run(['Singular','-q'],input=code,capture_output=True,text=True).stdout
            SQ[m] = [sp.sympify(l.split('=',1)[1].replace('^','**'), locals={str(v): v for v in q}) for l in out.splitlines() if '=' in l]
        gens = SQ[m]
    return Q, list(q), q[m], gens

def cases():
    L = []
    for e in EQ:
        L += [(e,1,1,dA,-1) for dA in range(1,5)] + [(e,1,2,dA,-1) for dA in range(1,5)]
        for m in (1,2): L += [(e,2,m,dA,dB) for dA in range(-1,5) for dB in range(-1,4) if not (dB==-1 and dA<=0)]
    blocks = [(e,3,4,3) for e in ('KdV','Boussinesq')]+[(e,4,3,2) for e in ('KdV','Boussinesq')]+\
             [(e,3,5,4) for e in ('mKdV','Gardner')]+[(e,4,3,2) for e in ('mKdV','Gardner')]+[('Burgers',m,4,3) for m in (3,4,5)]
    for e,m,mA,mB in blocks:
        L += [(e,2,m,dA,dB) for dA in range(-1,mA+1) for dB in range(-1,mB+1) if not (dB==-1 and dA<=0)]
    return L

def check_cert(gens, cof):
    return cof is not None and len(cof) == len(gens) and sp.expand(sum(ci*gi for ci, gi in zip(cof, gens)) - 1) == 0

fn = HERE / f'verify_{mode}.jsonl'
if len(sys.argv) > 2 and sys.argv[2] == 'fresh' and fn.exists(): fn.unlink()
done = set(tuple(json.loads(l)['case']) for l in fn.open()) if fn.exists() else set()
out = fn.open('a')
for cs in cases():
    if tuple(cs) in done: continue
    e, k, m, dA, dB = cs; t = time.time()
    Q, qv, qlead, sqg = (Qnorm if mode == 'norm' else Qgen)(k, m)
    a = sp.symbols(f'a0:{dA+1}') if dA >= 0 else ()
    b = sp.symbols(f'b0:{dB+1}') if dB >= 0 else ()
    A = sum(a[i]*p**i for i in range(len(a))); B = sum(b[i]*p**i for i in range(len(b)))
    f, pars, extra = EQ[e]
    eqs = system(f, A, B, Q, k)
    nz = sp.prod(([a[-1]] if a else []) + ([b[-1]] if b else [])) * qlead * extra
    V = list(qv)+list(a)+list(b)+pars+[s]
    lift = (mode == 'norm')
    G = eqs + [s*nz - 1]
    try:
        inc, cof, Gs = run(G, V, lift)
    except Exception as ex:
        out.write(json.dumps({'case': list(cs), 'status': 'UNDECIDED', 'detail': type(ex).__name__, 'seconds': round(time.time()-t, 2)})+'\n'); out.flush()
        print(cs, 'UNDECIDED', type(ex).__name__, flush=True); continue
    if inc:
        st = 'CERT_INCONSISTENT'
        if lift and not check_cert(Gs, cof): st = 'CERT_FAILED_CHECK'
    elif k == 2 and sqg:
        res = []
        for g in sqg:
            Gg = eqs + [s*nz*g - 1]; i2, c2, Ggs = run(Gg, V, lift)
            res.append(i2 and (not lift or check_cert(Ggs, c2)))
        st = 'CERT_SQUARE_ONLY' if all(res) else 'EXISTS_NONSQUARE'
    else:
        st = 'EXISTS'
    out.write(json.dumps({'case': list(cs), 'status': st, 'neqs': len(eqs), 'seconds': round(time.time()-t, 2)})+'\n'); out.flush()
    print(cs, st, round(time.time()-t, 1), flush=True)
