import subprocess, sympy as sp
class SingularError(RuntimeError): pass
def integral(e):
    """Clear denominators: same ideal, integer coefficients (avoids Singular reading g^2/4 as g^(2/4))."""
    e = sp.expand(e)
    den = sp.ilcm(1, 1, *[sp.fraction(t)[1] for t in sp.Add.make_args(e)]) if e != 0 else 1
    return sp.expand(e*den)
def to_s(e): return str(sp.expand(e)).replace('**', '^')
def run(gens, vars_, want_lift=False, timeout=1800):
    gens = [integral(g) for g in gens]
    names = [str(v) for v in vars_]
    code = f"ring r=0,({','.join(names)}),dp;\noption(redSB);\nideal I=" + ",\n".join(to_s(g) for g in gens) + ";\n"
    code += "ideal G=std(I);\nprint(\"RESULT \"+string(size(G)==1 && G[1]==1));\n"
    if want_lift:
        code += "if(size(G)==1 && G[1]==1){matrix T=lift(I,ideal(1));int j; for(j=1;j<=nrows(T);j++){print(\"COF \"+string(T[j,1]));}}\n"
    code += "quit;\n"
    r = subprocess.run(['Singular', '-q'], input=code, capture_output=True, text=True, timeout=timeout)
    out = r.stdout
    if '?' in out or r.stderr.strip() or 'RESULT' not in out:
        raise SingularError(out[:400] + r.stderr[:400])
    lines = out.splitlines()
    inconsistent = [l for l in lines if l.startswith('RESULT ')][0].split()[1] == '1'
    cof = None
    if want_lift and inconsistent:
        loc = {n: sp.Symbol(n) for n in names}
        cof = [sp.sympify(l[4:].replace('^', '**'), locals=loc) for l in lines if l.startswith('COF ')]
    return inconsistent, cof, gens
