import sympy as sp,sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from search import red,D
z,Y=sp.symbols('z Y')
def gaugeQ(k,m):
    if (k,m)==(1,2): q=sp.symbols('q0'); return z**2+q,[q]
    if (k,m)==(2,3): g2,g3=sp.symbols('g2 g3'); return 4*z**3-g2*z-g3,[g2,g3]
    if (k,m)==(2,4): q2,q1,q0=sp.symbols('q2 q1 q0'); return z**4+q2*z**2+q1*z+q0,[q2,q1,q0]
def run(Nfun,params,cases,label):
    for (k,m,dA,dB) in cases:
        Q,qv=gaugeQ(k,m)
        a=sp.symbols(f'a0:{dA+1}') if dA>=0 else ()
        b=sp.symbols(f'b0:{dB+1}') if dB>=0 else ()
        U=sum(a[i]*z**i for i in range(len(a)))+(sum(b[i]*z**i for i in range(len(b)))*Y if b else 0)
        Us=[U]
        for _ in range(3): Us.append(D(Us[-1],Q,k))
        R=red(Nfun(*Us,*params),Q,k)
        eqs=[]
        for j in range(k):
            Cj=sp.Poly(R,Y).coeff_monomial(Y**j) if k==2 else R
            eqs+=sp.Poly(Cj,z).coeffs()
        unk=list(qv)+list(a)+list(b)+[p for p in params if isinstance(p,sp.Symbol) and p.name=='c']
        lead=[a[-1]] if a and dA>=1 else []
        lead+= [b[-1]] if b else []
        G=sp.groebner(eqs+[sp.Symbol('s')*sp.prod(lead)-1],*unk,sp.Symbol('s'),order='grevlex')
        if list(G.exprs)==[1]:
            print(label,(k,m,dA,dB),'NO SOLUTION (Groebner basis = {1})'); continue
        sols=sp.solve(eqs,unk,dict=True)
        sols=[s for s in sols if all(sp.simplify(x.subs(s))!=0 for x in lead)]
        print(label,(k,m,dA,dB),'solutions:',sols)
