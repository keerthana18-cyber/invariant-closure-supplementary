import sympy as sp, itertools, sys
z,Y=sp.symbols('z Y')
def red(e,Q,k):
    e=sp.expand(e)
    if k==1: return sp.expand(e.subs(Y,Q))
    p=sp.Poly(e,Y); out=0
    for (d,),co in p.terms(): out+=co*Q**(d//2)*Y**(d%2)
    return sp.expand(out)
def D(e,Q,k):
    if k==1: return sp.expand(sp.diff(e,z)*Q)
    return red(sp.diff(e,z)*Y+sp.diff(e,Y)*sp.diff(Q,z)/2,Q,k)
def system(Nfun,k,m,dA,dB,params):
    q=sp.symbols(f'q0:{m+1}'); a=sp.symbols(f'a0:{dA+1}') if dA>=0 else ()
    b=sp.symbols(f'b0:{dB+1}') if dB>=0 else ()
    Q=sum(q[i]*z**i for i in range(m+1))
    U=sum(a[i]*z**i for i in range(len(a)))+ (sum(b[i]*z**i for i in range(len(b)))*Y if b else 0)
    Us=[U]
    for _ in range(4): Us.append(D(Us[-1],Q,k))
    R=red(Nfun(*Us,*params),Q,k)
    eqs=[]
    for j in range(k):
        Cj=sp.Poly(R,Y).coeff_monomial(Y**j) if k==2 else R
        eqs+= [e for e in sp.Poly(Cj,z).coeffs()]
    return eqs,q,a,b,Q,U
