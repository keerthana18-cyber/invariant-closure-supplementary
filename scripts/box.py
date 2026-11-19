import sympy as sp,sys,multiprocessing as mpc, json
import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from search import red,D
from run import gaugeQ
z,Y=sp.symbols('z Y')
def case(args):
    eqn,m,dA,dB=args
    c,nu,s=sp.symbols('c nu s')
    Q,qv=gaugeQ(2,m)
    a=sp.symbols(f'a0:{dA+1}') if dA>=0 else ()
    b=sp.symbols(f'b0:{dB+1}') if dB>=0 else ()
    U=sum(a[i]*z**i for i in range(len(a)))+(sum(b[i]*z**i for i in range(len(b)))*Y if b else 0)
    Us=[U]
    for _ in range(3): Us.append(D(Us[-1],Q,2))
    if eqn=='Fisher': N=Us[2]+c*Us[1]+Us[0]-Us[0]**2; extra=[c]; pars=[c]
    else: N=Us[3]-nu*Us[2]+Us[0]*Us[1]-c*Us[1]; extra=[nu]; pars=[c,nu]
    R=red(N,Q,2); eqs=[]
    for j in range(2): eqs+=sp.Poly(sp.Poly(R,Y).coeff_monomial(Y**j),z).coeffs()
    lead=([a[-1]] if dA>=0 else [])+([b[-1]] if b else [])
    G=sp.groebner(eqs+[s*sp.prod(extra)*sp.prod(lead)-1],*qv,*a,*b,*pars,s,order='grevlex')
    return 'NO_SOLUTION_CERTIFIED' if list(G.exprs)==[1] else 'CONSISTENT'
if __name__=='__main__':
    cases=[(e,m,dA,dB) for e in ['Fisher','KdV-Burgers'] for m in [3,4] for dA in range(-1,4) for dB in range(-1,3) if not (dB==-1 and dA<=0)]
    pool=mpc.Pool(8); res={}
    asyncs={cs:pool.apply_async(case,(cs,)) for cs in cases}
    for cs,ar in asyncs.items():
        try: res[cs]=ar.get(timeout=600)
        except mpc.TimeoutError: res[cs]='UNRESOLVED'
        print(cs,res[cs]); sys.stdout.flush()
    json.dump({str(k):v for k,v in res.items()},open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','results','box_results.json'),'w'),indent=1)
    pool.terminate()
