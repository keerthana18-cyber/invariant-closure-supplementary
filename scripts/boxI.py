import sympy as sp,sys
import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from search import red,D
z,Y=sp.symbols('z Y'); c,nu,s,q0=sp.symbols('c nu s q0')
Q=z**2+q0
for eqn in ['Fisher','KdV-Burgers']:
  for dA in range(1,5):
    a=sp.symbols(f'a0:{dA+1}'); U=sum(a[i]*z**i for i in range(dA+1))
    Us=[U]
    for _ in range(3): Us.append(D(Us[-1],Q,1))
    if eqn=='Fisher': N=Us[2]+c*Us[1]+Us[0]-Us[0]**2; pars=[c]; extra=1
    else: N=Us[3]-nu*Us[2]+Us[0]*Us[1]-c*Us[1]; pars=[c,nu]; extra=nu
    eqs=sp.Poly(red(N,Q,1),z).coeffs()
    G=sp.groebner(eqs+[s*extra*a[-1]-1],q0,*a,*pars,s,order='grevlex')
    if list(G.exprs)==[1]: print(eqn,dA,'NO_SOLUTION_CERTIFIED')
    else: print(eqn,dA,'FOUND',sp.solve(eqs,[q0,*a,c],dict=True))
