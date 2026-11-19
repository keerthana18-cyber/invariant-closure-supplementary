import sympy as sp
w,c,nu,al,k=sp.symbols('w c nu alpha k')
lam=nu/10; psi=-lam*sp.tanh(lam*w)
U=c+3*nu**2/25+12*nu/5*psi-12*psi**2
print('KdVB ODE:',sp.simplify((sp.diff(U,w,3)-nu*sp.diff(U,w,2)+U*sp.diff(U,w)-c*sp.diff(U,w)).rewrite(sp.exp)))
print('KdVB limits:',sp.simplify(U.subs(sp.tanh(lam*w),1)), sp.simplify(U.subs(sp.tanh(lam*w),-1)))
cc=5/sp.sqrt(6); U=sp.Rational(1,4)*(1-sp.tanh(w/sp.sqrt(24)))**2
print('Fisher AZ ODE:',sp.simplify((sp.diff(U,w,2)+cc*sp.diff(U,w)+U-U**2).rewrite(sp.exp)))
psi=-sp.tanh(w/sp.sqrt(24))/sp.sqrt(24)
print('Fisher closure form:',sp.simplify(sp.Rational(1,4)+sp.sqrt(6)*psi+6*psi**2-U))
# Gardner kink against its Riccati equation
g=sp.sqrt(-al); cg=-sp.Rational(3,2)/al-2*k**2
U=-1/(2*al)+k/g*sp.tanh(k*w); kq0=(1/g**2-cg)/(2*g)
print('Gardner kink (Riccati residual):',sp.simplify((sp.diff(U,w)-(-g*U**2+U/g+kq0)).subs({al:-3,k:sp.Rational(2,5)}).rewrite(sp.exp)))
