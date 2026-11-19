import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.special import ellipj
plt.rcParams.update({'font.family':'serif','font.size':10,'mathtext.fontset':'cm'})
fig,ax=plt.subplots(1,3,figsize=(10,2.9))
w=np.linspace(-30,30,800)
ax[0].plot(w,0.25*(1-np.tanh(w/np.sqrt(24)))**2,'k'); ax[0].set_title(r'(a) Fisher, $c=5/\sqrt{6}$')
nu,c=1.2,0.0; lam=nu/10; psi=-lam*np.tanh(lam*w)
ax[1].plot(w,c+3*nu**2/25+12*nu/5*psi-12*psi**2,'k'); ax[1].set_title(r'(b) KdV-Burgers, $\nu=6/5$, $c=0$')
e=np.array([-1.5,-0.2,1.2,2.5]); g=np.sqrt((e[2]-e[0])*(e[3]-e[1]))/2
m=(e[2]-e[1])*(e[3]-e[0])/((e[2]-e[0])*(e[3]-e[1])); w2=np.linspace(-6,6,800)
s=ellipj(g*w2,m)[0]**2
ax[2].plot(w2,(s*(e[2]-e[1])*e[0]-(e[2]-e[0])*e[1])/(s*(e[2]-e[1])-(e[2]-e[0])),'k')
ax[2].set_title(r'(c) Gardner, $\alpha=-1$, elliptic')
for a in ax: a.set_xlabel(r'$w$'); 
ax[0].set_ylabel(r'$U(w)$'); fig.tight_layout(); from pathlib import Path; HERE = Path(__file__).resolve().parent
fig.savefig(HERE / 'profiles.pdf'); fig.savefig(HERE / 'profiles.png', dpi=600)
