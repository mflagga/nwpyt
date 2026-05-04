# pyright: reportUndefinedVariable=false
import numpy as np
import matplotlib.pyplot as plt

def wczytaj(*nazwy): # funkcja do wczytywania plikow od clauda
    for nazwa in nazwy:
        globals()[nazwa] = np.loadtxt(nazwa + ".csv", delimiter=',')

# wczytanie danych
wczytaj("p1psi", "p1E", "p1V", "p2E", "p3E", "p3psi")
x = p1psi[:,0]
p1_psi = p1psi[:,1]
p1_E = p1E[:,0]
p1_T = p1E[:,1]
p1V0=np.max(p1V)
p2_E = p2E[:,0]
p2_T = p2E[:,1]
p3_E = p3E[:,0]
p3_T = p3E[:,1]
p3_psi = p3psi[:,0]
p3_V = p3psi[:,1]

# ==1==

# psi(x)
plt.figure(figsize=(8,8))
plt.plot(x,p1_psi,label=rf'$|\Psi|^2$')
for i in range(len(x)-1):
    plt.axvspan(x[i],x[i+1],color="black",alpha=0.5*p1V[i]/(np.max(p1V)),lw=0)
plt.legend()
plt.tight_layout()
plt.savefig("p1psi.png")
plt.close()

# T(E)
plt.figure(figsize=(8,8))
plt.axvline(p1V0, ls="--",color="black",label=rf'$V_0$')
plt.plot(p1_E,p1_T)
plt.xlabel(rf'$E\ [eV]$')
plt.ylabel(rf'T(E)')
plt.legend()
plt.grid(ls=':')
plt.tight_layout()
plt.savefig("p1E.png")
plt.close()

# ==2==

# T(E)
plt.figure(figsize=(8,8))
plt.plot(p2_E,p2_T)
plt.xlabel(rf'$E\ [eV]$')
plt.ylabel(rf'T(E)')
# plt.axvline(p1V0, ls="--",color="black",label=rf'$V_0$')
# plt.legend()
plt.grid(ls=':')
plt.tight_layout()
plt.savefig("p2E.png")
plt.close()

# ==3==

# T(E)
plt.figure(figsize=(8,8))
plt.axvline(p1V0, ls="--",color="black",label=rf'$V_0$')
plt.plot(p3_E,p3_T)
plt.xlabel(rf'$E\ [eV]$')
plt.ylabel(rf'T(E)')
plt.legend()
plt.tight_layout()
plt.grid(ls=':')
plt.savefig("p3E.png")
plt.close()

# psi(x)
plt.figure(figsize=(8,8))
plt.plot(x,p3_psi,label=rf'$|\Psi|^2$')
for i in range(len(x)-1):
    plt.axvspan(x[i],x[i+1],color="black",alpha=0.5*p3_V[i]/(np.max(p3_V)),lw=0)
plt.legend()
plt.tight_layout()
plt.savefig("p3psi.png")
plt.close()
