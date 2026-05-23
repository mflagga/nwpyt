from math import pi
import time
t0 =time.time()
# pyright: reportUndefinedVariable=false
import numpy as np
import matplotlib.pyplot as plt

def wczytaj(*nazwy): # funkcja do wczytywania plikow od clauda
    for nazwa in nazwy:
        globals()[nazwa] = np.loadtxt(nazwa + ".csv", delimiter=',')

# wczytanie plików
wczytaj("mp","misc","rsV","tbaV","tbaE","przerwy")
N=int(misc[0])
Nk=int(misc[1])
a=misc[2]
V0=misc[3]
mp=mp.reshape((Nk+1,N))
k=np.linspace(-pi/a,pi/a,Nk+1)
M=4

plt.figure(figsize=(5,5))
# plt.imshow(mp,origin='lower')
for i in range(M):
    plt.plot(k,mp[:,i],label=rf'Pasmo {i+1}')
plt.legend()
plt.xlabel(r"Wektor falowy $k\ [a.u.]$")
plt.ylabel(r"Energia $E\ [a.u.]$")
plt.title(rf"$V_0={-abs(V0)};\ a={a}$")
plt.grid(ls=":")
plt.tight_layout()
plt.savefig("mp.png")
plt.close()

plt.figure(figsize=(5,5))
plt.plot(rsV[:,0],rsV[:,1])
plt.xlabel(r"$x\ [a.u.]$")
plt.ylabel(r"$V\ [a.u.]$")
plt.tight_layout()
plt.savefig("rsV.png",dpi=150)
plt.close()

plt.figure(figsize=(5,5))
plt.plot(tbaV[:,0],tbaV[:,3])
plt.xlabel(r"$x\ [a.u.]$")
plt.ylabel(r"$V_{tot}\ [a.u.]$")
plt.tight_layout()
plt.savefig("tbaV.png",dpi=150)
plt.close()

plt.figure(figsize=(5,5))
plt.plot(k,mp[:,0],label=rf"RS",c='tab:blue',lw=5,alpha=0.5)
plt.plot(k,tbaE,ls="--",label=rf"TBA",c='tab:blue')
plt.grid(ls=":")
plt.legend()
plt.xlabel(r"Wektor falowy $k\ [a.u.]$")
plt.ylabel(r"Energia $E\ [a.u.]$")
plt.tight_layout()
plt.savefig("tbaE.png",dpi=150)
plt.close()

pM=2*M
plt.figure(figsize=(5,5))
plt.plot(przerwy[:pM,0],ls="-",alpha=0.5)
plt.scatter(range(len(przerwy[:pM,0])),przerwy[:pM,0],s=8)
plt.xlabel("Numer pasma")
plt.ylabel("Szerokość pasma $[a.u.]$")
plt.title(rf"$V_0={-abs(V0)};\ a={a}$")
plt.grid(ls=":")
plt.tight_layout()
plt.savefig("W.png")
plt.close()

plt.figure(figsize=(5,5))
plt.scatter(range(len(przerwy[:pM,1])),przerwy[:pM,1],s=8,c='tab:blue')
plt.plot(przerwy[:pM,1],c='tab:blue',alpha=0.5)
plt.xlabel("Numer przerwy")
plt.ylabel("Przerwa energetyczna $[a.u.]$")
plt.title(rf"$V_0={-abs(V0)};\ a={a}$")
plt.grid(ls=":")
plt.tight_layout()
plt.savefig("deltaE.png")
plt.close()
