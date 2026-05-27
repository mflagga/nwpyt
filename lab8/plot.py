# pyright: reportUndefinedVariable=false
import numpy as np
import matplotlib.pyplot as plt
from math import pi

def wczytaj(*nazwy): # funkcja do wczytywania plikow od clauda
    for nazwa in nazwy:
        globals()[nazwa] = np.loadtxt(nazwa + ".csv", delimiter=',')

wczytaj("jeden","H","Hk","misc","Fre","Fim")
N=int(misc[0])
a=misc[1]
eps=misc[2]
t=misc[3]

k = np.linspace(0.0,2.0*pi/a,N)
dysp=np.empty(N)
for i in range(N):
    dysp[i]=Hk[i*N+i]
teor = eps + 2.0*t*np.cos(k*a)

size = 4.5
Modcmap='binary'
Recmap='bwr'

plt.figure(figsize=(1.2*size,size))
plt.imshow(jeden.reshape((N,N)),cmap=Modcmap)
plt.title(rf"$|\hat{{F}}^\dagger \hat{{F}}|$")
plt.colorbar()
plt.tight_layout()
plt.savefig("jeden.png",dpi=150)
plt.close()

Hmax=np.abs(H).max()
plt.figure(figsize=(1.2*size,size))
plt.imshow(H.reshape((N,N)),cmap=Recmap,vmin=-Hmax,vmax=Hmax)
plt.title(rf"$\Re(\hat{{H}})$")
plt.colorbar()
plt.tight_layout()
plt.savefig("H.png",dpi=150)
plt.close()

Hkmax=np.abs(Hk).max()
plt.figure(figsize=(1.2*size,size))
plt.imshow(Hk.reshape((N,N)),cmap=Recmap,vmin=-Hkmax,vmax=Hkmax)
plt.title(rf"$\Re(\hat{{H}}_k)$")
plt.colorbar()
plt.tight_layout()
plt.savefig("Hk.png",dpi=150)
plt.close()

W=dysp.max()-dysp.min()
Esr=np.average(dysp)
plt.figure(figsize=(size,size))
plt.plot(k,dysp,ls='--',c='k',label=rf"Wyznaczona")
plt.plot(k,teor,c='tab:blue',alpha=0.75,label=rf"Analityczna")
plt.xlabel(rf"$k\ [a.u.]$")
plt.ylabel(rf"$E\ [a.u.]$")
plt.title(rf"$\epsilon={eps:.2f};\ t={t}$"+"\n"+rf"$\left< E \right>={Esr:.2f};\ W={W}$")
plt.grid(ls=":")
plt.legend()
plt.tight_layout()
plt.savefig("dysp.png",dpi=150)
plt.close()

Fremax=np.abs(Fre).max()
plt.figure(figsize=(1.2*size,size))
plt.imshow(Fre.reshape((N,N)),cmap=Recmap,vmin=-Fremax,vmax=Fremax)
plt.title(rf"$\Re(\hat{{F}})$")
plt.colorbar()
plt.tight_layout()
plt.savefig("Fre.png",dpi=150)
plt.close()

Fimmax=np.abs(Fim).max()
plt.figure(figsize=(1.2*size,size))
plt.imshow(Fim.reshape((N,N)),cmap=Recmap,vmin=-Fimmax,vmax=Fimmax)
plt.title(rf"$\Im(\hat{{F}})$")
plt.colorbar()
plt.tight_layout()
plt.savefig("Fim.png",dpi=150)
plt.close()
