# pyright: reportUndefinedVariable=false
import numpy as np
import matplotlib.pyplot as plt

def wczytaj(*nazwy): # funkcja do wczytywania plikow od clauda
    for nazwa in nazwy:
        globals()[nazwa] = np.loadtxt(nazwa + ".csv", delimiter=',')

# wczytanie plików
wczytaj("Vmap", "misc")
Lx=misc[0]
Ly=misc[1]
nx=int(misc[2])
ny=int(misc[3])
x=np.linspace(0,Lx,nx)
y=np.linspace(0,Ly,ny)
V=Vmap.reshape((nx+1,ny+1)).T

plt.figure(figsize=(9.6,5.4))
plt.imshow(V,origin='lower',extent=[0,Lx,0,Ly],cmap='hot')
plt.xlabel(rf'$x\ [nm]$')
plt.ylabel(rf'$y\ [nm]$')
# plt.colorbar()
plt.tight_layout()
plt.savefig("Vmap.png")
plt.close()
