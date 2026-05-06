# pyright: reportUndefinedVariable=false
import numpy as np
import matplotlib.pyplot as plt

def wczytaj(*nazwy): # funkcja do wczytywania plikow od clauda
    for nazwa in nazwy:
        globals()[nazwa] = np.loadtxt(nazwa + ".csv", delimiter=',')

# wczytanie plików
wczytaj("Vmap", "misc", "Vnmap", "Gfile")
Lx=misc[0]
Ly=misc[1]
nx=int(misc[2])
ny=int(misc[3])
# x=np.linspace(0,Lx,nx)
# y=np.linspace(0,Ly,ny)
V=Vmap.reshape((nx+1,ny+1)).T
Vn=Vnmap.reshape((nx+1,ny-1)).T
E=Gfile[:,0]
G=Gfile[:,1]

# wykres potencjalu
plt.figure(figsize=(9.6,5.4))
plt.imshow(V,origin='lower',extent=[0,Lx,0,Ly],cmap='hot')
plt.xlabel(rf'$x\ [nm]$')
plt.ylabel(rf'$y\ [nm]$')
plt.colorbar()
plt.tight_layout()
plt.savefig("Vmap.png")
plt.close()

# wykres roboczego potencjalu dla qtbm
plt.figure(figsize=(6,6))
plt.imshow(Vn,origin='lower',extent=[0,Lx,0,ny-1],cmap='hot')
plt.xlabel(rf'$x\ [nm]$')
plt.ylabel(rf'$n$')
plt.colorbar()
plt.tight_layout()
plt.savefig("Vnmap.png")
plt.close()

plt.figure(figsize=(6,6))
plt.plot(E,G)
plt.tight_layout()
plt.savefig("Gfile.png")
plt.close()
