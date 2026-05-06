# pyright: reportUndefinedVariable=false
import numpy as np
import matplotlib.pyplot as plt
import time

t0 =time.time()

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
plt.title('Potencjał')
plt.ylabel(rf'$y\ [nm]$')
plt.colorbar()
plt.tight_layout()
plt.savefig("Vmap.png")
plt.close()

# wykres roboczego potencjalu dla qtbm
plt.figure(figsize=(6,6))
plt.imshow(Vn,origin='lower',extent=[0,1,0,1],cmap='hot')
plt.xlabel(rf'$x\ [nm]$')
plt.title('Wartości własne dla poszczególnych przekrojów')
plt.ylabel(rf'$n$')
plt.colorbar()
plt.tight_layout()
plt.savefig("Vnmap.png")
plt.close()

plt.figure(figsize=(6,6))
plt.plot(E,G)
plt.xlabel(rf'$E\ [eV]$')
plt.ylabel(rf'Znormalizowna konduktancja $\frac{{G}}{{G_0}}$')
plt.grid(ls=":")
plt.tight_layout()
plt.savefig("Gfile.png")
plt.close()

t1=time.time()
print(f'Czas wykreslania: {t1-t0:.4f}s')
