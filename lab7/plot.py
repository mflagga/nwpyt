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
wczytaj("mp","misc","rsV")
N=int(misc[0])
Nk=int(misc[1])
a=misc[2]
mp=mp.reshape((Nk+1,N))
k=np.linspace(-pi/a,pi/a,Nk+1)
M=4

plt.figure(figsize=(8,8))
# plt.imshow(mp,origin='lower')
for i in range(M):
    plt.plot(k,mp[:,i],label=rf'Pasmo {i+1}')
plt.legend()
plt.xlabel(r"Wektor falowy $k$")
plt.ylabel(r"Energia $E$")
plt.grid(ls=":")
plt.tight_layout()
plt.savefig("mp.png")

plt.figure(figsize=(8,8))
plt.plot(rsV[:,0],rsV[:,1])
plt.tight_layout()
plt.savefig("rsV.png",dpi=150)
