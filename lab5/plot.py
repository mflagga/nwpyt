# pyright: reportUndefinedVariable=false
import numpy as np
import matplotlib.pyplot as plt

def wczytaj(*nazwy): # funkcja do wczytywania plikow od clauda
    for nazwa in nazwy:
        globals()[nazwa] = np.loadtxt(nazwa + ".csv", delimiter=',')

# wczytanie danych
wczytaj("p1psi", "p1E")
x = p1psi[:,0]
p1_psi = p1psi[:,1]
p1_E = p1E[:,0]
p1_T = p1E[:,1]

plt.figure(figsize=(8,8))
plt.plot(x,p1_psi)
plt.savefig("p1psi.png")
plt.close()

plt.figure(figsize=(8,8))
plt.plot(p1_E,p1_T)
plt.savefig("p1E.png")
plt.close()
