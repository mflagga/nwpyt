import numpy as np
import matplotlib.pyplot as plt
import os

os.makedirs("p1frames", exist_ok=True)

# wczytanie danych
p1misc = np.loadtxt('p1misc.csv',delimiter=',')
n=int(p1misc[0])
nt=int(p1misc[1])
L=p1misc[2]
N=n+1
Nt=nt+1
p1 = np.loadtxt('p1.csv',delimiter=',')
p1_psi1 = p1[0].reshape(N, Nt).T
p1_psi2 = p1[1].reshape(N, Nt).T
p1_psi3 = p1[2].reshape(N, Nt).T

x=np.linspace(0,L,N)

plt.figure(figsize=(8,8))
for p in range(Nt):
    plt.plot(x,p1_psi1[p])
    plt.plot(x,p1_psi2[p])
    plt.plot(x,p1_psi3[p])
    plt.savefig(f"./p1frames/frame_{p:04d}.png")
    plt.clf()
plt.close()
