import numpy as np
import matplotlib.pyplot as plt
import os
import time

os.makedirs("p1frames", exist_ok=True)
os.makedirs("p2frames", exist_ok=True)
os.makedirs("p3frames", exist_ok=True)

# wczytanie danych
p1misc = np.loadtxt('p1misc.csv',delimiter=',')
n=int(p1misc[0])
nt=int(p1misc[1])
L=p1misc[2]
N=n+1
Nt=nt+1
dt=p1misc[6]

p1 = np.loadtxt('p1.csv',delimiter=',')
p1_psi1 = p1[0].reshape(N, Nt).T
p1_psi2 = p1[1].reshape(N, Nt).T
p1_psi3 = p1[2].reshape(N, Nt).T

p2=np.loadtxt('p2.csv',delimiter=',')
p2_psi = p2.reshape(N,Nt).T

p2gamma=np.loadtxt("p2gamma.csv",delimiter=',')

p3 = np.loadtxt('p3.csv',delimiter=',')
p3_psi1 = p3[0].reshape(N, Nt).T
p3_psi2 = p3[1].reshape(N, Nt).T
p3_psi3 = p3[2].reshape(N, Nt).T

p3gamma=np.loadtxt("p3gamma.csv",delimiter=',')

p3misc=np.loadtxt("p3misc.csv",delimiter=',')

x=np.linspace(0,L,N)
psimax = max(p1_psi1.max(), p1_psi2.max(), p1_psi3.max())
psimax3 = max(p3_psi1.max(), p3_psi2.max(), p3_psi3.max())

# animacja pierwsza
t0=time.time()
plt.figure(figsize=(8,8))
for p in range(Nt):
    plt.plot(x,p1_psi1[p],label=rf"$\sigma_0={p1misc[3]*0.5}$")
    plt.plot(x,p1_psi2[p],label=rf"$\sigma_0={p1misc[4]*0.5}$")
    plt.plot(x,p1_psi3[p],label=rf"$\sigma_0={p1misc[5]*0.5}$")
    plt.legend()
    plt.title(rf"$t={p*dt:.3f}$")
    plt.ylim(0,psimax)
    plt.xlim(0,L)
    plt.savefig(f"./p1frames/frame_{p:04d}.png")
    plt.clf()
plt.close()
t1=time.time()
print(f"Animacja pierwsza: {Nt} klatek w czasie {t1-t0:.2f}s")

#animacja durga
plt.figure(figsize=(8,8))
for p in range(Nt):
    plt.plot(x,p2_psi[p],label=rf"$|\Psi|^2$")
    plt.plot(x,p2gamma[:,0],label=rf"$\Re(\Gamma)$")
    plt.plot(x,p2gamma[:,1],label=rf"$-\Im(\Gamma)$")
    plt.title(rf"$t={p*dt:.3f}$")
    plt.ylim(0,p2_psi.max())
    plt.xlim(0,L)
    plt.legend()
    plt.savefig(f"./p2frames/frame_{p:04d}.png")
    plt.clf()
plt.close()
t2=time.time()
print(f"Animacja druga: {Nt} klatek w czasie {t2-t1:.2f}s")

# animacja trzecia
plt.figure(figsize=(8,8))
for p in range(Nt):
    plt.plot(x,p3_psi1[p],label=rf"$E={p3misc[1]}$")
    plt.plot(x,p3_psi2[p],label=rf"$E={p3misc[2]}$")
    plt.plot(x,p3_psi3[p],label=rf"$E={p3misc[3]}$")
    plt.plot(x,p3gamma[:,0],label=rf"$\Re(\Gamma)$")
    plt.plot(x,p3gamma[:,1],label=rf"$\Im(\Gamma)$")
    plt.legend()
    plt.title(rf"$V_0={p3misc[0]}$; $t={p*dt:.3f}$")
    plt.ylim(0,psimax3)
    plt.xlim(0,L)
    plt.savefig(f"./p3frames/frame_{p:04d}.png")
    plt.clf()
plt.close()
t3=time.time()
print(f"Animacja trzecia: {Nt} klatek w czasie {t3-t2:.2f}s")
