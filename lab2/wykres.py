import numpy as np
import matplotlib.pyplot as plt

Aa = np.loadtxt("Aa.csv",delimiter=",")
Amisc = np.loadtxt("Amisc.csv",delimiter=",")

plt.figure(figsize=(12.8,7.2))
plt.plot(Aa[:,0],Aa[:,1],ls='--',color='black',label=rf'Analityczne; $E_1={np.pi**2/(2*Amisc[2]**2):.5f}$')
plt.plot(Aa[:,0],Aa[:,2],label=rf'Numeryczne; $E_1={Amisc[0]:.5f}$')
plt.xlabel(rf'$x$')
plt.ylabel(rf'$|\psi(x)|^2$')
plt.legend()
plt.grid(ls=":")
plt.tight_layout()
plt.savefig("Aa.png")
plt.close()