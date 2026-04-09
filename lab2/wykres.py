import numpy as np
import matplotlib.pyplot as plt

Aa = np.loadtxt("Aa.csv",delimiter=",")
Amisc = np.loadtxt("Amisc.csv",delimiter=",")

plt.figure(figsize=(9.6,5.4))
plt.plot(Aa[:,0],Aa[:,1],ls='--',color='black',label=rf'Analityczne; $E_1={np.pi**2/(2*Amisc[2]**2):.5f}$')
plt.plot(Aa[:,0],Aa[:,2],label=rf'Numeryczne; $E_1={Amisc[0]:.5f}$')
plt.title(rf'$|E_{{analityczne}}-E_{{numeryczne}}|={abs(np.pi**2/(2*Amisc[2]**2)-Amisc[0]):e}$')
plt.xlabel(rf'$x$')
plt.ylabel(rf'$|\psi(x)|^2$')
plt.legend()
plt.grid(ls=":")
plt.tight_layout()
plt.savefig("Aa.png")
plt.close()