import numpy as np
import matplotlib.pyplot as plt

# wczytanie plików
Aa = np.loadtxt("Aa.csv",delimiter=",")
Amisc = np.loadtxt("Amisc.csv",delimiter=",")
Ab = np.loadtxt("Ab.csv",delimiter=',')
Ba = np.loadtxt("Ba.csv",delimiter=',')
BEL = np.loadtxt("BEL.csv",delimiter=',')
Bmisc = np.loadtxt("Bmisc.csv",delimiter=',')
Ca = np.loadtxt("Ca.csv",delimiter=',')
Cmisc = np.loadtxt("Cmisc.csv",delimiter=',')
Da = np.loadtxt("Da.csv",delimiter=',')
Dmisc = np.loadtxt("Dmisc.csv",delimiter=',')

# wykres nieskonczonej studni normalne delta tau
plt.figure(figsize=(9.6,5.4))
plt.plot(Aa[:,0],Aa[:,1],ls='--',color='black',label=rf'Analityczne; $E_1={np.pi**2/(2*Amisc[2]**2):.5f}$')
plt.plot(Aa[:,0],Aa[:,2],label=rf'Numeryczne; $E_1={Amisc[0]:.5f}$')
plt.title(rf'$N={Amisc[1]:.0f};\ \epsilon={Amisc[3]:.0e};\ |E_{{analityczne}}-E_{{numeryczne}}|={abs(np.pi**2/(2*Amisc[2]**2)-Amisc[0]):e}$')
plt.xlabel(rf'$x$')
plt.ylabel(rf'$|\psi(x)|^2$')
plt.legend()
plt.grid(ls=":")
plt.tight_layout()
plt.savefig("Aa.png")
plt.close()

# wykres nieskonczonej studni nienormalna delta tau
plt.figure(figsize=(9.6,5.4))
plt.plot(Aa[:,0],Aa[:,1],ls='--',color='black',label=rf'Analityczne; $E_1={np.pi**2/(2*Amisc[2]**2):.5f}$')
plt.plot(Aa[:,0],Ab[:],label=rf'Numeryczne; $E_1={Amisc[4]:.5f}$')
plt.title(rf'$N={Amisc[1]:.0f};\ \epsilon={Amisc[3]:.0e};\ |E_{{analityczne}}-E_{{numeryczne}}|={abs(np.pi**2/(2*Amisc[2]**2)-Amisc[4])}$')
plt.xlabel(rf'$x$')
plt.ylabel(rf'$|\psi(x)|^2$')
plt.legend()
plt.grid(ls=":")
plt.tight_layout()
plt.savefig("Ab.png")
plt.close()

# wykres skończona studnia potencjalu
# plt.figure(figsize=(9.6,5.4))
# for i in range(int(Bmisc[2])):
#     plt.plot(Aa[:,0],Ba[:,i],label=rf'$L={BEL[i,1]};\ E={BEL[i,0]}$')
# plt.title(rf'Skończona studnia; $a={Bmisc[0]}$')
# plt.legend()
# plt.xlabel(rf'$x$')
# plt.ylabel(rf'$|\psi(x)|^2$')
# plt.xticks([0, Amisc[2]/2, Amisc[2]], [r'$0$', r'$\frac{L}{2}$', r'$L$'])
# plt.grid(ls=":")
# plt.tight_layout()
# plt.savefig("Ba.png")
# plt.close()
fig,axs = plt.subplots(2,2,figsize=(12,9))
for i, ax in enumerate(axs.flat):
    L_i = BEL[i, 1]
    x_i = np.linspace(0,L_i,len(Ba[:,i]))
    xa = L_i/2-Bmisc[0]/2
    xb = L_i/2+Bmisc[0]/2
    ax.plot(x_i,Ba[:,i])
    ax.axvline(xa,color='black', ls='--')
    ax.axvline(xb, color='black', ls='--')
    ax.set_title(rf'$L={L_i:.0f};\ E={BEL[i,0]:.5f}$')
    ax.set_xlabel(rf'$x$')
    ax.set_ylabel(rf'$|\psi(x)|^2$')
    ax.grid(ls=':')
fig.suptitle(rf'$a={Bmisc[0]:.0f},\ V_0={Bmisc[1]:.0f}$')
plt.tight_layout()
plt.savefig("Ba.png")
plt.close()

# wykres 2D rozwiazania
plt.figure(figsize=(7.5,6))
plt.imshow(Ca.reshape((int(Cmisc[0]+1),int(Cmisc[0]+1))).T,origin='lower',cmap='inferno',extent=[0,Cmisc[1],0,Cmisc[1]])
plt.title(rf"$E={Cmisc[2]}$")
plt.xlabel(rf'$x$')
plt.ylabel(rf'$y$')
plt.colorbar()
plt.tight_layout()
plt.savefig('Ca.png')
plt.close()

# wykres 2D rozwiazania z potencjalem
plt.figure(figsize=(7.5,6))
plt.imshow(Da.reshape((int(Cmisc[0]+1),int(Cmisc[0]+1))).T,origin='lower',cmap='inferno',extent=[0, Dmisc[1], 0, Dmisc[1]])
plt.title(rf"$E={Dmisc[0]}$")
plt.xlabel(rf'$x$')
plt.ylabel(rf'$y$')
plt.colorbar()
plt.tight_layout()
plt.savefig('Da.png')
plt.close()
