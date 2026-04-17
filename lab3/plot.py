import numpy as np
import matplotlib.pyplot as plt
import os

# wczytanie danych
S1a = np.loadtxt("S1a.csv",delimiter=',')
S1misc = np.loadtxt("S1misc.csv",delimiter=',')

# mapa układu
plt.figure(figsize=(8.5,6))
plt.imshow(S1a[:,1].reshape((int(S1misc[0])+1,int(S1misc[1])+1)).T,origin='lower',cmap='Dark2',vmin=-0.2,vmax=4.3)
plt.colorbar(ticks=[0,1,2,3,4]).ax.set_yticklabels(['Srodek',f'Dirichlet1={S1misc[2]:.2f}',f'Dirichlet2={S1misc[3]:.2f}','Neumann','Born-Karmann'])
plt.tight_layout()
plt.savefig('S1mask.png')
plt.close()

# mapa ładunku
plt.figure(figsize=(8,6))
plt.imshow(S1a[:,2].reshape((int(S1misc[0])+1,int(S1misc[1])+1)).T,origin='lower',cmap='bwr')
plt.colorbar()
plt.tight_layout()
plt.savefig('S1rho.png')
plt.close()
