import numpy as np
import matplotlib.pyplot as plt
import os

syscmap = 'Accent'
rhocmap = 'bwr'
vcmap = 'hot'

# wczytanie danych
S1a = np.loadtxt("S1a.csv",delimiter=',')
S1misc = np.loadtxt("S1misc.csv",delimiter=',')
S1rL = np.loadtxt("S1rL.csv",delimiter=',').T
S3a = np.loadtxt("S3a.csv",delimiter=',')
S3misc = np.loadtxt("S3misc.csv",delimiter=',')
S3r = np.loadtxt("S3r.csv",delimiter=',').T
S2a = np.loadtxt("S2a.csv",delimiter=',')
S2misc = np.loadtxt("S2misc.csv",delimiter=',')
S2rV = np.loadtxt("S2rV.csv", delimiter=',')

# Uklad pierwszy

# mapa układu
plt.figure(figsize=(7.5,6))
plt.imshow(S1a[:,1].reshape((int(S1misc[0])+1,int(S1misc[1])+1)).T,origin='lower',cmap=syscmap,vmin=-0.2,vmax=4.3,extent=[-S1misc[4]/2,S1misc[4]/2,-S1misc[5]/2,S1misc[5]/2])
plt.colorbar(ticks=[0,1,2,3,4]).ax.set_yticklabels(['Srodek',f'Dirichlet1={S1misc[2]:.2f}',f'Dirichlet2={S1misc[3]:.2f}','Neumann','Born-Karmann'])
# plt.tight_layout()
plt.xlabel(rf'$x$')
plt.ylabel(rf'$y$')
plt.savefig('S1mask.png')
plt.close()

# mapa ładunku
plt.figure(figsize=(7.5,6))
plt.imshow(S1a[:,2].reshape((int(S1misc[0])+1,int(S1misc[1])+1)).T,origin='lower',cmap=rhocmap,extent=[-S1misc[4]/2,S1misc[4]/2,-S1misc[5]/2,S1misc[5]/2])
plt.title(rf'$\rho(x,y)$')
plt.xlabel(rf'$x$')
plt.ylabel(rf'$y$')
plt.colorbar()
# plt.tight_layout()
plt.savefig('S1rho.png')
plt.close()

# mapa potencjału
plt.figure(figsize=(7.5,6))
plt.imshow(S1a[:,0].reshape((int(S1misc[0])+1,int(S1misc[1])+1)).T,origin='lower',cmap=vcmap,extent=[-S1misc[4]/2,S1misc[4]/2,-S1misc[5]/2,S1misc[5]/2])
plt.title(rf'$V(x,y)$')
plt.xlabel(rf'$x$')
plt.ylabel(rf'$y$')
plt.colorbar()
# plt.tight_layout()
plt.savefig('S1V.png')
plt.close()

# mapa potencjalu dla roznych wartosci L
Lvals = [5.0, 10.0, 15.0, 20.0]
nx = int(S1misc[0])
ny = int(S1misc[1])
L0 = S1misc[4]
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
# V0 = S1a[:,0].reshape((nx+1, ny+1)).T
# im = axes[0,0].imshow(V0, origin='lower', cmap=vcmap,
#                      extent=[-L0/2, L0/2, -L0/2, L0/2])
# axes[0,0].set_title(rf'$V(x,y)$, $L={L0:.1f}$')
# axes[0,0].set_xlabel(r'$x$')
# axes[0,0].set_ylabel(r'$y$')
# fig.colorbar(im, ax=axes[0,0])
positions = [(0,0), (0,1), (1,0), (1,1)]
for k, (L, pos) in enumerate(zip(Lvals, positions)):
    VL = S1rL[:,k].reshape((nx+1, ny+1)).T
    ax = axes[pos]
    im = ax.imshow(VL, origin='lower', cmap=vcmap,
                   extent=[-L/2, L/2, -L/2, L/2])
    ax.set_title(rf'$V(x,y)$, $L={L:.1f}$')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    fig.colorbar(im, ax=ax)
# plt.tight_layout()
plt.savefig('S1VrL.png')
plt.close()

# wykresy ukladu 3
# mapa układu
plt.figure(figsize=(7.5,6))
plt.imshow(S3a[:,1].reshape((int(S3misc[0])+1,int(S3misc[1])+1)).T,origin='lower',cmap=syscmap,vmin=-0.2,vmax=4.3,extent=[-S3misc[4]/2,S3misc[4]/2,-S3misc[5]/2,S3misc[5]/2])
plt.colorbar(ticks=[0,1,2,3,4]).ax.set_yticklabels(['Srodek',f'Dirichlet1={S3misc[2]:.2f}',f'Dirichlet2={S3misc[3]:.2f}','Neumann','Born-Karmann'])
# plt.tight_layout()
plt.xlabel(rf'$x$')
plt.ylabel(rf'$y$')
plt.savefig('S3mask.png')
plt.close()

# mapa ładunku
plt.figure(figsize=(7.5,6))
plt.imshow(S3a[:,2].reshape((int(S3misc[0])+1,int(S3misc[1])+1)).T,origin='lower',cmap=rhocmap,extent=[-S3misc[4]/2,S3misc[4]/2,-S3misc[5]/2,S3misc[5]/2])
plt.title(rf'$\rho(x,y)$')
plt.xlabel(rf'$x$')
plt.ylabel(rf'$y$')
plt.colorbar()
# plt.tight_layout()
plt.savefig('S3rho.png')
plt.close()

# mapa potencjału
plt.figure(figsize=(7.5,6))
plt.imshow(S3a[:,0].reshape((int(S3misc[0])+1,int(S3misc[1])+1)).T,origin='lower',cmap=vcmap,extent=[-S3misc[4]/2,S3misc[4]/2,-S3misc[5]/2,S3misc[5]/2])
plt.title(rf'$V(x,y)$')
plt.xlabel(rf'$x$')
plt.ylabel(rf'$y$')
plt.colorbar()
# plt.tight_layout()
plt.savefig('S3V.png')
plt.close()

# rozne rho i V
# mapa potencjalu dla roznych wartosci r i VB1 - uklad trzeci
params = [(7, -10), (14, -10), (7, -20), (14, -20)]
nx3 = int(S3misc[0])
ny3 = int(S3misc[1])
L3x = S3misc[4]
L3y = S3misc[5]

# dane dla kazdego panelu
V_data = []
for k in range(4):
    if k == 0:
        V_data.append(S3a[:,0].reshape((nx3+1, ny3+1)).T)
    else:
        V_data.append(S3r[:,k-1].reshape((nx3+1, ny3+1)).T)

# wspolne zakresy dla tego samego VB1
vmin_top = min(V_data[0].min(), V_data[1].min())
vmax_top = max(V_data[0].max(), V_data[1].max())
vmin_bot = min(V_data[2].min(), V_data[3].min())
vmax_bot = max(V_data[2].max(), V_data[3].max())
vranges = [(vmin_top, vmax_top), (vmin_top, vmax_top),
           (vmin_bot, vmax_bot), (vmin_bot, vmax_bot)]

fig, axes = plt.subplots(2, 2, figsize=(13, 10))
positions = [(0,0), (0,1), (1,0), (1,1)]
for k, ((r, VB1), pos) in enumerate(zip(params, positions)):
    vmin, vmax = vranges[k]
    ax = axes[pos]
    im = ax.imshow(V_data[k], origin='lower', cmap=vcmap,
                   extent=[-L3x/2, L3x/2, -L3y/2, L3y/2],
                   vmin=vmin, vmax=vmax)
    ax.set_title(rf'$V(x,y)$, $r={r}$, $V_{{B1}}={VB1}$')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    fig.colorbar(im, ax=ax)
# plt.tight_layout()
plt.savefig('S3Vr.png')
plt.close()

# wykresy dla ukladu drugiego
# mapa układu
plt.figure(figsize=(10,4))
plt.imshow(S2a[:,1].reshape((int(S2misc[0])+1,int(S2misc[1])+1)).T,origin='lower',cmap=syscmap,vmin=-0.2,vmax=4.3,extent=[-S2misc[4]/2,S2misc[4]/2,-S2misc[5]/2,S2misc[5]/2])
plt.colorbar(ticks=[0,1,2,3,4]).ax.set_yticklabels(['Srodek',f'Dirichlet1={S2misc[2]:.2f}',f'Dirichlet2={S2misc[3]:.2f}','Neumann','Born-Karmann'])
# plt.tight_layout()
plt.xlabel(rf'$x$')
plt.ylabel(rf'$y$')
plt.savefig('S2mask.png')
plt.close()

# mapa ładunku
plt.figure(figsize=(10,4))
plt.imshow(S2a[:,2].reshape((int(S2misc[0])+1,int(S2misc[1])+1)).T,origin='lower',cmap=rhocmap,extent=[-S2misc[4]/2,S2misc[4]/2,-S2misc[5]/2,S2misc[5]/2])
plt.title(rf'$\rho(x,y)$')
plt.xlabel(rf'$x$')
plt.ylabel(rf'$y$')
plt.colorbar()
# plt.tight_layout()
plt.savefig('S2rho.png')
plt.close()

# mapa potencjału
plt.figure(figsize=(10,4))
plt.imshow(S2a[:,0].reshape((int(S2misc[0])+1,int(S2misc[1])+1)).T,origin='lower',cmap=vcmap,extent=[-S2misc[4]/2,S2misc[4]/2,-S2misc[5]/2,S2misc[5]/2])
plt.title(rf'$V(x,y)$')
plt.xlabel(rf'$x$')
plt.ylabel(rf'$y$')
plt.colorbar()
# plt.tight_layout()
plt.savefig('S2V.png')
plt.close()

# przekroj dla roznych V1
VB1vals = [-10.0, -15.0, -20.0, -25.0]
ny2 = int(S2misc[1])
L2y = S2misc[5]
y = np.linspace(-L2y/2, L2y/2, ny2+1)

plt.figure(figsize=(8,6))
for k, VB1 in enumerate(VB1vals):
    plt.plot(y, S2rV[k], label=rf'$V_{{B1}}={VB1}$')
plt.xlabel(r'$y$')
plt.ylabel(r'$V(x=0, y)$')
plt.title(r'Przekroje potencjalu wzdluz $y$ przy $x=0$')
plt.legend()
plt.grid(alpha=0.3)
# plt.tight_layout()
plt.savefig('S2rV.png')
plt.close()
