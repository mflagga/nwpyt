import numpy as np
import matplotlib.pyplot as plt
import os

os.makedirs('./img',exist_ok=True)

# wczytanie danych
S1a = np.loadtxt("S1a.csv",delimiter=',')

plt.figure(figsize=(8,6))
plt.imshow(S1a[:,0].reshape((65,65)).T,origin='lower',cmap='jet')
plt.colorbar()
plt.tight_layout()
plt.savefig('./img/S1mask.png')
plt.close()
