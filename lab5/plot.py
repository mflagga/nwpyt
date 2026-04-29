import numpy as np
import matplotlib.pyplot as plt

psifile = np.loadtxt("psifile.csv",delimiter=',')

plt.figure(figsize=(8,8))
plt.plot(psifile)
plt.savefig("psi.png")
