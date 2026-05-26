# pyright: reportUndefinedVariable=false
import numpy as np
import matplotlib.pyplot as plt

def wczytaj(*nazwy): # funkcja do wczytywania plikow od clauda
    for nazwa in nazwy:
        globals()[nazwa] = np.loadtxt(nazwa + ".csv", delimiter=',')
        