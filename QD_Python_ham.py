# required libraries
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import math
import scipy
import itertools
import pandas as pd 
import time



site_size = 3000
site_elements = 130
PBC = False
V_max = 1 
r_0 = 1000
m_eff = 0.067
title =   "_" + str(site_elements) + "_" + str(site_size) + "_" + str(r_0) + "_" +  str(V_max) + "_" + str(PBC) + "_" + str(m_eff) + ".csv"

Hamiltonian = pd.read_csv("Hamiltonian" + title, header=None, float_precision='round_trip').to_numpy()

start = time.time()
print("start calculating")
E,psiT = np.linalg.eigh(Hamiltonian) # This computes the eigen values and eigenvectors
psi = np.transpose(psiT)   # We take the transpose of psiT to the wavefunction vectors can accessed as psi[n]
end = time.time()
print("stop calculation, time = ",end - start, "seconds")

df = pd.DataFrame(E)
df.to_csv("eigenvalues" + title, header=False, index=False)

df = pd.DataFrame(psi)
df.to_csv("psi" + title, header=False, index=False)