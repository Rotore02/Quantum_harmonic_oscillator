import numpy as np
import matplotlib.pyplot as mp
from scipy.linalg import eigh

mass = 1.0
omega = 1.0
hbar = 1.0

x_min = -10.0
x_max = 10.0
N_step = 2000
dx = (x_max - x_min)/N_step

x = np.linspace(-10.0, 10.0, N_step)
X = np.diag(x)

P = (-1j/(2*dx))*(np.diag(np.ones(N_step-1),1)-np.diag(np.ones(N_step-1),-1))

T = ((P @ P)/(2*mass)) 

V = 0.5 * mass * omega**2 * (X @ X)

H = T + V

E, psi = eigh(H)

file = open("Eigenvalues.txt", 'w')

for i in range(5): 
    file.write(str(E[i])+'\n')
    mp.plot(x, np.real(psi[:,i]))
file.close()
mp.title("Eigenfunctions")
mp.xlabel("Positions*{\displaystyle \hbar } (m)")
mp.savefig("Eigenfunctions.pdf")
    
