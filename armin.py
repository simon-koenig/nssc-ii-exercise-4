import numpy as np
import matplotlib.pyplot as plt


class Mesh(self):
    pass



    
_N = 100
H = np.zeros((_N,_N))
T = np.zeros((N))
P = np.zeros((N))

# splitting matrix and vectors

split_index = 90

T[90:] = ...
P[:90] = ...


LHS = H[:90][:90] # H 1...90 1...90 
RHS = P[:90] - np.dot(H[90:,:90],T[90:])

T[:90] = np.linalg.solve(LHS, RHS) 
P[90:] = np.dot(H[90:][:], T[:])