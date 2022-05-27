#! /usr/bin/env python
from lib2to3.pytree import Node
import numpy as np
from mesh_basic import elements_coords, L, elements, nx, ny

# Given data
k = 401.
h_z = 0.001
q = 1000000
c = 20.

# Imported geometry
coords = elements_coords
# coords = [[[0.,0.],[L,0.],[0.,L]],[[L,0.],[L,L],[0.,L]]]      #this was a test

## Local stiffness matrix
def local_stiffness_matrix(element):
    # Coordinates
    x1 = element[0][0]
    x2 = element[1][0]
    x3 = element[2][0]
    y1 = element[0][1]
    y2 = element[1][1]
    y3 = element[2][1]
    # Coefficients
    a1 = x2*y3-x3*y2
    a2 = x3*y1-x1*y3
    a3 = x1*y2-x2*y1
    b1 = y2-y3
    b2 = y3-y1
    b3 = y1-y2
    c1 = x3-x2
    c2 = x1-x3
    c3 = x2-x1
    # Area
    delta = (x1*b1 + x2*b2 + x3*b3)/2
    # As k_xx = k_yy and k_xy = 0 the stiffness matrix is easier
    ## Stiffness matrix
    w = k*h_z/(4*delta)
    H_e_out = [[(b1*b1+c1*c1)*w,(b1*b2+c1*c2)*w,(b1*b3+c1*c3)*w],
               [(b2*b1+c2*c1)*w,(b2*b2+c2*c2)*w,(b2*b3+c2*c3)*w],
               [(b3*b1+c3*c1)*w,(b3*b2+c3*c2)*w,(b3*b3+c3*c3)*w]]
    return H_e_out, delta, [b1,b2,b3], [c1,c2,c3]

# Save local stiffness matrices, geometry coefficients and areas for each element
H_e = {}
Area = {}
coeff_b = {}
coeff_c = {}

for i in range (len(coords)):
    var = i
    val = local_stiffness_matrix(coords[i])
    H_e[var] = val[0]
    Area[var] = val[1]
    coeff_b[var] = val[2]
    coeff_c[var] = val[3]

#print (coeff_b)

## Global stiffness matrix
def global_stiffness_matrix(dim_x,dim_y,local_matrices,All_elements):
    GSF = np.zeros((dim_x*dim_x,dim_y*dim_y))
    #print (H)
    for el in local_matrices:
        procEL = local_matrices[el]             # now processing this element
        nodesEL = All_elements[el]              # corresponding nodes
        # to allocate the nodes in the global stiffness matrix easily
        a = int(nodesEL[0])         
        b = int(nodesEL[1])
        c = int(nodesEL[2])
        # Assign the values from local matrix to global
        GSF[a,a] += procEL[0][0]
        GSF[a,b] += procEL[0][1]
        GSF[a,c] += procEL[0][2]
        GSF[b,a] += procEL[1][0]
        GSF[b,b] += procEL[1][1]
        GSF[b,c] += procEL[1][2]
        GSF[c,a] += procEL[2][0]
        GSF[c,b] += procEL[2][1]
        GSF[c,c] += procEL[2][2]
    return GSF

# Save global stiffness matrix
H = global_stiffness_matrix(nx, ny, H_e, elements)
#print(H)

##Solution and RHS vectors
T = np.zeros((nx*nx,1))
P = np.zeros((ny*ny,1))

## Boundary conditions
# Dirichlet on the upper boundary
T[nx*nx-nx:nx*nx]=293.
# Neumann on the lower boundary
P[0:ny] = q*L*h_z/(2*(nx-1))
P[1:ny-1] += q*L*h_z/(2*(nx-1))
#print (P)

## SOLVE
#Split the system
H1 = H[0:(nx*nx-(nx)),0:(ny*ny-(ny))]
H2 = H[0:(nx*nx-(nx)),ny*ny-ny:]
H3 = H[nx*nx-nx:,0:(ny*ny-(ny))]
H4 = H[nx*nx-nx:,ny*ny-ny:]
T1 = T[0:nx*nx-nx]
T2 = T[nx*nx-nx:]
P1 = P[0:ny*ny-ny]
P2 = P[ny*ny-ny:]

# First part
T1[0:nx*nx-(nx)] = np.linalg.solve(H1,P1-np.matmul(H2,T2))      #using H2 instead of H3 because of dimensions and symmetricity

# Second part
P2 = np.matmul(H3,T1)+np.matmul(H4,T2)                          #same as above

## Save whole T and P
T[0:nx*nx-nx] = T1
P[nx*nx-nx:] = P2

T_elems = {}
count = 0
for el in  elements:
    a = int(el[0])         
    b = int(el[1])
    c = int(el[2])
    T_elems[count] = [float(T[a]),float(T[b]),float(T[c])]
    count += 1
print (T_elems)



