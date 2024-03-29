#! /usr/bin/env python
import numpy as np
from mesh import make_elements, make_nodes
from plots import plot_graph
from print_HTP import print_HTP
import argparse


## Command line input to destinguish versions
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('variant', nargs='?', const=1, default=0)
    return parser.parse_args()
args = parse_arguments()

# Given data
k = 401.
h_z = 0.001
q = 1000000
c = 20.
elements_to_be_modified = [
                          40,41,42,43,44,45,46,47,
                          58,59,60,61,62,63,64,
                          76,77,78,79,80,
                          94,95,96
                          ]                             # -1 than in input file, as our numbering starts from 0

# Imported geometry
L = 0.01 # length of box
nx,ny = 10,10 # (n,n) box
variation = float(args.variant)                         # must be float as we use 4.1 and 4.2

nodes = make_nodes(L,nx,ny,variation)
elements = make_elements(nx,ny)

# Make list of elements with the entries not as nodes but as coordinates of nodes
element_coords = [[nodes[int(element[0])], nodes[int(element[1])], nodes[int(element[2])]] 
    for element in elements]
coords = element_coords

## V4
K = {}
for i in range (len(coords)):
    var = i
    val = k
    if variation == 4.1:
        if var in elements_to_be_modified:
            val = k*c
    elif variation == 4.2:
        if var in elements_to_be_modified:
            val = k/c
    K[var] = val

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
    delta = 1/2*np.abs(np.cross([x2-x1,y2-y1],[x3-x1,y3-y1]))
    # As k_xx = k_yy and k_xy = 0 the stiffness matrix is easier
    ## Stiffness matrix
    w = K[i]*h_z/(4*delta)
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

## Global stiffness matrix
def global_stiffness_matrix(dim_x,dim_y,local_matrices,All_elements):
    GSF = np.zeros((dim_x*dim_x,dim_y*dim_y))
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

##Solution and RHS vectors
T = np.zeros((nx*nx,1))
P = np.zeros((ny*ny,1))

## Boundary conditions
# Dirichlet on the upper boundary
T[nx*nx-nx:nx*nx]=293.
# Neumann on the lower boundary
P[0:ny] = q*L*h_z/(2*(nx-1))
P[1:ny-1] += q*L*h_z/(2*(nx-1))

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
T1[0:nx*nx-(nx)] = np.linalg.solve(H1,P1-np.matmul(np.transpose(H3),T2))
# Second part
P2 = np.matmul(np.transpose(H2),T1) + np.matmul(H4,T2)
## Save whole T and P
T[0:nx*nx-nx] = T1
P[nx*nx-nx:] = P2

## Save elementwise Temperatures
T_elems = {}
count = 0
for el in  elements:
    a = int(el[0])         
    b = int(el[1])
    c = int(el[2])
    T_elems[count] = [float(T[a]),float(T[b]),float(T[c])]
    count += 1

## Temperature gradient and heat flux
d_T = {}
dT = []
q_i = {}
qi = []
for el in T_elems:
    B = [coeff_b[el], coeff_c[el]]
    d_T[el] = 1/(2*Area[el])*np.matmul(B,np.transpose(T_elems[el]))
    dT.append(d_T[el])
    q_i[el] = -K[el]*d_T[el]
    qi.append(q_i[el])

## Save output
# Just to have nicer names
if variation % 1 == 0:
    VarNr = int(variation)
else:
    VarNr = variation
# Write output using given function
file = 'output_V'+str(VarNr)+'.txt'
print_HTP(H,T,P,file)
    
## Plot: Temperature Gradient and Fluxes
plot_graph(elements, nodes, T, dT, qi, nx, ny, VarNr)