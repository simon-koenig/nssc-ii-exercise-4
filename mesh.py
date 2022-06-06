#! /usr/bin/env python
import numpy as np 

# Create simple mesh of size L x L with nodes nx x ny
def make_nodes(L,nx,ny, variation = 0):
    x = np.linspace(0,L,nx)
    y = np.linspace(0,L,ny)

    # Make mesh
    xgrid,ygrid = np.meshgrid(x,y)
        
    #Print(xgrid,ygrid)
    nodes = [[xp,yp] for rows,cols in zip(xgrid,ygrid) for xp,yp in zip(rows,cols)]

    if variation == 1:
        def xp_trapez(xp,yp,L):
            return -(xp*yp)/(2*L) + xp
        nodes = [[xp_trapez(xp,yp,L),yp] for rows,cols in zip(xgrid,ygrid) for xp,yp in zip(rows,cols)]

    if variation == 2:
        def xp_quadratic(xp,yp,L):
            B = yp/(2*L)
            return xp*((xp*B/L) - B + 1)
        nodes = [[xp_quadratic(xp,yp,L),yp] for rows,cols in zip(xgrid,ygrid) for xp,yp in zip(rows,cols)]

    if variation == 3:
        def xp_annulus(xp,yp,L):
            r = 2*L - xp
            phi = np.pi - (np.pi * yp)/(4*L)
            return 2*L + r * np.cos(phi)
        def yp_annulus(xp,yp,L):
            r = 2*L - xp
            phi = np.pi - (np.pi * yp)/(4*L)
            return r * np.sin(phi)

        nodes = [[xp_annulus(xp,yp,L),yp_annulus(xp,yp,L)] for rows,cols in zip(xgrid,ygrid) for xp,yp in zip(rows,cols)]
    return nodes 


# Empty list of triangles, list will be of length (L*L*2), every list entry is a list of 
# 3 elements. With each element representing a node. 
def make_elements(nx,ny):
    elements = np.empty(((nx-1)*(ny-1)*2,3))
    for y in range(0,ny-1):
        for x in range(0,nx-1): 
            current_elem = (y*(ny-1) + x)*2
            v1 = y * ny + x
            v2 = v1 + 1
            v3 = v1 + ny
            v4 = v3 + 1
            elements[current_elem] = [v1,v2,v3]
            elements[current_elem+1] = [v2,v4,v3]
    return elements

