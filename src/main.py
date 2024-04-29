#import mesh_circle_triangles as mt
import numpy as np
import mesh_square_triangles as ms
from scipy.spatial import Delaunay
from element_matrix_triangle import element_matrix_triangle
from postprocess import  postprocessU
import matplotlib.pyplot as plt

element_type = 1 #1 for triangular

Lx = 0.05 #X length of the computational domain
Ly = 0.05 #Y length of the computational domain

delx = 0.0005 #Delta x
dely = delx   #Delta y

dx = 0.02 #X length of the dielectric and plates
dy = 0.01 #Y length of the dielectric

elements, outb_nodes, ptop_nodes, pbot_nodes, midpoints, dielectric_elements = ms.mesh_squareT(Lx,Ly,delx,dely)


M = len(elements.simplices) #Number of elements
N = len(elements.points)    #Number of nodes

x = elements.points[:,0] #x coordinate of points
y = elements.points[:,1] #y coordinate of points
conn = elements.simplices #connectivity matrix with element per row and node per column
#Input parameters
e0 = 8.85e-12 #Free space permittivity (F/m)
er = 1        #Dielectric constant inside the capacitor

#Boundary Conditions
Vtop = 1     #Dirichlet BC at the top plate
Vbottom = 0  #Dirichlet BC at the bottom plate

#Material and source parameters
eps = np.ones((M,1))*e0
eps[dielectric_elements] = er*e0

#Main Body
if element_type == 1:
    Ne = 3 #Number of nodes in each triangular element

A = np.zeros((N,N)) #Global matrix
b = np.zeros((N,1)) #Global b vector
area = np.zeros((M,1)) #Area of elements

for e in range(M):
    if element_type == 1:
        Ae, be, area[e], Jdet, x1, x2, x3, y1, y2, y3 = element_matrix_triangle(e, x, y, conn, eps[e][0], eps[e][0], 0, 0)
    pass

    for i in range(Ne):
        ig = conn[e,i]
        for j in range(Ne):
            jg = conn[e,j]
            A[ig,jg] += Ae[i,j]
        b[ig] += be[i]

#Impose Dirichlet Bounday Conditions

BCnodes = np.concatenate((outb_nodes,ptop_nodes,pbot_nodes))
BCvalues = np.concatenate((np.zeros(outb_nodes.size),Vtop*np.ones(ptop_nodes.size),Vbottom*np.ones(pbot_nodes.size)))

nodes = np.arange(N)
nodes[BCnodes] = 0
indexes = np.where(nodes != 0)[0]

A2 = np.diag(np.diag(np.ones((A.shape[0],A.shape[0]))))
A2[indexes,:] = 0.
A[BCnodes, :] = 0.
A = A + A2

b[BCnodes] = BCvalues.reshape(-1,1)

V = np.linalg.solve(A,b)

#Post-processing: map V to a 2D grid

x2d, y2d, V2d = postprocessU(Lx,Ly,delx,dely,N, elements, V)

cmap = plt.colormaps['plasma']
plt.pcolormesh(x2d,y2d,V2d,cmap=cmap)
plt.colorbar()
plt.show()


ms.plot_mesh_squareT(elements.points,elements.simplices,outb_nodes,ptop_nodes,pbot_nodes, midpoints, dielectric_elements)

#mt.mesh_plot(points,conn)