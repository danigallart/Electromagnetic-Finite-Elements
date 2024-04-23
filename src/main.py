#import mesh_circle_triangles as mt
import numpy as np
import mesh_square_triangles as ms
from scipy.spatial import Delaunay
from element_matrix_triangle import element_matrix_triangle
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
        Ae, be, area[e] = element_matrix_triangle(e, x, y, conn, eps[e], eps[e], 0, 0)
    pass



ms.plot_mesh_squareT(elements.points,elements.simplices,outb_nodes,ptop_nodes,pbot_nodes, midpoints, dielectric_elements)

#mt.mesh_plot(points,conn)