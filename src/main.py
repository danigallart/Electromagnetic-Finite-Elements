import mesh_circle_triangles as mt
import numpy as np
import mesh_square_triangles as ms
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

Lx = 0.05
Ly = 0.05

delx = 0.0005
dely = delx

dx = 0.02
dy = 0.01

elements, outb_nodes, ptop_nodes, pbot_nodes = ms.mesh_squareT(Lx,Ly,delx,dely)
ms.plot_mesh_squareT(elements.points,elements.simplices,outb_nodes,ptop_nodes,pbot_nodes)

#mt.mesh_plot(points,conn)