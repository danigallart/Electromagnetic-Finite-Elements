import numpy as np
from scipy.spatial import Delaunay
#Generate the mesh
import matplotlib.pyplot as plt


def mesh_squareT(Lx=1.0,Ly=1.0,delx=0.1,dely=0.1,dx=0.02,dy=0.01):
        
    Nx = round(Lx/delx) + 1
    Ny = round(Ly/dely) + 1

    xt = np.linspace(-Lx/2, Lx/2, Nx)
    yt = np.linspace(-Ly/2, Ly/2, Ny)

    points = []
    for i in xt:
        for j in yt:
            points.append([i,j])
    points = np.array(points)
    elements = Delaunay(points)
    
    #Boundary Nodes
    xLeftb = -Lx/2.
    xrightb = Lx/2.
    ybottomb = -Ly/2.
    ytopb = Ly/2.

    outb_nodes = np.where((elements.points[:,0]==xLeftb) | (elements.points[:,0]==xrightb) | (elements.points[:,1]==ybottomb) | (elements.points[:,1]==ytopb))

    #Nodes of the plates
    ptop_nodes = np.where((elements.points[:,1]<dy/2.+1e-8) & (elements.points[:,1]>dy/2.-1e-8) & (elements.points[:,0]<dx/2.+1e-8) & (elements.points[:,0]>-dx/2.-1e-8))
    pbot_nodes = np.where((elements.points[:,1]>-dy/2.-1e-8) & (elements.points[:,1]<-dy/2.+1e-8) & (elements.points[:,0]<dx/2.+1e-8) & (elements.points[:,0]>-dx/2.-1e-8))
    return elements, outb_nodes, ptop_nodes, pbot_nodes


def plot_mesh_squareT(points,simplices,outb_nodes, ptop_nodes, pbot_nodes):
    plt.triplot(points[:,0], points[:,1], simplices)
    plt.plot(points[:,0], points[:,1], 'ko',markersize=0.5)
    plt.scatter(points[outb_nodes][:,0],points[outb_nodes][:,1],s=2., c = 'orange')
    plt.plot(points[ptop_nodes][:,0],points[ptop_nodes][:,1],'r',markersize=1.)
    plt.plot(points[pbot_nodes][:,0],points[pbot_nodes][:,1],'r',markersize=1.)
    plt.show()