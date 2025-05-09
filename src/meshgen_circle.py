import numpy as np
from scipy.spatial import Delaunay

def meshgen_circle(obj, freq):
    c0 = 3*1e8        # m/sec, velocity of light in free space
    lmda = c0/freq  # meter, wavelength

    delh = obj.delh             # element size
    scat_type = obj.scat_type   # 'pec' or 'diel' (PEC or dielectric object)
    rc = obj.rc                 # radius of circular scatterer

    fsdim = 0.5*lmda      # distance btw scatterer and inner PML boundary
    pmldim = 0.5*lmda     # distance btw inner and outer PML boundaries
    
    rpmlin = rc+fsdim       # radius of inner PML boundary
    rpmlout = rpmlin+pmldim # radius of outer PML boundary

    #fit to delh incrementation

    rc = np.round(rc/delh)*delh
    rpmlout = np.round(rpmlout/delh)*delh
    rpmlin = np.round(rpmlin/delh)*delh

    # Create coordinates layer by layer
    x = np.array([0])
    y = np.array([0])

    rads = np.arange(delh,rpmlout+delh,delh)

    for r in rads:
        N = int(np.floor(np.pi/(2*np.arcsin(delh*0.5/r))))
        #print('N:', N)
        alfa_1 = np.linspace(0,np.pi,N)
        alfa_2 = np.linspace(-np.pi+np.pi/N,-np.pi/N,N-2)
        alfa = np.concatenate((alfa_1, alfa_2))
        x0 = np.array(r*np.cos(alfa)) # x-ccordinate of the circle
        y0 = np.array(r*np.sin(alfa)) # y-coordinate of the circle

        x = np.concatenate((x,x0))
        y = np.concatenate((y,y0))

        if (r < rc + 1e-8) and (r > rc - 1e-8):
            scatb_x = x0
            scatb_y = y0
        elif (r < rpmlin +1e-8) and (r > rpmlin - 1e-8):
            pmlbin_x = x0
            pmlbin_y = y0
        elif (r < rpmlout +1e-8) and (r > rpmlout - 1e-8):
            pmlbout_x = x0
            pmlbout_y = y0
    points = []
    for i,j in zip(x,y):
        points.append([i,j])
    elements = Delaunay(points)
    #Nodes
    scatb_nodes = []
    pmlbin_nodes = []
    pmlbout_nodes = []
    scatin_nodes = []
    scatin_points = []
    pmlin_nodes = []
    pmlin_points = []

    # Find boundary scattering nodes
    for i,j in zip(scatb_x,scatb_y):
        scatb_nodes.append(np.where((elements.points[:,0] == i) & (elements.points[:,1] == j))[0][0])
    
    # Find inside PML boundary nodes
    for i,j in zip(pmlbin_x,pmlbin_y):
        pmlbin_nodes.append(np.where((elements.points[:,0] == i) & (elements.points[:,1] == j))[0][0])
    
    # Find outside PML boundary nodes
    for i,j in zip(pmlbout_x,pmlbout_y):
        pmlbout_nodes.append(np.where((elements.points[:,0] == i) & (elements.points[:,1] == j))[0][0])

    #Find scattering and PML inside nodes
    for point in elements.points:
        if (np.sqrt(point[0]**2 + point[1]**2)) < rc - 1e-8:
            scatin_points.append(point)
        if ((np.sqrt(point[0]**2 + point[1]**2)) > rpmlin + 1e-8) & ((np.sqrt(point[0]**2 + point[1]**2)) < rpmlout - 1e-8) :
            pmlin_points.append(point)
    scatin_points = np.array(scatin_points)
    pmlin_points = np.array(pmlin_points)

    for i,j in zip(scatin_points[:,0],scatin_points[:,1]):
        scatin_nodes.append(np.where((elements.points[:,0] == i) & (elements.points[:,1] == j))[0][0])
    for i,j in zip(pmlin_points[:,0],pmlin_points[:,1]):
        pmlin_nodes.append(np.where((elements.points[:,0] == i) & (elements.points[:,1] == j))[0][0])

        

    return elements, scatb_nodes, scatin_nodes, pmlbin_nodes, pmlbout_nodes, pmlin_nodes
        
    







