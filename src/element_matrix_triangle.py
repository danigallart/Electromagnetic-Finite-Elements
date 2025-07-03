import numpy as np

def element_matrix_triangle(e, x, y, conn, pxe, pye, qe, fe):
    """Computation of the element A_e, b_e and area_e"""
    #e: element ID
    #x,y: x and y coordinates of all nodes in the mesh
    #conn: connectivity matrix
    #pxe: px value of the eth element
    #pye: py ""
    #qe: q ""
    #fe: f ""

    Ne = 3 #Number of nodes in an element
    Ae1 = np.zeros((Ne,Ne), dtype='complex_')
    Ae2 = np.zeros((Ne,Ne), dtype='complex_')
    be = np.zeros((Ne,1), dtype='complex_')

    n1 = conn[e, 0] #Node 1 of element e
    n2 = conn[e, 1] #Node 2 of element e
    n3 = conn[e, 2] #Node 3 of element e
    
    x1 = x[n1] #x coordinate of node 1
    x2 = x[n2] #x coordinate of node 2
    x3 = x[n3] #x coordinate of node 3

    y1 = y[n1] #y coordinate of node 1
    y2 = y[n2] #y coordinate of node 2
    y3 = y[n3] #y coordinate of node 3

    #Jacobian entries

    dxdpsi = x2 - x1
    dxdeta = x3 - x1
    dydpsi = y2 - y1
    dydeta = y3 - y1

    Jdet = dxdpsi*dydeta - dxdeta*dydpsi
    
    area_element = 0.5*Jdet

    #Entries of the inverse Jacobian matrix
    Jinv11 = (y3-y1)/Jdet
    Jinv12 = -(y2 - y1)/Jdet
    Jinv21 = -(x3 - x1)/Jdet
    Jinv22 = (x2 - x1)/Jdet

    #Derivative over the shape functions
    dNijdpsi = np.array([-1, 1, 0])
    dNijdeta = np.array([-1, 0, 1])

    #Compute the first part of the element matrix
    for i in range(Ne):
        for j in range(i,Ne):
            Ae1[i,j] = (pxe*(Jinv11*dNijdpsi[i] + Jinv12*dNijdeta[i]) * (Jinv11*dNijdpsi[j]+Jinv12*dNijdeta[j]) +\
                       pye*(Jinv21*dNijdpsi[i] + Jinv22*dNijdeta[i]) * (Jinv21*dNijdpsi[j]+Jinv22*dNijdeta[j])) *\
                       Jdet/2.0       
            Ae1[j,i] = Ae1[i,j]
            

    #Compute the second part of the element matrix
    gauss_weights = np.array([1./6., 1./6., 1./6.])
    psi = np.array([0.5, 0, 0.5])
    eta = np.array([0, 0.5, 0.5])

    for i in range(Ne):
        for j in range(i,Ne):
            sumg = 0.
            for k in range(len(gauss_weights)):
                if i == 0:
                    Ni = 1 - psi[k] - eta[k]
                elif i == 1:
                    Ni = psi[k]
                elif i == 2:
                    Ni = eta[k]
                if j == 0:
                    Nj = 1 - psi[k] - eta[k]
                elif j == 1:
                    Nj = psi[k]
                elif j == 2:
                    Nj = eta[k]
                sumg += gauss_weights[k]*qe*Ni*Nj*Jdet
            Ae2[i,j] = sumg
            Ae2[j,i] = Ae2[i,j]
    
    #Element matrix
    Ae = Ae1 + Ae2

    #Right hand side vector
    for i in range(Ne):
        sumg = 0.
        for k in range(len(gauss_weights)):
            if i == 0:
                Ni = 1 - psi[k] - eta[k]
            elif i == 1:
                Ni = psi[k]
            elif i == 2:
                Ni = eta[k]
            sumg += gauss_weights[k]*Ni*fe*Jdet
        be[i] = sumg
    return Ae, be, area_element, Jdet, x1, x2, x3, y1, y2, y3


