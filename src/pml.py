import numpy as np
from decimal import Decimal

def pml(x, y, k, pmlbin_x, pmlbin_y, pmlbout_x, pmlbout_y):
    alpha = 7.0*k
    alpha_complex = alpha/(1j*k)
    m = 3 #PML decay rate
    comptador = 0
    flag = False
    x = np.array(x)
    y = np.array(y)
    pmlbin_x = np.array(pmlbin_x)
    pmlbin_y = np.array(pmlbin_y)
    pmlbout_x = np.array(pmlbout_x)
    pmlbout_y = np.array(pmlbout_y)

    
    # Find the point in pmlbin which is closest to the point in the PML: (x,y)
    distance = (pmlbin_x - x)**2 + (pmlbin_y - y)**2
    ksi = np.sqrt(distance.min())
    minimum_distance_index = np.where(distance == distance.min())[0]
    if len(minimum_distance_index) > 1:
        minimum_distance_index = minimum_distance_index[0]
        flag = True
        comptador = 1
 
    x0 = pmlbin_x[minimum_distance_index]
    y0 = pmlbin_y[minimum_distance_index]

    # Find the point in pmlbout in the direction of the unit vector between (x,y) and (pmlbin_x,pmlbin_y)

    vpx = pmlbout_x - x0
    vpy = pmlbout_y - y0
    #print('len of vpx', len(vpx), '\n', 'len of vpy', len(vpy))
    lp = np.sqrt(vpx**2 + vpy**2)
    upx = vpx/lp
    upy = vpy/lp

    vx = x - x0
    vy = y - y0
    l = np.sqrt(vx**2 + vy**2)

    if l < 1.e8:
        xc = x
        yc = y
        return xc, yc, comptador, flag
    else:    
        uvx = vx/l
        uvy = vy/l

        scalar_product = (upx*uvx + upy*uvy)
        maximum_product_index = np.where(scalar_product == scalar_product.max())[0]

        x1 = x[maximum_product_index] #PML bout x
        y1 = y[maximum_product_index] #PML bout y

        dpml = np.sqrt((x1-x0)**2 + (y1-y0)**2) #Local PML thickness
        term = alpha_complex*((ksi**m) / (m*(dpml**(m-1))))

        xc = x + term*uvx #Complex coordinate x
        yc = y + term*uvy #Comples coordinate y

        return xc, yc, comptador, flag









