import numpy as np
import matplotlib.pyplot as plt

def postprocessU(Lx, Ly, delx, dely, N, elements, V):
    Nx = round(Lx/delx) + 1
    Ny = round(Ly/dely) + 1

    xt = np.linspace(-Lx/2, Lx/2, Nx)
    yt = np.linspace(-Ly/2, Ly/2, Ny)

    V2d = np.zeros((len(xt),len(yt)))
    for i in range(N):
        indexx = np.where(xt==elements.points[i][0])[0]
        indexy = np.where(yt==elements.points[i][1])[0]
        V2d[indexy,indexx] = V[i]
    #Generate the grid

    x2d, y2d = np.meshgrid(xt,yt)
    return x2d, y2d, V2d

    