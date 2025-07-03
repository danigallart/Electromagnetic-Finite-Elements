
import numpy as np
from meshgen_circle import meshgen_circle, plot_mesh_circleT
from matplotlib import pyplot as plt
from element_matrix_triangle import element_matrix_triangle
from pml import pml

# Constants
c0  = 3*1e8              # m/sec, velocity of light in free space
nu0 = 120*np.pi          # ohm, intrinsic impedance of the free space
e0  = (1e-9)/(36*np.pi)  # F/m, permittivity of free space
mu0 = 4*np.pi*1e-7       # H/m, permeability of free space

# Input parameters
freq = 300               # MHz frequency
freq = freq*1e6          # Hz frequency
lambda0 = c0/freq        # meter wavelength
k0 = 2*np.pi/lambda0     # 1/meter wavenumber
omg = 2*np.pi*freq       # rad/s radial frequency

# Polarization 
polrz = 'tm' # TM or TE polarization

mur = 1   # relative permittivity of anisotropic object
erxx = 9  # xx comp. of dielectric constant of anisotropic object
eryy = 4  # yy comp. of dielectric constant of anisotropic object
erzz = 2  # zzz comp. of dielectric constant of anisotropic object
erxy = 0  # xy comp. of dielectric constant of anisotropic object
erxz = 0  # xz comp. of dielectric constant of anisotropic object
eryx = 0  # yx comp. of dielectric constant of anisotropic object
eryz = 0  # yz comp. of dielectric constant of anisotropic object
erzx = 0  # zx comp. of dielectric constant of anisotropic object
erzy = 0  # zy comp. of dielectric constant of anisotropic object

er = np.array([[erxx, erxy, erxz], [eryx, eryy, eryz], [erzx, erzy, erzz]])
er = 5

phii = 0                # angle of incident field
phii = phii*np.pi/180   # radian of incident field
scat_type = 'diel'

# Mesh Generation
class structure:
    def __init__(self, rc, delh, scat_type, scat_shape):
        self.rc = rc
        self.delh = delh
        self.scat_type = scat_type
        self.scat_shape = scat_shape

delh = lambda0/40
scat_shape = 1

print('Creating the mesh')
rc = 1*lambda0    # radius of circular scatterer
obj = structure(rc, delh, scat_type, scat_shape)

elements, scatb_nodes, scatb_elm, scatin_nodes, scatin_elm, pmlbin_nodes, pmlbout_nodes, pmlin_nodes = meshgen_circle(obj,freq)


N = len(elements.points) # Number of nodes
M = len(elements.simplices) # Number of elements
conn = elements.simplices

# Material and Source parameters
epsr = np.ones(M)
epsr[scatin_elm] = er


#Relative permittivity of elements
epsrxx = np.ones(M)
epsryy = np.ones(M)
epsrzz = np.ones(M)
epsrxy = np.zeros(M)
epsryx = np.zeros(M)

epsrxx[scatin_elm] = erxx
epsryy[scatin_elm] = eryy
epsrzz[scatin_elm] = erzz
epsrxy[scatin_elm] = erxy
epsryx[scatin_elm] = eryx

x = elements.points[:,0]
x = x.astype(np.complex128)
y = elements.points[:,1]
y = y.astype(np.complex128)
#Implementation of the PML

pmlbin_x = x[pmlbin_nodes]
pmlbin_y = y[pmlbin_nodes]
pmlbout_x = x[pmlbout_nodes]
pmlbout_y = y[pmlbout_nodes]

print('The LC-PML is being created')
counter = 0
x_p = []
y_p = []
for i in pmlin_nodes:
    xc, yc, comptador, flag = pml(x[i],y[i],k0,pmlbin_x,pmlbin_y,pmlbout_x,pmlbout_y)
    if flag:
        x_p.append(x[i])
        y_p.append(y[i]) 
    x[i] = xc[0]
    y[i] = yc[0]
    counter += comptador
    #print(xc)
print('The LC-PML is finished', counter)

plot_mesh_circleT(elements, scatb_nodes, pmlbin_nodes, pmlbout_nodes, scatin_nodes, pmlin_nodes,x_p,y_p)

#Main Body
#Prepare the assembly of the matrix

Ne = 3 #Number of nodes in each triangular element (linear triangular elements)
A = np.zeros((N,N), dtype = 'complex_') #The global matrix A
b = np.zeros((N,1), dtype = 'complex_') #The global vector b

#Assemble the matrix

for e in range(M):
    if polrz == 'tm':
        pe = 1/mur
        qe = -k0**2*epsr[e]
    else:
        pe = 1/epsr[e]
        qe = -k0**2*mur
        
    Ae, be, area_element, Jdet, x1, x2, x3, y1, y2, y3 = element_matrix_triangle(e, x, y, conn, pe, pe, qe, 0) #Creation of element matrices
    #print('e:', e, 'Ae:', Ae)
#Assembly of the matrix
    for i in range(Ne):
        ig = conn[e,i]  #global node corresponding to i node in element e
        for j in range(Ne):
            jg = conn[e,j] #global node corresponding to j node in element e
            A[ig,jg] += Ae[i,j]
            
    
uinc = np.exp(1j*k0*(np.real(x)*np.cos(phii) + np.real(y)*np.sin(phii)))
uinc = uinc.reshape(-1,1)

Au = np.matmul(A,uinc).reshape(-1,1)
b[scatin_nodes] = -Au[scatin_nodes] 

u = np.linalg.solve(A,b)
utot = u + uinc
utot[pmlin_nodes] = 0.
#utot[pmlbout_nodes] = 0.

sc = plt.scatter(x, y, c=np.abs(utot), cmap='jet')#, vmin=0., vmax=5)  # Set range from 0 to 1
plt.colorbar(sc, label='u value')
plt.show()
