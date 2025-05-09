
import numpy as np
from meshgen_circle import meshgen_circle
from matplotlib import pyplot as plt

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

delh = lambda0/20
scat_shape = 1

print('Creating the mesh')
rc = 1*lambda0    # radius of circular scatterer
obj = structure(rc, delh, scat_type, scat_shape)

elements, scatb_nodes, scatin_nodes, pmlbin_nodes, pmlbout_nodes, pmlin_nodes = meshgen_circle(obj,freq)


N = len(elements.points) # Number of nodes
M = len(elements.simplices) # Number of elements

plt.scatter(elements.points[:,0],elements.points[:,1])
plt.scatter(elements.points[scatb_nodes][:,0], elements.points[scatb_nodes][:,1])
plt.scatter(elements.points[pmlbin_nodes][:,0],elements.points[pmlbin_nodes][:,1])
plt.scatter(elements.points[pmlbout_nodes][:,0],elements.points[pmlbout_nodes][:,1])
plt.scatter(elements.points[scatin_nodes][:,0], elements.points[scatin_nodes][:,1])
plt.scatter(elements.points[pmlin_nodes][:,0], elements.points[pmlin_nodes][:,1], c='black')
#plt.triplot(elements.points[:,0],elements.points[:,1], elements.simplices)
plt.show()