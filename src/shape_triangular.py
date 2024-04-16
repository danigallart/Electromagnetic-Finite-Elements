import scipy as sp
from scipy import io
import meshpy.triangle as triangle
#%matplotlib inline
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as pt
def round_trip_connect(start, end):
    result = []
    for i in range(start, end):
        result.append((i, i+1))
    result.append((end, start))
    return result

# Points list defines boundary
points = []
points.extend((np.cos(angle),  np.sin(angle)) for angle in np.linspace(0, 2*np.pi, 30, endpoint=False))


info = triangle.MeshInfo()
info.set_points(points)
info.set_facets(round_trip_connect(0, len(points)-1))

#Here you build mesh
mesh = triangle.build(info, max_volume=1e-2, min_angle=25)

#You can save it to file:
#write_dat(mesh,round_trip_connect(0, len(points)-1),'circle.dat')
#triangle.write_gnuplot_mesh("circle.dat", mesh)

mesh_points = np.array(mesh.points)
mesh_tris = np.array(mesh.elements)
mesh_facets = np.array(mesh.facets)


pt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris)
pt.show()