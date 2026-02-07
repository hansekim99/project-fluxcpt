from operator import index
import numpy as np
#Enable unfavourable CYs
#cytools.config.enable_experimental_features()
from cytools import fetch_polytopes, Polytope

from importlib import reload
import index_density_num as idn
reload(idn)

# |%%--%%| <Jn6Ef3is6D|Rd7f0EjWS0>

h21 = 3
#p = Polytope([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-6,-9]])
h21_polytope = fetch_polytopes(h11 = h21, lattice = "N", limit = 100)

# |%%--%%| <Rd7f0EjWS0|vNxZpeMhoy>

def index_density(p):
    #print(p)
    cy = p.triangulate().get_cy()
    dictK = cy.intersection_numbers(in_basis = True)
    print(dictK)

    mori_cone = cy.toric_mori_cone(in_basis = True)
    kahler_rays = np.linalg.inv(mori_cone.extremal_rays())
    arrayK = np.array([[[dictK.get(tuple(sorted((i,j,k))), 0) for i in range(h21)] for j in range(h21)] for k in range(h21)])
    
    cy_obj = idn.CalabiYau(h21, arrayK, kahler_rays,
                   moduli_interval = (1,100), qd3 = 50,
                   moduli_sample_no = int(5e6))

    _, udv_cmbn = cy_obj.uniform_eval(type = "cmbn", mrl = True)

    usv_cmbn = cy_obj.uniform_integrate(udv_cmbn)

    return usv_cmbn

# |%%--%%| <vNxZpeMhoy|kDkRakQzl4>

l = len(h21_polytope)

for i in range(5):
    print(index_density(h21_polytope[i]))
    #p = h21_polytope[i]
    #cy = p.triangulate().get_cy()

    #print(cy.h21())
