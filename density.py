from operator import index
import numpy as np
#Enable unfavourable CYs
#cytools.config.enable_experimental_features()
from cytools import fetch_polytopes, Polytope

from importlib import reload
import index_density_num as idn
reload(idn)

# |%%--%%| <Jn6Ef3is6D|yYbCQ6fuay>

h_s = 2
#p = Polytope([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-6,-9]])
h_s_polytope = fetch_polytopes(h11 = h_s, lattice = "N", limit = 100)

l = len(h_s_polytope)

#|%%--%%| <yYbCQ6fuay|qPUnwXUxOI>

for i in range(l):
    p = h_s_polytope[i]
    cy = p.triangulate().get_cy()
    dictK = cy.intersection_numbers(in_basis = True)
    print(dictK)
    
    cone_hyperplane = cy.toric_mori_cone(in_basis = True).extremal_rays()
    #print('kahler rays : \n', cy.toric_kahler_cone().extremal_rays())
    #print('kahler planes : \n', cy.toric_mori_cone(in_basis = True).extremal_rays())
    arrayK = np.array([[[dictK.get(tuple(sorted((i,j,k))), 0) for i in range(h_s)] for j in range(h_s)] for k in range(h_s)])
    
    cy_obj = idn.CalabiYau(h_s, arrayK, cone_hyperplane,
                   moduli_max = 10, moduli_cutoff = 1, qd3 = 50,
                   moduli_sample_no = int(5e5))
    
    _, udv_cmbn = cy_obj.uniform_eval(mrl = True)
    usv_cmbn = cy_obj.uniform_integrate(udv_cmbn)

    print(usv_cmbn)
