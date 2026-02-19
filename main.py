import numpy as np
#Enable unfavourable CYs
#cytools.config.enable_experimental_features()
from cytools import fetch_polytopes, Polytope

from importlib import reload
import num_index_density as idn
reload(idn)

# |%%--%%| <Jn6Ef3is6D|yYbCQ6fuay>

h_s = 2
#p = Polytope([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-6,-9]])
h_s_polytope = fetch_polytopes(h11 = h_s, lattice = "N", limit = 100)

l = len(h_s_polytope)

#|%%--%%| <yYbCQ6fuay|qPUnwXUxOI>

import matplotlib.pyplot as plt

for i in range(1):
    p = h_s_polytope[i]
    print(p.triangulate().get_cy().toric_kahler_cone().extremal_rays())
        
    cy_obj = idn.CalabiYau(p, moduli_max = 4, moduli_sample_factor = int(1e6), moduli_batch_no = int(1e3))

    moduli_distr, _, usv_cmbn, usv_se = cy_obj.integ_rho(sample_type = "projection", mrl = True)

    print(usv_cmbn)

