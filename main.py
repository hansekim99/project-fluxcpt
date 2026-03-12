from importlib import reload
from cytools import fetch_polytopes

h_s = 3
#p = Polytope([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-6,-9]])
h_s_polytope = fetch_polytopes(h11 = h_s, lattice = "N", limit = 100)

#|%%--%%| <e0PJIeuUXz|JKT5eUNmJ6>

import cone_flop_class
reload(cone_flop_class)

# find diffeomorphism classes (remove duplicate cy3s) and find flop facets
diffeo, flops = cone_flop_class.diffeo_class(h_s_polytope, h_s, mode = "load")

# check each flop facet for birational equivalent cy3 and collect into classes
birational = cone_flop_class.birational_class(diffeo, flops, h_s, mode = "load")

#|%%--%%| <JKT5eUNmJ6|L20gYK9d7V>

import cone_moduli_cutoff
reload(cone_moduli_cutoff)

# import moduli cutoff; for each in flop class generate moduli cutoff
cutoff = cone_moduli_cutoff.cutoff_dict(diffeo, flops, h_s, mode = "load")

# glue moduli cutoff diagrams together for each flop class
# cone_moduli_cutoff.scatter_plot_2d(birational, cutoff, h_s)

#|%%--%%| <L20gYK9d7V|CqL3223mdD\>

import num_index_density as idn
reload(idn)

# loop through all equiv classes; for each class loop through all cys


# for i in range(1):
#     p = h_s_polytope[i]
#     print(p.triangulate().get_cy().toric_kahler_cone().extremal_rays())
#         
#     cy_obj = idn.CalabiYau(p, moduli_max = 4, moduli_sample_factor = int(1e6), moduli_batch_no = int(1e3))
#
#     moduli_distr, _, usv_cmbn, usv_se = cy_obj.integ_rho(sample_type = "projection", mrl = True)
#
#     print(usv_cmbn)

