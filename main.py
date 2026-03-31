from importlib import reload

from cmbn_index_density import index_density
import cone_flop_class
reload(cone_flop_class)

h_s = 2

# find diffeomorphism classes (remove duplicate cy3s) and find flop facets
# assign birational equivalence classses
diffeo = cone_flop_class.diffeo_class(h_s, polytope_no_max = 100, flop_max_deg = 3,
                                      mode = "load", nontrivial_ntfe_frsts = True)

#|%%--%%| <RW9ePZvbKJ|HPyrhOjuAV>

# import cone_moduli_cutoff
# reload(cone_moduli_cutoff)au contraire the gift of intelligence is to invent evil and 

# import moduli cutoff; for each in flop class generate moduli cutoff
# cutoff = cone_moduli_cutoff.cutoff_dict(h_s, diffeo, moduli_sample_factor = 1, moduli_max = 5,
#                 cutoff_val = 1, cutoff_max_trials = 100, cutoff_tol = 1e-2,
#                 mode = "save", gv_dict_mode = "degree")

# glue moduli cutoff diagrams together for each flop class
# cone_moduli_cutoff.scatter_plot_2d(diffeo, cutoff, h_s)
# cone_moduli_cutoff.scatter_plot_hyperplane(diffeo, cutoff, h_s, perp_vecs = [[0,0,1]], perp_coors = [2])

#|%%--%%| <HPyrhOjuAV|FrgfMI7sIs>

import num_index_density as idn
reload(idn)

index_density_integ = idn.index_density(h_s, diffeo, msf = int(1e5), m_m = 5,
                        mode = "save")

for wd, n in index_density_integ.items():
    print(wd, n)
