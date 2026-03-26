from importlib import reload
from cytools import fetch_polytopes

h_s = 2
#p = Polytope([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-6,-9]])
h_s_polytope = fetch_polytopes(h11 = h_s, lattice = "N", limit = 10)

#|%%--%%| <BZ3AEmzjw0|dMTn4fSMHf>

import cone_flop_class
reload(cone_flop_class)

# find diffeomorphism classes (remove duplicate cy3s) and find flop facets
diffeo, flops = cone_flop_class.diffeo_class_flops(h_s_polytope, h_s, mode = "load")

# for h12, vv in diffeo.items():
#     for k, v in vv.items():
#         if h12 == 45:
#             print(h12, k, flops[h12][k][0])
#             print(flops[h12][k][1])

# check each flop facet for birational equivalent cy3 and collect into classes
birational = cone_flop_class.birational_class(diffeo, flops, h_s, mode = "load")

#|%%--%%| <dMTn4fSMHf|HPyrhOjuAV>

import cone_moduli_cutoff
reload(cone_moduli_cutoff)

# import moduli cutoff; for each in flop class generate moduli cutoff
cutoff = cone_moduli_cutoff.cutoff_dict(diffeo, flops, h_s, mode = "load")

# glue moduli cutoff diagrams together for each flop class
# cone_moduli_cutoff.scatter_plot_2d(birational, cutoff, h_s)
# cone_moduli_cutoff.scatter_plot_hyperplane(birational, cutoff, h_s, perp_vecs = [[0,0,1]], perp_coors = [1])

#|%%--%%| <HPyrhOjuAV|FrgfMI7sIs>

import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import num_index_density as idn
reload(idn)

# loop through all equiv classes; for each class loop through all cys

def scatter_2d_plot(h12, k):
    cy, points, rays = diffeo[h12][k][1], cutoff[k][0], cutoff[k][1]
    gv_flop = flops[h12][k][1][0] # multiple flops?
    cy_obj = idn.CalabiYau(cy, moduli_max = 10, moduli_sample_factor = int(5e5)) # 5e5

    plt.scatter(*np.array(points).T, s = 1, color = "red")
    plt.plot([0, 5*rays[0, 0]], [0, 5*rays[0, 1]], color="red")
    plt.plot([0, 5*rays[1, 0]], [0, 5*rays[1, 1]], color="red")

    return cy_obj, gv_flop

for h12, v in birational.items():
    for k, vv in v.items():
        if vv != [] and h12 == 74:
            # plt.cla()
            cy_obj, gv_flop = scatter_2d_plot(h12, k)
            moduli, scalar, npd_m = cy_obj.distr_rho(mrl = True, gv_flop = gv_flop)

            for fk in vv:
                flop_cy_obj, flop_gv_flop = scatter_2d_plot(h12, fk[0])
                flop_moduli, flop_scalar, flop_npd_m = flop_cy_obj.distr_rho(mrl = True, gv_flop = flop_gv_flop)

                moduli = np.concatenate((moduli, flop_moduli), axis=0)
                scalar = np.concatenate((scalar, flop_scalar), axis=0)
                npd_m = np.concatenate((npd_m, flop_npd_m), axis=0)

            mask = np.nonzero(np.abs(scalar))
            moduli, scalar = moduli[mask], scalar[mask]

            lim = 10

            plt.xlim((-lim,lim))
            plt.ylim((-lim,lim))

            grid_x, grid_y = np.linspace(-lim,lim,100), np.linspace(-lim,lim,100)
            grid_x, grid_y = np.meshgrid(grid_x, grid_y)
            grid_z = griddata((moduli[:,0], moduli[:,1]), np.log(np.abs(scalar)), (grid_x, grid_y), method='cubic')
            tree = cKDTree(np.c_[moduli[:,0], moduli[:,1]])
            distances, _ = tree.query(np.c_[grid_x.ravel(), grid_x.ravel()])
            distances = distances.reshape(grid_x.shape)
            grid_z_m = np.ma.masked_where(distances > 0.5, grid_z)

            contour = plt.contour(grid_x, grid_y, grid_z_m, levels=20, cmap='viridis')
            plt.colorbar(contour, label = "log rho")

            plt.title(f"Extended Kahler cone and GV cone for 2-face equiv. CY3s with h12 = {h12}\nContour map of index density")
            plt.gcf().set_size_inches(18.5, 10.5)
            plt.savefig(f"figures/num_index_density/h_s={h_s}_h12={h12}_ctr")
            plt.close()
            
            plt.xlim((-lim,lim))
            plt.ylim((-lim,lim))

            scatter = plt.scatter(moduli[:,0], moduli[:,1], c = np.log(np.abs(scalar)), s = 1)
            plt.colorbar(scatter, label = "log rho")

            plt.scatter(npd_m[:,0], npd_m[:,1], c = "black", s = 1)

            plt.title(f"Extended Kahler cone and GV cone for 2-face equiv. CY3s with h12 = {h12}\nColour map of index density")
            plt.gcf().set_size_inches(18.5, 10.5)

            _, _ = scatter_2d_plot(h12, k)

            for fk in vv:
                _, _ = scatter_2d_plot(h12, fk[0])

            plt.savefig(f"figures/num_index_density/h_s={h_s}_h12={h12}_col")
            plt.close()

            break


