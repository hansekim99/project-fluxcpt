import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

def scatter_plot_2d(diffeo, cutoff, h_s):
    for wd in diffeo.keys():
        plt.cla()
        plt.xlim((-5,5))
        plt.ylim((-5,5))
        
        cutoff_pts, cyd = cutoff[wd], diffeo[wd]
        h_l = cyd.h_l

        kahler_rays = cyd.kahler_extremal_rays.tolist()
        birational_wall_data = cyd.birational_class_wall_data_list

        plt.scatter(*np.array(cutoff_pts).T, s = 1)

        for flop_wd in birational_wall_data:
            flop_cutoff_pts = cutoff[flop_wd]
            plt.scatter(*np.array(flop_cutoff_pts).T, s = 1)

            flop_cyd = diffeo[flop_wd]
            flop_kahler_rays = flop_cyd.kahler_extremal_rays

            for flop_ray in flop_kahler_rays:
                kahler_rays.append(flop_ray)

        counter_rays = Counter(tuple(ray) for ray in kahler_rays)

        for ray, count in counter_rays.items():
            if count > 1:
                plt.plot([0, 5*ray[0]], [0, 5*ray[1]], color="red", ls="dashed")
            else:
                plt.plot([0, 5*ray[0]], [0, 5*ray[1]], color="red")
        
        plt.title(f"Extended Kahler cone and GV cone for 2-face equiv. CY3s with h12 = {h_l}")
        plt.savefig(f"figures/cone_moduli_cutoff/h_s={h_s}_h12={h_l}")

def scatter_plot_hyperplane(diffeo, cutoff, h_s, perp_vecs = None, perp_coors = None, crit = 0.1):
    perp_vecs, perp_coors = np.array(perp_vecs), np.array(perp_coors)
    def filter_points(pts):
        perp_ip = np.einsum('na,ba->nb', pts, perp_vecs)
        mask = np.all((perp_ip > perp_coors) & (perp_ip < perp_coors + crit), axis = 1)
        return pts[mask]

    for wd in diffeo.keys():
        plt.cla()
        plt.xlim((-5,5))
        plt.ylim((-5,5))

        cutoff_pts, cyd = cutoff[wd], diffeo[wd]
        h_l = cyd.h_l
        filtered_pts = filter_points(cutoff_pts)

        kahler_rays = cyd.kahler_extremal_rays
        birational_wall_data = cyd.birational_class_wall_data_list

        plt.scatter(*np.array(filtered_pts).T)
        
        for flop_wd in birational_wall_data:
            flop_cutoff_pts, flop_cyd = cutoff[flop_wd], diffeo[flop_wd]
            flop_filtered_pts = filter_points(flop_cutoff_pts)

            flop_rays = flop_cyd.kahler_extremal_rays
            plt.scatter(*np.array(flop_filtered_pts).T)
        
        plt.title(f"Extended Kahler cone and GV cone for 2-face equiv. CY3s with h11 = {h_l}\n 2-plane w. distance {perp_coors} along {perp_vecs}; tolerance {crit}")
        plt.savefig(f"figures/cone_moduli_cutoff/h_s={h_s}_h12={h_l}_vecs={str(perp_vecs)}_coors={str(perp_coors)}")

#|%%--%%| <QImUrKVOhc|BqO8y4VNTM>

from collections import Counter
from scipy.interpolate import griddata
from scipy.spatial import cKDTree

def contour_plot_2d(h_s, diffeo, cutoff, idd,
                    plot_size):
    for wd in diffeo.keys():
        plt.cla()
        plt.xlim((-plot_size, plot_size))
        plt.ylim((-plot_size, plot_size))
        
        cutoff_pts, cyd = cutoff[wd], diffeo[wd]
        scalar, moduli = idd[wd][1], idd[wd][2]
        # scalar_nc, moduli_nc = idd_non_corr[wd][1], idd_non_corr[wd][2]
        h_l = cyd.h_l

        kahler_rays = cyd.kahler_extremal_rays.tolist()
        birational_wall_data = cyd.birational_class_wall_data_list
        
        grid_x, grid_y = np.linspace(-plot_size,plot_size,100), np.linspace(-plot_size,plot_size,100)
        grid_x, grid_y = np.meshgrid(grid_x, grid_y)    
        
        for flop_wd in birational_wall_data:
            flop_cutoff_pts = cutoff[flop_wd]
            flop_scalar, flop_moduli = idd[flop_wd][1], idd[flop_wd][2]
            
            scalar = np.concatenate((scalar, flop_scalar))
            moduli = np.concatenate((moduli, flop_moduli))
            
            plt.scatter(*np.array(flop_cutoff_pts).T, s = 1)

            flop_cyd = diffeo[flop_wd]
            flop_kahler_rays = flop_cyd.kahler_extremal_rays

            for flop_ray in flop_kahler_rays:
                kahler_rays.append(flop_ray)
        
        plt.scatter(*np.array(cutoff_pts).T, s = 1)

        scalar_nz_mask = np.nonzero(scalar)
        scalar = scalar[scalar_nz_mask]
        moduli = moduli[scalar_nz_mask]
        log_abs_scalar = np.log(np.abs(scalar))
            
        grid_z = griddata((moduli[:,0], moduli[:,1]), log_abs_scalar, (grid_x, grid_y), method='cubic')
        # tree = cKDTree(np.c_[moduli[:,0], moduli[:,1]])
        # distances, _ = tree.query(np.c_[grid_x.ravel(), grid_x.ravel()])
        # distances = distances.reshape(grid_x.shape)
        # grid_z_m = np.ma.masked_where(distances > 0.5, grid_z)

        contour = plt.contour(grid_x, grid_y, grid_z, levels=20, cmap='viridis')
        plt.colorbar(contour, label = "log rho")

        counter_rays = Counter(tuple(ray) for ray in kahler_rays)

        for ray, count in counter_rays.items():
            if count > 1:
                plt.plot([0, 5*ray[0]], [0, 5*ray[1]], color="red", ls="dashed")
            else:
                plt.plot([0, 5*ray[0]], [0, 5*ray[1]], color="red")
        plt.scatter(*np.array(cutoff_pts).T, s = 1)
        
        plt.title(f"Extended Kahler cone and GV cone for 2-face equiv. CY3s with h12 = {h_l}")
        plt.gcf().set_size_inches(18.5, 10.5)
        plt.savefig(f"figures/num_index_density/h_s={h_s}_h12={h_l}_ctr")
        plt.close()

def scatter_plot_2d(h_s, diffeo, cutoff, idd,
                    plot_size):
    for wd in diffeo.keys():
        plt.cla()
        plt.xlim((-plot_size, plot_size))
        plt.ylim((-plot_size, plot_size))
        
        cutoff_pts, cyd = cutoff[wd], diffeo[wd]
        scalar, moduli = idd[wd][1], idd[wd][2]
        # scalar_nc, moduli_nc = idd_non_corr[wd][1], idd_non_corr[wd][2]
        h_l = cyd.h_l

        kahler_rays = cyd.kahler_extremal_rays.tolist()
        birational_wall_data = cyd.birational_class_wall_data_list
        
        for flop_wd in birational_wall_data:
            flop_cutoff_pts = cutoff[flop_wd]
            flop_scalar, flop_moduli = idd[flop_wd][1], idd[flop_wd][2]
            
            scalar = np.concatenate((scalar, flop_scalar))
            moduli = np.concatenate((moduli, flop_moduli))
            cutoff_pts = np.concatenate((cutoff_pts, flop_cutoff_pts))

            flop_cyd = diffeo[flop_wd]
            flop_kahler_rays = flop_cyd.kahler_extremal_rays

            for flop_ray in flop_kahler_rays:
                kahler_rays.append(flop_ray)
            
        scalar_nz_mask = np.nonzero(scalar)
        scalar = scalar[scalar_nz_mask]
        moduli = moduli[scalar_nz_mask]
        log_abs_scalar = np.log(np.abs(scalar))

        scatter = plt.scatter(moduli[:,0], moduli[:,1], c = log_abs_scalar, s = 1)
        plt.colorbar(scatter, label = "log rho")
        
        plt.scatter(*np.array(cutoff_pts).T, s = 1)

        counter_rays = Counter(tuple(ray) for ray in kahler_rays)

        for ray, count in counter_rays.items():
            if count > 1:
                plt.plot([0, 5*ray[0]], [0, 5*ray[1]], color="red", ls="dashed")
            else:
                plt.plot([0, 5*ray[0]], [0, 5*ray[1]], color="red")
        plt.scatter(*np.array(cutoff_pts).T, s = 1)
        
        plt.title(f"Extended Kahler cone and GV cone for 2-face equiv. CY3s with h12 = {h_l}")
        plt.gcf().set_size_inches(18.5, 10.5)
        plt.savefig(f"figures/num_index_density/h_s={h_s}_h12={h_l}_sct")
        plt.close()


