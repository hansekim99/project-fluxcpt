import numpy as np
import pickle
from tqdm import tqdm
import num_index_density as idn
import matplotlib.pyplot as plt

# |%%--%%| <Jn6Ef3is6D|a4OAbF0QDE>

# remove nilpotent rays
def _nilpotent_begone(gvs, ray_gv_list):
    for ray_gv in ray_gv_list:
        ray = ray_gv[0]
        d = 1
        while tuple(d * ray) in gvs:
            del gvs[tuple(d * ray)]

    return gvs

# given a sample of moduli, evalute the scaling where the instanton correction upto specified degree is less than specified cutoff
def _cutoff_ev_qvs_moduli(gvs, qvs, moduli, cutoff, max_trials, tol):
    gvs = np.asarray(gvs, dtype=np.float64)
    qvs = np.asarray(qvs, dtype=np.float64)
    moduli = np.asarray(moduli, dtype=np.float64)
    
    N = moduli.shape[0]
    exponent = -2*np.pi*np.einsum("Nmd,Nd->Nm", qvs, moduli)

    def log_inst(l):
        weighted_exponent = l[:,None] * exponent
        return np.log(np.einsum("Nm,Nm->N", np.abs(gvs), np.exp(weighted_exponent)))

    l_lo, l_hi = 1.0 * np.ones((N,)), 2.0 * np.ones((N,))
    iev_lo, iev_hi = log_inst(l_lo), log_inst(l_hi)

    log_cutoff = np.log(float(cutoff))
    
    for _ in range(max_trials):
        l_hi_bool = (iev_hi < log_cutoff)
        if (~l_hi_bool).all():
            break
        else:
            l_hi[l_hi_bool] /= 1.1
            iev_hi = log_inst(l_hi)
    for _ in range(max_trials):
        l_lo_bool = (iev_lo > log_cutoff)
        if (~l_lo_bool).all():
            break
        else:
            l_lo[l_lo_bool] *= 1.1
            iev_lo = log_inst(l_lo)
    
    # print(f"lo : {moduli * l_lo}, hi : {moduli * l_hi}")
    for _ in range(max_trials):
        l_mid = 0.5 * (l_lo + l_hi)
        # print(f"l_mid : {moduli * l_mid}")
        iev_mid = log_inst(l_mid)

        nan_filter = np.isnan(iev_mid)
        l_lo[nan_filter], l_mid[nan_filter], l_hi[nan_filter] = 0, 0, 0

        err = iev_mid - log_cutoff
        
        l_mid_bool = (np.abs(err) <= tol / np.abs(iev_mid))
        l_mid_hi_bool = (err > tol / np.abs(iev_mid))
        l_mid_lo_bool = (err < -tol/ np.abs(iev_mid))
        # print(l_mid_bool, l_mid_hi_bool, l_mid_lo_bool)
        # print(err, tol/iev_mid)

        if l_mid_bool.all():
            break
        else:
            l_lo[l_mid_bool] = l_hi[l_mid_bool] = l_mid[l_mid_bool]
            l_hi[l_mid_hi_bool] = l_mid[l_mid_hi_bool]
            l_lo[l_mid_lo_bool] = l_mid[l_mid_lo_bool]
        
    scaled_moduli = moduli * l_mid[:,None]

    return scaled_moduli

def _cutoff_ev(moduli, gv_dict, cutoff_val, cutoff_max_trials, cutoff_tol):
    moduli = moduli / np.linalg.norm(moduli, axis = 1).reshape(-1,1)
    scaled_moduli = np.zeros_like(moduli, dtype=np.float64)

    N, d = moduli.shape[0], moduli.shape[1]
    
    for i in tqdm(range(N)):
        gvs_i = np.array(list(gv_dict.values()), dtype=np.float64)
        qvs_i = np.array(list(gv_dict.keys()), dtype=np.float64)
        
        gvs_i, qvs_i, moduli_i = gvs_i.reshape(1, -1), qvs_i.reshape(1, -1, d), moduli[i].reshape(1, d)

        scaled_moduli[i] = _cutoff_ev_qvs_moduli(gvs_i, qvs_i, moduli_i, cutoff_val, cutoff_max_trials, cutoff_tol)

    return scaled_moduli

def cutoff_dict(h_s, diffeo, moduli_sample_factor, moduli_max,
                cutoff_val, cutoff_max_trials, cutoff_tol,
                mode, gv_dict_mode):
    if mode == "save" or mode == "run":
        cutoff = {}

        for wall_data, cy_data in diffeo.items():
            ray_gv_list = cy_data.nop_ray_gv_list
            ls_obj = idn.CalabiYau(cy_data, moduli_max,
                           moduli_sample_factor, moduli_batch_no = int(1e1*10**h_s))
            moduli = ls_obj._moduli_projection_sample()
            
            if gv_dict_mode == "degree":
                gv_dict_default = cy_data.cutoff_gv_dict_deg
            gv_dict = _nilpotent_begone(gv_dict_default.dok, ray_gv_list)

            scaled_moduli = _cutoff_ev(moduli, gv_dict, cutoff_val, cutoff_max_trials, cutoff_tol)
            cutoff[wall_data] = scaled_moduli

    if mode == "save":
        with open(f"data/cone_moduli_cutoff/h_s={h_s}_cutoff.json", "wb") as f:
            pickle.dump(cutoff, f)

    elif mode == "load":
        with open(f"data/cone_moduli_cutoff/h_s={h_s}_cutoff.json", "rb") as f:
            cutoff = pickle.load(f)

    return cutoff


#|%%--%%| <a4OAbF0QDE|EouizVw9ra>

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
