import numpy as np
import pickle
from tqdm import tqdm
import num_index_density as idn
import matplotlib.pyplot as plt

# |%%--%%| <Jn6Ef3is6D|hPiCUWr4nw>

# remove nilpotent rays
def nilpotent_begone(gvs, null_rays):
    for ray in null_rays:
        d = 1
        print(ray)
        while tuple(d * ray[0]) in gvs:
            del gvs[tuple(d * ray[0])]

    return gvs


# given a sample of moduli, evalute the scaling where the instanton correction upto specified degree is less than specified cutoff
    # if parallel, gv invariants are stored in memory
def cutoff_ev_qvs_moduli(gvs, qvs, moduli, cutoff, max_trials, tol):
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

def cutoff_ev(cy, null_rays, h_s, min_points = int(2e1), cutoff = 1, tol = 1e-3, max_trials = 100): 
    cy_obj = idn.CalabiYau(cy, moduli_max = 5,
                           moduli_sample_factor = int(1), moduli_batch_no = int(1e1*10**h_s))
    moduli = cy_obj._moduli_projection_sample()
    moduli = moduli / np.linalg.norm(moduli, axis = 1).reshape(-1,1)

    dual_rays = cy.mori_cone_cap(in_basis = True).extremal_rays()

    scaled_moduli = np.zeros_like(moduli, dtype=np.float64)

    N, d = moduli.shape[0], moduli.shape[1]

    gv_dict_default = cy.compute_gvs(min_points = round(min_points*1.2)).dok
    gv_dict = nilpotent_begone(gv_dict_default, null_rays)
    
    for i in tqdm(range(N)):
        gvs_i = np.array(list(gv_dict.values())[:min_points], dtype=np.float64)
        qvs_i = np.array(list(gv_dict.keys())[:min_points], dtype=np.float64)
        
        gvs_i, qvs_i, moduli_i = gvs_i.reshape(1, -1), qvs_i.reshape(1, -1, d), moduli[i].reshape(1, d)

        scaled_moduli[i] = cutoff_ev_qvs_moduli(gvs_i, qvs_i, moduli_i, cutoff, max_trials, tol)
        
    return scaled_moduli

def _cutoff_dict(diffeo, flops, h_s):
    cutoff = {}

    for h12, diffs in diffeo.items():
        for k, v in diffs.items():
            cy = v[1] # diffeo[h12][k][1]
            null_rays = flops[h12][k][1]
            scaled_moduli = cutoff_ev(cy, null_rays, h_s)
            rays = cy.toric_kahler_cone().extremal_rays()
            cutoff[k] = [scaled_moduli, rays]
            # cutoff[k] = scaled_moduli

    return cutoff

def cutoff_dict(diffeo, flops, h_s, mode):
    if mode == "save":
        cutoff = _cutoff_dict(diffeo, flops, h_s)
        with open(f"data/cone_moduli_cutoff/h_s={h_s}_cutoff.json", "wb") as f:
            pickle.dump(cutoff, f)

    elif mode == "load":
        with open(f"data/cone_moduli_cutoff/h_s={h_s}_cutoff.json", "rb") as f:
            cutoff = pickle.load(f)

    return cutoff

def scatter_plot_2d(birational, cutoff, h_s):
    for h12, v0 in birational.items():
        for k, v in v0.items():
            if v == []:
                continue
            else:
                plt.cla()
                plt.xlim((-5,5))
                plt.ylim((-5,5))
                
                cutoff_pts = cutoff[k][0]
                rays = cutoff[k][1]
                plt.scatter(*np.array(cutoff_pts).T, s = 1)
                plt.plot([0, 5*rays[0, 0]], [0, 5*rays[0, 1]], color="red")
                plt.plot([0, 5*rays[1, 0]], [0, 5*rays[1, 1]], color="red")
                
                for i in range(len(v)):
                    flop_cutoff_pts = cutoff[v[i][0]][0] # WIP : generalise for possible more
                    flop_rays = cutoff[v[i][0]][1] # WIP : generalise for possible more
                    plt.scatter(*np.array(flop_cutoff_pts).T, s = 1)
                    plt.plot([0, 5*flop_rays[0, 0]], [0, 5*flop_rays[0, 1]], color="red")
                    plt.plot([0, 5*flop_rays[1, 0]], [0, 5*flop_rays[1, 1]], color="red")
                
                plt.title(f"Extended Kahler cone and GV cone for 2-face equiv. CY3s with h12 = {h12}")
                plt.savefig(f"figures/cone_moduli_cutoff/h_s={h_s}_h12={h12}")

def scatter_plot_hyperplane(birational, cutoff, h_s, perp_vecs = None, perp_coors = None, crit = 0.1):
    perp_vecs, perp_coors = np.array(perp_vecs), np.array(perp_coors)
    def filter_points(pts):
        perp_ip = np.einsum('na,ba->nb', pts, perp_vecs)
        mask = np.all((perp_ip > perp_coors) & (perp_ip < perp_coors + crit), axis = 1)
        return pts[mask]

    for h12, v0 in birational.items():
        not_found = True
        for k, v in v0.items():
            if v == []:
                continue
            elif not_found:
                plt.cla()
                plt.xlim((-5,5))
                plt.ylim((-5,5))
                cutoff_pts = cutoff[k][0]
                filtered_pts = filter_points(cutoff_pts)
                rays = cutoff[k][1]
                plt.scatter(*np.array(filtered_pts).T)
                
                for i in range(len(v)):
                    flop_cutoff_pts = cutoff[v[i][0]][0]
                    flop_filtered_pts = filter_points(flop_cutoff_pts)
                    flop_rays = cutoff[v[i][0]][1]
                    plt.scatter(*np.array(flop_filtered_pts).T)
                
                plt.title(f"Extended Kahler cone and GV cone for 2-face equiv. CY3s with h12 = {h12}\n 2-plane w. distance {perp_coors} along {perp_vecs}; tolerance {crit}")
                plt.savefig(f"figures/cone_moduli_cutoff/h_s={h_s}_h12={h12}")

                not_found = False
