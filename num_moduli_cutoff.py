import numpy as np
#Enable unfavourable CYs
#cytools.config.enable_experimental_features()
from cytools import fetch_polytopes, Polytope
from tqdm import tqdm

# |%%--%%| <Jn6Ef3is6D|hPiCUWr4nw>

# given a sample of moduli, evaluate the closest relatively coprime integer
def coprime_integer_scaling(moduli, tol = 1e-2, denom_max = 100):
    err = np.ones_like(moduli)
    coprime_int = np.ones_like(moduli)

    i = 2
    updated = np.zeros((moduli.shape[0])).astype(bool)

    while i <= denom_max and not np.all(updated):
        mi = moduli[~updated] * i
        ri = np.round(mi)
        err = np.abs(mi - ri)

        filter = np.all(err < tol, axis=1)
        if np.any(filter):
            rows = np.where(~updated)[0][filter]
            coprime_int[rows] = ri[filter]
            updated[rows] = True

        i += 1
    
    coprime_int[~updated] = np.round(moduli[~updated] * denom_max)
    nonzero = ~np.any(coprime_int == 0, axis = 1)
    
    assoc_moduli = moduli[nonzero, :]
    coprime_int = coprime_int[nonzero, :].astype(int)

    return assoc_moduli, coprime_int

def rays_multiple(rays, vector, rtol = 1e-8, atol = 1e-12):
    for ray in rays:
        if np.std(ray/vector) == 0:
            return True 
    return False

# given a sample of moduli, evalute the scaling where the instanton correction upto specified degree is less than specified cutoff
    # if parallel, gv invariants are stored in memory
def cutoff_ev_qvs_moduli(gvs, qvs, moduli, cutoff, max_trials, tol):
    N = moduli.shape[0]
    exponent = -2*np.pi*np.einsum("Nmd,Nd->Nm", qvs, moduli)

    def log_inst(l):
        weighted_exponent = l[:,None] * exponent
        return np.log(np.einsum("Nm,Nm->N", gvs, np.exp(weighted_exponent)))

    l_lo, l_hi = 1.0 * np.ones((N,)), 2.0 * np.ones((N,))
    iev_lo, iev_hi = log_inst(l_lo), log_inst(l_hi) # iev_lo > iev_hi !
    
    for _ in range(max_trials):
        l_hi_bool = (iev_hi < np.log(cutoff))
        if (~l_hi_bool).all():
            break
        else:
            l_hi[l_hi_bool] /= 2.0
            iev_hi = log_inst(l_hi)
    for _ in range(max_trials):
        l_lo_bool = (iev_lo > np.log(cutoff))
        if (~l_lo_bool).all():
            break
        else:
            l_lo[l_lo_bool] *= 2.0
            iev_lo = log_inst(l_lo)
    
    #print(f"lo : {moduli[i] * l_lo}, hi : {moduli[i] * l_hi}")
    for _ in range(max_trials):
        l_mid = 0.5 * (l_lo + l_hi)
        iev_mid = log_inst(l_mid)
        
        l_mid_bool = (np.abs(iev_mid - np.log(cutoff)) <= tol / iev_mid)
        l_mid_hi_bool = (iev_mid - np.log(cutoff) > tol / iev_mid)
        l_mid_lo_bool = (iev_mid - np.log(cutoff) < -tol/ iev_mid)

        if l_mid_bool.all():
            break
        else:
            l_lo[l_mid_bool] = l_hi[l_mid_bool] = l_mid[l_mid_bool]
            l_hi[l_mid_hi_bool] = l_mid[l_mid_hi_bool]
            l_lo[l_mid_lo_bool] = l_mid[l_mid_lo_bool]
        
    scaled_moduli = moduli * l_mid[:,None]

    return scaled_moduli

def cutoff_ev(moduli, p, min_points = int(2e1), cutoff = 1, tol = 1e-3, max_trials = 100,
                   precise_grading = True, parallel = True): 

    cy = p.triangulate().get_cy()
    rays = cy.toric_kahler_cone().extremal_rays()

    # update moduli; exclude points where one of rounded integers is 0
    moduli, coprime_ints = coprime_integer_scaling(moduli)
    scaled_moduli = np.zeros_like(moduli)

    N, d = moduli.shape[0], moduli.shape[1]

    qvs = np.zeros((N, min_points, d))
    gvs = np.zeros((N, min_points))

    gv_dict_default = cy.compute_gvs(min_points = round(min_points*1.2))
    
    for i in tqdm(range(N)):
        if rays_multiple(rays, coprime_ints[i]) or not cy.toric_kahler_cone().contains(coprime_ints[i]):
            # update moduli; exclude points where gv invariant is not defined
            continue
        
        if precise_grading:
            # this step is the most time consuming
            gv_dict = cy.compute_gvs(min_points = round(min_points*1.2), grading_vec = coprime_ints[i]).dok
        else:
            gv_dict = gv_dict_default

        gvs_i, qvs_i = np.array(list(gv_dict.values())[:min_points]), np.array(list(gv_dict.keys())[:min_points])
        
        if parallel:
            gvs[i] = gvs_i
            qvs[i] = qvs_i
        else: # compute using gv invariants then discard
            gvs_i, qvs_i, moduli_i = gvs_i.reshape(1, min_points), qvs_i.reshape(1, min_points, d), moduli[i].reshape(1, d)
            scaled_moduli[i] = cutoff_ev_qvs_moduli(gvs_i, qvs_i, moduli_i, cutoff, max_trials, tol)

    if parallel:
        scaled_moduli = cutoff_ev_qvs_moduli(gvs, qvs, moduli, cutoff, max_trials, tol)
        
    return scaled_moduli

#|%%--%%| <hPiCUWr4nw|tYtbezb5mg>

h_s = 2
#p = Polytope([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-6,-9]])
h_s_polytope = fetch_polytopes(h11 = h_s, lattice = "N", limit = 100)

#|%%--%%| <tYtbezb5mg|lyY5tUZStl>

import matplotlib.pyplot as plt
import num_index_density as idn
import pickle

if __name__ == "__main__":
    #rho = index_density(h_s)
    moduli_max = 5
    cutoff = 1
    total_moduli = int(1e1*10**h_s)
    
    # for i in range(len(h_s_polytope)):
    for i in range(0,1):
        p = h_s_polytope[i]
        rays = p.triangulate().get_cy().toric_kahler_cone().extremal_rays()
        
        cy_obj = idn.CalabiYau(p, moduli_max = moduli_max,
                               moduli_sample_factor = int(1), moduli_batch_no = total_moduli)

        cy = p.triangulate().get_cy()
        gvs = cy.compute_gvs(min_points = 400)
        
        moduli = cy_obj._moduli_projection_sample()
        n1_moduli = moduli / np.linalg.norm(moduli, axis = 1).reshape(-1,1)

        scaled_moduli = cutoff_ev(n1_moduli, p, cutoff = 1, max_trials=300,
                                  precise_grading = True, parallel = False)

        with open(f"data/num_moduli_cutoff/3_num_moduli_cutoff={cutoff}_mm={moduli_max}_tm={total_moduli*10}_hs={h_s}_ind={i}.json", "wb") as f:
            pickle.dump(scaled_moduli, f)

        plt.close()
        plt.scatter(n1_moduli[:,0], n1_moduli[:,1], s = 1, color = "green")
        plt.scatter(scaled_moduli[:,0], scaled_moduli[:,1], s = 1)
        plt.xlim((-5, 5))
        plt.ylim((-5, 5))
        plt.plot([0,5*rays[0,0]], [0,5*rays[0,1]], color = "red")
        plt.plot([0,5*rays[1,0]], [0,5*rays[1,1]], color = "red")

        plt.savefig(f"figures/3_num_moduli_cutoff={cutoff}_mm={moduli_max}_tm={total_moduli*10}_hs={h_s}_ind={i}")
