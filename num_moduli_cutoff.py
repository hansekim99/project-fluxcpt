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
    coprime_int = coprime_int[nonzero, :].astype(int)

    return coprime_int

# given a dict of gv coeffs, and a moduli point compute instanton correction
def instanton_ev(moduli, gv_d):
    final = 0
    for qv, gv in gv_d.dok.items():
        final += gv * np.exp(-2*np.pi*(qv @ moduli))
    return final

# 
def rays_multiple(rays, vector):
    for ray in rays:
        if np.std(ray/vector) == 0:
            return True 
    return False

# given a sample of moduli, evalute the scaling where the instanton correction upto specified degree is less than specified cutoff
def cutoff_ev(moduli, p, min_points = int(2e1), cutoff = 1, tol = 1e-2, max_trials = 100):
    cy = p.triangulate().get_cy()
    rays = cy.toric_kahler_cone().extremal_rays()

    coprime_ints = coprime_integer_scaling(moduli)
    scaled_moduli = np.zeros_like(moduli)
    
    for i in tqdm(range(coprime_ints.shape[0])):
        #print(cy.toric_kahler_cone().contains(coprime_ints[i]))
        if rays_multiple(rays, coprime_ints[i]):
            continue
        print(coprime_ints[i])
        gvs = cy.compute_gvs(min_points = min_points, grading_vec = coprime_ints[i])
        #l_i = 1/(2*np.pi) * np.log(sum(gvs.dok.values()))
        l_i = 1
        l_curr, iev = l_i, instanton_ev(moduli[i] * l_i, gvs)
        trials = l_i
        while np.abs(iev - cutoff) > tol and trials < max_trials:
            iev = instanton_ev(moduli[i] * l_curr, gvs)
            l_del = np.clip((iev-cutoff),-(1.1)**(trials-1), (1.1)**(trials-1)) * l_curr * (1.1)**(-trials)
            l_curr += l_del
            trials += 1
        
        scaled_moduli[i] = moduli[i] * l_curr
        
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
    total_moduli = int(1e3)
    
    for i in range(3,len(h_s_polytope)):
        p = h_s_polytope[i]
        rays = p.triangulate().get_cy().toric_kahler_cone().extremal_rays()
        
        cy_obj = idn.CalabiYau(p, moduli_max = moduli_max,
                               moduli_sample_factor = int(1), moduli_batch_no = total_moduli)

        moduli = cy_obj._moduli_projection_sample()

        scaled_moduli = cutoff_ev(moduli, p, cutoff = 1)
        
        with open(f"data/num_moduli_cutoff/num_moduli_cutoff={cutoff}_mm={moduli_max}_tm={total_moduli*10}_hs={h_s}_ind={i}.json", "wb") as f:
            pickle.dump(scaled_moduli, f)

        plt.close()
        plt.scatter(scaled_moduli[:,0], scaled_moduli[:,1], s = 1)
        plt.xlim((-5, 5))
        plt.ylim((-5, 5))
        plt.plot([0,5*rays[0,0]], [0,5*rays[0,1]], color = "red")
        plt.plot([0,5*rays[1,0]], [0,5*rays[1,1]], color = "red")

        plt.savefig(f"figures/num_moduli_cutoff={cutoff}_mm={moduli_max}_tm={total_moduli*10}_hs={h_s}_ind={i}")
