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
def rays_multiple(rays, vector, rtol = 1e-8, atol = 1e-12):
    for ray in rays:
        if np.std(ray/vector) == 0:
            return True 
    return False

# given a sample of moduli, evalute the scaling where the instanton correction upto specified degree is less than specified cutoff
def cutoff_ev(moduli, p, min_points = int(2e1), cutoff = 1, tol = 1e-3, max_trials = 100):
    cy = p.triangulate().get_cy()
    rays = cy.toric_kahler_cone().extremal_rays()

    coprime_ints = coprime_integer_scaling(moduli)
    scaled_moduli = np.zeros_like(moduli)

    qv_gv_cache = {}
    
    for i in tqdm(range(coprime_ints.shape[0])):
        if rays_multiple(rays, coprime_ints[i]):
            continue 
        if not cy.toric_kahler_cone().contains(coprime_ints[i]):
            continue

        #gvs = cy.compute_gvs(min_points = min_points)
        gvs = cy.compute_gvs(min_points = min_points, grading_vec = coprime_ints[i])

        def inst(l):
            return instanton_ev(moduli[i] * l, gvs)

        l_lo, l_hi = 1.0, 2.0
        iev_lo, iev_hi = inst(l_lo), inst(l_hi) # iev_lo > iev_hi !

        for _ in range(max_trials):
            if iev_hi < cutoff + tol:
                break
            else:
                l_hi *= 2.0
                iev_hi = inst(l_hi)
        for _ in range(max_trials):
            if iev_lo > cutoff + tol:
                break
            else:
                l_lo /= 2.0
                iev_lo = inst(l_lo)
        
        #print(f"lo : {moduli[i] * l_lo}, hi : {moduli[i] * l_hi}")
        for _ in range(max_trials):
            l_mid = 0.5 * (l_lo + l_hi)
            iev_mid = inst(l_mid)

            if np.abs(iev_mid - cutoff) <= tol:
                l_lo = l_hi = l_mid
                break

            if iev_mid > cutoff:
                l_lo = l_mid
            else:
                l_hi = l_mid

            l_curr = l_hi
            #print(l_curr * moduli[i])
        
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
    
    #for i in range(len(h_s_polytope)):
    for i in range(13,14):
        p = h_s_polytope[i]
        rays = p.triangulate().get_cy().toric_kahler_cone().extremal_rays()
        
        cy_obj = idn.CalabiYau(p, moduli_max = moduli_max,
                               moduli_sample_factor = int(1), moduli_batch_no = total_moduli)

        cy = p.triangulate().get_cy()
        gvs = cy.compute_gvs(min_points = 400)
        print(gvs.dok.items())
        
        # moduli = cy_obj._moduli_projection_sample()
        # n1_moduli = moduli / np.linalg.norm(moduli, axis = 1).reshape(-1,1)
        #
        # scaled_moduli = cutoff_ev(n1_moduli, p, cutoff = 1, max_trials=300)
        # 
        # with open(f"data/num_moduli_cutoff/2_num_moduli_cutoff={cutoff}_mm={moduli_max}_tm={total_moduli*10}_hs={h_s}_ind={i}.json", "wb") as f:
        #     pickle.dump(scaled_moduli, f)
        #
        # plt.close()
        # plt.scatter(scaled_moduli[:,0], scaled_moduli[:,1], s = 1)
        # plt.xlim((-5, 5))
        # plt.ylim((-5, 5))
        # plt.plot([0,5*rays[0,0]], [0,5*rays[0,1]], color = "red")
        # plt.plot([0,5*rays[1,0]], [0,5*rays[1,1]], color = "red")
        #
        # plt.savefig(f"figures/2_num_moduli_cutoff={cutoff}_mm={moduli_max}_tm={total_moduli*10}_hs={h_s}_ind={i}")
