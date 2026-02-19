import numpy as np
#Enable unfavourable CYs
#cytools.config.enable_experimental_features()
from cytools import fetch_polytopes, Polytope

# |%%--%%| <Jn6Ef3is6D|hPiCUWr4nw>

def partition(d, n):
    if d == 1:
        yield (n,)
        return
    for k in range(n + 1):
        for rest in partition(d - 1, n - k):
            yield (k,) + rest

# given a sample of moduli, evaluate the instanton correction at a fixed degree (given gv invariants and maximum degree) along 
def instanton_ev(moduli, kahler_rays, gv, deg):
    h_s = moduli.shape[1]
    final = np.zeros(moduli.shape[0])
    #print(gv)
    for nums in partition(h_s, deg): # not all partitions are considered?
        qv = kahler_rays.T @ nums
        try:
            nq = gv[tuple(qv)]
            #print(qv, nq)
        except:
            continue
        instanton = np.exp(-2 * np.pi * (moduli @ qv))
        final += instanton * nq
        #print(instanton[:10], nq, final[:10])

    return final

#|%%--%%| <hPiCUWr4nw|tYtbezb5mg>

h_s = 2
#p = Polytope([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-6,-9]])
h_s_polytope = fetch_polytopes(h11 = h_s, lattice = "N", limit = 100)

#|%%--%%| <tYtbezb5mg|V89gww5Oac>

import matplotlib.pyplot as plt
import num_index_density as idn

moduli_max = 5
degree_max = 30

radial_moduli_sample = 50

if __name__ == "__main__":
    #rho = index_density(h_s)
    p = h_s_polytope[7]
    cy = p.triangulate().get_cy()
    dictK = cy.intersection_numbers(in_basis = True)
    print(dictK)
    
    cone_hyperplane = cy.toric_mori_cone(in_basis = True).extremal_rays() # kahler rays
    print(cy.toric_kahler_cone().extremal_rays())
    arrayK = np.array([[[dictK.get(tuple(sorted((i,j,k))), 0) for i in range(h_s)] for j in range(h_s)] for k in range(h_s)])
    
    cy_obj = idn.CalabiYau(h_s, arrayK, cone_hyperplane,
                   moduli_max, moduli_cutoff = 0,
                   moduli_sample_no = int(1e4) * moduli_max ** h_s, moduli_batch_no = int(1e4))

    # 1 or 2 for cy sampling, 3 for sampling entire circle
    moduli_sample = cy_obj._moduli_hitandrun_sample()
    #moduli_sample = cy_obj._moduli_uniform_sample()
    #moduli_sample = np.column_stack((np.cos(2*np.pi*np.linspace(0,1,radial_moduli_sample)), np.sin(2*np.pi*np.linspace(0,1,radial_moduli_sample)))) * moduli_max

    # 1 for radial sampling near origin, 2 for uniform sampling
    #moduli_sample_long = np.tensordot(moduli_sample, 1/np.linspace(1,radial_moduli_sample,radial_moduli_sample),axes=0).transpose(0,2,1).reshape(-1,h_s)
    moduli_sample_long = np.tensordot(moduli_sample, np.linspace(1/radial_moduli_sample,1,radial_moduli_sample),axes=0).transpose(0,2,1).reshape(-1,h_s)

    gv = cy.compute_gv(max_deg = degree_max)

    instanton = instanton_ev(moduli_sample_long, kahler_rays = cone_hyperplane, 
                             gv = gv.dok, deg = degree_max)
    
    m = np.isfinite(np.log(np.abs(instanton)))

    # 1 for h_s = 2, 2 for h_s = 1
    x, y, z = moduli_sample_long[:,0][m], moduli_sample_long[:,1][m], np.log(np.abs(instanton))[m]
    #x, z = moduli_sample_long[:,0][m], np.log(instanton)[m]
    #y = np.ones_like(x)

    plt.figure()
    # 1 for contour map, 2 for colour map
    cs = plt.tricontour(x, y, z, levels=20)
    plt.clabel(cs, inline=True, fontsize=8)
    #cs = plt.scatter(x, y, c=z, s=8)
    plt.colorbar(cs, label=f"sum_q N_q exp(-2 pi q * (t_1, t_2))for q degree = {degree_max}")
    plt.xlabel("t_1 = Im(z_1)")
    plt.ylabel("t_2 = Im(z_2)")
    plt.title(dictK)
    plt.tight_layout()

    plt.xlim(-moduli_max,moduli_max)
    plt.ylim(-moduli_max,moduli_max)

    plt.show()
