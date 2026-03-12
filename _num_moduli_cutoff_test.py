import numpy as np
from cytools import fetch_polytopes, Polytope

def multiples_subset(d, d_trunc, kmax=5):
    out = {}
    keys = set(d.keys())

    for (n, m) in d_trunc.keys():
        flag = True
        for k in range(2, kmax + 1):
            km = (k * n, k * m)
            if km not in keys:
                flag = False
        if flag:
            out[(n,m)] = d[(n,m)]
    return out

#|%%--%%| <OaDgb9jbfl|qFNo0QmB0U>

import pickle
import matplotlib.pyplot as plt

moduli_max = 5
cutoff = 1
total_moduli = int(1e3)

h_s = 2
#p = Polytope([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-6,-9]])
h_s_polytope = fetch_polytopes(h11 = h_s, lattice = "N", limit = 100)

#|%%--%%| <qFNo0QmB0U|7qJrRiUIyq>

def cutoff_check(moduli_sp, gvs_og):
    gvs_trunc = dict(list(gvs.dok.items())[:int(1e2)])

    gvs_subset = multiples_subset(gvs_og, gvs_trunc, kmax=5)

    exp_curves = np.zeros((moduli_sp.shape[0], len(gvs_subset), 5)).transpose(1,0,2)

    ind = 0
    for (n,m), v in gvs_subset.items():
        exp_curve = np.zeros((moduli_sp.shape[0], 5))
        for j in range(1,5):
            exp_curve[:,j] = gvs_og[(j*n,j*m)] * np.exp(-2*np.pi*(list((j*n,j*m)) @ moduli_sp.T))
        exp_curves[ind] = exp_curve
        ind += 1

    exp_curves = exp_curves.transpose(1,0,2)
    exp_curves[:,:,1:] = np.log(exp_curves[:,:,1:] + 1e-100)
    exp_curves[:,:,0] = 0

    return exp_curves

#|%%--%%| <7qJrRiUIyq|WjRoxEPkyl>

# check the exponential convergence of individual instanton corrections for points suggested to have converging instanton corrections

i = 0
with open(f"data/num_moduli_cutoff/num_moduli_cutoff={cutoff}_mm={moduli_max}_tm={total_moduli}_hs={h_s}/0_ind={i}.json", "rb") as f:
    moduli = pickle.load(f)

plt.close('all')

p = h_s_polytope[i]
rays = p.triangulate().get_cy().toric_kahler_cone().extremal_rays()
cy = p.triangulate().get_cy()
gvs = cy.compute_gvs(min_points = int(2e3))

#|%%--%%| <WjRoxEPkyl|aIc2c3x4PS>

ref_dir = [0.8,0.8]
print(str(ref_dir))
moduli_sp = np.array([0.6,1])[None,:] * np.linspace(0.6,1.2,5)[:,None]

for tar in range(moduli_sp.shape[0]):
    exp_curves = cutoff_check(moduli_sp[tar], gvs.dok)
    plt.close()
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(6, 9), constrained_layout=True)

    ax0.scatter(moduli[:, 0], moduli[:, 1], s=1)
    ax0.scatter(moduli_sp[tar, 0], moduli_sp[tar, 1], color="red")
    ax0.set_xlim((-5, 5))
    ax0.set_ylim((-5, 5))
    ax0.plot([0, 5*rays[0, 0]], [0, 5*rays[0, 1]], color="red")
    ax0.plot([0, 5*rays[1, 0]], [0, 5*rays[1, 1]], color="red")

    ax1.plot(exp_curves[0].T)
    
    plt.savefig(f"figures/check_moduli_cutoff={cutoff}_mm={moduli_max}_tm={total_moduli}_hs={h_s}/0_ind={i}_refdir={str(ref_dir).replace(".","-")}_tar={tar}")

