import pickle
import numpy as np
import matplotlib.pyplot as plt
from cytools import fetch_polytopes, Polytope

#|%%--%%| <IRAmQwrIo5|GzXC3Hgzzg>

moduli_max = 5
cutoff = 1
total_moduli = int(1e3)

h_s = 2
#p = Polytope([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-6,-9]])
h_s_polytope = fetch_polytopes(h11 = h_s, lattice = "N", limit = 100)

# explicitly check infinity cone of each cym

# for i in range(2,3):
for i in range(len(h_s_polytope)):
    with open(f"data/num_moduli_cutoff/num_moduli_cutoff={cutoff}_mm={moduli_max}_tm={total_moduli}_hs={h_s}/3_ind={i}.json", "rb") as f:
        moduli = pickle.load(f)

    plt.close('all')

    p = h_s_polytope[i]
    dual_rays = p.triangulate().get_cy().mori_cone_cap(in_basis=True).extremal_rays()
    rays = p.triangulate().get_cy().mori_cone_cap(in_basis=True).hyperplanes()
    cy = p.triangulate().get_cy()
    gvs = cy.compute_gvs(min_points = int(1e3))

    cone_infty = np.array(list(gvs.dok.keys())).T

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(6, 9), constrained_layout=True)
    fig.suptitle(f"{i}th CY3 in KS database for h_11 = {h_s}")
    ax0.scatter(cone_infty[0], cone_infty[1])
    ax0.set_xlim((-10.5,10.5))
    ax0.set_ylim((-10.5,10.5))
    ax0.plot([0, 5*dual_rays[0, 0]], [0, 5*dual_rays[0, 1]], color="red")
    ax0.plot([0, 5*dual_rays[1, 0]], [0, 5*dual_rays[1, 1]], color="red")
    ax0.set_title("Mori cone and integer sites with nonvanishing GVs")

    ax1.scatter(moduli[:,0],moduli[:,1],s=1)
    ax1.set_xlim((-5, 5))
    ax1.set_ylim((-5, 5))
    ax1.plot([0, 5*rays[0, 0]], [0, 5*rays[0, 1]], color="red")
    ax1.plot([0, 5*rays[1, 0]], [0, 5*rays[1, 1]], color="red")
    ax1.set_title("Kahler cone and points where instanton correction = 1")

    plt.savefig(f"figures/check_flop_cutoff={cutoff}_mm={moduli_max}_tm={total_moduli}_hs={h_s}/3_ind={i}")
