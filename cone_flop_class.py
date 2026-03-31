import numpy as np
import pickle
from cytools import fetch_polytopes
from tqdm import tqdm

#|%%--%%| <VyiaPz74Fb|uQF5UqC7hp>

class CYData:
    def __init__(self, polytope, polytope_index, triangulation, cy,
                 cutoff_gv_dict_deg = 2, cutoff_gv_dict_facet_no = 10):
        # self.p = polytope
        # self.t = triangulation
        # self.cy = cy
        self.p_i = polytope_index # index : [i,j,k], ith polytope, jth triangulation out of k

        self.h_s, self.h_l = cy.h11(), cy.h12()
        self.scc = cy.second_chern_class(in_basis = True)
        self.kijk = cy.intersection_numbers(in_basis = True)

        self.mori_extremal_rays = cy.mori_cone_cap(in_basis = True).extremal_rays()
        self.kahler_extremal_rays = cy.mori_cone_cap(in_basis = True).dual_cone().extremal_rays()

        self.cutoff_gv_dict_deg = cy.compute_gvs(max_deg = cutoff_gv_dict_deg)

        self.nop_ray_gv_list = []
        self.diffeo_class_wall_data_list = []
        self.birational_class_wall_data_list = []
        self.birational_class_ray_gv_list = []

    def assign_nop_curves(self, nop_ray_gv_list):
        self.nop_ray_gv_list = nop_ray_gv_list

    def assign_diffeo_duplicate(self, dup_polytope, dup_polytope_index, dup_triangulation, dup_cy):
        wall_data = (dup_polytope, dup_polytope_index, dup_triangulation, dup_cy)
        self.diffeo_class_wall_data_list.append(wall_data)

    def assign_birational_duplicate(self, wall_data, ray_gv):
        self.birational_class_wall_data_list.append(wall_data)
        self.birational_class_ray_gv_list.append(ray_gv)

def identify_nilpotent(cy, max_deg = 3):
    gvs = cy.compute_gvs(max_deg = max_deg)

    nonzero_charges = gvs.dok # frak(C)
    grading_vector = gvs.grading_vec

    app_nilpotent_charges = [] # frak(N)
    
    for charge in nonzero_charges.keys():
        charge = np.array(charge)
        if np.gcd.reduce(charge) != 1: # check for coprime
            continue

        k = 1
        while tuple(k * charge) in nonzero_charges:
            k += 1
        
        # 2212.10573; equation A.1
        if (k * charge) @ grading_vector <= max_deg:
            cond = 0
            for j in range(1, k): # 1 to k-1
                cond += j ** 2 * nonzero_charges[tuple(j * charge)]
            if cond > 0:
                charge_gv = (charge, nonzero_charges[tuple(charge)])
                app_nilpotent_charges.append(charge_gv)

    nop_charges = app_nilpotent_charges # implement rest of algorithm

    return nop_charges

def wall_data_cy_obj(t):
    cy = t.get_cy()
    h_l = cy.h12()
    second_chern_class = tuple(cy.second_chern_class(in_basis = True).tolist())
    intersection_numbers = tuple(sorted(cy.intersection_numbers(in_basis = True).items()))
    wall_data = (h_l, second_chern_class, intersection_numbers)

    return wall_data, cy

def birational_flop_wall_data(wall_data, ray_gv):
    flop_ray, flop_gv = ray_gv[0], ray_gv[1]
    second_chern_class, intersection_numbers = np.array(wall_data[1]), dict(wall_data[2])
    h_s = len(second_chern_class)
    
    flop_second_chern_class = second_chern_class.copy() + 2 * flop_gv * flop_ray
    flop_intersection_numbers = intersection_numbers.copy()

    for a in range(h_s):
        for b in range(h_s):
            for c in range(h_s):
                if (a,b,c) != tuple(sorted((a,b,c))):
                    continue
                if (a,b,c) in flop_intersection_numbers:
                    flop_intersection_numbers[(a,b,c)] -= int(flop_gv*flop_ray[a]*flop_ray[b]*flop_ray[c])
                else:
                    flop_intersection_numbers[(a,b,c)] = int(-flop_gv*flop_ray[a]*flop_ray[b]*flop_ray[c])
                if flop_intersection_numbers[(a,b,c)] == 0:
                    del flop_intersection_numbers[(a,b,c)]
    
    flop_second_chern_class = tuple(flop_second_chern_class.tolist())
    flop_intersection_numbers = tuple(sorted(flop_intersection_numbers.items()))

    flop_wall_data = (wall_data[0], flop_second_chern_class, flop_intersection_numbers)

    return flop_wall_data

def _nontrivial_ntfe_frsts_diffeo_class(diffeo_class):
    final_dict = {}

    for k, v in diffeo_class.items():
        p_i = v.p_i
        if p_i[2] != 1:
            final_dict[k] = v

    return final_dict

def diffeo_class(h_s, polytope_no_max, flop_max_deg = 3, 
                 mode = "load", nontrivial_ntfe_frsts = False):
    if mode == "save" or mode == "run":
        polytopes = fetch_polytopes(h11 = h_s, lattice = "N", limit = polytope_no_max)
        diffeo = {}

        for i, p in tqdm(enumerate(polytopes)):
            polytope_frst = p.ntfe_frsts()

            for j, t in enumerate(polytope_frst):
                wall_data, cy = wall_data_cy_obj(t)

                if wall_data not in diffeo:
                    p_i = [i,j,len(polytope_frst)]
                    cy_data = CYData(p, p_i, t, cy)
                    cy_data.assign_nop_curves(identify_nilpotent(cy, flop_max_deg))
                    diffeo[wall_data] = cy_data

                else:
                    diffeo[wall_data].assign_diffeo_duplicate(p, p_i, t, cy)

        for wall_data, cydata in diffeo.items():
            ray_gv_list = cydata.nop_ray_gv_list
            for ray_gv in ray_gv_list:
                flop_wall_data = birational_flop_wall_data(wall_data, ray_gv)
                if flop_wall_data in diffeo: # toric found
                    cydata.assign_birational_duplicate(flop_wall_data, ray_gv)

    if mode == "save":       
        with open(f"data/cone_diffeo_class/h_s={h_s}_pnomax={polytope_no_max}.json", "wb") as f:
            pickle.dump(diffeo, f)

    elif mode == "load":
        with open(f"data/cone_diffeo_class/h_s={h_s}_pnomax={polytope_no_max}.json", "rb") as f:
            diffeo = pickle.load(f)

    if nontrivial_ntfe_frsts and (mode == "run" or mode == "load"):
        diffeo = _nontrivial_ntfe_frsts_diffeo_class(diffeo)

    return diffeo

