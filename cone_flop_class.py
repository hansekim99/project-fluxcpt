import numpy as np
import pickle

#|%%--%%| <VyiaPz74Fb|uQF5UqC7hp>

from collections import defaultdict

def identify_nilpotent(cy, min_points = 1e3, criter = 3):
    dual_rays = cy.mori_cone_cap(in_basis = True).extremal_rays()
    null_rays = []
    gvs = cy.compute_gvs(min_points = int(min_points)).dok
    # WIP : find flop facets
    for i in range(dual_rays.shape[0]):
        ray = dual_rays[i]
        if tuple(ray) in gvs:
            if not tuple(criter * ray) in gvs:
                null_rays.append([ray, gvs[(tuple(ray))]])
                # ray_gvs = []
                # for i in range(criter):
                #     if tuple((i+1) * ray) in gvs:
                #         print(ray, i)
                #         ray_gvs.append(gvs[tuple((i+1) * ray)])
                #     else:
                #         break
                # null_rays.append([ray, ray_gvs])

    return null_rays

def _diffeo_class_flops(polytopes):
    diffeo = defaultdict(dict)
    flops = defaultdict(dict)

    for i, p in enumerate(polytopes):
        # print(i, p.points())
        polytope_frst = p.ntfe_frsts()
        for j, t in enumerate(polytope_frst):
            if j != 0:
                print(i,j)
            cy = t.get_cy()
            wall = (tuple(cy.second_chern_class(in_basis = True).tolist()), tuple(sorted(cy.intersection_numbers(in_basis = True).items())))
            if wall not in diffeo:
                diffeo[cy.h12()][wall] = [[i, j], cy]
                flop_ray_gv0_list = identify_nilpotent(cy)
                flops[cy.h12()][wall] = [[i, j], flop_ray_gv0_list]

            else:
                print([i,j], " is duplicate of ", diffeo[wall][0])
                #duplicates[wall] = [[i,j], cy]

    return diffeo, flops

def diffeo_class_flops(polytopes, h_s, mode = "load"):
    if mode == "save" or mode == "run":
        diffeo, flops = _diffeo_class_flops(polytopes)

    if mode == "save":       
        with open(f"data/cone_diffeo_class/h_s={h_s}.json", "wb") as f:
            pickle.dump(diffeo, f)
        with open(f"data/cone_diffeo_class/h_s={h_s}_flops.json", "wb") as f:
            pickle.dump(flops, f)

    elif mode == "load":
        with open(f"data/cone_diffeo_class/h_s={h_s}.json", "rb") as f:
            diffeo = pickle.load(f)
        with open(f"data/cone_diffeo_class/h_s={h_s}_flops.json", "rb") as f:
            flops = pickle.load(f)

    return diffeo, flops

def _birational_class(diffeo, flops):
    birational = defaultdict(dict)

    for h12, diffs in diffeo.items():
        for k, v in diffs.items():
            h_s = len(k[0])
            flop_ray_gv0_list = flops[h12][k][1]

            equiv_class = []
            
            if len(flop_ray_gv0_list) > 0:
                print("before flop : ", k)

            for flop_ray_gv0 in flop_ray_gv0_list:
                flop_ray, flop_gv0 = flop_ray_gv0[0], flop_ray_gv0[1]
                flop_second_chern_class = np.array(k[0])
                flop_second_chern_class += 2 * flop_gv0 * flop_ray
            
                flop_intersection_numbers = dict(k[1])
                for a in range(h_s):
                    for b in range(h_s):
                        for c in range(h_s):
                            if (a,b,c) in flop_intersection_numbers:
                                flop_intersection_numbers[(a,b,c)] -= int(flop_gv0*flop_ray[a]*flop_ray[b]*flop_ray[c])
                            else:
                                flop_intersection_numbers[(a,b,c)] = int(-flop_gv0*flop_ray[a]*flop_ray[b]*flop_ray[c])
                            if flop_intersection_numbers[(a,b,c)] == 0:
                                del flop_intersection_numbers[(a,b,c)]
                
                flop_key = (tuple(flop_second_chern_class.tolist()), tuple(sorted(flop_intersection_numbers.items())))

                if flop_key in diffs:
                    equiv_class.append([flop_key, flop_ray])
                else:
                    # WIP : find integer transformations
                    cy = v[1] # cy = diffeo[h12][k][1]
                    # second_chern_class = list(k[0])
                    # intersection_numbers = dict(k[1])

            birational[h12][k] = equiv_class

    return birational

def birational_class(diffeo, flops, h_s, mode = "load"):
    if mode == "save" or mode == "run":
        birational = _birational_class(diffeo, flops)

    if mode == "save":
        with open(f"data/cone_diffeo_class/h_s={h_s}_birational.json", "wb") as f:
            pickle.dump(birational, f)

    elif mode == "load":
        with open(f"data/cone_diffeo_class/h_s={h_s}_birational.json", "rb") as f:
            birational = pickle.load(f)

    return birational
