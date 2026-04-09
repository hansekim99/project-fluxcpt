import compression
import json, hashlib
import itertools
import numpy as np
from cytools import fetch_polytopes

def identify_nilpotent(gvs, grading_vector, max_deg):
    # follows algorithm of 2212.10573
    nonzero_charges = gvs # frak(C)

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
                charge_gv = (charge.tolist(), nonzero_charges[tuple(charge)])
                app_nilpotent_charges.append(charge_gv)

    nop_charges = app_nilpotent_charges # TODO : implement rest of algorithm

    return nop_charges

def birational_flop_scc_kijk_dict(h_s, scc, kijk_dict, ray_gv):
    flop_ray, flop_gv = np.array(ray_gv[0]), ray_gv[1]
    
    flop_scc = scc.copy() + 2 * flop_gv * flop_ray
    flop_kijk_dict = kijk_dict.copy()

    for a in range(h_s):
        for b in range(h_s):
            for c in range(h_s):
                if (a,b,c) != tuple(sorted((a,b,c))):
                    continue
                if (a,b,c) in flop_kijk_dict:
                    flop_kijk_dict[(a,b,c)] -= int(flop_gv*flop_ray[a]*flop_ray[b]*flop_ray[c])
                else:
                    flop_kijk_dict[(a,b,c)] = int(-flop_gv*flop_ray[a]*flop_ray[b]*flop_ray[c])
                if flop_kijk_dict[(a,b,c)] == 0:
                    del flop_kijk_dict[(a,b,c)]

    return flop_scc, flop_kijk_dict

class CYData:
    def __init__(self, triangulation, cutoff_gv_dict_deg, flop_gv_dict_deg,
                 wall_data_set = None, polytope_index = None, _internal = False): # wall_data_set : dataset of wall data to identify birationally equiv. cy3s as toric flops
        cy = triangulation.get_cy() # TODO : filter out nonfavourable cys
        self.p_i = polytope_index # index : [i,j,k], ith polytope, jth triangulation out of k

        self.h_s = cy.h11() # h21, h11 of type iib theory respectively
        self.h_l = np.array(cy.h12())
        self.scc = np.array(cy.second_chern_class(in_basis = True)) # second chern class; np array

        kijk_dict = cy.intersection_numbers(in_basis = True) # intersection numbers; dict
        self.kijk = self._kijk_dict_to_kijk(kijk_dict) # intersection numbers; np array

        self.wall_data = self._compute_wall_data_hash(self.h_l, self.scc, self.kijk)

        if not _internal:
            self.mori_extremal_rays = cy.mori_cone_cap(in_basis = True).extremal_rays()
            self.kahler_extremal_rays = cy.mori_cone_cap(in_basis = True).dual_cone().extremal_rays()
            
            kahler_hyperplane = self.mori_extremal_rays
            self.hplane_n = kahler_hyperplane.T / (np.linalg.norm(kahler_hyperplane, axis = 1))

            self.cutoff_gv_dict = cy.compute_gvs(max_deg = cutoff_gv_dict_deg).dok # TODO : evaluate gvs along facets of mori cone
            flop_gv = cy.compute_gvs(max_deg = flop_gv_dict_deg)
            self.flop_gv_dict = flop_gv.dok
            self.flop_gv_grading_vec = flop_gv.grading_vec

            self.flop_facet_ray_gv_list = identify_nilpotent(self.flop_gv_dict, self.flop_gv_grading_vec, flop_gv_dict_deg)
            self.birational_equiv_ray_gv_list = []
            self.birational_equiv_wall_data_list = []

            if wall_data_set == None:
                wall_data_set = set()
                print("CYData object initialised without list of toric phases; birationally equivalent toric flops will not be calculated.")

            for ray_gv in self.flop_facet_ray_gv_list:
                flop_scc, flop_kijk_dict = birational_flop_scc_kijk_dict(self.h_s, self.scc, kijk_dict, ray_gv)
                flop_kijk = self._kijk_dict_to_kijk(flop_kijk_dict)
                flop_wall_data = self._compute_wall_data_hash(self.h_l, flop_scc, flop_kijk)
                if flop_wall_data in wall_data_set:
                    self.birational_equiv_ray_gv_list.append(ray_gv)
                    self.birational_equiv_wall_data_list.append(self.wall_data)

    @staticmethod
    def _compute_wall_data_hash(*arrays) -> str:
        hasher = hashlib.sha256()
        for arr in arrays:
            hasher.update(np.ascontiguousarray(arr).tobytes())
        return hasher.hexdigest()

    def _kijk_dict_to_kijk(self, kijk_dict):
        kijk = np.zeros((self.h_s, self.h_s, self.h_s), dtype=int) # intersection numbers; np array
        for ind, val in kijk_dict.items():
            for perm in itertools.permutations(ind):
                kijk[perm] = val
        return kijk

    @classmethod
    def wall_data_set_from_KS(cls, polytope_list, cutoff_gv_dict_deg, flop_gv_dict_deg):
        wall_data_set = set()
        for p in polytope_list:
            polytope_frst = p.ntfe_frsts()
            for t in polytope_frst:
                cy_data = cls(t, cutoff_gv_dict_deg, flop_gv_dict_deg, _internal = True)
                wall_data_set.add(cy_data.wall_data)
        return wall_data_set

    def save_to_hdf5(self, h5_group):
        h5_group.attrs['p_i'] = self.p_i if self.p_i is not None else [-1, -1, -1]
        h5_group.attrs['h_s'] = self.h_s
        h5_group.attrs['h_l'] = self.h_l
        
        h5_group.create_dataset('scc', data=np.array(self.scc), compression="gzip")
        h5_group.create_dataset('kijk', data=np.array(self.kijk), compression="gzip")
        h5_group.create_dataset('mori_extremal_rays', data=np.array(self.mori_extremal_rays), compression="gzip")
        h5_group.create_dataset('kahler_extremal_rays', data=np.array(self.kahler_extremal_rays), compression="gzip")
        h5_group.create_dataset('hplane_n', data=np.array(self.hplane_n), compression="gzip")
        
        def save_dok_dict(name, d):
            if not d: return
            h5_group.create_dataset(f'{name}_keys', data=np.array(list(d.keys())), compression="gzip")
            h5_group.create_dataset(f'{name}_vals', data=np.array(list(d.values())), compression="gzip")

        save_dok_dict('cutoff_gv_dict', self.cutoff_gv_dict)
        save_dok_dict('flop_gv_dict', self.flop_gv_dict)

        str_dtype = h5py.string_dtype(encoding='utf-8')
        h5_group.create_dataset('flop_facet_ray_gv_list', data=json.dumps(self.flop_facet_ray_gv_list), dtype=str_dtype)
        h5_group.create_dataset('birational_equiv_ray_gv_list', data=json.dumps(self.birational_equiv_ray_gv_list), dtype=str_dtype)
        h5_group.create_dataset('birational_equiv_wall_data_list', data=json.dumps(self.birational_equiv_wall_data_list), dtype=str_dtype)

    @classmethod
    def from_hdf5(cls, h5_group):
        obj = cls.__new__(cls)
        
        obj.p_i = list(h5_group.attrs['p_i'])
        obj.h_s = h5_group.attrs['h_s']
        obj.h_l = h5_group.attrs['h_l']
        
        obj.scc = tuple(h5_group['scc'][:])
        obj.kijk = tuple(h5_group['kijk'][:])
        obj.mori_extremal_rays = list(h5_group['mori_extremal_rays'][:])
        obj.kahler_extremal_rays = list(h5_group['kahler_extremal_rays'][:])
        obj.hplane_n = list(h5_group['hplane_n'][:])
        
        obj.wall_data = cls._compute_wall_data_hash(obj.h_l, obj.scc, obj.kijk)
        
        def load_dict(name):
            if f'{name}_keys' in h5_group and f'{name}_vals' in h5_group:
                keys = h5_group[f'{name}_keys'][:]
                vals = h5_group[f'{name}_vals'][:]
                return {tuple(k): v for k, v in zip(keys, vals)}
            return {}

        obj.cutoff_gv_dict = load_dict('cutoff_gv_dict')
        obj.flop_gv_dict = load_dict('flop_gv_dict')

        def load_json_dataset(name):
            if name in h5_group:
                val = h5_group[name][()]
                if isinstance(val, bytes):
                    val = val.decode('utf-8')
                return json.loads(val)
            return []

        obj.flop_facet_ray_gv_list = load_json_dataset('flop_facet_ray_gv_list')
        obj.birational_equiv_ray_gv_list = load_json_dataset('birational_equiv_ray_gv_list')
        
        raw_wall_data_list = load_json_dataset('birational_equiv_wall_data_list')
        obj.birational_equiv_wall_data_list = [tuple(w) for w in raw_wall_data_list]

        return obj

#|%%--%%| <l4k7YG8zLq|D1L36ySRP5>

import h5py

def save_cy_data_from_KS(h_s, polytope_no_start = 0, polytope_no_max = 100, cutoff_gv_dict_deg = 2, flop_gv_dict_deg = 3):
    with h5py.File(f"data/cydata/cy_data_h_s={h_s}.h5", "a") as db:
        polytope_list = fetch_polytopes(h11 = h_s, lattice = "N", limit = polytope_no_max)
        wall_data_set = CYData.wall_data_set_from_KS(polytope_list, cutoff_gv_dict_deg, flop_gv_dict_deg)

        for p_index, p in enumerate(polytope_list):
            if p_index < polytope_no_start:
                continue

            polytope_frst = p.ntfe_frsts()

            for t_index, t in enumerate(polytope_frst):
                t_max = len(polytope_frst)
                p_i = [p_index, t_index+1, t_max]

                cy_data = CYData(t, cutoff_gv_dict_deg, flop_gv_dict_deg, wall_data_set, p_i)
                wd_str = f"wd={cy_data.wall_data}"

                if wd_str in db:
                    del db[wd_str]
                grp = db.create_group(wd_str)
                cy_data.save_to_hdf5(grp)

def load_cy_data_from_KS(h_s):
    with h5py.File(f"data/cydata/cy_data_h_s={h_s}.h5", "r") as db:
        for wd_str in db.keys():
            cy_data = CYData.from_hdf5(db[wd_str])
            yield cy_data

#|%%--%%| <D1L36ySRP5|CS7uSa0LAP>

if __name__ == "__main__":
    h_s = 3

    # save_cy_data_from_KS(h_s, polytope_no_max = 100)

    cy_data_gen = load_cy_data_from_KS(h_s)

    for cy_data in cy_data_gen:
        print(cy_data.birational_equiv_ray_gv_list)
