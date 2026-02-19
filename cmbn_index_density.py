import itertools, tqdm, math
from typing import List, Optional
from dataclasses import dataclass

from importlib import reload
import cmbn_utils, cmbn_graph_iso

reload(cmbn_utils)
reload(cmbn_graph_iso)

@dataclass()
class Block:
    species: List[int]
    pre_coeff : object = 1
    cb_ph : Optional = None
    
    # total integers
    @property
    def int_points_num(self):
        return sum(cmbn_utils.int_pts_num(s) for s in self.species)
    
    @property
    def power_D(self):
        return sum(2*(1-cmbn_utils.int_pts_num(s)) for s in self.species)
    
    def insert(self, new_species):
        return Block(self.species + [new_species], self.pre_coeff)

    def incr_pre_coeff(self):
        return Block(self.species, self.pre_coeff + 1)

def prelist_from_prev_prelist(block_list : list[Block]):
    if len(block_list) == 0:
        return [Block([s]) for s in range(4)] 
    else:
        return [new_block for s in range(4) for block in block_list
            if (new_block := block.insert(s)).power_D >= -4]

def collect_permutations(block_list : list[Block]):
    # sort list by sorted species
    keyf = lambda b: tuple(sorted(b.species))
    block_list.sort(key=keyf)

    new_block_list: list[Block] = []
    for block in block_list:
        if not new_block_list or keyf(new_block_list[-1]) != keyf(block):
            new_block_list.append(block)
        else:
            temp_block = new_block_list[-1]
            new_block_list[-1] = temp_block.incr_pre_coeff()

    return new_block_list

def weight_coeffs(block_list : list[Block], n : int, d : int, tqdm_hide=False):
    out = {}

    for block in tqdm.tqdm(block_list, bar_format=cmbn_graph_iso.bar_fmt, desc=f"d={d}", disable=tqdm_hide):
        ext_conn, int_conn, excluded_pts = [], [], []
        int_pts_count = 0
        for i, s in enumerate(block.species):
            ext_conn += cmbn_utils.pre_ext_contraction_ext_points(s, 3*i, int_pts_count)
            int_conn += cmbn_utils.pre_ext_contraction_int_points(s, 3*i, int_pts_count)
            excluded_pts += cmbn_utils.excluded_pre_ext_contr_points(s, 3*i, int_pts_count)
            int_pts_count += cmbn_utils.int_pts_num(s)

        base = cmbn_graph_iso.block_coeff(n, d, int_pts_count, Block(block.species, 1))

        for perm_ext in itertools.permutations(range(d)):
            conn_1 = ext_conn[:]
            for i in range(len(ext_conn) // 3):
                conn_1[3*i+1] = excluded_pts[perm_ext[i]]

            coeff_ext = base * block.pre_coeff * cmbn_utils.one_ov_fact(perm_ext, d)

            for perm_int in itertools.permutations(range(int_pts_count)):
                conn_2 = []
                for conn in conn_1:
                    if conn[0] == 1:
                        conn_2.append(int_conn[perm_int[conn[1]]][1])
                    else:
                        conn_2.append(conn[1])

                w = [[0] * d for _ in range(d)]
                for i in range(d):
                    b = 3 * i
                    w[i][conn_2[b]     // 3] += 1
                    w[i][conn_2[b + 1] // 3] += 1
                    w[i][conn_2[b + 2] // 3] += 1

                key = tuple(tuple(row) for row in w)
                canon = cmbn_graph_iso.canonical_weight(key)
                out[canon] = out.get(canon, 0) + coeff_ext

    return out

#|%%--%%| <dtYypXKuDv|dq53fCWSkj>

def index_density(n): 
    weight_coeff_dict_total, _, _, block_list = cmbn_utils.init_prep(n)

    for d in range(1,n+1):
        block_list_precollection = prelist_from_prev_prelist(block_list)
        block_list_ = collect_permutations(block_list_precollection)

        weight_coeff_dict = weight_coeffs(block_list_, n, d, tqdm_hide=(d<=3))
        
        weight_coeff_dict_total.update(weight_coeff_dict)
        
        block_list = block_list_precollection

    return weight_coeff_dict_total

#|%%--%%| <dq53fCWSkj|MNbnT2PYKr>

import numpy as np 
import pickle, string

def partition(A):
    A = np.asarray(A)
    n = A.shape[0]
    seen_r, seen_c = np.zeros(n, bool), np.zeros(n, bool)
    out = []

    for r0 in range(n):
        if seen_r[r0]: 
            continue
        R, C = [], []
        stack = [('r', r0)]
        seen_r[r0] = True

        while stack:
            side, k = stack.pop()
            if side == 'r':
                R.append(k)
                for j in np.flatnonzero(A[k] > 0):
                    if not seen_c[j]:
                        seen_c[j] = True
                        stack.append(('c', j))
            else:
                C.append(k)
                for i in np.flatnonzero(A[:, k] > 0):
                    if not seen_r[i]:
                        seen_r[i] = True
                        stack.append(('r', i))
        out.append(A[np.ix_(sorted(R), sorted(C))])
    return out

def block_to_str(block):
    m = block.shape[0]
    assert block.ndim == 2 and block.shape[1] == m

    letters = iter(string.ascii_lowercase) #+ string.ascii_uppercase)
    L = [[] for _ in range(m)]
    R = [[] for _ in range(m)]

    for i in range(m):
        for j in range(m):
            k = block[i, j]
            labs = [next(letters) for _ in range(k)]
            L[i].extend(labs)
            R[j].extend(labs)

    def term(idx_list):
        core = ''.join(idx_list)
        return ('Z' + core)

    ins = [term(L[i]) for i in range(m)] + [term(R[j]) for j in range(m)]
    return ','.join(ins) + '->Z'

def h_s_to_invariant(h_s):
    with open(f"index_cmbn/index_cmbn_{h_s}.json", "rb") as f:
        rho = pickle.load(f)

    tot_str = dict()

    for k, v in rho.items():
        if v == 0:
            continue
        elif k == 0:
            tot_str[0] = v 
            continue
        else:
            invariant_blocks = partition(k)
            invariant_str = []
            for block in invariant_blocks:
                invariant_str.append(block_to_str(block))
            tot_str[tuple(invariant_str)] = v

    return tot_str

#|%%--%%| <MNbnT2PYKr|XhZdbuJo7u>

if __name__ == "__main__":
    h_s = 3
    #rho = index_density(h_s)
    
    #cmbn_utils.dict_print(rho)

    #with open(f"index_cmbn/index_cmbn_{h_s}.json", "wb") as f:
    #    pickle.dump(rho, f)
    
    tot_str = h_s_to_invariant(h_s)
    print(tot_str)
