import numpy as np
import igraph as ig
import itertools, tqdm
from typing import List, Optional
from dataclasses import dataclass

from utils import *
from graph_iso import *

@dataclass()
class Block:
    species: List[int]
    pre_coeff : object = 1
    cb_ph : Optional = None
    
    # total integers
    @property
    def int_points_num(self):
        return sum(int_pts_num(s) for s in self.species)
    
    @property
    def power_D(self):
        return sum(2*(1-int_pts_num(s)) for s in self.species)
    
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

def pre_contr_conn(block_list : list[Block]):
    new_block_list: list[Block] = []

    for block in block_list:
        int_pts_count = 0
        ext_conn, int_conn, excluded_pts = [], [], []

        for i, s in enumerate(block.species):
            ext_conn += pre_ext_contraction_ext_points(s, 3*i, int_pts_count)
            int_conn += pre_ext_contraction_int_points(s, 3*i, int_pts_count)
            excluded_pts += excluded_pre_ext_contr_points(s, 3*i, int_pts_count)
            int_pts_count += int_pts_num(s)

        new_block_list.append(Block(block.species, block.pre_coeff,
                                (ext_conn, int_conn, excluded_pts)))
    return new_block_list

def ext_contr_conn(block_list : list[Block]):
    new_block_list: list[Block] = []
    
    for block in block_list:
        ext_conn, int_conn, excluded_pts = block.cb_ph
        perms_ext = list(itertools.permutations(range(len(ext_conn) // 3)))
        contr_conn_list = []
        
        for perm in perms_ext:
            new_ext_conn = ext_conn[:]
            for i in range(len(ext_conn) // 3):
                new_ext_conn[3*i+1] = excluded_pts[perm[i]]

            one_ov_fact_ = one_ov_fact(perm, len(block.species))

            new_block_list.append(Block(block.species, block.pre_coeff * one_ov_fact_,
                                    (new_ext_conn, int_conn)))

    return new_block_list

def int_contr_conn(block_list : list[Block], tqdm_hide = False):
    new_block_list = []

    for block in tqdm.tqdm(block_list, bar_format = bar_fmt, desc = "Int. Contr.", disable = tqdm_hide):
        ext_conn, int_conn = block.cb_ph
        perms_int = list(itertools.permutations(range(len(int_conn))))

        for perm in perms_int:
            new_ext_conn = []
            for conn in ext_conn:
                if conn[0] == 1:
                    new_ext_conn.append(int_conn[perm[conn[1]]][1])
                else:
                    new_ext_conn.append(conn[1])
            new_block_list.append(Block(block.species, block.pre_coeff, new_ext_conn))
    
    return new_block_list

def contr_conn_to_weight(block_list: list[Block], tqdm_hide = False):
    new_block_list = []
    
    for block in tqdm.tqdm(block_list, bar_format = bar_fmt, desc = "Graph Conv.", disable = tqdm_hide):
        perm = block.cb_ph
        d = len(perm) // 3
        w = [[0]*d for _ in range(d)]

        for i in range(d):
            b = 3*i
            w[i][perm[b]   // 3] += 1
            w[i][perm[b+1] // 3] += 1
            w[i][perm[b+2] // 3] += 1

        new_block_list.append(Block(block.species, block.pre_coeff, w))

    return new_block_list

#|%%--%%| <dtYypXKuDv|dq53fCWSkj>

def index_density(n): # n = h21
    weight_coeff_dict_total, weight_list_total, coeff_list_total, block_list = init_prep(n)

    for d in range(1,n+1):
        #print("Order : ", d)
        with timer("Construct block list", hide = True):
            block_list_precollection = prelist_from_prev_prelist(block_list)
            #block_print(block_list_precollection)

        with timer("Collect permtuations", hide = True):
            block_list_ = collect_permutations(block_list_precollection)
            #block_print(block_list_)

        with timer("Pre-contraction", hide = True):
            pre_contr_conn_list = pre_contr_conn(block_list_)
            #block_print(pre_contr_conn_list)
            
        with timer("External indices contraction", hide = True):
            ext_contr_conn_list = ext_contr_conn(pre_contr_conn_list)
            #block_print(ext_contr_conn_list)

        with timer("Internal indices contraction", hide = True):
            int_contr_conn_list = int_contr_conn(ext_contr_conn_list, tqdm_hide = (d <= 4))
            #block_print(int_contr_conn_list)

        with timer("Graph construction", hide = True):
            weight_list = contr_conn_to_weight(int_contr_conn_list, tqdm_hide = (d <= 3))
            #block_print(weight_list)

        with timer("Graph isomorphism check", hide = True):
            #weight_coeff_dict, weight_singular_coeff_dict = coeff_list(weight_list, n, d, tqdm_hide = (d <= 3))
            weight_coeff_dict = coeff_list(weight_list, n, d, tqdm_hide = (d <= 3))

            #print(weight_coeff_dict)
        
        weight_coeff_dict_total.update(weight_coeff_dict)
        #weight_singular_coeff_dict_total.update(weight_singular_coeff_dict)
        
        block_list = block_list_precollection

    #print("Coefficients")
    #dict_print(weight_coeff_dict_total)
    #print("--")

    #print("Singular X^-2 terms")
    #dict_print(weight_singular_coeff_dict_total)
    return weight_coeff_dict_total

#|%%--%%| <dq53fCWSkj|42VGVABTt5>

if __name__ == "__main__":
    n = 5
    rho = index_density(n)
    
    dict_print(rho)
    
