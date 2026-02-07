from collections import defaultdict
from functools import lru_cache
import igraph as ig
import tqdm

from math import factorial

bar_fmt = "{l_bar}{bar}| n={n:.2e} / total={total:.2e} [{remaining}]"

def block_coeff(n, d, m, block):
    if block.power_D == 0:
       fact_coeff = [1,-2*n-2*m,
                    n**2+n+(2*n+2)*m+m*(m-1)]
    elif block.power_D == -2:
        fact_coeff = [0,1,-n-m]
    elif block.power_D == -4:
        fact_coeff = [0,0,1]
    
    coeff = fact_coeff[0] * factorial(n - d + 1) + fact_coeff[1] * factorial(n - d)
    if n != d:
        coeff += fact_coeff[2] * factorial(n - d - 1)
        #singular_coeff = 0
    #else:
        #singular_coeff = fact_coeff[2]
    #return block.pre_coeff * coeff, block.pre_coeff * singular_coeff
    return block.pre_coeff * coeff


def coeff_list(block_list : list, n : int, d : int, tqdm_hide = False):
    @lru_cache(maxsize = None)
    def _weight_to_graph(weight):
        edges_perm = [(i, d+j) for i in range(d) for j in range(d) if weight[i][j] > 0]
        
        graph_perm = ig.Graph(n = 2*d, edges = edges_perm, directed = True)
        graph_perm.vs["color"] = [0]*d + [1]*d
        graph_perm.es["weight"] = [weight[i][j-d] for i, j in edges_perm]
        #graph_perm.simplify(loops=false, combine_edges=dict(weight="sum"))

        return graph_perm
    
    @lru_cache(maxsize = None)
    def _row_col_invariant(key):
        d_ = len(key)
        row_ms = tuple(sorted(tuple(sorted(r)) for r in key))
        col_ms = tuple(sorted(tuple(sorted(key[i][j] for i in range(d_))) for j in range(d_)))
        return (row_ms, col_ms)

    graphs_dict_ = defaultdict(list)
    weight_coeff_dict = dict()
    #weight_singular_coeff_dict = dict()

    for block in tqdm.tqdm(block_list, bar_format = bar_fmt, desc = "graph iso.", disable = tqdm_hide):
        key = tuple(tuple(r) for r in block.cb_ph)
        graph = _weight_to_graph(key)
        #coeff, singular_coeff = block_coeff(n, d, block.int_points_num, block)
        coeff = block_coeff(n, d, block.int_points_num, block)

        graph_key = _row_col_invariant(key)

        if not graphs_dict_[graph_key]:
            weight_coeff_dict[key] = coeff
            #weight_singular_coeff_dict[key] = singular_coeff
            graphs_dict_[graph_key].append(key)
            continue

        has_iso = False
        for i, key_ in enumerate(graphs_dict_[graph_key]):
            graph_ = _weight_to_graph(key_)
            is_iso = graph.isomorphic_vf2(graph_, 
                            color1 = graph.vs["color"], edge_color1=graph.es["weight"],
                            color2 = graph_.vs["color"], edge_color2=graph_.es["weight"])
            if is_iso:
                weight_coeff_dict[key_] += coeff
                #weight_singular_coeff_dict[key_] += singular_coeff
                has_iso = True
                break

        if not has_iso:
            weight_coeff_dict[key] = coeff
            #weight_singular_coeff_dict[key] = singular_coeff
            graphs_dict_[graph_key].append(key)

    #return weight_coeff_dict, weight_singular_coeff_dict
    return weight_coeff_dict

