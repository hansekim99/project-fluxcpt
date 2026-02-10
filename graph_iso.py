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

@lru_cache(maxsize=None)
def canonical_weight(weight):
    d = len(weight)
    if d == 1:
        return weight 
    
    cols, rows = tuple(range(d)), tuple(range(d))
    best_mat, best_flat = None, None

    row_heur = sorted(rows, key=lambda r: (sum(weight[r]), tuple(sorted(weight[r])))) # heuristic; pruning

    def prefix_cols_for_rows(row_order):
        k = len(row_order)
        if k == 0:
            return cols
        col_keys = [tuple(weight[r][c] for r in row_order) for c in cols]
        return tuple(sorted(cols, key=lambda c: col_keys[c]))

    def flat_prefix(row_order, col_order):
        flat = []
        for r in row_order:
            flat.extend(weight[r][c] for c in col_order)
        return tuple(flat)

    def full_mat_for_rows(row_order):
        # fixed row order : lexicographically minimal column order wrt.
        # row-major flattening obtained by sorting columns by their full column vectors
        col_order = tuple(sorted(cols, key=lambda c: tuple(weight[r][c] for r in row_order)))
        mat = tuple(tuple(weight[r][c] for c in col_order) for r in row_order)
        flat = tuple(x for row in mat for x in row)
        return mat, flat

    def rec(row_order, remaining):
        nonlocal best_mat, best_flat

        k = len(row_order)
        if best_flat is not None and k > 0:
            col_order = prefix_cols_for_rows(row_order)
            pref = flat_prefix(row_order, col_order)
            if pref > best_flat[:k * d]:
                return

        if not remaining:
            mat, flat = full_mat_for_rows(row_order)
            if best_flat is None or flat < best_flat:
                best_mat, best_flat = mat, flat
            return

        for r in row_heur:
            if r in remaining:
                rec(row_order + (r,), tuple(x for x in remaining if x != r))

    rec((), rows)
    return best_mat

def coeff_list_canon(block_list : list, n : int, d : int, tqdm_hide = False):
    weight_coeff_dict = {}

    for block in tqdm.tqdm(block_list, bar_format=bar_fmt, desc="Canon.", disable=tqdm_hide):
        key = tuple(tuple(r) for r in block.cb_ph)
        canon = canonical_weight(key)
        coeff = block_coeff(n, d, block.int_points_num, block)
        weight_coeff_dict[canon] = weight_coeff_dict.get(canon, 0) + coeff

    return weight_coeff_dict

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

    for block in tqdm.tqdm(block_list, bar_format = bar_fmt, desc = "Graph Iso.", disable = tqdm_hide):
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


