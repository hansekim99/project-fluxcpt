import time
from fractions import Fraction
from math import factorial
from contextlib import contextmanager

@contextmanager
def timer(label="", hide = False):
    t0 = time.time()
    yield
    if not hide:
        print(f"{label} : {time.time() - t0:.3f}s")

def int_pts_num(sp):
    if sp == 0:
        return 1
    elif sp == 1 or sp == 2:
        return 2
    elif sp == 3:
        return 3

def pre_ext_contraction_ext_points(sp, e_p_c, i_p_c):
    if sp == 0 or sp == 1:
        return [(0,e_p_c),None,(1,i_p_c)]
    elif sp == 2:
        return [(1,i_p_c),None,(1,i_p_c+1)]
    elif sp == 3:
        return [(1,i_p_c),None,(1,i_p_c+2)]

def pre_ext_contraction_int_points(sp, e_p_c, i_p_c):
    if sp == 0:
        return [(0,e_p_c+2)]
    elif sp == 1:
        return [(0,e_p_c+1),(0,e_p_c+2)]
    elif sp == 2:
        return [(0,e_p_c),(0,e_p_c+1)]
    elif sp == 3:
        return [(0,e_p_c),(0,e_p_c+1),(0,e_p_c+2)]

def excluded_pre_ext_contr_points(sp, e_p_c, i_p_c):
    if sp == 0:
        return [(0,e_p_c+1)]
    elif sp == 1 or sp == 3:
        return [(1,i_p_c+1)]
    elif sp == 2:
        return [(0,e_p_c+2)]

def one_ov_fact(p, n):
    seen = [0]*len(p)
    c = 0
    for i in range(len(p)):
        if not seen[i]:
            while not seen[i]:
                seen[i] = 1
                i = p[i]
            c += 1
    return Fraction((-1)**(c), factorial(n))

def block_print(block_list : list):
    print("--")
    for block in block_list:
        print(block)

def list_print(list_a : list, list_b : list):
    print("--")
    for i, block in enumerate(list_a):
        print(block, " : ", list_b[i])

def dict_print(dict_a):
    for k, v in dict_a.items():
        print(k, " : ", v)

def init_prep(n):
    weight_coeff_dict_total = {0 : 2 * factorial(n)}
    weight_list_total, coeff_list_total = [0], [2 * factorial(n)]
    block_list = []

    return weight_coeff_dict_total, weight_list_total, coeff_list_total, block_list

def dict_print_tex(dict_a):
    inv_str = ""
    fin_str = ""
    for k, v in dict_a.items():
        if k == 0:
            fin_str += str(v)
            s = 0
        else:
            if len(k) > s:
                s = len(k)
                i = 1
            if v == 0:
                continue
            elif v > 0:
                fin_str += "+"

            if v != 1 and v != -1:
                if v.denominator == 1:
                    fin_str += str(v)
                else:
                    if v.numerator > 0:
                        fin_str += f"\\tfrac{{{v.numerator}}}{{{v.denominator}}}"
                    else:
                        fin_str += f"-\\tfrac{{{-v.numerator}}}{{{v.denominator}}}"
            elif v == -1:
                fin_str += "-"
            inv_str += f"\mathcal{{F}}^{{({s})}}_{{{i}}}&={k}\\\\\n"
            fin_str += f"\mathcal{{F}}^{{({s})}}_{{{i}}}"
            i += 1
    print(fin_str)
    print(inv_str)

def sub_matr_col_row_sum(sub_matr, indices):
    col_tot = [0] * len(indices)
    row_tot = [0] * len(indices)
    for vi, i in enumerate(indices):
        for j in indices:
            col_tot[vi] += sub_matr[i][j]
            row_tot[vi] += sub_matr[j][i]
    return (all(x == 3 for x in col_tot) & all(x == 3 for x in row_tot))

def decomp_matr(matr):
    def restr_matr(matr, indices):
        return tuple(tuple(matr[j][k] for j in indices) for k in indices)
    res = []
    nonres = list(range(len(matr[0])))
    for i in range(len(matr[0])):
        if i in nonres:
            sub_matr = restr_matr(matr, nonres)
            temp_res = [nonres.pop(0)]
            while not sub_matr_col_row_sum(matr, temp_res):
                temp_res.append(nonres.pop(0))
            res.append(temp_res)
    
    fin_dict = {}
    for ind in res:
        key = restr_matr(matr, ind)
        fin_dict[key] = fin_dict.get(key, 0) + 1

    return fin_dict

#def dict_print_res(dict_a):
#    prev_degs, temp_degs = [], []
#    final_prods, final_vals = [], []
#    inv_str = ""
#    
#    for k, v in dict_a.items():
#        if k == 0:
#            s = 0
#            final_prods.append(1)
#            final_vals.append(v)
#        else:
#            if len(k) > s:
#                s = len(k)
#                i = 1
#                prev_degs.append(temp_degs)
#                temp_degs = []
#            
#            decomp = decomp_matr(k)

#            if len(next(iter(decomp.items()))[0]) != s:
#                for kn, vn in decomp.items():
#                    sind = len(kn)
#                    ind = prev_degs[sind].index(kn) + 1
#                    if vn != 1:
#                        if sind == 1 and ind == 1:
#                            inv_str += f"|\\mathcal{{F}}|^{2*vn}"
#                        elif sind == 2 and ind == 2:
#                            inv_str += f"(\\mathcal{{F}}_{{(2,1)}})^{vn}"
#                        else:
#                            inv_str += f"(\\mathcal{{F}}^{{({sind})}}_{{{ind}}})^{{{vn}}}"
#                    else:
#                        if sind == 1 and ind == 1:
#                            inv_str += r"|\mathcal{F}|^2"
#                        elif sind == 2 and ind == 2:
#                            inv_str += f"\\mathcal{{F}}_{{(2,1)}}"
#                        else:
#                            inv_str += f"\\mathcal{{F}}^{{({sind})}}_{{{ind}}}"
#            else:
#                inv_str += f"\\mathcal{{F}}^{{({s})}}_{{{i}}}"
#            inv_str += f"&={k}\\\\\n"

#            temp_degs.append(k)
#            i += 1

#            final_prods.append(decomp)
#            final_vals.append(v.numerator/v.denominator)

#    return inv_str, final_prods, final_vals

