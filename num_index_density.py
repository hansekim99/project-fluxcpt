import numpy as np
import vegas
from tqdm import tqdm

#|%%--%%| <TW544k8r37|6WU9KONeFg>

from cmbn_index_density import h_s_to_invariant

def compute_invariant(invariant, specgeo):
    msno = specgeo.shape[0]
    total_invariant = np.zeros(msno)
        
    for k, v in invariant.items():
        if k == 0:
            total_invariant += np.ones(msno) * np.float64(v)
        else:
            term_invariant = np.ones(msno) * np.float64(v)
            for term_str in k:
                order = (len(term_str) - 2)//10
                paras = []
                for _ in range(order):
                    paras.append(specgeo)
                    #paras.append(specgeo.real)
                for _ in range(order):
                    paras.append(np.conj(specgeo))
                    #paras.append(specgeo.real)
                term_invariant = (term_invariant * np.einsum(term_str, *paras)).real
            total_invariant += term_invariant

    return total_invariant

def polylog(s, z, co = 20):
    res = np.zeros_like(z)

    for k in range(1, co):
        res += z ** k / k ** s

    return res

# |%%--%%| <6WU9KONeFg|NixyaaHNsh>

from dataclasses import dataclass
from cydata import CYData

@dataclass
class MCSampleParas:
    moduli_max: float
    sample_number: int = int(1e4)

@dataclass
class MCIntegParas:
    instanton_corr_kijk_mode: bool = True
    instanton_corr_metr_mode: bool = True
    instanton_cutoff_mode : float = 1
    graph_data_mode: bool = False

def kahler_num(kijk, moduli_sample, ray_gv_list = None):
    polyK = 8/6 * np.einsum('abc,na,nb,nc->n', kijk, moduli_sample, moduli_sample, moduli_sample)
    scalarK = -np.log(polyK)

    polyK_a = 8/2 * np.einsum('abc,nb,nc->na', kijk, moduli_sample, moduli_sample)
    polyK_ab = 8 * np.einsum('abc,nc->nab', kijk, moduli_sample)

    matrG_num = 1/4 * (polyK_a[:,:,None] * polyK_a[:,None,:] - polyK_ab[:,:,:] * polyK[:,None,None]) / polyK[:,None,None] ** 2
    matrG_inv_num = np.linalg.inv(matrG_num)
    
    eigs = np.linalg.eigvalsh(matrG_inv_num)
    pos_def_mask = np.all(eigs > 0, axis = -1)

    vb_num = np.linalg.cholesky(matrG_inv_num[pos_def_mask])
    scalarK = scalarK[pos_def_mask]
    
    if ray_gv_list == None:
        specgeo_num = np.einsum('nai,nbj,nck,abc->nijk', vb_num, vb_num, vb_num, kijk) * np.exp(scalarK[:,None,None,None])
    else:
        k_ijk = np.stack([kijk] * moduli_sample.shape[0]).astype('float64')
        for ray_gv in ray_gv_list:
            ray, gv = np.array(ray_gv[0]), ray_gv[1]
            k_ijk = k_ijk + gv * ray[None,:,None,None] * ray[None,None,:,None] * ray[None,None,None,:] * polylog(0, np.exp(-2*np.pi*(moduli_sample @ ray)))[:,None,None,None]

        specgeo_num = np.einsum('nai,nbj,nck,nabc->nijk', vb_num, vb_num, vb_num, k_ijk) * np.exp(scalarK[:,None,None,None])
    
    _, logdetG_num = np.linalg.slogdet(matrG_num)
    logdetG_num = logdetG_num[pos_def_mask]

    return logdetG_num, specgeo_num, pos_def_mask

#|%%--%%| <NixyaaHNsh|yMQdCgjdif>

def _nilpotent_begone(gvs, ray_gv_list):
    for ray_gv in ray_gv_list:
        ray = ray_gv[0]
        d = 1
        while tuple(d * ray) in gvs:
            del gvs[tuple(d * ray)]
    return gvs

def moduli_check(ray_gv_list, gv_dict_default, moduli):
    gv_dict = _nilpotent_begone(gv_dict_default, ray_gv_list)

    gvs = np.array(list(gv_dict.values()), dtype=np.float64)
    qvs = np.array(list(gv_dict.keys()), dtype=np.float64)

    exponent = np.exp(-2*np.pi*np.einsum("md,Nd->Nm", qvs, moduli))
    instanton = np.einsum("m,Nm->N", np.abs(gvs), polylog(3,exponent))
    mask = instanton < 1

    return mask

def integ_rho(cy_data : CYData, sample_paras : MCSampleParas, config_paras : MCIntegParas):
    ray_gv_list = cy_data.flop_facet_ray_gv_list
    gv_dict_default = cy_data.cutoff_gv_dict

    invariant = h_s_to_invariant(cy_data.h_s)
    
    @vegas.batchintegrand
    def distr_rho(ms):
        rho_final = np.zeros(ms.shape[0])

        kc_mask = (ms @ cy_data.hplane_n > 0).all(axis = 1) # kc : inside kahler cone
        if not np.any(kc_mask):
            return rho_final
        ms_kc = ms[kc_mask]

        ic_mask = moduli_check(ray_gv_list, gv_dict_default, ms_kc) # ic : instanton convergence
        if not np.any(ic_mask):
            return rho_final
        ms_ic = ms_kc[ic_mask]

        logdetG, specgeo, pd_mask = kahler_num(cy_data.kijk, ms_ic, ray_gv_list) # pd : positive definitive
        if not np.any(pd_mask):
            return rho_final
        logdetG_pd, specgeo_pd = logdetG[pd_mask], specgeo[pd_mask]

        rho = compute_invariant(invariant, specgeo_pd) * np.exp(logdetG_pd) * np.pi**(-(cy_data.h_s + 1))

        rho_pd = np.zeros(len(ms_ic))
        rho_pd[pd_mask] = rho

        rho_ic = np.zeros(len(ms_kc))
        rho_ic[ic_mask] = rho_pd

        rho_final[kc_mask] = rho_ic 

        return rho_final

    domain = [[-sample_paras.moduli_max, sample_paras.moduli_max]] * cy_data.h_s
    integrand = vegas.Integrator(domain)
    
    # grid training
    integrand(distr_rho, nitn=10, neval=sample_paras.sample_number)
    
    # evaluation
    index_no = integrand(distr_rho, nitn=20, neval=sample_paras.sample_number * 2)

    return index_no.mean, index_no.sdev


#|%%--%%| <yMQdCgjdif|U9bfRiIKZf>

from cydata import load_cy_data_from_KS

if __name__ == "__main__":
    # h_s = 2 : sno ~ 1e2~3
    # h_s = 3 : sno ~ 1e4
    h_s = 3

    cy_data_gen = load_cy_data_from_KS(h_s)

    sample_paras = MCSampleParas(moduli_max = 20, 
                                 sample_number = int(1e4))

    config_paras = MCIntegParas(instanton_corr_kijk_mode = True,
                                instanton_corr_metr_mode = True,
                                instanton_cutoff_mode = 1,
                                graph_data_mode = False)

    for i, cy_data in enumerate(cy_data_gen):
        if i < 1:
            index_no, index_no_sdev = integ_rho(cy_data, sample_paras, config_paras)
        
            print(f"{i} : {index_no:.2E} , {index_no_sdev:.2E}")
            # TODO : save data to files
            # TODO : measure elapsed time
            # TODO : compare with old algo

            # TODO : include prefactor
            
            # TODO : implement faster tensor calculations
            # TODO : implement conformal + conifold integration
