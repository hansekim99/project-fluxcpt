import numpy as np
import math
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
    sample_batch_number: int = int(1e2)
    seed: int = 1
    
def kahler_cone_volume(h_s, hplane_n, sample_paras : MCSampleParas, sample_batch_no_custom = 0):
    rng = np.random.default_rng(seed = sample_paras.seed)

    sbno = max(10_000, sample_paras.sample_batch_number * 5, sample_batch_no_custom)

    sample = rng.normal(size = (sbno, h_s))
    sample /= np.linalg.norm(sample, axis=1).reshape(-1,1)

    sample = (sample @ hplane_n > 0).all(axis=1)
    volume = (2.0 * np.pi**(h_s/2.0) / math.gamma(h_s/2.0)) * sample.mean()

    return volume

def moduli_projection_sample(h_s, kahler_rays, sample_paras : MCSampleParas):
    rng = np.random.default_rng(seed = sample_paras.seed)
    rays_num = len(kahler_rays)

    raw_moduli_im_samples = rng.uniform(0, sample_paras.moduli_max, size = (sample_paras.sample_batch_number, rays_num))
    concat_krays = np.zeros((rays_num, rays_num))
    concat_krays[:, :h_s] = kahler_rays
    concat_krays[:, h_s:] = np.eye(rays_num, rays_num - h_s)
    
    moduli_im_samples_new = (raw_moduli_im_samples @ concat_krays)[:,:h_s]
    
    return moduli_im_samples_new

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

def distr_rho(cy_data : CYData, sample_paras : MCSampleParas):
    sno, sbno = sample_paras.sample_number, sample_paras.sample_batch_number
    moduli_distr = np.zeros((sno, cy_data.h_s))
    integrand_distr = np.zeros((sno,))
    
    ray_gv_list = cy_data.flop_facet_ray_gv_list

    sample_repeat_no = round(sno/sbno)

    invariant = h_s_to_invariant(cy_data.h_s)

    curr_pos = 0
    
    for _ in tqdm(range(sample_repeat_no)):
        # part 1 : sample points inside kahler cone
        ms_sample = moduli_projection_sample(cy_data.h_s, cy_data.kahler_extremal_rays, sample_paras)
            
        # part 2 : compute scalar density
        logdetG, specgeo, pdmask = kahler_num(cy_data.kijk, ms_sample, ray_gv_list)
        scalar = compute_invariant(invariant, specgeo) * np.pi**(-(cy_data.h_s+1)) # index vacua density = rho
        integrand = scalar * np.exp(logdetG)
        
        moduli_distr[curr_pos:curr_pos + pdmask.shape[0]] = ms_sample[pdmask]
        integrand_distr[curr_pos:curr_pos + pdmask.shape[0]] = integrand[pdmask]

        curr_pos += pdmask.shape[0]

    return moduli_distr[:curr_pos], integrand_distr[:curr_pos]

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

def integ_rho(cy_data : CYData, sample_paras : MCSampleParas):
    ray_gv_list = cy_data.flop_facet_ray_gv_list
    gv_dict_default = cy_data.cutoff_gv_dict

    moduli_distr, scalar_distr = distr_rho(cy_data, sample_paras)
    mask = moduli_check(ray_gv_list, gv_dict_default, moduli_distr)
    scalar_distr = scalar_distr[mask]

    scalar_mean = np.mean(scalar_distr)
    volume = kahler_cone_volume(cy_data.h_s, cy_data.hplane_n, sample_paras) * np.sum(mask) / scalar_distr.shape[0]

    scalar_integ = scalar_mean * volume

    return moduli_distr, scalar_distr, scalar_integ

#|%%--%%| <yMQdCgjdif|U9bfRiIKZf>

from cydata import load_cy_data_from_KS

if __name__ == "__main__":
    h_s = 3

    cy_data_gen = load_cy_data_from_KS(h_s)

    sample_paras = MCSampleParas(moduli_max = 10, sample_number = int(1e5), seed = None)

    for i, cy_data in enumerate(cy_data_gen):
        if i < 1:
            _, _, index_no = integ_rho(cy_data, sample_paras)
            
            print(i, " : ", index_no)
