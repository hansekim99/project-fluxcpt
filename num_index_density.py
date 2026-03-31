import numpy as np
import math
import pickle
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

class CalabiYau:
    def __init__(self, cy_data, moduli_max = 10.0, 
                moduli_sample_factor = int(1e4), moduli_batch_no = int(1e2)):

        self.h_s = cy_data.h_s
        
        dictK = cy_data.kijk
        self.arrayK = np.array([[[dictK.get(tuple(sorted((i,j,k))), 0) for i in range(self.h_s)] for j in range(self.h_s)] for k in range(self.h_s)])

        # cone_hyperplane = cy.toric_mori_cone(in_basis = True).extremal_rays()
        cone_hyperplane = cy_data.mori_extremal_rays
        self.hplane_n = cone_hyperplane.T / (np.linalg.norm(cone_hyperplane, axis = 1))
        self.krays = cy_data.kahler_extremal_rays

        self.m_m = moduli_max

        self.msno = moduli_batch_no                          # size of paralellised calculations
        self.mrno = int(moduli_sample_factor / self.msno)    # size of for loop

        self.rng = np.random.default_rng()

    def _section_area(self, msno = None):
        if msno is None:
            msno = max(10_000, self.msno * 5)
        sample = self.rng.normal(size=(msno, self.h_s))
        sample /= np.linalg.norm(sample, axis=1).reshape(-1,1)
        sample = (sample @ self.hplane_n > 0).all(axis=1)
        return (2.0 * np.pi**(self.h_s/2.0) / math.gamma(self.h_s/2.0)) * sample.mean()

    def _moduli_projection_sample(self):
        rays_num = self.krays.shape[0]
        raw_moduli_im_samples = self.rng.uniform(0, self.m_m, size = (self.msno, rays_num))
        concat_krays = np.zeros((rays_num, rays_num))
        concat_krays[:, :self.h_s] = self.krays
        concat_krays[:, self.h_s:] = np.eye(rays_num, rays_num - self.h_s)
        
        moduli_im_samples_new = (raw_moduli_im_samples @ concat_krays)[:,:self.h_s]
        
        return moduli_im_samples_new

    def _ms_num(self, ms, gv_flop = None):
        polyK = 8/6 * np.einsum('abc,na,nb,nc->n', self.arrayK, ms, ms, ms)
        scalarK = -np.log(polyK)

        polyK_a = 8/2 * np.einsum('abc,nb,nc->na', self.arrayK, ms, ms)
        polyK_ab = 8 * np.einsum('abc,nc->nab', self.arrayK, ms)

        if gv_flop != None:
            gv, flop = np.array(gv_flop[1]), np.array(gv_flop[0])
            # multiple flops?
            
            z = np.exp(-2*np.pi*(ms @ flop))
            # polyK_inst = 1/(2*np.pi)**3 * gv * np.exp(-2*np.pi*(ms @ flop))
            polyK_inst = 1/(2*np.pi)**3 * gv * polylog(3, z)
            polyK += polyK_inst

            # polyK_a_inst = -1/(2*np.pi)**2 * gv * np.exp(-2*np.pi*(ms @ flop))[:, None] * flop[None, :]
            polyK_a_inst = - 1/(2*np.pi)**2 * gv * polylog(2, z)[:, None] * flop[None, :]
            polyK_a += polyK_a_inst

            # polyK_ab_inst = 1/(2*np.pi) * gv * np.exp(-2*np.pi*(ms @ flop))[:, None, None] * flop[None, :, None] * flop[None, None, :]
            polyK_ab_inst = 1/(2*np.pi) * gv * polylog(1, z)[:, None, None] * flop[None, :, None] * flop[None, None, :]
            polyK_ab += polyK_ab_inst
        
        matrG_num = 1/4 * (polyK_a[:,:,None] * polyK_a[:,None,:] - polyK_ab[:,:,:] * polyK[:,None,None]) / polyK[:,None,None] ** 2
        matrG_inv_num = np.linalg.inv(matrG_num)
        
        eigs = np.linalg.eigvalsh(matrG_inv_num)
        pos_def_mask = np.all(eigs > 0, axis = -1)

        vb_num = np.linalg.cholesky(matrG_inv_num[pos_def_mask])
        scalarK = scalarK[pos_def_mask]
        specgeo_num = np.einsum('nai,nbj,nck,abc->nijk', vb_num, vb_num, vb_num, self.arrayK) * np.exp(scalarK[:,None,None,None])
        
        _, logdetG_num = np.linalg.slogdet(matrG_num)
        logdetG_num = logdetG_num[pos_def_mask]

        return logdetG_num, specgeo_num, pos_def_mask

    def distr_rho(self, gv_flop = None):
        moduli_distr = []
        integrand_distr = []
        # npd_distr = []

        invariant = h_s_to_invariant(self.h_s)
        
        for _ in tqdm(range(self.mrno)):
            # part 1 : sample points inside kahler cone
            ms_sample = self._moduli_projection_sample()
                
            # part 2 : compute scalar density
            logdetG, specgeo, pdmask = self._ms_num(ms_sample, gv_flop)
            scalar = compute_invariant(invariant, specgeo) * np.pi**(-(self.h_s+1)) # index vacua density = rho
            integrand = scalar * np.exp(logdetG)
            
            moduli_distr.append(ms_sample[pdmask])
            integrand_distr.append(integrand)
            # npd_distr.append(ms_sample[np.logical_not(pdmask)])

        moduli_distr = np.concatenate(moduli_distr, axis=0)
        scalar_distr = np.concatenate(integrand_distr, axis = 0)
        # npd_distr = np.concatenate(npd_distr, axis = 0)

        return moduli_distr, scalar_distr

#|%%--%%| <NixyaaHNsh|qHQxDZINw6>

def _nilpotent_begone(gvs, ray_gv_list):
    for ray_gv in ray_gv_list:
        ray = ray_gv[0]
        d = 1
        while tuple(d * ray) in gvs:
            del gvs[tuple(d * ray)]
    return gvs

def moduli_check(cyd, moduli, gv_dict_mode = "degree"):
    ray_gv_list = cyd.nop_ray_gv_list

    if gv_dict_mode == "degree":
        gv_dict_default = cyd.cutoff_gv_dict_deg
    gv_dict = _nilpotent_begone(gv_dict_default.dok, ray_gv_list)

    gvs = np.array(list(gv_dict.values()), dtype=np.float64)
    qvs = np.array(list(gv_dict.keys()), dtype=np.float64)

    exponent = np.exp(-2*np.pi*np.einsum("md,Nd->Nm", qvs, moduli))
    instanton = np.einsum("m,Nm->N", gvs, exponent)
    mask = instanton < 1

    return mask

def index_density(h_s, diffeo, msf = int(1e3), m_m = 5,
                  mode = "load"):
    if mode == "save":
        index_density_distr = {}

        for wd, cyd in diffeo.items():
            cy_obj = CalabiYau(cyd, moduli_max = m_m, moduli_sample_factor = m_m ** h_s * msf)
            moduli_distr, scalar_distr = cy_obj.distr_rho()

            mask = moduli_check(cyd, moduli_distr)
            volume = cy_obj._section_area() * np.sum(mask) / scalar_distr.shape[0]
            scalar_distr = scalar_distr[mask]
            
            scalar_mean = np.mean(scalar_distr)
            scalar_integ = scalar_mean * volume

            index_density_distr[wd] = scalar_integ

        with open(f"data/num_index_density/h_s={h_s}_density.json", "wb") as f:
            pickle.dump(index_density_distr, f)

    elif mode == "load":
        with open(f"data/num_index_density/h_s={h_s}_density.json", "rb") as f:
            index_density_distr = pickle.load(f)

    return index_density_distr
