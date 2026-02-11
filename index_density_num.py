import numpy as np
import math
from tqdm import tqdm

#|%%--%%| <TW544k8r37|6WU9KONeFg>

from index_density import h_s_to_invariant

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

# |%%--%%| <6WU9KONeFg|TW544k8r37>

class CalabiYau:
    def __init__(self, h21, arrayK, cone_hyperplane = np.array([[1.]]),
                moduli_max = 10, moduli_cutoff = 1, qd3 = 10,
                moduli_sample_no = int(5e6)):
        self.n = h21
        self.arrayK = arrayK
        self.hyperplane = cone_hyperplane

        self.m_m, self.m_c = moduli_max, moduli_cutoff
        
        self.qd3 = qd3

        self.msno = int(1e1)                          # size of paralellised calculations
        self.mrno = int(moduli_sample_no / self.msno) # size of for loop

        self.rng = np.random.default_rng()
    
    def _moduli_uniform_sample(self):
        # rejection sampling within the kahler cone
        # for each choice of n-1 rays, obtain one hyperplane

        moduli_im_samples = np.empty((self.msno, self.n))
        n_hplane = self.hyperplane.T / (np.linalg.norm(self.hyperplane, axis = 1))
        filled = 0
        
        while filled < self.msno:
            raw_moduli_im_samples = self.rng.uniform(-self.m_m, self.m_m, size = (self.msno, self.n))
            dist_from_hyperplane = raw_moduli_im_samples @ n_hplane
            filter = (dist_from_hyperplane > self.m_c).all(axis = 1)
            moduli_im_samples_new = raw_moduli_im_samples[filter]

            k = min(moduli_im_samples_new.shape[0], self.msno - filled)
            moduli_im_samples[filled:filled+k] = moduli_im_samples_new[:k]
            filled += k
            
            self.accepted_n += moduli_im_samples_new.shape[0]
            self.sampled_n += self.msno
        
        return moduli_im_samples

    def _ms_num(self, ms):
        polyK = 1/6 * np.einsum('abc,na,nb,nc->n', self.arrayK, ms, ms, ms)
        scalarK = -np.log(2*polyK)

        polyK_a = 1/2 * np.einsum('abc,nb,nc->na', self.arrayK, ms, ms)
        polyK_ab = np.einsum('abc,nc->nab', self.arrayK, ms)
        
        matrG_num = 1/4 * (polyK_a[:,:,None] * polyK_a[:,None,:] - polyK_ab[:,:,:] * polyK[:,None,None]) / polyK[:,None,None] ** 2
        matrG_inv_num = np.linalg.inv(matrG_num)

        _, logdetG_num = np.linalg.slogdet(matrG_num)
        
        vb_num = np.linalg.cholesky(matrG_inv_num)
        specgeo_num = -1j * np.einsum('nai,nbj,nck,abc->nijk', vb_num, vb_num, vb_num, self.arrayK) * np.exp(scalarK)[:,None,None,None]

        return logdetG_num, specgeo_num.imag * 1j
    
    def uniform_eval(self, mrl = False):
        moduli_distr = []
        weighted_ind_vac_den_distr = []
        self.accepted_n, self.sampled_n = 0, 0

        invariant = h_s_to_invariant(self.n)

        for _ in tqdm(range(self.mrno), disable = not mrl):
            ms_uniform = self._moduli_uniform_sample()
            logdetG, specgeo = self._ms_num(ms_uniform)
            
            scalar = compute_invariant(invariant, specgeo)
            index_vacua_density = scalar * np.pi**(-(self.n+1))
            
            weighted_ind_vac_den = index_vacua_density * np.exp(logdetG)

            moduli_distr.append(ms_uniform)
            weighted_ind_vac_den_distr.append(weighted_ind_vac_den)
        
        moduli_distr = np.array(moduli_distr)
        weighted_ind_vac_den_distr = np.array(weighted_ind_vac_den_distr)

        return moduli_distr, weighted_ind_vac_den_distr

    def uniform_integrate(self, scalar_distr):
        averaged_scalar = np.mean(scalar_distr)
        
        volume_uniform = (2*self.m_m)**self.n * self.accepted_n/self.sampled_n

        prefactor = (2*np.pi) ** (2*(self.n+1)) / math.factorial(2*(self.n+1))
        axiodilaton = np.pi / 12

        return averaged_scalar * prefactor * axiodilaton * volume_uniform
