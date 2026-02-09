import numpy as np
import math
from tqdm import tqdm

# |%%--%%| <YwStxGLPTl|e8noAh2Wcz>

class CalabiYau:
    def __init__(self, h21, arrayK, kahler_rays = np.array([1.]),
                moduli_interval = (2,5), qd3 = 10,
                moduli_sample_no = int(5e6)):
        self.n = h21
        self.arrayK = arrayK
        self.kahler_rays = kahler_rays
        
        self.mi = moduli_interval
        self.qd3 = qd3

        self.msno = int(1e1)                          # size of paralellised calculations
        self.mrno = int(moduli_sample_no / self.msno) # size of for loop

        self.rng = np.random.default_rng()
    
    def _moduli_uniform_sample(self):
        raw_moduli_im_samples = self.rng.uniform(*self.mi, size = (self.msno, self.n))
        
        moduli_im_samples = raw_moduli_im_samples @ self.kahler_rays.T

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
    
    def calc_index_vacua_density_cmbn(self, specgeo): # specgeo : single np array (rank 4 array; size : msno x h21 x h21 x h21)
        scalar_1_r = np.einsum('nabc,nabc->n', specgeo, np.conj(specgeo))
            
        scalar_30_03 = scalar_1_r ** 2
        scalar_21_12 = np.einsum('nabc,ndef,nabd,ncef->n', specgeo, specgeo, np.conj(specgeo), np.conj(specgeo))
        
        scalar_2_r = 1/2 * scalar_30_03 - 1/2 * scalar_21_12

        scalar_300_030_003 = scalar_1_r ** 3
        scalar_300_021_012 = scalar_1_r * scalar_21_12
        scalar_210_021_102 = np.einsum('nabc,ndef,nghl,ndbc,ngef,nahl->n', specgeo, specgeo, specgeo, np.conj(specgeo), np.conj(specgeo), np.conj(specgeo))
        scalar_210_111_012 = np.einsum('nabc,ndef,nghl,nabd,nceg,nfhl->n', specgeo, specgeo, specgeo, np.conj(specgeo), np.conj(specgeo), np.conj(specgeo))
        scalar_111_111_111 = np.einsum('nabc,ndef,nghl,ngbf,nael,ndhc->n', specgeo, specgeo, specgeo, np.conj(specgeo), np.conj(specgeo), np.conj(specgeo))
        
        scalar_3_r = - 1/6 * scalar_300_030_003 + 1/2 * scalar_300_021_012 - scalar_210_021_102 + scalar_210_111_012 - 1/3 * scalar_111_111_111

        #scalar = 2 - scalar_1_r
        #scalar = 4 - 2 * scalar_1_r + scalar_2_r
        scalar = 12 - 4 * scalar_1_r + 2 * scalar_2_r + scalar_3_r
        index_vacua_density = scalar * np.pi**(-(self.n+1))

        return index_vacua_density # single np array (rank 1 array; size : msno)
    
    def uniform_eval(self, mrl = False): # type : flux or cmbn
        moduli_distr = []
        weighted_ind_vac_den_distr = []

        for _ in tqdm(range(self.mrno), disable = not mrl):
            ms_uniform = self._moduli_uniform_sample() # change for new sampling
            logdetG, specgeo = self._ms_num(ms_uniform)
            
            index_vacua_density = self.calc_index_vacua_density_cmbn(specgeo)
            
            weighted_ind_vac_den = index_vacua_density * np.exp(logdetG)

            moduli_distr.append(ms_uniform)
            weighted_ind_vac_den_distr.append(weighted_ind_vac_den)
        
        moduli_distr = np.array(moduli_distr)
        weighted_ind_vac_den_distr = np.array(weighted_ind_vac_den_distr)

        return moduli_distr, weighted_ind_vac_den_distr

    def uniform_integrate(self, scalar_distr, mrl = False):
        averaged_scalar = np.mean(scalar_distr)
        
        volume_uniform = (self.mi[1] - self.mi[0]) ** self.n * np.linalg.det(self.kahler_rays) # change for new sampling

        prefactor = (2*np.pi) ** (2*(self.n+1)) / math.factorial(2*(self.n+1))
        axiodilaton = np.pi / 12

        return averaged_scalar * volume_uniform * prefactor * axiodilaton
