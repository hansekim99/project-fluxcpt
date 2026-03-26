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

# |%%--%%| <6WU9KONeFg|ytmUumXbRo>

class CalabiYau:
    def __init__(self, cy, moduli_max = 10.0, 
                moduli_sample_factor = int(1e4), moduli_batch_no = int(1e2)):

        self.h_s = cy.h11()
        self.kahler_cone = cy.toric_kahler_cone()
        
        dictK = cy.intersection_numbers(in_basis = True)
        self.arrayK = np.array([[[dictK.get(tuple(sorted((i,j,k))), 0) for i in range(self.h_s)] for j in range(self.h_s)] for k in range(self.h_s)])

        # cone_hyperplane = cy.toric_mori_cone(in_basis = True).extremal_rays()
        cone_hyperplane = cy.mori_cone_cap(in_basis = True).extremal_rays()
        self.hplane_n = cone_hyperplane.T / (np.linalg.norm(cone_hyperplane, axis = 1))
        # self.krays = cy.toric_kahler_cone().extremal_rays()
        self.krays = cy.mori_cone_cap(in_basis = True).hyperplanes()

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
    
    def _moduli_discrete_sample(self):
        moduli_im_samples = self.kahler_cone.find_lattice_points(min_points = self.msno * self.mrno)

        size = moduli_im_samples.max()

        return moduli_im_samples / size * self.m_m

    def _moduli_uniform_sample(self):
        # rejection sampling within the kahler cone
        # for each choice of n-1 rays, obtain one hyperplane

        moduli_im_samples = np.empty((self.msno, self.h_s))
        filled = 0
        
        while filled < self.msno:
            raw_moduli_im_samples = self.rng.uniform(-self.m_m, self.m_m, size = (self.msno, self.h_s)).astype(np.float128)
            dist_from_hyperplane = raw_moduli_im_samples @ self.hplane_n
            filter = (dist_from_hyperplane > 0).all(axis = 1)
            moduli_im_samples_new = raw_moduli_im_samples[filter]

            k = min(moduli_im_samples_new.shape[0], self.msno - filled)
            moduli_im_samples[filled:filled+k] = moduli_im_samples_new[:k]
            filled += k
        
        return moduli_im_samples

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
            
            z = gv * np.exp(-2*np.pi*(ms @ flop))
            # polyK_inst = 1/(2*np.pi)**3 * gv * np.exp(-2*np.pi*(ms @ flop))
            polyK_inst = 1/(2*np.pi)**3 * polylog(3, z)
            polyK += polyK_inst

            # polyK_a_inst = -1/(2*np.pi)**2 * gv * np.exp(-2*np.pi*(ms @ flop))[:, None] * flop[None, :]
            polyK_a_inst = -1/(2*np.pi)**2 * polylog(2, z)[:, None] * flop[None, :]
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

    def distr_rho(self, mrl = False, gv_flop = None):
        moduli_distr = []
        integrand_distr = []
        npd_distr = []

        invariant = h_s_to_invariant(self.h_s)
        
        for _ in tqdm(range(self.mrno), disable = not mrl):
            # part 1 : sample points inside kahler cone
            ms_sample = self._moduli_projection_sample()
                
            # part 2 : compute scalar density
            logdetG, specgeo, pdmask = self._ms_num(ms_sample, gv_flop)
            scalar = compute_invariant(invariant, specgeo) * np.pi**(-(self.h_s+1)) # index vacua density = rho
            integrand = scalar * np.exp(logdetG)
            
            moduli_distr.append(ms_sample[pdmask])
            integrand_distr.append(integrand)
            npd_distr.append(ms_sample[np.logical_not(pdmask)])

        moduli_distr = np.concatenate(moduli_distr, axis=0)
        scalar_distr = np.concatenate(integrand_distr, axis = 0)
        npd_distr = np.concatenate(npd_distr, axis = 0)

        return moduli_distr, scalar_distr, npd_distr

    
    def integ_rho(self, test_func = None, mrl = False, gv_flop = None, k = 20):
        moduli_distr = []
        integrand_distr = []

        invariant = h_s_to_invariant(self.h_s)

        section_area = self._section_area()

        for _ in tqdm(range(self.mrno), disable = not mrl):
            # part 1 : sample points inside kahler cone
            ms_sample = self._moduli_projection_sample()

            # rescale moduli points to spherical shell
            ms_norms = np.linalg.norm(ms_sample, axis=1).reshape(-1,1)
            ms_sample = ms_sample / ms_norms * self.m_m

            # part 2 : compute scalar density
            if test_func == None:
                logdetG, specgeo, posdef = self._ms_num(ms_sample)
                scalar = compute_invariant(invariant, specgeo) * np.pi**(-(self.h_s+1)) # index vacua density = rho
                integrand = scalar * np.exp(logdetG)
                p = -self.h_s

            elif test_func == "one":
                integrand = np.ones_like(ms_sample[:,0])
                p = self.h_s
            
            # part 3 : integrate; perform radial integral analytically
            # WIP : include cutoff as lower bound for analytic integral
            th = ms_sample / self.m_m
            r_b = self.m_m
            r_a = np.ones(th.shape[0])
            #r_a = np.asarray(cutoff(th), dtype=float)

            filter = np.isfinite(r_a) & (r_a < r_b)
            rf = np.zeros_like(r_a)
            rf[filter] = (r_b**p - r_a[filter]**p) / p
            integrand *= rf
            
            # part 4 : integrate; perform angular integral numerically; importance sampling via knn
            angular_dist = np.arccos(np.clip(th @ th.T,-1,1))
            np.fill_diagonal(angular_dist, np.inf)
            alpha = np.partition(angular_dist, k-1, axis=1)[:, k-1]
            impt = np.maximum(alpha, 1e-15) ** (self.h_s - 1)
            integrand = (integrand * impt) / np.mean(impt)

            moduli_distr.append(ms_sample)
            integrand_distr.append(integrand)
                    
        moduli_distr, scalar_distr = np.array(moduli_distr), np.array(integrand_distr)
        moduli_distr_lin = moduli_distr.reshape(-1,self.h_s)
        scalar_distr_lin = scalar_distr.reshape(-1)
            
        batch_means = scalar_distr.mean(axis=1) * section_area
        scalar_mean = batch_means.mean()
        scalar_se  = batch_means.std(ddof=1) / np.sqrt(batch_means.size)

        if test_func == None:
            prefactor = (2*np.pi) ** (2*(self.h_s+1)) / math.factorial(2*(self.h_s+1))
            axiodilaton = np.pi / 12

            const = prefactor * axiodilaton 
            scalar_mean *= const
            scalar_se *= const

        return moduli_distr_lin, scalar_distr_lin, scalar_mean, scalar_se

