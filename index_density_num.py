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
                moduli_max = 10, moduli_cutoff = 1,
                moduli_sample_no = int(5e6), moduli_batch_no = int(1e2)):
        self.n = h21
        self.arrayK = arrayK
        self.hplane_n = cone_hyperplane.T / (np.linalg.norm(cone_hyperplane, axis = 1))

        self.m_m, self.m_c = moduli_max, moduli_cutoff

        self.msno = moduli_batch_no                          # size of paralellised calculations
        self.mrno = int(moduli_sample_no / self.msno) # size of for loop

        self.rng = np.random.default_rng()

    def volume_uniform(self):
        raw = self.rng.uniform(-self.m_m, self.m_m, size=(self.msno*self.mrno, self.n))
        ok = (raw @ self.hplane_n > self.m_c).all(axis=1)
        return (2*self.m_m)**self.n * ok.mean()
    
    def _hr_seed_point(self):
        while True:
            x = self.rng.uniform(-self.m_m, self.m_m, size=self.n)
            if (x @ self.hplane_n > self.m_c).all():
                return x

    def _moduli_hitandrun_sample(self, thin=1, eps=1e-14):
        M, c, H = float(self.m_m), float(self.m_c), self.hplane_n  # H: (n, L)
        x = getattr(self, "_hr_x", None)
        if x is None:
            x = self._hr_seed_point()

        out = np.empty((self.msno, self.n), dtype=np.float64)

        for k in range(self.msno):
            for _ in range(thin):
                while True:
                    d = self.rng.normal(size=self.n)
                    d /= np.linalg.norm(d)

                    # box chord: -M <= x + t d <= M
                    lo = np.where(d > 0, (-M - x) / d, np.where(d < 0, ( M - x) / d, -np.inf))
                    hi = np.where(d > 0, ( M - x) / d, np.where(d < 0, (-M - x) / d,  np.inf))
                    t_lo, t_hi = lo.max(), hi.min()

                    # halfspaces: (x + t d)·h_j > c
                    a = d @ H
                    b = x @ H
                    z = np.abs(a) < 1e-15
                    if np.any(z & (b <= c + eps)):
                        continue

                    pos = a >  1e-15
                    neg = a < -1e-15
                    if np.any(pos):
                        t_lo = max(t_lo, np.max((c + eps - b[pos]) / a[pos]))
                    if np.any(neg):
                        t_hi = min(t_hi, np.min((c + eps - b[neg]) / a[neg]))

                    if t_hi > t_lo:
                        x = x + self.rng.uniform(t_lo, t_hi) * d
                        break
            out[k] = x
        self._hr_x = x
        return out

    def _moduli_uniform_sample(self):
        # rejection sampling within the kahler cone
        # for each choice of n-1 rays, obtain one hyperplane

        moduli_im_samples = np.empty((self.msno, self.n))
        filled = 0
        
        while filled < self.msno:
            raw_moduli_im_samples = self.rng.uniform(-self.m_m, self.m_m, size = (self.msno, self.n))
            dist_from_hyperplane = raw_moduli_im_samples @ self.hplane_n
            filter = (dist_from_hyperplane > self.m_c).all(axis = 1)
            moduli_im_samples_new = raw_moduli_im_samples[filter]

            k = min(moduli_im_samples_new.shape[0], self.msno - filled)
            moduli_im_samples[filled:filled+k] = moduli_im_samples_new[:k]
            filled += k
        
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
        integrand_distr = []
        self.accepted_n, self.sampled_n = 0, 0

        invariant = h_s_to_invariant(self.n)

        for _ in tqdm(range(self.mrno), disable = not mrl):
            #ms_uniform = self._moduli_uniform_sample()
            ms_uniform = self._moduli_hitandrun_sample()

            logdetG, specgeo = self._ms_num(ms_uniform)
            
            scalar = compute_invariant(invariant, specgeo)
            index_vacua_density = scalar * np.pi**(-(self.n+1))
            
            integrand_surface = index_vacua_density * np.exp(logdetG)
            
            #dist_at_cutoff = self.m_c * np.linalg.norm(ms_uniform,axis=1) / (np.min(ms_uniform @ self.hplane_n, axis = 1))
            dist_at_cutoff = self.m_c * np.linalg.norm(ms_uniform,axis=1) / (np.min(ms_uniform @ self.hplane_n, axis = 1))
            r_b, r_a = np.linalg.norm(ms_uniform,axis=1), dist_at_cutoff
            r_max = self.m_m * r_b / np.max(np.abs(ms_uniform), axis=1)

            g_u = integrand_surface * (r_b ** (2 * self.n))
            #integrand = integrand_surface * 1/(-2*self.n+1) * (1-(r_a/r_b)**(-2*self.n+1))
            integrand = g_u * (r_a**(-self.n) - r_max**(-self.n)) / (r_max**(self.n) - r_a**(self.n))

            #integrand = integrand_surface

            moduli_distr.append(ms_uniform)
            integrand_distr.append(integrand)
                    
        return np.array(moduli_distr), np.array(integrand_distr)

    def uniform_integrate(self, scalar_distr):
        batch_means = scalar_distr.mean(axis=1)
        scalar_mean = batch_means.mean()
        scalar_se  = batch_means.std(ddof=1) / np.sqrt(batch_means.size)

        prefactor = (2*np.pi) ** (2*(self.n+1)) / math.factorial(2*(self.n+1))
        axiodilaton = np.pi / 12

        const = prefactor * axiodilaton * self.volume_uniform()

        integ = const * scalar_mean
        integ_se = const * scalar_se

        return integ, integ_se
