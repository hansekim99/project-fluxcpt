import sympy as sym
import numpy as np
import math, datetime, os
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.spatial import ConvexHull
from scipy.stats import binned_statistic_2d

from scipy.special import erf, erfinv
from scipy.stats import qmc

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

    def _set_flux_mc_paras(self, quanta_sample_no = int(1e3), quanta_range = 5, stdev = np.sqrt(2/3)):
        self.qr = quanta_range

        self.qsno = int(1e3)                          # size of paralellised calculations
        self.qrno = int(quanta_sample_no / self.qsno) # size of for loop

        self.stdev = stdev
    
    def _moduli_uniform_sample(self):
        raw_moduli_im_samples = self.rng.uniform(*self.mi, size = (self.msno, self.n))
        
        moduli_im_samples = raw_moduli_im_samples @ self.kahler_rays.T

        return moduli_im_samples

    def _quanta_uniform_sample(self, sample_mode, u = None):
        if sample_mode == "uniform":
            zrs, zis = self.rng.uniform(-self.qr, self.qr, size = (self.qsno, self.n)), self.rng.uniform(-self.qr, self.qr, size = (self.qsno, self.n))
            xrs, xis = self.rng.uniform(-self.qr, self.qr, size = self.qsno), self.rng.uniform(-self.qr, self.qr, size = self.qsno)
        elif sample_mode == "exp":
            u = self.rng.random((self.qsno, 2*(self.n+1)))

            t = self.qr * self.stdev * erfinv((2.0 * u - 1.0) * erf(1/self.stdev))

            zrs, zis = t[:, :self.n], t[:, self.n:2 * self.n]
            xrs, xis = t[:, 2 * self.n], t[:, 2 * self.n + 1]
        
        zzs, xxs = zrs + 1j*zis, xrs + 1j*xis
        
        return zzs, xxs

    def _flux_volume(self, length, sampler):
        dim = 2*(self.n+1)
        if sampler == "exp":
            Z1 = self.qr * np.sqrt(np.pi) * self.stdev * erf(1/self.stdev)
            flux_volume = Z1 ** dim
        else:
            flux_volume = (2 * length) ** dim
        return flux_volume

    def _ms_num(self, ms):
        polyK = 1/6 * np.einsum('abc,na,nb,nc->n', self.arrayK, ms, ms, ms)
        scalarK = -np.log(2*polyK)

        polyK_a = 1/2 * np.einsum('abc,nb,nc->na', self.arrayK, ms, ms)
        polyK_ab = np.einsum('abc,nc->nab', self.arrayK, ms)
        #kappa_inv_ab = np.linalg.inv(polyK_ab)
        
        matrG_num = 1/4 * (polyK_a[:,:,None] * polyK_a[:,None,:] - polyK_ab[:,:,:] * polyK[:,None,None]) / polyK[:,None,None] ** 2
        #matrG_inv_num = (-np.exp(-scalarK)[:,None,None] * kappa_inv_ab / 2 + 2*np.einsum('na,nb->nab', ms, ms))
        matrG_inv_num = np.linalg.inv(matrG_num)

        # cy.compute_inverse_metric; polyK = volume
        #print(16 * matrG_num * polyK ** 2)

        _, logdetG_num = np.linalg.slogdet(matrG_num)
        
        vb_num = np.linalg.cholesky(matrG_inv_num)
        specgeo_num = -1j * np.einsum('nai,nbj,nck,abc->nijk', vb_num, vb_num, vb_num, self.arrayK) * np.exp(scalarK)[:,None,None,None]

        return logdetG_num, specgeo_num.imag * 1j

    def _flux_num(self, zzs, xxs, specgeo):
        zcs, xcs = np.conj(zzs), np.conj(xxs)
        zas_sr, xas_sr = zzs*zcs, xxs*xcs 

        I = np.eye(self.n)[None, :, :]

        matrB = np.einsum('nabc,mc->nmab', specgeo, zcs, optimize = True)
        matrC = np.einsum('nabc,mc->nmab', np.conj(specgeo), zzs, optimize = True)
        
        matrD = xxs[:,None,None]*I - zzs[:,:,None]*zcs[:,None,:]/xcs[:,None,None]
        
        # matrA = X*I - U*V^T; X scalar, U, V vector
        vecU, vecV = zcs, zzs / xxs[:,None]
        sclX, den = xcs, xcs - np.einsum('ma,ma->m', vecV, vecU)
        detA = (sclX ** (self.n-1)) * den

        # A^-1 B = 1/X * B + 1/(X * (X - V^T U)) * U * (V^T * B)
        VTB = np.einsum('ma,nmab->nmb', vecV, matrB, optimize = True)
        corr = np.einsum('ma,nmb->nmab', vecU, VTB, optimize = True)
        matrY = matrB / sclX[None, :, None, None] + corr / (sclX[None, :, None, None] * den[None, :, None, None])
        
        matrS = matrD - matrC @ matrY
        signS, logdetS = np.linalg.slogdet(matrS)
        integrand_log = - xas_sr - np.sum(zas_sr, axis=1) + np.log(xas_sr) + logdetS + np.log(np.abs(detA))
        
        integrand_sign = (signS * detA/np.abs(detA)).real

        exp_imp_pdist_log = (xas_sr + np.sum(zas_sr, axis=1))/((self.qr * self.stdev) ** 2)

        return integrand_log, integrand_sign, exp_imp_pdist_log
    
    def calc_index_vacua_density_cmbn(self, specgeo): # specgeo : single np array (rank 4 array; size : msno x h21 x h21 x h21)
        scalar_1_r = np.einsum('nabc,nabc->n', specgeo, np.conj(specgeo))
            
        scalar_30_03 = scalar_1_r ** 2
        scalar_21_12 = np.einsum('nabc,ndef,nabd,ncef->n', specgeo, specgeo, np.conj(specgeo), np.conj(specgeo))
        
        scalar_2_r = scalar_30_03 - scalar_21_12
        
        scalar = 4 - 2 * scalar_1_r + 1/2 * scalar_2_r
        #scalar = 2 - scalar_1_r
        index_vacua_density = scalar * np.pi**(-(self.n+1))

        return index_vacua_density # single np array (rank 1 array; size : msno)

    def calc_vacua_density_flux(self, specgeo, type, qrl = False, 
                                sampler = "exp"):
        index_vacua_density = np.zeros(specgeo.shape[0])
        
        dim, qrno = 2*(self.n+1), self.qrno
        pi_pref = np.pi**(-dim)

        flux_volume = self._flux_volume(self.qr, sampler)

        for _ in tqdm(range(qrno), disable = not qrl):
            zzs, xxs = self._quanta_uniform_sample(sample_mode = sampler)
            integrand_log, sign, imp_pdist  = self._flux_num(zzs, xxs, specgeo)
            if sampler == "exp":
                integrand_log += imp_pdist

            integrand = np.exp(integrand_log).real
            if type == "index":
                integrand *= sign
            
            averaged_scalar = np.mean(integrand, axis = 1)

            scalar = averaged_scalar * flux_volume * pi_pref

            index_vacua_density += scalar / qrno
        
        return index_vacua_density
    
    def uniform_eval(self, type, mrl = False): # type : flux or cmbn
        moduli_distr = []
        weighted_ind_vac_den_distr = []

        for _ in tqdm(range(self.mrno), disable = not mrl):
            ms_uniform = self._moduli_uniform_sample() # change for new sampling
            logdetG, specgeo = self._ms_num(ms_uniform)
            
            if type == "flux_index" or type == "flux_normal":
                if type == "flux_index":
                    index_vacua_density = self.calc_vacua_density_flux(specgeo, type = "index")
                elif type == "flux_normal":
                    index_vacua_density = self.calc_vacua_density_flux(specgeo, type = "normal")
            elif type == "cmbn":
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
