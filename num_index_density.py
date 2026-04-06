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

    def _ms_num(self, ms, ray_gv_list = None):
        polyK = 8/6 * np.einsum('abc,na,nb,nc->n', self.arrayK, ms, ms, ms)
        scalarK = -np.log(polyK)

        polyK_a = 8/2 * np.einsum('abc,nb,nc->na', self.arrayK, ms, ms)
        polyK_ab = 8 * np.einsum('abc,nc->nab', self.arrayK, ms)

        matrG_num = 1/4 * (polyK_a[:,:,None] * polyK_a[:,None,:] - polyK_ab[:,:,:] * polyK[:,None,None]) / polyK[:,None,None] ** 2
        matrG_inv_num = np.linalg.inv(matrG_num)
        
        eigs = np.linalg.eigvalsh(matrG_inv_num)
        pos_def_mask = np.all(eigs > 0, axis = -1)

        vb_num = np.linalg.cholesky(matrG_inv_num[pos_def_mask])
        scalarK = scalarK[pos_def_mask]
        
        if ray_gv_list == None:
            specgeo_num = np.einsum('nai,nbj,nck,abc->nijk', vb_num, vb_num, vb_num, self.arrayK) * np.exp(scalarK[:,None,None,None])
        else:
            k_ijk = np.stack([self.arrayK.copy()] * ms.shape[0]).astype('float64')
            for ray_gv in ray_gv_list:
                ray, gv = ray_gv[0], ray_gv[1]
                k_ijk = k_ijk + gv * ray[None,:,None,None] * ray[None,None,:,None] * ray[None,None,None,:] * polylog(0, np.exp(-2*np.pi*(ms @ ray)))[:,None,None,None]

            specgeo_num = np.einsum('nai,nbj,nck,nabc->nijk', vb_num, vb_num, vb_num, k_ijk) * np.exp(scalarK[:,None,None,None])
        
        _, logdetG_num = np.linalg.slogdet(matrG_num)
        logdetG_num = logdetG_num[pos_def_mask]

        return logdetG_num, specgeo_num, pos_def_mask

    def distr_rho(self, ray_gv_list = None):
        moduli_distr = []
        integrand_distr = []
        # npd_distr = []

        invariant = h_s_to_invariant(self.h_s)
        
        for _ in tqdm(range(self.mrno)):
            # part 1 : sample points inside kahler cone
            ms_sample = self._moduli_projection_sample()
                
            # part 2 : compute scalar density
            logdetG, specgeo, pdmask = self._ms_num(ms_sample, ray_gv_list)
            scalar = compute_invariant(invariant, specgeo) * np.pi**(-(self.h_s+1)) # index vacua density = rho
            integrand = scalar * np.exp(logdetG)
            
            moduli_distr.append(ms_sample[pdmask])
            integrand_distr.append(integrand)
            # npd_distr.append(ms_sample[np.logical_not(pdmask)])

        moduli_distr = np.concatenate(moduli_distr, axis=0)
        scalar_distr = np.concatenate(integrand_distr, axis = 0)
        # npd_distr = np.concatenate(npd_distr, axis = 0)

        return moduli_distr, scalar_distr

#|%%--%%| <NixyaaHNsh|0i0ujtyS8B>

def _nilpotent_begone(gvs, ray_gv_list):
    for ray_gv in ray_gv_list:
        ray = ray_gv[0]
        d = 1
        while tuple(d * ray) in gvs:
            del gvs[tuple(d * ray)]
    return gvs

def moduli_check(ray_gv_list, gv_dict_default, moduli):
    gv_dict = _nilpotent_begone(gv_dict_default.dok, ray_gv_list)

    gvs = np.array(list(gv_dict.values()), dtype=np.float64)
    qvs = np.array(list(gv_dict.keys()), dtype=np.float64)

    exponent = np.exp(-2*np.pi*np.einsum("md,Nd->Nm", qvs, moduli))
    instanton = np.einsum("m,Nm->N", np.abs(gvs), polylog(3,exponent))
    mask = instanton < 1

    return mask

def index_density(h_s, diffeo, msf = int(1e3), m_m = 5,
                  mode = "load", instanton_corr = True, cutoff_corr = True, gv_dict_mode = "degree"):
    if mode == "save":
        index_density_distr = {}

        for wd, cyd in diffeo.items():
            cy_obj = CalabiYau(cyd, moduli_max = m_m, moduli_sample_factor = m_m ** h_s * msf)
            ray_gv_list = cyd.nop_ray_gv_list

            if gv_dict_mode == "degree":
                gv_dict_default = cyd.cutoff_gv_dict_deg
            
            if instanton_corr:
                moduli_distr, scalar_distr = cy_obj.distr_rho(ray_gv_list)
            else:
                moduli_distr, scalar_distr = cy_obj.distr_rho()
            
            if cutoff_corr:
                mask = moduli_check(ray_gv_list, gv_dict_default, moduli_distr)
                masked_scalar_distr = scalar_distr[mask]

                scalar_mean = np.mean(masked_scalar_distr)
                volume = cy_obj._section_area() * np.sum(mask) / scalar_distr.shape[0]
                scalar_distr, moduli_distr = scalar_distr[mask], moduli_distr[mask]
            else:
                scalar_mean = np.mean(scalar_distr)
                volume = cy_obj._section_area()

            scalar_integ = scalar_mean * volume


            index_density_distr[wd] = (scalar_integ, scalar_distr, moduli_distr)

        with open(f"data/num_index_density/h_s={h_s}_density.json", "wb") as f:
            pickle.dump(index_density_distr, f)

    elif mode == "load":
        with open(f"data/num_index_density/h_s={h_s}_density.json", "rb") as f:
            index_density_distr = pickle.load(f)

    return index_density_distr

#|%%--%%| <0i0ujtyS8B|Hzevv0vjyw>

import matplotlib.pyplot as plt
from collections import Counter
from scipy.interpolate import griddata
from scipy.spatial import cKDTree

def contour_plot_2d(h_s, diffeo, cutoff, idd,
                    plot_size):
    for wd in diffeo.keys():
        plt.cla()
        plt.xlim((-plot_size, plot_size))
        plt.ylim((-plot_size, plot_size))
        
        cutoff_pts, cyd = cutoff[wd], diffeo[wd]
        scalar, moduli = idd[wd][1], idd[wd][2]
        # scalar_nc, moduli_nc = idd_non_corr[wd][1], idd_non_corr[wd][2]
        h_l = cyd.h_l

        kahler_rays = cyd.kahler_extremal_rays.tolist()
        birational_wall_data = cyd.birational_class_wall_data_list
        
        grid_x, grid_y = np.linspace(-plot_size,plot_size,100), np.linspace(-plot_size,plot_size,100)
        grid_x, grid_y = np.meshgrid(grid_x, grid_y)    
        
        for flop_wd in birational_wall_data:
            flop_cutoff_pts = cutoff[flop_wd]
            flop_scalar, flop_moduli = idd[flop_wd][1], idd[flop_wd][2]
            
            scalar = np.concatenate((scalar, flop_scalar))
            moduli = np.concatenate((moduli, flop_moduli))
            
            plt.scatter(*np.array(flop_cutoff_pts).T, s = 1)

            flop_cyd = diffeo[flop_wd]
            flop_kahler_rays = flop_cyd.kahler_extremal_rays

            for flop_ray in flop_kahler_rays:
                kahler_rays.append(flop_ray)
        
        plt.scatter(*np.array(cutoff_pts).T, s = 1)

        scalar_nz_mask = np.nonzero(scalar)
        scalar = scalar[scalar_nz_mask]
        moduli = moduli[scalar_nz_mask]
        log_abs_scalar = np.log(np.abs(scalar))
            
        grid_z = griddata((moduli[:,0], moduli[:,1]), log_abs_scalar, (grid_x, grid_y), method='cubic')
        # tree = cKDTree(np.c_[moduli[:,0], moduli[:,1]])
        # distances, _ = tree.query(np.c_[grid_x.ravel(), grid_x.ravel()])
        # distances = distances.reshape(grid_x.shape)
        # grid_z_m = np.ma.masked_where(distances > 0.5, grid_z)

        contour = plt.contour(grid_x, grid_y, grid_z, levels=20, cmap='viridis')
        plt.colorbar(contour, label = "log rho")

        counter_rays = Counter(tuple(ray) for ray in kahler_rays)

        for ray, count in counter_rays.items():
            if count > 1:
                plt.plot([0, 5*ray[0]], [0, 5*ray[1]], color="red", ls="dashed")
            else:
                plt.plot([0, 5*ray[0]], [0, 5*ray[1]], color="red")
        plt.scatter(*np.array(cutoff_pts).T, s = 1)
        
        plt.title(f"Extended Kahler cone and GV cone for 2-face equiv. CY3s with h12 = {h_l}")
        plt.gcf().set_size_inches(18.5, 10.5)
        plt.savefig(f"figures/num_index_density/h_s={h_s}_h12={h_l}_ctr")
        plt.close()

def scatter_plot_2d(h_s, diffeo, cutoff, idd,
                    plot_size):
    for wd in diffeo.keys():
        plt.cla()
        plt.xlim((-plot_size, plot_size))
        plt.ylim((-plot_size, plot_size))
        
        cutoff_pts, cyd = cutoff[wd], diffeo[wd]
        scalar, moduli = idd[wd][1], idd[wd][2]
        # scalar_nc, moduli_nc = idd_non_corr[wd][1], idd_non_corr[wd][2]
        h_l = cyd.h_l

        kahler_rays = cyd.kahler_extremal_rays.tolist()
        birational_wall_data = cyd.birational_class_wall_data_list
        
        for flop_wd in birational_wall_data:
            flop_cutoff_pts = cutoff[flop_wd]
            flop_scalar, flop_moduli = idd[flop_wd][1], idd[flop_wd][2]
            
            scalar = np.concatenate((scalar, flop_scalar))
            moduli = np.concatenate((moduli, flop_moduli))
            cutoff_pts = np.concatenate((cutoff_pts, flop_cutoff_pts))

            flop_cyd = diffeo[flop_wd]
            flop_kahler_rays = flop_cyd.kahler_extremal_rays

            for flop_ray in flop_kahler_rays:
                kahler_rays.append(flop_ray)
            
        scalar_nz_mask = np.nonzero(scalar)
        scalar = scalar[scalar_nz_mask]
        moduli = moduli[scalar_nz_mask]
        log_abs_scalar = np.log(np.abs(scalar))

        scatter = plt.scatter(moduli[:,0], moduli[:,1], c = log_abs_scalar, s = 1)
        plt.colorbar(scatter, label = "log rho")
        
        plt.scatter(*np.array(cutoff_pts).T, s = 1)

        counter_rays = Counter(tuple(ray) for ray in kahler_rays)

        for ray, count in counter_rays.items():
            if count > 1:
                plt.plot([0, 5*ray[0]], [0, 5*ray[1]], color="red", ls="dashed")
            else:
                plt.plot([0, 5*ray[0]], [0, 5*ray[1]], color="red")
        plt.scatter(*np.array(cutoff_pts).T, s = 1)
        
        plt.title(f"Extended Kahler cone and GV cone for 2-face equiv. CY3s with h12 = {h_l}")
        plt.gcf().set_size_inches(18.5, 10.5)
        plt.savefig(f"figures/num_index_density/h_s={h_s}_h12={h_l}_sct")
        plt.close()

