import numpy as np
import cytools, inspect
#Enable unfavourable CYs
#cytools.config.enable_experimental_features()
from cytools import fetch_polytopes, Polytope

import matplotlib.pyplot as plt

# |%%--%%| <Jn6Ef3is6D|Rd7f0EjWS0>

h21 = 1

#p = Polytope([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-6,-9]])
h21_polytope = fetch_polytopes(h11 = h21, lattice = "N", limit = 100)

# |%%--%%| <Rd7f0EjWS0|vNxZpeMhoy>

p = h21_polytope[2]
#print(p)
cy = p.triangulate().get_cy()
dictK = cy.intersection_numbers(in_basis = True)

mori_cone = cy.toric_mori_cone(in_basis = True)
kahler_rays = np.linalg.inv(mori_cone.extremal_rays())
arrayK = np.array([[[dictK.get(tuple(sorted((i,j,k))), 0) for i in range(h21)] for j in range(h21)] for k in range(h21)])

# |%%--%%| <vNxZpeMhoy|kDkRakQzl4>

from importlib import reload
import index_density_num as idn
reload(idn)

# |%%--%%| <kDkRakQzl4|KPBh1spmQT>

cy_obj = idn.CalabiYau(h21, arrayK, kahler_rays,
                   moduli_interval = (1,3), qd3 = 50,
                   moduli_sample_no = int(5e4))
cy_obj._set_flux_mc_paras(quanta_sample_no = int(1e4), quanta_range = 5)

#|%%--%%| <KPBh1spmQT|D7PqYPM7Af>

_, specgeo = cy_obj._ms_num([[1]])

print(specgeo)

cy_obj.qsno = int(1e3)
zzs, xxs = cy_obj._quanta_uniform_sample(sample_mode = "exp")
integrand_log, sign, _ = cy_obj._flux_num(zzs, xxs, specgeo)
integrand = np.exp(integrand_log).real * sign

x, y, z = np.abs(zzs).T[0], np.abs(xxs), integrand[0]

print(z[:10])

#|%%--%%| <D7PqYPM7Af|guW8i5X1jE>

from scipy.interpolate import griddata
xi = np.linspace(x.min(), x.max(), 200)
yi = np.linspace(y.min(), y.max(), 200)
Xi, Yi = np.meshgrid(xi, yi)

# interpolate scattered data onto grid
Zi = griddata((x, y), z, (Xi, Yi), method="linear")

# contour plot
plt.contour(Xi, Yi, Zi, levels=20)
plt.xlabel("x")
plt.ylabel("y")
plt.colorbar(label="z")
plt.title(f'{specgeo[0,0,0,0]}')
plt.show()

# |%%--%%| <guW8i5X1jE|PZaS237ODP>

_, udv_cmbn        = cy_obj.uniform_eval(type = "cmbn",        mrl = True)
_, udv_flux_index  = cy_obj.uniform_eval(type = "flux_index",  mrl = True)
_, udv_flux_normal = cy_obj.uniform_eval(type = "flux_normal", mrl = True)

#|%%--%%| <PZaS237ODP|EJXmZfsibM>

usv_cmbn = cy_obj.uniform_integrate(udv_cmbn)
usv_flux_index = cy_obj.uniform_integrate(udv_flux_index)
usv_flux_normal = cy_obj.uniform_integrate(udv_flux_normal)

#|%%--%%| <EJXmZfsibM|x6N9eL8Ntt>

print(usv_cmbn)
print(usv_flux_index)
print(usv_flux_normal)
