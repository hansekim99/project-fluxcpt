import numpy as np
import cytools, inspect
#Enable unfavourable CYs
#cytools.config.enable_experimental_features()
from cytools import fetch_polytopes, Polytope

import matplotlib.pyplot as plt

# |%%--%%| <Jn6Ef3is6D|Rd7f0EjWS0>

h21 = 2
#p = Polytope([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-6,-9]])
h21_polytope = fetch_polytopes(h11 = h21, lattice = "N", limit = 100)

# |%%--%%| <Rd7f0EjWS0|vNxZpeMhoy>

p = h21_polytope[2]
#print(p)
cy = p.triangulate().get_cy()
dictK = cy.intersection_numbers(in_basis = True)
print(dictK)

mori_cone = cy.toric_mori_cone(in_basis = True)
kahler_rays = np.linalg.inv(mori_cone.extremal_rays())
arrayK = np.array([[[dictK.get(tuple(sorted((i,j,k))), 0) for i in range(h21)] for j in range(h21)] for k in range(h21)])

# |%%--%%| <vNxZpeMhoy|kDkRakQzl4>

from importlib import reload
import index_density_num as idn
reload(idn)

# |%%--%%| <kDkRakQzl4|KPBh1spmQT>

qr = 10
cy_obj = idn.CalabiYau(h21, arrayK, kahler_rays,
                   moduli_interval = (2,2.1), qd3 = 50,
                   moduli_sample_no = int(1e4))
cy_obj._set_flux_mc_paras(quanta_sample_no = int(1e8), quanta_range = qr, stdev = np.sqrt(1./1.))

#|%%--%%| <KPBh1spmQT|5f390XU7Pm>

_, udv_cmbn        = cy_obj.uniform_eval(type = "cmbn",        mrl = True)
_, udv_flux_index  = cy_obj.uniform_eval(type = "flux_index",  mrl = True)
_, udv_flux_normal = cy_obj.uniform_eval(type = "flux_normal", mrl = True)

#|%%--%%| <5f390XU7Pm|kdbFTt1Z9l>

usv_cmbn = cy_obj.uniform_integrate(udv_cmbn)
usv_flux_index = cy_obj.uniform_integrate(udv_flux_index)
usv_flux_normal = cy_obj.uniform_integrate(udv_flux_normal)

#|%%--%%| <kdbFTt1Z9l|JbpHLs7LA8>

print(usv_cmbn)
print(usv_flux_index)
print(usv_flux_normal)

#|%%--%%| <JbpHLs7LA8|gpC1EF8mDS>

print(cy.compute_inverse_kahler_metric([2,8]))
logdet, specgeo = cy_obj._ms_num([[2,8]])

print(specgeo)
print(cy_obj.calc_index_vacua_density_cmbn(specgeo))
print(cy_obj.calc_vacua_density_flux(specgeo, type="index", sampler="exp"))
print(cy_obj.calc_vacua_density_flux(specgeo, type="normal", sampler="exp"))

#|%%--%%| <gpC1EF8mDS|SuLPslMJlo>

cy_obj.qsno = int(1e4)
zzs, xxs = cy_obj._quanta_uniform_sample(sample_mode = "exp")
integrand_log, sign, _ = cy_obj._flux_num(zzs, xxs, specgeo)
integrand = np.exp(integrand_log).real * sign

x, y, z = np.abs(zzs), np.abs(xxs), integrand[0]
pts = np.column_stack([x, y]) 

print(pts[:10])
print(z[:10])

#|%%--%%| <SuLPslMJlo|guW8i5X1jE>

from scipy.interpolate import griddata

grid_n = 200
d = x.shape[1]

axes = [np.linspace(0, qr) for _ in range(d)]
y_axis = np.linspace(0, qr, grid_n)
axes.append(y_axis)
grid = np.meshgrid(*axes, indexing="xy")

z_grid = griddata(pts, z, tuple(grid), method="linear")

print(z_grid)

plt.contour(axes[0], axes[1], z_grid, levels=20)
plt.colorbar(label="z")
plt.show()

