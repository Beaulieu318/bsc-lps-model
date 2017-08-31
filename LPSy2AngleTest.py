import LPSLayer as ll
import LPSArgoAnalysis as laa
import numpy as np
import ArgoLayer as al
import matplotlib.pyplot as plt

TRegion = 1

density = np.array([24.7, 25.5, 26.2, 27.0, 27.5])
regions = [10, 30, 35, 50]
gammas = [0]+list((density[1:]-density[:-1])/np.mean(density+1000)*9.8)[:3]
lamda_e = 343
lamda_minmax = [280, 343]
H_0 = 500
H_4 = 1100
phi = [40, 25, 10]
we_max = 30 / float(365*24*60*60)
ekmanpumping = 'Actual'

LdataA = ll.LPSLayer(TRegion, regions, gammas, lamda_e, lamda_minmax, H_0, H_4, phi, we_max, ekmanpumping)
LdataA.Layers()

depths = LdataA.Depths()
pvs = LdataA.PV()
    
d_line = np.linspace(45, 35, depths[0].shape[1])

y_args = (np.abs(LdataA.Depths()[1] - d_line)).argmin(0)
x_args = np.arange(len(y_args))

lons = LdataA.Depths()[0][y_args, x_args]
lats = LdataA.Depths()[1][y_args, x_args]
d_depths = LdataA.Depths()[2][0][y_args, x_args]
d_pv = LdataA.PV()[2][0][y_args, x_args]

d_depths_sort_arg = d_depths.argsort()
d_depths = d_depths[d_depths_sort_arg]
d_depths_arg_max = np.nanargmax(d_depths)
d_depths = d_depths[:d_depths_arg_max]
d_pv = d_pv[d_depths_sort_arg]
d_pv = d_pv[:d_depths_arg_max]

plt.plot(-d_depths, d_pv, label='Slanted')
plt.legend(loc='best')
plt.ylabel('Potential Vorticity')
plt.xlabel('Depth (m)')
plt.title(r'PV against depth for slanted $y_2$')

poly = np.polyfit(-1/d_depths, d_pv, deg=1)

def G(x, poly):
    return poly[-1]+poly[-2]/x#+ poly[-3]*x**2# + poly[-4]*x**3
    
x_range = np.linspace(max(-d_depths), min(-d_depths), 50)
y_range = [G(x,poly) for x in x_range]
#plt.plot(x_range, y_range)

plt.show()