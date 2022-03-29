from scipy import special
from matplotlib import pyplot as plt
import numpy as np

ORDER = 3

def psi_array(sjn, z):
    return sjn * z

def xi_array(sjn, syn, z):
    return z * (sjn + complex(0, 1) * syn)

# MAIN
x = np.linspace(0.01, 10, 500)

sjn = []
syn = []
psin = []
xin = []
for i in range(ORDER):
    sjn.append(special.spherical_jn(i, x))
    syn.append(special.spherical_yn(i, x))
    psin.append(psi_array(sjn[i], x))
    xin.append(xi_array(sjn[i], syn[i], x))
    plt.plot(x, xin[i].imag)

plt.show()