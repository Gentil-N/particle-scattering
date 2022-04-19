from scipy import special
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
import math

#def sy_array_single(z, n):
#    y0 = -math.cos(z) / z
#    if n == 0:
#        return [y0]
#    y1 = -math.cos(z) / (z * z) - math.sin(z) / z
#    if n == 1:
#        return [y0, y1]
#    yn = [complex(0.0, 0.0)] * (n + 1)
#    yn[0] = y0
#    yn[1] = y1
#    for i in range(2, (n + 1)):
#        yn[i] = (2 * i + 1) / z * yn[i - 1] - yn[i - 2]
#    return yn
#
#def sy_array(z, n):
#    res = np.zeros(len(z), dtype=complex)
#    for i in range(0, len(z)):
#        a = sy_array_single(z[i], n + 1)
#        res[i] = a[n]
#    return res
#
#def sj_upward_reccurence(z, n):
#    count = n + 1
#    jn = [complex(0, 0)] * count
#    jn[0] = math.sin(z) / z
#    if n < 1:
#        return jn
#    jn[1] = math.sin(z) / z**2 - math.cos(z) / z
#    for i in range(1, count - 1):
#        jn[i + 1] = (2 * i + 1) / z * jn[i] - jn[i - 1]
#    return jn
#
#def sj_downward_reccurence(z, nl, nu):
#    count = nu - nl + 1
#    jn = [complex(0.0, 0.0)] * (count + 2);
#    jn[count + 1] = complex(0.0, 0.0)
#    jn[count] = complex(1.0, 0.0)
#    for i in range(count, 0, -1):
#        jn[i - 1] = (2 * i + 1) / z * jn[i] - jn[i + 1]
#    del jn[count:]
#    return jn
#
#def sj_array_single(z, n):
#    if abs(z) > n / 2:
#        return sj_upward_reccurence(z, n)
#    else:
#        PACK = 50
#        num_big_loop = int(n / PACK)
#        jn_all = []
#        for i in range(0, num_big_loop):
#            jn_all.append(sj_downward_reccurence(z, i * PACK, (i + 1) * PACK))
#        rest = n % PACK
#        if rest != 0:
#            jn_all.append(sj_downward_reccurence(z, n - rest, n))
#
#        jn = []
#        norm = math.sin(z) / z / jn_all[0][0]
#        for i in range(0, len(jn_all[0])):
#            jn.append(jn_all[0][i] * norm)
#        for i in range(1, len(jn_all)):
#            norm = jn[-1] / jn_all[i][0]
#            for k in range(1, len(jn_all[i])):
#                jn.append(jn_all[i][k] * norm);
#        return jn
#
#def sj_array(z, n):
#    res = np.zeros(len(z), dtype=complex)
#    for i in range(0, len(z)):
#        a = sj_array_single(z[i], n + 1)
#        res[i] = a[n]
#    return res

ORDER = 3
ORDER_MAX = ORDER + 1

def psi_array(sjn, z):
    return sjn * z

def xi_array(sjn, syn, z):
    return z * (sjn + complex(0, 1) * syn)

def d_array(psin, z):
    dn = [complex(0.0, 0.0)] * len(psin)
    #for i in range(1, len(dn)):
    #    dn[i] = psin[i - 1] / psin[i] - i / z
    psi_last = psin[-1]
    psi_before_last = psin[-2]
    dn[-1] = psi_before_last / psi_last - (len(psin) - 1) / z
    for i in range(len(psin) - 1, 0, -1):
        ndivz = i / z
        dn[i - 1] = ndivz - 1 / (dn[i] + ndivz)
    return dn

def a_array(m, psin, xin, dn, z):
    an = [complex(0.0, 0.0)] * len(psin)
    for i in range(1, len(an)):
        coeff_ai = dn[i] / m + i / z
        an[i] = (coeff_ai * psin[i] - psin[i - 1]) / (coeff_ai * xin[i] - xin[i - 1])
    return an

def b_array(m, psin, xin, dn, z):
    bn = [complex(0.0, 0.0)] * len(psin)
    for i in range(1, len(bn)):
        coeff_bi = dn[i] * m + i / z
        bn[i] = (coeff_bi * psin[i] - psin[i - 1]) / (coeff_bi * xin[i] - xin[i - 1])
    return bn

def load_ref_index(filename):
    ref_file = open(filename, 'r')
    lines = ref_file.readlines()
    data = []
    for line in lines:
        numbers = line.split(',')
        data.append((float(numbers[0].strip()) * 1e-6, float(numbers[1].strip()), float(numbers[2].strip())))
    ref_file.close()
    return data

def get_ref_index(data, wavelength):
    for i in range(len(data) - 1):
        if data[i][0] < wavelength and data[i + 1][0] > wavelength:
            alpha_re = (data[i + 1][1] - data[i][1]) / (data[i + 1][0] - data[i][0])
            gamma_re = data[i][1] - alpha_re * data[i][0]
            ref_re = alpha_re * wavelength + gamma_re
            alpha_im = (data[i + 1][2] - data[i][2]) / (data[i + 1][0] - data[i][0])
            gamma_im = data[i][2] - alpha_im * data[i][0]
            ref_im = alpha_im * wavelength + gamma_im
            return complex(ref_re, ref_im)
    return complex(data[-1][1], data[-1][2])

def compute_cross_sections(ref_indices_raw, wavelengths, particle_size):
    medium_n = 1.0
    upper_x = 2 * math.pi * medium_n * particle_size
    x = upper_x / wavelengths
    m = np.zeros(len(wavelengths), dtype=complex)
    mx = np.zeros(len(wavelengths), dtype=complex)
    for i in range(len(wavelengths)):
        m[i] = get_ref_index(ref_indices_raw, wavelengths[i]) / medium_n
        mx[i] = m[i] * x[i]

    sjn_x = []
    syn_x = []
    psin_x = []
    xin_x = []
    sjn_mx = []
    psin_mx = []
    dn_mx = []
    for i in range(ORDER_MAX):
        sjn_x.append(special.spherical_jn(i, x))
        syn_x.append(special.spherical_yn(i, x))
        psin_x.append(psi_array(sjn_x[i], x))
        xin_x.append(xi_array(sjn_x[i], syn_x[i], x))

        sjn_mx.append(special.spherical_jn(i, mx))
        psin_mx.append(psi_array(sjn_mx[i], mx))

    dn_mx = d_array(psin_mx, mx)
    an = a_array(m, psin_x, xin_x, dn_mx, x)
    bn = b_array(m, psin_x, xin_x, dn_mx, x)
    mul = (wavelengths / medium_n)**2 / (2 * math.pi) / (particle_size**2 * math.pi)

    part_res_csa = [0] * ORDER_MAX
    part_res_ext = [0] * ORDER_MAX
    part_res_csa_an = [0] * ORDER_MAX
    part_res_csa_bn = [0] * ORDER_MAX
    part_res_ext_an = [0] * ORDER_MAX
    part_res_ext_bn = [0] * ORDER_MAX
    for i in range(1, ORDER_MAX):
        part_res_csa_an[i] = (2 * i + 1) * an[i].real**2 + an[i].imag**2
        part_res_csa_bn[i] = (2 * i + 1) * bn[i].real**2 + bn[i].imag**2
        part_res_ext_an[i] = (2 * i + 1) * an[i].real
        part_res_ext_bn[i] = (2 * i + 1) * bn[i].real
        part_res_csa[i] = part_res_csa_an[i] + part_res_csa_bn[i]
        part_res_ext[i] = part_res_ext_an[i] + part_res_ext_bn[i]
    res_csa_an = mul * sum(part_res_csa_an)
    res_csa_bn = mul * sum(part_res_csa_bn)
    res_ext_an = mul * sum(part_res_ext_an)
    res_ext_bn = mul * sum(part_res_ext_bn)
    res_csa = mul * sum(part_res_csa)
    res_ext = mul * sum(part_res_ext)

    return (res_csa, res_ext, res_csa_an, res_csa_bn, res_ext_an, res_ext_bn)

def plot_surface_sca_ext():
    PARTSIZE_LOWER = 40e-9
    PARTSIZE_UPPER = 160e-9
    partsizes = np.linspace(PARTSIZE_LOWER, PARTSIZE_UPPER, DIV)
    scattering_cross_section = np.zeros((len(partsizes), len(WAVELENGTHS)))
    extinction_cross_section = np.zeros((len(partsizes), len(WAVELENGTHS)))
    for i in range(len(partsizes)):
        print(partsizes[i])
        res = compute_cross_sections(REF_INDICES_RAW, WAVELENGTHS, partsizes[i])
        scattering_cross_section[i] = res[0]
        extinction_cross_section[i] = res[1]
    fig = plt.figure(num=0, figsize=(12, 5))
    axs = fig.subplots(nrows=1, ncols=2)
    axs[0].set_title("Scattering Cross Section")
    axs[0].contourf(WAVELENGTHS, partsizes, scattering_cross_section, cmap='inferno', levels=70)
    axs[0].set(xlabel="wavelength", ylabel="particle radius")
    axs[0].grid()
    axs[1].set_title("Extinction Cross Section")
    axs[1].contourf(WAVELENGTHS, partsizes, extinction_cross_section, cmap='inferno', levels=70)
    axs[1].set(xlabel="wavelength", ylabel="particle radius")
    axs[1].grid()
    fig.colorbar(mappable=ScalarMappable(norm=Normalize(vmin=0, vmax=10), cmap='inferno'), ax=axs[0])
    fig.colorbar(mappable=ScalarMappable(norm=Normalize(vmin=0, vmax=10), cmap='inferno'), ax=axs[1])
    plt.show()

def plot_coeff_sca_ext(particle_size):
    res = compute_cross_sections(REF_INDICES_RAW, WAVELENGTHS, particle_size)
    fig = plt.figure(num=0, figsize=(13, 4))
    axs = fig.subplots(nrows=1, ncols=2)
    #axs.plot(WAVELENGTHS, res[0])
    #axs.plot(WAVELENGTHS, res[1])
    axs[0].set_title("Scattering Cross Section with Coeff")
    axs[0].plot(WAVELENGTHS, res[0], label="Sca Total")
    axs[0].plot(WAVELENGTHS, res[2], label="\'an\'")
    axs[0].plot(WAVELENGTHS, res[3], label="\'bn\'")
    axs[0].set(xlabel="wavelength")
    axs[0].legend()
    axs[0].grid()
    axs[1].set_title("Extinction Cross Section with coeff")
    axs[1].plot(WAVELENGTHS, res[1], label="Ext Total")
    axs[1].plot(WAVELENGTHS, res[4], label="\'an\'")
    axs[1].plot(WAVELENGTHS, res[5], label="\'bn\'")
    axs[1].set(xlabel="wavelength")
    axs[1].legend()
    axs[1].grid()
    plt.show()

####################### MAIN #######################

DIV = 500
REF_INDICES_RAW = load_ref_index("./res/refractive-index-silicon.csv")
WAVELENGTHS = np.linspace(REF_INDICES_RAW[0][0], REF_INDICES_RAW[-1][0], DIV)

#plot_surface_sca_ext()
#plot_coeff_sca_ext(90e-9)

x = np.linspace(0, 30, 1000)
sjn_x = []
syn_x = []
for i in range(ORDER_MAX):
    sjn_x.append(special.spherical_jn(i, x))
    syn_x.append(special.spherical_yn(i, x))
plt.plot(x, sjn_x[0], label="j0")
plt.plot(x, sjn_x[1], label="j1")
plt.plot(x, syn_x[0], label="y0")
plt.plot(x, syn_x[1], label="y1")
plt.ylim(top=1.2, bottom=-1)
plt.xlim(right=20)
plt.xlim(left=0)
plt.grid()
plt.legend()
plt.show()