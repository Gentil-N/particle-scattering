from scipy import special
from scipy import integrate
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
import math
import cmath
import plyfile as pf

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
#
#def sj_downward_reccurence(z, nu):
#    count = nu + 1
#    jn = [0.0] * (count + 2);
#    jn[count + 1] = 0.0
#    jn[count] = 1.0
#    for i in range(count, 0, -1):
#        jn[i - 1] = (2 * i + 1) / z * jn[i] - jn[i + 1]
#        if (i - 1) < 3 and z > 8 and z < 10:
#            print(z, " ", i - 1, " ", jn[i - 1])
#    norm = math.sin(z) / z / jn[0]
#    for i in range(count):
#        jn[i] *= norm;
#    return jn

def sy_array(z, n):
    y0 = -cmath.cos(z) / z;
    if n == 0:
        return [y0]
    y1 = -cmath.cos(z) / (z * z) - cmath.sin(z) / z
    if n == 1:
        return [y0, y1]
    yn = [complex(0.0, 0.0)] * (n + 1)
    yn[0] = y0
    yn[1] = y1
    for i in range(2, n + 1):
        yn[i] = (2 * i + 1) / z * yn[i - 1] - yn[i - 2]
    return yn

def sj_upward_reccurence(z, n):
    count = n + 1
    jn = [complex(0.0, 0.0)] * count
    jn[0] = cmath.sin(z) / z
    jn[1] = cmath.sin(z) / (z**2) - cmath.cos(z) / z
    for i in range(1, count - 1):
        jn[i + 1] = (2 * i + 1) / z * jn[i] - jn[i - 1]
    return jn

def sj_downward_reccurence(z, nl, nu):
    count = nu - nl + 1
    jn = [complex(0.0, 0.0)] * (count + 2)
    jn[count + 1] = complex(0.0, 0.0)
    jn[count] = complex(1.0, 0.0)
    for i in range(count, 0, -1):
        jn[i - 1] = (2 * i + 1) / z * jn[i] - jn[i + 1]
    jn = jn[:count]
    return jn

def sj_array(z, n):
    if abs(z) > n / 2.0:
        return sj_upward_reccurence(z, n)
    else:
        PACK = 50
        num_big_loop = int(n / PACK)
        jn_all = []
        for i in range(0,num_big_loop):
            jn_all.append(sj_downward_reccurence(z, i * PACK, (i + 1) * PACK))
        rest = n % PACK
        if rest != 0:
            jn_all.append(sj_downward_reccurence(z, n - rest, n))

        jn = []
        norm = cmath.sin(z) / z / jn_all[0][0]
        for i in range(len(jn_all[0])):
            jn.append(jn_all[0][i] * norm)
        for i in range(1,len(jn_all)):
            norm = jn[-1] / jn_all[i][0]
            for k in range(1, len(jn_all[i])):
                jn.append(jn_all[i][k] * norm)
        return jn

ORDER = 3
ORDER_LEN = ORDER + 1

def psi_array(sjn, z):
    return sjn * z

def xi_array(sjn, syn, z):
    return z * (sjn + complex(0, 1) * syn)

def xi_der_array(xin, z):
    xin_der = [complex(0.0)] * len(xin)
    for i in range(1, len(xin_der)):
        xin_der[i] = xin[i - 1] - i * xin[i] / z
    return xin_der
    

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

def pi_tau_array(n_len, theta):
    assert n_len >= 3
    mu = math.cos(theta)
    pi_tau_n = [(0.0, 0.0)] * n_len
    pi_tau_n[0] = (0.0, 0.0)
    pi_tau_n[1] = (1.0, mu * 1.0 - 2 * 0.0)
    for i in range(2, len(pi_tau_n)):
        curr_pi = (2 * i - 1) / (i - 1) * mu * pi_tau_n[i - 1][0] - i / (i - 1) * pi_tau_n[i - 2][0]
        pi_tau_n[i] = (curr_pi, i * mu * curr_pi - (i + 1) * pi_tau_n[i - 1][0])
    return pi_tau_n

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

def en(n):
    return complex(0.0, 1.0)**n * (2 * n + 1) / (n * (n + 1))

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
        if data[i][0] <= wavelength and data[i + 1][0] >= wavelength:
            alpha_re = (data[i + 1][1] - data[i][1]) / (data[i + 1][0] - data[i][0])
            gamma_re = data[i][1] - alpha_re * data[i][0]
            ref_re = alpha_re * wavelength + gamma_re
            alpha_im = (data[i + 1][2] - data[i][2]) / (data[i + 1][0] - data[i][0])
            gamma_im = data[i][2] - alpha_im * data[i][0]
            ref_im = alpha_im * wavelength + gamma_im
            return complex(ref_re, ref_im)
    return complex(data[0][1], data[0][2])

def theta_func_first(an, bn, xin, xin_der, theta):
    pi_tau_n = pi_tau_array(len(an), theta)
    sum_first = 0
    sum_second = 0
    for i in range(1, len(an)):
        sum_first += en(i) * (complex(0.0, 1.0) * an[i] * xin_der[i] * pi_tau_n[i][1] - bn[i] * xin[i] * pi_tau_n[i][0])
        sum_second += en(i) * (complex(0.0, 1.0) * bn[i] * xin_der[i] * pi_tau_n[i][0] - an[i] * xin[i] * pi_tau_n[i][1])
    return (sum_first * sum_second.conjugate()).real * math.sin(theta)

def theta_func_second(an, bn, xin, xin_der, theta):
    pi_tau_n = pi_tau_array(len(an), theta)
    sum_first = 0
    sum_second = 0
    for i in range(1, len(an)):
        sum_first += en(i) * (bn[i] * xin[i] * pi_tau_n[i][1] - complex(0.0, 1.0) * an[i] * xin_der[i] * pi_tau_n[i][0])
        sum_second += en(i) * (complex(0.0, 1.0) * bn[i] * xin_der[i] * pi_tau_n[i][1] - an[i] * xin[i] * pi_tau_n[i][0])
    return (sum_first * sum_second.conjugate()).real * math.sin(theta)

def theta_func(an, bn, xin, xin_der, theta):
    pi_tau_n = pi_tau_array(len(an), theta)
    sum_first = 0
    sum_second = 0
    sum_third = 0
    sum_fourth = 0
    for i in range(1, len(an)):
        sum_first += en(i) * (complex(0.0, 1.0) * an[i] * xin_der[i] * pi_tau_n[i][1] - bn[i] * xin[i] * pi_tau_n[i][0])
        sum_second += en(i) * (complex(0.0, 1.0) * bn[i] * xin_der[i] * pi_tau_n[i][0] - an[i] * xin[i] * pi_tau_n[i][1])
        sum_third += en(i) * (bn[i] * xin[i] * pi_tau_n[i][1] - complex(0.0, 1.0) * an[i] * xin_der[i] * pi_tau_n[i][0])
        sum_fourth += en(i) * (complex(0.0, 1.0) * bn[i] * xin_der[i] * pi_tau_n[i][1] - an[i] * xin[i] * pi_tau_n[i][0])
    return (sum_first * sum_second.conjugate()).real * math.sin(theta) - (sum_third * sum_fourth.conjugate()).real * math.sin(theta)

def trapz(func, inf, sup, div):
    step = (sup - inf) / div
    res = 0
    for i in range(div):
        a = i * step + inf
        b = a + step
        res += (func(a) + func(b)) * (step) / 2
    return res

def compute_integrated_scattering_cross_section(phi_inf, phi_sup, theta_inf, theta_sup, ref_indices_raw, wavelengths, particle_size):
    medium_n = 1.0
    upper_x = 2 * math.pi * medium_n * particle_size
    res = np.zeros(len(wavelengths), dtype=float)
    for j in range(len(wavelengths)):
        #print(wavelengths[j])
        x = upper_x / wavelengths[j]
        m = get_ref_index(ref_indices_raw, wavelengths[j]) / medium_n
        mx = m * x
    
        sjn_x = []
        syn_x = []
        psin_x = []
        xin_x = []
        xin_der_x = []
        sjn_mx = []
        psin_mx = []
        dn_mx = []
        for i in range(ORDER_LEN):
            sjn_x.append(special.spherical_jn(i, x))
            syn_x.append(special.spherical_yn(i, x))
            psin_x.append(psi_array(sjn_x[i], x))
            xin_x.append(xi_array(sjn_x[i], syn_x[i], x))

            sjn_mx.append(special.spherical_jn(i, mx))
            psin_mx.append(psi_array(sjn_mx[i], mx))

        xin_der_x = xi_der_array(xin_x, x)
        dn_mx = d_array(psin_mx, mx)
        an = a_array(m, psin_x, xin_x, dn_mx, x)
        bn = b_array(m, psin_x, xin_x, dn_mx, x)

        #phi_mul = wavelengths[j]**2 / (4 * math.pi**2 * 3*10e8 * 4 * math.pi * 10e-7)
        integ_phi_first = trapz(lambda phi: math.cos(phi)**2, phi_inf, phi_sup, 50)
        integ_phi_second = trapz(lambda phi: math.sin(phi)**2, phi_inf, phi_sup, 50)
        mul = 0.5 * wavelengths[j]**2 * 0.318 / (2 * math.pi) / (particle_size**2 * math.pi)
        integ_theta_first = trapz(lambda theta: theta_func_first(an, bn, xin_x, xin_der_x, theta), theta_inf, theta_sup, 50)
        integ_theta_second = trapz(lambda theta: theta_func_second(an, bn, xin_x, xin_der_x, theta), theta_inf, theta_sup, 50)
        #integ_theta = trapz(lambda theta: theta_func(an, bn, xin_x, xin_der_x, theta), theta_inf, theta_sup, 100)

        res[j] = (integ_phi_first * integ_theta_first - integ_phi_second * integ_theta_second) * mul
        #mul = wavelengths[j]**2 / (2 * math.pi) / (particle_size**2 * math.pi)
        #res_an[j] = mul * (3 * an[1].real**2 + an[1].imag**2 + 5 * an[2].real**2 + an[2].imag**2 + 7 * an[3].real**2 + an[3].imag**2)
        #res_bn[j] = mul * (3 * bn[1].real**2 + bn[1].imag**2 + 5 * bn[2].real**2 + bn[2].imag**2 + 7 * bn[3].real**2 + bn[3].imag**2)

    return res
    
def distance(point1, point2):
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

def clamp(val, mmin, mmax):
    return max(min(val, mmax), mmin)

def function_value(an, bn, xin_x, xin_der_x, phi_theta):
    return math.cos(phi_theta[0])**2 * theta_func_first(an, bn, xin_x, xin_der_x, phi_theta[1]) - math.sin(phi_theta[0])**2 * theta_func_second(an, bn, xin_x, xin_der_x, phi_theta[1])

def compute_cross_sections_triangle(mul, phi_theta_1, phi_theta_2, phi_theta_3, res1, res2, res3):
    a = 0
    b = 0
    c = 0
    d = 0
    e = 0
    f = 0
    if res1 < res2:
        if res3 < res1:
            d = res3
            e = res1
            f = res2
            a = distance(phi_theta_3, phi_theta_1)
            b = distance(phi_theta_3, phi_theta_2)
            c = distance(phi_theta_1, phi_theta_2)
        elif res2 < res3:
            d = res1
            e = res2
            f = res3
            a = distance(phi_theta_1, phi_theta_2)
            b = distance(phi_theta_1, phi_theta_3)
            c = distance(phi_theta_3, phi_theta_2)
        else:
            d = res1
            e = res3
            f = res2
            a = distance(phi_theta_1, phi_theta_3)
            b = distance(phi_theta_1, phi_theta_2)
            c = distance(phi_theta_3, phi_theta_2)
    else:
        if res3 < res2:
            d = res3
            e = res2
            f = res1
            a = distance(phi_theta_3, phi_theta_2)
            b = distance(phi_theta_3, phi_theta_1)
            c = distance(phi_theta_1, phi_theta_2)
        elif res1 < res3:
            d = res2
            e = res1
            f = res3
            a = distance(phi_theta_2, phi_theta_1)
            b = distance(phi_theta_2, phi_theta_3)
            c = distance(phi_theta_1, phi_theta_3)
        else:
            d = res2
            e = res3
            f = res1
            a = distance(phi_theta_2, phi_theta_3)
            b = distance(phi_theta_2, phi_theta_1)
            c = distance(phi_theta_3, phi_theta_2)
    if a == 0 or c == 0:
        return (0, 0)
    h = a * math.sin(math.acos(clamp((a**2 + c**2 - b**2) / (2 * a * c), -1, 1)))
    v1 = 0.5 * c * h * d
    v2 = (1/6) * (f + e) * c * h
    return (mul * (v1 + v2), 0.5 * c * h)

def compute_cross_sections_triangle_simple(mul, phi_theta_1, phi_theta_2, phi_theta_3, res1, res2, res3):
    a = 0
    b = 0
    c = 0
    d = 0
    e = 0
    f = 0
    if res1 < res2:
        if res3 < res1:
            d = res3
            e = res1
            f = res2
            a = distance(phi_theta_3, phi_theta_1)
            b = distance(phi_theta_3, phi_theta_2)
            c = distance(phi_theta_1, phi_theta_2)
        elif res2 < res3:
            d = res1
            e = res2
            f = res3
            a = distance(phi_theta_1, phi_theta_2)
            b = distance(phi_theta_1, phi_theta_3)
            c = distance(phi_theta_3, phi_theta_2)
        else:
            d = res1
            e = res3
            f = res2
            a = distance(phi_theta_1, phi_theta_3)
            b = distance(phi_theta_1, phi_theta_2)
            c = distance(phi_theta_3, phi_theta_2)
    else:
        if res3 < res2:
            d = res3
            e = res2
            f = res1
            a = distance(phi_theta_3, phi_theta_2)
            b = distance(phi_theta_3, phi_theta_1)
            c = distance(phi_theta_1, phi_theta_2)
        elif res1 < res3:
            d = res2
            e = res1
            f = res3
            a = distance(phi_theta_2, phi_theta_1)
            b = distance(phi_theta_2, phi_theta_3)
            c = distance(phi_theta_1, phi_theta_3)
        else:
            d = res2
            e = res3
            f = res1
            a = distance(phi_theta_2, phi_theta_3)
            b = distance(phi_theta_2, phi_theta_1)
            c = distance(phi_theta_3, phi_theta_2)
    if a == 0 or c == 0:
        return 0
    h = a * math.sin(math.acos(clamp((a**2 + c**2 - b**2) / (2 * a * c), -1, 1)))
    return mul * 0.5 * c * h * ((d + e + f) / 3)

def compute_cross_sections_square(mul, c, d, e, f, phi_step, theta_step):
    h = (c + d + e + f) / 4
    return phi_step * theta_step * h

def compute_cross_sections_whole_by_triangle(ref_indices_raw, wavelengths, particle_size):
    medium_n = 1.0
    upper_x = 2 * math.pi * medium_n * particle_size
    res = np.zeros(len(wavelengths), dtype=float)
    print(particle_size, ": 0.0 %", end='\r', flush=True)
    num_wavelengths = len(wavelengths)
    area = 0
    for j in range(num_wavelengths):
        #print(wavelengths[j])
        x = upper_x / wavelengths[j]
        m = get_ref_index(ref_indices_raw, wavelengths[j]) / medium_n
        mx = m * x
    
        sjn_x = []
        syn_x = []
        psin_x = []
        xin_x = []
        xin_der_x = []
        sjn_mx = []
        psin_mx = []
        dn_mx = []
        for i in range(ORDER_LEN):
            sjn_x.append(special.spherical_jn(i, x))
            syn_x.append(special.spherical_yn(i, x))
            psin_x.append(psi_array(sjn_x[i], x))
            xin_x.append(xi_array(sjn_x[i], syn_x[i], x))

            sjn_mx.append(special.spherical_jn(i, mx))
            psin_mx.append(psi_array(sjn_mx[i], mx))

        xin_der_x = xi_der_array(xin_x, x)
        dn_mx = d_array(psin_mx, mx)
        an = a_array(m, psin_x, xin_x, dn_mx, x)
        bn = b_array(m, psin_x, xin_x, dn_mx, x)
        mul = 0.5 * wavelengths[j]**2 * 0.318 / (2 * math.pi) / (particle_size**2 * math.pi)

        curr_res = 0
        div = 50
        phi_step = 2 * math.pi / div
        theta_step = math.pi / div
        fn_values = []
        for i in range(div + 1):
            fn_row = []
            for k in range(div + 1):
                curr_phi = i * phi_step
                curr_theta = k * theta_step
                fn_row.append(function_value(an, bn, xin_x, xin_der_x, [curr_phi, curr_theta]))
            fn_values.append(fn_row)

        area = 0
        for i in range(div):
            for k in range(div):
                curr_phi = i * phi_step
                curr_theta = k * theta_step
                a = [curr_phi, curr_theta]
                b = [curr_phi + phi_step, curr_theta]
                c = [curr_phi + phi_step, curr_theta + theta_step]
                d = [curr_phi, curr_theta + theta_step]
                #resa = function_value(an, bn, xin_x, xin_der_x, a)
                #resb = function_value(an, bn, xin_x, xin_der_x, b)
                #resc = function_value(an, bn, xin_x, xin_der_x, c)
                #resd = function_value(an, bn, xin_x, xin_der_x, d)
                comp = compute_cross_sections_triangle(mul, a, b, c, abs(fn_values[i][k]), abs(fn_values[i + 1][k]), abs(fn_values[i + 1][k + 1]))
                curr_res += comp[0]
                area += comp[1]
                comp = compute_cross_sections_triangle(mul, a, c, d, abs(fn_values[i][k]), abs(fn_values[i + 1][k + 1]), abs(fn_values[i][k + 1]))
                curr_res += comp[0]
                area += comp[1]
                #curr_res += compute_cross_sections_square(mul, fn_values[i][k], fn_values[i + 1][k], fn_values[i + 1][k + 1], fn_values[i][k + 1], phi_step, theta_step)
        res[j] = curr_res #* 0.581
        print(particle_size, ":", float(int(j / num_wavelengths * 1000)) / 10, " %", end='\r', flush=True)
    print(particle_size, ": 100.0 %", flush=True)
    print(area)
    # area ne correspond pas du tout ! Theorie : 19.74, WholeTri : 21., SphereTri : 42.
    return res

def compute_cross_sections(ref_indices_raw, wavelengths, particle_size):
    medium_n = 1.0
    upper_x = 2 * math.pi * medium_n * particle_size
    res_csa = np.zeros(len(wavelengths), dtype=float)
    res_ext = np.zeros(len(wavelengths), dtype=float)
    res_csa_an = np.zeros(len(wavelengths), dtype=float)
    res_csa_bn = np.zeros(len(wavelengths), dtype=float)
    res_ext_an = np.zeros(len(wavelengths), dtype=float)
    res_ext_bn = np.zeros(len(wavelengths), dtype=float)
    for j in range(len(wavelengths)):
        x = upper_x / wavelengths[j]
        m = get_ref_index(ref_indices_raw, wavelengths[j]) / medium_n
        mx = m * x

        sjn_x = []
        syn_x = []
        psin_x = []
        xin_x = []
        sjn_mx = []
        psin_mx = []
        dn_mx = []
        for i in range(ORDER_LEN):
            sjn_x.append(special.spherical_jn(i, x))
            syn_x.append(special.spherical_yn(i, x))
            psin_x.append(psi_array(sjn_x[i], x))
            xin_x.append(xi_array(sjn_x[i], syn_x[i], x))

            sjn_mx.append(special.spherical_jn(i, mx))
            psin_mx.append(psi_array(sjn_mx[i], mx))

        dn_mx = d_array(psin_mx, mx)
        an = a_array(m, psin_x, xin_x, dn_mx, x)
        bn = b_array(m, psin_x, xin_x, dn_mx, x)
        mul = (wavelengths[j] / medium_n)**2 / (2 * math.pi) / (particle_size**2 * math.pi)

        part_res_csa = [0] * ORDER_LEN
        part_res_ext = [0] * ORDER_LEN
        part_res_csa_an = [0] * ORDER_LEN
        part_res_csa_bn = [0] * ORDER_LEN
        part_res_ext_an = [0] * ORDER_LEN
        part_res_ext_bn = [0] * ORDER_LEN
        for i in range(1, ORDER_LEN):
            part_res_csa_an[i] = (2 * i + 1) * (an[i].real**2 + an[i].imag**2)
            part_res_csa_bn[i] = (2 * i + 1) * (bn[i].real**2 + bn[i].imag**2)
            part_res_ext_an[i] = (2 * i + 1) * an[i].real
            part_res_ext_bn[i] = (2 * i + 1) * bn[i].real
            part_res_csa[i] = part_res_csa_an[i] + part_res_csa_bn[i]
            part_res_ext[i] = part_res_ext_an[i] + part_res_ext_bn[i]
        res_csa_an[j] = mul * sum(part_res_csa_an)
        res_csa_bn[j] = mul * sum(part_res_csa_bn)
        res_ext_an[j] = mul * sum(part_res_ext_an)
        res_ext_bn[j] = mul * sum(part_res_ext_bn)
        res_csa[j] = mul * sum(part_res_csa)
        res_ext[j] = mul * sum(part_res_ext)
    return (res_csa, res_ext, res_csa_an, res_csa_bn, res_ext_an, res_ext_bn)

def transform_cartesian_to_spherical_angles(x, y, z):
    r = math.sqrt(x**2 + y**2 + z**2)
    cos_theta = z / r
    theta = math.acos(cos_theta)
    if x > 0:
        return (math.atan(y / x), theta)
    if x < 0 and y >= 0:
        return (math.atan(y / x) + math.pi, theta)
    if x < 0 and y < 0:
        return (math.atan(y / x) - math.pi, theta)
    if x == 0 and y > 0:
        return (math.pi / 2, theta)
    if x == 0 and y < 0:
        return (-math.pi / 2, theta)
    else:
        return (0.0, theta)

def is_white(r, g, b):
    return (r == 255 and g == 255 and b == 255)

def triangle_has_no_white(vert_0, vert_1, vert_2):
    if not(is_white(vert_0[3], vert_0[4], vert_0[5])) and not(is_white(vert_1[3], vert_1[4], vert_1[5])) and not(is_white(vert_2[3],vert_2[4], vert_2[5])):
        return True
    return False

def load_selected_triangle(filename):
    data = pf.PlyData.read(filename)
    point_coords = []
    for vertex in data.elements[0]:
        if not(is_white(vertex[3], vertex[4], vertex[5])):
            point_coords.append(transform_cartesian_to_spherical_angles(vertex[0], vertex[1], vertex[2]))
        else:
            point_coords.append((7.7, 7.7))
    indices = []
    #print(len(data.elements[1]))
    for face in data.elements[1]:
        if len(face[0]) == 3:
            vert_0 = data.elements[0][face[0][0]]
            vert_1 = data.elements[0][face[0][1]]
            vert_2 = data.elements[0][face[0][2]]
            if triangle_has_no_white(vert_0, vert_1, vert_2):
                indices.append(face[0][0])
                indices.append(face[0][1])
                indices.append(face[0][2])
        elif len(face[0]) == 4:
            vert_0 = data.elements[0][face[0][0]]
            vert_1 = data.elements[0][face[0][1]]
            vert_2 = data.elements[0][face[0][2]]
            vert_3 = data.elements[0][face[0][3]]
            if triangle_has_no_white(vert_0, vert_1, vert_2):
                indices.append(face[0][0])
                indices.append(face[0][1])
                indices.append(face[0][2])
            if triangle_has_no_white(vert_2, vert_3, vert_0):
                indices.append(face[0][2])
                indices.append(face[0][3])
                indices.append(face[0][0])
    print(len(point_coords))
    print(len(indices))
    return (point_coords, indices)

def compute_cross_section_by_triangle(ref_indices_raw, wavelengths, particle_size, point_coords, indices):
    medium_n = 1.0
    upper_x = 2 * math.pi * medium_n * particle_size
    res = np.zeros(len(wavelengths), dtype=float)
    print(particle_size, ": 0.0 %", end='\r', flush=True)
    num_wavelengths = len(wavelengths)
    area = 0
    for j in range(num_wavelengths):
        #print(wavelengths[j])
        x = upper_x / wavelengths[j]
        m = get_ref_index(ref_indices_raw, wavelengths[j]) / medium_n
        mx = m * x
    
        sjn_x = []
        syn_x = []
        psin_x = []
        xin_x = []
        xin_der_x = []
        sjn_mx = []
        psin_mx = []
        dn_mx = []
        for i in range(ORDER_LEN):
            sjn_x.append(special.spherical_jn(i, x))
            syn_x.append(special.spherical_yn(i, x))
            psin_x.append(psi_array(sjn_x[i], x))
            xin_x.append(xi_array(sjn_x[i], syn_x[i], x))

            sjn_mx.append(special.spherical_jn(i, mx))
            psin_mx.append(psi_array(sjn_mx[i], mx))

        xin_der_x = xi_der_array(xin_x, x)
        dn_mx = d_array(psin_mx, mx)
        an = a_array(m, psin_x, xin_x, dn_mx, x)
        bn = b_array(m, psin_x, xin_x, dn_mx, x)
        mul = 0.5 * wavelengths[j]**2 * 0.318 / (2 * math.pi) / (particle_size**2 * math.pi)

        curr_res = 0
        area = 0
        fn_values = []
        for coord in point_coords:
            fn_values.append(function_value(an, bn, xin_x, xin_der_x, coord))

        for i in range(0, int(len(indices)), 3):
            a = point_coords[indices[i]]
            b = point_coords[indices[i + 1]]
            c = point_coords[indices[i + 2]]
            comp = compute_cross_sections_triangle(mul, a, b, c, abs(fn_values[indices[i]]), abs(fn_values[indices[i + 1]]), abs(fn_values[indices[i + 2]]))
            curr_res += comp[0]
            area += comp[1]
        print("\n", area)
        exit()

        res[j] = curr_res
        print(particle_size, ":", float(int(j / num_wavelengths * 1000)) / 10, " %", end='\r', flush=True)
    print(particle_size, ": 100.0 %", flush=True)
    print(area)
    return res

def plot_surface_sca_ext():
    PARTSIZE_LOWER = 50e-9
    PARTSIZE_UPPER = 100e-9
    partsizes = np.linspace(PARTSIZE_LOWER, PARTSIZE_UPPER, DIV)
    scattering_cross_section = np.zeros((len(partsizes), len(WAVELENGTHS)))
    extinction_cross_section = np.zeros((len(partsizes), len(WAVELENGTHS)))
    for i in range(len(partsizes)):
        print(partsizes[i])
        res = compute_cross_sections(REF_INDICES_RAW, WAVELENGTHS, partsizes[i])
        scattering_cross_section[i] = res[0]
        extinction_cross_section[i] = res[1]
    fig0 = plt.figure(num=0)
    fig1 = plt.figure(num=1)
    ax0 = fig0.subplots(nrows=1, ncols=1)
    ax1 = fig1.subplots(nrows=1, ncols=1)
    ax0.set_title("Scattering Cross Section")
    ax0.contourf(WAVELENGTHS, partsizes, scattering_cross_section, cmap='inferno', levels=70)
    ax0.set(xlabel="wavelength", ylabel="particle radius")
    ax1.set_title("Extinction Cross Section")
    ax1.contourf(WAVELENGTHS, partsizes, extinction_cross_section, cmap='inferno', levels=70)
    ax1.set(xlabel="wavelength", ylabel="particle radius")
    fig0.colorbar(mappable=ScalarMappable(norm=Normalize(vmin=0, vmax=10), cmap='inferno'), ax=ax0)
    fig1.colorbar(mappable=ScalarMappable(norm=Normalize(vmin=0, vmax=10), cmap='inferno'), ax=ax1)
    plt.show()

def plot_coeff_sca_ext(particle_size):
    res = compute_cross_sections(REF_INDICES_RAW, WAVELENGTHS, particle_size)
    fig0 = plt.figure(num=0)
    fig1 = plt.figure(num=1)
    ax0 = fig0.subplots(nrows=1, ncols=1)
    ax1 = fig1.subplots(nrows=1, ncols=1)
    #axs.plot(WAVELENGTHS, res[0])
    #axs.plot(WAVELENGTHS, res[1])
    ax0.set_title("Scattering Cross Section with Coeff")
    ax0.plot(WAVELENGTHS, res[0], label="Sca Total")
    ax0.plot(WAVELENGTHS, res[2], label="\'an\'", color='red')
    ax0.plot(WAVELENGTHS, res[3], label="\'bn\'", color='green')
    ax0.set(xlabel="wavelength")
    ax0.legend()
    ax0.grid()
    ax1.set_title("Extinction Cross Section with coeff")
    ax1.plot(WAVELENGTHS, res[1], label="Ext Total")
    ax1.plot(WAVELENGTHS, res[4], label="\'an\'", color='red')
    ax1.plot(WAVELENGTHS, res[5], label="\'bn\'", color='green')
    ax1.set(xlabel="wavelength")
    ax1.legend()
    ax1.grid()
    plt.show()

def plot_sca_ext(particle_size):
    res = compute_cross_sections(REF_INDICES_RAW, WAVELENGTHS, particle_size)
    fig0 = plt.figure(num=0)
    ax0 = fig0.subplots(nrows=1, ncols=1)
    ax0.set_title("Scattering/Extinction Cross Sections")
    ax0.plot(WAVELENGTHS, res[0], label="Sca")
    ax0.plot(WAVELENGTHS, res[1], label="Ext")
    ax0.plot(WAVELENGTHS, res[1] - res[0], label="Abs")
    ax0.set(xlabel="wavelength")
    ax0.legend()
    ax0.grid()
    plt.show()

def plot_ref_indices():
    fig0 = plt.figure(num=0)
    ax0 = fig0.subplots(nrows=1, ncols=1)
    ref_real = []
    ref_imag = []
    for wavelength in WAVELENGTHS:
        ref_real.append(get_ref_index(REF_INDICES_RAW, wavelength).real)
        ref_imag.append(get_ref_index(REF_INDICES_RAW, wavelength).imag)
    ax0.set_title("Refractive Indices")
    ax0.plot(WAVELENGTHS, ref_real, color='cyan')
    ax0.plot(WAVELENGTHS, ref_imag, color='#EBC23A')
    ax0.grid()
    plt.show()

def plot_integ_sca(particle_size):
    res0 = compute_integrated_scattering_cross_section(0, 2 * math.pi, math.pi / 2, math.pi, REF_INDICES_RAW, WAVELENGTHS, particle_size)
    #res1 = compute_integrated_scattering_cross_section(0, 2 * math.pi, 0, math.pi * 0.57, REF_INDICES_RAW, WAVELENGTHS, particle_size)
    #res_exact = compute_cross_sections(REF_INDICES_RAW, WAVELENGTHS, particle_size)
    fig0 = plt.figure(num=0)
    ax0 = fig0.subplots(nrows=1, ncols=1)
    ax0.set_title("Scattering Cross Sections")
    ax0.plot(WAVELENGTHS, res0, label="Integrated")
    #ax0.plot(WAVELENGTHS, res1, label="Integrated")
    #ax0.plot(WAVELENGTHS, res_exact[0], label="Exact")
    ax0.set(xlabel="wavelength")
    ax0.legend()
    ax0.grid()
    plt.show()

def plot_integ_sca_surface():
    PARTSIZE_LOWER = 50e-9
    PARTSIZE_UPPER = 100e-9
    partsizes = np.linspace(PARTSIZE_LOWER, PARTSIZE_UPPER, 50)
    scattering_cross_section = np.zeros((len(partsizes), len(WAVELENGTHS)))
    for i in range(len(partsizes)):
        print(partsizes[i])
        scattering_cross_section[i] = compute_integrated_scattering_cross_section(0, 2 * math.pi, 0.0 * math.pi, 1.0 * math.pi, REF_INDICES_RAW, WAVELENGTHS, partsizes[i])
    fig0 = plt.figure(num=0)
    ax0 = fig0.subplots(nrows=1, ncols=1)
    ax0.set_title("Scattering Cross Section")
    ax0.contourf(WAVELENGTHS, partsizes, scattering_cross_section, cmap='inferno', levels=70)
    ax0.set(xlabel="wavelength", ylabel="particle radius")
    fig0.colorbar(mappable=ScalarMappable(norm=Normalize(vmin=0, vmax=10), cmap='inferno'), ax=ax0)
    plt.show()

def plot_integ_sca_by_triangle(particle_size):
    res0 = compute_cross_sections_whole_by_triangle(REF_INDICES_RAW, WAVELENGTHS, particle_size)
    #res1 = compute_integrated_scattering_cross_section(0, 2 * math.pi, 0, math.pi * 0.57, REF_INDICES_RAW, WAVELENGTHS, particle_size)
    #res_exact = compute_cross_sections(REF_INDICES_RAW, WAVELENGTHS, particle_size)
    fig0 = plt.figure(num=0)
    ax0 = fig0.subplots(nrows=1, ncols=1)
    ax0.set_title("Scattering Cross Sections")
    ax0.plot(WAVELENGTHS, res0, label="Integrated")
    #ax0.plot(WAVELENGTHS, res1, label="Integrated")
    #ax0.plot(WAVELENGTHS, res_exact[0], label="Exact")
    ax0.set(xlabel="wavelength")
    ax0.legend()
    ax0.grid()
    plt.show()

def plot_integ_sca_by_triangle_file(particle_size, filename):
    data = load_selected_triangle("./res/sphere.ply")
    res0 = compute_cross_section_by_triangle(REF_INDICES_RAW, WAVELENGTHS, particle_size, data[0], data[1])
    #res1 = compute_integrated_scattering_cross_section(0, 2 * math.pi, 0, math.pi * 0.57, REF_INDICES_RAW, WAVELENGTHS, particle_size)
    #res_exact = compute_cross_sections(REF_INDICES_RAW, WAVELENGTHS, particle_size)
    fig0 = plt.figure(num=0)
    ax0 = fig0.subplots(nrows=1, ncols=1)
    ax0.set_title("Scattering Cross Sections")
    ax0.plot(WAVELENGTHS, res0, label="Integrated")
    #ax0.plot(WAVELENGTHS, res1, label="Integrated")
    #ax0.plot(WAVELENGTHS, res_exact[0], label="Exact")
    ax0.set(xlabel="wavelength")
    ax0.legend()
    ax0.grid()
    plt.show()

def plot_integ_sca_surface_by_triangle():
    PARTSIZE_LOWER = 50e-9
    PARTSIZE_UPPER = 100e-9
    partsizes = np.linspace(PARTSIZE_LOWER, PARTSIZE_UPPER, 50)
    scattering_cross_section = np.zeros((len(partsizes), len(WAVELENGTHS)))
    for i in range(len(partsizes)):
        scattering_cross_section[i] = compute_cross_sections_whole_by_triangle(REF_INDICES_RAW, WAVELENGTHS, partsizes[i])
    print("plotting...")
    fig0 = plt.figure(num=0)
    ax0 = fig0.subplots(nrows=1, ncols=1)
    ax0.set_title("Scattering Cross Section")
    ax0.contourf(WAVELENGTHS, partsizes, scattering_cross_section, cmap='inferno', levels=70)
    ax0.set(xlabel="wavelength", ylabel="particle radius")
    fig0.colorbar(mappable=ScalarMappable(norm=Normalize(vmin=0, vmax=10), cmap='inferno'), ax=ax0)
    plt.show()
    print("DONE!!!")

def plot_integ_sca_surface_by_triangle():
    PARTSIZE_LOWER = 50e-9
    PARTSIZE_UPPER = 100e-9
    partsizes = np.linspace(PARTSIZE_LOWER, PARTSIZE_UPPER, 50)
    scattering_cross_section = np.zeros((len(partsizes), len(WAVELENGTHS)))
    data = load_selected_triangle("./res/sphere.ply")
    for i in range(len(partsizes)):
        scattering_cross_section[i] = compute_cross_section_by_triangle(REF_INDICES_RAW, WAVELENGTHS, particle_size, data[0], data[1])
    print("plotting...")
    fig0 = plt.figure(num=0)
    ax0 = fig0.subplots(nrows=1, ncols=1)
    ax0.set_title("Scattering Cross Section")
    ax0.contourf(WAVELENGTHS, partsizes, scattering_cross_section, cmap='inferno', levels=70)
    ax0.set(xlabel="wavelength", ylabel="particle radius")
    fig0.colorbar(mappable=ScalarMappable(norm=Normalize(vmin=0, vmax=10), cmap='inferno'), ax=ax0)
    plt.show()
    print("DONE!!!")
####################### MAIN #######################

DIV = 500
REF_INDICES_RAW = load_ref_index("./res/refractive-index-silicon-2.csv")
WAVELENGTHS = np.linspace(REF_INDICES_RAW[0][0], REF_INDICES_RAW[-1][0], DIV)

#plot_surface_sca_ext()
#plot_coeff_sca_ext(80e-9)
#plot_sca_ext(80e-9)
#plot_ref_indices()
#plot_integ_sca(90e-9)
#plot_integ_sca_surface()
plot_integ_sca_by_triangle(80e-9)
#plot_integ_sca_surface_by_triangle()
#plot_integ_sca_by_triangle_file(80e-9, "./res/sphere.ply")

#x = np.linspace(0, 2 * math.pi, 100)
#pi2 = []
#pi3 = []
#pi4 = []
#tau1 = []
#tau2 = []
#tau3 = []
#for coord in x:
#    pitau = pi_tau_array(10, coord)
#    pi2.append(pitau[2][0])
#    pi3.append(pitau[3][0])
#    pi4.append(pitau[4][0])
#    tau1.append(pitau[1][1])
#    tau2.append(pitau[2][1])
#    tau3.append(pitau[3][1])
#plt.plot(x, pi2, label="pi2")
#plt.plot(x, pi3, label="pi3")
#plt.plot(x, pi4, label="pi4")
#plt.plot(x, tau1, label="tau1")
#plt.plot(x, tau2, label="tau2")
#plt.plot(x, tau3, label="tau3")
#plt.legend()
#plt.show()

#x = np.linspace(0.1, 30, 1000)
#sjn_0 = []
#sjn_1 = []
#sjn_2 = []
#for coord in x:
#    curr_sjn = sj_downward_reccurence(coord, 240)
#    sjn_0.append(curr_sjn[0])
#    sjn_1.append(curr_sjn[1])
#    sjn_2.append(curr_sjn[2])
#syn_x = []
#for i in range(ORDER_MAX):
    #sjn_x.append(special.spherical_jn(i, x))
    #syn_x.append(special.spherical_yn(i, x))
#plt.plot(x, sjn_0, label="j0")
#plt.plot(x, sjn_1, label="j1")
#plt.plot(x, sjn_2, label="j2")
#plt.legend()
#plt.xlim(left=0.0, right=30.0)
#plt.ylim(bottom=-1.0, top=1.0)

#plt.plot(x, syn_x[0], label="y0")
#plt.plot(x, syn_x[1], label="y1")
#plt.ylim(top=1.2, bottom=-1)
#plt.xlim(right=20)
#plt.xlim(left=0)
#plt.grid()
#plt.legend()
#plt.show()