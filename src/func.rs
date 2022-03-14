use scilib::math::complex::Complex;
use std::f64::consts::{PI};
use std::cmp::*;

pub fn sy_array<T: Into<Complex>, U: Into<usize>>(x: T, n: U) -> Vec<Complex> {
    let z: Complex = x.into();
    let order: usize = n.into();
    let y0 = -z.cos() / z;
    if order == 0 {
        return vec![y0];
    }
    let y1 = -z.cos() / (z * z) - z.sin() / z;
    if order == 1 {
        return vec![y0, y1];
    }
    let mut yn: Vec<Complex> = vec![Complex::from(0.0, 0.0); order + 1];
    yn[0] = y0;
    yn[1] = y1;
    for i in 2..(order + 1) {
        yn[i] = ((2 * i + 1) as f64) / z * yn[i - 1] - yn[i - 2];
    }
    return yn;
}

fn sj_upward_reccurence<T: Into<Complex>, U: Into<usize>>(z: T, n: U) -> Vec<Complex> {
    let order: usize = n.into();
    let count = order + 1;
    let x: Complex = z.into();
    let mut jn = vec![Complex::from(0.0, 0.0); count];
    jn[0] = x.sin() / x;
    jn[1] = x.sin() / x.powi(2) - x.cos() / x;
    for i in 1..=count - 2 {
        jn[i + 1] = (2 * i + 1) as f64 / x * jn[i] - jn[i - 1];
    }
    jn
}

fn sj_downward_reccurence<T: Into<Complex>, U: Into<usize>>(z: T, nl: U, nu: U) -> Vec<Complex> {
    let lower_order: usize = nl.into();
    let upper_order: usize = nu.into();
    let count = upper_order - lower_order + 1;
    let x: Complex = z.into();
    let mut jn= vec![Complex::from(0.0, 0.0); count + 2];
    jn[count + 1] = Complex::from(0.0, 0.0);
    jn[count] = Complex::from(1.0, 0.0);
    for i in (1..=count).rev() {
        jn[i - 1] = (2 * i + 1) as f64 / x * jn[i] - jn[i + 1];
    }
    jn.resize(count, Complex::from(0.0, 0.0));
    jn
}

pub fn sj_array<T: Into<Complex>, U: Into<usize>>(z: T, n: U) -> Vec<Complex> {
    let order: usize = n.into();
    let x: Complex = z.into();
    if x.modulus() > order as f64 / 2.0 {
        sj_upward_reccurence(x, order)
    } else {
        const PACK: usize = 50;
        let num_big_loop = order / PACK;
        let mut jn_all = Vec::<Vec<Complex>>::new();
        for i in 0..num_big_loop {
            jn_all.push(sj_downward_reccurence(x, i * PACK, (i + 1) * PACK));
        }
        let rest = order % PACK;
        if rest != 0 {
            jn_all.push(sj_downward_reccurence(x, order - rest, order));
        }

        let mut jn = Vec::<Complex>::with_capacity(order);
        let mut norm = x.sin() / x / jn_all[0][0];
        for i in 0..jn_all[0].len() {
            jn.push(jn_all[0][i] * norm);
        }
        for i in 1..jn_all.len() {
            norm = *jn.last().unwrap() / jn_all[i][0];
            for k in 1..jn_all[i].len() {
                jn.push(jn_all[i][k] * norm);
            }
        }
        jn
    }
}

pub fn pitau_array<U: Into<usize>>(theta: f64, n: U) -> (Vec<f64>, Vec<f64>) {
    let order: usize = n.into();
    assert!(order > 0);
    let mut pi_n = vec![0.0; order + 1];
    let mut tau_n = vec![0.0; order + 1];
    pi_n[0] = 0.0;
    pi_n[1] = 1.0;
    let nu = theta.cos();
    tau_n[1] = nu * pi_n[1] - 2.0 * pi_n[0];
    for i in 2..=order {
        let n = i as f64;
        pi_n[i] = (2.0 * n - 1.0) / (n - 1.0) * nu * pi_n[i - 1] - n / (n - 1.0) * pi_n[i - 2];
        tau_n[i] = n * nu * pi_n[i] - (n + 1.0) * pi_n[i - 1]
    }
    return (pi_n, tau_n);
}

pub fn psi_array<T: Into<Complex>>(x: T, sjn: &Vec<Complex>) -> Vec<Complex> {
    let mut psin: Vec<Complex> = Vec::with_capacity(sjn.len());
    let z: Complex = x.into();
    for j in sjn {
        psin.push(*j * z);
    }
    return psin;
}

pub fn xi_array<T: Into<Complex>>(x: T, jn: &Vec<Complex>, yn: &Vec<Complex>) -> Vec<Complex> {
    let z: Complex = x.into();
    let min = min(jn.len(), yn.len());
    let mut xinz: Vec<Complex> = Vec::with_capacity(min);
    for i in 0..min {
        xinz.push(z * (jn[i] + Complex::from(0.0, 1.0) * yn[i]));
    }
    return xinz;
}

pub fn xi_deriv_array<T: Into<Complex>>(x: T, xin: &Vec<Complex>) -> Vec<Complex> {
    let z: Complex = x.into();
    let mut xidn = Vec::<Complex>::with_capacity(xin.len() - 1);
    xidn.push(Complex::new());
    for i in 1..xidn.len() {
        let n: f64 = i as f64;
        xidn.push((n * xin[i - 1] - (n + 1.0) * xin[i + 1]) / (2.0 * n + 1.0));
    }
    xidn
}

pub fn d_array<T: Into<Complex>>(x: T, psin: &Vec<Complex>) -> Vec<Complex> {
    assert!(psin.len() > 1);
    let z: Complex = x.into();
    let mut dn = vec![Complex::from(0.0, 0.0); psin.len()];
    {
        let psi_last_n = psin.len() - 1;
        let psi_last = psin[psi_last_n];
        let psi_before_last = psin[psi_last_n - 1];
        dn[psin.len() - 1] = psi_before_last / psi_last - (psi_last_n as f64) / z;
    }
    for i in (1..dn.len()).rev() {
        let ndivz = (i as f64) / z;
        dn[i - 1] = ndivz - 1.0 / (dn[i] + ndivz);
    }
    return dn;
}

pub fn a_array<T: Into<Complex>, U: Into<usize>>(x: f64, m: T, n: U,  psin: &Vec<Complex>, xin: &Vec<Complex>, dn: &Vec<Complex>) -> Vec<Complex> {
    let order: usize = n.into();
    assert!(order < psin.len() && order < xin.len() && order < dn.len());
    let coeff_m: Complex = m.into();
    let mut an: Vec<Complex> = Vec::with_capacity(psin.len());
    an.push(Complex::from(0.0, 0.0));
    for j in 1..=order {
        let coeff_aj = dn[j] / coeff_m + (j as f64) / x;
        let aj: Complex = (coeff_aj * psin[j] - psin[j - 1]) / (coeff_aj * xin[j] - xin[j - 1]);
        an.push(aj);
    }
    return an;
}

pub fn b_array<T: Into<Complex>, U: Into<usize>>(x: f64, m: T, n: U,  psin: &Vec<Complex>, xin: &Vec<Complex>, dn: &Vec<Complex>) -> Vec<Complex> {
    let order: usize = n.into();
    assert!(order < psin.len() && order < xin.len() && order < dn.len());
    let coeff_m: Complex = m.into();
    let mut bn: Vec<Complex> = Vec::with_capacity(psin.len());
    bn.push(Complex::from(0.0, 0.0));
    for j in 1..=order {
        let coeff_bj = dn[j] * coeff_m + (j as f64) / x;
        let bj: Complex = (coeff_bj * psin[j] - psin[j - 1]) / (coeff_bj * xin[j] - xin[j - 1]);
        bn.push(bj);
    }
    return bn;
}