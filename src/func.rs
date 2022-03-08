use scilib::math::complex::Complex;
use std::f64::consts::{PI};

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
        yn[i] = Complex::from((2 * i + 1) as f64, 0.0) / z * yn[i - 1] - yn[i - 2];
    }
    return yn;
}

pub fn sj_array<T: Into<Complex>, U: Into<usize>>(x: T, n: U) -> Vec<Complex> {
    let z: Complex = x.into();
    let order: usize = n.into();
    let max_len = (order + 1) * 10; // ten times more to have the limit
    let mut jn: Vec<Complex> = vec![Complex::from(0.0, 0.0); max_len];
    jn[max_len - 1] = Complex::from(0.0, 0.0);
    jn[max_len - 2] = Complex::from(1.0, 0.0);
    for i in (1..=max_len - 2).rev() {
        jn[i - 1] = ((2 * i + 1) as f64) / z * jn[i] - jn[i + 1];
    }
    jn.resize(order + 1, Complex::from(0.0, 0.0));
    let coeff = z.sin() / (z * jn[0]);
    for i in 0..jn.len() {
        jn[i] *= coeff;
    }
    return jn;
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
    assert!(jn.len() == yn.len());
    let z: Complex = x.into();
    let mut xinz: Vec<Complex> = Vec::with_capacity(jn.len());
    for i in 0..jn.len() {
        xinz.push(z * (jn[i] + Complex::from(0.0, 1.0) * yn[i]));
    }
    return xinz;
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
    assert!(psin.len() == xin.len() && xin.len() == dn.len());
    let coeff_m: Complex = m.into();
    let order: usize = n.into();
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
    assert!(psin.len() == xin.len() && xin.len() == dn.len());
    let coeff_m: Complex = m.into();
    let order: usize = n.into();
    let mut bn: Vec<Complex> = Vec::with_capacity(psin.len());
    bn.push(Complex::from(0.0, 0.0));
    for j in 1..=order {
        let coeff_bj = dn[j] * coeff_m + (j as f64) / x;
        let bj: Complex = (coeff_bj * psin[j] - psin[j - 1]) / (coeff_bj * xin[j] - xin[j - 1]);
        bn.push(bj);
    }
    return bn;
}