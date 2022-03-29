use plotters::prelude::*;
use scilib::constant::*;
use scilib::math::complex::Complex;
use std::cmp::*;
use std::f64::consts::*;

mod csv;
mod func;
mod utils;

/*
use num_traits::{Num, Zero, One, NumCast};

pub trait Trigo {
    fn sin(&self) -> Self;
    fn cos(&self) -> Self;
    fn tan(&self) -> Self;
}

macro_rules! trigo_impl {
    ($t:ty) => {
        impl Trigo for $t {
            #[inline]
            fn sin(&self) -> $t {
                <$t>::sin(*self)
            }
            #[inline]
            fn cos(&self) -> $t {
                <$t>::cos(*self)
            }
            #[inline]
            fn tan(&self) -> $t {
                <$t>::tan(*self)
            }
        }
    };
}

trigo_impl!(f64);

impl Trigo for Complex {
    #[inline]
    fn sin(&self) -> Complex {
        Complex::sin(self)
    }
    #[inline]
    fn cos(&self) -> Complex {
        Complex::cos(self)
    }
    #[inline]
    fn tan(&self) -> Complex {
        Complex::tan(self)
    }
}*/

/*
xmodmap -e "keycode 105 = less greater less greater bar brokenbar lessthanequal greaterthanequal"
 */
fn y_n(order: usize, x: f64) -> Vec<f64> {
    /*let mut yn: Vec<f64> = Vec::with_capacity(max(order + 1, 2));
    yn.push(-x.cos() / x); // y0
    yn.push(-x.cos() / (x * x) - x.sin() / x); // y1
    for i in 2..(order+1) {
        yn.push((2 * i + 1) as f64 / x * yn.last().unwrap() - yn[yn.len() - 2]);
    }*/
    let y_0 = -x.cos() / x;
    if order == 0 {
        return vec![y_0];
    }
    let y_1 = -x.cos() / (x * x) - x.sin() / x;
    if order == 1 {
        return vec![y_1];
    }
    let mut yn: Vec<f64> = vec![0.0; order + 1];
    yn[0] = y_0;
    yn[1] = y_1;
    for i in 2..(order + 1) {
        yn[i] = (2 * i + 1) as f64 / x * yn[i - 1] - yn[i - 2];
    }
    return yn;
}

fn y_n_z(order: usize, z: Complex) -> Vec<Complex> {
    let y_0 = -z.cos() / z;
    if order == 0 {
        return vec![y_0];
    }
    let y_1 = -z.cos() / (z * z) - z.sin() / z;
    if order == 1 {
        return vec![y_1];
    }
    let mut yn: Vec<Complex> = vec![Complex::from(0.0, 0.0); order + 1];
    yn[0] = y_0;
    yn[1] = y_1;
    for i in 2..(order + 1) {
        yn[i] = Complex::from((2 * i + 1) as f64, 0.0) / z * yn[i - 1] - yn[i - 2];
    }
    return yn;
}

/*
fn main() -> Result<(), Box<dyn std::error::Error>> {

    let mut y0: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut y1: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut y2: Vec<(f64, f64)> = Vec::with_capacity(400);
    for i in 0..399 {
        let coord_x = (i + 1) as f64 / 20.0;
        let y_n = y_n(2, coord_x);
        y0.push((coord_x, y_n[0]));
        y1.push((coord_x, y_n[1]));
        y2.push((coord_x, y_n[2]));
    }

    plot_png("./results/bessel-yn.png", (800, 800), "bessel yn", (0.0, 20.0), (-1.0, 0.5), &vec![y0, y1, y2], &vec!["y0", "y1", "y2"], &vec![RED, BLUE, GREEN])?;

    Ok(())
}*/

fn j_n(order: usize, x: f64) -> Vec<f64> {
    let max_len = (order + 1) * 10; // ten times more to have the limit
    let mut jn: Vec<f64> = vec![0.0; max_len];
    jn[max_len - 1] = 0.0;
    jn[max_len - 2] = 1.0;
    for i in (1..=max_len - 2).rev() {
        jn[i - 1] = (2 * i + 1) as f64 / x * jn[i] - jn[i + 1];
    }
    jn.resize(order + 1, 0.0);
    let coeff = x.sin() / (x * jn[0]);
    for i in 0..jn.len() {
        jn[i] *= coeff;
    }
    return jn;
}

fn j_n_z(order: usize, z: Complex) -> Vec<Complex> {
    let max_len = order + 1; //(order + 1) * 5; // ten times more to have the limit
    let mut jn: Vec<Complex> = vec![Complex::from(0.0, 0.0); max_len];
    jn[0] = z.sin() / z;
    jn[1] = z.sin() / z.powi(2) - z.cos() / z;
    for i in 1..=max_len - 2 {
        jn[i + 1] = Complex::from((2 * i + 1) as f64, 0.0) / z * jn[i] - jn[i - 1];
        if jn[i + 1].re.abs() > 1e10 {
            jn[i + 1].re = 0.4;
        }
        //println!("{} /// {} /// {}", z.re, i, jn[i + 1]);
    }
    /*jn[max_len - 1] = Complex::from(0.0, 0.0);
    jn[max_len - 2] = Complex::from(1.0, 0.0);
    for i in (1..=max_len - 2).rev() {
        jn[i - 1] = Complex::from((2 * i + 1) as f64, 0.0) / z * jn[i] - jn[i + 1];
        //println!("{} /// {} /// {}", z.re, i, jn[i - 1]);
    }
    //println!("");
    jn.resize(order + 1, Complex::from(0.0, 0.0));
    let coeff = z.sin() / (z * jn[0]);
    for i in 0..jn.len() {
        jn[i] *= coeff;
        //println!("{} /// {} /// {}", z.re, i, jn[i]);
    }*/
    //println!("");
    return jn;
}

fn upward_reccurence<T: Into<Complex>>(z: T, n: usize) -> Vec<Complex> {
    let count = n + 1;
    let x: Complex = z.into();
    let mut jn = vec![Complex::from(0.0, 0.0); count];
    jn[0] = x.sin() / x;
    jn[1] = x.sin() / x.powi(2) - x.cos() / x;
    for i in 1..=count - 2 {
        jn[i + 1] = (2 * i + 1) as f64 / x * jn[i] - jn[i - 1];
    }
    jn
}

fn downward_reccurence<T: Into<Complex>>(z: T, nl: usize, nu: usize) -> Vec<Complex> {
    let count = nu - nl + 1;
    let x: Complex = z.into();
    let mut jn = vec![Complex::from(0.0, 0.0); count + 2];
    jn[count + 1] = Complex::from(0.0, 0.0);
    jn[count] = Complex::from(1.0, 0.0);
    for i in (1..=count).rev() {
        jn[i - 1] = (2 * i + 1) as f64 / x * jn[i] - jn[i + 1];
    }
    jn.resize(count, Complex::from(0.0, 0.0));
    jn
}

const LOOP: usize = 50;

fn fucking_jn<T: Into<Complex>>(z: T, n: usize) -> Vec<Complex> {
    let x: Complex = z.into();
    if x.modulus() > n as f64 / 2.0 {
        upward_reccurence(x, n)
    } else {
        let num_big_loop = n / LOOP;
        let mut jn_all = Vec::<Vec<Complex>>::new();
        for i in 0..num_big_loop {
            jn_all.push(downward_reccurence(x, i * LOOP, (i + 1) * LOOP));
        }
        let rest = n % LOOP;
        if rest != 0 {
            jn_all.push(downward_reccurence(x, n - rest, n));
        }

        let mut jn = Vec::<Complex>::with_capacity(n);
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

fn miller(x: f64, order: usize) -> Vec<f64> {
    let mut jn: Vec<f64> = vec![0.0; order + 2];
    jn[order + 1] = 0.0;
    jn[order] = 1.0;
    for i in (1..order + 1).rev() {
        jn[i - 1] = (2 * i + 1) as f64 / x * jn[i] - jn[i + 1];
    }

    let norm = x.sin() / x / jn[0];
    for i in 0..jn.len() {
        jn[i] *= norm;
        println!("{} /// {}", i, jn[i]);
    }
    return jn;
}

/*
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut j0: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut j1: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut j2: Vec<(f64, f64)> = Vec::with_capacity(400);
    for i in 0..399 {
        let coord_x = (i + 1) as f64 / 10.0;
        let jn_x = fucking_jn(Complex::from(coord_x, 0.0), 500usize); //j_n_z(100usize, Complex::from(coord_x, 0.0));
        j0.push((coord_x, jn_x[0].re));
        j1.push((coord_x, jn_x[1].re));
        j2.push((coord_x, jn_x[2].re));
    }

    plot_png(
        "./results/bessel-jn-z.png",
        (800, 800),
        "bessel jn",
        (0.0, 40.0),
        (-0.5, 1.2),
        &vec![j0, j1, j2],
        &vec!["j0", "j1", "j2"],
        &vec![RED, BLUE, GREEN],
    )?;

    Ok(())
}*/

fn pitau_n(order: usize, theta: f64) -> (Vec<f64>, Vec<f64>) {
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

fn psi_n(jn: &Vec<f64>, x: f64) -> Vec<f64> {
    let mut psin: Vec<f64> = Vec::with_capacity(jn.len());
    for j in jn {
        psin.push(j * x);
    }
    return psin;
}

fn psi_n_z(jnz: &Vec<Complex>, z: Complex) -> Vec<Complex> {
    let mut psinz: Vec<Complex> = Vec::with_capacity(jnz.len());
    for j in jnz {
        psinz.push(*j * z);
    }
    return psinz;
}

fn d_n(psin: &Vec<f64>, x: f64) -> Vec<f64> {
    assert!(psin.len() > 1);
    let mut dn = vec![0.0; psin.len()];
    {
        let psi_last_n = psin.len() - 1;
        let psi_last = psin[psi_last_n];
        let psi_before_last = psin[psi_last_n - 1];
        dn[psin.len() - 1] = psi_before_last / psi_last - psi_last_n as f64 / x;
    }
    for i in (1..dn.len()).rev() {
        let ndivx = i as f64 / x;
        dn[i - 1] = ndivx - 1.0 / (dn[i] + ndivx);
    }
    return dn;
}

fn d_n_z(psinz: &Vec<Complex>, z: Complex) -> Vec<Complex> {
    assert!(psinz.len() > 1);
    let mut dnz = vec![Complex::from(0.0, 0.0); psinz.len()];
    {
        let psi_last_n = psinz.len() - 1;
        let psi_last = psinz[psi_last_n];
        let psi_before_last = psinz[psi_last_n - 1];
        dnz[psinz.len() - 1] =
            psi_before_last / psi_last - Complex::from(psi_last_n as f64, 0.0) / z;
    }
    for i in (1..dnz.len()).rev() {
        let ndivz = Complex::from(i as f64, 0.0) / z;
        dnz[i - 1] = ndivz - Complex::from(1.0, 0.0) / (dnz[i] + ndivz);
    }
    return dnz;
}

/*
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut psi0: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut psi1: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut psi2: Vec<(f64, f64)> = Vec::with_capacity(400);
    for i in 0..399 {
        let coord_x = (i + 1) as f64 / 20.0;
        let jn = j_n(2, coord_x);
        let psin = psi_n(&jn, coord_x);
        psi0.push((coord_x, psin[0]));
        psi1.push((coord_x, psin[1]));
        psi2.push((coord_x, psin[2]));
    }

    plot_png("./results/riccati-bessel-psin.png", (800, 800), "riccati-bessel psin", (0.0, 20.0), (-2.0, 2.0), &vec![psi0, psi1, psi2], &vec!["psi0", "psi1", "psi2"], &vec![RED, BLUE, GREEN])?;

    Ok(())
}*/

/*
fn main() -> Result<(), Box<dyn std::error::Error>> {

    let mut d0: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut d1: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut d2: Vec<(f64, f64)> = Vec::with_capacity(400);
    for i in 0..399 {
        let coord_x = (i + 1) as f64 / 20.0;
        let jn = j_n(2, coord_x);
        let psin = psi_n(&jn, coord_x);
        let dn = d_n(&psin, coord_x);
        d0.push((coord_x, dn[0]));
        d1.push((coord_x, dn[1]));
        d2.push((coord_x, dn[2]));
    }

    plot_png("./results/log-deriv-dn.png", (800, 800), "log derivative dn", (0.0, 20.0), (-10.0, 10.0), &vec![d0, d1, d2], &vec!["d0", "d1", "d2"], &vec![RED, BLUE, GREEN])?;

    Ok(())
}*/

fn xi_n(jn: &Vec<f64>, yn: &Vec<f64>, x: f64) -> Vec<Complex> {
    assert!(jn.len() == yn.len());
    let mut xin: Vec<Complex> = Vec::with_capacity(jn.len());
    for i in 0..jn.len() {
        xin.push(Complex::from(jn[i], yn[i]) * x);
    }
    return xin;
}

fn xi_n_z(jnz: &Vec<Complex>, ynz: &Vec<Complex>, z: Complex) -> Vec<Complex> {
    assert!(jnz.len() == ynz.len());
    let mut xinz: Vec<Complex> = Vec::with_capacity(jnz.len());
    for i in 0..jnz.len() {
        xinz.push(z * (jnz[i] + Complex::from(0.0, 1.0) * ynz[i]));
    }
    return xinz;
}

/*
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut xid0im: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut xid1im: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut xid2im: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut xid0re: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut xid1re: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut xid2re: Vec<(f64, f64)> = Vec::with_capacity(400);
    for i in 0..399 {
        let coord_x = (i + 1) as f64 / 20.0;
        let jn = func::sj_array(coord_x, 2usize);
        let yn = func::sy_array(coord_x, 2usize);
        let xin = func::xi_array(coord_x, &jn, &yn);
        let xidn = func::xi_deriv_array(coord_x, &xin);
        xid0im.push((coord_x, xidn[0].im));
        xid1im.push((coord_x, xidn[1].im));
        xid2im.push((coord_x, xidn[2].im));
        xid0re.push((coord_x, xidn[0].re));
        xid1re.push((coord_x, xidn[1].re));
        xid2re.push((coord_x, xidn[2].re));
    }

    plot_png(
        "./results/riccati-bessel-xin-imaginary.png",
        (800, 800),
        "riccati-bessel xin (i)",
        (0.0, 20.0),
        (-10.0, 10.0),
        &vec![xid0im, xid1im, xid2im],
        &vec!["xi0im", "xi1im", "xi2im"],
        &vec![RED, BLUE, GREEN],
    )?;

    plot_png(
        "./results/riccati-bessel-xin-real.png",
        (800, 800),
        "riccati-bessel xin (r)",
        (0.0, 20.0),
        (-10.0, 10.0),
        &vec![xid0re, xid1re, xid2re],
        &vec!["xi0re", "xi1re", "xi2re"],
        &vec![RED, BLUE, GREEN],
    )?;

    Ok(())
}*/

/*
fn main() -> Result<(), Box<dyn std::error::Error>> {
    {
        let mut pi2: Vec<(f64, f64)> = Vec::with_capacity(400);
        let mut pi3: Vec<(f64, f64)> = Vec::with_capacity(400);
        let mut pi4: Vec<(f64, f64)> = Vec::with_capacity(400);
        for i in 0..399 {
            let coord_x = (i + 1) as f64 / 20.0;
            let pin = pitau_n(4, coord_x).0;
            pi2.push((coord_x, pin[2]));
            pi3.push((coord_x, pin[3]));
            pi4.push((coord_x, pin[4]));
        }

        plot_png("./results/angle-dependent-pin.png", (800, 800), "angle dependent pin", (0.0, 7.0), (-11.0, 11.0), &vec![pi2, pi3, pi4], &vec!["pi2", "pi3", "pi4"], &vec![RED, BLUE, GREEN])?;
    }
    {
        let mut tau1: Vec<(f64, f64)> = Vec::with_capacity(400);
        let mut tau2: Vec<(f64, f64)> = Vec::with_capacity(400);
        let mut tau3: Vec<(f64, f64)> = Vec::with_capacity(400);
        for i in 0..399 {
            let coord_x = (i + 1) as f64 / 20.0;
            let taun = pitau_n(3, coord_x).1;
            tau1.push((coord_x, taun[1]));
            tau2.push((coord_x, taun[2]));
            tau3.push((coord_x, taun[3]));
        }

        plot_png("./results/angle-dependent-taun.png", (800, 800), "angle dependent taun", (0.0, 7.0), (-7.0, 7.0), &vec![tau1, tau2, tau3], &vec!["tau1", "tau2", "tau3"], &vec![RED, BLUE, GREEN])?;
    }
    Ok(())
}*/

fn plot_png(
    path: &str,
    size: (u32, u32),
    title: &str,
    boundaries_x: (f64, f64),
    boundaries_y: (f64, f64),
    functions: &Vec<Vec<(f64, f64)>>,
    labels: &Vec<&str>,
    colors: &Vec<RGBColor>,
) -> Result<(), Box<dyn std::error::Error>> {
    assert!(functions.len() == labels.len() && labels.len() == colors.len());
    let root = BitMapBackend::new(path, size).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(
            boundaries_x.0..boundaries_x.1,
            boundaries_y.0..boundaries_y.1,
        )?;

    chart.configure_mesh().draw()?;
    println!("pushing data...");
    for i in 0..functions.len() {
        let closure =
            move |(x, y): (i32, i32)| PathElement::new(vec![(x, y), (x + 20, y)], &colors[i]);
        chart
            .draw_series(LineSeries::new(functions[i].clone(), &colors[i]))?
            .label(labels[i])
            .legend(closure);
    }
    println!("drawing...");
    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}

fn a_n(
    order: usize,
    psin: &Vec<f64>,
    xin: &Vec<Complex>,
    dn: &Vec<f64>,
    m: f64,
    x: f64,
) -> Vec<Complex> {
    assert!(psin.len() == xin.len() && xin.len() == dn.len());
    let mut an: Vec<Complex> = Vec::with_capacity(psin.len());
    an.push(Complex::from(0.0, 0.0));
    for j in 1..=order {
        let n = j as f64;
        let coeff_aj = dn[j] / m + n / x;
        let aj: Complex =
            Complex::from(coeff_aj * psin[j] - psin[j - 1], 0) / (coeff_aj * xin[j] - xin[j - 1]);
        an.push(aj);
    }
    return an;
}

fn b_n(
    order: usize,
    psin: &Vec<f64>,
    xin: &Vec<Complex>,
    dn: &Vec<f64>,
    m: f64,
    x: f64,
) -> Vec<Complex> {
    assert!(psin.len() == xin.len() && xin.len() == dn.len());
    let mut bn: Vec<Complex> = Vec::with_capacity(psin.len());
    bn.push(Complex::from(0.0, 0.0));
    for j in 1..=order {
        let n = j as f64;
        let coeff_bj = dn[j] * m + n / x;
        let bj: Complex =
            Complex::from(coeff_bj * psin[j] - psin[j - 1], 0) / (coeff_bj * xin[j] - xin[j - 1]);
        bn.push(bj);
    }
    return bn;
}

fn a_n_z(
    z: f64,
    m: Complex,
    order: usize,
    psin: &Vec<f64>,
    xin: &Vec<Complex>,
    dnz: &Vec<Complex>,
) -> Vec<Complex> {
    assert!(psin.len() == xin.len() && xin.len() == dnz.len());
    let mut anz: Vec<Complex> = Vec::with_capacity(psin.len());
    anz.push(Complex::from(0.0, 0.0));
    for j in 1..=order {
        let n = Complex::from(j as f64, 0.0);
        let coeff_aj = dnz[j] / m + n / z;
        let aj: Complex = (coeff_aj * psin[j] - psin[j - 1]) / (coeff_aj * xin[j] - xin[j - 1]);
        anz.push(aj);
    }
    return anz;
}

fn b_n_z(
    z: f64,
    m: Complex,
    order: usize,
    psin: &Vec<f64>,
    xin: &Vec<Complex>,
    dnz: &Vec<Complex>,
) -> Vec<Complex> {
    assert!(psin.len() == xin.len() && xin.len() == dnz.len());
    let mut bnz: Vec<Complex> = Vec::with_capacity(psin.len());
    bnz.push(Complex::from(0.0, 0.0));
    for j in 1..=order {
        let n = Complex::from(j as f64, 0.0);
        let coeff_bj = dnz[j] * m + n / z;
        let bj: Complex = (coeff_bj * psin[j] - psin[j - 1]) / (coeff_bj * xin[j] - xin[j - 1]);
        bnz.push(bj);
    }
    return bnz;
}

/*
fn main() -> Result<(), Box<dyn std::error::Error>> {
    //let upper_x = 2.0 * std::f64::consts::PI * MED_N * PART_SIZE;
    let mut an_r: Vec<Vec<(f64, f64)>> = vec![Vec::with_capacity(500); 4];
    let mut an_i: Vec<Vec<(f64, f64)>> = vec![Vec::with_capacity(500); 4];
    let mut bn_r: Vec<Vec<(f64, f64)>> = vec![Vec::with_capacity(500); 4];
    let mut bn_i: Vec<Vec<(f64, f64)>> = vec![Vec::with_capacity(500); 4];

    let wavelenght_boundaries = (400e-9, 600e-9);
    let medium_n = 1.0;
    let m = 3.5 / medium_n; // 1.0 : just to remember the refractive index of the air
    let upper_x = 2.0 * std::f64::consts::PI * medium_n * 200e-9;
    let step = (wavelenght_boundaries.1 - wavelenght_boundaries.0) / 500.0;//((upper_limit - lower_limit) / 500.0).abs();

    for i in 0..=499 {
        let current_wavelenght = wavelenght_boundaries.0 + step * (i + 1) as f64;
        let coord_x = upper_x / current_wavelenght;

        let jn_x = j_n(3, coord_x);
        let yn_x = y_n(3, coord_x);
        let psin_x = psi_n(&jn_x, coord_x);
        let xin_x = xi_n(&jn_x, &yn_x, coord_x);

        let jn_mx = j_n(3, coord_x * m);
        let psin_mx = psi_n(&jn_mx, coord_x * m);
        let dn_mx = d_n(&psin_mx, coord_x * m);

        let an = a_n(3, &psin_x, &xin_x, &dn_mx, m, coord_x);
        let bn = b_n(3, &psin_x, &xin_x, &dn_mx, m, coord_x);

        for j in 1..4 {
            an_r[j].push((current_wavelenght, an[j].re));
            an_i[j].push((current_wavelenght, an[j].im));
            bn_r[j].push((current_wavelenght, bn[j].re));
            bn_i[j].push((current_wavelenght, bn[j].im));
        }
    }

    println!("calculating...");
    plot_png(
        "./results/scattering-coeff-anr.png",
        (2000, 400),
        "scattering coeff anr",
        (wavelenght_boundaries.0, wavelenght_boundaries.1),
        (-0.2, 1.2),
        &an_r[1..].to_vec(),
        &vec!["anr1", "anr2", "anr3"],
        &vec![RED, BLUE, GREEN],
    )?;

    plot_png(
        "./results/scattering-coeff-ani.png",
        (2000, 400),
        "scattering coeff ani",
        (wavelenght_boundaries.0, wavelenght_boundaries.1),
        (-1.0, 1.0),
        &an_i[1..].to_vec(),
        &vec!["ani1", "ani2", "ani3"],
        &vec![RED, BLUE, GREEN],
    )?;

    plot_png(
        "./results/scattering-coeff-bnr.png",
        (2000, 400),
        "scattering coeff bnr",
        (wavelenght_boundaries.0, wavelenght_boundaries.1),
        (-0.2, 1.2),
        &bn_r[1..].to_vec(),
        &vec!["bnr1", "bnr2", "bnr3"],
        &vec![RED, BLUE, GREEN],
    )?;

    plot_png(
        "./results/scattering-coeff-bni.png",
        (2000, 400),
        "scattering coeff bni",
        (wavelenght_boundaries.0, wavelenght_boundaries.1),
        (-1.0, 1.0),
        &bn_i[1..].to_vec(),
        &vec!["bni1", "bni2", "bni3"],
        &vec![RED, BLUE, GREEN],
    )?;

    Ok(())
}*/

/*
fn main() -> Result<(), Box<dyn std::error::Error>> {

    let mut c_sca: Vec<(f64, f64)> = Vec::with_capacity(500);
    let mut c_ext: Vec<(f64, f64)> = Vec::with_capacity(500);

    let wavelenght_boundaries = (400e-9, 1800e-9);
    let medium_n = 1.0;
    let m = Complex::from(3.5, 0.05) / medium_n; // 1.0 : just to remember the refractive index of the air
    let upper_x = 2.0 * std::f64::consts::PI * medium_n * 200e-9;
    let step = (wavelenght_boundaries.1 - wavelenght_boundaries.0) / 500.0;//((upper_limit - lower_limit) / 500.0).abs();

    for i in 0..=499 {
        let current_wavelenght = wavelenght_boundaries.0 + step * (i + 1) as f64;
        let coord_x = upper_x / current_wavelenght;

        let jn_x = func::sj_array(coord_x, 3usize);
        let yn_x = func::sy_array(coord_x, 3usize);
        let psin_x = func::psi_array(coord_x, &jn_x);
        let xin_x = func::xi_array(coord_x, &jn_x, &yn_x);

        let jn_mx = func::sj_array(coord_x * m, 3usize);
        let psin_mx = func::psi_array(coord_x * m, &jn_mx);
        let dn_mx = func::d_array(coord_x * m, &psin_mx);

        let an = func::a_array(coord_x, m, 3usize, &psin_x, &xin_x, &dn_mx);
        let bn = func::b_array(coord_x, m, 3usize, &psin_x, &xin_x, &dn_mx);

        let mut sum_sca = 0.0;
        let mut sum_ext = 0.0;
        for j in 1..4 {
            sum_sca += (2.0 * j as f64 + 1.0) * (an[j].re.powi(2) + an[j].im.powi(2) + bn[j].re.powi(2) + bn[j].im.powi(2));
            sum_ext += (2.0 * j as f64 + 1.0) * (an[j].re + bn[j].re);
        }
        let surface = (200.0e-9) * (200.0e-9) * std::f64::consts::PI;
        let mul = (current_wavelenght / medium_n).powi(2) / (2.0 * std::f64::consts::PI) / surface;
        sum_sca *= mul;
        sum_ext *= mul;
        c_sca.push((current_wavelenght, sum_sca));
        c_ext.push((current_wavelenght, sum_ext));
    }

    println!("calculating...");
    plot_png(
        "./results/cross-section-test.png",
        (2000, 400),
        "cross section c_sca & c_ext",
        (wavelenght_boundaries.0, wavelenght_boundaries.1),
        (0.0, 10.0),
        &vec![c_sca, c_ext],
        &vec!["c_sca", "c_ext"],
        &vec![RED, BLUE],
    )?;

    Ok(())
}*/

/*
fn main() -> Result<(), Box<dyn std::error::Error>> {

    let mut c_sca_coeff: Vec<Vec<(f64, f64)>> = vec![Vec::with_capacity(500); 3];
    let mut c_ext_coeff: Vec<Vec<(f64, f64)>> = vec![Vec::with_capacity(500); 3];

    let wavelenght_boundaries = (400e-9, 1800e-9);
    let medium_n = 1.0;
    let m = Complex::from(3.5, 0.05) / medium_n; // 1.0 : just to remember the refractive index of the air
    let upper_x = 2.0 * std::f64::consts::PI * medium_n * 200e-9;
    let step = (wavelenght_boundaries.1 - wavelenght_boundaries.0) / 500.0;

    for i in 0..=499 {
        let current_wavelenght = wavelenght_boundaries.0 + step * (i + 1) as f64;
        let coord_x = upper_x / current_wavelenght;

        let jn_x = func::sj_array(coord_x, 3usize);
        let yn_x = func::sy_array(coord_x, 3usize);
        let psin_x = func::psi_array(coord_x, &jn_x);
        let xin_x = func::xi_array(coord_x, &jn_x, &yn_x);

        let jn_mx = func::sj_array(coord_x * m, 3usize);
        let psin_mx = func::psi_array(coord_x * m, &jn_mx);
        let dn_mx = func::d_array(coord_x * m, &psin_mx);

        let an = func::a_array(coord_x, m, 3usize, &psin_x, &xin_x, &dn_mx);
        let bn = func::b_array(coord_x, m, 3usize, &psin_x, &xin_x, &dn_mx);

        let mul = (current_wavelenght / medium_n).powi(2) / (2.0 * std::f64::consts::PI) / ((200.0e-9) * (200.0e-9) * std::f64::consts::PI);
        for j in 1..4 {
            c_sca_coeff[j - 1].push((current_wavelenght, mul * (2.0 * j as f64 + 1.0) * (an[j].re.powi(2) + an[j].im.powi(2) + bn[j].re.powi(2) + bn[j].im.powi(2))));
            c_ext_coeff[j - 1].push((current_wavelenght, mul * (2.0 * j as f64 + 1.0) * (an[j].re + bn[j].re)));
        }
    }

    println!("calculating...");
    plot_png(
        "./results/c-sca-coeffs.png",
        (2000, 400),
        "c_sca coeffs",
        (wavelenght_boundaries.0, wavelenght_boundaries.1),
        (0.0, 10.0),
        &c_sca_coeff,
        &vec!["coeff 1", "coeff 2", "coeff 3"],
        &vec![RED, BLUE, GREEN],
    )?;
    plot_png(
        "./results/c-ext-coeffs.png",
        (2000, 400),
        "c_ext coeffs",
        (wavelenght_boundaries.0, wavelenght_boundaries.1),
        (0.0, 10.0),
        &c_ext_coeff,
        &vec!["coeff 1", "coeff 2", "coeff 3"],
        &vec![RED, BLUE, GREEN],
    )?;

    Ok(())
}*/

trait Cast {
    fn to_usize(&self) -> usize;
}

impl Cast for f64 {
    fn to_usize(&self) -> usize {
        *self as usize
    }
}

use std::ops::{Add, Sub};

fn approx_eq<T: PartialOrd + Add<Output = T> + Sub<Output = T> + Copy>(
    exact: T,
    value: T,
    approx: T,
) -> bool {
    if (exact - approx) < value && (exact + approx) > value {
        return true;
    }
    false
}

fn get_ref_index(csv_file: &Vec<(f64, f64, f64)>, wavelength: f64) -> Complex {
    for i in 0..csv_file.len() - 1 {
        if approx_eq(csv_file[i].0, wavelength, 1e-3) {
            return Complex::from(csv_file[i].1, csv_file[i].2);
        } else if csv_file[i].0 < wavelength && csv_file[i + 1].0 > wavelength {
            let point_re = utils::linear_interpolation_2D(
                (csv_file[i].0, csv_file[i].1),
                (csv_file[i + 1].0, csv_file[i + 1].1),
                wavelength,
            );
            let point_im = utils::linear_interpolation_2D(
                (csv_file[i].0, csv_file[i].2),
                (csv_file[i + 1].0, csv_file[i + 1].2),
                wavelength,
            );
            return Complex::from(point_re.1, point_im.1);
        }
    }
    Complex::from(csv_file.last().unwrap().1, csv_file.last().unwrap().2)
}

#[macro_export]
macro_rules! summation {
    ($type:ty, $var:pat, $range:expr, $form:expr) => {{
        let mut result: $type = 0 as $type;
        for $var in $range {
            result += $form;
        }
        result
    }};
}

#[macro_export]
macro_rules! product {
    ($type:ty, $var:pat, $range:expr, $form:expr) => {{
        let mut result: $type = 1 as $type;
        for $var in $range {
            result *= $form;
        }
        result
    }};
}

mod integration {

    fn chunk_trapez(a: f64, b: f64, offset: f64) -> f64 {
        (a + b) * offset
    }

    pub fn trapez_fn(
        function: impl Fn(f64) -> f64,
        lower_bound: f64,
        upper_bound: f64,
        div: usize,
    ) -> f64 {
        assert!(div >= 1);
        let mut area = 0.0;
        let offset = upper_bound - lower_bound;
        let mut a = function(lower_bound);
        let step = offset / div as f64;
        for i in 0..div {
            let b = function(step * (i as f64 + 1.0));
            area += chunk_trapez(a, b, step);
            a = b;
        }
        area
    }

    pub fn trapez_dt(data: &Vec<(f64, f64)>) -> f64 {
        let mut area = 0.0;
        for i in 0..data.len() - 1 {
            area += chunk_trapez(data[i].1, data[i + 1].1, data[i + 1].0 - data[i].0);
        }
        area
    }
}

fn get_en(e0: f64, n: i32) -> Complex {
    let nf = n as f64;
    Complex::i().powi(n as i32) * e0 * (2.0 * nf + 1.0) / (nf * (nf + 1.0))
}

fn get_comp(
    theta: f64,
    e0: f64,
    last_mode: usize,
    rho: f64,
    x: f64,
    m: Complex,
) -> (Complex, Complex, Complex, Complex) {
    let jn_x = func::sj_array(x, last_mode as usize);
    let yn_x = func::sy_array(x, last_mode as usize);
    let pitaun = func::pitau_array(theta, last_mode as usize);

    let psin_x = func::psi_array(x, &jn_x);
    let xin_x = func::xi_array(x, &jn_x, &yn_x);

    let jn_rho = func::sj_array(rho, last_mode as usize);
    let yn_rho = func::sy_array(rho, last_mode as usize);
    let xin_rho = func::xi_array(rho, &jn_rho, &yn_rho);
    let xidn_rho = func::xi_deriv_array(rho, &xin_rho);

    let jn_mx = func::sj_array(x * m, last_mode as usize);
    let psin_mx = func::psi_array(x * m, &jn_mx);
    let dn_mx = func::d_array(x * m, &psin_mx);

    let an = func::a_array(x, m, last_mode as usize, &psin_x, &xin_x, &dn_mx);
    let bn = func::b_array(x, m, last_mode as usize, &psin_x, &xin_x, &dn_mx);

    let mut sum_es_theta = Complex::new();
    let mut sum_hs_theta = Complex::new();
    let mut sum_es_phi = Complex::new();
    let mut sum_hs_phi = Complex::new();
    for n in 1..=last_mode {
        let en = get_en(e0, n as i32);
        sum_es_theta += en
            * (Complex::i() * an[n] * xidn_rho[n] * pitaun.1[n] - bn[n] * xin_rho[n] * pitaun.0[n]);
        sum_hs_theta += en
            * (Complex::i() * bn[n] * xidn_rho[n] * pitaun.1[n] - an[n] * xin_rho[n] * pitaun.0[n]);
        sum_es_phi += en
            * (bn[n] * xin_rho[n] * pitaun.1[n] - Complex::i() * an[n] * xidn_rho[n] * pitaun.0[n]);
        sum_hs_phi += en
            * (Complex::i() * bn[n] * xidn_rho[n] * pitaun.0[n] - an[n] * xin_rho[n] * pitaun.1[n]);
    }
    (sum_es_theta, sum_hs_theta, sum_es_phi, sum_hs_phi)
}

fn integ_whole_particle(
    wavelength: f64,
    medium_n: f64,
    medium_mu: f64,
    ref_index: Complex,
    particle_size: f64,
    e0: f64,
    last_mode: usize,
    div: usize,
) -> f64 {
    assert!(div >= 1);
    let mul = wavelength.powi(2) / (8.0 * PI * C * medium_mu * std::f64::consts::PI * particle_size.powi(2));
    let x = 2.0 * PI * medium_n * particle_size / wavelength;
    let m = ref_index / medium_n;
    let step = PI / div as f64;

    let mut data = Vec::<(f64, f64)>::with_capacity(div + 1);
    for i in 0..=div {
        let theta = step * i as f64;
        let res = get_comp(theta, e0, last_mode, x, x, m);
        data.push((
            theta,
            (res.0 * res.3.conjugate()).re * theta.sin()
                - (res.2 * res.1.conjugate()).re * theta.sin(),
        ));
        //println!("{} {}", (res.0 * res.3.conjugate()).re * theta.sin(), (res.2 * res.1.conjugate()).re * theta.sin());
    }
    let integ = integration::trapez_dt(&data);
    mul * integ * 300.0
}


fn main() -> Result<(), Box<dyn std::error::Error>> {

    let mut c_sca: Vec<(f64, f64)> = Vec::with_capacity(500);
    let mut c_ext: Vec<(f64, f64)> = Vec::with_capacity(500);

    let ref_indices = csv::parse("./res/refractive-index-silicon.csv");
    let wavelength_boundaries = (
        ref_indices[0].0 * 1e-6,
        ref_indices.last().unwrap().0 * 1e-6,
    );
    let medium_n = 1.0;
    let particle_size = 85e-9;
    let upper_x = 2.0 * std::f64::consts::PI * medium_n * particle_size;
    let step = (wavelength_boundaries.1 - wavelength_boundaries.0) / 500.0; //((upper_limit - lower_limit) / 500.0).abs();
    let max_n: usize = 10;

    for i in 0..=499 {
        let current_wavelength = wavelength_boundaries.0 + step * (i + 1) as f64;
        let ref_index = get_ref_index(&ref_indices, current_wavelength * 1e6);
        let m = ref_index / medium_n;
        let coord_x = upper_x / current_wavelength;

        let jn_x = func::sj_array(coord_x, max_n);
        let yn_x = func::sy_array(coord_x, max_n);
        let psin_x = func::psi_array(coord_x, &jn_x);
        let xin_x = func::xi_array(coord_x, &jn_x, &yn_x);

        let jn_mx = func::sj_array(coord_x * m, max_n);
        let psin_mx = func::psi_array(coord_x * m, &jn_mx);
        let dn_mx = func::d_array(coord_x * m, &psin_mx);

        let an = func::a_array(coord_x, m, max_n, &psin_x, &xin_x, &dn_mx);
        let bn = func::b_array(coord_x, m, max_n, &psin_x, &xin_x, &dn_mx);

        let mul = (current_wavelength / medium_n).powi(2)
            / (2.0 * std::f64::consts::PI)
            / ((particle_size) * (particle_size) * std::f64::consts::PI);
        let sum_sca = mul
            * summation!(
                f64,
                j,
                1..=max_n,
                (2.0 * j as f64 + 1.0)
                    * (an[j].re.powi(2) + an[j].im.powi(2) + bn[j].re.powi(2) + bn[j].im.powi(2))
            );
        let sum_ext = mul
            * summation!(
                f64,
                j,
                1..=max_n,
                (2.0 * j as f64 + 1.0) * (an[j].re + bn[j].re)
            );

        c_sca.push((current_wavelength, sum_sca));
        c_ext.push((current_wavelength, sum_ext));
    }

    println!("calculating...");
    plot_png(
        "./results/cross-section-z-3.png",
        (2000, 400),
        "cross section c_sca & c_ext - 3",
        (wavelength_boundaries.0, wavelength_boundaries.1),
        (0.0, 10.0),
        &vec![c_sca, c_ext],
        &vec!["c_sca", "c_ext"],
        &vec![RED, BLUE],
    )?;

    Ok(())
}

/*
fn main() -> Result<(), Box<dyn std::error::Error>> {

    let mut c_sca: Vec<(f64, f64)> = Vec::with_capacity(500);

    let ref_indices = csv::parse("./res/refractive-index-silicon.csv");
    let wavelength_boundaries = (
        ref_indices[0].0 * 1e-6,
        ref_indices.last().unwrap().0 * 1e-6,
    );
    let medium_n = 1.0;
    let step = (wavelength_boundaries.1 - wavelength_boundaries.0) / 500.0; //((upper_limit - lower_limit) / 500.0).abs();
    let max_n: usize = 3;

    for i in 0..=499 {
        let current_wavelength = wavelength_boundaries.0 + step * (i + 1) as f64;
        let ref_index = get_ref_index(&ref_indices, current_wavelength * 1e6);
        let c_csa_local = integ_whole_particle(current_wavelength, medium_n, MU_0, ref_index, 85e-9, 1.0, max_n, 1000);
        println!("{} /// {}", i, c_csa_local);
        c_sca.push((current_wavelength, c_csa_local));
    }

    println!("calculating...");
    plot_png(
        "./results/integ-cross-section-z.png",
        (2000, 400),
        "integ cross section c_sca",
        (wavelength_boundaries.0, wavelength_boundaries.1),
        (0.0, 10.0),
        &vec![c_sca],
        &vec!["c_sca"],
        &vec![RED],
    )?;

    Ok(())
}*/
