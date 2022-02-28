use plotters::prelude::*;
use scilib::math::complex::Complex;
use std::cmp::*;

const MED_N: f64 = 1.0;
const PART_N: f64 = 2.0;
const PART_SIZE: f64 = 1.0;
const WAVE_LENGHT: f64 = 0.1;

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

/*
fn main() -> Result<(), Box<dyn std::error::Error>> {

    let mut j0: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut j1: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut j2: Vec<(f64, f64)> = Vec::with_capacity(400);
    for i in 0..399 {
        let coord_x = (i + 1) as f64 / 20.0;
        let j_n = j_n(2, coord_x);
        j0.push((coord_x, j_n[0]));
        j1.push((coord_x, j_n[1]));
        j2.push((coord_x, j_n[2]));
    }

    plot_png("./results/bessel-jn.png", (800, 800), "bessel jn", (0.0, 20.0), (-0.5, 1.2), &vec![j0, j1, j2], &vec!["j0", "j1", "j2"], &vec![RED, BLUE, GREEN])?;

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

/*
fn main() -> Result<(), Box<dyn std::error::Error>> {

    let mut xi0: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut xi1: Vec<(f64, f64)> = Vec::with_capacity(400);
    let mut xi2: Vec<(f64, f64)> = Vec::with_capacity(400);
    for i in 0..399 {
        let coord_x = (i + 1) as f64 / 20.0;
        let jn = j_n(2, coord_x);
        let yn = y_n(2, coord_x);
        let xin = xi_n(&jn, &yn, coord_x);
        xi0.push((coord_x, xin[0].im));
        xi1.push((coord_x, xin[1].im));
        xi2.push((coord_x, xin[2].im));
    }

    plot_png("./results/riccati-bessel-xin-imaginary.png", (800, 800), "riccati-bessel xin (i)", (0.0, 20.0), (-10.0, 10.0), &vec![xi0, xi1, xi2], &vec!["xi0", "xi1", "xi2"], &vec![RED, BLUE, GREEN])?;

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

    for i in 0..functions.len() {
        chart
            .draw_series(LineSeries::new(functions[i].clone(), &colors[i]))?
            .label(labels[i])
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    }

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}

fn a_n(order: usize,  psin: &Vec<f64>, xin: &Vec<Complex>, dn: &Vec<f64>, m: f64, x: f64) -> Vec<Complex> {
    assert!(psin.len() == xin.len() && xin.len() == dn.len());
    let mut an: Vec<Complex> = Vec::with_capacity(psin.len());
    an.push(Complex::from(0.0, 0.0));
    for j in 1..=order {
        let n = j as f64;
        let coeff_aj = dn[j] / m + n / x;
        let aj: Complex = Complex::from(coeff_aj * psin[j] - psin[j - 1], 0)
            / (coeff_aj * xin[j] - xin[j - 1]);
        an.push(aj);
    }
    return an;
}

fn b_n(order: usize,  psin: &Vec<f64>, xin: &Vec<Complex>, dn: &Vec<f64>, m: f64, x: f64) -> Vec<Complex> {
    assert!(psin.len() == xin.len() && xin.len() == dn.len());
    let mut bn: Vec<Complex> = Vec::with_capacity(psin.len());
    bn.push(Complex::from(0.0, 0.0));
    for j in 1..=order {
        let n = j as f64;
        let coeff_bj = dn[j] * m + n / x;
            let bj: Complex = Complex::from(coeff_bj * psin[j] - psin[j - 1], 0)
                / (coeff_bj * xin[j] - xin[j - 1]);
                bn.push(bj);
    }
    return bn;
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    //let upper_x = 2.0 * std::f64::consts::PI * MED_N * PART_SIZE;
    let m = PART_N / MED_N;
    let mut an_r: Vec<Vec<(f64, f64)>> = vec![Vec::with_capacity(500); 4];
    let mut an_i: Vec<Vec<(f64, f64)>> = vec![Vec::with_capacity(500); 4];
    let mut bn_r: Vec<Vec<(f64, f64)>> = vec![Vec::with_capacity(500); 4];
    let mut bn_i: Vec<Vec<(f64, f64)>> = vec![Vec::with_capacity(500); 4];

    for i in 0..499 {
        let coord_x = (i + 1) as f64 / 10.0;

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
            an_r[j].push((coord_x, an[j].re));
            an_i[j].push((coord_x, an[j].im));
            bn_r[j].push((coord_x, bn[j].re));
            bn_i[j].push((coord_x, bn[j].im));
        }
    }

    plot_png(
        "./results/scattering-coeff-anr.png",
        (2000, 400),
        "scattering coeff anr",
        (0.0, 50.0),
        (-0.2, 1.2),
        &an_r[1..].to_vec(),
        &vec!["anr1", "anr2", "anr3"],
        &vec![RED, BLUE, GREEN],
    )?;

    plot_png(
        "./results/scattering-coeff-ani.png",
        (2000, 400),
        "scattering coeff ani",
        (0.0, 50.0),
        (-1.0, 1.0),
        &an_i[1..].to_vec(),
        &vec!["ani1", "ani2", "ani3"],
        &vec![RED, BLUE, GREEN],
    )?;

    plot_png(
        "./results/scattering-coeff-bnr.png",
        (2000, 400),
        "scattering coeff bnr",
        (0.0, 50.0),
        (-0.2, 1.2),
        &bn_r[1..].to_vec(),
        &vec!["bnr1", "bnr2", "bnr3"],
        &vec![RED, BLUE, GREEN],
    )?;

    plot_png(
        "./results/scattering-coeff-bni.png",
        (2000, 400),
        "scattering coeff bni",
        (0.0, 50.0),
        (-1.0, 1.0),
        &bn_i[1..].to_vec(),
        &vec!["bni1", "bni2", "bni3"],
        &vec![RED, BLUE, GREEN],
    )?;

    Ok(())
}