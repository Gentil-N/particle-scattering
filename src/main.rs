use plotters::prelude::*;
use std::cmp::*;

const MED_N: f64 = 1.0;
const PART_N: f64 = 2.0;
const PART_SIZE: f64 = 0.1;

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
    for i in 2..(order+1) {
        yn[i] = (2 * i + 1) as f64 / x * yn[i - 1] - yn[i - 2];
    }
    return yn;
}

 /*SHOW BESSEL yn*/
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("./results/bessel_yn.png", (800, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("bessel yn", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0.0..20.0, -1.0..0.5)?;

    chart.configure_mesh().draw()?;

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

    chart
        .draw_series(LineSeries::new(
            y0, &RED,
        ))?
        .label("y0")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    chart
        .draw_series(LineSeries::new(
            y1, &BLUE,
        ))?
        .label("y1")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    chart
        .draw_series(LineSeries::new(
            y2, &GREEN,
        ))?
        .label("y2")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}

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
    let root = BitMapBackend::new("./results/bessel_jn.png", (800, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("bessel jn", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0.0..20.0, -0.5..1.2)?;

    chart.configure_mesh().draw()?;

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

    chart
        .draw_series(LineSeries::new(j0, &RED))?
        .label("j0")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    chart
        .draw_series(LineSeries::new(j1, &BLUE))?
        .label("j1")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    chart
        .draw_series(LineSeries::new(j2, &GREEN))?
        .label("j2")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}
*/

fn pitau_n(order: usize, theta : f64) -> (Vec<f64>, Vec<f64>) {
    assert!(order > 0);
    let mut pi_n = vec![0.0; order + 1];
    let mut tau_n = vec![0.0; order + 1];
    pi_n[0] = 0.0;
    pi_n[1] = 1.0;
    let nu = theta.cos();
    tau_n[1] = nu * pi_n[0] - 2.0 * pi_n[1];
    for i in 2..=order {
        let n = i as f64;
        pi_n[i] = (2.0 * n + 1.0) / (n - 1.0) * nu * pi_n[i - 1] - n / (n - 1.0)  * pi_n[i - 2];
        tau_n[i] = n * nu * pi_n[i] - (n + 1.0) * pi_n[i - 1]
    }
    return (pi_n, tau_n);
}