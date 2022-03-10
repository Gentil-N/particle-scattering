
pub fn linear_interpolation_2D_by_coeff(a: (f64, f64), b: (f64, f64), coeff: f64) -> (f64, f64) {
    (a.0 + (b.0 - a.0) * coeff, a.1 + (b.1 - a.1) * coeff)
}

pub fn linear_interpolation_2D(a: (f64, f64), b: (f64, f64), c: f64) -> (f64, f64) {
    let coeff = (c - a.0) / (b.0 - a.0);
    linear_interpolation_2D_by_coeff(a, b, coeff)
}