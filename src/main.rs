use plotters::prelude::*;

const MED_N: f64 = 1;
const PART_N: f64 = 2;
const PART_SIZE: f64 = 0.1;
/*
xmodmap -e "keycode 105 = less greater less greater bar brokenbar lessthanequal greaterthanequal"
 */
fn yn_array(order: usize, x: f64) -> Vec<f64> {
    let mut yn: Vec<f64> = vec![0.0; order];
    let y0 = -x.cos() / x;
    let y1 = -x.cos() / (x * x) - x.sin() / x;

    return Vec::new();
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("0.png", (800, 500)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("y=x^2", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f32..1f32, 0f32..1f32)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            (-50..=50).map(|x| x as f32 / 50.0).map(|x| (x, x * x)),
            &RED,
        ))?
        .label("y = x^2")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}
