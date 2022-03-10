use std::fs;

pub fn parse(filename: &str) -> Vec<(f64, f64, f64)> {
    let contents = fs::read_to_string(filename).unwrap();
    let lines = contents.lines().collect::<Vec<&str>>();
    let mut data: Vec<(f64, f64, f64)> = Vec::with_capacity(lines.len());
    for line in lines {
        let triple_str = line.split(",").collect::<Vec<&str>>();
        data.push((
            triple_str[0].trim().parse::<f64>().unwrap(),
            triple_str[1].trim().parse::<f64>().unwrap(),
            triple_str[2].trim().parse::<f64>().unwrap(),
        ));
    }
    data
}