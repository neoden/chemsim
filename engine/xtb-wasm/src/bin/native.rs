//! Native reference run: `native <file.xyz> [reps]` prints the result buffer as JSON
//! (full precision) plus timings, for comparison with the WebAssembly build.
use std::time::Instant;
use xtb_wasm::Engine;

const BOHR: f64 = 0.529177210903;

fn jitter(base: &[[f64; 3]], step: usize) -> Vec<[f64; 3]> {
    // deterministic small displacements (same formula as the JS harness), ~0.01 bohr
    base.iter()
        .enumerate()
        .map(|(i, p)| {
            let mut q = *p;
            for k in 0..3 {
                q[k] += 0.01 * ((1.7 * step as f64 + 2.3 * i as f64 + 0.9 * k as f64).sin());
            }
            q
        })
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let text = std::fs::read_to_string(&args[1]).unwrap();
    let reps: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(20);
    let mut lines = text.lines();
    let nat: usize = lines.next().unwrap().trim().parse().unwrap();
    lines.next();
    let mut numbers = Vec::new();
    let mut pos = Vec::new();
    for line in lines.take(nat) {
        let mut it = line.split_whitespace();
        let z = match it.next().unwrap() {
            "H" => 1,
            "B" => 5,
            "C" => 6,
            "N" => 7,
            "O" => 8,
            s => panic!("element {s}"),
        };
        let v: Vec<f64> = it.map(|x| x.parse::<f64>().unwrap() / BOHR).collect();
        numbers.push(z);
        pos.push([v[0], v[1], v[2]]);
    }
    let t = Instant::now();
    let mut eng = Engine::new(&numbers, &pos, 0.0, 300.0, None).unwrap();
    let status = eng.singlepoint(&pos, true, true);
    let cold = t.elapsed().as_secs_f64() * 1e3;
    assert_eq!(status, 0);
    let out = eng.out.clone();
    let mut warm = Vec::new();
    for step in 1..=reps {
        let p = jitter(&pos, step);
        let t = Instant::now();
        assert_eq!(eng.singlepoint(&p, true, true), 0);
        warm.push(t.elapsed().as_secs_f64() * 1e3);
    }
    warm.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let last = eng.out.clone();
    println!(
        "{{\"cold_ms\":{cold:?},\"warm_median_ms\":{:?},\"out\":{out:?},\"last\":{last:?}}}",
        warm[warm.len() / 2]
    );
}
