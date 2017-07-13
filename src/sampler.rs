extern crate cgmath;
extern crate rand;

use rand::distributions::{IndependentSample, Range};

pub trait Sampler {
    fn sample_unit_square(&self) -> Vec<(f32, f32)>;
}

pub struct SingleSample;

impl SingleSample {
    pub fn new() -> SingleSample {
        SingleSample
    }
}

impl Sampler for SingleSample {
    fn sample_unit_square(&self) -> Vec<(f32, f32)> {
        vec![ (0.0, 0.0) ]
    }
}

pub struct MultiJittered {
    n: usize,
}

impl MultiJittered {
    pub fn new(n: usize) -> MultiJittered {
        MultiJittered{
            n: n,
        }
    }
}

impl Sampler for MultiJittered {
    fn sample_unit_square(&self) -> Vec<(f32, f32)> {
        let num_samples = self.n * self.n;
        let subcell_width = 1.0 / num_samples as f32;
        let range = Range::new(0f32, subcell_width);
        let mut rng = rand::thread_rng();
        let mut samples: Vec<(f32, f32)>  = Vec::new();
        for i in 0..num_samples {
            samples.push((0.0, 0.0));
        }
        for i in 0..self.n {
            for j in 0..self.n {
                samples[i * self.n + j] =
                    ((i as f32 * self.n as f32 + j as f32) * subcell_width + range.ind_sample(&mut rng),
                     (j as f32 * self.n as f32 + i as f32) * subcell_width + range.ind_sample(&mut rng));
            }
        }
        for i in 0..self.n {
            for j in 0..self.n {
                let range = Range::new(j, self.n);
                let k = range.ind_sample(&mut rng);
                let (x1, y1) = samples[i * self.n + j];
                let (x2, y2) = samples[i * self.n + k];
                samples[i * self.n + j] = (x2, y1);
                samples[i * self.n + k] = (x1, y2);
            }
        }
        for i in 0..self.n {
            for j in 0..self.n {
                let range = Range::new(j, self.n);
                let k = range.ind_sample(&mut rng);
                let (x1, y1) = samples[i * self.n + j];
                let (x2, y2) = samples[i * self.n + k];
                samples[i * self.n + j] = (x1, y2);
                samples[i * self.n + k] = (x2, y1);
            }
        }
        samples
    }
}
