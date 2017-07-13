extern crate cgmath;

use cgmath::{Point3, Vector3};

#[derive(Copy, Clone)]
pub struct Shade {
    pub local:  Point3<f32>,
    pub normal: Vector3<f32>,
    pub mat: usize,
//    pub color: [u8; 3],
    pub tmin: f32,
    pub depth: usize,
}
