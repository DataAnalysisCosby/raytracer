extern crate cgmath;

use cgmath::{Point3, Vector3};

#[derive(Copy, Clone)]
pub struct Ray {
    pub o: Point3<f32>,
    pub d: Vector3<f32>,
}

