extern crate cgmath;
extern crate image;

use std::f32;
use std::num;
use std::ops::*;
use std::vec::Vec;
use cgmath::{Point, Vector, EuclideanVector, Point3, Vector3};

pub trait Texture {
    fn get_color(&self, local: Point3<f32>) -> Vector3<f32>;
//    fn get_texel_coordinates(local: Point3<f32>, w: u32, h: u32) -> (u32, u32);
}

pub struct SphereTexture {
    pub img: image::RgbImage,
}

impl Texture for SphereTexture {
    fn get_color(&self, local: Point3<f32>) -> Vector3<f32> {
        let theta = local.y.acos();
        let t_phi = local.x.atan2(local.z);
        let phi = if t_phi < 0.0 {
            t_phi + 2.0 * f32::consts::PI
        } else {
            t_phi
        };
        let u = phi * (1.0 / (2.0 * f32::consts::PI));
        let v = 1.0 - theta * f32::consts::FRAC_1_PI;
        let (w, h) = self.img.dimensions();
        let loc = (((w - 1) as f32 * u) as u32,
                   ((h - 1) as f32 * v) as u32);
        let pixel = self.img[loc];
        Vector3::<f32>::new( pixel.data[0] as f32 / 255.0,
                             pixel.data[1] as f32 / 255.0,
                             pixel.data[2] as f32 / 255.0)
    }
}
