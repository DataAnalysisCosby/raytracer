extern crate cgmath;

use std::f32;
use std::num;
use std::ops::*;
use cgmath::{Point, Vector, EuclideanVector, Point3, Vector3};

use bvh::BVH;
use geom::{Object, Ray};
use shade::Shade;
use sampler::Sampler;
use material::Material;

pub trait Light {
    fn get_direction(&mut self, hit_point: Point3<f32>) -> Vector3<f32>;
    fn get_color(&self, hit_point: Point3<f32>) -> Vector3<f32>;
    fn in_shadow<'a>(&self, ray: Ray, id: usize, bvh: &BVH<'a>) -> bool;
}

#[derive(Copy, Clone)]
pub struct Ambient {
    pub color: Vector3<f32>,
}

impl Light for Ambient {
    fn get_direction(&mut self, _: Point3<f32>) -> Vector3<f32> {
        Vector3::<f32>::new(0.0, 0.0, 0.0)
    }

    fn get_color(&self, _: Point3<f32>) -> Vector3<f32> {
        self.color
    }

    fn in_shadow<'a>(&self, ray: Ray, id: usize, bvh: &BVH<'a>) -> bool {
        false
    }
}


#[derive(Copy, Clone)]
pub struct PointLight {
    pub p: Point3<f32>,
    pub color: Vector3<f32>,
}

impl Light for PointLight {
    fn get_direction(&mut self, q: Point3<f32>) -> Vector3<f32> {
        (self.p - q).normalize()
    }

    fn get_color(&self,  q: Point3<f32>) -> Vector3<f32> {
        self.color
    }

    fn in_shadow<'a>(&self, ray: Ray, id: usize, bvh: &BVH<'a>) -> bool {
        let d = (ray.o - self.p).length();
        // If we exit prematurely, the area is in shadow.
        !bvh.raytrace(&ray, |i: usize, s: Shade| -> bool {
            id == i || s.tmin <= 0.0001 || s.tmin >= d
        })
    }
}

// An area light is both a light source and a geometric object.
#[derive(Copy, Clone)]
pub struct AreaLight<'a> {
    pub obj: &'a Object,
    pub color: Vector3<f32>,
    pub sample: (Point3<f32>, Vector3<f32>),
}

impl<'a> Light for AreaLight<'a> {
    fn get_direction(&mut self, q: Point3<f32>) -> Vector3<f32> {
        self.sample = self.obj.sample();
        (self.sample.0 - q).normalize()
    }

    fn in_shadow<'b>(&self, ray: Ray, id: usize, bvh: &BVH<'b>) -> bool {
        let d = (self.sample.0 - ray.o).length();
        // If we exit prematurely, the area is in shadow.
        !bvh.raytrace(&ray, |i: usize, s: Shade| -> bool {
            i == id || s.tmin <= 0.0001 || s.tmin >= d
        })
    }

    fn get_color(&self,  q: Point3<f32>) -> Vector3<f32> {
        let (p, n) = self.sample;
        let ndotd = n.dot((p - q).normalize());
        if ndotd > 0.0 {
            self.color
        } else {
            return Vector3::new(0.0f32, 0.0, 0.)
        }
    }
}

