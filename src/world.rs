extern crate cgmath;
extern crate image;

use std::cmp;
use std::num;
use std::f32;
use std::ops::*;
use cgmath::{Point, Vector, EuclideanVector, Point3, Vector3};

use sampler::Sampler;
use light::Light;
use camera::Camera;
use geom::{Object, Ray};
use bvh::BVH;
use shade::Shade;
use material::Material;
// use tracer::Tracer;

pub struct World<'a, T> where
    T: Camera + Iterator<Item=(u32, u32, Box<Fn(f32, f32) -> Ray>)> {
        bg_color: Vector3<f32>,
        lights: Vec<&'a mut  Light>,
        materials: Vec<&'a Material>,
        bvh: BVH<'a>,
        rays: Box<T>,
        sampler: Box<Sampler>,
}

impl<'a, T> World<'a, T> where
    T: Camera + Iterator<Item=(u32, u32, Box<Fn(f32, f32) -> Ray>)> {
    pub fn new<ST: Sampler + 'static>(rays: T, sampler: ST) -> Self {
        World{
            bg_color: Vector3::<f32>::new(1.0, 1.0, 1.0),
            lights: vec![],
            materials: vec![],
            bvh: BVH::new(),
            rays: Box::new(rays),
            sampler: Box::new(sampler),
        }
    }


    pub fn add_object(&mut self, obj: &'a Object) {
        self.bvh.insert_object(obj);
    }

    pub fn add_light(&mut self, light:&'a mut  Light) -> usize {
        self.lights.push(light);
        self.lights.len() - 1
    }

    pub fn add_material(&mut self, mat:&'a Material) -> usize {
        self.materials.push(mat);
        self.materials.len() - 1
    }

    fn collide_ray(&mut self, ray: &Ray) -> (Vector3<f32>, usize) {
        let mut shade = Shade{
            local: Point3::<f32>::new(0.0, 0.0, 0.0),
            normal: Vector3::<f32>::new(0.0, 1.0, 0.0),
            mat: 0,
            tmin: f32::INFINITY,
            depth: 0,
        };
        let mut id = 0;
        self.bvh.raytrace(ray, |i: usize, s: Shade| {
            if s.tmin < shade.tmin {
                shade = s;
                id = i;
            }
            true
        });
        if shade.tmin == f32::INFINITY {
            (self.bg_color, 1)
        } else {
            (self.materials[shade.mat].shade(shade, ray, &mut self.lights, id,
                                             &self.bvh, &self.materials),
             shade.depth)
        }
    }

    pub fn render_scene(&mut self, filename: &str) {
        let (w, h) = self.rays.dims();
        let size = (w * h) as f32;
        let mut count = 0;
        let mut rays = 0;
        let mut img = image::ImageBuffer::<image::Rgb<u8>>::new(w, h);
        while let Some((i, j, ray_gen)) = self.rays.next() {
            count += 1;
            let completion = (((count as f32) / size) * 100.0) as i32;
            print!("Rendering...%{}\r", completion);
            let samples = self.sampler.sample_unit_square();
            let mut color = [0f32, 0.0, 0.0];
            for &(x, y) in samples.iter() {
                let ray = ray_gen(x, y);
                let (c, num_cast) = self.collide_ray(&ray);
                rays += num_cast;
                color[0] += c.x;
                color[1] += c.y;
                color[2] += c.z;
            }
            for k in 0..3 {
                color[k] /= samples.len() as f32;
            }
            img.get_pixel_mut(i, j).data = [ (color[0] * 255.0).min(255.0) as u8,
                                            (color[1] * 255.0).min(255.0) as u8,
                                            (color[2] * 255.0).min(255.0) as u8];
        }
        println!("Rendering complete. {} rays cast.", rays);
        img.save(filename).unwrap()
    }
}


