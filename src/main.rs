#[macro_use]
extern crate cgmath;
extern crate rand;
extern crate image;

mod bvh;
mod world;
mod geom;
mod shade;
mod light;
mod camera;
mod sampler;
mod material;
mod mesh;
mod texture;

use std::env;
use std::io::{self, Write};
use std::cell::Cell;
use cgmath::{Point3, Vector3, Rotation3, EuclideanVector, Vector, Decomposed, Matrix3, Quaternion};

use bvh::BVH;
use world::World;
use light::{Ambient, PointLight, AreaLight};
use sampler::{Sampler, SingleSample, MultiJittered};
use geom::{Plane, Sphere, Triangle, Instance};
use material::{GlossySpecular, Matte, MatteText, Phong, Reflective, Emissive, Transparent, GlossyReflector };
use camera::{OrthogonalProj, PerspectiveProj};
use texture::{Texture, SphereTexture};

fn make_spheres(w: f32, h: f32, n: usize) -> Vec<Sphere>{
    let mut spheres = Vec::<Sphere>::new();
    let w_step = w / (n as f32);
    let h_step = h / (n as f32);
    for i in 0..n {
        for j in 0..n {
            spheres.push(Sphere{
                c: Point3::<f32>::new(w_step * (i as f32) - w / 2.0,
                                      h_step * (j as f32) - h / 2.0,
                                      0.0),
                r: 30.0,
                mat: 4,
            });
        }
    }
    spheres
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let texture = SphereTexture{
        img: image::open("EarthHighRes.jpg").unwrap().to_rgb(),
    };
    {
        let left_sphere = Sphere{
            c: Point3::<f32>::new(-50.0, 0.0, 0.0),
            r: 25.0,
            mat: 0,
        };
        let center_sphere = Instance{
            geom: Box::new(Sphere{
                c: Point3::<f32>::new(0.0, 0.0, -30.0),
                r: 25.0,
                mat: 1,
            }),
            inv_transform: Matrix3::<f32>::new(
                1.0 / 2.0, 0.0, 0.0,
                0.0, 1.0 / 1.5, 0.0,
                0.0, 0.0, 1.0,
                ),
        };
        let right_sphere = Sphere{
                c: Point3::<f32>::new(40.0, 0.0, 30.0),
                r: 25.0,
                mat: 4,
        };
        let plane = Plane{
            p: Point3::<f32>::new(0.0, 0.0, -100.0),
            n: Vector3::<f32>::new(0.0, 0.0, 1.0),
            mat: 3,
        };
        let plane2 = Plane{
            p: Point3::<f32>::new(0.0, 0.0, 500.0),
            n: Vector3::<f32>::new(0.0, 0.0, 1.0),
            mat: 3,
        };
        let mut ambient_light = Ambient{
            color: Vector3::<f32>::new(1.0, 1.0, 1.0),
        };
        let mut point_light =  PointLight{
            p: Point3::<f32>::new(500.0, 500.0, 500.0),
            color: Vector3::<f32>::new(0.5, 0.5, 0.5),
        };
        /*
        let mut area_light = AreaLight{
            obj: &right_sphere,
            color: Vector3::<f32>::new(1.0, 1.0, 1.0),
            sample: (Point3::new(1.0f32, 1.0, 1.0),
                     Vector3::new(1.0f32, 1.0, 1.0)),
        };*/
        let mat0 = MatteText{
            ka: 1.0,
            kd: 1.0,
            tex: &texture,
        };
        let sampler = MultiJittered::new(32);
        let mat1 = GlossyReflector{
            samples: sampler.sample_unit_square(),
            curr_s: Cell::new(0),
            gs: GlossySpecular{
                exp: 50.0,
                ks: Vector3::new(0.3f32, 0.5f32, 0.1f32),
            },
            ph: Phong::new(0.25, 0.5, 0.15, Vector3::<f32>::new(0.25, 0.25, 0.25)),
        };
        let mat2 = Phong::new(0.25, 0.65, 8.0,
                              Vector3::<f32>::new(1.0, 0.0, 1.0));
        let mat3 = Phong::new(0.25, 0.65, 8.0,
                              Vector3::<f32>::new(1.0, 1.0, 0.0));
        let mat4 = Transparent::new(0.1, Vector3::<f32>::new(1.0, 5.0, 0.0),
                                    0.1, 0.9,
                                    Phong::new(0.5, 0.5, 2000.0,
                                               Vector3::<f32>::new(1.0, 1.0, 0.0)));
        let mut world = World::new(
            PerspectiveProj::new(2000, 2000, 850.0, 1.0 / 2.0,
                                 Point3::<f32>::new(0.0, 100.0, 500.0),
                                 Point3::<f32>::new(0.0, 0.0, 0.0),
                                 Vector3::<f32>::new(1.0, 0.0, 0.0)),
            MultiJittered::new(6));
        world.add_light(&mut ambient_light);
        world.add_light(&mut point_light);
//        world.add_light(&mut area_light);
        world.add_material(&mat0);
        world.add_material(&mat1);
        world.add_material(&mat2);
        world.add_material(&mat3);
        world.add_material(&mat4);
        world.add_object(&left_sphere);
        world.add_object(&center_sphere);
        world.add_object(&right_sphere);
        world.add_object(&plane);
        world.render_scene("nicescene10.png");
    }
}
