extern crate cgmath;

use std::f32;
use std::num;
use std::ops::*;
use std::vec::Vec;
use std::cell::Cell;
use cgmath::{Point, Vector, EuclideanVector, Point3, Vector3};

use bvh::BVH;
use geom::Ray;
use light::Light;
use shade::Shade;
use texture::Texture;

#[derive(Copy, Clone)]
struct Lambertian {
    kd: f32,
    cd: Vector3<f32>,
}

impl Lambertian {
    fn f(&self) -> Vector3<f32> {
        self.cd * self.kd * f32::consts::FRAC_1_PI
    }

    fn rho(&self) -> Vector3<f32> {
        self.cd * self.kd
    }
}

#[derive(Copy, Clone)]
struct Specular {
    exp: f32,
    ks: Vector3<f32>,
}

impl Specular {
    fn f(&self, normal: Vector3<f32>, wo: Vector3<f32>, wi: Vector3<f32>) -> Vector3<f32> {
        let ndotwo = normal.dot(wo);
        let r = -wo + normal * ndotwo * 2.0;
        let rdotwo = r.dot(wo);
        if rdotwo > 0.0 {
            self.ks * rdotwo.powf(self.exp)
        } else {
            Vector3::<f32>::new(0.0, 0.0, 0.0)
        }
    }
}

#[derive(Copy, Clone)]
pub struct GlossySpecular {
    pub exp: f32,
    pub ks: Vector3<f32>,
}

impl GlossySpecular {
    fn f(&self, normal: Vector3<f32>, sp: Point3<f32>, wo: Vector3<f32>) -> (Vector3<f32>, Vector3<f32>, f32) {
        let ndotwo = normal.dot(wo);
        let r = -wo + normal * ndotwo * 2.0;
        let u = Vector3::new(0.00424f32, 1.0, 0.00764).cross(r).normalize();
        let v = u.cross(r);
        let wi_t = u * sp.x + v * sp.y + r * sp.z;
        let wi = if normal.dot(wi_t) < 0.0 {
            u * -sp.x - v * sp.y + r * sp.z
        } else {
            wi_t
        };
        let phong_lobe = r.dot(wi).powf(self.exp);
        let pdf = phong_lobe * normal.dot(wi);
        (wi, self.ks * phong_lobe, pdf)
    }
}

#[derive(Copy, Clone)]
struct PerfectSpecular {
    kr: f32,
    cr: Vector3<f32>,
}

impl PerfectSpecular {
    fn f(&self, normal: Vector3<f32>, wo: Vector3<f32>) -> (Vector3<f32>, Vector3<f32>) {
        let ndotwo = normal.dot(wo);
        let wi = normal * 2.0 * ndotwo - wo;
        (wi, self.cr * self.kr / normal.dot(wi).abs())
    }
}

#[derive(Copy, Clone)]
struct PerfectTransmitter {
    kt: f32,
    ior: f32,
}

impl PerfectTransmitter {
    fn tir(&self, normal: Vector3<f32>, wo: Vector3<f32>) -> bool {
        let cos_thetai = normal.dot(wo);
        let eta = if cos_thetai < 0.0 {
            1.0 / self.ior
        } else {
            self.ior
        };
        ((1.0 - (1.0 - cos_thetai * cos_thetai) / (eta * eta)) < 0.0)
    }

    fn f(&self, normal: Vector3<f32>, wo: Vector3<f32>) -> (Vector3<f32>, Vector3<f32>) {
        let ndotwo = normal.dot(wo);
        let (n, cos_thetai, eta) = if ndotwo < 0.0 {
            (-normal, -ndotwo, 1.0/self.ior)
        } else {
            (normal, ndotwo, self.ior)
        };
        let t = 1.0 - (1.0 - cos_thetai * cos_thetai) / (eta * eta);
        let wt = -wo / eta - n * (t.sqrt() - cos_thetai / eta);
        (wt, Vector3::new(1.0f32, 1.0, 1.0) * self.kt / ((eta * eta) * normal.dot(wt).abs()))
    }
}

pub trait Material {
    fn shade<'a>(&self, s: Shade, ray: &Ray, lights: &mut Vec<&'a mut Light>,
                 id: usize, bvh: &BVH<'a>, materials: &Vec<&'a Material>) -> Vector3<f32>;
}


#[derive(Copy, Clone)]
pub struct Emissive {
    ls: f32,
    ce: Vector3<f32>,
}

impl Emissive {
    pub fn new(ls: f32, ce: Vector3<f32>) -> Emissive {
        Emissive{
            ls: ls,
            ce: ce,
        }
    }
}

impl Material for Emissive {
    fn shade<'a>(&self, s: Shade, ray: &Ray, lights: &mut Vec<&'a mut Light>,
                 id: usize, bvh: &BVH<'a>, materials: &Vec<&'a Material>) -> Vector3<f32> {
        if -s.normal.dot(ray.d) > 0.0 {
            self.ce * self.ls
        } else {
            Vector3::<f32>::new(0.0, 0.0, 0.0)
        }
    }
}

#[derive(Copy, Clone)]
pub struct Matte {
    ambient_brdf: Lambertian,
    diffuse_brdf: Lambertian,
}

impl Matte {
    pub fn new(ka: f32, kd: f32, c: Vector3<f32>) -> Matte {
        Matte{
            ambient_brdf: Lambertian{
                kd: ka,
                cd: c,
            },
            diffuse_brdf: Lambertian{
                kd: kd,
                cd: c,
            },
        }
    }
}

impl Material for Matte {
    fn shade<'a>(&self, s: Shade, ray: &Ray, lights: &mut Vec<&'a mut Light>,
                 id: usize, bvh: &BVH<'a>, materials: &Vec<&'a Material>) -> Vector3<f32> {
        let wo = -ray.d;
        let mut l = self.ambient_brdf.rho() * lights[0].get_color(s.local);
        for i in 1..lights.len() {
            let wi = lights[i].get_direction(s.local);
            let ndotwi = s.normal.dot(wi);
            if ndotwi > 0.0 {
                l = l + self.diffuse_brdf.f() * lights[i].get_color(s.local) * ndotwi;
            }
        }
        l
    }
}

#[derive(Copy, Clone)]
pub struct Phong {
    ambient_brdf: Lambertian,
    diffuse_brdf: Lambertian,
    specular_brdf: Specular,
}


impl Phong {
    pub fn new(ka: f32, kd: f32, shine: f32, c: Vector3<f32>) -> Phong {
        Phong{
            ambient_brdf: Lambertian{
                kd: ka,
                cd: c,
            },
            diffuse_brdf: Lambertian{
                kd: kd,
                cd: c,
            },
            specular_brdf: Specular{
                exp: shine,
                ks: c,
            },
        }
    }
}

impl Material for Phong {
    fn shade<'a>(&self, s: Shade, ray: &Ray, lights: &mut Vec<&'a mut Light>,
                 id: usize, bvh: &BVH<'a>, materials: &Vec<&'a Material>) -> Vector3<f32> {
        let wo = -ray.d;
        let mut l = self.ambient_brdf.rho() * lights[0].get_color(s.local);
        for i in 1..lights.len() {
            let wi = lights[i].get_direction(s.local);
            let ndotwi = s.normal.dot(wi);
            if ndotwi > 0.0 {
                if !lights[i].in_shadow(Ray::new(s.local, wi),id, bvh) {
                    l = l + (self.diffuse_brdf.f() +
                             self.specular_brdf.f(s.normal, wo, wi)) * lights[i].get_color(s.local) * ndotwi;
                }
            }
        }
        l
    }
}

#[derive(Copy, Clone)]
pub struct Reflective {
    kr: f32,
    cr: Vector3<f32>,
    ph: Phong,
}

impl Reflective {
    pub fn new(kr: f32, cr: Vector3<f32>, ph: Phong) -> Reflective {
        Reflective{
            kr: kr,
            cr: cr,
            ph: ph,
        }
    }
}

const MAX_DEPTH: usize = 500;

impl Material for Reflective {
    fn shade<'a>(&self, s: Shade, ray: &Ray, lights: &mut Vec<&'a mut Light>,
                 id: usize, bvh: &BVH<'a>, materials: &Vec<&'a Material>) -> Vector3<f32> {
        let wo = -ray.d;
        let ndotwo = s.normal.dot(wo);
        let wi = -wo + s.normal * 2.0 * ndotwo;
        let fr = self.cr * self.kr / s.normal.dot(wi);
        let reflected = Ray::new(s.local, wi);
        let mut shade = Shade{
            local: Point3::<f32>::new(0.0, 0.0, 0.0),
            normal: Vector3::<f32>::new(0.0, 1.0, 0.0),
            mat: 0,
            tmin: f32::INFINITY,
            depth: 1,
        };
        let mut reflected_id = 0;
        bvh.raytrace(&reflected, |i: usize, s: Shade| {
            if i != id && s.tmin < shade.tmin {
                shade = s;
                reflected_id = i;
            }
            true
        });
        self.ph.shade(s, ray, lights, id, bvh, materials) + fr * s.normal.dot(wi) *
            if shade.depth + s.depth > MAX_DEPTH {
                Vector3::<f32>::new(0.0, 0.0, 0.0)
            } else if shade.tmin == f32::INFINITY {
                Vector3::<f32>::new(1.0, 1.0, 1.0)
            } else {
                materials[shade.mat].shade(shade, &reflected, lights,
                                           reflected_id, bvh, materials)
            }
    }
}

#[derive(Copy, Clone)]
pub struct Transparent {
    reflective: PerfectSpecular,
    specular: PerfectTransmitter,
    ph: Phong,
}

impl Transparent {
    pub fn new(kr: f32, cr: Vector3<f32>, kt: f32, ior: f32, ph: Phong) -> Transparent {
        Transparent{
            reflective: PerfectSpecular{
                kr: kr,
                cr: cr,
            },
            specular: PerfectTransmitter{
                kt: kt,
                ior: ior,
            },
            ph: ph,
        }
    }
}

impl Material for Transparent {
    fn shade<'a>(&self, s: Shade, ray: &Ray, lights: &mut Vec<&'a mut Light>,
                 id: usize, bvh: &BVH<'a>, materials: &Vec<&'a Material>) -> Vector3<f32> {
        let wo = -ray.d;
        let ndotwo = s.normal.dot(wo);
        let (wi, fr) = self.reflective.f(s.normal, wo);
        let reflected = Ray::new(s.local, wi);
        let mut shade = Shade{
            local: Point3::<f32>::new(0.0, 0.0, 0.0),
            normal: Vector3::<f32>::new(0.0, 1.0, 0.0),
            mat: 0,
            tmin: f32::INFINITY,
            depth: 1,
        };
        let mut reflected_id = 0;
        bvh.raytrace(&reflected, |i: usize, s: Shade| {
            if i != id && s.tmin < shade.tmin {
                shade = s;
                reflected_id = i;
            }
            true
        });
        let reflected_shade = if shade.depth + s.depth > MAX_DEPTH {
            Vector3::new(0.0f32, 0.0, 0.0)
        } else if shade.tmin == f32::INFINITY {
            Vector3::new(1.0f32, 1.0, 1.0)
        } else {
            materials[shade.mat].shade(shade, &reflected, lights,
                                       reflected_id, bvh, materials)
        };
        self.ph.shade(s, ray, lights, id, bvh, materials) +
            if self.specular.tir(s.normal, wo) {
                reflected_shade
            } else {
                let (wt, ft) = self.specular.f(s.normal, wo);
                let transmitted = Ray::new(s.local, wt);
                let mut shade = Shade{
                    local: Point3::<f32>::new(0.0, 0.0, 0.0),
                    normal: Vector3::<f32>::new(0.0, 1.0, 0.0),
                    mat: 0,
                    tmin: f32::INFINITY,
                    depth: 1,
                };
                let mut transmitted_id = 0;
                bvh.raytrace(&transmitted, |i: usize, s: Shade| {
                    if i != id && s.tmin < shade.tmin {
                        shade = s;
                        reflected_id = i;
                    }
                    true
                });
                let transmitted_shade = if shade.depth + s.depth > MAX_DEPTH {
                    Vector3::new(0.0f32, 0.0, 0.0)
                } else if shade.tmin == f32::INFINITY {
                    Vector3::new(1.0f32, 1.0, 1.0)
                } else {
                    materials[shade.mat].shade(shade, &reflected, lights,
                                               transmitted_id, bvh, materials)
                };
                reflected_shade * fr * s.normal.dot(wi).abs() +
                    transmitted_shade * ft * s.normal.dot(wt).abs()
            }
    }
}

#[derive(Copy, Clone)]
pub struct MatteText<'a> {
    pub ka: f32,
    pub kd: f32,
    pub tex: &'a Texture,
}

impl<'a> Material for MatteText<'a> {
    fn shade<'b>(&self, s: Shade, ray: &Ray, lights: &mut Vec<&'b mut Light>,
                 id: usize, bvh: &BVH<'b>, materials: &Vec<&'b Material>) -> Vector3<f32> {
        let wo = -ray.d;
        let ndotwo = s.normal.dot(wo);
        let cd = self.tex.get_color(Point3::from_vec(s.normal));
        let mut l = cd * self.ka * lights[0].get_color(s.local);
        let diffuse = Lambertian{
            kd: self.kd,
            cd: cd,
        };
        for i in 1..lights.len() {
            let wi = lights[i].get_direction(s.local);
            let ndotwi = s.normal.dot(wi);
            let ndotwo = s.normal.dot(wo);
            if ndotwi > 0.0 && ndotwo > 0.0 {
                if !lights[i].in_shadow(Ray::new(s.local, wi),id, bvh) {
                    l = l + diffuse.f() * lights[i].get_color(s.local) * ndotwi;
                }
            }
        }
        l
    }
}


pub struct GlossyReflector {
    pub samples: Vec<(f32, f32)>,
    pub curr_s: Cell<usize>,
    pub gs: GlossySpecular,
    pub ph: Phong,
}

impl Material for GlossyReflector {
    fn shade<'a>(&self, s: Shade, ray: &Ray, lights: &mut Vec<&'a mut Light>,
                 id: usize, bvh: &BVH<'a>, materials: &Vec<&'a Material>) -> Vector3<f32> {
        let wo = -ray.d;
        let ndotwo = s.normal.dot(wo);
        let l = self.ph.shade(s, ray, lights, id, bvh, materials);
        let (x, y) = self.samples[self.curr_s.get()];
        if self.curr_s.get() == self.samples.len() - 1 {
            self.curr_s.set(0);
        } else {
            self.curr_s.set(self.curr_s.get() + 1);
        }
        let phi = 2.0 * f32::consts::PI * x;
        let cos_theta = (1.0 - y).powf(1.0 / (self.gs.exp + 1.0));
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        let p = Point3::new(sin_theta * phi.cos(),
                            sin_theta * phi.sin(),
                            cos_theta);
        let (wi, fr, pdf) = self.gs.f(s.normal, p, wo);
        let reflected = Ray::new(s.local, wi);
        let mut shade = Shade{
            local: Point3::<f32>::new(0.0, 0.0, 0.0),
            normal: Vector3::<f32>::new(0.0, 1.0, 0.0),
            mat: 0,
            tmin: f32::INFINITY,
            depth: 1,
        };
        let mut reflected_id = 0;
        bvh.raytrace(&reflected, |i: usize, s: Shade| {
            if i != id && s.tmin < shade.tmin {
                shade = s;
                reflected_id = i;
            }
            true
        });
        let reflected_shade = if shade.depth + s.depth > MAX_DEPTH {
            Vector3::new(1.0f32, 1.0, 1.0)
        } else if shade.tmin == f32::INFINITY {
            Vector3::new(0.0f32, 0.0, 0.0)
        } else {
            materials[shade.mat].shade(shade, &reflected, lights,
                                       reflected_id, bvh, materials)
        };
        l + fr * reflected_shade * s.normal.dot(wi) / pdf
    }
}
