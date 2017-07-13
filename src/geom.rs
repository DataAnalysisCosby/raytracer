extern crate cgmath;
extern crate rand;

use std::f32;
use std::num;
use std::ops::*;
use rand::distributions::{IndependentSample};
use cgmath::{Point, EuclideanVector, Vector, SquareMatrix, Point3, Vector3, Matrix3};

use material::Material;
use shade::Shade;

#[derive(Copy, Clone)]
pub struct Ray {
    pub o: Point3<f32>,
    pub d: Vector3<f32>,
}

impl Ray {
    pub fn new(p: Point3<f32>, d: Vector3<f32>) -> Ray {
        Ray{
            o: p,
            d: d,
        }
    }
}

#[derive(Copy, Clone)]
pub struct AABB {
    min: Point3<f32>,
    max: Point3<f32>
}

impl AABB {
    pub fn empty() -> AABB {
        AABB{
            min: Point3::<f32>::new(0.0, 0.0, 0.0),
            max: Point3::<f32>::new(0.0, 0.0, 0.0),
        }
    }

    pub fn collides(a: AABB, b:  AABB) -> bool {
        a.max.x >= b.min.x && a.min.x <= b.max.x &&
            a.max.y >= b.min.y && a.min.y <= b.max.y &&
            a.max.z >= b.min.z && a.min.z <= b.max.z
    }

    pub fn combine(a: AABB, b: AABB) -> AABB {
        AABB{
            min: Point3::<f32>::new(a.min.x.min(b.min.x),
                                  a.min.y.min(b.min.y),
                                  a.min.z.min(b.min.z)),
            max: Point3::<f32>::new(a.max.x.max(b.max.x),
                                  a.max.y.max(b.max.y),
                                  a.max.z.max(b.max.z)),
        }
    }

    pub fn surface_area(&self) -> f32 {
        let w = self.max.x - self.min.x;
        let h = self.max.y - self.min.y;
        let d = self.max.z - self.min.z;
        2.0 * w * h + 2.0 * w * d + 2.0 * h * d
    }

    pub fn collides_ray(&self, ray: &Ray) -> bool {
        let mut tmin: f32 = 0.0;
        let mut tmax: f32 = f32::INFINITY;
        if ray.d.x == 0.0 && (ray.o.x < self.min.x || ray.o.x > self.max.x) {
            return false;
        } else {
            let ood = 1.0 / ray.d.x;
            let t1 = (self.min.x - ray.o.x) * ood;
            let t2 = (self.max.x - ray.o.x) * ood;
            let (t1, t2) = if t1 > t2 {
                (t2, t1)
            } else {
                (t1, t2)
            };
            tmin = tmin.max(t1);
            tmax = tmax.min(t2);
            if tmin > tmax {
                return false;
            }
        }
        if ray.d.y == 0.0 && (ray.o.y < self.min.y || ray.o.y > self.max.y) {
            return false;
        } else {
            let ood = 1.0 / ray.d.y;
            let t1 = (self.min.y - ray.o.y) * ood;
            let t2 = (self.max.y - ray.o.y) * ood;
            let (t1, t2) = if t1 > t2 {
                (t2, t1)
            } else {
                (t1, t2)
            };
            tmin = tmin.max(t1);
            tmax = tmax.min(t2);
            if tmin > tmax {
                return false;
            }
        }
        if ray.d.z == 0.0 && (ray.o.z < self.min.z || ray.o.z > self.max.z) {
            return false;
        } else {
            let ood = 1.0 / ray.d.z;
            let t1 = (self.min.z - ray.o.z) * ood;
            let t2 = (self.max.z - ray.o.z) * ood;
            let (t1, t2) = if t1 > t2 {
                (t2, t1)
            } else {
                (t1, t2)
            };
            tmin = tmin.max(t1);
            tmax = tmax.min(t2);
            if tmin > tmax {
                return false;
            }
        }
        true
    }
}


pub trait Object {
    fn collides(&self, ray: &Ray) -> Option<Shade>;
    fn make_bounds(&self) -> AABB;
    fn sample(&self) -> (Point3<f32>, Vector3<f32>);
}

#[derive(Copy, Clone)]
pub struct Plane {
    pub p: Point3<f32>,
    pub n: Vector3<f32>,
    pub mat: usize,
}

impl Object for Plane {
    fn collides(&self, ray: &Ray) -> Option<Shade> {
        let t = (self.p - ray.o).dot(self.n) / ray.d.dot(self.n);
        if t > 0.0 {
            Some(Shade{
                local: ray.o + ray.d * t,
                normal: self.n,
                mat: self.mat,
                tmin: t,
                depth: 1,
            })
        } else {
            None
        }
    }

    fn make_bounds(&self) -> AABB {
        AABB {
            min: Point3::<f32>::new(-f32::INFINITY,
                                  -f32::INFINITY,
                                  -f32::INFINITY),
            max: Point3::<f32>::new(f32::INFINITY,
                                  f32::INFINITY,
                                  f32::INFINITY),
        }
    }

    fn sample(&self) -> (Point3<f32>, Vector3<f32>) {
        // We won't have area light planes
        (self.p, self.n)
    }
}

#[derive(Copy, Clone)]
pub struct Triangle {
    pub a: Point3<f32>,
    pub b: Point3<f32>,
    pub c: Point3<f32>,
    pub mat: usize,
}

const EPSILON: f32 = 0.0000001;

impl Object for Triangle {
    fn collides(&self, ray: &Ray) -> Option<Shade> {
        // We shouldn't be calculating the normal every time, but whatever.
        let ab = self.b - self.a;
        let ac = self.c - self.a;
        let n = ab.cross(ac).normalize();
        let t = (self.a - ray.o).dot(n) / (ray.d.dot(n));
        if t > EPSILON {
            let p = ray.o + ray.d * t;
            let v = p - self.a;
            let dot1 = ac.dot(ac);
            let dot2 = ac.dot(ab);
            let dot3 = ac.dot(v);
            let dot4 = ab.dot(ab);
            let dot5 = ab.dot(v);
            let invd = 1.0 / (dot1 * dot4 - dot2 * dot2);
            let u = (dot4 * dot3 - dot2 * dot5) * invd;
            let v = (dot1 * dot5 - dot2 * dot3) * invd;
            if u >= 0.0 && v >= 0.0 && (u + v) < 1.0 {
                Some(Shade{
                    local: p,
                    normal: n,
                    mat: self.mat,
                    tmin: t,
                    depth: 1,
                })
            } else {
                None
            }
        } else {
            None
        }
    }

    fn make_bounds(&self) -> AABB {
        AABB {
            min: Point3::<f32>::new(self.a.x.min(self.b.x.min(self.c.x)),
                                  self.a.y.min(self.b.y.min(self.c.y)),
                                  self.a.z.min(self.b.z.min(self.c.z))),
            max: Point3::<f32>::new(self.a.x.max(self.b.x.max(self.c.x)),
                                  self.a.y.max(self.b.y.max(self.c.y)),
                                  self.a.z.max(self.b.z.max(self.c.z))),
        }
    }

    fn sample(&self) -> (Point3<f32>, Vector3<f32>) {
        let range = rand::distributions::Range::new(0f32, 1.0);
        let mut rng = rand::thread_rng();
        let v1 = self.b - self.a;
        let v2 = self.c - self.a;
        let x = range.ind_sample(&mut rng);
        let y = range.ind_sample(&mut rng);
        let p = self.a.to_vec() + self.b.to_vec() * x + self.c.to_vec() * y;
        (Point3::from_vec(p), v2.cross(v1))
    }
}

#[derive(Copy, Clone)]
pub struct Sphere {
    pub c: Point3<f32>,
    pub r: f32,
    pub mat: usize,
}

impl Object for Sphere {
    fn collides(&self, ray: &Ray) -> Option<Shade> {
        let temp = ray.o - self.c;
        let a = ray.d.dot(ray.d);
        let b = 2.0 * temp.dot(ray.d);
        let c = temp.dot(temp) - self.r * self.r;
        let d = b * b - 4.0 * a * c;
        if d < 0.0 {
            None
        } else {
            let e: f32 = d.sqrt();
            let denom = 2.0 * a;
            let t = (-b - e) / denom;
            if t > EPSILON {
                let l = ray.o + ray.d * t;
                let n = (l - self.c).normalize();
                Some(Shade{
                    local: l,
                    normal: n,
                    mat: self.mat,
                    tmin: t,
                    depth: 1,
                })
            } else {
                let t = (-b + e) / denom;
                let l = ray.o + ray.d * t;
                let n = (l - self.c).normalize();
                if t > EPSILON {
                    Some(Shade{
                        local: l,
                        normal: n,
                        mat: self.mat,
                        tmin: t,
                        depth: 1,
                    })
                } else {
                    None
                }
            }
        }
    }

    fn make_bounds(&self) -> AABB {
        AABB {
            min: Point3::<f32>::new(self.c.x - self.r,
                                  self.c.y - self.r,
                                  self.c.z - self.r),
            max: Point3::<f32>::new(self.c.x + self.r,
                                  self.c.y + self.r,
                                  self.c.z + self.r),
        }
    }

    fn sample(&self) -> (Point3<f32>, Vector3<f32>) {
        let range = rand::distributions::Range::new(0f32, 1.0);
        let mut rng = rand::thread_rng();
        let phi = 2.0 * f32::consts::PI * range.ind_sample(&mut rng);
        let theta = (2.0 * range.ind_sample(&mut rng) - 1.0).acos();
        let nv = Vector3::<f32>::new(self.r * theta.cos() * phi.sin(),
                                     self.r * theta.sin() * phi.cos(),
                                     self.r * phi.cos());
        let np = Point3::<f32>::new(nv.x + self.c.x,
                                    nv.y + self.c.y,
                                    nv.z + self.c.z);
        (np, nv.normalize())
    }
}

pub struct Instance {
    pub geom: Box<Object>,
    pub inv_transform: Matrix3<f32>,
}

impl Object for Instance {
    fn collides(&self, ray: &Ray) -> Option<Shade> {
        let inv_ray = Ray{
            o: Point3::from_vec(self.inv_transform * ray.o.to_vec()),
            d: self.inv_transform * ray.d,
        };
        if let Some(s) = self.geom.collides(&inv_ray) {
            Some(Shade{
                local: ray.o + ray.d * s.tmin,
                normal: (self.inv_transform * s.normal).normalize(),
                mat: s.mat,
                tmin: s.tmin,
                depth: 1,
            })
        } else {
            None
        }
    }

    fn make_bounds(&self) -> AABB {
        let bounds = self.geom.make_bounds();
        if let Some(trans) = self.inv_transform.invert() {
            AABB{
                min: Point3::from_vec(trans * bounds.min.to_vec()),
                max: Point3::from_vec(trans * bounds.max.to_vec()),
            }
        } else {
            bounds
        }
    }

    fn sample(&self) -> (Point3<f32>, Vector3<f32>) {
        self.geom.sample()
    }
}
