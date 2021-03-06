extern crate cgmath;

use std::ops::*;
use cgmath::{Point, Vector, EuclideanVector, Point3, Vector3};

use geom::Ray;

#[derive(Copy, Clone)]
pub struct ViewPlane {
    pub w: u32,
    pub h: u32,
    pub s: f32,
    pub gamma: f32,
}


#[derive(Copy, Clone)]
pub struct OrthogonalProj {
    i: u32,
    j: u32,
    p: Point3<f32>,
    vp: ViewPlane,
}

pub trait Camera {
    fn dims(&self) -> (u32, u32);
}

impl OrthogonalProj {
    pub fn new(w: u32, h: u32, p: Point3<f32>) -> OrthogonalProj {
        OrthogonalProj{
            i: 0,
            j: 0,
            p: p,
            vp: ViewPlane{
                w: w,
                h: h,
                s: 0.2,
                gamma: 1.0,
            },
        }
    }
}

impl Camera for OrthogonalProj {
    fn dims(&self) -> (u32, u32) {
        (self.vp.w, self.vp.h)
    }
}

impl Iterator for OrthogonalProj {
    type Item = (u32, u32, Box<Fn(f32, f32) -> Ray>);

    fn next(&mut self) -> Option<(u32, u32, Box<Fn(f32, f32) -> Ray>)> {
        if self.j == self.vp.w {
            self.i += 1;
            if self.i == self.vp.h {
                return None;
            }
            self.j = 0;
        }
        let li = self.i;
        let lj = self.j;
        let lp = self.p;
        let lvp = self.vp;
        self.j += 1;
        Some((li, lj, Box::new(move |x: f32, y: f32| -> Ray {
            let xp = lvp.s * (lj as f32 - 0.5 * (lvp.w as f32 + x));
            let yp = lvp.s * (li as f32 - 0.5 * (lvp.h as f32 + y));
            Ray{
                o: Point3::<f32>::new(xp, yp, 0.0) + lp.to_vec(),
                d: Vector3::<f32>::new(0.0, 0.0, -1.0),
            }
        })))
    }
}


#[derive(Copy, Clone)]
pub struct PerspectiveProj {
    i: u32,
    j: u32,
    eye: Point3<f32>,
    u: Vector3<f32>,
    v: Vector3<f32>,
    w: Vector3<f32>,
    d: f32,
    vp: ViewPlane,
}

impl PerspectiveProj {
    pub fn new(width: u32, height: u32, d: f32, s: f32, eye: Point3<f32>, lookat: Point3<f32>, up: Vector3<f32>) -> PerspectiveProj {
        let w = (eye - lookat).normalize();
        let u = up.cross(w).normalize();
        let v = w.cross(u);
        PerspectiveProj{
            i: 0,
            j: 0,
            eye: eye,
            u: u,
            v: v,
            w: w,
            d: d,
            vp: ViewPlane{
                w: width,
                h: height,
                s: s,
                gamma: 1.0,
            },
        }
    }
}

impl Camera for PerspectiveProj {
    fn dims(&self) -> (u32, u32) {
        (self.vp.w, self.vp.h)
    }
}

impl Iterator for PerspectiveProj {
    type Item = (u32, u32, Box<Fn(f32, f32) -> Ray>);

    fn next(&mut self) -> Option<(u32, u32, Box<Fn(f32, f32) -> Ray>)> {
        if self.j == self.vp.w {
            self.i += 1;
            if self.i == self.vp.h {
                return None;
            }
            self.j = 0;
        }
        let li = self.i;
        let lj = self.j;
        let leye = self.eye;
        let lu = self.u;
        let lv = self.v;
        let lw = self.w;
        let ld = self.d;
        let lvp = self.vp;
        self.j += 1;
        Some((li, lj, Box::new(move |x: f32, y: f32| -> Ray {
            let xp = lvp.s * (lj as f32 - 0.5 * (lvp.w as f32 + x));
            let yp = lvp.s * (li as f32 - 0.5 * (lvp.h as f32 + y));
            let d = (lu * xp + lv * yp - lw * ld).normalize();
            Ray{
                o: leye,
                d: d,
            }
        })))
    }
}
