extern crate cgmath;

use std::f32;
use std::num;
use std::ops::*;
use std::vec::Vec;
use std::io::{ Result, BufRead, BufReader };
use std::fs::File;
use cgmath::{Point, Vector, Transform, EuclideanVector, Decomposed, Quaternion,
             Point3, Vector3};

use geom::Triangle;

pub fn load_mesh(file: &str, scale: f32) -> Vec<Triangle> {
    let f = File::open(file).unwrap();
    let file = BufReader::new(&f);
    let mut nverts = 0;
    let mut nfaces = 0;
    let mut lines = file.lines();
    while let Some(line) = lines.next() {
        let curr_line = line.unwrap();
        let words: Vec<&str> = curr_line.split_whitespace().collect();
        match words[0] {
            // We ignore all properties, comments, and format directives.
            // We assume that vertices are three floats in x y z order.
            // We assume that faces are lists of three vertices.
            "element" => {
                match words[1] {
                    "vertex" => {
                        nverts = words[2].parse::<usize>().unwrap();
                    },
                    "face" => {
                        nfaces = words[2].parse::<usize>().unwrap();
                    },
                    _ => ()
                };
            },
            "end_header" => {
                break;
            },
            // Everything else we ignore.
            _  => ()
        }
    }
    // We also assume that vertices come before faces :-P
    let mut verts = Vec::<Point3<f32>>::new();
    for i in 0..nverts {
        if let Some(l) = lines.next() {
            let curr_line = l.unwrap();
            let words: Vec<&str> = curr_line.split_whitespace().collect();
            let p = Point3::<f32>::new(words[0].parse::<f32>().unwrap(),
                                       words[1].parse::<f32>().unwrap(),
                                       words[2].parse::<f32>().unwrap());
            verts.push(p * scale);//transform.transform_point(p));
        }
    }
    let mut faces = Vec::<Triangle>::new();
    for i in 0..nfaces {
        if let Some(l) = lines.next() {
            let curr_line = l.unwrap();
            let words: Vec<&str> = curr_line.split_whitespace().collect();
            let i1 = words[1].parse::<usize>().unwrap();
            let i2 = words[2].parse::<usize>().unwrap();
            let i3 = words[3].parse::<usize>().unwrap();
            let p1 = verts[i1];
            let p2 = verts[i2];
            let p3 = verts[i3];
            faces.push(Triangle{
                a: p1,
                b: p2,
                c: p3,
                mat: 0,
            });
        }
    }
    faces
}
