extern crate cgmath;

use std::cmp;
use std::vec::Vec;
use cgmath::{Point, Vector, Point3, Vector3};

use geom::{Object, AABB, Ray};
use shade::Shade;

#[derive(Copy, Clone)]
enum BVHNodeType<'a> {
    Leaf(&'a (Object + 'a)),
    Parent(usize, usize),
}

#[derive(Copy, Clone)]
struct BVHNode<'a> {
    bounds: AABB,
    height: i32,
    parent: usize,
    node_type: BVHNodeType<'a>,
}

pub struct BVH<'a> {
    root: usize,
    free: Vec<usize>,
    pool: Vec<BVHNode<'a>>,
}

impl<'a> BVH<'a> {
    pub fn new() -> BVH<'a>{
        BVH{
            root: 0,
            free: Vec::new(),
            pool: Vec::new(),
        }
    }

    pub fn empty(&self) -> bool {
        self.free.len() == self.pool.len()
    }

    fn allocate_node(&mut self) -> usize {
        if let Some(free_elem) = self.free.pop() {
            self.pool[free_elem].height = -1;
            free_elem
        } else {
            let elem = self.pool.len();
            self.pool.push(BVHNode{
                bounds: AABB::empty(),
                height: -1,
                parent: 0,
                node_type: BVHNodeType::Parent(0, 0),
            });
            elem
        }
    }

    pub fn remove_leaf(&mut self, leaf: usize) {
        self.free.push(leaf);
        let parent = self.pool[leaf].parent;
        if let BVHNodeType::Parent(child1, child2) = self.pool[parent].node_type {
            let sibling = if child1 == leaf {
                child2
            } else {
                child1
            };
            if self.root != parent {
                let grand_parent = self.pool[parent].parent;
                if let BVHNodeType::Parent(child1, child2) =
                    self.pool[grand_parent].node_type {
                        self.pool[grand_parent].node_type =
                            if child1 == parent {
                                BVHNodeType::Parent(sibling, child2)
                            } else {
                                BVHNodeType::Parent(child1, sibling)
                            }
                    }
                self.pool[sibling].parent = grand_parent;
                self.free.push(parent);
                let mut i = grand_parent;
                loop {
                    i = self.balance(i);
                    if let BVHNodeType::Parent(child1, child2)
                        = self.pool[i].node_type {
                            self.pool[i].bounds = AABB::combine(
                                self.pool[child1].bounds, self.pool[child2].bounds);
                            self.pool[i].height = 1 + cmp::max(
                                self.pool[child1].height,
                                self.pool[child2].height);
                            if self.root == i {
                                break;
                            }

                        i = self.pool[i].parent;
                    }
                }
            } else {
                self.root = sibling;
                self.free.push(parent);
            }
        }
    }

    // Call a given function f for every object ray collides with. If f returns
    // false, prematurely exit and return false. Otherwise continue with every
    // object and return true.
    pub fn raytrace<F>(&self, ray: &Ray, mut callback: F) -> bool
        where F: FnMut(usize, Shade) -> bool {
            if self.empty() {
                return true;
            }
            let mut stack = vec![ self.root ];
            while !stack.is_empty() {
                if let Some(top) = stack.pop() {
                    if self.pool[top].bounds.collides_ray(ray) {
                        match self.pool[top].node_type {
                            BVHNodeType::Leaf(obj) => {
                                if let Some(s) = obj.collides(ray) {
                                    if !callback(top, s) {
                                        return false;
                                    }
                                }
                            },
                            BVHNodeType::Parent(c1, c2) => {
                                stack.push(c1);
                                stack.push(c2);
                            }
                        }
                    }
                }
            }
            true
        }

    pub fn insert_object(&mut self, obj: &'a Object) -> usize {
        let leaf = self.allocate_node();
        let bounds = obj.make_bounds();
        self.pool[leaf].node_type = BVHNodeType::Leaf(obj);
        self.pool[leaf].bounds = bounds;
        if self.free.len() == self.pool.len() - 1 {
            self.root = leaf;
            return leaf;
        }
        let mut best = self.root;
        loop {
            if let BVHNodeType::Parent(child1, child2) = self.pool[best].node_type {
                let curr_bounds = self.pool[best].bounds;
                let area = curr_bounds.surface_area();
                let combined_bounds = AABB::combine(curr_bounds, bounds);
                let combined_area = combined_bounds.surface_area();
                let cost = combined_area * 2.0;
                let inheritance_cost = (combined_area - area) * 2.0;
                let child1_cost =
                    if let BVHNodeType::Parent(_, _) = self.pool[child1].node_type {
                        let combined = AABB::combine(bounds,
                                                     self.pool[child1].bounds);
                        let old_area = self.pool[child1].bounds.surface_area();
                        let new_area = combined.surface_area();
                        new_area - old_area + inheritance_cost
                    } else {
                        let combined = AABB::combine(bounds,
                                                     self.pool[child1].bounds);
                        combined.surface_area() + inheritance_cost
                    };
                let child2_cost =
                    if let BVHNodeType::Parent(_, _) = self.pool[child2].node_type {
                        let combined = AABB::combine(bounds,
                                                     self.pool[child2].bounds);
                        let old_area = self.pool[child2].bounds.surface_area();
                        let new_area = combined.surface_area();
                        new_area - old_area + inheritance_cost
                    } else {
                        let combined = AABB::combine(bounds,
                                                     self.pool[child2].bounds);
                        combined.surface_area() + inheritance_cost
                    };
                // Descend according to minimum cost
                if cost < child1_cost && cost < child2_cost {
                    break;
                }

                best = if child1_cost < child2_cost {
                    child1
                } else {
                    child2
                };
            } else {
                break;
            }
        }

        // Create a new parent
        let old_parent = self.pool[best].parent;
        let new_parent = self.allocate_node();
        self.pool[new_parent].parent = old_parent;
        self.pool[new_parent].bounds = AABB::combine(bounds,
                                                     self.pool[best].bounds);
        self.pool[new_parent].height = self.pool[best].height + 1;

        if best != self.root {
            if let BVHNodeType::Parent(child1, child2) = self.pool[old_parent].node_type {
                self.pool[old_parent].node_type = if child1 == best {
                    BVHNodeType::Parent(new_parent, child2)
                } else {
                    BVHNodeType::Parent(child1, new_parent)
                };
            }
        } else {
            self.root = new_parent;
        }
        self.pool[new_parent].node_type = BVHNodeType::Parent(best, leaf);
        self.pool[best].parent = new_parent;
        self.pool[leaf].parent = new_parent;

        // Walk up the tree fixing the heights and bounds.
        let mut i = self.pool[leaf].parent;
        loop {
            i = self.balance(i);

            // Gauranteed to be a parent.
            if let BVHNodeType::Parent(child1, child2) = self.pool[i].node_type {
                self.pool[i].height = 1 + cmp::max(self.pool[child1].height,
                                                   self.pool[child2].height);
                self.pool[i].bounds = AABB::combine(self.pool[child1].bounds,
                                                    self.pool[child2].bounds);
                if i == self.root {
                    break;
                };
            }

            i = self.pool[i].parent;
        }

        leaf
    }

    // This could be really cleaned up with using pointers instead of array
    // references everywhere.
    fn balance(&mut self, a: usize) -> usize {
        if self.pool[a].height < 2 {
            return a
        }
        if let BVHNodeType::Parent(b, c) = self.pool[a].node_type {
            if self.pool[c].height > self.pool[b].height {
                if let BVHNodeType::Parent(f, g) = self.pool[c].node_type {
                    // Swap A and C
                    self.pool[c].parent = self.pool[a].parent;
                    self.pool[a].parent = c;

                    // A's old parent should point to C
                    if self.root == a {
                        self.root = c;
                    } else if let BVHNodeType::Parent(pchild1, pchild2)
                        = self.pool[self.pool[c].parent].node_type {
                            let parent = self.pool[c].parent;
                            self.pool[parent].node_type =
                                if pchild1 == a {
                                    BVHNodeType::Parent(c, pchild2)
                                } else {
                                    BVHNodeType::Parent(pchild1, c)
                                };
                        }

                    // Rotate and readjust
                    // This really needs to be extracted into its own function
                    if self.pool[f].height > self.pool[g].height {
                        self.pool[c].node_type = BVHNodeType::Parent(a, f);
                        self.pool[a].node_type = BVHNodeType::Parent(b, g);
                        self.pool[g].parent = a;
                        self.pool[a].bounds = AABB::combine(
                            self.pool[b].bounds, self.pool[g].bounds);
                        self.pool[c].bounds = AABB::combine(
                            self.pool[a].bounds, self.pool[f].bounds);
                        self.pool[a].height = 1 + cmp::max(
                            self.pool[b].height, self.pool[g].height);
                        self.pool[c].height = 1 + cmp::max(
                            self.pool[a].height, self.pool[f].height);
                    } else {
                        self.pool[c].node_type = BVHNodeType::Parent(a, g);
                        self.pool[a].node_type = BVHNodeType::Parent(b, f);
                        self.pool[f].parent = a;
                        self.pool[a].bounds = AABB::combine(
                            self.pool[b].bounds, self.pool[f].bounds);
                        self.pool[c].bounds = AABB::combine(
                            self.pool[a].bounds, self.pool[g].bounds);
                        self.pool[a].height = 1 + cmp::max(
                            self.pool[b].height, self.pool[f].height);
                        self.pool[c].height = 1 + cmp::max(
                            self.pool[a].height, self.pool[g].height);
                    }
                }
                return c;
            }
            if self.pool[b].height > self.pool[c].height {
                if let BVHNodeType::Parent(d, e) = self.pool[b].node_type {
                    // Swap A and B
                    self.pool[b].parent = self.pool[a].parent;
                    self.pool[a].parent = b;

                    // A's old parent should point to B
                    if self.root == a {
                        self.root = b;
                    } else if let BVHNodeType::Parent(pchild1, pchild2)
                        = self.pool[self.pool[b].parent].node_type {
                            let parent = self.pool[b].parent;
                            self.pool[parent].node_type =
                                if pchild1 == a {
                                    BVHNodeType::Parent(b, pchild2)
                                } else {
                                    BVHNodeType::Parent(pchild1, b)
                                };
                        }

                    // Rotate and readjust
                    // This really needs to be extracted into its own function
                    if self.pool[d].height > self.pool[e].height {
                        self.pool[b].node_type = BVHNodeType::Parent(a, d);
                        self.pool[a].node_type = BVHNodeType::Parent(e, c);
                        self.pool[e].parent = a;
                        self.pool[a].bounds = AABB::combine(
                            self.pool[c].bounds, self.pool[e].bounds);
                        self.pool[b].bounds = AABB::combine(
                            self.pool[a].bounds, self.pool[d].bounds);
                        self.pool[a].height = 1 + cmp::max(
                            self.pool[c].height, self.pool[e].height);
                        self.pool[b].height = 1 + cmp::max(
                            self.pool[a].height, self.pool[d].height);
                    } else {
                        self.pool[b].node_type = BVHNodeType::Parent(a, e);
                        self.pool[a].node_type = BVHNodeType::Parent(d, c);
                        self.pool[d].parent = a;
                        self.pool[a].bounds = AABB::combine(
                            self.pool[c].bounds, self.pool[d].bounds);
                        self.pool[b].bounds = AABB::combine(
                            self.pool[a].bounds, self.pool[e].bounds);
                        self.pool[a].height = 1 + cmp::max(
                            self.pool[c].height, self.pool[d].height);
                        self.pool[b].height = 1 + cmp::max(
                            self.pool[a].height, self.pool[e].height);
                    }
                }
                return b
            }
        }
        a
    }
}
