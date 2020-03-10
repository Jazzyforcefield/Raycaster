// Copyright 2020 Michael Ung

#ifndef TRIANGLETYPE_H_
#define TRIANGLETYPE_H_

#include "ray.h"
#include "shape.h"
#include "vectortype.h"

extern std::vector<VectorType> vb;
extern std::vector<VectorType> normals;
extern std::vector<VectorType> tex_coords;

class TriangleType : public ShapeType {
 public:
  TriangleType() {
    vertices_ = new unsigned int[3];
    normals_ = new unsigned int[3];
    tex_coords_ = new unsigned int[3];
    m_ = mtlc - 1;
    t_ = texc - 1;
  }

  TriangleType(unsigned int one, unsigned int two, unsigned int three,
               unsigned int n1 = 0, unsigned int n2 = 0, unsigned int n3 = 0,
               unsigned int t1 = 0, unsigned int t2 = 0, unsigned int t3 = 0) {
    vertices_ = new unsigned int[3];
    normals_ = new unsigned int[3];
    tex_coords_ = new unsigned int[3];

    vertices_[0] = one;
    vertices_[1] = two;
    vertices_[2] = three;

    normals_[0] = n1;
    normals_[1] = n2;
    normals_[2] = n3;

    tex_coords_[0] = n1;
    tex_coords_[1] = n2;
    tex_coords_[2] = n3;

    m_ = mtlc - 1;
    t_ = texc - 1;
  }

  ~TriangleType() { delete vertices_; }


  bool check_point(VectorType intersection, float & alpha, float & beta, float & gamma) {
    VectorType e1, e2, e3;
    float a, b, c;
    float sum;
    float actual;

    actual = (vb[vertices_[1] - 1] - vb[vertices_[0] - 1])
                   .cross(vb[vertices_[2] - 1] - vb[vertices_[0] - 1]).length() / 2.f;

    e1 = vb[vertices_[0] - 1] - intersection;
    e2 = vb[vertices_[1] - 1] - intersection;
    e3 = vb[vertices_[2] - 1] - intersection;

    c = e1.cross(e2).length();
    a = e2.cross(e3).length();
    b = e3.cross(e1).length();

    sum = (a + b + c) / 2.f;

    if (fabs(actual - sum) < 0.0001) {
      alpha = 0.5f * a / actual;
      beta = 0.5f * b / actual;
      gamma = 0.5f * c / actual;

      return true;
    } return false;
  }

  VectorType normal(VectorType bay = VectorType(1.f, 1.f, 1.f)) {
    if (n_def_ && !(bay.x == bay.y && bay.y == bay.z && bay.z == 1.f)) {
      return (normals[normals_[0] - 1].scalar(bay.x) +
              normals[normals_[1] - 1].scalar(bay.y) +
              normals[normals_[2] - 1].scalar(bay.z)).normalize();
    } else {
      return (vb[vertices_[1] - 1] - vb[vertices_[0] - 1]).cross(vb[vertices_[2] - 1] -
              vb[vertices_[0] - 1]).normalize();
    }
  }

  VectorType texture(VectorType bay) {
    if (textured_ && !(bay.x == bay.y && bay.y == bay.z && bay.z == 1.f)) {    
  return (tex_coords[tex_coords_[0] - 1].scalar(bay.x) +
              tex_coords[tex_coords_[1] - 1].scalar(bay.y) +
              tex_coords[tex_coords_[2] - 1].scalar(bay.z));
    } else {
      return VectorType(0, 0, 0);
    }
  }

  VectorType position() {
    return (vb[vertices_[0] - 1] + vb[vertices_[1] - 1] + vb[vertices_[2] - 1])
           .scalar(1.f / 3);
  }

  unsigned int * vertices_;
  unsigned int * normals_;
  unsigned int * tex_coords_;
  int m_;
  int t_;
  bool textured_;
  bool n_def_;
};

#endif  // TRIANGLETYPE_H_
