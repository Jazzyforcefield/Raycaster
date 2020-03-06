// Copyright 2020 Michael Ung

#ifndef TRIANGLETYPE_H_
#define TRIANGLETYPE_H_

#include "ray.h"
#include "shape.h"
#include "vectortype.h"

extern std::vector<VectorType> vb;

class TriangleType : public ShapeType {
 public:
  TriangleType() {
    vertices_ = new unsigned int[3];
    m_ = 0;
  }

  TriangleType(unsigned int one, unsigned int two, unsigned int three) {
    vertices_ = new unsigned int[3];

    vertices_[0] = one;
    vertices_[1] = two;
    vertices_[2] = three;

    m_ = 0;
  }

  ~TriangleType() { delete vertices_; }

  VectorType normal() {
    return (vb[vertices_[1] - 1] - vb[vertices_[0] - 1]).cross(vb[vertices_[2] - 1] -
            vb[vertices_[0] - 1]);
  }

  VectorType position() {
    return (vb[vertices_[0] - 1] + vb[vertices_[1] - 1] + vb[vertices_[2] - 1])
           .scalar(1.f / 3);
  }

  unsigned int * vertices_;
  int m_;
};

#endif  // TRIANGLETYPE_H_
