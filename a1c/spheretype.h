// Copyright 2020 Michael Ung

#ifndef SPHERETYPE_H_
#define SPHERETYPE_H_

#include "shape.h"

int mtlc = 0;
int texc = 0;

class SphereType : public ShapeType {
 public:
  SphereType(float x = 0, float y = 0, float z = 0, float r = 0) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->r = r;

    if (mtlc == 0) {
      this->m = 0;
    } else {
      this->m = mtlc - 1;
    }

    textured = false;
    t = texc - 1;
  }

  float x;
  float y;
  float z;
  float r;
  int m;
  int t;
  bool textured;
};

#endif  // SPHERETYPE_H_
