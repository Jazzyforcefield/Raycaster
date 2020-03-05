// Copyright 2020 Michael Ung

#ifndef SPHERETYPE_H_
#define SPHERETYPE_H_

int mtlc = 0;
int cc = 0;

class SphereType {
 public:
  SphereType(float x = 0, float y = 0, float z = 0, float r = 0) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->r = r;

    if (mtlc == 0) {
      this->m = 0;
    } else {
      this->m = cc++;
      cc %= mtlc;
    }
  }

  float x;
  float y;
  float z;
  float r;
  int m;
};

#endif  // SPHERETYPE_H_
