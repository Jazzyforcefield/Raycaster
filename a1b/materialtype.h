// Copyright 2020 Michael Ung

#ifndef MATERIALTYPE_H_
#define MATERIALTYPE_H_

#include "colortype.h"

class MaterialType {
 public:
  float ambient_, diffuse_, specular_;
  ColorType albedo_, highlight_;
  float n_;

  MaterialType(float odr = 0, float odg = 0, float odb = 0, float osr = 0,
               float osg = 0, float osb = 0, float ka = 0, float kd = 0,
               float ks = 0, float n = 0) {
    ambient_ = ka; diffuse_ = kd; specular_ = ks;
    albedo_ = ColorType(odr, odg, odb);
    highlight_ = ColorType(osr, osg, osb);
    n_ = n;
  }

};

#endif  // MATERIALTYPE_H_
