// Copyright 2020 Michael Ung

#ifndef LIGHTTYPE_H_
#define LIGHTTYPE_H_

class LightType {
 public:
  VectorType components_;
  int type_;
  ColorType color_;

  LightType(float x = 0, float y = 0, float z = 0, float w = 0, float r = 0,
            float g = 0, float b = 0) {
    components_ = VectorType(x, y, z);
    type_ = w;
    color_ = ColorType(r, g, b);
  }
};

#endif  // LIGHTTYPE_H_
