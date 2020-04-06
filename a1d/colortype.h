// Copyright 2020 Michael Ung

#ifndef COLORTYPE_H_
#define COLORTYPE_H_

class ColorType {
 public:
  float r;
  float g;
  float b;  

  ColorType(float r = 0, float g = 0, float b = 0) {
    this->r = r;
    this->g = g;
    this->b = b;
  }

  ColorType scalar(float n) {
    return ColorType(r * n, g * n, b * n);
  }

  ColorType operator+(ColorType other) {
    return ColorType(r + other.r, g + other.g, b + other.b);
  }

  ColorType operator-(ColorType other) {
    return ColorType(r - other.r, g - other.g, b - other.b);
  }

  ColorType operator*(ColorType other) {
    return ColorType(r * other.r, g * other.g, b * other.b);
  }

  ColorType normalize() {
    return ColorType(r / 3.f, g / 3.f, b / 3.f);
  }

  ColorType clamp(float min, float max) {
    float red, green, blue;

    red = std::max(min, r);
    red = std::min(red, max);
    green = std::max(min, g);
    green = std::min(green, max);
    blue = std::max(min, b);
    blue = std::min(blue, max);

    return ColorType(red, green, blue);
  }

  void Print() {
    std::cout << "Red: " << r << " Green: " << g << " Blue: " << b << std::endl;
  }
};

#endif  // COLORTYPE_H_
