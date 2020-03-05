// Copyright 2020 Michael Ung

#include <algorithm>

int mtlc = 0;
int cc = 0;

class RayType {
 public:
  RayType(float x = 0, float y = 0, float z = 0, float dx = 0, float dy = 0, float dz = 0) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
  }

  float x;
  float y;
  float z;
  float dx;
  float dy;
  float dz;
};

class SphereType {
 public:
  SphereType(float x = 0, float y = 0, float z = 0, float r = 0) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->r = r;
    this->m = cc++;
    cc %= mtlc;
  }

  float x;
  float y;
  float z;
  float r;
  int m;
};

class ColorType {
 public:
  ColorType(float r = 0, float g = 0, float b = 0) {
    this->r = r;
    this->g = g;
    this->b = b;
  }

  float r;
  float g;
  float b;
};

class VectorType {
 public:
  VectorType(float x = 0, float y = 0, float z = 0) {
    this->x = x;
    this->y = y;
    this->z = z;
  }
  
  VectorType cross(VectorType other) {
    return VectorType(y * other.z - z * other.y, z * other.x - x * other.z,
                          x * other.y - y * other.x);
  }

  VectorType normalize() {
    float length = sqrt(x * x + y * y + z * z);
    return VectorType(x / length, y / length, z / length);
  }

  VectorType scalar(float s) {
    return VectorType(x * s, y * s, z * s);
  }

  float dot(VectorType other) {
    return x * other.x + y * other.y + z * other.z;
  }

  VectorType operator-(VectorType other) {
    return VectorType(x - other.x, y - other.y, z - other.z);
  }

  VectorType operator+(VectorType other) {
    return VectorType(x + other.x, y + other.y, z + other.z);
  }

  float x;
  float y;
  float z;
};

// Raycast related variables
VectorType eye;
VectorType viewdir;
VectorType updir;
float hfov;
int width;
int height;
ColorType bkgcolor;
ColorType mtlcolor[10];
ColorType ** pixels;
std::vector<SphereType *> objects;

float aspect;


ColorType Shade_Ray(float x, float y, float z, SphereType & s) {
  return ColorType(mtlcolor[s.m].r, mtlcolor[s.m].g, mtlcolor[s.m].b);
}

ColorType Trace_Ray(RayType ray) {
  float B;
  float C;
  int sindex = -1;
  float mint = 99999;

  for (int s = 0; s < objects.size(); s++) {
    B = 2 * (ray.dx * (eye.x - objects[s]->x) + ray.dy * (eye.y - objects[s]->y) +
        ray.dz * (eye.z - objects[s]->z));
    C = pow(eye.x - objects[s]->x, 2) + pow(eye.y - objects[s]->y, 2) +
        pow(eye.z - objects[s]->z, 2) - pow(objects[s]->r, 2);
    float disc = pow(B, 2) - 4 * C;
    if (disc >= 0) {
      float t1 = (-B + sqrt(disc)) / 2;
      float t2 = (-B - sqrt(disc)) / 2;
      float ct = std::min(t1, t2);
      if (t1 < 0 || t2 < 0) {
        continue;
      }

      if (ct <= mint) {
        sindex = s;
      }
    } else {
      continue;
    }
  }

  if (sindex < 0) {
    return ColorType(0, 0, 0);
  } 

  VectorType dir = VectorType(ray.dx, ray.dy, ray.dz);
  VectorType intersect = eye + dir.scalar(mint);

  return Shade_Ray(intersect.x, intersect.y, intersect.z, *objects[sindex]);
}

