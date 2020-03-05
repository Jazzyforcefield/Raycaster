// Copyright 2020 Michael Ung

#ifndef VECTORTYPE_H_
#define VECTORTYPE_H_

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

  bool operator==(VectorType other) {
    return (x == other.x && y == other.y && z == other.z);
  }

  void Print() {
    std::cout << "X: " << x << " Y: " << y << " Z: " << z << std::endl;
  }

  float length() {
    return sqrt(x * x + y * y + z * z);
  }

  float x;
  float y;
  float z;
};

#endif  // VECTORTYPE_H_
