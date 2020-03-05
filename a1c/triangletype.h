// Copyright 2020 Michael Ung

#ifndef TRIANGLETYPE_H_
#define TRIANGLETYPE_H_

#include "vectortype.h"

class TriangleType {
 public:
  TriangleType() { vertices_ = new VectorType[3]; m_ = 0; }

  VectorType vertices_;
  int m_;
};

#endif  // TRIANGLETYPE_H_
