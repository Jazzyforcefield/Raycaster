// Copyright 2020 Michael Ung

#ifndef TEXTURETYPE_H_
#define TEXTURETYPE_H_

class TextureType {
 public:
  TextureType() {};
  TextureType(int width, int height) {
    map_ = new ColorType * [height];
    for (int i = 0; i < height; i++) {
      map_[i] = new ColorType[width];
    }
  }

  int width_;
  int height_;
  ColorType ** map_;
};

#endif  // TEXTURETYPE_H_
