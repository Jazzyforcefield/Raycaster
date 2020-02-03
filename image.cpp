// Copyright 2020 Michael Ung

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

int image_write(std::ostream & fout, int width, int height) {
  int size = width * height;
  int r = 0, g = 0, b = 0;

  float rad = 0;

  fout << "P3" << std::endl;
  fout << "# A masterful piece of art" << std::endl;
  fout << width << " " << height << std::endl;
  fout << "255" << std::endl;

  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      int res = i > height * 0.05 && i < height * 0.95 && j > width * 0.05 && j < width * 0.95;
      r = 100 + 155 * ((float) (i * j) / size);
      g = 100 + 155 * ((float) (i * j) / (4 * size));
      b = 100 + 155 * ((float) (i * j) / (8 * size));
      
      if (abs((height / 4) * sin(j * 3.14 / (width / 2)) + (height / 2) - i) < 5) {
        fout << r << " " << 255 - g << " " << "0" << std::endl;
      } else {
        if (i > width / 2) {
          if (j > height / 2) {
            fout << 255 - b << " " << 255 - g << " " << 255 - r << std::endl;            
          } else {
            fout << b << " " << g << " " << r << std::endl;
          }
        } else {
          if (j < height / 2) {
            fout << 255 - b << " " << 255 - g << " " << 255 - r << std::endl;            
          } else {
            fout << b << " " << g << " " << r << std::endl;
          }
        }
      }
      rad = rad;
    }
  }
}

int main(int argc, char ** argv) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " " << "input_file" << std::endl;
    exit(0);
  }

  std::ifstream fin;
  std::ofstream fout;
  std::string key;
  char buf[50];
  int width, height;

  fin.open(argv[1]);
  if (!fin) {
    std::cerr << "Invalid input file!" << std::endl;
    exit(1);
  }

  fin >> key;
  fin >> width;
  fin >> height;

  if (key.compare("imsize") != 0) {
    std::cerr << "Invalid file contents!" << std::endl;
    exit(1);
  }

  if (width <= 0 || height <= 0) {
    std::cerr << "Invalid dimensions!" << std::endl;
    exit(1);
  }

  sprintf(buf, "%s.ppm", argv[1]);

  fout.open(buf);

  image_write(fout, width, height);

  fin.close();
  fout.close();
}
