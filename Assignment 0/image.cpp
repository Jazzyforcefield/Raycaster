// Copyright 2020 Michael Ung

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

// Function that writes all contents to the .ppm file
int image_write(std::ostream & fout, int width, int height) {
  int size = width * height;
  int r = 0, g = 0, b = 0;

  float rad = 0;

  std::cout << "Creating image..." << std::endl;

  fout << "P3" << std::endl;
  fout << "# A masterful piece of art" << std::endl;
  fout << width << " " << height << std::endl;
  fout << "255" << std::endl;

  // Constraints to have a sine curve printed
  for (int i = 0; i < height; i++) {    // Loop through every row
    for (int j = 0; j < width; j++) {   // Loop through ever column
      // Varying colors
      r = 100 + 155 * ((float) (i * j) / size);
      g = 100 + 155 * ((float) (i * j) / (4 * size));
      b = 100 + 155 * ((float) (i * j) / (8 * size));
      
      // Checks if row is part of curve
      if (abs((height / 4) * sin(j * 3.14 / (width / 2)) + (height / 2) - i) < 5) {
        fout << r << " " << 255 - g << " " << "0" << std::endl;
      } else {
        if (i > height / 2) {
          if (j > width / 2) {
            fout << 255 - b << " " << 255 - g << " " << 255 - r << std::endl;            
          } else {
            fout << b << " " << g << " " << r << std::endl;
          }
        } else {
          if (j < width / 2) {
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

// Main function that handles file input
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
  if (!fin) {   // Checks for valid file
    std::cerr << "Invalid input file!" << std::endl;
    exit(1);
  }

  fin >> key;
  fin >> width;
  fin >> height;

  if (key.compare("imsize") != 0) {   // Checks for keyword
    std::cerr << "Invalid file contents!" << std::endl;
    exit(1);
  }

  if (width <= 0 || height <= 0) {    // Checks for proper dimensions
    std::cerr << "Invalid dimensions!" << std::endl;
    exit(1);
  }

  if  (width > 3840 && height > 2160) {   // Additional check because bigger = longer
    std::cerr << "Size is too big! I'm gonna have to limit you to 4k resolution." << std::endl; 
  }

  sprintf(buf, "%s.ppm", argv[1]);    // Output file name

  fout.open(buf);

  image_write(fout, width, height);   // Actually writing to file

  fin.close();
  fout.close();
}
