// Copyright 2020 Michael Ung

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include "ray.h"

#define MAX_SHAPES 20

// File related variables
std::ifstream fin;
std::ofstream fout;
std::string name;
std::string buffer;


// Initialize defaults and ifstream
int setup(std::ifstream & fin, char * name) {
  fin.open(name);

  if (!fin) {
    return -1;
  }

  hfov = 60;
  width = 800;
  height = 600;

  return 0;
}

// Aaaa yikes, reads file
int read_file(std::ifstream & fin) {
  while (fin >> buffer) {
    if (buffer == "v") {
      VectorType temp_vector;
      fin >> temp_vector.x;
      fin >> temp_vector.y;
      fin >> temp_vector.z;
      vb.push_back(temp_vector);
    } else if (buffer == "f") {
      TriangleType * temp_tri = new TriangleType();
      std::string b1;
      std::string b2;
      for (int i = 0; i < 3; i++) {
        fin >> b1;
        if (std::count(b1.begin(), b1.end(), '/') == 1) {
          b2 = b1.substr(0, b1.find("/"));
          temp_tri->vertices_[i] = stoi(b2);
          b2 = b1.substr(b1.find("/") + 1, b1.length());
          temp_tri->tex_coords_[i] = stoi(b2);
          temp_tri->n_def_ = false;
          temp_tri->textured_ = true;
        } else if (std::count(b1.begin(), b1.end(), '/') == 2) {
          if (b1.find("//") == std::string::npos) {
            b2 = b1.substr(0, b1.find("/"));
            temp_tri->vertices_[i] = stoi(b2);
            b2 = b1.substr(b1.find("/") + 1, b1.length() - b1.find("/"));
            temp_tri->tex_coords_[i] = stoi(b2);
            b2 = b1.substr(b1.find("/") + 1, b1.length());
            temp_tri->normals_[i] = stoi(b2);
            temp_tri->n_def_ = true;
            temp_tri->textured_ = true;
          } else {
            b2 = b1.substr(0, b1.find("//"));
            temp_tri->vertices_[i] = stoi(b2);
            b2 = b1.substr(b1.find("//") + 2, b1.length());
            temp_tri->normals_[i] = stoi(b2);
            temp_tri->n_def_ = true;
          }
        } else {
          temp_tri->vertices_[i] = stoi(b1);
          temp_tri->n_def_ = false;
        }
      }
      tri.push_back(temp_tri);

    } else if (buffer == "vn") {
      float xn, yn, zn;
      fin >> xn;
      fin >> yn;
      fin >> zn;
      normals.push_back(VectorType(xn, yn, zn));

    } else if (buffer == "vt") {
      float xt, yt;
      fin >> xt;
      fin >> yt;
      tex_coords.push_back(VectorType(xt, yt, 0));

    } else if (buffer == "eye") {
      if (!fin.eof())
        fin >> eye.x;
      else return -1;
      if (!fin.eof())
        fin >> eye.y;
      else return -1;
      if (!fin.eof())
        fin >> eye.z;
      else return -1;

    } else if (buffer == "viewdir") {
      if (!fin.eof())
        fin >> viewdir.x;
      else return -1;
      if (!fin.eof())
        fin >> viewdir.y;
      else return -1;
      if (!fin.eof())
        fin >> viewdir.z;
      else return -1;

    } else if (buffer == "updir") {
      if (!fin.eof())
        fin >> updir.x;
      else return -1;
      if (!fin.eof())
        fin >> updir.y;
      else return -1;
      if (!fin.eof())
        fin >> updir.z;
      else return -1;

    } else if (buffer == "hfov") {
      if (!fin.eof())
        fin >> hfov;
      else return -1;

    } else if (buffer == "imsize") {
      float w, h;
      if (!fin.eof()) {
        fin >> w;
        width = w;
      } else return -1;
      if (!fin.eof()) {
        fin >> h;
        height = h;
      } else return -1;

    } else if (buffer == "bkgcolor") {
      if (!fin.eof())
        fin >> bkgcolor.r;
      else return -1;
      if (!fin.eof())
        fin >> bkgcolor.g;
      else return -1;
      if (!fin.eof())
        fin >> bkgcolor.b;
      else return -1;

    } else if (buffer == "mtlcolor") {
      if (mtlc > 9) {
        std::cout << "Max materials reached!" << std::endl;
        exit(1);
      }
      if (!fin.eof())
        fin >> mtlcolor[mtlc].albedo_.r;
      else return -1;
      if (!fin.eof())
        fin >> mtlcolor[mtlc].albedo_.g;
      else return -1;
      if (!fin.eof())
        fin >> mtlcolor[mtlc].albedo_.b;
      else return -1;
      if (!fin.eof())
        fin >> mtlcolor[mtlc].highlight_.r;
      else return -1;
      if (!fin.eof())
        fin >> mtlcolor[mtlc].highlight_.g;
      else return -1;
      if (!fin.eof())
        fin >> mtlcolor[mtlc].highlight_.b;
      else return -1;
      if (!fin.eof())
        fin >> mtlcolor[mtlc].ambient_;
      else return -1;
      if (!fin.eof())
        fin >> mtlcolor[mtlc].diffuse_;
      else return -1;
      if (!fin.eof())
        fin >> mtlcolor[mtlc].specular_;
      else return -1;
      if (!fin.eof()) {
        fin >> mtlcolor[mtlc].n_;
        mtlc++;
      } else return -1;

    } else if (buffer == "texture") {
      if (texc > 9) {
        std::cout << "Max textures reached!" << std::endl;
        exit(1);
      }
      std::ifstream tfin;
      fin >> buffer;
      tfin.open(buffer);

      if (!tfin) {
        std::cerr << "Invalid .ppm file!" << std::endl;
        exit(1);
      }
      float max_c;
      tfin >> buffer;
      tfin >> textures[texc].width_;
      tfin >> textures[texc].height_;
      tfin >> max_c;
      textures[texc].map_ = new ColorType * [textures[texc].height_];
      for (int i = 0; i < textures[texc].height_; i++) {
        textures[texc].map_[i] = new ColorType[textures[texc].width_];
      }

      for (int row = 0; row < textures[texc].height_; row++) {
        for (int col = 0; col < textures[texc].width_; col++) {
          tfin >> textures[texc].map_[row][col].r;
          tfin >> textures[texc].map_[row][col].g;
          tfin >> textures[texc].map_[row][col].b;
          textures[texc].map_[row][col] = textures[texc].map_[row][col].scalar(1.f / max_c);
          //textures[texc].map_[row][col].Print();
        }
      }
      tfin.close();
      texc++;

    } else if (buffer == "light") {
      LightType * temp = new LightType();
      if (!fin.eof())
        fin >> temp->components_.x;
      else return -1;
      if (!fin.eof())
        fin >> temp->components_.y;
      else return -1;
      if (!fin.eof())
        fin >> temp->components_.z;
      else return -1;
      if (!fin.eof())
        fin >> temp->type_;
      else return -1;
      if (!fin.eof())
        fin >> temp->color_.r;
      else return -1;
      if (!fin.eof())
        fin >> temp->color_.g;
      else return -1;
      if (!fin.eof())
        fin >> temp->color_.b;
      else return -1;
      lights.push_back(temp);

    } else if (buffer == "attlight") {
      for (int i = 0; i < 7; i++) {
        if (!fin.eof()) {
          fin >> buffer;
        } else return -1;
      }
      if (!fin.eof())
      fin >> buffer;
      else return -1;
      if (!fin.eof())
      fin >> attenuation;
      else return -1;
      if (!fin.eof())
      fin >> buffer;
      else return -1;

    } else if (buffer == "depthcueing") {
      if (!fin.eof())
      fin >> depthcue.r;
      else return -1;
      if (!fin.eof())
      fin >> depthcue.g;
      else return -1;
      if (!fin.eof())
      fin >> depthcue.b;
      else return -1;
      if (!fin.eof())
      fin >> dcmax;
      else return -1;
      if (!fin.eof())
      fin >> dcmin;
      else return -1;
      if (!fin.eof())
      fin >> dc_far;
      else return -1;
      if (!fin.eof())
      fin >> dc_near;
      else return -1;

    } else if (buffer == "sphere") {
      SphereType * temp = new SphereType();
      if (!fin.eof())
        fin >> temp->x;
      else return -1;
      if (!fin.eof())
        fin >> temp->y;
      else return -1;
      if (!fin.eof())
        fin >> temp->z;
      else return -1;
      if (!fin.eof())
        fin >> temp->r;
      else return -1;
      objects.push_back(temp);

    } else if (buffer == "triangle") {
      std::cout << "Triangle not implemented" << std::endl;

    } else if (buffer == "cylinder") {
      std::cout << "Cylinder not implemented" << std::endl;
      for (int i = 0; i < 7; i++) {
        if (!fin.eof())
          fin >> buffer;
        else return -1;
      }

    } else {
      std::cout << "Error on: " << buffer << std::endl;
      return -1;

    }
  } return 0;
}

// Prints values of inputs for debugging
int print_values() {
  std::cout << "eye: " << eye.x << " " << eye.y << " " << eye.z << std::endl;
  std::cout << "viewdir: " << viewdir.x << " " << viewdir.y << " " << viewdir.z << std::endl;
  std::cout << "updir: " << updir.x << " " << updir.y << " " << updir.z << std::endl;
  std::cout << "hfov: " << hfov << std::endl;
  std::cout << "width: " << width << " height: " << height << std::endl;

  std::cout << "bkgcolor: " << bkgcolor.r << " " << bkgcolor.g << " " << bkgcolor.b;
  std::cout << std::endl;

  for (int i = 0; i < mtlc; i++) {
    std::cout << "mtlcolor: " << mtlcolor[i].albedo_.r << " " << mtlcolor[i].albedo_.g;
    std::cout << " " << mtlcolor[i].albedo_.b << " " << mtlcolor[i].highlight_.r;
    std::cout << " " << mtlcolor[i].highlight_.g << " ";
    std::cout << mtlcolor[i].highlight_.b << " " << mtlcolor[i].n_;
    std::cout << std::endl;
  }

  for (int i = 0; i < objects.size(); i++) {
    std::cout << "sphere: " << objects[i]->x << " " << objects[i]->y << " " << objects[i]->z;
    std::cout << " " << objects[i]->r << std::endl;
  }

  for (int i = 0; i < lights.size(); i++) {
    std::cout << "light: type_: " << lights[i]->type_ << std::endl;
    lights[i]->components_.Print();
    lights[i]->color_.Print();
    std::cout << std::endl;
  }

  for (int i = 0; i < vb.size(); i++) {
    std::cout << "v: " << vb[i].x << " " << vb[i].y << " " << vb[i].z;
    std::cout << std::endl;
  }

  for (int i = 0; i < normals.size(); i++) {
    std::cout << "vn: " << normals[i].x << " " << normals[i].y << " " << normals[i].z;
    std::cout << std::endl;
  }

  for (int i = 0; i < tex_coords.size(); i++) {
    std::cout << "vt: " << tex_coords[i].x << " " << tex_coords[i].y << " " << tex_coords[i].z;
    std::cout << std::endl;
  }

  for (int i = 0; i < tri.size(); i++) {
    std::cout << "f: " << tri[i]->vertices_[0] << " " << tri[i]->vertices_[1];
    std::cout << " " << tri[i]->vertices_[2];
    std::cout << std::endl;
  }

  return 0;
}

// Constrains inputs to a max
int check_values() {
  int status = 0;
  if (eye.x > 10000 || eye.y > 10000 || eye.z > 10000) {
    std::cout << "Eye is too far." << std::endl;
    status = -1;
  }

  if (viewdir.x == 0 && viewdir.y == 0 && viewdir.z == 0) {
    std::cout << "Bad viewdir." << std::endl;
    status = -1;
  }

  if (updir.x == 0 && updir.y == 0 && updir.z == 0) {
    std::cout << "Bad updir." << std::endl;
    status = -1;
  }

  if (updir.x == -viewdir.x && updir.y == -viewdir.y && updir.z == -viewdir.z) {
    std::cout << "Updir opposite of viewdir." << std::endl;
    status = -1;
  }

  if (height > 5000 && width > 5000) {
    std::cout << "Image too big, limit yourself to 5000x5000." << std::endl;
    status = -1;
  }

  if (height <= 0 || width <= 0) {
    std::cout << "Bad dimensions." << std::endl;
    status = -1;
  }

  return status;
}

int create_image(std::ofstream& fout) {
  fout << "P3" << std::endl;
  fout << "# A masterful piece of art" << std::endl;
  fout << width << " " << height << std::endl;
  fout << "255" << std::endl;

  for (int row = 0; row < height; row++) {
    for (int col = 0; col < width; col++) {
      fout << (int)(255 * pixels[row][col].r)  << " ";
      fout << (int)(255 * pixels[row][col].g) << " ";
      fout << (int)(255 * pixels[row][col].b) << std::endl;
    }
  }
  return 0;
}

void Update_Normals() {
  std::vector<VectorType> normals_new;
  for (int i = 0; i < normals.size(); i++) {
    int usage = 0;
    VectorType usagev = VectorType();

    for (int j = 0; j < tri.size(); j++) {
      if (tri[j]->normals_[0] == i + 1 || tri[j]->normals_[1] == i + 1 ||
          tri[j]->normals_[2] == i + 1) {
        usage++;
        usagev = usagev + tri[j]->normal();
      }
    }

    if (usage != 0) {
      normals_new.push_back(usagev.normalize());
    }
  }

  normals = normals_new;
}

int main(int argc, char ** argv) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " input_file" << std::endl;
    exit(0);
  }

  int ret = setup(fin, argv[1]);
  if (ret < 0) {
    std::cerr << "Invalid file!" << std::endl;
    exit(1);
  }

  ret = read_file(fin);
  if (ret < 0) {
    std::cerr << "Invalid file contents!" << std::endl;
    exit(1);
  }

  print_values();
  //ret = check_values();
  if (ret < 0) {
  //  std::cerr << "Please fix your errors." << std::endl;
   // exit(1);
  }
  
  name = argv[1];
  name = name.substr(0, name.length() - 4) + ".ppm";
  pixels = new ColorType * [height];
  for (int i = 0; i < height; i++) {  
    pixels[i] = new ColorType[width];
  }

  // Rays
  VectorType u = viewdir.cross(updir).normalize();
  VectorType v = u.cross(viewdir).normalize();
  VectorType n = viewdir.normalize();
  float focal = 10.0;


  // All points
  VectorType center = eye;
  aspect = (float)width / height;
  float swidth = 2 * focal * tan(hfov * 3.14 / 360);
  float sheight = swidth / aspect;
  VectorType un = u.scalar((float)swidth / 2);
  VectorType vn = v.scalar((float)sheight / 2);

  VectorType ul = VectorType(center.x + n.scalar(focal).x - un.x + vn.x,
                             center.y + n.scalar(focal).y - un.y + vn.y,
                             center.z + n.scalar(focal).z - un.z + vn.z);

  VectorType ur = VectorType(center.x + n.scalar(focal).x + un.x + vn.x,
                             center.y + n.scalar(focal).y + un.y + vn.y,
                             center.z + n.scalar(focal).z + un.z + vn.z);

  VectorType ll = VectorType(center.x + n.scalar(focal).x - un.x - vn.x,
                             center.y + n.scalar(focal).y - un.y - vn.y,
                             center.z + n.scalar(focal).z - un.z - vn.z);

  VectorType lr = VectorType(center.x + n.scalar(focal).x + un.x - vn.x,
                             center.y + n.scalar(focal).y + un.y - vn.y,
                             center.z + n.scalar(focal).z + un.z - vn.z);

  //Update_Normals();   // ????? I commentted this out and it worked...

  for (int row = 0; row < height; row++) {
    for (int col = 0; col < width; col++) {

      VectorType dh = (ur - ul).scalar((float)1.0 / (float)width);
      VectorType dv = (ll - ul).scalar((float)1.0 / (float)height);

      VectorType ch = dh.scalar(0.5); //(ur - ul).scalar(2.f * width);
      VectorType cv = dh.scalar(0.5); //(ll - ul).scalar(2.f * height);
      VectorType point = ((ul + dh.scalar(col) + dv.scalar(row) + ch + cv) -
                          eye).normalize();
      RayType ray = RayType(eye.x, eye.y, eye.z, point.x, point.y, point.z);
      pixels[row][col] = Trace_Ray(ray);
    }
  }
  
  fout.open(name);
  create_image(fout);

  return 0;
}








