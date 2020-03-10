// Copyright 2020 Michael Ung

#ifndef RAY_H_
#define RAY_H_

#include <algorithm>
#include <iomanip>

#include "vectortype.h"
#include "colortype.h"
#include "materialtype.h"
#include "spheretype.h"
#include "triangletype.h"
#include "lighttype.h"
#include "texturetype.h"


// Raycast related variables
VectorType eye;
VectorType viewdir;
VectorType updir;
float hfov;
float aspect;
int width;
int height;
ColorType bkgcolor;
ColorType depthcue = ColorType(0, 0, 0);
ColorType ** pixels;
MaterialType mtlcolor[10];
TextureType textures[10];
std::vector<SphereType *> objects;

// Vertex buffer
std::vector<VectorType> vb;
std::vector<VectorType> normals;
std::vector<VectorType> tex_coords;

std::vector<TriangleType *> tri;

std::vector<LightType *> lights;
float attenuation = 1;

float dc_near = 10000.f;
float dc_far = 999999.f;
float dcmax = 1.f;
float dcmin = 0.f;

class RayType {
 public:
  RayType(float x = 0, float y = 0, float z = 0,
          float dx = 0, float dy = 0, float dz = 0) {
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

ColorType Shade_Ray(float x, float y, float z, SphereType & s) {
  float epsilon = 0.01;  // Sparse ray threshold
  float ka, kd, ks;      // Ambient, diffuse, and specular constants
  float n;               // Specular highlight fall-off
  float fshadow;
  ColorType od, os;      // Albedo and highlight color
  ColorType ambient;
  ColorType diffspec;
  ColorType result;
  VectorType N, V, L, H; // Intersection normal, surface-light direction, 
  VectorType Lbefore;

  VectorType intersection = VectorType(x, y, z);
  VectorType sphere_center = VectorType(s.x, s.y, s.z);

  N = (VectorType(x, y, z) - sphere_center).normalize();

  int u, v;
  float negative = atan2(N.y, N.x);
  VectorType t_coords;

  if (negative < 0.0001) {
    t_coords = VectorType((2 * 3.1415926 + atan2(N.y, N.x)) / (2 * 3.1415926),
                                   acos(N.z) / 3.1415926, 0);
  } else {
    t_coords = VectorType(atan2(N.y, N.x)/ (2 * 3.1415926),
                                   acos(N.z) / 3.1415926, 0);
  }

  if (s.textured_ && texc > 0) {
    u = t_coords.x * textures[s.t_].width_;
    v = t_coords.y * textures[s.t_].height_;
  }

  // Assigning equation variables
  ka = mtlcolor[s.m].ambient_;
  kd = mtlcolor[s.m].diffuse_;
  ks = mtlcolor[s.m].specular_;
  od = (s.textured_ && texc > 0) ? textures[s.t_].map_[v][u] : mtlcolor[s.m].albedo_;
  os = mtlcolor[s.m].highlight_;
  n = mtlcolor[s.m].n_;
  V = viewdir.scalar(-1).normalize();

  ambient = od.scalar(ka);

  // Loop through lights
  for (int i = 0; i < lights.size(); i++) {
    if (lights[i]->type_ == 1) {
      Lbefore = lights[i]->components_ - intersection;
    } else {
      Lbefore = lights[i]->components_.scalar(-1);
    }

    float distance;
    float B, C, discrim;
    float t1, t2, ct;

    distance = Lbefore.length();
    L = Lbefore.normalize();
    fshadow = 200.f;

    // Avoid shadow bleeding
    bool in_shadow = false;

    // Shadow calculations
    for (int j = 0; j < objects.size(); j++) {
      VectorType otherpos = VectorType(objects[j]->x, objects[j]->y, objects[j]->z);
      VectorType otherdir = (intersection - otherpos).normalize();
      float length = (intersection - otherpos - otherdir.scalar(objects[j]->r)).length();

      // Checks if the length between surfaces larger than light distance
      // Also checks if it's a directional light or if it's in a shadow
      // Or if it's in the same direction as the light
       //if (in_shadow) break;
      if (length >= distance && lights[i]->type_ != 0 ||
          otherdir.dot(L.scalar(-1)) <= 0 || in_shadow) {
        continue;
      }

      // Cast 50 rays for softer shadows
      for (int k = 0; k < 50; k++) {
        float random[3];

        // Max 20% deviation from original direction
        for (int r = 0; r < 3; r++) {
          random[r] = (rand() % 200) / 1000.f;
        }   // r

        VectorType rand_vec = VectorType(random[0], random[1], random[2]).normalize();
        if (rand() % 1) rand_vec = rand_vec.scalar(-1.f);


        // Calculates B, C, and the discriminant for the ray-sphere intersection equation
        B = 2 * (L.dot(intersection - otherpos));
        C = (intersection - otherpos).dot(intersection - otherpos) - pow(objects[j]->r, 2);
        discrim = pow(B, 2) - 4 * C;

        // If there are multiple solutions, calculate the values for t
        if (discrim >= 0) {
          t1 = (-B + sqrt(discrim)) / 2;
          t2 = (-B - sqrt(discrim)) / 2;
          ct = std::max(t1, t2);

          if (t1 >= epsilon && t2 >= epsilon) {
            ct = std::min(t1, t2);
          }

          // For when intersection is behind
          if (ct <= epsilon) {
            continue;           // Go to next ray
          } else if (ct > epsilon) {
            // If already calculated, then don't add more than the shadow for one object
            in_shadow = true;

            // Subtracting from the shadow factor for darker shadows, alterable
            // Controls some of shadow fall-off, scalable by distance
            fshadow -= otherdir.dot(N.scalar(-1)) * (lights.size() / (length + 1));

            // Change shadow behavior if a directional light
            // Side note: if only i had read the canvas page I would've saved
            // An entire day wondering why my shadows looked weird, and that
            // Was because my directional light fell off. Still don't know
            // Why it took so long to figure out tho...         
            if (lights[i]->type_ == 0) {
              fshadow = 0.f;
              break;
            }
          } else {
            std::cerr << "Problem has occurred in ShadeRay!" << std::endl;
            exit(1);
          }
        }
      }     // k
    }       // j

    // Set the shadow factor and the attentuation factors
    // The attenuation used only has linear scaling with distance
    // Increase numerator to increase contrast/decreases fall-off rate
    // Basically a = 0, b = 0.125, c = 0 | a = 1
    attenuation = (lights[i]->type_ == 1) ? 8.f / distance : 2.f;
    fshadow = std::min(std::max(0.f, fshadow / 200.f), 1.f);

    // Calculate H, which is the normalized sum of the light direction and the
    // Negative normalzed viewing direction
    H = (L + V).normalize();

    // A running total of the diffusion and specular components
    // Multiplied by the attenuation factor and the shadow factor
    diffspec = diffspec + lights[i]->color_
               .scalar(attenuation * fshadow * (1)) *
               (od.scalar(kd * std::max(0.f, N.dot(L))) +
               os.scalar(ks * pow(std::max(0.f, N.dot(H)), n)));
    fshadow = 200.f;
  }           // i

  // Reset shadow factor

  // Clamp the result of the addition between the ambient light and the other
  // Light components using a vector component clamp

  float dc_dist = (intersection - eye).length();
  float dca = 0.f;

  if (dc_dist >= dc_far) {
    dca = dcmin;
  } else if (dc_dist <= dc_near) {
    dca = dcmax;
  } else {
    dca = (dcmax - dcmin) * (dc_far - dc_dist) / (dc_far - dc_near);
  }

  result = (ambient + diffspec).clamp(0.f, 1.f);
  result = result.scalar(dca) + depthcue.scalar(1.f - dca);
  return result.clamp(0.f, 1.f);
}

ColorType Shade_RayT(float x, float y, float z, TriangleType & s, VectorType bay) {
  float epsilon = 0.01;  // Sparse ray threshold
  float ka, kd, ks;      // Ambient, diffuse, and specular constants
  float n;               // Specular highlight fall-off
  float fshadow;
  ColorType od, os;      // Albedo and highlight color
  ColorType ambient;
  ColorType diffspec;
  ColorType result;
  VectorType N, V, L, H; // Intersection normal, surface-light direction, 
  VectorType Lbefore;

  VectorType intersection = VectorType(x, y, z);

  int u, v;
  VectorType t_coords = s.texture(bay);

  // Check to see if texture is provided
  if (s.textured_ && texc > 0) {
    u = t_coords.x * textures[s.t_].width_;
    v = t_coords.y * textures[s.t_].height_;
  }

  // Assigning equation variables
  ka = mtlcolor[s.m_].ambient_;
  kd = mtlcolor[s.m_].diffuse_;
  ks = mtlcolor[s.m_].specular_;
  od = (s.textured_ && texc > 0) ? textures[s.t_].map_[v][u] : mtlcolor[s.m_].albedo_;
  os = mtlcolor[s.m_].highlight_;
  n = mtlcolor[s.m_].n_;
  N = s.normal(bay).normalize();
  V = viewdir.scalar(-1).normalize();

  ambient = od.scalar(ka);

  // Loop through lights
  for (int i = 0; i < lights.size(); i++) {
    if (lights[i]->type_ == 1) {
      Lbefore = lights[i]->components_ - intersection;
    } else {
      Lbefore = lights[i]->components_.scalar(-1);
    }

    float distance;
    float B, C, discrim;
    float t1, t2, ct;

    distance = Lbefore.length();
    L = Lbefore.normalize();
    fshadow = 200.f;

    // Avoid shadow bleeding
    bool in_shadow = false;

    // Shadow calculations
    for (int j = 0; j < objects.size(); j++) {
      VectorType otherpos = VectorType(objects[j]->x, objects[j]->y, objects[j]->z);
      VectorType otherdir = (intersection - otherpos).normalize();
      float length = (intersection - otherpos - otherdir.scalar(objects[j]->r)).length();

      // Checks if the length between surfaces larger than light distance
      // Also checks if it's a directional light or if it's in a shadow
      // Or if it's in the same direction as the light
       //if (in_shadow) break;
      if (length >= distance && lights[i]->type_ != 0 ||
          otherdir.dot(L.scalar(-1)) <= 0 || in_shadow) {
        continue;
      }

      // Cast 50 rays for softer shadows
      for (int k = 0; k < 50; k++) {
        float random[3];

        // Max 20% deviation from original direction
        for (int r = 0; r < 3; r++) {
          random[r] = (rand() % 200) / 1000.f;
        }   // r

        VectorType rand_vec = VectorType(random[0], random[1], random[2]).normalize();
        if (rand() % 1) rand_vec = rand_vec.scalar(-1.f);


        // Calculates B, C, and the discriminant for the ray-sphere intersection equation
        B = 2 * (L.dot(intersection - otherpos));
        C = (intersection - otherpos).dot(intersection - otherpos) - pow(objects[j]->r, 2);
        discrim = pow(B, 2) - 4 * C;

        // If there are multiple solutions, calculate the values for t
        if (discrim >= 0) {
          t1 = (-B + sqrt(discrim)) / 2;
          t2 = (-B - sqrt(discrim)) / 2;
          ct = std::max(t1, t2);

          if (t1 >= epsilon && t2 >= epsilon) {
            ct = std::min(t1, t2);
          }

          // For when intersection is behind
          if (ct <= epsilon) {
            continue;           // Go to next ray
          } else if (ct > epsilon) {
            // If already calculated, then don't add more than the shadow for one object
            in_shadow = true;

            // Subtracting from the shadow factor for darker shadows, alterable
            // Controls some of shadow fall-off, scalable by distance
            fshadow -= otherdir.dot(N.scalar(-4)) * ((lights.size()) / (length + 1));

            // Change shadow behavior if a directional lightcc
            // Side note: if only i had read the canvas page I would've saved
            // An entire day wondering why my shadows looked weird, and that
            // Was because my directional light fell off. Still don't know
            // Why it took so long to figure out tho...         
            if (lights[i]->type_ == 0) {
              fshadow = 0.f;
              break;
            }
          } else {
            std::cerr << "Problem has occurred in ShadeRay!" << std::endl;
            exit(1);
          }
        }
      }     // k
    }       // j




    for (int j = 0; j < tri.size(); j++) {
      // Checks if the length between surfaces larger than light distance
      // Also checks if it's a directional light or if it's in a shadow
      // Or if it's in the same direction as the light
       //if (in_shadow) break;
      if (in_shadow) {
        continue;
      }

      // Cast 50 rays for softer shadows
      for (int k = 0; k < 50; k++) {
        float random[3];

        // Max 20% deviation from original direction
        for (int r = 0; r < 3; r++) {
          random[r] = (rand() % 200) / 1000.f;
        }   // r

        VectorType rand_vec = VectorType(random[0], random[1], random[2]).normalize();
        if (rand() % 1) rand_vec = rand_vec.scalar(-1.f);

        float alpha = 1.f, beta = 1.f, gamma = 1.f;
        VectorType point = vb[tri[j]->vertices_[0] - 1];
        VectorType normal = tri[j]->normal();
        float tD = -normal.dot(point);
        float Bt = normal.dot(L);

        if (fabs(Bt) < 0.0001) {
          return bkgcolor;
        }

        float t = -(normal.dot(intersection) + tD) / Bt;
        VectorType intersect = intersection + L.scalar(t);
        bool res = tri[j]->check_point(intersect, alpha, beta, gamma);
        if (t > 0 && res) {
          in_shadow = true;
          fshadow -= L.dot(N.scalar(-4)) * ((1.f / lights.size()) / (Lbefore.length() + 1));
        } else if (t <= 0 || !res) {
          continue;
        } else {
          std::cerr << "Problem has occurred in ShadeRay!" << std::endl;
          exit(1);
        }
      }
    }     // k






    // Set the shadow factor and the attentuation factors
    // The attenuation used only has linear scaling with distance
    // Increase numerator to increase contrast/decreases fall-off rate
    // Basically a = 0, b = 0.125, c = 0 | a = 1
    attenuation = (lights[i]->type_ == 1) ? 8.f / Lbefore.length() : 2.f;
    fshadow = std::max(std::min(1.f, fshadow / 200.f), 0.f);

    // Calculate H, which is the normalized sum of the light direction and the
    // Negative normalzed viewing direction
    H = (L + V).normalize();

    // A running total of the diffusion and specular components
    // Multiplied by the attenuation factor and the shadow factor
    diffspec = diffspec + lights[i]->color_
               .scalar(attenuation * fshadow * (1)) *
               (od.scalar(kd * std::max(0.f, N.dot(L))) +
               os.scalar(ks * pow(std::max(0.f, N.dot(H)), n)));
    fshadow = 200.f;
  }           // i

  // Reset shadow factor

  // Clamp the result of the addition between the ambient light and the other
  // Light components using a vector component clamp

  float dc_dist = (intersection - eye).length();
  float dca = 0.f;

  if (dc_dist >= dc_far) {
    dca = dcmin;
  } else if (dc_dist <= dc_near) {
    dca = dcmax;
  } else {
    dca = (dcmax - dcmin) * (dc_far - dc_dist) / (dc_far - dc_near);
  }

  result = (ambient + diffspec).clamp(0.f, 1.f);
  result = result.scalar(dca) + depthcue.scalar(1.f - dca);
  return result.clamp(0.f, 1.f);
}

ColorType Trace_Ray(RayType ray) {
  float B, C, discrim;
  float t1, t2, ct;
  float min_t = 99999;
  int sindex = -1;
  bool triangle = false;

  VectorType ray_dir = VectorType(ray.dx, ray.dy, ray.dz);

  for (int s = 0; s < objects.size(); s++) {
    VectorType object_pos = VectorType(objects[s]->x, objects[s]->y, objects[s]->z);
    B = 2 * (ray_dir.dot(eye - object_pos));
    C = (eye - object_pos).dot(eye - object_pos) - pow(objects[s]->r, 2);
    discrim = pow(B, 2) - 4 * C;

    // Check if solutions exist and take the nearest one
    if (discrim >= 0) {
      t1 = (-B + sqrt(discrim)) / 2;
      t2 = (-B - sqrt(discrim)) / 2;
      ct = std::min(t1, t2);

      // Check if out of view/behind
      if (ct <= 0) {
        continue;   // Go to check next object
      }

      // Update new closest value and index
      if (ct <= min_t) {
        min_t = ct;
        sindex = s;
        triangle = false;
      }
    }
  }


  // Check triangles
  float alpha = 1.f, beta = 1.f, gamma = 1.f;
  for (int s = 0; s < tri.size(); s++) {
    VectorType point = vb[tri[s]->vertices_[0] - 1];
    VectorType normal = tri[s]->normal();
    float tD = -normal.dot(point);
    float Bt = normal.dot(ray_dir);

    if (fabs(Bt) < 0.0001) {
      return bkgcolor;
    }

    float t = -(normal.dot(eye) + tD) / Bt;
    float temp_min_t = std::min(t, min_t);
    VectorType intersect = eye + ray_dir.scalar(temp_min_t);

    if (t < min_t && tri[s]->check_point(intersect, alpha, beta, gamma)) {
      triangle = true;
      sindex = s;
      min_t = temp_min_t;
    }
  }



  // Once checked all objects, check if index was not updated
  if (sindex < 0) {
    return bkgcolor;
  } 

  // Calculate intersection point and pass it and the sphere to ShadeRay()
  VectorType intersect = eye + ray_dir.scalar(min_t);
  if (triangle) {
    return Shade_RayT(intersect.x, intersect.y, intersect.z, *tri[sindex],
                      VectorType(alpha, beta, gamma));
  } return Shade_Ray(intersect.x, intersect.y, intersect.z, *objects[sindex]);
}

#endif  // RAY_H_


