## Assignment 1C
# Introduction
In this portion of the assignment we extend the functionality of our raycaster to be able to
draw triangles. Additionally, smooth shading and textures are implemented.

# How to Use
There is a makefile inside the folder. You should run: `make` in order to create an executable
called `gen_image`. Then you can do `./gen_image <location_of_input_file>` in order to run the
program. The .ppm image will be generated in the folder that the input is in.

# Implemented Rubric Criteria
## Note
This is what I think I implemented.

## Criteria
* ___The program robustly accepts extended scene description files that include texture images, texture coordinates and surface normal vectors.  The program is able to robustly handle triangle definitions that include per-vertex normal directions and/or per-vertex texture coordinates in addition to vertex locations. The implementation is done in a way that enables easily working with triangle mesh models originally defined in .obj format. (5pts)\
    * Able to accept extended scene description, yes. Robustly, no.\
    * Easy to work with, yes. Code structure, yikes.\
* ___The program correctly computes ray/plane intersections, and correctly performs point-in-triangle testing using Barycentric coordinates, enabling the rendering of scenes containing triangles as well as spheres. (20pts)\
    * Yes.\
* ___The program is capable of rendering triangles using flat shading, in which every pixel in a triangle is assigned the same color, obtained by correctly evaluating the Phong illumination equation using the unit length normal of the plane in which the triangle lies. (10pts)\
    * Yes, although the inputs provided for the sphere_square.txt confused me. It specified smooth shading, but the provided image looked like it was flat shaded. Changing the input file to actually use flat shading made it look like the provided image.\
* ___The program is capable of rendering triangles using smooth shading, in which every pixel within a triangle is assigned a unique color, obtained by correctly evaluating the Phong illumination equation using a unit length normal direction interpolated from the three normal directions defined at the three triangle vertices. (15pts)\
    * Yes.\
* ___The program is capable of rendering textured spheres. An appropriate texture coordinate is computed at each ray/sphere intersection point using a pre-defined, hard-coded mapping. That texture coordinate is used to retrieve a correctly corresponding color from the texture map, which specifies the object’s diffuse color in the Phong illumination model at that point. (20pts)\
    * Yes.\
* ___The program is capable of rendering textured triangles.  The texture coordinate at the ray/triangle intersection point is correctly interpolated from the texture coordinates defined at each of the three triangle vertices. The interpolated texture coordinate is used to retrieve a correctly corresponding color from the texture map, which specifies the object’s diffuse color in the Phong illumination model at that point. (25pts)\
    * Yes.\
* ___The student has submitted one or more scene description files and accompanying rendered images to successfully demonstrate all of the capabilities of their program, including the ability to render: at least one flat-shaded triangle, multiple smooth-shaded triangles, one or more textured triangles, and one or more textured spheres. (5pts)\
    * Yes. Flat shaded triangles include the black pin on the map and the right two triangles. The smooth shaded triangles are the ones in the middle. The textured triangles are the folded map. The texture sphere is the globe. Untextured sphere is also included as the ball on the pin.\
* ___Extra credit: The program is capable of reading a normal map from a file and correctly using the values in that map to vary the surface normal direction used in calculating the illumination equation at each point across a bump-mapped surface. (7pts)\
    * Nope.
