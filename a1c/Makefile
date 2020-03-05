all: raycaster.cpp ray.h vectortype.h colortype.h spheretype.h materialtype.h lighttype.h
	g++ raycaster.cpp -o gen_image

.PHONY: clean

clean:
	rm -f *.ppm gen_image scenes/*.ppm
