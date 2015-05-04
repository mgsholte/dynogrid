#ifndef VECTOR_H
#define VECTOR_H

// a 2D vector
typedef struct {
	double x,y;
} vec2;

// a 3D vector
typedef struct {
	double x,y,z;
} vec3;

// scale the vector v (in-place) by the scalar factor
void vec3_scale(vec3 *v, double factor);

// add the vectors v,w and store result in v
void vec3_add(vec3 *v, vec3 w);

#endif //VECTOR_H
