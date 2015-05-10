#include "vector.h"
#include "math.h"

double vec3_norm(const vec3 v) {
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

void vec3_scale(vec3 *v, double factor) {
	v->x *= factor;	v->y *= factor;	v->z *= factor;
}

void vec3_add(vec3 *v, vec3 w) {
	v->x += w.x; v->y += w.y; v->z += w.z;
}
