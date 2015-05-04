#include "vector.h"

void vec3_scale(vec3 *v, double factor) {
	v->x *= factor;	v->y *= factor;	v->z *= factor;
}

void vec3_add(vec3 *v, vec3 w) {
	v->x += w.x; v->y += w.y; v->z += w.z;
}
