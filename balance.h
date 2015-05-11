#ifndef BALANCE_H
#define BALANCE_H

#include "decs.h"
#include "list.h"

typedef struct {
	List cells;
	int neighbor;
	int dir; //dir for "direction"...always either negatitve or postive for left or right, up or down, in or out
	int layer;
} surface;



#endif //BALANCE_H
