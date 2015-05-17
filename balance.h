#ifndef BALANCE_H
#define BALANCE_H

#include "decs.h"
#include "list.h"
#include "tree.h"

typedef struct {
	List *cells;
	int neighbor;
	int dir; //dir for "direction"...always either negatitve or postive for left or right, up or down, in or out
	int layer;
} surface;

void Balance(tree ****grid);
void give_take_surface(tree**** base_grid, List* list_u, int id_u, List* list_d, int id_d);
void determine_neighbor_matchings(List* ne_matchings[], char dim, tree**** grid);
void resize_allocation(tree**** base_grid);
void convert_ghost2real_and_reghost(tree**** base_grid, tree* new_tree, char dir6);

#endif //BALANCE_H
