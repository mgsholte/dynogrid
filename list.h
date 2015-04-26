#ifndef LIST_H
#define LIST_H

#include "decs.h"

typedef struct Node {
	particle *payload;
	struct Node *next;
} Node;

// singly linked circular list
typedef struct {
	Node *sentinel;  // dummy node which marks the beginning and end of the list
	Node *iter;      // pointer to next node returned when iterating over the list
	Node *prev;      // allows for popping the current element during an iteration
} List;

// list modification
List list_init();
void list_free(List list);
void list_add(List list, particle *payload);
void list_pop(List *list);

// list iteration
void list_reset_iter(List *l);
bool list_has_next(List l);
particle* list_get_next(List *l);

// list query
int list_length(List list);

#endif //LIST_H
