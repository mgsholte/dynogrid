#ifndef LIST_H
#define LIST_H

#include "decs.h"

typedef struct Node {
	void *payload;
	struct Node *next;
} Node;

// singly linked circular list
typedef struct {
	Node *sentinel;  // dummy node which marks the beginning and end of the list
	Node *iter;      // pointer to next node returned when iterating over the list
	Node *prev;      // allows for popping the current element during an iteration
	int length;
} List;

// list creation/deletion
List list_init();
void list_free(List list);

// list modification
void list_add(List *list, void *payload);
void list_pop(List *list);

// list iteration
void list_reset_iter(List *l);
bool list_has_next(List l);
void* list_get_next(List *l);

// list query
int list_length(List list);

// particle passing
void list_pass(List *recip, List *donor, void *payload);
void list_combine(List *recip, List *donor);

#endif //LIST_H
