#include <stdlib.h>
#include <string.h> // for memcpy
#include <stdio.h> // for printf in debugging

#include "list.h"

List* list_init() {
	Node *sentinel = (Node*) malloc(sizeof(Node));
	sentinel->payload = NULL;
	sentinel->next = sentinel;

	List *list = (List*) malloc(sizeof(List));
	// the following allows us to set sentinel even though it is a const pointer
	List tmp = (List) { sentinel, sentinel, sentinel, 0 };
	memcpy(list, &tmp, sizeof(List));

	return list;
}

void list_free(List *list, bool free_payloads) {
	list_reset_iter(list);
	while(list_has_next(list)) {
		list_get_next(list);
		list_pop(list, free_payloads);
	}
	//TODO: is this necessary?
	list_reset_iter(list);
	if (list_has_next(list)) {
		list_get_next(list);
		list_pop(list, free_payloads);
	}
	free(list->sentinel);
	//TODO: should be unnecessary
	memset(list, (size_t) NULL, sizeof(Node*)); // <==> list->sentinel = NULL;
	

	free(list);
}

void list_add(List *list, void *payload) {
	Node *new_node = (Node*) malloc(sizeof(Node));

	new_node->payload = payload;
	new_node->next = list->sentinel->next;
	list->sentinel->next = new_node;
	++(list->length);
}

// remove the node currently being iterated (the one last returned by list_get_next)
void list_pop(List *l, bool free_payload) {
	Node *x = l->prev->next;  // get node to be removed
	if (x == l->sentinel) {
		printf("warning! popping the sentinel!\n");
	}
	l->prev->next = x->next;  // unlink it
	l->iter = l->prev;        // reset iterator
	if (free_payload) {
		free(x->payload);         // free the payload at the deleted node
	}
	free(x);                  // free node itself
	--(l->length);
}

// Moves an elem (the one last returned by list_get_next(donor)) from a donor list to a recipient list. donated item will be 1st elem in recip list
// no memory is alloced or freed in the process
void list_pass(List *recip, List *donor) {
	// get the node to be passed
	Node *x = donor->iter; 
	if (x == donor->sentinel) {
		printf("warning! passing the sentinel!\n");
	}
	// unlink it from donor list
	donor->prev->next = x->next;
	// reset donor iter and elem count
	donor->iter = donor->prev;
	--(donor->length);
	// link node into recip list
	x->next = recip->sentinel->next;
	recip->sentinel->next = x;
	++(recip->length);
}

// Add everithing in the donor list to the recipient list
// You have to traverse the list to connect the last node no matter what, so this probably isn't too slow
void list_combine(List *recip, List *donor) {
	list_reset_iter(donor);
	while (list_has_next(donor)) {
		list_get_next(donor);
		list_pass(recip, donor);
	}
	//TODO: debugging
	if (list_length(donor) != 0) {
		printf("warning, found a bug in list_combine\n");
	}
}

void list_reset_iter(List *l) {
	l->iter = l->sentinel;
	l->prev = l->sentinel;
}

bool list_has_next(List *l) {
	return (bool) (l->iter->next != l->sentinel);
}

//NB: iter should point to the node whose payload it returns. it must start off as the sentinel after a reset
void* list_get_next(List *l) {
	l->prev = l->iter;
	l->iter = l->iter->next;

	return l->iter->payload;
}

static int node_length(Node *sentinel, Node *cur, int acc) {
	return (cur == sentinel)
		? acc
		: node_length(sentinel, cur->next, acc+1);
}

int list_length(List *list) {
	return list->length;
	//return node_length(list.sentinel, list.sentinel->next, 0);
}
