#include <stdlib.h>

#include "list.h"

List list_init() {
	Node *sentinel = (Node*) malloc(sizeof(Node));
	sentinel->next = sentinel;

	return (List) { sentinel, sentinel, sentinel };
}

void list_free(List list) {
	list_reset_iter(&list);
	while(list_has_next(list)) {
		list_pop(&list);
	}
	free(list.iter);
	free(list.prev);
	free(list.sentinel);
}

void list_add(List *list, particle *payload) {
	Node *new_node = (Node*) malloc(sizeof(Node));

	new_node->payload = payload;
	new_node->next = (list->sentinel)->next;
	(list->sentinel)->next = new_node;
}

// remove the node currently being iterated
void list_pop(List *l) {
	Node *x = l->prev->next;  // get node to be removed
	l->prev->next = x->next;  // unlink it
	l->iter = x->next;        // reset iterator
	free(x->payload);         // free the particle at the deleted node
	free(x);                  // free node itself
}

void list_reset_iter(List *l) {
	l->iter = l->sentinel;
	l->prev = l->sentinel;
}

bool list_has_next(List l) {
	return l.iter->next != l.sentinel;
}

//NB: should point to the node whose payload it returns. it must start off as the sentinel after a reset
particle* list_get_next(List *l) {
	l->prev = l->iter;
	l->iter = l->iter->next;

	return l->iter->payload;
}

static int node_length(Node *sentinel, Node *cur, int acc) {
	return (cur == sentinel)
		? acc
		: node_length(sentinel, cur->next, acc+1);
}

int list_length(List list) {
	return node_length(list.sentinel, list.sentinel->next, 0);
}
