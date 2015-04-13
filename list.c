#include "list.h"

List list_init() {
	Node *sentinel = (Node*) malloc(sizeof(Node));
	sentinel->next = sentinel;

	return (List) { sentinel, sentinel };
}

void list_add(List *list, particle *payload) {
	Node *new_node = (Node*) malloc(sizeof(Node));

	new_node->payload = payload;
	new_node->next = (list->sentinel)->next;
	(list->sentinel)->next = new_node;
}

// remove the node that comes after prev
void list_pop(List *list, Node *prev) {
	Node *x = prev->next;
	list->iter = prev;
	prev->next = x->next;
	free(x);
}

void list_reset_iter(List *l) {
	l->iter = l->sentinel->next;
}

bool list_has_next(List l) {
	return l.iter != l.sentinel;
}

particle* list_get_next(List *l) {
	particle *ans = l->iter->payload;
	l->iter = l->iter->next;

	return ans;
}

static int node_length(Node *sentinel, Node *cur, int acc) {
	if (cur == sentinel)
		return acc;
	else
		return node_length(sentinel, cur->next, acc+1);
}

int list_length(List list) {
	return node_length(list.sentinel, list.sentinel->next, 0);
}
