#include "list.h"

List list_init() {
	Node *sentinel = (Node*) malloc(sizeof(Node));
	sentinel->next = sentinel;

	return (List) { sentinel, sentinel };
}

void list_add(List list, particle *payload) {
	Node *new_node = (Node*) malloc(sizeof(Node));

	new_node->payload = payload;
	new_node->next = (list.sentinel)->next;
	(list.sentinel)->next = new_node;
}

void list_reset_iter(List l) {
	l.iter = l.sentinel->next;
}

bool list_has_next(List l) {
	return l.iter != l.sentinel;
}

particle* list_get_next(List l) {
	particle *ans = l.iter->payload;
	l.iter = l.iter->next;

	return ans;
}

