#include "list.h"

List list_init() {
	Node *sentinel = (Node*) malloc(sizeof(Node));
	sentinel->next = sentinel;

	return (List) { sentinel, sentinel };
}

void list_add(List list, void *payload) {
	Node *new_node = (Node*) malloc(sizeof(Node));

	new_node->payload = payload;
	new_node->next = (list.sentinel)->next;
	(list.sentinel)->next = new_node;
}

void list_reset_iter(List l) {
	l.iter = sentinel->next;
}

bool list_has_next(List l) {
	return l.iter != l.sentinel;
}

void* list_get_next(List l) {
	void *ans = l.iter->payload;
	l.iter = l.iter->next;

	return ans;
}

