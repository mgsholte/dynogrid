#include "list.h"

List list_init(int size) {
	Node *sentinel = (Node *)malloc(sizeof(Node));
	sentinel->next = sentinel;

	return (List) { size, sentinel };
}

void list_add(List list, void *payload) {
	Node *new_node = (Node *)malloc(sizeof(Node));

	new_node->payload = payload;
	new_node->next = (list.sentinel)->next;
	(list.sentinel)->next = new_node;
}
