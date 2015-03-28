#include "list.h"

List* list_init(int size) {
	List *ans = (List *)malloc(sizeof(List));

	ans->size = size;

	Node *sentinel = (Node *)malloc(sizeof(Node));
	sentinel->next = sentinel;

	ans->sentinel = sentinel;

	return ans;
}
