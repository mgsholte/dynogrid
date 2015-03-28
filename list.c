#include "list.h"

List list_init(int size) {
	Node *sentinel = (Node *)malloc(sizeof(Node));
	sentinel->next = sentinel;

	return (List) { size, sentinel };
}
