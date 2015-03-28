#ifndef LIST_H
#define LIST_H

struct Node;
typedef struct {
	void *payload;
	Node *next;
} Node;

// singly linked circular list
typedef struct {
	int size;  // size (in bytes) of type stored in list
	Node *sentinel;  // dummy node which marks the beginning and end of the list
} List;

List list_init(int size);
void list_add(List list, void *payload);

#endif //LIST_H
