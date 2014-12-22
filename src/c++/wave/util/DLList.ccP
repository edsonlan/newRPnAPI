#ifdef __GNUC__
#pragma implementation
#endif
#include <stdlib.h>
#include <<T>DLList.h>

void
<T>DLList::append(const <T>DLList & list)
{
	<T>DLListNode::reference(list.head);
	if (head == NULL)
		head = list.head;
	else {
		tail->next = list.head;
		if (list.head != NULL)
			list.head->prev = tail;
	}
	tail = list.tail;
}

void
<T>DLList::prepend(const <T>DLList & list)
{
	<T>DLListNode *h = list.head;
	if (h != NULL) {
		<T>DLListNode *node = new <T>DLListNode(h->obj);
		node->prev = NULL;
		<T>DLListNode *trailer = node;
		for (h = h->next; h != NULL; h = h->next) {
			<T>DLListNode *n = new <T>DLListNode(h->obj);
			trailer->next = n;
			n->prev = trailer;
			trailer = n;
		}
		trailer->next = head;
		if (head != NULL)
			head->prev = trailer;
		head = node;
	}
}

void
<T>DLList::del(const <T> & obj)
{
	<T>DLListNode *h = head;

	for ( ; ; ) {
		if (h == NULL) {
			tail = head = h;
			return;
		}
		else if (&h->obj == &obj) {
			<T>DLListNode *n = h->next;
			if (n != NULL)
				n->prev = NULL;
			if (--h->ref == 0)
				delete h;
			h = n;
		}
		else
			break;
	}

	<T>DLListNode *prev = h;
	<T>DLListNode *node = h->next;
	while (node != NULL) {
		if (&node->obj == &obj) {
			<T>DLListNode *n = node->next;
			if (n != NULL)
				n->prev = node->prev;
			if (--node->ref == 0)
				delete node;
			prev->next = n;
			node = n;
		}
		else {
			prev = node;
			node = node->next;
		}
	}
	head = h;
	tail = prev;
}

void
<T>DLList::del_last(const <T> & obj)
{
	<T>DLListNode *t = tail;

	if (t == NULL)
		return;

	if (&t->obj == &obj) {
		<T>DLListNode *p = t->prev;
		if (p != NULL)
			p->next = NULL;
		else
			head = NULL;
		tail = p;
		if (--t->ref == 0)
			delete t;
		return;
	}

	<T>DLListNode *next = t;
	<T>DLListNode *node = t->prev;
	while (node != NULL) {
		if (&node->obj == &obj) {
			<T>DLListNode *p = node->prev;
			if (p != NULL)
				p->next = node->next;
			else
				head = next;
			next->prev = p;
			if (--node->ref == 0)
				delete node;
			return;
		}
		else {
			next = node;
			node = node->prev;
		}
	}
}
