
#ifndef __LIST_H
#define __LIST_H


template <class TObject, class TKey = int>
class List;


template <class TObject, class TKey = int>
class ListItem {
	TObject *m_object;
	TKey m_key;
	ListItem *m_next;
	ListItem *m_prev;

public:
	ListItem();
	ListItem(TObject *obj, TKey key = 0);
	ListItem(const ListItem &item);
	~ListItem();

	TObject *object() const { return m_object; }
	TKey key() const { return m_key; }
	ListItem *next() const { return m_next; }
	ListItem *prev() const { return m_prev; }

	void setObject(TObject *obj) { m_object = obj; }
	void setKey(TKey key) { m_key = key; }

	void removeObject();
	void releaseObject();

	friend class List<TObject, TKey>;
};


template <class TObject, class TKey>
class List {
protected:
	ListItem<TObject, TKey> *m_first;
	ListItem<TObject, TKey> *m_last;
	ListItem<TObject, TKey> *m_current;
	int m_size;
	bool m_is_copy;

public:
	List();
	List(const List &list);
	~List();

	ListItem<TObject, TKey> *first() const { return m_first; }
	ListItem<TObject, TKey> *last() const { return m_last; }
	ListItem<TObject, TKey> *current() const { return m_current; }
	int size() const { return m_size; }
	bool is_copy() const { return m_is_copy; }

	void setIsCopy(bool is_copy) { m_is_copy = is_copy; }

	List<TObject, TKey> &operator =(const List &list) { return assign(list); }
	List<TObject, TKey> &assign(const List &list);

	ListItem<TObject, TKey> *moveFirst();
	ListItem<TObject, TKey> *moveLast();
	ListItem<TObject, TKey> *moveNext();
	ListItem<TObject, TKey> *movePrev();
	ListItem<TObject, TKey> *findItem(TKey key);

	TObject *getFirst();
	TObject *getLast();
	TObject *getNext();
	TObject *getPrev();
	TObject *getCurrent() const;
	TObject *getItem(ListItem<TObject, TKey> *item);
	TObject *get(TKey key);

	int addFirstItem(ListItem<TObject, TKey> *item);
	int addLastItem(ListItem<TObject, TKey> *item);
	int addNextItem(ListItem<TObject, TKey> *item);
	int addPrevItem(ListItem<TObject, TKey> *item);
	int addAscentItem(ListItem<TObject, TKey> *item);
	int addDescentItem(ListItem<TObject, TKey> *item);

	int addFirst(TObject *obj, TKey key = 0);
	int addLast(TObject *obj, TKey key = 0);
	int addNext(TObject *obj, TKey key = 0);
	int addPrev(TObject *obj, TKey key = 0);
	int addAscent(TObject *obj, TKey key = 0);
	int addDescent(TObject *obj, TKey key = 0);

	int removeFirst();
	int removeLast();
	int removeCurrent();
	int removeAll();
	int removeItem(ListItem<TObject, TKey> *item);
	int remove(TKey key);

	int releaseFirst();
	int releaseLast();
	int releaseCurrent();
	int releaseAll();
	int releaseItem(ListItem<TObject, TKey> *item);
	int release(TKey key);

	int toArray(TObject **&array);
	int fromArray(TObject **&array, int size);

protected:
	void deleteItem(ListItem<TObject, TKey> *item);
};


////////////////////////////////
//
// implementation of ListItem


template <class TObject, class TKey>
inline ListItem<TObject, TKey>::ListItem()
{
	m_object = NULL;
	m_key = 0;
	m_next = m_prev = NULL;
}

template <class TObject, class TKey>
inline ListItem<TObject, TKey>::ListItem(TObject *obj, TKey key)
{
	m_object = obj;
	m_key = key;
	m_next = m_prev = NULL;
}

template <class TObject, class TKey>
inline ListItem<TObject, TKey>::ListItem(const ListItem &item)
{
	if (item.m_object != NULL) {
		m_object = new TObject(*item.m_object);
	}
	else {
		m_object = NULL;
	}
	m_key = item.m_key;
	m_next = m_prev = NULL;
}

template <class TObject, class TKey>
inline ListItem<TObject, TKey>::~ListItem()
{
	removeObject();
}

template <class TObject, class TKey>
inline void ListItem<TObject, TKey>::removeObject()
{
	delete m_object;
	m_object = NULL;
}

template <class TObject, class TKey>
inline void ListItem<TObject, TKey>::releaseObject()
{
	m_object = NULL;
}


////////////////////////////////
//
// implementation of List


template <class TObject, class TKey>
inline List<TObject, TKey>::List()
{
	m_first = m_last = m_current = NULL;
	m_size = 0;
	m_is_copy = false;
}

template <class TObject, class TKey>
inline List<TObject, TKey>::List(const List &list)
{
	ListItem<TObject, TKey> *temp = list.m_first;
	while (temp != NULL) {
		addLastItem(new ListItem<TObject, TKey>(*temp));
		temp = temp->m_next;
	}
}

template <class TObject, class TKey>
inline List<TObject, TKey>::~List()
{
	removeAll();
}

template <class TObject, class TKey>
inline List<TObject, TKey> &List<TObject, TKey>::assign(const List &list)
{
	ListItem<TObject, TKey> *temp;
	removeAll();
	temp = list.m_first;
	while (temp != NULL) {
		addLastItem(new ListItem<TObject, TKey>(*temp));
		temp = temp->m_next;
	}
	return *this;
}

template <class TObject, class TKey>
inline ListItem<TObject, TKey> *List<TObject, TKey>::moveFirst()
{
	m_current = m_first;
	return m_current;
}

template <class TObject, class TKey>
inline ListItem<TObject, TKey> *List<TObject, TKey>::moveLast()
{
	m_current = m_last;
	return m_current;
}

template <class TObject, class TKey>
inline ListItem<TObject, TKey> *List<TObject, TKey>::moveNext()
{
	if (m_current != NULL) {
		m_current = m_current->m_next;
	}
	else {
		m_current = m_first;
	}
	return m_current;
}

template <class TObject, class TKey>
inline ListItem<TObject, TKey> *List<TObject, TKey>::movePrev()
{
	if (m_current != NULL) {
		m_current = m_current->m_prev;
	}
	else {
		m_current = m_last;
	}
	return m_current;
}

template <class TObject, class TKey>
inline ListItem<TObject, TKey> *List<TObject, TKey>::findItem(TKey key)
{
	ListItem<TObject, TKey> *temp = m_first;
	while (temp != NULL && temp->m_key != key) {
		temp = temp->m_next;
	}
	return temp;
}

template <class TObject, class TKey>
inline TObject *List<TObject, TKey>::getFirst()
{
	moveFirst();
	return getCurrent();
}

template <class TObject, class TKey>
inline TObject *List<TObject, TKey>::getLast()
{
	moveLast();
	return getCurrent();
}

template <class TObject, class TKey>
inline TObject *List<TObject, TKey>::getNext()
{
	moveNext();
	return getCurrent();
}

template <class TObject, class TKey>
inline TObject *List<TObject, TKey>::getPrev()
{
	movePrev();
	return getCurrent();
}

template <class TObject, class TKey>
inline TObject *List<TObject, TKey>::getCurrent() const
{
	if (m_current != NULL) {
		return m_current->m_object;
	}
	else {
		return NULL;
	}
}

template <class TObject, class TKey>
inline TObject *List<TObject, TKey>::getItem(ListItem<TObject, TKey> *item)
{
	m_current = item;
	return getCurrent();
}

template <class TObject, class TKey>
inline TObject *List<TObject, TKey>::get(TKey key)
{
	ListItem<TObject, TKey> *item = findItem(key);
	if (item != NULL) {
		m_current = item;
		return m_current->m_object;
	}
	else {
		return NULL;
	}
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::addFirstItem(ListItem<TObject, TKey> *item)
{
	if (m_first != NULL) {
		item->m_next = m_first;
		item->m_prev = NULL;
		m_first->m_prev = item;
		m_first = item;
		m_size++;
	}
	else {
		m_first = m_last = item;
		item->m_next = item->m_prev = NULL;
		m_size = 1;
	}
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::addLastItem(ListItem<TObject, TKey> *item)
{
	if (m_last != NULL) {
		item->m_prev = m_last;
		item->m_next = NULL;
		m_last->m_next = item;
		m_last = item;
		m_size++;
	}
	else {
		m_first = m_last = item;
		item->m_prev = item->m_next = NULL;
		m_size = 1;
	}
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::addNextItem(ListItem<TObject, TKey> *item)
{
	if (m_current == m_last) {
		return addLastItem(item);
	}
	else if (m_current != NULL) {
		item->m_prev = m_current;
		item->m_next = m_current->m_next;
		m_current->m_next->m_prev = item;
		m_current->m_next = item;
		m_size++;
		return m_size;
	}
	else {
		return addFirstItem(item);
	}
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::addPrevItem(ListItem<TObject, TKey> *item)
{
	if (m_current == m_first) {
		return addFirstItem(item);
	}
	else if (m_current != NULL) {
		item->m_next = m_current;
		item->m_prev = m_current->m_prev;
		m_current->m_prev->m_next = item;
		m_current->m_prev = item;
		m_size++;
		return m_size;
	}
	else {
		return addLastItem(item);
	}
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::addAscentItem(ListItem<TObject, TKey> *item)
{
	ListItem<TObject, TKey> *temp = m_current;
	m_current = m_first;
	while (m_current != NULL && m_current->m_key < item->m_key) {
		m_current = m_current->m_next;
	}
	addPrevItem(item);
	m_current = temp;
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::addDescentItem(ListItem<TObject, TKey> *item)
{
	ListItem<TObject, TKey> *temp = m_current;
	m_current = m_last;
	while (m_current != NULL && m_current->m_key < item->m_key) {
		m_current = m_current->m_prev;
	}
	addNextItem(item);
	m_current = temp;
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::addFirst(TObject *obj, TKey key)
{
	return addFirstItem(new ListItem<TObject, TKey> (obj, key));
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::addLast(TObject *obj, TKey key)
{
	return addLastItem(new ListItem<TObject, TKey> (obj, key));
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::addNext(TObject *obj, TKey key)
{
	return addNextItem(new ListItem<TObject, TKey> (obj, key));
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::addPrev(TObject *obj, TKey key)
{
	return addPrevItem(new ListItem<TObject, TKey> (obj, key));
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::addAscent(TObject *obj, TKey key)
{
	return addAscentItem(new ListItem<TObject, TKey> (obj, key));
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::addDescent(TObject *obj, TKey key)
{
	return addDescentItem(new ListItem<TObject, TKey> (obj, key));
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::removeFirst()
{
	ListItem<TObject, TKey> *temp;
	if (m_first == m_last) {
		deleteItem(m_first);
		m_first = m_last = m_current = NULL;
		m_size = 0;
	}
	else {
		if (m_first == m_current) {
			m_current = m_first->m_next;
		}
		m_first->m_next->m_prev = NULL;
		temp = m_first;
		m_first = m_first->m_next;
		deleteItem(temp);
		m_size--;
	}
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::removeLast()
{
	ListItem<TObject, TKey> *temp;
	if (m_last == m_first) {
		deleteItem(m_last);
		m_first = m_last = m_current = NULL;
		m_size = 0;
	}
	else {
		if (m_last == m_current) {
			m_current = m_last->m_next;
		}
		m_last->m_prev->m_next = NULL;
		temp = m_last;
		m_last = m_last->m_prev;
		deleteItem(temp);
		m_size--;
	}
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::removeCurrent()
{
	ListItem<TObject, TKey> *temp;
	if (m_current == m_first) {
		removeFirst();
	}
	else if (m_current == m_last) {
		removeLast();
	}
	else if (m_current != NULL) {
		m_current->m_prev->m_next = m_current->m_next;
		m_current->m_next->m_prev = m_current->m_prev;
		temp = m_current;
		m_current = m_current->m_next;
		deleteItem(temp);
		m_size--;
	}
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::removeAll()
{
	ListItem<TObject, TKey> *temp;
	m_current = m_first;
	while (m_current != NULL) {
		temp = m_current;
		m_current = m_current->m_next;
		deleteItem(temp);
	}
	m_first = m_last = m_current = NULL;
	m_size = 0;
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::removeItem(ListItem<TObject, TKey> *item)
{
	ListItem<TObject, TKey> *temp;
	if (item != NULL) {
		if (item == m_current) {
			removeCurrent();
		}
		else {
			temp = m_current;
			m_current = item;
			removeCurrent();
			m_current = temp;
		}
	}
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::remove(TKey key)
{
	return removeItem(findItem(key));
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::releaseFirst()
{
	if (m_first != NULL) {
		m_first->releaseObject();
		removeFirst();
	}
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::releaseLast()
{
	if (m_last != NULL) {
		m_last->releaseObject();
		removeLast();
	}
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::releaseCurrent()
{
	if (m_current != NULL) {
		m_current->releaseObject();
		removeCurrent();
	}
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::releaseAll()
{
	ListItem<TObject, TKey> *temp;
	temp = m_first;
	while (temp != NULL) {
		temp->releaseObject();
		temp = temp->m_next;
	}
	return removeAll();
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::releaseItem(ListItem<TObject, TKey> *item)
{
	if (item != NULL) {
		item->releaseObject();
		removeItem(item);
	}
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::release(TKey key)
{
	return releaseItem(findItem(key));
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::toArray(TObject **&array)
{
	int i;
	ListItem<TObject, TKey> *temp = m_current;
	array = new TObject * [m_size];
	moveFirst();
	for (i=0; i<m_size; i++) {
		array[i] = getCurrent();
		moveNext();
	}
	m_current = temp;
	return m_size;
}

template <class TObject, class TKey>
inline int List<TObject, TKey>::fromArray(TObject **&array, int size)
{
	int i;
	removeAll();
	for (i=0; i<size; i++) {
		addLast(array[i]);
	}
	return m_size;
}

template <class TObject, class TKey>
inline void List<TObject, TKey>::deleteItem(ListItem<TObject, TKey> *item)
{
	if (m_is_copy) item->releaseObject();
	delete item;
}


#endif // __LIST_H
