
#ifndef __TREE_H
#define __TREE_H


#include <stdlib.h>


template <class TObject, class TKey = int>
class Tree;


template <class TObject, class TKey = int>
class TreeNode {
	TObject *m_object;
	TKey m_key;
	TreeNode *m_parent;
	TreeNode **m_children;
	int m_children_num;

public:
	TreeNode(int child_num = 2);
	TreeNode(TObject *obj, TKey key = 0, int child_num = 2);
	TreeNode(const TreeNode &node);
	virtual ~TreeNode();

	TObject *object() const { return m_object; }
	int children_num() const { return m_children_num; }

	void setObject(TObject *obj, TKey key = 0) { m_object = obj; m_key = key; }
	void setKey(TKey key = 0) { m_key = key; }

	void setChildrenNum(int child_num);

	void addChildNode(int i, TreeNode *node) { m_children[i] = node; }
	TreeNode *getChildNode(int i) const { return m_children[i]; }
	void removeChildNode(int i) { delete m_children[i]; m_children[i] = NULL; }

	virtual void addChild(int i, TObject *obj = NULL, TKey key = 0);
	virtual TObject *getChild(int i) const;

	int getSize() const;
	void release();

	friend class Tree<TObject, TKey>;
};


template <class TObject, class TKey>
class Tree {
protected:
	TreeNode<TObject, TKey> m_root;
	int m_children_num;

public:
	Tree(int child_num = 2);
	Tree(const Tree &tree);
	~Tree();

	TreeNode<TObject, TKey> &root() { return m_root; }
	int children_num() const { return m_children_num; }

	void setChildrenNum(int child_num);

	int getSize() const { return (m_root.getSize()-1); }
	void release() { m_root.release(); }

	TreeNode<TObject, TKey> *newTreeNode(TObject *obj, TKey key = 0);
};


////////////////////////////////
//
// implementation of TreeNode

template <class TObject, class TKey>
inline TreeNode<TObject, TKey>::TreeNode(int child_num)
{
	int i;
	m_children_num = child_num;
	m_children = new TreeNode * [m_children_num];
	for (i=0; i<m_children_num; i++) {
		m_children[i] = NULL;
	}
	m_parent = NULL;
	m_object = NULL;
	m_key = 0;
}

template <class TObject, class TKey>
inline TreeNode<TObject, TKey>::TreeNode(TObject *obj, TKey key, int child_num)
{
	int i;
	m_children_num = child_num;
	m_children = new TreeNode * [m_children_num];
	for (i=0; i<m_children_num; i++) {
		m_children[i] = NULL;
	}
	m_parent = NULL;
	m_object = obj;
	m_key = key;
}

template <class TObject, class TKey>
inline TreeNode<TObject, TKey>::TreeNode(const TreeNode &node)
{
	int i;
	m_children_num = node.m_children_num;
	m_children = new TreeNode * [m_children_num];
	for (i=0; i<m_children_num; i++) {
		if (node.m_children[i] != NULL) {
			m_children[i] = new TreeNode (*node.m_children[i]);
			m_children[i]->m_parent = this;
		}
		else {
			m_children[i] = NULL;
		}
	}
	m_parent = node.m_parent;
	if (node.m_object != NULL) {
		m_object = new TObject (*node.m_object);
	}
	else {
		m_object = NULL;
	}
	m_key = node.m_key;
}

template <class TObject, class TKey>
inline TreeNode<TObject, TKey>::~TreeNode()
{
	delete m_object;
	m_key = 0;
	m_parent = NULL;
	delete[] m_children;
	m_children = NULL;
	m_children_num = 0;
}

template <class TObject, class TKey>
inline void TreeNode<TObject, TKey>::setChildrenNum(int child_num)
{
	int i;
	if (m_children_num != child_num) {
		delete[] m_children;
		m_children_num = child_num;
		m_children = new TreeNode * [m_children_num];
		for (i=0; i<m_children_num; i++) {
			m_children[i] = NULL;
		}
	}
}

template <class TObject, class TKey>
inline void TreeNode<TObject, TKey>::addChild(int i, TObject *obj, TKey key)
{
	m_children[i] = new TreeNode (obj, key, m_children_num);
}

template <class TObject, class TKey>
inline TObject *TreeNode<TObject, TKey>::getChild(int i) const
{
	if (m_children[i] != NULL) {
		return m_children[i]->m_object;
	}
	else {
		return NULL;
	}
}

template <class TObject, class TKey>
inline int TreeNode<TObject, TKey>::getSize() const
{
	int i, s;
	s = 1;
	for (i=0; i<m_children_num; i++) {
		if (m_children[i] != NULL) {
			s += m_children[i]->getSize();
		}
	}
	return s;
}

template <class TObject, class TKey>
inline void TreeNode<TObject, TKey>::release()
{
	int i;
	m_object = NULL;
	for (i=0; i<m_children_num; i++) {
		if (m_children[i] != NULL) {
			m_children[i]->release();
		}
	}
}


////////////////////////////////
//
// implementation of Tree

template <class TObject, class TKey>
inline Tree<TObject, TKey>::Tree(int child_num)
	: m_root(child_num)
{
	m_children_num = child_num;
}

template <class TObject, class TKey>
inline Tree<TObject, TKey>::Tree(const Tree &tree)
	: m_root(tree.m_root)
{
	m_children_num = tree.m_children_num;
}

template <class TObject, class TKey>
inline Tree<TObject, TKey>::~Tree()
{
}

template <class TObject, class TKey>
inline void Tree<TObject, TKey>::setChildrenNum(int child_num)
{
	if (m_children_num != child_num) {
		m_children_num = child_num;
		m_root.setChildrenNum(m_children_num);
	}
}

template <class TObject, class TKey>
inline TreeNode<TObject, TKey> *Tree<TObject, TKey>::newTreeNode(TObject *obj, TKey key)
{
	return (new TreeNode<TObject, TKey> (obj, key, m_children_num));
}


#endif //__TREE_H
