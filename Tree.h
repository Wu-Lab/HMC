
#ifndef __TREE_H
#define __TREE_H


#include <vector>


template <class T>
class TreeNode {
	TreeNode *m_parent;
	vector<TreeNode*> m_children;
	T m_data;

public:
	explicit TreeNode(int n = 0, T d = T()) : m_parent(0), m_children(n), m_data(d) { }
	TreeNode(const TreeNode &node);
	~TreeNode();

	T data() const { return m_data; }
	int size() const { return m_children.size(); }
	void resize(int n);

	TreeNode *addChild(int i, T d = T());
	TreeNode *setChild(int i, T d);
	TreeNode *getChild(int i) const { return m_children[i]; }
};


template <class T>
inline TreeNode<T>::TreeNode(const TreeNode &node)
: m_parent(node.m_parent),
  m_children(node.m_children.size()),
  m_data(node.m_data)
{
	for (int i=0; i<size(); ++i) {
		if (node.m_children[i]) {
			m_children[i] = new TreeNode (*node.m_children[i]);
			m_children[i]->m_parent = this;
		}
	}
}

template <class T>
inline TreeNode<T>::~TreeNode()
{
	for (int i=0; i<size(); ++i) {
		delete m_children[i];
	}
}

template <class T>
inline void TreeNode<T>::resize(int n)
{
	int i;
	for (i=n; i<size(); ++i) {
		delete m_children[i];
	}
	m_children.resize(n);
	for (i=0; i<size(); ++i) {
		if (m_children[i]) {
			m_children[i]->resize(n);
		}
	}
}

template <class T>
inline TreeNode<T> *TreeNode<T>::addChild(int i, T d)
{
	if (!m_children[i]) {
		m_children[i] = new TreeNode(size(), d);
		m_children[i]->m_parent = this;
	}
	return m_children[i];
}

template <class T>
inline TreeNode<T> *TreeNode<T>::setChild(int i, T d)
{
	if (m_children[i]) {
		m_children[i]->m_data = d;
	}
	else {
		m_children[i] = new TreeNode(size(), d);
		m_children[i]->m_parent = this;
	}
	return m_children[i];
}


#endif //__TREE_H
