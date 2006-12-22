
#ifndef __MATRIX_H
#define __MATRIX_H


template <class T>
class Matrix {
public:
	static T **alloc(int row, int col);
	static T **alloc(int row_l, int row_h, int col_l, int col_h);
	static void free(T **m, int row_l = 0, int col_l = 0);
};


////////////////////////////////
//
// implementations of Matrix


template <class T>
inline T **Matrix<T>::alloc(int row, int col)
{
	T **m =	alloc(0, row-1, 0, col-1);
	return m;
}

template <class T>
inline T **Matrix<T>::alloc(int row_l, int row_h, int col_l, int col_h)
{
	T **m;
	long i, nrow=row_h-row_l+1, ncol=col_h-col_l+1;
	m = new T * [nrow];
	m -= row_l;
	m[row_l] = new T [nrow*ncol];
	m[row_l] -= col_l;
	for(i=row_l+1; i<=row_h; i++) m[i] = m[i-1] + ncol;
	return m;
}

template <class T>
inline void Matrix<T>::free(T **m, int row_l, int col_l)
{
	if (m != NULL) {
		delete[] (m[row_l]+col_l);
		delete[] (m+row_l);
	}
	m = NULL;
}


#endif // __MATRIX_H
