
#ifndef __MEMLEAK_H
#define __MEMLEAK_H


// detect memory leak

#ifdef _DEBUG
#include <cstdlib>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK, __FILE__, __LINE__)
#endif // _DEBUG


inline void EnableMemLeakCheck()
{
#ifdef _DEBUG
	_CrtSetDbgFlag(_CrtSetDbgFlag(_CRTDBG_REPORT_FLAG) | _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif // _DEBUG
}


#endif // __MEMLEAK_H
