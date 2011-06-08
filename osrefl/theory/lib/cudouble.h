#ifndef _CUDOUBLE_H_
#define _CUDOUBLE_H_
#ifdef _CUFLOAT_H_
# warning cufloat.h and cudouble.h cannot coexist
#endif

typedef double Real;
const double eps = 2.22044604925e-16;
const double RealMax = 1.79769313486+308;

#include <lib/_cucomplex.h>

#endif /* _CUDOUBLE_H_ */
