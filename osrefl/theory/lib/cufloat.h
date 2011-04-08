#ifndef _CUFLOAT_H_
#define _CUFLOAT_H_
#ifdef _CUDOUBLE_H_
# warning cufloat.h and cudouble.h cannot coexist
#endif

typedef float Real;
const float eps = 1.19209289551e-07f;
const float RealMax = 3.402823466e38f;

#include <_cucomplex.h>

#endif /* _CUFLOAT_H_ */
