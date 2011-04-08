#if !defined(_CUFLOAT_H_) && !defined(_CUDOUBLE_H_)
# warning _cucomplex.h should not be directly included
#endif

#ifdef HAVE_CUDA
# include <cuda.h>
# include <math_functions.h>
# include <cuComplex.h>
# define INLINE __device__ __inline__
# define CUDA_KERNEL extern "C" __global__ void
#else
# include <cmath>
# define INLINE inline
# define CUDA_KERNEL extern "C" void
#endif

class Cplx
{
private:
  Real _real, _imag;
public:
  INLINE Cplx(const Real& r= 0, const Real& i= 0) : _real(r), _imag(i) { }
  INLINE Cplx(const Cplx& c) : _real(c._real), _imag(c._imag) { }
  INLINE const Real& real() const { return _real; }
  INLINE const Real& imag() const { return _imag; }
  INLINE Real& real() { return _real; }
  INLINE Real& imag() { return _imag; }
  INLINE void real(const Real& r) { _real = r; }
  INLINE void imag(const Real& r) { _imag = r; }


  // Real args
  INLINE Cplx& operator*=(const Real& r)
    { _real*=r; _imag*=r; return *this; }
  INLINE Cplx& operator/=(const Real& r)
    { _real/=r; _imag/=r; return *this; }
  INLINE Cplx& operator+=(const Real& r)
    { _real+=r; return *this; }
  INLINE Cplx& operator-=(const Real& r)
    { _real-=r; return *this; }
  INLINE Cplx& operator=(const Real& r)
    { _real=r; _imag=0; return *this; }

  // Complex args
  INLINE Cplx& operator+=(const Cplx& r)
    { _real += r._real; _imag += r._imag; return *this; }
  INLINE Cplx& operator-=(const Cplx& r)
    { _real -= r._real; _imag -= r._imag; return *this; }
  INLINE Cplx& operator*=(const Cplx& r)
    {
      const Real u = _real*r._real - _imag*r._imag;
      const Real v = _real*r._imag + _imag*r._real;
//std::cout<<"*"<<u<<" "<<v<<std::endl;
      _real = u; _imag = v; return *this;
    }
  INLINE Cplx& operator/=(const Cplx& r)
    {
      /* Robert L. Smith, Algorithm 116: Complex division,
       * Communications of the ACM, v.5 n.8, p.435, Aug. 1962
       * */
      if (fabs(r._real) >= fabs(r._imag)) {
        const Real t = r._imag/r._real;
        const Real den = r._real + r._imag*t;
        const Real u = (_real + _imag*t)/den;
        const Real v = (_imag - _real*t)/den;
//std::cout<<"/"<<u<<" "<<v<<std::endl;
        _real = u; _imag = v;

      }else {
          const Real t = r._real/r._imag;
          const Real den = r._real*t + r._imag;
          const Real u = (_real*t + _imag)/den;
          const Real v = (_imag*t - _real)/den;
  //std::cout<<"/"<<u<<" "<<v<<std::endl;
          _real = u; _imag = v;
        }

      return *this;
    }
  INLINE Cplx& operator=(const Cplx& r)
    { _real = r._real; _imag = r._imag; return *this; }
} ;

/*
#ifndef HAVE_SINCOS
INLINE void sincos(const Real v, Real *s, Real *c)
{
  *s = sin(v);
  *c = cos(v);
}
#endif
*/

// Operators
INLINE Cplx operator-(const Cplx& x)
  { return Cplx(-x.real(),-x.imag()); }

INLINE Cplx operator+(const Cplx& x, const Cplx& y)
  { Cplx r = x; r += y; return r; }
INLINE Cplx operator-(const Cplx& x, const Cplx& y)
  { Cplx r = x; r -= y; return r; }
INLINE Cplx operator*(const Cplx& x, const Cplx& y)
  { Cplx r = x; r *= y; return r; }
INLINE Cplx operator/(const Cplx& x, const Cplx& y)
  { Cplx r = x; r /= y; return r; }

INLINE Cplx operator+(const Cplx& x, const Real& y)
  { Cplx r = x; r += y; return r; }
INLINE Cplx operator-(const Cplx& x, const Real& y)
  { Cplx r = x; r -= y; return r; }
INLINE Cplx operator*(const Cplx& x, const Real& y)
  { Cplx r = x; r *= y; return r; }
INLINE Cplx operator/(const Cplx& x, const Real& y)
  { Cplx r = x; r /= y; return r; }

INLINE Cplx operator+(const Real& x, const Cplx& y)
  { Cplx r = y; r += x; return r; }
INLINE Cplx operator-(const Real& x, const Cplx& y)
  { Cplx r = -y; r += x; return r; }
INLINE Cplx operator*(const Real& x, const Cplx& y)
  { Cplx r = x; r *= y; return r; }
INLINE Cplx operator/(const Real& x, const Cplx& y)
  { Cplx r = x; r /= y; return r; }


// Functions
INLINE Real abs(const Cplx& Q)
{
  const Real a = fabs(Q.real());
  const Real b = fabs(Q.imag());
  Real v, w, t;
  if (a > b) { v = a; w = b; }
  else { v = b; w = a; }
  t = w/v;
  t = Real(1) + t*t;
  t = v * sqrt(t);
  if ((v == 0) || (v > RealMax) || (w > RealMax)) { t = v+w; }
  return t;
}
INLINE Real norm(const Cplx& z)
{
  const Real a = z.real();
  const Real b = z.imag();
  return a*a + b*b;
}
INLINE Cplx conj(const Cplx& z)
  { return Cplx(z.real(), -z.imag()); }
INLINE Cplx polar(const Real& rho, const Real& theta)
{
  Real c,s;
  sincos(theta,&c,&s);
  return Cplx(rho*c, rho*s);
}

INLINE Cplx sqrt(const Cplx& z)
{
  const Real a = z.real();
  const Real b = z.imag();
//std::cout<<"sqrt"<<a<<" "<<b<<std::endl;
  if (a == 0) {
    const Real real = sqrt(fabs(b) / 2);
    const Real imag = a < 0 ? -real : real;
    return Cplx(real, imag);
  } else {
    const Real t = sqrt(2 * (abs(z) + fabs(a)));
    const Real u = t/2;
    if (a > 0) {
      return Cplx(u,z.imag()/t);
    } else {
      const Real real = fabs(b)/t;
      const Real imag = (b<0 ? -u : u);
      return Cplx(real,imag);
    }
  }
}

INLINE Cplx pow(const Cplx& z, int n)
{
	Cplx zn(z), result(1,0);
	while (n) {
		if (n&1) result *= zn;
		zn *= zn;
		n >>= 1;
	}
	return result;
}

INLINE Cplx exp(const Cplx& z)
{
  const Real r = exp(z.real());
  Real s, c;
  sincos(z.imag(),&s,&c);
  return Cplx(r*c, r*s);
}
INLINE Cplx sin(const Cplx& z)
{
  const Real a = z.real();
  const Real b = z.imag();
  return Cplx(sin(a)*cosh(b), -cos(a)*sinh(b));
}
INLINE Cplx sinh(const Cplx& z)
{
  const Real a = z.real();
  const Real b = z.imag();
  return Cplx(sinh(a)*cos(b), -cosh(a)*sin(b));
}
INLINE Cplx cos(const Cplx& z)
{
  const Real a = z.real();
  const Real b = z.imag();
  return Cplx(cos(a)*cosh(b), -sin(a)*sinh(b));
}
INLINE Cplx cosh(const Cplx& z)
{
  const Real a = z.real();
  const Real b = z.imag();
  return Cplx(cosh(a)*cos(b), -sinh(a)*sin(b));
}
INLINE Cplx tan(const Cplx& z)
  { return sin(z)/cos(z); }
INLINE Cplx tanh(const Cplx& z)
  { return sinh(z)/cosh(z); }
