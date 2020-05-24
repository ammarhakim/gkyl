// Gkyl ------------------------------------------------------------------------
//
// Forward-mode AD using HyperReal numbers and operator overloading
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------


#pragma once

// std includes
#include <cmath>
#include <iostream>

namespace Gkyl {
  
  template <typename RT, typename AT> class HyperReal;

  // Private types to extract real and adjoint parts from a
  // number. The number can be a POD (double/float) or a HyperReal

  /* Fetch real part of POD number */
  template <typename T>
  struct _R {
      static T g(T r) { return r; }
  };
  /* Fetch infinitesimal part of POD number */
  template <typename T>
  struct _I {
      static T g(T r) { return 0; }
  };
  
  /* Fetch real part of HyperReal number */
  template <typename RT, typename AT>
  struct _R<HyperReal<RT, AT> > {
      static RT g(const HyperReal<RT, AT>& r) { return r.real(); }
  };
  /* Fetch infinitesimal part of HyperReal number */
  template <typename RT, typename AT>
  struct _I<HyperReal<RT, AT> > {
      static AT g(const HyperReal<RT, AT>& r) { return r.inf(); }
  };

  /* Hyperreal number: real + infinitesimal (adjoint). RT is type of
   * the real-part and AT the type of the adoint part */
  template <typename RT, typename AT=RT>
  class HyperReal {
    public:
      // various ctors
      HyperReal() : rp(0), ip(0) { }
      HyperReal(const RT& rel) : rp(rel), ip(0) { }
      HyperReal(const RT& rel, const AT& inf) : rp(rel), ip(inf) { }

      // real and infinitesimal parts of number
      RT real() const { return rp; }
      AT inf() const { return ip; }

      // assignment =
      template<typename RHT>
      HyperReal& operator=(const RHT& rv) {
        RT y0 = _R<RHT>::g(rv); AT y1 = _I<RHT>::g(rv);
        rp = y0; ip = y1;
        return *this;    
      }

      // compound assignment +=
      template <typename RHT>
      HyperReal& operator+=(const RHT& rv) {
        RT y0 = _R<RHT>::g(rv); AT y1 = _I<RHT>::g(rv);
        rp += y0; ip += y1;
        return *this;
      }
      // compound assignment -=
      template <typename RHT>
      HyperReal& operator-=(const RHT& rv) {
        RT y0 = _R<RHT>::g(rv); AT y1 = _I<RHT>::g(rv);
        rp -= y0; ip -= y1;
        return *this;
      }
      // compound assignment *=
      template <typename RHT>
      HyperReal& operator*=(const RHT& rv) {
        RT y0 = _R<RHT>::g(rv); AT y1 = _I<RHT>::g(rv);
        rp *= y0; ip *= y1;
        return *this;
      }
      // compound assignment /=
      template <typename RHT>
      HyperReal& operator/=(const RHT& rv) {
        RT y0 = _R<RHT>::g(rv); AT y1 = _I<RHT>::g(rv);
        rp /= y0; ip /= y1;
        return *this;
      }

      // binary +
      template <typename LHT, typename RHT>
      friend HyperReal operator+(const LHT& lv, const RHT& rv) {
        RT x0 = _R<LHT>::g(lv), y0 = _R<RHT>::g(rv);
        AT x1 = _I<LHT>::g(lv), y1 = _I<RHT>::g(rv);
        return HyperReal<RT,AT>(x0+y0, x1+y1);
      }
      // binary -
      template <typename LHT, typename RHT>
      friend HyperReal operator-(const LHT& lv, const RHT& rv) {
        RT x0 = _R<LHT>::g(lv), y0 = _R<RHT>::g(rv);
        AT x1 = _I<LHT>::g(lv), y1 = _I<RHT>::g(rv);
        return HyperReal<RT,AT>(x0-y0, x1-y1);
      }
      // binary *
      template <typename LHT, typename RHT>
      friend HyperReal operator*(const LHT& lv, const RHT& rv) {
        RT x0 = _R<LHT>::g(lv), y0 = _R<RHT>::g(rv);
        AT x1 = _I<LHT>::g(lv), y1 = _I<RHT>::g(rv);
        return HyperReal<RT,AT>(x0*y0, x0*y1+x1*y0);
      }
      // binary /
      template <typename LHT, typename RHT>
      friend HyperReal operator/(const LHT& lv, const RHT& rv) {
        RT x0 = _R<LHT>::g(lv), y0 = _R<RHT>::g(rv);
        AT x1 = _I<LHT>::g(lv), y1 = _I<RHT>::g(rv);
        return HyperReal<RT,AT>(x0/y0, -(x0*y1-x1*y0)/(y0*y0));
      }

      // unary -, +
      HyperReal operator-() { return HyperReal<RT,AT>(-rp, -ip); }
      HyperReal operator+() { return HyperReal<RT,AT>(rp, ip); }

      // relational <
      template<typename LHT, typename RHT>
      friend bool operator<(const LHT& lv, const RHT& rv) {
        return _R<LHT>::g(lv) < _R<RHT>::g(rv);
      }
      // relational >
      template<typename LHT, typename RHT>
      friend bool operator>(const LHT& lv, const RHT& rv) { return rv < lv; }
      // relational <=
      template<typename LHT, typename RHT>
      friend bool operator<=(const LHT& lv, const RHT& rv) { return !(lv > rv); }
      // relational >=
      template<typename LHT, typename RHT>
      friend bool operator>=(const LHT& lv, const RHT& rv) { return !(lv < rv); }

      // equality == 
      template<typename LHT, typename RHT>
      friend bool operator==(const LHT& lv, const RHT& rv) {
        return _R<LHT>::g(lv) == _R<RHT>::g(rv);
      }
      // inequality != 
      template<typename LHT, typename RHT>
      friend bool operator!=(const LHT& lv, const RHT& rv) { return !(lv == rv); }
      
    private:
      RT rp; /* Real part */
      AT ip; /* Infinitesimal part */
  };

  // Predefined types
  using HyperDouble = HyperReal<double>;
  using HyperFloat = HyperReal<float>;

  /* Derivatives of functions from std::math library */

  // sign of value
  template <typename T>
  int sgn(T val) { return (T(0) < val) - (val < T(0)); }
    
  // this default private struct supplies methods for use with POD
  // types (double and float)
  template <typename T>
  struct _m {
      static T sqrt(const T& x) { return std::sqrt(x); }      
      static T cos(const T& x) { return std::cos(x); }
      static T sin(const T& x) { return std::sin(x); }
      static T tan(const T& x) { return std::tan(x); }
      static T asin(const T& x) { return std::asin(x); }
      static T acos(const T& x) { return std::acos(x); }
      static T atan(const T& x) { return std::atan(x); }
      static T sinh(const T& x) { return std::sinh(x); }
      static T cosh(const T& x) { return std::cosh(x); }
      static T tanh(const T& x) { return std::tanh(x); }
      static T exp(const T& x) { return std::exp(x); }
      static T log(const T& x) { return std::log(x); }
      static T abs(const T& x) { return std::abs(x); }
      static T floor(const T& x) { return std::floor(x); }
      static T ceil(const T& x) { return std::ceil(x); }
  };
  
  // specialization to HyperReal number
  template <typename RT, typename AT>
  struct _m<HyperReal<RT,AT> > {
      
      static HyperReal<RT,AT> sqrt(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        double y0 = std::sqrt(x0);
        return HyperReal<RT,AT>(y0, 0.5*x1/y0);
      }
      
      static HyperReal<RT,AT> cos(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        return HyperReal<RT,AT>(std::cos(x0), -x1*std::sin(x0));
      }
        
      static HyperReal<RT,AT> sin(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        return HyperReal<RT,AT>(std::sin(x0), x1*std::cos(x0));
      }

      static HyperReal<RT,AT> tan(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        double tx0 = std::tan(x0);
        return HyperReal<RT,AT>(tx0, x1*(1+tx0*tx0));
      }

      static HyperReal<RT,AT> asin(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        return HyperReal<RT,AT>(std::asin(x0), x1/std::sqrt(1-x0*x0));
      }

      static HyperReal<RT,AT> acos(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        return HyperReal<RT,AT>(std::acos(x0), -x1/std::sqrt(1-x0*x0));
      }

      static HyperReal<RT,AT> atan(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        return HyperReal<RT,AT>(std::atan(x0), x1/(1+x0*x0));
      }

      static HyperReal<RT,AT> sinh(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        return HyperReal<RT,AT>(std::sinh(x0), x1*std::cosh(x0));
      }

      static HyperReal<RT,AT> cosh(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        return HyperReal<RT,AT>(std::cosh(x0), x1*std::sinh(x0));
      }

      static HyperReal<RT,AT> tanh(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        double tx0 = std::tanh(x0);
        return HyperReal<RT,AT>(tx0, x1*(1-tx0*tx0));
      }

      static HyperReal<RT,AT> exp(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        double ex0 = std::exp(x0);
        return HyperReal<RT,AT>(ex0, x1*ex0);
      }

      static HyperReal<RT,AT> log(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        return HyperReal<RT,AT>(std::log(x0), x1/x0);
      }

      static HyperReal<RT,AT> abs(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        return HyperReal<RT,AT>(std::abs(x0), x1*sgn(x0));
      }

      static HyperReal<RT,AT> floor(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        return HyperReal<RT,AT>(std::floor(x0), 0.0);
      }

      static HyperReal<RT,AT> ceil(const HyperReal<RT,AT>& x) {
        RT x0 = x.real(); AT x1 = x.inf();
        return HyperReal<RT,AT>(std::ceil(x0), 0.0);
      }
  };

  // Math functions
  template <typename T> inline T sqrt(const T& x) { return _m<T>::sqrt(x); }
  template <typename T> inline T cos(const T& x) { return _m<T>::cos(x); }
  template <typename T> inline T sin(const T& x) { return _m<T>::sin(x); }
  template <typename T> inline T tan(const T& x) { return _m<T>::tan(x); }
  template <typename T> inline T asin(const T& x) { return _m<T>::asin(x); }
  template <typename T> inline T acos(const T& x) { return _m<T>::acos(x); }
  template <typename T> inline T atan(const T& x) { return _m<T>::atan(x); }
  template <typename T> inline T sinh(const T& x) { return _m<T>::sinh(x); }
  template <typename T> inline T cosh(const T& x) { return _m<T>::cosh(x); }
  template <typename T> inline T tanh(const T& x) { return _m<T>::tanh(x); }
  template <typename T> inline T exp(const T& x) { return _m<T>::exp(x); }
  template <typename T> inline T log(const T& x) { return _m<T>::log(x); }
  template <typename T> inline T abs(const T& x) { return _m<T>::abs(x); }
  template <typename T> inline T floor(const T& x) { return _m<T>::floor(x); }
  template <typename T> inline T ceil(const T& x) { return _m<T>::ceil(x); }
}
