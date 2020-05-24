#define CATCH_CONFIG_MAIN

#include <GkylForwardAutoDiff.h>

#include <catch.hpp>
#include <cmath>
#include <vector>

TEST_CASE("Tests for relational operators", "[relops]") {
  double x = 5.0, y = 6.0;
  Gkyl::HyperDouble hx(5.0);
  Gkyl::HyperDouble hy(6.0);

  REQUIRE( x==hx );
  REQUIRE( hx==x );

  REQUIRE( x<hy );
  REQUIRE( hy>x );

  REQUIRE( y>hx );
  REQUIRE( hx<y );

  REQUIRE( hx<=hx );
  REQUIRE( hx>=hx );

  REQUIRE( hx<=5.0 );
  REQUIRE( hx>=5.0 );
}

// these set of test ensure that we can call the Gkyl::m math
// functions with regular POD types (double and float)
TEST_CASE("Tests for calling spec funcs with PODs", "[pod-func]") {
  double x = 5.0;
  double z;

  // f(x) = x*cos(x)
  z = x*Gkyl::cos(x);
  REQUIRE( z == Approx(5.0*std::cos(5.0)) );

  // f(x) = x*sin(x)
  z = x*Gkyl::sin(x);
  REQUIRE( z == Approx(5.0*std::sin(5.0)) );

  // f(x) = cos(x)^2 + sin(x)^2
  z = Gkyl::sin(x)*Gkyl::sin(x) + Gkyl::cos(x)*Gkyl::cos(x);
  REQUIRE( z == Approx(1.0) );

  // f(x) = cos(x*sin(x))
  z = Gkyl::cos(x*Gkyl::sin(x));
  REQUIRE( z == Approx(std::cos(5*std::sin(5))) );

  // f(x) = sqrt(x)
  z = Gkyl::sqrt(x);
  REQUIRE( z == Approx(std::sqrt(5)) );

  // f(x) = tan(x)
  z = Gkyl::tan(x);
  REQUIRE( z == Approx(std::tan(5)) );

  // f(x) = asin(t)
  double t = 0.5;
  z = Gkyl::asin(t);
  REQUIRE( z == Approx(std::asin(0.5)) );

  // f(x) = acos(t)
  t = 0.5;
  z = Gkyl::acos(t);
  REQUIRE( z == Approx(std::acos(0.5)) );

  // f(x) = atan(t)
  t = 0.5;
  z = Gkyl::atan(t);
  REQUIRE( z == Approx(std::atan(0.5)) );

  // sinh, cosh, tanh
  z = Gkyl::sinh(x);
  REQUIRE( z == Approx(std::sinh(5.0)) );
  z = Gkyl::cosh(x);
  REQUIRE( z == Approx(std::cosh(5.0)) );
  z = Gkyl::tanh(x);
  REQUIRE( z == Approx(std::tanh(5.0)) );

  // f(x) = exp(x)
  z = Gkyl::exp(x);
  REQUIRE( z == Approx(std::exp(5.0)) );

  // f(x) = log(x)
  z = Gkyl::log(x);
  REQUIRE( z == Approx(std::log(5.0)) );

  // f(x) = abs(x)
  z = Gkyl::abs(-x);
  REQUIRE( z == 5.0 );
}

TEST_CASE("Compound assignment tests", "[assignment]") {
  Gkyl::HyperDouble x(5.0, 1.0);
  Gkyl::HyperDouble z;

  // +=
  z = 0.0;
  z += 3;
  REQUIRE( z.real() == 3 );
  REQUIRE( z.inf() == 0.0 );

  z = 1.0;
  z += x;
  REQUIRE( z.real() == 6.0 );
  REQUIRE( z.inf() == 1.0 );

  // -=
  z = 0.0;
  z -= 3;
  REQUIRE( z.real() == -3 );
  REQUIRE( z.inf() == 0.0 );

  z = 1.0;
  z -= x;
  REQUIRE( z.real() == -4.0 );
  REQUIRE( z.inf() == -1.0 );

  // *=
  z = 1.5;
  z *= 3;
  REQUIRE( z.real() == 4.5 );
  REQUIRE( z.inf() == 0.0 );

  z = 1.0;
  z *= x;
  REQUIRE( z.real() == 5.0 );
  REQUIRE( z.inf() == 0.0 );

  // /=
  z = 1.0;
  z /= x;
  REQUIRE( z.real() == 0.2 );
  REQUIRE( z.inf() == 0.0 );

}

TEST_CASE("Basic first derivative tests", "[simple-first-diff]") {
  Gkyl::HyperDouble x(5.0, 1.0);
  Gkyl::HyperDouble z;

  // f(x) = 3
  z = 3;
  REQUIRE( z.real() == 3 );
  REQUIRE( z.inf() == 0.0 );

  // f(x) = -x+3
  z = -x+3;
  REQUIRE( z.real() == -5.0+3 );
  REQUIRE( z.inf() == -1.0 );

  // f(x) = 3-x
  z = 3-x;
  REQUIRE( z.real() == -5.0+3 );
  REQUIRE( z.inf() == -1.0 );

  // f(x) = x+3;
  z = 3; z += x;
  REQUIRE( z.real() == 5.0+3 );
  REQUIRE( z.inf() == 1.0 );

  // f(x) = 2*x^2
  z = 2*x*x;
  REQUIRE( z.real() == 2*25.0 );
  REQUIRE( z.inf() == 20.0 );

  // f(x) = x/(1+x)
  z = x/(1+x);
  REQUIRE( z.real() == Approx(5./(1+5)) );
  REQUIRE( z.inf() == Approx(1./(5*5+2*5+1)) );

  // f(x) = x*cos(x)
  z = x*Gkyl::cos(x);
  REQUIRE( z.real() == Approx(5.0*std::cos(5.0)) );
  REQUIRE( z.inf() == Approx(std::cos(5.0)-5.0*std::sin(5.0)) );

  // f(x) = x*sin(x)
  z = x*Gkyl::sin(x);
  REQUIRE( z.real() == Approx(5.0*std::sin(5.0)) );
  REQUIRE( z.inf() == Approx(std::sin(5.0)+5.0*std::cos(5.0)) );

  // f(x) = cos(x)^2 + sin(x)^2
  z = Gkyl::sin(x)*Gkyl::sin(x) + Gkyl::cos(x)*Gkyl::cos(x);
  REQUIRE( z.real() == Approx(1.0) );
  REQUIRE( z.inf() == Approx(0.0) );

  // f(x) = cos(x*sin(x))
  z = Gkyl::cos(x*Gkyl::sin(x));
  REQUIRE( z.real() == Approx(std::cos(5*std::sin(5))) );
  REQUIRE( z.inf() == Approx(-0.4578343032148585) );

  // f(x) = cos(w); w = x*sin(x)
  Gkyl::HyperDouble w = x*Gkyl::sin(x);
  z = Gkyl::cos(w);
  REQUIRE( z.real() == Approx(std::cos(5*std::sin(5))) );
  REQUIRE( z.inf() == Approx(-0.4578343032148585) );

  // f(x) = sqrt(x)
  z = Gkyl::sqrt(x);
  REQUIRE( z.real() == Approx(std::sqrt(5)) );
  REQUIRE( z.inf() == Approx(0.223606797749979) );

  // f(x) = tan(x)
  z = Gkyl::tan(x);
  REQUIRE( z.real() == Approx(std::tan(5)) );
  REQUIRE( z.inf() == Approx(12.42788170745835) );

  // f(x) = asin(x)
  Gkyl::HyperDouble t(0.5, 1.0);
  z = Gkyl::asin(t);
  REQUIRE( z.real() == Approx(std::asin(0.5)) );
  REQUIRE( z.inf() == Approx(1.154700538379252) );

  // f(x) = acos(x)
  z = Gkyl::acos(t);
  REQUIRE( z.real() == Approx(std::acos(0.5)) );
  REQUIRE( z.inf() == Approx(-1.154700538379252) );

  // f(x) = atan(x)
  z = Gkyl::atan(t);
  REQUIRE( z.real() == Approx(std::atan(0.5)) );
  REQUIRE( z.inf() == Approx(0.8) );

  // f(x) = sinh(x)
  z = Gkyl::sinh(x);
  REQUIRE( z.real() == Approx(std::sinh(5.0)) );
  REQUIRE( z.inf() == Approx(74.20994852478785) );

  // f(x) = cosh(x)
  z = Gkyl::cosh(x);
  REQUIRE( z.real() == Approx(std::cosh(5.0)) );
  REQUIRE( z.inf() == Approx(74.20321057778875) );

  // f(x) = tanh(x)
  z = Gkyl::tanh(x);
  REQUIRE( z.real() == Approx(std::tanh(5.0)) );
  REQUIRE( z.inf() == Approx(1.815832309438067e-4) );

  // f(x) = exp(x)
  z = Gkyl::exp(x);
  REQUIRE( z.real() == Approx(std::exp(5.0)) );
  REQUIRE( z.inf() == Approx(std::exp(5.0)) );

  // f(x) = log(x)
  z = Gkyl::log(x);
  REQUIRE( z.real() == Approx(std::log(5.0)) );
  REQUIRE( z.inf() == Approx(1/5.0) );

  // f(x) = abs(x)
  z = Gkyl::abs(Gkyl::HyperDouble(-5.0, 1.0));
  REQUIRE( z.real() == 5.0 );
  REQUIRE( z.inf() == -1.0 );

  // f(x) = x>5 then x*x else x*x*x
  z = x>5 ? x*x : x*x*x;
  REQUIRE( z.real() == 125.0 );
  REQUIRE( z.inf() == 75.0 );
}
