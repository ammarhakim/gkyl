#define CATCH_CONFIG_MAIN

#include <catch.hpp>
#include <GkylRange.h>
#include <cmath>
#include <vector>

TEST_CASE("Tests for range object", "[range]") {
  GkylRange_t crange;
  crange.ndim = 2;
  crange.lower[0] = 1; crange.upper[0] = 5;
  crange.lower[1] = 2; crange.upper[1] = 10;

  Gkyl::Range range(&crange);

  REQUIRE( range.ndim() == 2 );
  for (unsigned i=0; i<range.ndim(); ++i) {
    REQUIRE( range.lower(i) == crange.lower[i] );
    REQUIRE( range.upper(i) == crange.upper[i] );
    REQUIRE( range.shape(i) == (crange.upper[i]-crange.lower[i]+1) );
  }
  REQUIRE( range.volume() == 5*9 );
}

