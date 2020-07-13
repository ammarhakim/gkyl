#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include <GkylBasisTypes.h>

template <int N> class CompileTimeCheck {
  public:
    static int val() { return N; }
};


TEST_CASE("Tests for BasisTypes[m]", "[max-order]") {

  // 1D
  REQUIRE( Gkyl::BasisCount<1,0,Gkyl::G_MAX_ORDER_C>::numBasis() == 1 );
  REQUIRE( Gkyl::BasisCount<1,1,Gkyl::G_MAX_ORDER_C>::numBasis() == 2 );
  REQUIRE( Gkyl::BasisCount<1,2,Gkyl::G_MAX_ORDER_C>::numBasis() == 3 );

  REQUIRE( Gkyl::BasisCount<1,2,Gkyl::G_MAX_ORDER_C>::id() == "maximal-order" );

  // 2D
  REQUIRE( Gkyl::BasisCount<2,0,Gkyl::G_MAX_ORDER_C>::numBasis() == 1 );
  REQUIRE( Gkyl::BasisCount<2,1,Gkyl::G_MAX_ORDER_C>::numBasis() == 3 );
  REQUIRE( Gkyl::BasisCount<2,2,Gkyl::G_MAX_ORDER_C>::numBasis() == 6 );
}

TEST_CASE("Tests for BasisTypes[s]", "[serendipity]") {

  // 1D
  REQUIRE( Gkyl::BasisCount<1,0,Gkyl::G_SERENDIPITY_C>::numBasis() == 1 );
  REQUIRE( Gkyl::BasisCount<1,1,Gkyl::G_SERENDIPITY_C>::numBasis() == 2 );
  REQUIRE( Gkyl::BasisCount<1,2,Gkyl::G_SERENDIPITY_C>::numBasis() == 3 );

  // 2D
  REQUIRE( Gkyl::BasisCount<2,0,Gkyl::G_SERENDIPITY_C>::numBasis() == 1 );
  REQUIRE( Gkyl::BasisCount<2,1,Gkyl::G_SERENDIPITY_C>::numBasis() == 4 );
  REQUIRE( Gkyl::BasisCount<2,2,Gkyl::G_SERENDIPITY_C>::numBasis() == 8 );

  REQUIRE( Gkyl::BasisCount<1,2,Gkyl::G_SERENDIPITY_C>::id() == "serendipity" );
}

TEST_CASE("Tests if BasisTypes is compile-time", "[compile-time]") {
  
  //std::cout << Gkyl::BasisCount<1,2,Gkyl::G_SERENDIPITY_C>::numBasis() << std::endl;
  REQUIRE(
    CompileTimeCheck<Gkyl::BasisCount<1,2,Gkyl::G_SERENDIPITY_C>::numBasis()>::val()
    ==
    Gkyl::BasisCount<1,2,Gkyl::G_SERENDIPITY_C>::numBasis()
  );
}
