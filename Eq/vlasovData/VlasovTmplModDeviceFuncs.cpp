// Gkyl ------------------------------------------------------------------------
//
// Code to copy various kernerls from device to host
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <VlasovTmplModDecl.h>

using namespace Gkyl;



__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x1vMaxP1 =
  &Gkyl::VlasovModDecl<1,1,1,G_MAX_ORDER_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_1x1vMaxP1 =
  &Gkyl::VlasovModDecl<1,1,1,G_MAX_ORDER_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x1vMaxP1 =
  &Gkyl::VlasovModDecl<1,1,1,G_MAX_ORDER_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x1vMaxP1 =
  &Gkyl::VlasovModDecl<1,1,1,G_MAX_ORDER_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x1vMaxP2 =
  &Gkyl::VlasovModDecl<1,1,2,G_MAX_ORDER_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_1x1vMaxP2 =
  &Gkyl::VlasovModDecl<1,1,2,G_MAX_ORDER_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x1vMaxP2 =
  &Gkyl::VlasovModDecl<1,1,2,G_MAX_ORDER_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x1vMaxP2 =
  &Gkyl::VlasovModDecl<1,1,2,G_MAX_ORDER_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x2vMaxP1 =
  &Gkyl::VlasovModDecl<1,2,1,G_MAX_ORDER_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_1x2vMaxP1 =
  &Gkyl::VlasovModDecl<1,2,1,G_MAX_ORDER_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x2vMaxP1 =
  &Gkyl::VlasovModDecl<1,2,1,G_MAX_ORDER_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x2vMaxP1 =
  &Gkyl::VlasovModDecl<1,2,1,G_MAX_ORDER_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x2vMaxP2 =
  &Gkyl::VlasovModDecl<1,2,2,G_MAX_ORDER_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_1x2vMaxP2 =
  &Gkyl::VlasovModDecl<1,2,2,G_MAX_ORDER_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x2vMaxP2 =
  &Gkyl::VlasovModDecl<1,2,2,G_MAX_ORDER_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x2vMaxP2 =
  &Gkyl::VlasovModDecl<1,2,2,G_MAX_ORDER_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x3vMaxP1 =
  &Gkyl::VlasovModDecl<1,3,1,G_MAX_ORDER_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_1x3vMaxP1 =
  &Gkyl::VlasovModDecl<1,3,1,G_MAX_ORDER_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x3vMaxP1 =
  &Gkyl::VlasovModDecl<1,3,1,G_MAX_ORDER_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x3vMaxP1 =
  &Gkyl::VlasovModDecl<1,3,1,G_MAX_ORDER_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x3vMaxP2 =
  &Gkyl::VlasovModDecl<1,3,2,G_MAX_ORDER_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_1x3vMaxP2 =
  &Gkyl::VlasovModDecl<1,3,2,G_MAX_ORDER_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x3vMaxP2 =
  &Gkyl::VlasovModDecl<1,3,2,G_MAX_ORDER_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x3vMaxP2 =
  &Gkyl::VlasovModDecl<1,3,2,G_MAX_ORDER_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_2x2vMaxP1 =
  &Gkyl::VlasovModDecl<2,2,1,G_MAX_ORDER_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_2x2vMaxP1 =
  &Gkyl::VlasovModDecl<2,2,1,G_MAX_ORDER_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_2x2vMaxP1 =
  &Gkyl::VlasovModDecl<2,2,1,G_MAX_ORDER_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_2x2vMaxP1 =
  &Gkyl::VlasovModDecl<2,2,1,G_MAX_ORDER_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_2x2vMaxP2 =
  &Gkyl::VlasovModDecl<2,2,2,G_MAX_ORDER_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_2x2vMaxP2 =
  &Gkyl::VlasovModDecl<2,2,2,G_MAX_ORDER_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_2x2vMaxP2 =
  &Gkyl::VlasovModDecl<2,2,2,G_MAX_ORDER_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_2x2vMaxP2 =
  &Gkyl::VlasovModDecl<2,2,2,G_MAX_ORDER_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_2x3vMaxP1 =
  &Gkyl::VlasovModDecl<2,3,1,G_MAX_ORDER_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_2x3vMaxP1 =
  &Gkyl::VlasovModDecl<2,3,1,G_MAX_ORDER_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_2x3vMaxP1 =
  &Gkyl::VlasovModDecl<2,3,1,G_MAX_ORDER_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_2x3vMaxP1 =
  &Gkyl::VlasovModDecl<2,3,1,G_MAX_ORDER_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_2x3vMaxP2 =
  &Gkyl::VlasovModDecl<2,3,2,G_MAX_ORDER_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_2x3vMaxP2 =
  &Gkyl::VlasovModDecl<2,3,2,G_MAX_ORDER_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_2x3vMaxP2 =
  &Gkyl::VlasovModDecl<2,3,2,G_MAX_ORDER_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_2x3vMaxP2 =
  &Gkyl::VlasovModDecl<2,3,2,G_MAX_ORDER_C>::surfElcMagTerm;




__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_3x3vMaxP1 =
  &Gkyl::VlasovModDecl<3,3,1,G_MAX_ORDER_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_3x3vMaxP1 =
  &Gkyl::VlasovModDecl<3,3,1,G_MAX_ORDER_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_3x3vMaxP1 =
  &Gkyl::VlasovModDecl<3,3,1,G_MAX_ORDER_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_3x3vMaxP1 =
  &Gkyl::VlasovModDecl<3,3,1,G_MAX_ORDER_C>::surfElcMagTerm;




__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x1vSerP1 =
  &Gkyl::VlasovModDecl<1,1,1,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_1x1vSerP1 =
  &Gkyl::VlasovModDecl<1,1,1,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x1vSerP1 =
  &Gkyl::VlasovModDecl<1,1,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x1vSerP1 =
  &Gkyl::VlasovModDecl<1,1,1,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x1vSerP2 =
  &Gkyl::VlasovModDecl<1,1,2,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_1x1vSerP2 =
  &Gkyl::VlasovModDecl<1,1,2,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x1vSerP2 =
  &Gkyl::VlasovModDecl<1,1,2,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x1vSerP2 =
  &Gkyl::VlasovModDecl<1,1,2,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x2vSerP1 =
  &Gkyl::VlasovModDecl<1,2,1,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_1x2vSerP1 =
  &Gkyl::VlasovModDecl<1,2,1,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x2vSerP1 =
  &Gkyl::VlasovModDecl<1,2,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x2vSerP1 =
  &Gkyl::VlasovModDecl<1,2,1,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x2vSerP2 =
  &Gkyl::VlasovModDecl<1,2,2,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_1x2vSerP2 =
  &Gkyl::VlasovModDecl<1,2,2,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x2vSerP2 =
  &Gkyl::VlasovModDecl<1,2,2,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x2vSerP2 =
  &Gkyl::VlasovModDecl<1,2,2,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x3vSerP1 =
  &Gkyl::VlasovModDecl<1,3,1,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_1x3vSerP1 =
  &Gkyl::VlasovModDecl<1,3,1,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x3vSerP1 =
  &Gkyl::VlasovModDecl<1,3,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x3vSerP1 =
  &Gkyl::VlasovModDecl<1,3,1,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x3vSerP2 =
  &Gkyl::VlasovModDecl<1,3,2,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_1x3vSerP2 =
  &Gkyl::VlasovModDecl<1,3,2,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x3vSerP2 =
  &Gkyl::VlasovModDecl<1,3,2,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x3vSerP2 =
  &Gkyl::VlasovModDecl<1,3,2,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_2x2vSerP1 =
  &Gkyl::VlasovModDecl<2,2,1,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_2x2vSerP1 =
  &Gkyl::VlasovModDecl<2,2,1,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_2x2vSerP1 =
  &Gkyl::VlasovModDecl<2,2,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_2x2vSerP1 =
  &Gkyl::VlasovModDecl<2,2,1,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_2x2vSerP2 =
  &Gkyl::VlasovModDecl<2,2,2,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_2x2vSerP2 =
  &Gkyl::VlasovModDecl<2,2,2,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_2x2vSerP2 =
  &Gkyl::VlasovModDecl<2,2,2,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_2x2vSerP2 =
  &Gkyl::VlasovModDecl<2,2,2,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_2x3vSerP1 =
  &Gkyl::VlasovModDecl<2,3,1,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_2x3vSerP1 =
  &Gkyl::VlasovModDecl<2,3,1,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_2x3vSerP1 =
  &Gkyl::VlasovModDecl<2,3,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_2x3vSerP1 =
  &Gkyl::VlasovModDecl<2,3,1,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_2x3vSerP2 =
  &Gkyl::VlasovModDecl<2,3,2,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_2x3vSerP2 =
  &Gkyl::VlasovModDecl<2,3,2,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_2x3vSerP2 =
  &Gkyl::VlasovModDecl<2,3,2,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_2x3vSerP2 =
  &Gkyl::VlasovModDecl<2,3,2,G_SERENDIPITY_C>::surfElcMagTerm;




__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_3x3vSerP1 =
  &Gkyl::VlasovModDecl<3,3,1,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_3x3vSerP1 =
  &Gkyl::VlasovModDecl<3,3,1,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_3x3vSerP1 =
  &Gkyl::VlasovModDecl<3,3,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_3x3vSerP1 =
  &Gkyl::VlasovModDecl<3,3,1,G_SERENDIPITY_C>::surfElcMagTerm;

