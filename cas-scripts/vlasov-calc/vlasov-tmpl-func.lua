local xsys = require "xsys"

local vlasovHeaderTopString = [[
// Gkyl ------------------------------------------------------------------------
//
// Code to copy various kernerls from device to host
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <VlasovTmplModDecl.h>

using namespace Gkyl;
]]

vlasovHeaderTemplateString = [[


|for ci = CMIN, CMAX do
| for vi = ci, VMAX do
|  for pi = PMIN, PMAX do

__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_${ci}x${vi}v${basisShortNm}P${pi} =
  &Gkyl::VlasovModDecl<${ci},${vi},${pi},${basisNm}>::volumeStreamTerm;

__device__ Vlasov_surfSreamTerm_t p_Vlasov_surfSreamTerm_${ci}x${vi}v${basisShortNm}P${pi} =
  &Gkyl::VlasovModDecl<${ci},${vi},${pi},${basisNm}>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_${ci}x${vi}v${basisShortNm}P${pi} =
  &Gkyl::VlasovModDecl<${ci},${vi},${pi},${basisNm}>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_${ci}x${vi}v${basisShortNm}P${pi} =
  &Gkyl::VlasovModDecl<${ci},${vi},${pi},${basisNm}>::surfElcMagTerm;

|    end
|  end
|end
]]

-- instantiate template
local vlasovHeaderTemplate = xsys.template ( vlasovHeaderTemplateString )

-- concatinate various pieces to generate final header
vlasovHeaderOut =
   vlasovHeaderTopString
   ..
   vlasovHeaderTemplate {
      CMIN = 1, CMAX = 2,
      VMAX = 3,
      PMIN = 1, PMAX = 2,
      basisNm = 'G_MAX_ORDER_C',
      basisShortNm = 'Max',
   }
   ..
   vlasovHeaderTemplate {
      CMIN = 3, CMAX = 3,
      VMAX = 3,
      PMIN = 1, PMAX = 1,
      basisNm = 'G_MAX_ORDER_C',
      basisShortNm = 'Max',
   }
   ..
   vlasovHeaderTemplate {
      CMIN = 1, CMAX = 2,
      VMAX = 3,
      PMIN = 1, PMAX = 2,
      basisNm = 'G_SERENDIPITY_C',
      basisShortNm = 'Ser',
   }
   ..
   vlasovHeaderTemplate {
      CMIN = 3, CMAX = 3,
      VMAX = 3,
      PMIN = 1, PMAX = 1,
      basisNm = 'G_SERENDIPITY_C',
      basisShortNm = 'Ser',
   }

-- write out string
io.write( vlasovHeaderOut )
