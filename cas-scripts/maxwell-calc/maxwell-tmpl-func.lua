local xsys = require "xsys"

local maxwellHeaderTopString = [[
// Gkyl ------------------------------------------------------------------------
//
// Code to copy various kernels from device to host
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <MaxwellTmplModDecl.h>

using namespace Gkyl;
]]

maxwellHeaderTemplateString = [[


|for ci = CMIN, CMAX do
|  for pi = PMIN, PMAX do

__device__ Maxwell_volumeTerm_t p_Maxwell_volumeTerm_${ci}x${basisShortNm}P${pi} =
  &Gkyl::MaxwellModDecl<${ci},${pi},${basisNm}>::volumeTerm;

__device__ Maxwell_surfTerm_t p_Maxwell_surfTerm_${ci}x${basisShortNm}P${pi} =
  &Gkyl::MaxwellModDecl<${ci},${pi},${basisNm}>::surfTerm;

|    end
|end
]]

-- instantiate template
local maxwellHeaderTemplate = xsys.template ( maxwellHeaderTemplateString )

-- concatinate various pieces to generate final header
maxwellHeaderOut =
   maxwellHeaderTopString
   ..
   maxwellHeaderTemplate {
      CMIN = 1, CMAX = 2,
      PMIN = 1, PMAX = 2,
      basisNm = 'G_MAX_ORDER_C',
      basisShortNm = 'Max',
   }
   ..
   maxwellHeaderTemplate {
      CMIN = 3, CMAX = 3,
      PMIN = 1, PMAX = 1,
      basisNm = 'G_MAX_ORDER_C',
      basisShortNm = 'Max',
   }
   ..
   maxwellHeaderTemplate {
      CMIN = 1, CMAX = 2,
      PMIN = 1, PMAX = 2,
      basisNm = 'G_SERENDIPITY_C',
      basisShortNm = 'Ser',
   }
   ..
   maxwellHeaderTemplate {
      CMIN = 3, CMAX = 3,
      PMIN = 1, PMAX = 1,
      basisNm = 'G_SERENDIPITY_C',
      basisShortNm = 'Ser',
   }


maxwellFuncCopyTemplateString = [[

Maxwell_volumeTerm_t
Maxwell_getVolumeTerm(int cdim, int polyOrder, int basisType) {
  Maxwell_volumeTerm_t func;
  
  if (cdim == 1) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_volumeTerm_1xSerP1, sizeof(Maxwell_volumeTerm_t));
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_volumeTerm_1xSerP2, sizeof(Maxwell_volumeTerm_t));
      }
  }
  else if (cdim == 2) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_volumeTerm_2xSerP1, sizeof(Maxwell_volumeTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_volumeTerm_2xSerP2, sizeof(Maxwell_volumeTerm_t));
      }      
  }
  else if (cdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_volumeTerm_3xSerP1, sizeof(Maxwell_volumeTerm_t));
      }
  }
  return func;
}

Maxwell_surfTerm_t
Maxwell_getSurfTerm(int cdim, int polyOrder, int basisType) {
  Maxwell_surfTerm_t func;
  
  if (cdim == 1) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_surfTerm_1xSerP1, sizeof(Maxwell_surfTerm_t));
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_surfTerm_1xSerP2, sizeof(Maxwell_surfTerm_t));
      }
  }
  else if (cdim == 2) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_surfTerm_2xSerP1, sizeof(Maxwell_surfTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_surfTerm_2xSerP2, sizeof(Maxwell_surfTerm_t));
      }      
  }
  else if (cdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_surfTerm_3xSerP1, sizeof(Maxwell_surfTerm_t));
      }
  }
  return func;
}
]]


-- write out string
io.write( maxwellHeaderOut .. maxwellFuncCopyTemplateString )
