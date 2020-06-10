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

__device__ Vlasov_surfStreamTerm_t p_Vlasov_surfStreamTerm_${ci}x${vi}v${basisShortNm}P${pi} =
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


vlasovFuncCopyTemplateString = [[

Vlasov_volumeStreamTerm_t
Vlasov_getVolumeStreamTerm(int cdim, int vdim, int polyOrder, int basisType) {
  Vlasov_volumeStreamTerm_t func;
  
  if (cdim == 1) {
    if (vdim == 1) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeStreamTerm_1x1vSerP1, sizeof(Vlasov_volumeStreamTerm_t));
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeStreamTerm_1x1vSerP2, sizeof(Vlasov_volumeStreamTerm_t));
      }
    }
    else if (vdim == 2) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeStreamTerm_1x2vSerP1, sizeof(Vlasov_volumeStreamTerm_t));
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeStreamTerm_1x2vSerP2, sizeof(Vlasov_volumeStreamTerm_t));        
      }      
    }
    else if (vdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeStreamTerm_1x3vSerP1, sizeof(Vlasov_volumeStreamTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeStreamTerm_1x3vSerP2, sizeof(Vlasov_volumeStreamTerm_t));        
      }
    }
  }
  else if (cdim == 2) {
    if (vdim == 2) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeStreamTerm_2x2vSerP1, sizeof(Vlasov_volumeStreamTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeStreamTerm_2x2vSerP2, sizeof(Vlasov_volumeStreamTerm_t));        
      }      
    }
    else if (vdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeStreamTerm_2x3vSerP1, sizeof(Vlasov_volumeStreamTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeStreamTerm_2x3vSerP2, sizeof(Vlasov_volumeStreamTerm_t));
      }
    }    
  }
  else if (cdim == 3) {
    if (vdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeStreamTerm_3x3vSerP1, sizeof(Vlasov_volumeStreamTerm_t));
      }
    }
  }
  return func;
}

Vlasov_surfStreamTerm_t
Vlasov_getSurfStreamTerm(int cdim, int vdim, int polyOrder, int basisType) {
  Vlasov_surfStreamTerm_t func;
  
  if (cdim == 1) {
    if (vdim == 1) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfStreamTerm_1x1vSerP1, sizeof(Vlasov_surfStreamTerm_t));
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfStreamTerm_1x1vSerP2, sizeof(Vlasov_surfStreamTerm_t));
      }
    }
    else if (vdim == 2) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfStreamTerm_1x2vSerP1, sizeof(Vlasov_surfStreamTerm_t));
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfStreamTerm_1x2vSerP2, sizeof(Vlasov_surfStreamTerm_t));        
      }      
    }
    else if (vdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfStreamTerm_1x3vSerP1, sizeof(Vlasov_surfStreamTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfStreamTerm_1x3vSerP2, sizeof(Vlasov_surfStreamTerm_t));        
      }
    }
  }
  else if (cdim == 2) {
    if (vdim == 2) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfStreamTerm_2x2vSerP1, sizeof(Vlasov_surfStreamTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfStreamTerm_2x2vSerP2, sizeof(Vlasov_surfStreamTerm_t));        
      }      
    }
    else if (vdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfStreamTerm_2x3vSerP1, sizeof(Vlasov_surfStreamTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfStreamTerm_2x3vSerP2, sizeof(Vlasov_surfStreamTerm_t));
      }
    }    
  }
  else if (cdim == 3) {
    if (vdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfStreamTerm_3x3vSerP1, sizeof(Vlasov_surfStreamTerm_t));
      }
    }
  }
  return func;  
}

Vlasov_volumeTerm_t
Vlasov_getVolumeTerm(int cdim, int vdim, int polyOrder, int basisType) {
  Vlasov_volumeTerm_t func;
  
  if (cdim == 1) {
    if (vdim == 1) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeTerm_1x1vSerP1, sizeof(Vlasov_volumeTerm_t));
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeTerm_1x1vSerP2, sizeof(Vlasov_volumeTerm_t));
      }
    }
    else if (vdim == 2) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeTerm_1x2vSerP1, sizeof(Vlasov_volumeTerm_t));
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeTerm_1x2vSerP2, sizeof(Vlasov_volumeTerm_t));        
      }      
    }
    else if (vdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeTerm_1x3vSerP1, sizeof(Vlasov_volumeTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeTerm_1x3vSerP2, sizeof(Vlasov_volumeTerm_t));        
      }
    }
  }
  else if (cdim == 2) {
    if (vdim == 2) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeTerm_2x2vSerP1, sizeof(Vlasov_volumeTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeTerm_2x2vSerP2, sizeof(Vlasov_volumeTerm_t));        
      }      
    }
    else if (vdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeTerm_2x3vSerP1, sizeof(Vlasov_volumeTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeTerm_2x3vSerP2, sizeof(Vlasov_volumeTerm_t));
      }
    }    
  }
  else if (cdim == 3) {
    if (vdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_volumeTerm_3x3vSerP1, sizeof(Vlasov_volumeTerm_t));
      }
    }
  }
  return func;
}

Vlasov_surfElcMagTerm_t
Vlasov_getSurfElcMagTerm(int cdim, int vdim, int polyOrder, int basisType) {
  Vlasov_surfElcMagTerm_t func;
  
  if (cdim == 1) {
    if (vdim == 1) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfElcMagTerm_1x1vSerP1, sizeof(Vlasov_surfElcMagTerm_t));
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfElcMagTerm_1x1vSerP2, sizeof(Vlasov_surfElcMagTerm_t));
      }
    }
    else if (vdim == 2) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfElcMagTerm_1x2vSerP1, sizeof(Vlasov_surfElcMagTerm_t));
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfElcMagTerm_1x2vSerP2, sizeof(Vlasov_surfElcMagTerm_t));        
      }      
    }
    else if (vdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfElcMagTerm_1x3vSerP1, sizeof(Vlasov_surfElcMagTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfElcMagTerm_1x3vSerP2, sizeof(Vlasov_surfElcMagTerm_t));        
      }
    }
  }
  else if (cdim == 2) {
    if (vdim == 2) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfElcMagTerm_2x2vSerP1, sizeof(Vlasov_surfElcMagTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfElcMagTerm_2x2vSerP2, sizeof(Vlasov_surfElcMagTerm_t));        
      }      
    }
    else if (vdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfElcMagTerm_2x3vSerP1, sizeof(Vlasov_surfElcMagTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfElcMagTerm_2x3vSerP2, sizeof(Vlasov_surfElcMagTerm_t));
      }
    }    
  }
  else if (cdim == 3) {
    if (vdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Vlasov_surfElcMagTerm_3x3vSerP1, sizeof(Vlasov_surfElcMagTerm_t));
      }
    }
  }
  return func;
}
]]


-- write out string
io.write( vlasovHeaderOut .. vlasovFuncCopyTemplateString )
