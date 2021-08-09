// Gkyl ------------------------------------------------------------------------
//
// Code to copy various kernerls from device to host
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <VlasovTmplModDecl.h>

using namespace Gkyl;



__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x1vSerP1 =
  &Gkyl::VlasovModDecl<1,1,1,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfStreamTerm_t p_Vlasov_surfStreamTerm_1x1vSerP1 =
  &Gkyl::VlasovModDecl<1,1,1,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x1vSerP1 =
  &Gkyl::VlasovModDecl<1,1,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x1vSerP1 =
  &Gkyl::VlasovModDecl<1,1,1,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x1vSerP2 =
  &Gkyl::VlasovModDecl<1,1,2,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfStreamTerm_t p_Vlasov_surfStreamTerm_1x1vSerP2 =
  &Gkyl::VlasovModDecl<1,1,2,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x1vSerP2 =
  &Gkyl::VlasovModDecl<1,1,2,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x1vSerP2 =
  &Gkyl::VlasovModDecl<1,1,2,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x2vSerP1 =
  &Gkyl::VlasovModDecl<1,2,1,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfStreamTerm_t p_Vlasov_surfStreamTerm_1x2vSerP1 =
  &Gkyl::VlasovModDecl<1,2,1,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x2vSerP1 =
  &Gkyl::VlasovModDecl<1,2,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x2vSerP1 =
  &Gkyl::VlasovModDecl<1,2,1,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x2vSerP2 =
  &Gkyl::VlasovModDecl<1,2,2,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfStreamTerm_t p_Vlasov_surfStreamTerm_1x2vSerP2 =
  &Gkyl::VlasovModDecl<1,2,2,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x2vSerP2 =
  &Gkyl::VlasovModDecl<1,2,2,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x2vSerP2 =
  &Gkyl::VlasovModDecl<1,2,2,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x3vSerP1 =
  &Gkyl::VlasovModDecl<1,3,1,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfStreamTerm_t p_Vlasov_surfStreamTerm_1x3vSerP1 =
  &Gkyl::VlasovModDecl<1,3,1,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x3vSerP1 =
  &Gkyl::VlasovModDecl<1,3,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x3vSerP1 =
  &Gkyl::VlasovModDecl<1,3,1,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_1x3vSerP2 =
  &Gkyl::VlasovModDecl<1,3,2,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfStreamTerm_t p_Vlasov_surfStreamTerm_1x3vSerP2 =
  &Gkyl::VlasovModDecl<1,3,2,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_1x3vSerP2 =
  &Gkyl::VlasovModDecl<1,3,2,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_1x3vSerP2 =
  &Gkyl::VlasovModDecl<1,3,2,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_2x2vSerP1 =
  &Gkyl::VlasovModDecl<2,2,1,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfStreamTerm_t p_Vlasov_surfStreamTerm_2x2vSerP1 =
  &Gkyl::VlasovModDecl<2,2,1,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_2x2vSerP1 =
  &Gkyl::VlasovModDecl<2,2,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_2x2vSerP1 =
  &Gkyl::VlasovModDecl<2,2,1,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_2x2vSerP2 =
  &Gkyl::VlasovModDecl<2,2,2,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfStreamTerm_t p_Vlasov_surfStreamTerm_2x2vSerP2 =
  &Gkyl::VlasovModDecl<2,2,2,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_2x2vSerP2 =
  &Gkyl::VlasovModDecl<2,2,2,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_2x2vSerP2 =
  &Gkyl::VlasovModDecl<2,2,2,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_2x3vSerP1 =
  &Gkyl::VlasovModDecl<2,3,1,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfStreamTerm_t p_Vlasov_surfStreamTerm_2x3vSerP1 =
  &Gkyl::VlasovModDecl<2,3,1,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_2x3vSerP1 =
  &Gkyl::VlasovModDecl<2,3,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_2x3vSerP1 =
  &Gkyl::VlasovModDecl<2,3,1,G_SERENDIPITY_C>::surfElcMagTerm;


__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_2x3vSerP2 =
  &Gkyl::VlasovModDecl<2,3,2,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfStreamTerm_t p_Vlasov_surfStreamTerm_2x3vSerP2 =
  &Gkyl::VlasovModDecl<2,3,2,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_2x3vSerP2 =
  &Gkyl::VlasovModDecl<2,3,2,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_2x3vSerP2 =
  &Gkyl::VlasovModDecl<2,3,2,G_SERENDIPITY_C>::surfElcMagTerm;




__device__ Vlasov_volumeStreamTerm_t p_Vlasov_volumeStreamTerm_3x3vSerP1 =
  &Gkyl::VlasovModDecl<3,3,1,G_SERENDIPITY_C>::volumeStreamTerm;

__device__ Vlasov_surfStreamTerm_t p_Vlasov_surfStreamTerm_3x3vSerP1 =
  &Gkyl::VlasovModDecl<3,3,1,G_SERENDIPITY_C>::surfStreamTerm;;

__device__ Vlasov_volumeTerm_t p_Vlasov_volumeTerm_3x3vSerP1 =
  &Gkyl::VlasovModDecl<3,3,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Vlasov_surfElcMagTerm_t p_Vlasov_surfElcMagTerm_3x3vSerP1 =
  &Gkyl::VlasovModDecl<3,3,1,G_SERENDIPITY_C>::surfElcMagTerm;


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
