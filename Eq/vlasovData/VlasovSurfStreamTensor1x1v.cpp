#include <VlasovModDecl.h> 
void VlasovSurfStream1x1vTensor_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[1]; 
  double wxl = wl[1]; 

  double dvxr = dxvr[1]; 
  double wxr = wr[1]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[4]; 

  if (wxr>0) { 
  incr[0] = (1.732050807568877*fl[1]+fl[0])*wxl+(0.5*fl[3]+0.2886751345948129*fl[2])*dvxl; 
  incr[1] = ((-3.0*fl[1])-1.732050807568877*fl[0])*wxl+((-0.8660254037844386*fl[3])-0.5*fl[2])*dvxl; 
  incr[2] = (1.732050807568877*fl[3]+fl[2])*wxl+(0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[3] = ((-3.0*fl[3])-1.732050807568877*fl[2])*wxl+((-0.8660254037844386*fl[1])-0.5*fl[0])*dvxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += incr[3]*dxl; 
  } else { 
  incr[0] = (fr[0]-1.732050807568877*fr[1])*wxr+(0.2886751345948129*fr[2]-0.5*fr[3])*dvxr; 
  incr[1] = (3.0*fr[1]-1.732050807568877*fr[0])*wxr+(0.8660254037844386*fr[3]-0.5*fr[2])*dvxr; 
  incr[2] = (fr[2]-1.732050807568877*fr[3])*wxr+(0.2886751345948129*fr[0]-0.5*fr[1])*dvxr; 
  incr[3] = (3.0*fr[3]-1.732050807568877*fr[2])*wxr+(0.8660254037844386*fr[1]-0.5*fr[0])*dvxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += incr[3]*dxl; 
  } 
} 
void VlasovSurfStream1x1vTensor_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[1]; 
  double wxl = wl[1]; 

  double dvxr = dxvr[1]; 
  double wxr = wr[1]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[9]; 

  if (wxr>0) { 
  incr[0] = (2.23606797749979*fl[4]+1.732050807568877*fl[1]+fl[0])*wxl+(0.6454972243679028*fl[6]+0.5*fl[3]+0.2886751345948129*fl[2])*dvxl; 
  incr[1] = ((-3.872983346207417*fl[4])-3.0*fl[1]-1.732050807568877*fl[0])*wxl+((-1.118033988749895*fl[6])-0.8660254037844386*fl[3]-0.5*fl[2])*dvxl; 
  incr[2] = (2.23606797749979*fl[6]+1.732050807568877*fl[3]+fl[2])*wxl+(0.5773502691896258*fl[8]+0.447213595499958*fl[7]+0.2581988897471612*fl[5]+0.6454972243679029*fl[4]+0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[3] = ((-3.872983346207417*fl[6])-3.0*fl[3]-1.732050807568877*fl[2])*wxl+((-1.0*fl[8])-0.7745966692414834*fl[7]-0.4472135954999579*fl[5]-1.118033988749895*fl[4]-0.8660254037844386*fl[1]-0.5*fl[0])*dvxl; 
  incr[4] = (5.0*fl[4]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*wxl+(1.443375672974065*fl[6]+1.118033988749895*fl[3]+0.6454972243679029*fl[2])*dvxl; 
  incr[5] = (2.23606797749979*fl[8]+1.732050807568877*fl[7]+fl[5])*wxl+(0.5773502691896257*fl[6]+0.4472135954999579*fl[3]+0.2581988897471612*fl[2])*dvxl; 
  incr[6] = (5.0*fl[6]+3.872983346207417*fl[3]+2.23606797749979*fl[2])*wxl+(1.290994448735806*fl[8]+fl[7]+0.5773502691896257*fl[5]+1.443375672974065*fl[4]+1.118033988749895*fl[1]+0.6454972243679028*fl[0])*dvxl; 
  incr[7] = ((-3.872983346207417*fl[8])-3.0*fl[7]-1.732050807568877*fl[5])*wxl+((-1.0*fl[6])-0.7745966692414834*fl[3]-0.447213595499958*fl[2])*dvxl; 
  incr[8] = (5.0*fl[8]+3.872983346207417*fl[7]+2.23606797749979*fl[5])*wxl+(1.290994448735806*fl[6]+fl[3]+0.5773502691896258*fl[2])*dvxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += incr[7]*dxl; 
  outl[8] += -1.0*incr[8]*dxl; 
  } else { 
  incr[0] = (2.23606797749979*fr[4]-1.732050807568877*fr[1]+fr[0])*wxr+(0.6454972243679028*fr[6]-0.5*fr[3]+0.2886751345948129*fr[2])*dvxr; 
  incr[1] = ((-3.872983346207417*fr[4])+3.0*fr[1]-1.732050807568877*fr[0])*wxr+((-1.118033988749895*fr[6])+0.8660254037844386*fr[3]-0.5*fr[2])*dvxr; 
  incr[2] = (2.23606797749979*fr[6]-1.732050807568877*fr[3]+fr[2])*wxr+(0.5773502691896258*fr[8]-0.447213595499958*fr[7]+0.2581988897471612*fr[5]+0.6454972243679029*fr[4]-0.5*fr[1]+0.2886751345948129*fr[0])*dvxr; 
  incr[3] = ((-3.872983346207417*fr[6])+3.0*fr[3]-1.732050807568877*fr[2])*wxr+((-1.0*fr[8])+0.7745966692414834*fr[7]-0.4472135954999579*fr[5]-1.118033988749895*fr[4]+0.8660254037844386*fr[1]-0.5*fr[0])*dvxr; 
  incr[4] = (5.0*fr[4]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*wxr+(1.443375672974065*fr[6]-1.118033988749895*fr[3]+0.6454972243679029*fr[2])*dvxr; 
  incr[5] = (2.23606797749979*fr[8]-1.732050807568877*fr[7]+fr[5])*wxr+(0.5773502691896257*fr[6]-0.4472135954999579*fr[3]+0.2581988897471612*fr[2])*dvxr; 
  incr[6] = (5.0*fr[6]-3.872983346207417*fr[3]+2.23606797749979*fr[2])*wxr+(1.290994448735806*fr[8]-1.0*fr[7]+0.5773502691896257*fr[5]+1.443375672974065*fr[4]-1.118033988749895*fr[1]+0.6454972243679028*fr[0])*dvxr; 
  incr[7] = ((-3.872983346207417*fr[8])+3.0*fr[7]-1.732050807568877*fr[5])*wxr+((-1.0*fr[6])+0.7745966692414834*fr[3]-0.447213595499958*fr[2])*dvxr; 
  incr[8] = (5.0*fr[8]-3.872983346207417*fr[7]+2.23606797749979*fr[5])*wxr+(1.290994448735806*fr[6]-1.0*fr[3]+0.5773502691896258*fr[2])*dvxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += incr[7]*dxl; 
  outl[8] += -1.0*incr[8]*dxl; 
  } 
} 
void VlasovSurfStream1x1vTensor_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[1]; 
  double wxl = wl[1]; 

  double dvxr = dxvr[1]; 
  double wxr = wr[1]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[16]; 

  if (wxr>0) { 
  incr[0] = (2.645751311064591*fl[8]+2.23606797749979*fl[4]+1.732050807568877*fl[1]+fl[0])*wxl+(0.7637626158259733*fl[11]+0.6454972243679028*fl[6]+0.5*fl[3]+0.2886751345948129*fl[2])*dvxl; 
  incr[1] = ((-4.58257569495584*fl[8])-3.872983346207417*fl[4]-3.0*fl[1]-1.732050807568877*fl[0])*wxl+((-1.322875655532295*fl[11])-1.118033988749895*fl[6]-0.8660254037844386*fl[3]-0.5*fl[2])*dvxl; 
  incr[2] = (2.645751311064591*fl[11]+2.23606797749979*fl[6]+1.732050807568877*fl[3]+fl[2])*wxl+(0.6831300510639734*fl[13]+0.5773502691896258*fl[10]+0.7637626158259735*fl[8]+0.447213595499958*fl[7]+0.2581988897471612*fl[5]+0.6454972243679029*fl[4]+0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[3] = ((-4.58257569495584*fl[11])-3.872983346207417*fl[6]-3.0*fl[3]-1.732050807568877*fl[2])*wxl+((-1.183215956619923*fl[13])-1.0*fl[10]-1.322875655532295*fl[8]-0.7745966692414834*fl[7]-0.4472135954999579*fl[5]-1.118033988749895*fl[4]-0.8660254037844386*fl[1]-0.5*fl[0])*dvxl; 
  incr[4] = (5.916079783099617*fl[8]+5.0*fl[4]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*wxl+(1.707825127659933*fl[11]+1.443375672974065*fl[6]+1.118033988749895*fl[3]+0.6454972243679029*fl[2])*dvxl; 
  incr[5] = (2.645751311064591*fl[13]+2.23606797749979*fl[10]+1.732050807568877*fl[7]+fl[5])*wxl+(0.6708203932499369*fl[15]+0.5669467095138409*fl[14]+0.4391550328268399*fl[12]+0.6831300510639731*fl[11]+0.253546276418555*fl[9]+0.5773502691896257*fl[6]+0.4472135954999579*fl[3]+0.2581988897471612*fl[2])*dvxl; 
  incr[6] = (5.916079783099615*fl[11]+5.0*fl[6]+3.872983346207417*fl[3]+2.23606797749979*fl[2])*wxl+(1.527525231651947*fl[13]+1.290994448735806*fl[10]+1.707825127659933*fl[8]+fl[7]+0.5773502691896257*fl[5]+1.443375672974065*fl[4]+1.118033988749895*fl[1]+0.6454972243679028*fl[0])*dvxl; 
  incr[7] = ((-4.58257569495584*fl[13])-3.872983346207417*fl[10]-3.0*fl[7]-1.732050807568877*fl[5])*wxl+((-1.161895003862225*fl[15])-0.9819805060619657*fl[14]-0.760638829255665*fl[12]-1.183215956619923*fl[11]-0.4391550328268399*fl[9]-1.0*fl[6]-0.7745966692414834*fl[3]-0.447213595499958*fl[2])*dvxl; 
  incr[8] = ((-7.0*fl[8])-5.916079783099617*fl[4]-4.58257569495584*fl[1]-2.645751311064591*fl[0])*wxl+((-2.02072594216369*fl[11])-1.707825127659933*fl[6]-1.322875655532295*fl[3]-0.7637626158259735*fl[2])*dvxl; 
  incr[9] = (2.645751311064591*fl[15]+2.236067977499789*fl[14]+1.732050807568877*fl[12]+fl[9])*wxl+(0.6708203932499368*fl[13]+0.5669467095138407*fl[10]+0.4391550328268399*fl[7]+0.253546276418555*fl[5])*dvxl; 
  incr[10] = (5.916079783099616*fl[13]+5.0*fl[10]+3.872983346207417*fl[7]+2.23606797749979*fl[5])*wxl+(1.5*fl[15]+1.267731382092775*fl[14]+0.9819805060619656*fl[12]+1.527525231651947*fl[11]+0.5669467095138407*fl[9]+1.290994448735806*fl[6]+fl[3]+0.5773502691896258*fl[2])*dvxl; 
  incr[11] = ((-7.0*fl[11])-5.916079783099615*fl[6]-4.58257569495584*fl[3]-2.645751311064591*fl[2])*wxl+((-1.807392228230128*fl[13])-1.527525231651947*fl[10]-2.02072594216369*fl[8]-1.183215956619923*fl[7]-0.6831300510639731*fl[5]-1.707825127659933*fl[4]-1.322875655532295*fl[1]-0.7637626158259733*fl[0])*dvxl; 
  incr[12] = ((-4.58257569495584*fl[15])-3.872983346207417*fl[14]-3.0*fl[12]-1.732050807568877*fl[9])*wxl+((-1.161895003862225*fl[13])-0.9819805060619656*fl[10]-0.760638829255665*fl[7]-0.4391550328268399*fl[5])*dvxl; 
  incr[13] = ((-7.0*fl[13])-5.916079783099616*fl[10]-4.58257569495584*fl[7]-2.645751311064591*fl[5])*wxl+((-1.774823934929885*fl[15])-1.5*fl[14]-1.161895003862225*fl[12]-1.807392228230128*fl[11]-0.6708203932499368*fl[9]-1.527525231651947*fl[6]-1.183215956619923*fl[3]-0.6831300510639734*fl[2])*dvxl; 
  incr[14] = (5.916079783099616*fl[15]+5.0*fl[14]+3.872983346207417*fl[12]+2.236067977499789*fl[9])*wxl+(1.5*fl[13]+1.267731382092775*fl[10]+0.9819805060619657*fl[7]+0.5669467095138409*fl[5])*dvxl; 
  incr[15] = ((-7.0*fl[15])-5.916079783099616*fl[14]-4.58257569495584*fl[12]-2.645751311064591*fl[9])*wxl+((-1.774823934929885*fl[13])-1.5*fl[10]-1.161895003862225*fl[7]-0.6708203932499369*fl[5])*dvxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 
  outr[15] += incr[15]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += incr[7]*dxl; 
  outl[8] += incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += incr[11]*dxl; 
  outl[12] += incr[12]*dxl; 
  outl[13] += incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += incr[15]*dxl; 
  } else { 
  incr[0] = ((-2.645751311064591*fr[8])+2.23606797749979*fr[4]-1.732050807568877*fr[1]+fr[0])*wxr+((-0.7637626158259733*fr[11])+0.6454972243679028*fr[6]-0.5*fr[3]+0.2886751345948129*fr[2])*dvxr; 
  incr[1] = (4.58257569495584*fr[8]-3.872983346207417*fr[4]+3.0*fr[1]-1.732050807568877*fr[0])*wxr+(1.322875655532295*fr[11]-1.118033988749895*fr[6]+0.8660254037844386*fr[3]-0.5*fr[2])*dvxr; 
  incr[2] = ((-2.645751311064591*fr[11])+2.23606797749979*fr[6]-1.732050807568877*fr[3]+fr[2])*wxr+((-0.6831300510639734*fr[13])+0.5773502691896258*fr[10]-0.7637626158259735*fr[8]-0.447213595499958*fr[7]+0.2581988897471612*fr[5]+0.6454972243679029*fr[4]-0.5*fr[1]+0.2886751345948129*fr[0])*dvxr; 
  incr[3] = (4.58257569495584*fr[11]-3.872983346207417*fr[6]+3.0*fr[3]-1.732050807568877*fr[2])*wxr+(1.183215956619923*fr[13]-1.0*fr[10]+1.322875655532295*fr[8]+0.7745966692414834*fr[7]-0.4472135954999579*fr[5]-1.118033988749895*fr[4]+0.8660254037844386*fr[1]-0.5*fr[0])*dvxr; 
  incr[4] = ((-5.916079783099617*fr[8])+5.0*fr[4]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*wxr+((-1.707825127659933*fr[11])+1.443375672974065*fr[6]-1.118033988749895*fr[3]+0.6454972243679029*fr[2])*dvxr; 
  incr[5] = ((-2.645751311064591*fr[13])+2.23606797749979*fr[10]-1.732050807568877*fr[7]+fr[5])*wxr+((-0.6708203932499369*fr[15])+0.5669467095138409*fr[14]-0.4391550328268399*fr[12]-0.6831300510639731*fr[11]+0.253546276418555*fr[9]+0.5773502691896257*fr[6]-0.4472135954999579*fr[3]+0.2581988897471612*fr[2])*dvxr; 
  incr[6] = ((-5.916079783099615*fr[11])+5.0*fr[6]-3.872983346207417*fr[3]+2.23606797749979*fr[2])*wxr+((-1.527525231651947*fr[13])+1.290994448735806*fr[10]-1.707825127659933*fr[8]-1.0*fr[7]+0.5773502691896257*fr[5]+1.443375672974065*fr[4]-1.118033988749895*fr[1]+0.6454972243679028*fr[0])*dvxr; 
  incr[7] = (4.58257569495584*fr[13]-3.872983346207417*fr[10]+3.0*fr[7]-1.732050807568877*fr[5])*wxr+(1.161895003862225*fr[15]-0.9819805060619657*fr[14]+0.760638829255665*fr[12]+1.183215956619923*fr[11]-0.4391550328268399*fr[9]-1.0*fr[6]+0.7745966692414834*fr[3]-0.447213595499958*fr[2])*dvxr; 
  incr[8] = (7.0*fr[8]-5.916079783099617*fr[4]+4.58257569495584*fr[1]-2.645751311064591*fr[0])*wxr+(2.02072594216369*fr[11]-1.707825127659933*fr[6]+1.322875655532295*fr[3]-0.7637626158259735*fr[2])*dvxr; 
  incr[9] = ((-2.645751311064591*fr[15])+2.236067977499789*fr[14]-1.732050807568877*fr[12]+fr[9])*wxr+((-0.6708203932499368*fr[13])+0.5669467095138407*fr[10]-0.4391550328268399*fr[7]+0.253546276418555*fr[5])*dvxr; 
  incr[10] = ((-5.916079783099616*fr[13])+5.0*fr[10]-3.872983346207417*fr[7]+2.23606797749979*fr[5])*wxr+((-1.5*fr[15])+1.267731382092775*fr[14]-0.9819805060619656*fr[12]-1.527525231651947*fr[11]+0.5669467095138407*fr[9]+1.290994448735806*fr[6]-1.0*fr[3]+0.5773502691896258*fr[2])*dvxr; 
  incr[11] = (7.0*fr[11]-5.916079783099615*fr[6]+4.58257569495584*fr[3]-2.645751311064591*fr[2])*wxr+(1.807392228230128*fr[13]-1.527525231651947*fr[10]+2.02072594216369*fr[8]+1.183215956619923*fr[7]-0.6831300510639731*fr[5]-1.707825127659933*fr[4]+1.322875655532295*fr[1]-0.7637626158259733*fr[0])*dvxr; 
  incr[12] = (4.58257569495584*fr[15]-3.872983346207417*fr[14]+3.0*fr[12]-1.732050807568877*fr[9])*wxr+(1.161895003862225*fr[13]-0.9819805060619656*fr[10]+0.760638829255665*fr[7]-0.4391550328268399*fr[5])*dvxr; 
  incr[13] = (7.0*fr[13]-5.916079783099616*fr[10]+4.58257569495584*fr[7]-2.645751311064591*fr[5])*wxr+(1.774823934929885*fr[15]-1.5*fr[14]+1.161895003862225*fr[12]+1.807392228230128*fr[11]-0.6708203932499368*fr[9]-1.527525231651947*fr[6]+1.183215956619923*fr[3]-0.6831300510639734*fr[2])*dvxr; 
  incr[14] = ((-5.916079783099616*fr[15])+5.0*fr[14]-3.872983346207417*fr[12]+2.236067977499789*fr[9])*wxr+((-1.5*fr[13])+1.267731382092775*fr[10]-0.9819805060619657*fr[7]+0.5669467095138409*fr[5])*dvxr; 
  incr[15] = (7.0*fr[15]-5.916079783099616*fr[14]+4.58257569495584*fr[12]-2.645751311064591*fr[9])*wxr+(1.774823934929885*fr[13]-1.5*fr[10]+1.161895003862225*fr[7]-0.6708203932499369*fr[5])*dvxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 
  outr[15] += incr[15]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += incr[7]*dxl; 
  outl[8] += incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += incr[11]*dxl; 
  outl[12] += incr[12]*dxl; 
  outl[13] += incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += incr[15]*dxl; 
  } 
} 
