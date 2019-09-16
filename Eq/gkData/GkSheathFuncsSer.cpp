#include <cmath> 
#include <fenv.h> 
#include <GyrokineticModDecl.h> 
// approximation for inverse Langevin function 
long double invL(long double x) { 
  // from Kroger 
  return (3.*x-x*x*x*(6. + x*x - 2.*x*x*x*x)/5.)/(1.-x*x); 
}

void calcSheathReflection3x2vSer_P1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq_i; long double xc, b, rVal, fac; 
  double fReflXYQuad[4][8]; 
  double fReflXYMuQuad[2][4]; 
  feclearexcept(FE_ALL_EXCEPT); 
  

// quadrature node (x,y)_i=1 
  vcutSq_i = (0.1666666666666667*q_*((2.449489742783178*phiWall[7]-2.449489742783178*phi[7]-4.242640687119286*phiWall[6]+4.242640687119286*phi[6]-4.242640687119286*phiWall[5]+4.242640687119286*phi[5]+7.348469228349534*phiWall[3]-7.348469228349534*phi[3])*zVal+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]-2.449489742783178*phiWall[2]+2.449489742783178*phi[2]-2.449489742783178*phiWall[1]+2.449489742783178*phi[1]+4.242640687119286*phiWall[0]-4.242640687119286*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[0][0] = 0.0; 
  fReflXYQuad[0][1] = 0.0; 
  fReflXYQuad[0][2] = 0.0; 
  fReflXYQuad[0][3] = 0.0; 
  fReflXYQuad[0][4] = 0.0; 
  fReflXYQuad[0][5] = 0.0; 
  fReflXYQuad[0][6] = 0.0; 
  fReflXYQuad[0][7] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[0][0] = 0.1666666666666667*(f[6]-1.732050807568877*(f[2]+f[1])+3.0*f[0]); 
  fReflXYQuad[0][1] = 0.1666666666666667*(f[16]-1.732050807568877*(f[8]+f[7])+3.0*f[3]); 
  fReflXYQuad[0][2] = 0.1666666666666667*(f[17]-1.732050807568877*(f[10]+f[9])+3.0*f[4]); 
  fReflXYQuad[0][3] = 0.1666666666666667*(f[20]-1.732050807568877*(f[13]+f[12])+3.0*f[5]); 
  fReflXYQuad[0][4] = 0.1666666666666667*(f[26]-1.732050807568877*(f[19]+f[18])+3.0*f[11]); 
  fReflXYQuad[0][5] = 0.1666666666666667*(f[27]-1.732050807568877*(f[22]+f[21])+3.0*f[14]); 
  fReflXYQuad[0][6] = 0.1666666666666667*(f[28]-1.732050807568877*(f[24]+f[23])+3.0*f[15]); 
  fReflXYQuad[0][7] = 0.1666666666666667*(f[31]-1.732050807568877*(f[30]+f[29])+3.0*f[25]); 
  } else { // partial reflection 
  rVal = (0.5773502691896258*(1.732050807568877*f[27]-3.0*(f[22]+f[21]+f[16])+5.196152422706631*(f[14]+f[8]+f[7])-9.0*f[3]))/(1.732050807568877*f[20]-3.0*(f[13]+f[12]+f[6])+5.196152422706631*(f[5]+f[2]+f[1])-9.0*f[0]); 
  // if f is not realizable, no reflection 
  if(std::abs(rVal)>=.99 || -0.0392837100659193*(1.732050807568877*f[20]-3.0*(f[13]+f[12]+f[6])+5.196152422706631*(f[5]+f[2]+f[1])-9.0*f[0]) < 0.) { 
  fReflXYMuQuad[0][0] = 0.0; 
  fReflXYMuQuad[0][1] = 0.0; 
  fReflXYMuQuad[0][2] = 0.0; 
  fReflXYMuQuad[0][3] = 0.0; 
  } else {
   b = invL(rVal); 
   if(wv > 0) {
    xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
    fac = b>500? 0. : b<-500? 1. : abs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,1,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][0] = (-0.0392837100659193*(1.732050807568877*f[20]-3.0*(f[13]+f[12]+f[6])+5.196152422706631*(f[5]+f[2]+f[1])-9.0*f[0]))*fac; 
    fac = (b>500 || abs(b)<1e-10)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,2,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][1] = (-0.0392837100659193*(1.732050807568877*f[27]-3.0*(f[22]+f[21]+f[16])+5.196152422706631*(f[14]+f[8]+f[7])-9.0*f[3]))*fac; 
   } else { 
    xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
    fac = b>500? 1. : b<-500? 0. : abs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,1,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][0] = (-0.0392837100659193*(1.732050807568877*f[20]-3.0*(f[13]+f[12]+f[6])+5.196152422706631*(f[5]+f[2]+f[1])-9.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || abs(b)<1e-10)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,2,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][1] = (-0.0392837100659193*(1.732050807568877*f[27]-3.0*(f[22]+f[21]+f[16])+5.196152422706631*(f[14]+f[8]+f[7])-9.0*f[3]))*fac; 
   } 
  } 
  rVal = (0.5773502691896258*(1.732050807568877*f[27]+3.0*(f[16]-1.0*(f[22]+f[21]))+5.196152422706631*(f[14]-1.0*(f[8]+f[7]))+9.0*f[3]))/(1.732050807568877*f[20]+3.0*(f[6]-1.0*(f[13]+f[12]))+5.196152422706631*(f[5]-1.0*(f[2]+f[1]))+9.0*f[0]); 
  // if f is not realizable, no reflection 
  if(std::abs(rVal)>=.99 || 0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[6]-1.0*(f[13]+f[12]))+5.196152422706631*(f[5]-1.0*(f[2]+f[1]))+9.0*f[0]) < 0.) { 
  fReflXYMuQuad[1][0] = 0.0; 
  fReflXYMuQuad[1][1] = 0.0; 
  fReflXYMuQuad[1][2] = 0.0; 
  fReflXYMuQuad[1][3] = 0.0; 
  } else {
   b = invL(rVal); 
   if(wv > 0) {
    xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
    fac = b>500? 0. : b<-500? 1. : abs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,1,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][0] = (0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[6]-1.0*(f[13]+f[12]))+5.196152422706631*(f[5]-1.0*(f[2]+f[1]))+9.0*f[0]))*fac; 
    fac = (b>500 || abs(b)<1e-10)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,2,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][1] = (0.0392837100659193*(1.732050807568877*f[27]+3.0*(f[16]-1.0*(f[22]+f[21]))+5.196152422706631*(f[14]-1.0*(f[8]+f[7]))+9.0*f[3]))*fac; 
   } else { 
    xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
    fac = b>500? 1. : b<-500? 0. : abs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,1,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][0] = (0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[6]-1.0*(f[13]+f[12]))+5.196152422706631*(f[5]-1.0*(f[2]+f[1]))+9.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || abs(b)<1e-10)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,2,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][1] = (0.0392837100659193*(1.732050807568877*f[27]+3.0*(f[16]-1.0*(f[22]+f[21]))+5.196152422706631*(f[14]-1.0*(f[8]+f[7]))+9.0*f[3]))*fac; 
   } 
  } 
  fReflXYQuad[0][0] = 0.7071067811865468*(fReflXYMuQuad[1][0]+fReflXYMuQuad[0][0]); 
  fReflXYQuad[0][1] = 0.7071067811865468*(fReflXYMuQuad[1][1]+fReflXYMuQuad[0][1]); 
  fReflXYQuad[0][2] = 0.7071067811865468*(fReflXYMuQuad[1][2]+fReflXYMuQuad[0][2]); 
  fReflXYQuad[0][3] = 1.224744871391589*(fReflXYMuQuad[1][0]-1.0*fReflXYMuQuad[0][0]); 
  fReflXYQuad[0][4] = 0.7071067811865468*(fReflXYMuQuad[1][3]+fReflXYMuQuad[0][3]); 
  fReflXYQuad[0][5] = 1.224744871391589*(fReflXYMuQuad[1][1]-1.0*fReflXYMuQuad[0][1]); 
  fReflXYQuad[0][6] = 1.224744871391589*(fReflXYMuQuad[1][2]-1.0*fReflXYMuQuad[0][2]); 
  fReflXYQuad[0][7] = 1.224744871391589*(fReflXYMuQuad[1][3]-1.0*fReflXYMuQuad[0][3]); 
  } 

 
// quadrature node (x,y)_i=2 
  vcutSq_i = -(0.1666666666666667*q_*((2.449489742783178*phiWall[7]-2.449489742783178*phi[7]+4.242640687119286*phiWall[6]-4.242640687119286*(phi[6]+phiWall[5])+4.242640687119286*phi[5]-7.348469228349534*phiWall[3]+7.348469228349534*phi[3])*zVal+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]+2.449489742783178*phiWall[2]-2.449489742783178*(phi[2]+phiWall[1])+2.449489742783178*phi[1]-4.242640687119286*phiWall[0]+4.242640687119286*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[1][0] = 0.0; 
  fReflXYQuad[1][1] = 0.0; 
  fReflXYQuad[1][2] = 0.0; 
  fReflXYQuad[1][3] = 0.0; 
  fReflXYQuad[1][4] = 0.0; 
  fReflXYQuad[1][5] = 0.0; 
  fReflXYQuad[1][6] = 0.0; 
  fReflXYQuad[1][7] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[1][0] = -0.1666666666666667*(f[6]+1.732050807568877*f[2]-1.0*(1.732050807568877*f[1]+3.0*f[0])); 
  fReflXYQuad[1][1] = -0.1666666666666667*(f[16]+1.732050807568877*f[8]-1.0*(1.732050807568877*f[7]+3.0*f[3])); 
  fReflXYQuad[1][2] = -0.1666666666666667*(f[17]+1.732050807568877*f[10]-1.0*(1.732050807568877*f[9]+3.0*f[4])); 
  fReflXYQuad[1][3] = -0.1666666666666667*(f[20]+1.732050807568877*f[13]-1.0*(1.732050807568877*f[12]+3.0*f[5])); 
  fReflXYQuad[1][4] = -0.1666666666666667*(f[26]+1.732050807568877*f[19]-1.0*(1.732050807568877*f[18]+3.0*f[11])); 
  fReflXYQuad[1][5] = -0.1666666666666667*(f[27]+1.732050807568877*f[22]-1.0*(1.732050807568877*f[21]+3.0*f[14])); 
  fReflXYQuad[1][6] = -0.1666666666666667*(f[28]+1.732050807568877*f[24]-1.0*(1.732050807568877*f[23]+3.0*f[15])); 
  fReflXYQuad[1][7] = -0.1666666666666667*(f[31]+1.732050807568877*f[30]-1.0*(1.732050807568877*f[29]+3.0*f[25])); 
  } else { // partial reflection 
  rVal = (0.5773502691896258*(1.732050807568877*f[27]+3.0*(f[22]-1.0*(f[21]+f[16]))+5.196152422706631*(f[7]-1.0*(f[14]+f[8]))+9.0*f[3]))/(1.732050807568877*f[20]+3.0*(f[13]-1.0*(f[12]+f[6]))+5.196152422706631*(f[1]-1.0*(f[5]+f[2]))+9.0*f[0]); 
  // if f is not realizable, no reflection 
  if(std::abs(rVal)>=.99 || 0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[13]-1.0*(f[12]+f[6]))+5.196152422706631*(f[1]-1.0*(f[5]+f[2]))+9.0*f[0]) < 0.) { 
  fReflXYMuQuad[0][0] = 0.0; 
  fReflXYMuQuad[0][1] = 0.0; 
  fReflXYMuQuad[0][2] = 0.0; 
  fReflXYMuQuad[0][3] = 0.0; 
  } else {
   b = invL(rVal); 
   if(wv > 0) {
    xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
    fac = b>500? 0. : b<-500? 1. : abs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,1,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][0] = (0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[13]-1.0*(f[12]+f[6]))+5.196152422706631*(f[1]-1.0*(f[5]+f[2]))+9.0*f[0]))*fac; 
    fac = (b>500 || abs(b)<1e-10)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,2,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][1] = (0.0392837100659193*(1.732050807568877*f[27]+3.0*(f[22]-1.0*(f[21]+f[16]))+5.196152422706631*(f[7]-1.0*(f[14]+f[8]))+9.0*f[3]))*fac; 
   } else { 
    xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
    fac = b>500? 1. : b<-500? 0. : abs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,1,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][0] = (0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[13]-1.0*(f[12]+f[6]))+5.196152422706631*(f[1]-1.0*(f[5]+f[2]))+9.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || abs(b)<1e-10)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,2,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][1] = (0.0392837100659193*(1.732050807568877*f[27]+3.0*(f[22]-1.0*(f[21]+f[16]))+5.196152422706631*(f[7]-1.0*(f[14]+f[8]))+9.0*f[3]))*fac; 
   } 
  } 
  rVal = (0.5773502691896258*(1.732050807568877*f[27]+3.0*(f[22]-1.0*f[21]+f[16])+5.196152422706631*(f[8]-1.0*f[14])-1.0*(5.196152422706631*f[7]+9.0*f[3])))/(1.732050807568877*f[20]+3.0*(f[13]-1.0*f[12]+f[6])+5.196152422706631*(f[2]-1.0*f[5])-1.0*(5.196152422706631*f[1]+9.0*f[0])); 
  // if f is not realizable, no reflection 
  if(std::abs(rVal)>=.99 || -0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[13]-1.0*f[12]+f[6])+5.196152422706631*(f[2]-1.0*f[5])-1.0*(5.196152422706631*f[1]+9.0*f[0])) < 0.) { 
  fReflXYMuQuad[1][0] = 0.0; 
  fReflXYMuQuad[1][1] = 0.0; 
  fReflXYMuQuad[1][2] = 0.0; 
  fReflXYMuQuad[1][3] = 0.0; 
  } else {
   b = invL(rVal); 
   if(wv > 0) {
    xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
    fac = b>500? 0. : b<-500? 1. : abs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,1,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][0] = (-0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[13]-1.0*f[12]+f[6])+5.196152422706631*(f[2]-1.0*f[5])-1.0*(5.196152422706631*f[1]+9.0*f[0])))*fac; 
    fac = (b>500 || abs(b)<1e-10)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,2,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][1] = (-0.0392837100659193*(1.732050807568877*f[27]+3.0*(f[22]-1.0*f[21]+f[16])+5.196152422706631*(f[8]-1.0*f[14])-1.0*(5.196152422706631*f[7]+9.0*f[3])))*fac; 
   } else { 
    xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
    fac = b>500? 1. : b<-500? 0. : abs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,1,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][0] = (-0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[13]-1.0*f[12]+f[6])+5.196152422706631*(f[2]-1.0*f[5])-1.0*(5.196152422706631*f[1]+9.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || abs(b)<1e-10)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,2,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][1] = (-0.0392837100659193*(1.732050807568877*f[27]+3.0*(f[22]-1.0*f[21]+f[16])+5.196152422706631*(f[8]-1.0*f[14])-1.0*(5.196152422706631*f[7]+9.0*f[3])))*fac; 
   } 
  } 
  fReflXYQuad[1][0] = 0.7071067811865468*(fReflXYMuQuad[1][0]+fReflXYMuQuad[0][0]); 
  fReflXYQuad[1][1] = 0.7071067811865468*(fReflXYMuQuad[1][1]+fReflXYMuQuad[0][1]); 
  fReflXYQuad[1][2] = 0.7071067811865468*(fReflXYMuQuad[1][2]+fReflXYMuQuad[0][2]); 
  fReflXYQuad[1][3] = 1.224744871391589*(fReflXYMuQuad[1][0]-1.0*fReflXYMuQuad[0][0]); 
  fReflXYQuad[1][4] = 0.7071067811865468*(fReflXYMuQuad[1][3]+fReflXYMuQuad[0][3]); 
  fReflXYQuad[1][5] = 1.224744871391589*(fReflXYMuQuad[1][1]-1.0*fReflXYMuQuad[0][1]); 
  fReflXYQuad[1][6] = 1.224744871391589*(fReflXYMuQuad[1][2]-1.0*fReflXYMuQuad[0][2]); 
  fReflXYQuad[1][7] = 1.224744871391589*(fReflXYMuQuad[1][3]-1.0*fReflXYMuQuad[0][3]); 
  } 

 
// quadrature node (x,y)_i=3 
  vcutSq_i = -(0.1666666666666667*q_*((2.449489742783178*phiWall[7]-2.449489742783178*phi[7]-4.242640687119286*phiWall[6]+4.242640687119286*(phi[6]+phiWall[5])-4.242640687119286*phi[5]-7.348469228349534*phiWall[3]+7.348469228349534*phi[3])*zVal+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]-2.449489742783178*phiWall[2]+2.449489742783178*(phi[2]+phiWall[1])-2.449489742783178*phi[1]-4.242640687119286*phiWall[0]+4.242640687119286*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[2][0] = 0.0; 
  fReflXYQuad[2][1] = 0.0; 
  fReflXYQuad[2][2] = 0.0; 
  fReflXYQuad[2][3] = 0.0; 
  fReflXYQuad[2][4] = 0.0; 
  fReflXYQuad[2][5] = 0.0; 
  fReflXYQuad[2][6] = 0.0; 
  fReflXYQuad[2][7] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[2][0] = -0.1666666666666667*(f[6]+1.732050807568877*(f[1]-1.0*f[2])-3.0*f[0]); 
  fReflXYQuad[2][1] = -0.1666666666666667*(f[16]+1.732050807568877*(f[7]-1.0*f[8])-3.0*f[3]); 
  fReflXYQuad[2][2] = -0.1666666666666667*(f[17]+1.732050807568877*(f[9]-1.0*f[10])-3.0*f[4]); 
  fReflXYQuad[2][3] = -0.1666666666666667*(f[20]+1.732050807568877*(f[12]-1.0*f[13])-3.0*f[5]); 
  fReflXYQuad[2][4] = -0.1666666666666667*(f[26]+1.732050807568877*(f[18]-1.0*f[19])-3.0*f[11]); 
  fReflXYQuad[2][5] = -0.1666666666666667*(f[27]+1.732050807568877*(f[21]-1.0*f[22])-3.0*f[14]); 
  fReflXYQuad[2][6] = -0.1666666666666667*(f[28]+1.732050807568877*(f[23]-1.0*f[24])-3.0*f[15]); 
  fReflXYQuad[2][7] = -0.1666666666666667*(f[31]+1.732050807568877*(f[29]-1.0*f[30])-3.0*f[25]); 
  } else { // partial reflection 
  rVal = (0.5773502691896258*(1.732050807568877*f[27]+3.0*((-1.0*f[22])+f[21]-1.0*f[16])+5.196152422706631*((-1.0*f[14])+f[8]-1.0*f[7])+9.0*f[3]))/(1.732050807568877*f[20]+3.0*((-1.0*f[13])+f[12]-1.0*f[6])+5.196152422706631*((-1.0*f[5])+f[2]-1.0*f[1])+9.0*f[0]); 
  // if f is not realizable, no reflection 
  if(std::abs(rVal)>=.99 || 0.0392837100659193*(1.732050807568877*f[20]+3.0*((-1.0*f[13])+f[12]-1.0*f[6])+5.196152422706631*((-1.0*f[5])+f[2]-1.0*f[1])+9.0*f[0]) < 0.) { 
  fReflXYMuQuad[0][0] = 0.0; 
  fReflXYMuQuad[0][1] = 0.0; 
  fReflXYMuQuad[0][2] = 0.0; 
  fReflXYMuQuad[0][3] = 0.0; 
  } else {
   b = invL(rVal); 
   if(wv > 0) {
    xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
    fac = b>500? 0. : b<-500? 1. : abs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,1,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][0] = (0.0392837100659193*(1.732050807568877*f[20]+3.0*((-1.0*f[13])+f[12]-1.0*f[6])+5.196152422706631*((-1.0*f[5])+f[2]-1.0*f[1])+9.0*f[0]))*fac; 
    fac = (b>500 || abs(b)<1e-10)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,2,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][1] = (0.0392837100659193*(1.732050807568877*f[27]+3.0*((-1.0*f[22])+f[21]-1.0*f[16])+5.196152422706631*((-1.0*f[14])+f[8]-1.0*f[7])+9.0*f[3]))*fac; 
   } else { 
    xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
    fac = b>500? 1. : b<-500? 0. : abs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,1,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][0] = (0.0392837100659193*(1.732050807568877*f[20]+3.0*((-1.0*f[13])+f[12]-1.0*f[6])+5.196152422706631*((-1.0*f[5])+f[2]-1.0*f[1])+9.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || abs(b)<1e-10)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,2,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][1] = (0.0392837100659193*(1.732050807568877*f[27]+3.0*((-1.0*f[22])+f[21]-1.0*f[16])+5.196152422706631*((-1.0*f[14])+f[8]-1.0*f[7])+9.0*f[3]))*fac; 
   } 
  } 
  rVal = (0.5773502691896258*(1.732050807568877*f[27]+3.0*((-1.0*f[22])+f[21]+f[16])+5.196152422706631*(f[7]-1.0*(f[14]+f[8]))-9.0*f[3]))/(1.732050807568877*f[20]+3.0*((-1.0*f[13])+f[12]+f[6])+5.196152422706631*(f[1]-1.0*(f[5]+f[2]))-9.0*f[0]); 
  // if f is not realizable, no reflection 
  if(std::abs(rVal)>=.99 || -0.0392837100659193*(1.732050807568877*f[20]+3.0*((-1.0*f[13])+f[12]+f[6])+5.196152422706631*(f[1]-1.0*(f[5]+f[2]))-9.0*f[0]) < 0.) { 
  fReflXYMuQuad[1][0] = 0.0; 
  fReflXYMuQuad[1][1] = 0.0; 
  fReflXYMuQuad[1][2] = 0.0; 
  fReflXYMuQuad[1][3] = 0.0; 
  } else {
   b = invL(rVal); 
   if(wv > 0) {
    xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
    fac = b>500? 0. : b<-500? 1. : abs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,1,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][0] = (-0.0392837100659193*(1.732050807568877*f[20]+3.0*((-1.0*f[13])+f[12]+f[6])+5.196152422706631*(f[1]-1.0*(f[5]+f[2]))-9.0*f[0]))*fac; 
    fac = (b>500 || abs(b)<1e-10)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,2,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][1] = (-0.0392837100659193*(1.732050807568877*f[27]+3.0*((-1.0*f[22])+f[21]+f[16])+5.196152422706631*(f[7]-1.0*(f[14]+f[8]))-9.0*f[3]))*fac; 
   } else { 
    xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
    fac = b>500? 1. : b<-500? 0. : abs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,1,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][0] = (-0.0392837100659193*(1.732050807568877*f[20]+3.0*((-1.0*f[13])+f[12]+f[6])+5.196152422706631*(f[1]-1.0*(f[5]+f[2]))-9.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || abs(b)<1e-10)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,2,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][1] = (-0.0392837100659193*(1.732050807568877*f[27]+3.0*((-1.0*f[22])+f[21]+f[16])+5.196152422706631*(f[7]-1.0*(f[14]+f[8]))-9.0*f[3]))*fac; 
   } 
  } 
  fReflXYQuad[2][0] = 0.7071067811865468*(fReflXYMuQuad[1][0]+fReflXYMuQuad[0][0]); 
  fReflXYQuad[2][1] = 0.7071067811865468*(fReflXYMuQuad[1][1]+fReflXYMuQuad[0][1]); 
  fReflXYQuad[2][2] = 0.7071067811865468*(fReflXYMuQuad[1][2]+fReflXYMuQuad[0][2]); 
  fReflXYQuad[2][3] = 1.224744871391589*(fReflXYMuQuad[1][0]-1.0*fReflXYMuQuad[0][0]); 
  fReflXYQuad[2][4] = 0.7071067811865468*(fReflXYMuQuad[1][3]+fReflXYMuQuad[0][3]); 
  fReflXYQuad[2][5] = 1.224744871391589*(fReflXYMuQuad[1][1]-1.0*fReflXYMuQuad[0][1]); 
  fReflXYQuad[2][6] = 1.224744871391589*(fReflXYMuQuad[1][2]-1.0*fReflXYMuQuad[0][2]); 
  fReflXYQuad[2][7] = 1.224744871391589*(fReflXYMuQuad[1][3]-1.0*fReflXYMuQuad[0][3]); 
  } 

 
// quadrature node (x,y)_i=4 
  vcutSq_i = (0.1666666666666667*q_*((2.449489742783178*phiWall[7]-2.449489742783178*phi[7]+4.242640687119286*phiWall[6]-4.242640687119286*phi[6]+4.242640687119286*phiWall[5]-4.242640687119286*phi[5]+7.348469228349534*phiWall[3]-7.348469228349534*phi[3])*zVal+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]+2.449489742783178*phiWall[2]-2.449489742783178*phi[2]+2.449489742783178*phiWall[1]-2.449489742783178*phi[1]+4.242640687119286*phiWall[0]-4.242640687119286*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[3][0] = 0.0; 
  fReflXYQuad[3][1] = 0.0; 
  fReflXYQuad[3][2] = 0.0; 
  fReflXYQuad[3][3] = 0.0; 
  fReflXYQuad[3][4] = 0.0; 
  fReflXYQuad[3][5] = 0.0; 
  fReflXYQuad[3][6] = 0.0; 
  fReflXYQuad[3][7] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[3][0] = 0.1666666666666667*(f[6]+1.732050807568877*(f[2]+f[1])+3.0*f[0]); 
  fReflXYQuad[3][1] = 0.1666666666666667*(f[16]+1.732050807568877*(f[8]+f[7])+3.0*f[3]); 
  fReflXYQuad[3][2] = 0.1666666666666667*(f[17]+1.732050807568877*(f[10]+f[9])+3.0*f[4]); 
  fReflXYQuad[3][3] = 0.1666666666666667*(f[20]+1.732050807568877*(f[13]+f[12])+3.0*f[5]); 
  fReflXYQuad[3][4] = 0.1666666666666667*(f[26]+1.732050807568877*(f[19]+f[18])+3.0*f[11]); 
  fReflXYQuad[3][5] = 0.1666666666666667*(f[27]+1.732050807568877*(f[22]+f[21])+3.0*f[14]); 
  fReflXYQuad[3][6] = 0.1666666666666667*(f[28]+1.732050807568877*(f[24]+f[23])+3.0*f[15]); 
  fReflXYQuad[3][7] = 0.1666666666666667*(f[31]+1.732050807568877*(f[30]+f[29])+3.0*f[25]); 
  } else { // partial reflection 
  rVal = (0.5773502691896258*(1.732050807568877*f[27]+3.0*(f[22]+f[21]-1.0*f[16])+5.196152422706631*f[14]-1.0*(5.196152422706631*(f[8]+f[7])+9.0*f[3])))/(1.732050807568877*f[20]+3.0*(f[13]+f[12]-1.0*f[6])+5.196152422706631*f[5]-1.0*(5.196152422706631*(f[2]+f[1])+9.0*f[0])); 
  // if f is not realizable, no reflection 
  if(std::abs(rVal)>=.99 || -0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[13]+f[12]-1.0*f[6])+5.196152422706631*f[5]-1.0*(5.196152422706631*(f[2]+f[1])+9.0*f[0])) < 0.) { 
  fReflXYMuQuad[0][0] = 0.0; 
  fReflXYMuQuad[0][1] = 0.0; 
  fReflXYMuQuad[0][2] = 0.0; 
  fReflXYMuQuad[0][3] = 0.0; 
  } else {
   b = invL(rVal); 
   if(wv > 0) {
    xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
    fac = b>500? 0. : b<-500? 1. : abs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,1,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][0] = (-0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[13]+f[12]-1.0*f[6])+5.196152422706631*f[5]-1.0*(5.196152422706631*(f[2]+f[1])+9.0*f[0])))*fac; 
    fac = (b>500 || abs(b)<1e-10)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,2,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][1] = (-0.0392837100659193*(1.732050807568877*f[27]+3.0*(f[22]+f[21]-1.0*f[16])+5.196152422706631*f[14]-1.0*(5.196152422706631*(f[8]+f[7])+9.0*f[3])))*fac; 
   } else { 
    xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
    fac = b>500? 1. : b<-500? 0. : abs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,1,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][0] = (-0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[13]+f[12]-1.0*f[6])+5.196152422706631*f[5]-1.0*(5.196152422706631*(f[2]+f[1])+9.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || abs(b)<1e-10)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("1,2,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[0][1] = (-0.0392837100659193*(1.732050807568877*f[27]+3.0*(f[22]+f[21]-1.0*f[16])+5.196152422706631*f[14]-1.0*(5.196152422706631*(f[8]+f[7])+9.0*f[3])))*fac; 
   } 
  } 
  rVal = (0.5773502691896258*(1.732050807568877*f[27]+3.0*(f[22]+f[21]+f[16])+5.196152422706631*(f[14]+f[8]+f[7])+9.0*f[3]))/(1.732050807568877*f[20]+3.0*(f[13]+f[12]+f[6])+5.196152422706631*(f[5]+f[2]+f[1])+9.0*f[0]); 
  // if f is not realizable, no reflection 
  if(std::abs(rVal)>=.99 || 0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[13]+f[12]+f[6])+5.196152422706631*(f[5]+f[2]+f[1])+9.0*f[0]) < 0.) { 
  fReflXYMuQuad[1][0] = 0.0; 
  fReflXYMuQuad[1][1] = 0.0; 
  fReflXYMuQuad[1][2] = 0.0; 
  fReflXYMuQuad[1][3] = 0.0; 
  } else {
   b = invL(rVal); 
   if(wv > 0) {
    xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
    fac = b>500? 0. : b<-500? 1. : abs(b)<1e-10? (1.+xc)/2. : (exp(b*xc)-exp(-b))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,1,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][0] = (0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[13]+f[12]+f[6])+5.196152422706631*(f[5]+f[2]+f[1])+9.0*f[0]))*fac; 
    fac = (b>500 || abs(b)<1e-10)? 0. : b<-500? 1. : ((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,2,1: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][1] = (0.0392837100659193*(1.732050807568877*f[27]+3.0*(f[22]+f[21]+f[16])+5.196152422706631*(f[14]+f[8]+f[7])+9.0*f[3]))*fac; 
   } else { 
    xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
    fac = b>500? 1. : b<-500? 0. : abs(b)<1e-10? (1.-xc)/2. : (exp(b)-exp(b*xc))/(2.*sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,1,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][0] = (0.0392837100659193*(1.732050807568877*f[20]+3.0*(f[13]+f[12]+f[6])+5.196152422706631*(f[5]+f[2]+f[1])+9.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || abs(b)<1e-10)? 0. : ((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/2./(b*cosh(b)-sinh(b)); 
    if(fetestexcept(FE_OVERFLOW | FE_DIVBYZERO)) {printf("2,2,2: b = %Le, sinh(b) = %Le\n", b, sinh(b)); throw "Overflow!";}
    fReflXYMuQuad[1][1] = (0.0392837100659193*(1.732050807568877*f[27]+3.0*(f[22]+f[21]+f[16])+5.196152422706631*(f[14]+f[8]+f[7])+9.0*f[3]))*fac; 
   } 
  } 
  fReflXYQuad[3][0] = 0.7071067811865468*(fReflXYMuQuad[1][0]+fReflXYMuQuad[0][0]); 
  fReflXYQuad[3][1] = 0.7071067811865468*(fReflXYMuQuad[1][1]+fReflXYMuQuad[0][1]); 
  fReflXYQuad[3][2] = 0.7071067811865468*(fReflXYMuQuad[1][2]+fReflXYMuQuad[0][2]); 
  fReflXYQuad[3][3] = 1.224744871391589*(fReflXYMuQuad[1][0]-1.0*fReflXYMuQuad[0][0]); 
  fReflXYQuad[3][4] = 0.7071067811865468*(fReflXYMuQuad[1][3]+fReflXYMuQuad[0][3]); 
  fReflXYQuad[3][5] = 1.224744871391589*(fReflXYMuQuad[1][1]-1.0*fReflXYMuQuad[0][1]); 
  fReflXYQuad[3][6] = 1.224744871391589*(fReflXYMuQuad[1][2]-1.0*fReflXYMuQuad[0][2]); 
  fReflXYQuad[3][7] = 1.224744871391589*(fReflXYMuQuad[1][3]-1.0*fReflXYMuQuad[0][3]); 
  } 

 
  fRefl[0] = 0.5*(fReflXYQuad[3][0]+fReflXYQuad[2][0]+fReflXYQuad[1][0]+fReflXYQuad[0][0]); 
  fRefl[1] = 0.8660254037844385*(fReflXYQuad[3][0]-1.0*fReflXYQuad[2][0]+fReflXYQuad[1][0]-1.0*fReflXYQuad[0][0]); 
  fRefl[2] = 0.8660254037844385*(fReflXYQuad[3][0]+fReflXYQuad[2][0]-1.0*(fReflXYQuad[1][0]+fReflXYQuad[0][0])); 
  fRefl[3] = 0.5*(fReflXYQuad[3][1]+fReflXYQuad[2][1]+fReflXYQuad[1][1]+fReflXYQuad[0][1]); 
  fRefl[4] = 0.5*(fReflXYQuad[3][2]+fReflXYQuad[2][2]+fReflXYQuad[1][2]+fReflXYQuad[0][2]); 
  fRefl[5] = 0.5*(fReflXYQuad[3][3]+fReflXYQuad[2][3]+fReflXYQuad[1][3]+fReflXYQuad[0][3]); 
  fRefl[6] = 1.5*(fReflXYQuad[3][0]-1.0*(fReflXYQuad[2][0]+fReflXYQuad[1][0])+fReflXYQuad[0][0]); 
  fRefl[7] = 0.8660254037844385*(fReflXYQuad[3][1]-1.0*fReflXYQuad[2][1]+fReflXYQuad[1][1]-1.0*fReflXYQuad[0][1]); 
  fRefl[8] = 0.8660254037844385*(fReflXYQuad[3][1]+fReflXYQuad[2][1]-1.0*(fReflXYQuad[1][1]+fReflXYQuad[0][1])); 
  fRefl[9] = 0.8660254037844385*(fReflXYQuad[3][2]-1.0*fReflXYQuad[2][2]+fReflXYQuad[1][2]-1.0*fReflXYQuad[0][2]); 
  fRefl[10] = 0.8660254037844385*(fReflXYQuad[3][2]+fReflXYQuad[2][2]-1.0*(fReflXYQuad[1][2]+fReflXYQuad[0][2])); 
  fRefl[11] = 0.5*(fReflXYQuad[3][4]+fReflXYQuad[2][4]+fReflXYQuad[1][4]+fReflXYQuad[0][4]); 
  fRefl[12] = 0.8660254037844385*(fReflXYQuad[3][3]-1.0*fReflXYQuad[2][3]+fReflXYQuad[1][3]-1.0*fReflXYQuad[0][3]); 
  fRefl[13] = 0.8660254037844385*(fReflXYQuad[3][3]+fReflXYQuad[2][3]-1.0*(fReflXYQuad[1][3]+fReflXYQuad[0][3])); 
  fRefl[14] = 0.5*(fReflXYQuad[3][5]+fReflXYQuad[2][5]+fReflXYQuad[1][5]+fReflXYQuad[0][5]); 
  fRefl[15] = 0.5*(fReflXYQuad[3][6]+fReflXYQuad[2][6]+fReflXYQuad[1][6]+fReflXYQuad[0][6]); 
  fRefl[16] = 1.5*(fReflXYQuad[3][1]-1.0*(fReflXYQuad[2][1]+fReflXYQuad[1][1])+fReflXYQuad[0][1]); 
  fRefl[17] = 1.5*(fReflXYQuad[3][2]-1.0*(fReflXYQuad[2][2]+fReflXYQuad[1][2])+fReflXYQuad[0][2]); 
  fRefl[18] = 0.8660254037844385*(fReflXYQuad[3][4]-1.0*fReflXYQuad[2][4]+fReflXYQuad[1][4]-1.0*fReflXYQuad[0][4]); 
  fRefl[19] = 0.8660254037844385*(fReflXYQuad[3][4]+fReflXYQuad[2][4]-1.0*(fReflXYQuad[1][4]+fReflXYQuad[0][4])); 
  fRefl[20] = 1.5*(fReflXYQuad[3][3]-1.0*(fReflXYQuad[2][3]+fReflXYQuad[1][3])+fReflXYQuad[0][3]); 
  fRefl[21] = 0.8660254037844385*(fReflXYQuad[3][5]-1.0*fReflXYQuad[2][5]+fReflXYQuad[1][5]-1.0*fReflXYQuad[0][5]); 
  fRefl[22] = 0.8660254037844385*(fReflXYQuad[3][5]+fReflXYQuad[2][5]-1.0*(fReflXYQuad[1][5]+fReflXYQuad[0][5])); 
  fRefl[23] = 0.8660254037844385*(fReflXYQuad[3][6]-1.0*fReflXYQuad[2][6]+fReflXYQuad[1][6]-1.0*fReflXYQuad[0][6]); 
  fRefl[24] = 0.8660254037844385*(fReflXYQuad[3][6]+fReflXYQuad[2][6]-1.0*(fReflXYQuad[1][6]+fReflXYQuad[0][6])); 
  fRefl[25] = 0.5*(fReflXYQuad[3][7]+fReflXYQuad[2][7]+fReflXYQuad[1][7]+fReflXYQuad[0][7]); 
  fRefl[26] = 1.5*(fReflXYQuad[3][4]-1.0*(fReflXYQuad[2][4]+fReflXYQuad[1][4])+fReflXYQuad[0][4]); 
  fRefl[27] = 1.5*(fReflXYQuad[3][5]-1.0*(fReflXYQuad[2][5]+fReflXYQuad[1][5])+fReflXYQuad[0][5]); 
  fRefl[28] = 1.5*(fReflXYQuad[3][6]-1.0*(fReflXYQuad[2][6]+fReflXYQuad[1][6])+fReflXYQuad[0][6]); 
  fRefl[29] = 0.8660254037844385*(fReflXYQuad[3][7]-1.0*fReflXYQuad[2][7]+fReflXYQuad[1][7]-1.0*fReflXYQuad[0][7]); 
  fRefl[30] = 0.8660254037844385*(fReflXYQuad[3][7]+fReflXYQuad[2][7]-1.0*(fReflXYQuad[1][7]+fReflXYQuad[0][7])); 
  fRefl[31] = 1.5*(fReflXYQuad[3][7]-1.0*(fReflXYQuad[2][7]+fReflXYQuad[1][7])+fReflXYQuad[0][7]); 
}
