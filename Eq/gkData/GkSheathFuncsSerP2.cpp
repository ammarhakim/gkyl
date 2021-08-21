#include <GkSheathModDecl.h> 


void calcSheathReflection1x1vSer_P2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
  double fReflZQuad[3][8]; 
  

  vcutSq = (0.5*q_*(zVal*((9.48683298050514*phiWall[2]-9.48683298050514*phi[2])*zVal+4.898979485566357*phiWall[1]-4.898979485566357*phi[1])-3.16227766016838*phiWall[2]+3.16227766016838*phi[2]+2.828427124746191*phiWall[0]-2.828427124746191*phi[0]))/m_; 
  if(vcutSq <= vlowerSq) { // absorb (no reflection) 
  fRefl[0] = 0.0; 
  fRefl[1] = 0.0; 
  fRefl[2] = 0.0; 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.0; 
  fRefl[5] = 0.0; 
  fRefl[6] = 0.0; 
  fRefl[7] = 0.0; 
  } else if(vcutSq > vupperSq) { // full reflection 
  fRefl[0] = f[0]; 
  fRefl[1] = f[1]; 
  fRefl[2] = f[2]; 
  fRefl[3] = f[3]; 
  fRefl[4] = f[4]; 
  fRefl[5] = f[5]; 
  fRefl[6] = f[6]; 
  fRefl[7] = f[7]; 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))/(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflZQuad[0][0] = 0.0; 
  fReflZQuad[0][1] = 0.0; 
  fReflZQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.1414213562373095*(6.708203932499369*f[7]-5.0*f[5]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]-3.0*f[1])+5.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(5.0*f[2]-6.708203932499369*f[3])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[0][2] = (-0.1414213562373095*(6.708203932499369*f[7]-5.0*f[5]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(6.708203932499369*f[6]-6.0*f[2]))/(2.23606797749979*f[4]-2.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflZQuad[1][0] = 0.0; 
  fReflZQuad[1][1] = 0.0; 
  fReflZQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (-0.1178511301977579*(6.708203932499369*f[6]-6.0*f[2]))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (0.7071067811865475*f[5])*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][0] = (-0.3535533905932737*(2.23606797749979*f[4]-2.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][1] = (-0.1178511301977579*(6.708203932499369*f[6]-6.0*f[2]))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[1][2] = (0.7071067811865475*f[5])*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))/(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflZQuad[2][0] = 0.0; 
  fReflZQuad[2][1] = 0.0; 
  fReflZQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][2] = (0.1414213562373095*(6.708203932499369*f[7]+5.0*f[5]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][0] = (0.1414213562373095*(2.23606797749979*(2.0*f[4]+3.0*f[1])+5.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][1] = (0.04714045207910316*(13.41640786499874*f[6]+3.0*(6.708203932499369*f[3]+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZQuad[2][2] = (0.1414213562373095*(6.708203932499369*f[7]+5.0*f[5]))*fac; 
   } 
  } 
  fRefl[0] = 4.189452333721185e-16*(9.37681275180578e+14*fReflZQuad[2][0]+1.500290040288925e+15*fReflZQuad[1][0]+9.37681275180578e+14*fReflZQuad[0][0]); 
  fRefl[1] = 0.5270462766947305*(fReflZQuad[2][0]-1.0*fReflZQuad[0][0]); 
  fRefl[2] = 4.189452333721185e-16*(9.37681275180578e+14*fReflZQuad[2][1]+1.500290040288925e+15*fReflZQuad[1][1]+9.37681275180578e+14*fReflZQuad[0][1]); 
  fRefl[3] = 0.5270462766947305*(fReflZQuad[2][1]-1.0*fReflZQuad[0][1]); 
  fRefl[4] = 0.3513641844631534*(fReflZQuad[2][0]-2.0*fReflZQuad[1][0]+fReflZQuad[0][0]); 
  fRefl[5] = 4.189452333721185e-16*(9.37681275180578e+14*fReflZQuad[2][2]+1.500290040288925e+15*fReflZQuad[1][2]+9.37681275180578e+14*fReflZQuad[0][2]); 
  fRefl[6] = 0.3513641844631535*(fReflZQuad[2][1]-2.0*fReflZQuad[1][1]+fReflZQuad[0][1]); 
  fRefl[7] = 0.5270462766947306*(fReflZQuad[2][2]-1.0*fReflZQuad[0][2]); 
  } 

 
}

void calcSheathReflection1x2vSer_P2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq; long double xc, b, xbarVal, fac; 
  double fReflZMuQuad[8][8]; 
  

  vcutSq = (0.5*q_*(zVal*((9.48683298050514*phiWall[2]-9.48683298050514*phi[2])*zVal+4.898979485566357*phiWall[1]-4.898979485566357*phi[1])-3.16227766016838*phiWall[2]+3.16227766016838*phi[2]+2.828427124746191*phiWall[0]-2.828427124746191*phi[0]))/m_; 
  if(vcutSq <= vlowerSq) { // absorb (no reflection) 
  fRefl[0] = 0.0; 
  fRefl[1] = 0.0; 
  fRefl[2] = 0.0; 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.0; 
  fRefl[5] = 0.0; 
  fRefl[6] = 0.0; 
  fRefl[7] = 0.0; 
  fRefl[8] = 0.0; 
  fRefl[9] = 0.0; 
  fRefl[10] = 0.0; 
  fRefl[11] = 0.0; 
  fRefl[12] = 0.0; 
  fRefl[13] = 0.0; 
  fRefl[14] = 0.0; 
  fRefl[15] = 0.0; 
  fRefl[16] = 0.0; 
  fRefl[17] = 0.0; 
  fRefl[18] = 0.0; 
  fRefl[19] = 0.0; 
  } else if(vcutSq > vupperSq) { // full reflection 
  fRefl[0] = f[0]; 
  fRefl[1] = f[1]; 
  fRefl[2] = f[2]; 
  fRefl[3] = f[3]; 
  fRefl[4] = f[4]; 
  fRefl[5] = f[5]; 
  fRefl[6] = f[6]; 
  fRefl[7] = f[7]; 
  fRefl[8] = f[8]; 
  fRefl[9] = f[9]; 
  fRefl[10] = f[10]; 
  fRefl[11] = f[11]; 
  fRefl[12] = f[12]; 
  fRefl[13] = f[13]; 
  fRefl[14] = f[14]; 
  fRefl[15] = f[15]; 
  fRefl[16] = f[16]; 
  fRefl[17] = f[17]; 
  fRefl[18] = f[18]; 
  fRefl[19] = f[19]; 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))/(4.47213595499958*(6.708203932499369*(f[15]+f[13])-5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[5])-25.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])-5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[5])-25.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflZMuQuad[0][0] = 0.0; 
  fReflZMuQuad[0][1] = 0.0; 
  fReflZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (-0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])-5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[5])-25.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (-0.03333333333333333*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][0] = (-0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])-5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[5])-25.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][1] = (-0.03333333333333333*(2.0*(9.0*(f[19]+f[17])-6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[6]+f[4])-3.0*f[10])-5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[0][2] = (0.1*(9.0*f[18]-6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))/(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-5.0*f[7])+2.0*(5.0*f[0]-6.708203932499369*f[3])); 
  // if f is not realizable, no reflection from this node 
  if(0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-5.0*f[7])+2.0*(5.0*f[0]-6.708203932499369*f[3])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflZMuQuad[1][0] = 0.0; 
  fReflZMuQuad[1][1] = 0.0; 
  fReflZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-5.0*f[7])+2.0*(5.0*f[0]-6.708203932499369*f[3])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (0.01666666666666667*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.1*(6.708203932499369*f[14]-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][0] = (0.05*(2.23606797749979*(6.708203932499369*f[13]+4.0*f[9]-5.0*f[7])+2.0*(5.0*f[0]-6.708203932499369*f[3])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][1] = (0.01666666666666667*(45.0*f[17]+6.708203932499369*(4.0*f[16]-5.0*f[11])+6.0*(5.0*f[2]-6.708203932499369*f[6])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[1][2] = (-0.1*(6.708203932499369*f[14]-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))/(4.47213595499958*(6.708203932499369*(f[15]-1.0*f[13])+5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[5])+25.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.02*(4.47213595499958*(6.708203932499369*(f[15]-1.0*f[13])+5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[5])+25.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflZMuQuad[2][0] = 0.0; 
  fReflZMuQuad[2][1] = 0.0; 
  fReflZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (0.02*(4.47213595499958*(6.708203932499369*(f[15]-1.0*f[13])+5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[5])+25.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (0.03333333333333333*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.1*(9.0*f[18]+6.708203932499369*f[14]-1.0*(6.708203932499369*f[12]+5.0*f[8])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][0] = (0.02*(4.47213595499958*(6.708203932499369*(f[15]-1.0*f[13])+5.0*(f[9]+f[7]))+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[5])+25.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][1] = (0.03333333333333333*(2.0*(9.0*(f[19]-1.0*f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(2.23606797749979*(f[4]-1.0*f[6])-3.0*f[10])+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[2][2] = (-0.1*(9.0*f[18]+6.708203932499369*f[14]-1.0*(6.708203932499369*f[12]+5.0*f[8])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(5.0*(9.0*f[19]-6.708203932499369*f[16])+2.0*(13.41640786499874*f[11]+3.0*(5.0*f[2]-6.708203932499369*f[4]))))/(2.23606797749979*(6.708203932499369*f[15]-5.0*f[9])+2.0*(2.23606797749979*(2.0*f[7]-3.0*f[1])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.05*(2.23606797749979*(6.708203932499369*f[15]-5.0*f[9])+2.0*(2.23606797749979*(2.0*f[7]-3.0*f[1])+5.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflZMuQuad[3][0] = 0.0; 
  fReflZMuQuad[3][1] = 0.0; 
  fReflZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.05*(2.23606797749979*(6.708203932499369*f[15]-5.0*f[9])+2.0*(2.23606797749979*(2.0*f[7]-3.0*f[1])+5.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.01666666666666667*(5.0*(9.0*f[19]-6.708203932499369*f[16])+2.0*(13.41640786499874*f[11]+3.0*(5.0*f[2]-6.708203932499369*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.1*(6.708203932499369*f[12]-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][0] = (0.05*(2.23606797749979*(6.708203932499369*f[15]-5.0*f[9])+2.0*(2.23606797749979*(2.0*f[7]-3.0*f[1])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][1] = (0.01666666666666667*(5.0*(9.0*f[19]-6.708203932499369*f[16])+2.0*(13.41640786499874*f[11]+3.0*(5.0*f[2]-6.708203932499369*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[3][2] = (-0.1*(6.708203932499369*f[12]-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(5.0*(9.0*f[19]+6.708203932499369*f[16])-2.0*(13.41640786499874*f[11]+3.0*(6.708203932499369*f[4]+5.0*f[2]))))/(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9])-2.0*(2.23606797749979*(2.0*f[7]+3.0*f[1])+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9])-2.0*(2.23606797749979*(2.0*f[7]+3.0*f[1])+5.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflZMuQuad[4][0] = 0.0; 
  fReflZMuQuad[4][1] = 0.0; 
  fReflZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9])-2.0*(2.23606797749979*(2.0*f[7]+3.0*f[1])+5.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (-0.01666666666666667*(5.0*(9.0*f[19]+6.708203932499369*f[16])-2.0*(13.41640786499874*f[11]+3.0*(6.708203932499369*f[4]+5.0*f[2]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.1*(6.708203932499369*f[12]+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[15]+5.0*f[9])-2.0*(2.23606797749979*(2.0*f[7]+3.0*f[1])+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][1] = (-0.01666666666666667*(5.0*(9.0*f[19]+6.708203932499369*f[16])-2.0*(13.41640786499874*f[11]+3.0*(6.708203932499369*f[4]+5.0*f[2]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[4][2] = (0.1*(6.708203932499369*f[12]+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[4]-1.0*f[6]))-5.0*f[2])))/(4.47213595499958*(6.708203932499369*f[15]-1.0*(6.708203932499369*f[13]+5.0*(f[9]+f[7])))+3.0*(15.0*f[5]+11.18033988749895*(f[1]-1.0*f[3]))-25.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.02*(4.47213595499958*(6.708203932499369*f[15]-1.0*(6.708203932499369*f[13]+5.0*(f[9]+f[7])))+3.0*(15.0*f[5]+11.18033988749895*(f[1]-1.0*f[3]))-25.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflZMuQuad[5][0] = 0.0; 
  fReflZMuQuad[5][1] = 0.0; 
  fReflZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.02*(4.47213595499958*(6.708203932499369*f[15]-1.0*(6.708203932499369*f[13]+5.0*(f[9]+f[7])))+3.0*(15.0*f[5]+11.18033988749895*(f[1]-1.0*f[3]))-25.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.03333333333333333*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[4]-1.0*f[6]))-5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (-0.1*(9.0*f[18]+6.708203932499369*(f[12]-1.0*f[14])-5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][0] = (-0.02*(4.47213595499958*(6.708203932499369*f[15]-1.0*(6.708203932499369*f[13]+5.0*(f[9]+f[7])))+3.0*(15.0*f[5]+11.18033988749895*(f[1]-1.0*f[3]))-25.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][1] = (-0.03333333333333333*(2.0*(9.0*f[19]-1.0*(9.0*f[17]+6.708203932499369*(f[16]+f[11])))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[4]-1.0*f[6]))-5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[5][2] = (-0.1*(9.0*f[18]+6.708203932499369*(f[12]-1.0*f[14])-5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(45.0*f[17]+6.708203932499369*(5.0*f[11]-4.0*f[16])-6.0*(6.708203932499369*f[6]+5.0*f[2])))/(2.23606797749979*(6.708203932499369*f[13]-4.0*f[9]+5.0*f[7])-2.0*(6.708203932499369*f[3]+5.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.05*(2.23606797749979*(6.708203932499369*f[13]-4.0*f[9]+5.0*f[7])-2.0*(6.708203932499369*f[3]+5.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflZMuQuad[6][0] = 0.0; 
  fReflZMuQuad[6][1] = 0.0; 
  fReflZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[13]-4.0*f[9]+5.0*f[7])-2.0*(6.708203932499369*f[3]+5.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (-0.01666666666666667*(45.0*f[17]+6.708203932499369*(5.0*f[11]-4.0*f[16])-6.0*(6.708203932499369*f[6]+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (0.1*(6.708203932499369*f[14]+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][0] = (-0.05*(2.23606797749979*(6.708203932499369*f[13]-4.0*f[9]+5.0*f[7])-2.0*(6.708203932499369*f[3]+5.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][1] = (-0.01666666666666667*(45.0*f[17]+6.708203932499369*(5.0*f[11]-4.0*f[16])-6.0*(6.708203932499369*f[6]+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[6][2] = (0.1*(6.708203932499369*f[14]+5.0*f[8]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))/(4.47213595499958*(6.708203932499369*(f[15]+f[13])+5.0*(f[9]+f[7]))+3.0*(15.0*f[5]+11.18033988749895*(f[3]+f[1]))+25.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])+5.0*(f[9]+f[7]))+3.0*(15.0*f[5]+11.18033988749895*(f[3]+f[1]))+25.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflZMuQuad[7][0] = 0.0; 
  fReflZMuQuad[7][1] = 0.0; 
  fReflZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])+5.0*(f[9]+f[7]))+3.0*(15.0*f[5]+11.18033988749895*(f[3]+f[1]))+25.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (0.03333333333333333*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (0.1*(9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][0] = (0.02*(4.47213595499958*(6.708203932499369*(f[15]+f[13])+5.0*(f[9]+f[7]))+3.0*(15.0*f[5]+11.18033988749895*(f[3]+f[1]))+25.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][1] = (0.03333333333333333*(2.0*(9.0*(f[19]+f[17])+6.708203932499369*(f[16]+f[11]))+3.0*(3.0*(3.0*f[10]+2.23606797749979*(f[6]+f[4]))+5.0*f[2])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflZMuQuad[7][2] = (0.1*(9.0*f[18]+6.708203932499369*(f[14]+f[12])+5.0*f[8]))*fac; 
   } 
  } 
  fRefl[0] = 0.05555555555555555*(fReflZMuQuad[7][0]+8.0*fReflZMuQuad[6][0]+fReflZMuQuad[5][0]+8.0*(fReflZMuQuad[4][0]+fReflZMuQuad[3][0])+fReflZMuQuad[2][0]+8.0*fReflZMuQuad[1][0]+fReflZMuQuad[0][0]); 
  fRefl[1] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflZMuQuad[7][0]-4188761.0*fReflZMuQuad[6][0])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflZMuQuad[4][0]-1.0*fReflZMuQuad[3][0])-4.63256860547201e+14*fReflZMuQuad[5][0])+1.6692641e+7*(2.266096151179001e+23*fReflZMuQuad[2][0]-1.0*(8377522.0*fReflZMuQuad[1][0]+2.266096151179001e+23*fReflZMuQuad[0][0])))); 
  fRefl[2] = 0.05555555555555555*(fReflZMuQuad[7][1]+8.0*fReflZMuQuad[6][1]+fReflZMuQuad[5][1]+8.0*(fReflZMuQuad[4][1]+fReflZMuQuad[3][1])+fReflZMuQuad[2][1]+8.0*fReflZMuQuad[1][1]+fReflZMuQuad[0][1]); 
  fRefl[3] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflZMuQuad[7][0]+7.4121097687552e+14*fReflZMuQuad[6][0]+4.63256860547201e+14*fReflZMuQuad[5][0])-1.0*(9.988783372543001e+12*(fReflZMuQuad[4][0]+fReflZMuQuad[3][0])+4.80816459958139e+14*(4.63256860547201e+14*fReflZMuQuad[2][0]+7.4121097687552e+14*fReflZMuQuad[1][0]+4.63256860547201e+14*fReflZMuQuad[0][0]))); 
  fRefl[4] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflZMuQuad[7][1]-4188761.0*fReflZMuQuad[6][1])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflZMuQuad[4][1]-1.0*fReflZMuQuad[3][1])-4.63256860547201e+14*fReflZMuQuad[5][1])+1.6692641e+7*(2.266096151179001e+23*fReflZMuQuad[2][1]-1.0*(8377522.0*fReflZMuQuad[1][1]+2.266096151179001e+23*fReflZMuQuad[0][1])))); 
  fRefl[5] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflZMuQuad[7][0]+9.0*fReflZMuQuad[6][0])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflZMuQuad[4][0]-1.346286087882789e+17*fReflZMuQuad[5][0])+1.346286087882789e+16*(9.93730136036331e+15*fReflZMuQuad[0][0]-1.0*(9.93730136036331e+15*fReflZMuQuad[2][0]+fReflZMuQuad[1][0])))); 
  fRefl[6] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflZMuQuad[7][1]+7.4121097687552e+14*fReflZMuQuad[6][1]+4.63256860547201e+14*fReflZMuQuad[5][1])-1.0*(9.988783372543001e+12*(fReflZMuQuad[4][1]+fReflZMuQuad[3][1])+4.80816459958139e+14*(4.63256860547201e+14*fReflZMuQuad[2][1]+7.4121097687552e+14*fReflZMuQuad[1][1]+4.63256860547201e+14*fReflZMuQuad[0][1]))); 
  fRefl[7] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflZMuQuad[7][0]-4.17103711299121e+14*fReflZMuQuad[6][0])+5.028593299999999e+7*(3.190559553141742e+22*fReflZMuQuad[5][0]+2384663.0*(fReflZMuQuad[4][0]-1.0*fReflZMuQuad[3][0])+3.190559553141742e+22*(fReflZMuQuad[2][0]-2.0*fReflZMuQuad[1][0]+fReflZMuQuad[0][0]))); 
  fRefl[8] = 0.05555555555555555*(fReflZMuQuad[7][2]+8.0*fReflZMuQuad[6][2]+fReflZMuQuad[5][2]+8.0*(fReflZMuQuad[4][2]+fReflZMuQuad[3][2])+fReflZMuQuad[2][2]+8.0*fReflZMuQuad[1][2]+fReflZMuQuad[0][2]); 
  fRefl[9] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflZMuQuad[7][0]+1.6692641e+7*fReflZMuQuad[6][0])+2.087265252531456e+15*(1.881394845173929e+38*fReflZMuQuad[5][0]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflZMuQuad[2][0]+2384663.0*fReflZMuQuad[1][0]+3.190559553141742e+22*fReflZMuQuad[0][0])-3.20880527843592e+30*(fReflZMuQuad[4][0]+fReflZMuQuad[3][0])))); 
  fRefl[10] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflZMuQuad[7][1]+9.0*fReflZMuQuad[6][1])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflZMuQuad[4][1]-1.346286087882789e+17*fReflZMuQuad[5][1])+1.346286087882789e+16*(9.93730136036331e+15*fReflZMuQuad[0][1]-1.0*(9.93730136036331e+15*fReflZMuQuad[2][1]+fReflZMuQuad[1][1])))); 
  fRefl[11] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflZMuQuad[7][1]-4.17103711299121e+14*fReflZMuQuad[6][1])+5.028593299999999e+7*(3.190559553141743e+22*fReflZMuQuad[5][1]+2384663.0*(fReflZMuQuad[4][1]-1.0*fReflZMuQuad[3][1])+3.190559553141743e+22*(fReflZMuQuad[2][1]+fReflZMuQuad[0][1]))); 
  fRefl[12] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflZMuQuad[7][2]-4188761.0*fReflZMuQuad[6][2])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflZMuQuad[4][2]-4.63256860547201e+14*fReflZMuQuad[5][2])+1.6692641e+7*(2.266096151179001e+23*fReflZMuQuad[2][2]-1.0*(8377522.0*fReflZMuQuad[1][2]+2.266096151179001e+23*fReflZMuQuad[0][2])))); 
  fRefl[13] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflZMuQuad[7][0]-2.0*fReflZMuQuad[6][0]+fReflZMuQuad[5][0])-27.0*fReflZMuQuad[3][0]+1.077028870306231e+18*((-1.0*fReflZMuQuad[2][0])+2.0*fReflZMuQuad[1][0]-1.0*fReflZMuQuad[0][0])); 
  fRefl[14] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflZMuQuad[7][2]+7.4121097687552e+14*fReflZMuQuad[6][2]+4.63256860547201e+14*fReflZMuQuad[5][2])-1.0*(9.988783372543001e+12*(fReflZMuQuad[4][2]+fReflZMuQuad[3][2])+4.80816459958139e+14*(4.63256860547201e+14*fReflZMuQuad[2][2]+7.4121097687552e+14*fReflZMuQuad[1][2]+4.63256860547201e+14*fReflZMuQuad[0][2]))); 
  fRefl[15] = 1.719407810605221e-19*(1.077028870306231e+18*fReflZMuQuad[7][0]+27.0*fReflZMuQuad[6][0]+1.077028870306231e+18*((-1.0*fReflZMuQuad[5][0])+2.0*(fReflZMuQuad[3][0]-1.0*fReflZMuQuad[4][0])+fReflZMuQuad[2][0]-1.0*fReflZMuQuad[0][0])); 
  fRefl[16] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflZMuQuad[7][1]+1.6692641e+7*fReflZMuQuad[6][1])+2.087265252531456e+15*(1.881394845173929e+38*fReflZMuQuad[5][1]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflZMuQuad[2][1]+2384663.0*fReflZMuQuad[1][1]+3.190559553141743e+22*fReflZMuQuad[0][1])-3.20880527843592e+30*(fReflZMuQuad[4][1]+fReflZMuQuad[3][1])))); 
  fRefl[17] = 1.719407810605222e-19*(1.077028870306231e+18*(fReflZMuQuad[7][1]+fReflZMuQuad[5][1])-27.0*fReflZMuQuad[3][1]+2.154057740612463e+17*((-5.0*fReflZMuQuad[2][1])+10.0*fReflZMuQuad[1][1]-5.0*fReflZMuQuad[0][1])); 
  fRefl[18] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflZMuQuad[7][2]+9.0*fReflZMuQuad[6][2])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflZMuQuad[4][2]-1.346286087882789e+17*fReflZMuQuad[5][2])+1.346286087882789e+16*(9.93730136036331e+15*fReflZMuQuad[0][2]-1.0*(9.93730136036331e+15*fReflZMuQuad[2][2]+fReflZMuQuad[1][2])))); 
  fRefl[19] = 1.719407810605222e-19*(1.077028870306231e+18*fReflZMuQuad[7][1]+27.0*fReflZMuQuad[6][1]+2.154057740612463e+17*((-5.0*fReflZMuQuad[5][1])+2.0*(5.0*fReflZMuQuad[3][1]-5.0*fReflZMuQuad[4][1])+5.0*fReflZMuQuad[2][1])); 
  } 

 
}

void calcSheathReflection3x2vSer_P2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq_i; long double xc, b, xbarVal, fac; 
  double fReflXYQuad[8][20]; 
  double fReflXYZMuQuad[8][8]; 
  

// node (x,y)_1 
  vcutSq_i = (0.01*q_*(zVal*((426.9074841227313*phiWall[19]-426.9074841227313*phi[19]-318.1980515339465*phiWall[16]+318.1980515339465*phi[16]-318.1980515339465*phiWall[15]+318.1980515339465*phi[15]+237.1708245126285*phiWall[9]-237.1708245126285*phi[9])*zVal-146.9693845669907*phiWall[18]+146.9693845669907*phi[18]-146.9693845669907*phiWall[17]+146.9693845669907*phi[17]+109.5445115010333*phiWall[14]-109.5445115010333*phi[14]+109.5445115010333*phiWall[13]-109.5445115010333*phi[13]+220.454076850486*phiWall[10]-220.454076850486*phi[10]-164.3167672515499*phiWall[6]+164.3167672515499*phi[6]-164.3167672515499*phiWall[5]+164.3167672515499*phi[5]+122.4744871391589*phiWall[3]-122.4744871391589*phi[3])-142.3024947075771*phiWall[19]+142.3024947075771*phi[19]+106.0660171779822*phiWall[16]-106.0660171779822*phi[16]+106.0660171779822*phiWall[15]-106.0660171779822*phi[15]-84.85281374238573*phiWall[12]+84.85281374238573*phi[12]-84.85281374238573*phiWall[11]+84.85281374238573*phi[11]-79.0569415042095*phiWall[9]+79.0569415042095*phi[9]+63.24555320336762*phiWall[8]-63.24555320336762*phi[8]+63.24555320336762*phiWall[7]-63.24555320336762*phi[7]+127.2792206135786*phiWall[4]-127.2792206135786*phi[4]-94.86832980505142*phiWall[2]+94.86832980505142*phi[2]-94.86832980505142*phiWall[1]+94.86832980505142*phi[1]+70.71067811865477*phiWall[0]-70.71067811865477*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[0][0] = 0.0; 
  fReflXYQuad[0][1] = 0.0; 
  fReflXYQuad[0][2] = 0.0; 
  fReflXYQuad[0][3] = 0.0; 
  fReflXYQuad[0][4] = 0.0; 
  fReflXYQuad[0][5] = 0.0; 
  fReflXYQuad[0][6] = 0.0; 
  fReflXYQuad[0][7] = 0.0; 
  fReflXYQuad[0][8] = 0.0; 
  fReflXYQuad[0][9] = 0.0; 
  fReflXYQuad[0][10] = 0.0; 
  fReflXYQuad[0][11] = 0.0; 
  fReflXYQuad[0][12] = 0.0; 
  fReflXYQuad[0][13] = 0.0; 
  fReflXYQuad[0][14] = 0.0; 
  fReflXYQuad[0][15] = 0.0; 
  fReflXYQuad[0][16] = 0.0; 
  fReflXYQuad[0][17] = 0.0; 
  fReflXYQuad[0][18] = 0.0; 
  fReflXYQuad[0][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[0][0] = -0.02*(4.47213595499958*(6.708203932499369*(f[32]+f[31])-5.0*(f[17]+f[16]))+3.0*(11.18033988749895*(f[2]+f[1])-15.0*f[6])-25.0*f[0]); 
  fReflXYQuad[0][1] = -0.03333333333333333*(2.0*(9.0*(f[57]+f[56])-6.708203932499369*(f[34]+f[33]))+3.0*(3.0*(2.23606797749979*(f[8]+f[7])-3.0*f[21])-5.0*f[3])); 
  fReflXYQuad[0][2] = -0.03333333333333333*(2.0*(9.0*(f[60]+f[59])-6.708203932499369*(f[38]+f[37]))+3.0*(3.0*(2.23606797749979*(f[10]+f[9])-3.0*f[22])-5.0*f[4])); 
  fReflXYQuad[0][3] = -0.03333333333333333*(2.0*(9.0*(f[69]+f[68])-6.708203932499369*(f[44]+f[43]))+3.0*(3.0*(2.23606797749979*(f[13]+f[12])-3.0*f[25])-5.0*f[5])); 
  fReflXYQuad[0][4] = -0.02*(4.47213595499958*(6.708203932499369*(f[88]+f[87])-5.0*(f[62]+f[61]))+3.0*(11.18033988749895*(f[24]+f[23])-15.0*f[51])-25.0*f[11]); 
  fReflXYQuad[0][5] = -0.02*(4.47213595499958*(6.708203932499369*(f[92]+f[91])-5.0*(f[71]+f[70]))+3.0*(11.18033988749895*(f[27]+f[26])-15.0*f[52])-25.0*f[14]); 
  fReflXYQuad[0][6] = -0.02*(4.47213595499958*(6.708203932499369*(f[95]+f[94])-5.0*(f[75]+f[74]))+3.0*(11.18033988749895*(f[29]+f[28])-15.0*f[53])-25.0*f[15]); 
  fReflXYQuad[0][7] = 0.1*(9.0*f[58]-6.708203932499369*(f[36]+f[35])+5.0*f[18]); 
  fReflXYQuad[0][8] = 0.1*(9.0*f[65]-6.708203932499369*(f[41]+f[40])+5.0*f[19]); 
  fReflXYQuad[0][9] = 0.1*(9.0*f[80]-6.708203932499369*(f[48]+f[47])+5.0*f[20]); 
  fReflXYQuad[0][10] = -0.03333333333333333*(2.0*(9.0*(f[108]+f[107])-6.708203932499369*(f[97]+f[96]))+3.0*(3.0*(2.23606797749979*(f[55]+f[54])-3.0*f[86])-5.0*f[30])); 
  fReflXYQuad[0][11] = 0.1*(9.0*f[89]-6.708203932499369*(f[64]+f[63])+5.0*f[39]); 
  fReflXYQuad[0][12] = 0.1*(9.0*f[90]-6.708203932499369*(f[67]+f[66])+5.0*f[42]); 
  fReflXYQuad[0][13] = 0.1*(9.0*f[93]-6.708203932499369*(f[73]+f[72])+5.0*f[45]); 
  fReflXYQuad[0][14] = 0.1*(9.0*f[100]-6.708203932499369*(f[78]+f[77])+5.0*f[46]); 
  fReflXYQuad[0][15] = 0.1*(9.0*f[103]-6.708203932499369*(f[82]+f[81])+5.0*f[49]); 
  fReflXYQuad[0][16] = 0.1*(9.0*f[104]-6.708203932499369*(f[84]+f[83])+5.0*f[50]); 
  fReflXYQuad[0][17] = 0.1*(9.0*f[109]-6.708203932499369*(f[99]+f[98])+5.0*f[76]); 
  fReflXYQuad[0][18] = 0.1*(9.0*f[110]-6.708203932499369*(f[102]+f[101])+5.0*f[79]); 
  fReflXYQuad[0][19] = 0.1*(9.0*f[111]-6.708203932499369*(f[106]+f[105])+5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])-60.37383539249431*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87]))+9.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+20.12461179749811*(f[55]+f[54]+f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(11.18033988749895*(f[15]+f[11]+f[10]+f[9])-15.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22]))-25.0*f[4])))/(269.9999999999999*(f[103]+f[93]+f[92]+f[91])-9.0*(22.3606797749979*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*(f[27]+f[26]+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(55.90169943749476*(f[5]+f[3]+f[2]+f[1])-75.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6]))-125.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(269.9999999999999*(f[103]+f[93]+f[92]+f[91])-9.0*(22.3606797749979*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*(f[27]+f[26]+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(55.90169943749476*(f[5]+f[3]+f[2]+f[1])-75.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6]))-125.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.002*(269.9999999999999*(f[103]+f[93]+f[92]+f[91])-9.0*(22.3606797749979*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*(f[27]+f[26]+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(55.90169943749476*(f[5]+f[3]+f[2]+f[1])-75.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6]))-125.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])-60.37383539249431*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87]))+9.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+20.12461179749811*(f[55]+f[54]+f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(11.18033988749895*(f[15]+f[11]+f[10]+f[9])-15.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22]))-25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(81.0*f[110]-60.37383539249431*(f[102]+f[101]+f[100]+f[90])+5.0*(9.0*(f[79]+f[78]+f[77]+f[67]+f[66]+f[65])-6.708203932499369*(f[46]+f[42]+f[41]+f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.002*(269.9999999999999*(f[103]+f[93]+f[92]+f[91])-9.0*(22.3606797749979*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*(f[27]+f[26]+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(55.90169943749476*(f[5]+f[3]+f[2]+f[1])-75.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6]))-125.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])-60.37383539249431*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87]))+9.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+20.12461179749811*(f[55]+f[54]+f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(11.18033988749895*(f[15]+f[11]+f[10]+f[9])-15.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22]))-25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(81.0*f[110]-60.37383539249431*(f[102]+f[101]+f[100]+f[90])+5.0*(9.0*(f[79]+f[78]+f[77]+f[67]+f[66]+f[65])-6.708203932499369*(f[46]+f[42]+f[41]+f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+60.37383539249431*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94])-5.0*f[89])+9.0*(5.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))-2.0*(10.0*(f[60]+f[59])+20.12461179749811*f[53]))+33.54101966249684*(4.0*f[50]-5.0*f[39])+2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*(f[29]+f[28]+f[22])-11.18033988749895*(f[15]+f[10]+f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-5.0*f[58])+6.708203932499369*((-4.0*(f[48]+f[47]))+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))-2.0*(13.41640786499874*(f[32]+f[31])+27.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[13]+f[12]+f[6])-11.18033988749895*(f[5]+f[2]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-5.0*f[58])+6.708203932499369*((-4.0*(f[48]+f[47]))+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))-2.0*(13.41640786499874*(f[32]+f[31])+27.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[13]+f[12]+f[6])-11.18033988749895*(f[5]+f[2]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-5.0*f[58])+6.708203932499369*((-4.0*(f[48]+f[47]))+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))-2.0*(13.41640786499874*(f[32]+f[31])+27.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[13]+f[12]+f[6])-11.18033988749895*(f[5]+f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.001666666666666667*(405.0*f[109]+60.37383539249431*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94])-5.0*f[89])+9.0*(5.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))-2.0*(10.0*(f[60]+f[59])+20.12461179749811*f[53]))+33.54101966249684*(4.0*f[50]-5.0*f[39])+2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*(f[29]+f[28]+f[22])-11.18033988749895*(f[15]+f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01*(60.37383539249431*f[100]+5.0*((-9.0*(f[78]+f[77]+f[65]))+6.708203932499369*(f[46]+f[41]+f[40])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]-5.0*(f[73]+f[72])+4.0*(f[69]+f[68])-5.0*f[58])+6.708203932499369*((-4.0*(f[48]+f[47]))+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))-2.0*(13.41640786499874*(f[32]+f[31])+27.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[13]+f[12]+f[6])-11.18033988749895*(f[5]+f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.001666666666666667*(405.0*f[109]+60.37383539249431*(4.0*f[104]-5.0*(f[99]+f[98])+4.0*(f[95]+f[94])-5.0*f[89])+9.0*(5.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))-2.0*(10.0*(f[60]+f[59])+20.12461179749811*f[53]))+33.54101966249684*(4.0*f[50]-5.0*f[39])+2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*(f[29]+f[28]+f[22])-11.18033988749895*(f[15]+f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01*(60.37383539249431*f[100]+5.0*((-9.0*(f[78]+f[77]+f[65]))+6.708203932499369*(f[46]+f[41]+f[40])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])+60.37383539249431*((-1.0*(f[106]+f[105]))+f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94]+f[89])+f[88]+f[87])))+9.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+20.12461179749811*(f[55]+f[54]-1.0*f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[30])+f[29]+f[28]-1.0*(f[24]+f[23]-1.0*f[22]))+11.18033988749895*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9])))+25.0*f[4])))/(269.9999999999999*(f[103]-1.0*f[93]+f[92]+f[91])+9.0*(22.3606797749979*((-1.0*(f[82]+f[81]))+f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68]+f[58])))-1.0*(22.3606797749979*(f[57]+f[56])+45.0*f[52]))+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+27.0*(f[27]+f[26]-1.0*f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]+f[12]-1.0*(f[8]+f[7]-1.0*f[6]))+55.90169943749476*((-1.0*f[5])+f[3]-1.0*(f[2]+f[1])))+125.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(269.9999999999999*(f[103]-1.0*f[93]+f[92]+f[91])+9.0*(22.3606797749979*((-1.0*(f[82]+f[81]))+f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68]+f[58])))-1.0*(22.3606797749979*(f[57]+f[56])+45.0*f[52]))+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+27.0*(f[27]+f[26]-1.0*f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]+f[12]-1.0*(f[8]+f[7]-1.0*f[6]))+55.90169943749476*((-1.0*f[5])+f[3]-1.0*(f[2]+f[1])))+125.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.002*(269.9999999999999*(f[103]-1.0*f[93]+f[92]+f[91])+9.0*(22.3606797749979*((-1.0*(f[82]+f[81]))+f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68]+f[58])))-1.0*(22.3606797749979*(f[57]+f[56])+45.0*f[52]))+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+27.0*(f[27]+f[26]-1.0*f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]+f[12]-1.0*(f[8]+f[7]-1.0*f[6]))+55.90169943749476*((-1.0*f[5])+f[3]-1.0*(f[2]+f[1])))+125.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])+60.37383539249431*((-1.0*(f[106]+f[105]))+f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94]+f[89])+f[88]+f[87])))+9.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+20.12461179749811*(f[55]+f[54]-1.0*f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[30])+f[29]+f[28]-1.0*(f[24]+f[23]-1.0*f[22]))+11.18033988749895*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9])))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(81.0*f[110]+60.37383539249431*((-1.0*(f[102]+f[101]))+f[100]-1.0*f[90])+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]-1.0*f[67])+f[66]-1.0*f[65])+6.708203932499369*(f[46]-1.0*f[42]+f[41]+f[40])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.002*(269.9999999999999*(f[103]-1.0*f[93]+f[92]+f[91])+9.0*(22.3606797749979*((-1.0*(f[82]+f[81]))+f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68]+f[58])))-1.0*(22.3606797749979*(f[57]+f[56])+45.0*f[52]))+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+27.0*(f[27]+f[26]-1.0*f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]+f[12]-1.0*(f[8]+f[7]-1.0*f[6]))+55.90169943749476*((-1.0*f[5])+f[3]-1.0*(f[2]+f[1])))+125.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])+60.37383539249431*((-1.0*(f[106]+f[105]))+f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94]+f[89])+f[88]+f[87])))+9.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+20.12461179749811*(f[55]+f[54]-1.0*f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[30])+f[29]+f[28]-1.0*(f[24]+f[23]-1.0*f[22]))+11.18033988749895*((-1.0*f[15])+f[11]-1.0*(f[10]+f[9])))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(81.0*f[110]+60.37383539249431*((-1.0*(f[102]+f[101]))+f[100]-1.0*f[90])+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]-1.0*f[67])+f[66]-1.0*f[65])+6.708203932499369*(f[46]-1.0*f[42]+f[41]+f[40])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+60.37383539249431*(4.0*(f[89]+f[88]+f[87])-5.0*(f[106]+f[105]+f[104]))+225.0*(f[85]+f[84]+f[83])-1.0*(18.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+20.12461179749811*f[51])+167.7050983124842*f[50])+2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[24]+f[23]+f[22])-11.18033988749895*(f[11]+f[10]+f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[103]+9.0*(4.0*(f[58]+f[57]+f[56])-5.0*(f[82]+f[81]+f[80]))+33.54101966249684*(f[49]+f[48]+f[47])-1.0*(2.0*(13.41640786499874*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*f[21])+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[8]+f[7]+f[6])-11.18033988749895*(f[3]+f[2]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(4.0*(f[58]+f[57]+f[56])-5.0*(f[82]+f[81]+f[80]))+33.54101966249684*(f[49]+f[48]+f[47])-1.0*(2.0*(13.41640786499874*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*f[21])+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[8]+f[7]+f[6])-11.18033988749895*(f[3]+f[2]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(4.0*(f[58]+f[57]+f[56])-5.0*(f[82]+f[81]+f[80]))+33.54101966249684*(f[49]+f[48]+f[47])-1.0*(2.0*(13.41640786499874*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*f[21])+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[8]+f[7]+f[6])-11.18033988749895*(f[3]+f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.001666666666666667*(405.0*f[111]+60.37383539249431*(4.0*(f[89]+f[88]+f[87])-5.0*(f[106]+f[105]+f[104]))+225.0*(f[85]+f[84]+f[83])-1.0*(18.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+20.12461179749811*f[51])+167.7050983124842*f[50])+2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[24]+f[23]+f[22])-11.18033988749895*(f[11]+f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01*(60.37383539249431*f[90]+5.0*((-9.0*(f[67]+f[66]+f[65]))+6.708203932499369*(f[42]+f[41]+f[40])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(4.0*(f[58]+f[57]+f[56])-5.0*(f[82]+f[81]+f[80]))+33.54101966249684*(f[49]+f[48]+f[47])-1.0*(2.0*(13.41640786499874*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*f[21])+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[8]+f[7]+f[6])-11.18033988749895*(f[3]+f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.001666666666666667*(405.0*f[111]+60.37383539249431*(4.0*(f[89]+f[88]+f[87])-5.0*(f[106]+f[105]+f[104]))+225.0*(f[85]+f[84]+f[83])-1.0*(18.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+20.12461179749811*f[51])+167.7050983124842*f[50])+2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[24]+f[23]+f[22])-11.18033988749895*(f[11]+f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01*(60.37383539249431*f[90]+5.0*((-9.0*(f[67]+f[66]+f[65]))+6.708203932499369*(f[42]+f[41]+f[40])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+60.37383539249431*(5.0*(f[104]-1.0*(f[106]+f[105]))+4.0*((-1.0*f[89])+f[88]+f[87]))+9.0*(25.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-20.12461179749811*f[51]))+167.7050983124842*f[50]+2.0*(3.0*(3.0*(15.0*(f[24]+f[23]-1.0*f[22])+11.18033988749895*((-1.0*f[11])+f[10]+f[9]))-25.0*f[4])-67.08203932499369*(f[39]+f[38]+f[37]))))/(2.23606797749979*(60.37383539249431*f[103]+9.0*((-5.0*(f[82]+f[81]))+5.0*f[80]+4.0*((-1.0*f[58])+f[57]+f[56]))+6.708203932499369*(5.0*f[49]-5.0*(f[48]+f[47]))+2.0*(13.41640786499874*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))-27.0*f[21])+25.0*f[20])+2.0*((-22.3606797749979*(f[18]+f[17]+f[16]))+3.0*(15.0*(f[8]+f[7]-1.0*f[6])+11.18033988749895*((-1.0*f[3])+f[2]+f[1]))-25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*((-5.0*(f[82]+f[81]))+5.0*f[80]+4.0*((-1.0*f[58])+f[57]+f[56]))+6.708203932499369*(5.0*f[49]-5.0*(f[48]+f[47]))+2.0*(13.41640786499874*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))-27.0*f[21])+25.0*f[20])+2.0*((-22.3606797749979*(f[18]+f[17]+f[16]))+3.0*(15.0*(f[8]+f[7]-1.0*f[6])+11.18033988749895*((-1.0*f[3])+f[2]+f[1]))-25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*((-5.0*(f[82]+f[81]))+5.0*f[80]+4.0*((-1.0*f[58])+f[57]+f[56]))+6.708203932499369*(5.0*f[49]-5.0*(f[48]+f[47]))+2.0*(13.41640786499874*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))-27.0*f[21])+25.0*f[20])+2.0*((-22.3606797749979*(f[18]+f[17]+f[16]))+3.0*(15.0*(f[8]+f[7]-1.0*f[6])+11.18033988749895*((-1.0*f[3])+f[2]+f[1]))-25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*(f[104]-1.0*(f[106]+f[105]))+4.0*((-1.0*f[89])+f[88]+f[87]))+9.0*(25.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-20.12461179749811*f[51]))+167.7050983124842*f[50]+2.0*(3.0*(3.0*(15.0*(f[24]+f[23]-1.0*f[22])+11.18033988749895*((-1.0*f[11])+f[10]+f[9]))-25.0*f[4])-67.08203932499369*(f[39]+f[38]+f[37]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.01*(60.37383539249431*f[90]+5.0*(9.0*(f[65]-1.0*(f[67]+f[66]))+6.708203932499369*(f[42]-1.0*(f[41]+f[40]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*((-5.0*(f[82]+f[81]))+5.0*f[80]+4.0*((-1.0*f[58])+f[57]+f[56]))+6.708203932499369*(5.0*f[49]-5.0*(f[48]+f[47]))+2.0*(13.41640786499874*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31])))-27.0*f[21])+25.0*f[20])+2.0*((-22.3606797749979*(f[18]+f[17]+f[16]))+3.0*(15.0*(f[8]+f[7]-1.0*f[6])+11.18033988749895*((-1.0*f[3])+f[2]+f[1]))-25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*(f[104]-1.0*(f[106]+f[105]))+4.0*((-1.0*f[89])+f[88]+f[87]))+9.0*(25.0*(f[85]-1.0*(f[84]+f[83]))+2.0*(10.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59])))-20.12461179749811*f[51]))+167.7050983124842*f[50]+2.0*(3.0*(3.0*(15.0*(f[24]+f[23]-1.0*f[22])+11.18033988749895*((-1.0*f[11])+f[10]+f[9]))-25.0*f[4])-67.08203932499369*(f[39]+f[38]+f[37]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.01*(60.37383539249431*f[90]+5.0*(9.0*(f[65]-1.0*(f[67]+f[66]))+6.708203932499369*(f[42]-1.0*(f[41]+f[40]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]+f[107]))+60.37383539249431*((-1.0*(f[106]+f[105]+f[104]-1.0*f[99]))+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+9.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-20.12461179749811*(f[55]+f[54]+f[53]-1.0*f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]+f[29]+f[28]-1.0*(f[24]+f[23]+f[22]))+11.18033988749895*((-1.0*f[15])+f[11]+f[10]+f[9]))-25.0*f[4])))/(269.9999999999999*(f[103]-1.0*(f[93]+f[92]+f[91]))+9.0*(22.3606797749979*((-1.0*(f[82]+f[81]+f[80]-1.0*f[73]))+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31])-1.0*(27.0*(f[27]+f[26]+f[25]-1.0*f[21])+10.0*(f[20]+f[18]+f[17]+f[16])))+3.0*(75.0*(f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6]))+55.90169943749476*((-1.0*f[5])+f[3]+f[2]+f[1]))-125.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(269.9999999999999*(f[103]-1.0*(f[93]+f[92]+f[91]))+9.0*(22.3606797749979*((-1.0*(f[82]+f[81]+f[80]-1.0*f[73]))+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31])-1.0*(27.0*(f[27]+f[26]+f[25]-1.0*f[21])+10.0*(f[20]+f[18]+f[17]+f[16])))+3.0*(75.0*(f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6]))+55.90169943749476*((-1.0*f[5])+f[3]+f[2]+f[1]))-125.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.002*(269.9999999999999*(f[103]-1.0*(f[93]+f[92]+f[91]))+9.0*(22.3606797749979*((-1.0*(f[82]+f[81]+f[80]-1.0*f[73]))+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31])-1.0*(27.0*(f[27]+f[26]+f[25]-1.0*f[21])+10.0*(f[20]+f[18]+f[17]+f[16])))+3.0*(75.0*(f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6]))+55.90169943749476*((-1.0*f[5])+f[3]+f[2]+f[1]))-125.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]+f[107]))+60.37383539249431*((-1.0*(f[106]+f[105]+f[104]-1.0*f[99]))+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+9.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-20.12461179749811*(f[55]+f[54]+f[53]-1.0*f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]+f[29]+f[28]-1.0*(f[24]+f[23]+f[22]))+11.18033988749895*((-1.0*f[15])+f[11]+f[10]+f[9]))-25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(81.0*f[110]-60.37383539249431*(f[102]+f[101]+f[100]-1.0*f[90])+5.0*(9.0*(f[79]+f[78]+f[77]-1.0*(f[67]+f[66]+f[65]))+6.708203932499369*((-1.0*f[46])+f[42]+f[41]+f[40])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.002*(269.9999999999999*(f[103]-1.0*(f[93]+f[92]+f[91]))+9.0*(22.3606797749979*((-1.0*(f[82]+f[81]+f[80]-1.0*f[73]))+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31])-1.0*(27.0*(f[27]+f[26]+f[25]-1.0*f[21])+10.0*(f[20]+f[18]+f[17]+f[16])))+3.0*(75.0*(f[14]+f[13]+f[12]-1.0*(f[8]+f[7]+f[6]))+55.90169943749476*((-1.0*f[5])+f[3]+f[2]+f[1]))-125.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]+f[107]))+60.37383539249431*((-1.0*(f[106]+f[105]+f[104]-1.0*f[99]))+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87])))+9.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-20.12461179749811*(f[55]+f[54]+f[53]-1.0*f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]+f[29]+f[28]-1.0*(f[24]+f[23]+f[22]))+11.18033988749895*((-1.0*f[15])+f[11]+f[10]+f[9]))-25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(81.0*f[110]-60.37383539249431*(f[102]+f[101]+f[100]-1.0*f[90])+5.0*(9.0*(f[79]+f[78]+f[77]-1.0*(f[67]+f[66]+f[65]))+6.708203932499369*((-1.0*f[46])+f[42]+f[41]+f[40])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+60.37383539249431*((-1.0*(4.0*f[104]+5.0*(f[99]+f[98])))+4.0*(f[95]+f[94])+5.0*f[89])+9.0*(5.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(10.0*(f[60]+f[59])-20.12461179749811*f[53]))+33.54101966249684*(5.0*f[39]-4.0*f[50])+2.0*(3.0*(3.0*(15.0*(f[29]+f[28]-1.0*f[22])+11.18033988749895*((-1.0*f[15])+f[10]+f[9]))-25.0*f[4])-67.08203932499369*(f[38]+f[37]))))/(2.23606797749979*(60.37383539249431*f[93]+9.0*((-1.0*(4.0*f[80]+5.0*(f[73]+f[72])))+4.0*(f[69]+f[68])+5.0*f[58])+6.708203932499369*(4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35])))+2.0*(13.41640786499874*(f[32]+f[31])-27.0*f[25])+5.0*(5.0*f[18]-4.0*f[20]))+2.0*((-22.3606797749979*(f[17]+f[16]))+3.0*(15.0*(f[13]+f[12]-1.0*f[6])+11.18033988749895*((-1.0*f[5])+f[2]+f[1]))-25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*((-1.0*(4.0*f[80]+5.0*(f[73]+f[72])))+4.0*(f[69]+f[68])+5.0*f[58])+6.708203932499369*(4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35])))+2.0*(13.41640786499874*(f[32]+f[31])-27.0*f[25])+5.0*(5.0*f[18]-4.0*f[20]))+2.0*((-22.3606797749979*(f[17]+f[16]))+3.0*(15.0*(f[13]+f[12]-1.0*f[6])+11.18033988749895*((-1.0*f[5])+f[2]+f[1]))-25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*((-1.0*(4.0*f[80]+5.0*(f[73]+f[72])))+4.0*(f[69]+f[68])+5.0*f[58])+6.708203932499369*(4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35])))+2.0*(13.41640786499874*(f[32]+f[31])-27.0*f[25])+5.0*(5.0*f[18]-4.0*f[20]))+2.0*((-22.3606797749979*(f[17]+f[16]))+3.0*(15.0*(f[13]+f[12]-1.0*f[6])+11.18033988749895*((-1.0*f[5])+f[2]+f[1]))-25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.001666666666666667*(405.0*f[109]+60.37383539249431*((-1.0*(4.0*f[104]+5.0*(f[99]+f[98])))+4.0*(f[95]+f[94])+5.0*f[89])+9.0*(5.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(10.0*(f[60]+f[59])-20.12461179749811*f[53]))+33.54101966249684*(5.0*f[39]-4.0*f[50])+2.0*(3.0*(3.0*(15.0*(f[29]+f[28]-1.0*f[22])+11.18033988749895*((-1.0*f[15])+f[10]+f[9]))-25.0*f[4])-67.08203932499369*(f[38]+f[37]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(60.37383539249431*f[100]+5.0*(9.0*(f[65]-1.0*(f[78]+f[77]))+6.708203932499369*(f[46]-1.0*(f[41]+f[40]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*((-1.0*(4.0*f[80]+5.0*(f[73]+f[72])))+4.0*(f[69]+f[68])+5.0*f[58])+6.708203932499369*(4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35])))+2.0*(13.41640786499874*(f[32]+f[31])-27.0*f[25])+5.0*(5.0*f[18]-4.0*f[20]))+2.0*((-22.3606797749979*(f[17]+f[16]))+3.0*(15.0*(f[13]+f[12]-1.0*f[6])+11.18033988749895*((-1.0*f[5])+f[2]+f[1]))-25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.001666666666666667*(405.0*f[109]+60.37383539249431*((-1.0*(4.0*f[104]+5.0*(f[99]+f[98])))+4.0*(f[95]+f[94])+5.0*f[89])+9.0*(5.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(10.0*(f[60]+f[59])-20.12461179749811*f[53]))+33.54101966249684*(5.0*f[39]-4.0*f[50])+2.0*(3.0*(3.0*(15.0*(f[29]+f[28]-1.0*f[22])+11.18033988749895*((-1.0*f[15])+f[10]+f[9]))-25.0*f[4])-67.08203932499369*(f[38]+f[37]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(60.37383539249431*f[100]+5.0*(9.0*(f[65]-1.0*(f[78]+f[77]))+6.708203932499369*(f[46]-1.0*(f[41]+f[40]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]-1.0*(f[108]+f[107]))+60.37383539249431*((-1.0*(f[106]+f[105]))+f[104]-1.0*(f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]-1.0*f[89]+f[88]+f[87])))+9.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])+f[60]+f[59]))+20.12461179749811*((-1.0*(f[55]+f[54]))+f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]-1.0*(f[29]+f[28]+f[24]+f[23]-1.0*f[22]))+11.18033988749895*(f[15]+f[11]-1.0*(f[10]+f[9])))+25.0*f[4])))/(269.9999999999999*(f[103]+f[93]-1.0*(f[92]+f[91]))+9.0*(22.3606797749979*((-1.0*(f[82]+f[81]))+f[80]-1.0*(f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]-1.0*f[58]+f[57]+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31]))+27.0*((-1.0*(f[27]+f[26]))+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*(f[13]+f[12]+f[8]+f[7]-1.0*f[6]))+55.90169943749476*(f[5]+f[3]-1.0*(f[2]+f[1])))+125.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(269.9999999999999*(f[103]+f[93]-1.0*(f[92]+f[91]))+9.0*(22.3606797749979*((-1.0*(f[82]+f[81]))+f[80]-1.0*(f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]-1.0*f[58]+f[57]+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31]))+27.0*((-1.0*(f[27]+f[26]))+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*(f[13]+f[12]+f[8]+f[7]-1.0*f[6]))+55.90169943749476*(f[5]+f[3]-1.0*(f[2]+f[1])))+125.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.002*(269.9999999999999*(f[103]+f[93]-1.0*(f[92]+f[91]))+9.0*(22.3606797749979*((-1.0*(f[82]+f[81]))+f[80]-1.0*(f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]-1.0*f[58]+f[57]+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31]))+27.0*((-1.0*(f[27]+f[26]))+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*(f[13]+f[12]+f[8]+f[7]-1.0*f[6]))+55.90169943749476*(f[5]+f[3]-1.0*(f[2]+f[1])))+125.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*(f[108]+f[107]))+60.37383539249431*((-1.0*(f[106]+f[105]))+f[104]-1.0*(f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]-1.0*f[89]+f[88]+f[87])))+9.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])+f[60]+f[59]))+20.12461179749811*((-1.0*(f[55]+f[54]))+f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]-1.0*(f[29]+f[28]+f[24]+f[23]-1.0*f[22]))+11.18033988749895*(f[15]+f[11]-1.0*(f[10]+f[9])))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(81.0*f[110]+60.37383539249431*((-1.0*(f[102]+f[101]))+f[100]+f[90])+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]+f[67]+f[66]-1.0*f[65]))+6.708203932499369*(f[46]+f[42]-1.0*(f[41]+f[40]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.002*(269.9999999999999*(f[103]+f[93]-1.0*(f[92]+f[91]))+9.0*(22.3606797749979*((-1.0*(f[82]+f[81]))+f[80]-1.0*(f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]-1.0*f[58]+f[57]+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31]))+27.0*((-1.0*(f[27]+f[26]))+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*(f[13]+f[12]+f[8]+f[7]-1.0*f[6]))+55.90169943749476*(f[5]+f[3]-1.0*(f[2]+f[1])))+125.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*(f[108]+f[107]))+60.37383539249431*((-1.0*(f[106]+f[105]))+f[104]-1.0*(f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]-1.0*f[89]+f[88]+f[87])))+9.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])+f[60]+f[59]))+20.12461179749811*((-1.0*(f[55]+f[54]))+f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]-1.0*(f[29]+f[28]+f[24]+f[23]-1.0*f[22]))+11.18033988749895*(f[15]+f[11]-1.0*(f[10]+f[9])))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(81.0*f[110]+60.37383539249431*((-1.0*(f[102]+f[101]))+f[100]+f[90])+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]+f[67]+f[66]-1.0*f[65]))+6.708203932499369*(f[46]+f[42]-1.0*(f[41]+f[40]))+5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[0][0] = 0.05555555555555555*(fReflXYZMuQuad[7][0]+8.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0]+8.0*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[0][1] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][0]-4188761.0*fReflXYZMuQuad[6][0])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])-4.63256860547201e+14*fReflXYZMuQuad[5][0])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][0]-1.0*(8377522.0*fReflXYZMuQuad[1][0]+2.266096151179001e+23*fReflXYZMuQuad[0][0])))); 
  fReflXYQuad[0][2] = 0.05555555555555555*(fReflXYZMuQuad[7][1]+8.0*fReflXYZMuQuad[6][1]+fReflXYZMuQuad[5][1]+8.0*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[0][3] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][0]+7.4121097687552e+14*fReflXYZMuQuad[6][0]+4.63256860547201e+14*fReflXYZMuQuad[5][0])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][0]+7.4121097687552e+14*fReflXYZMuQuad[1][0]+4.63256860547201e+14*fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[0][4] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][1]-4188761.0*fReflXYZMuQuad[6][1])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])-4.63256860547201e+14*fReflXYZMuQuad[5][1])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][1]-1.0*(8377522.0*fReflXYZMuQuad[1][1]+2.266096151179001e+23*fReflXYZMuQuad[0][1])))); 
  fReflXYQuad[0][5] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][0]+9.0*fReflXYZMuQuad[6][0])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][0]-1.346286087882789e+17*fReflXYZMuQuad[5][0])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][0]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])))); 
  fReflXYQuad[0][6] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][1]+7.4121097687552e+14*fReflXYZMuQuad[6][1]+4.63256860547201e+14*fReflXYZMuQuad[5][1])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][1]+7.4121097687552e+14*fReflXYZMuQuad[1][1]+4.63256860547201e+14*fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[0][7] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][0]-4.17103711299121e+14*fReflXYZMuQuad[6][0])+5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[5][0]+2384663.0*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])+3.190559553141742e+22*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[0][8] = 0.05555555555555555*(fReflXYZMuQuad[7][2]+8.0*fReflXYZMuQuad[6][2]+fReflXYZMuQuad[5][2]+8.0*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[0][9] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][0]+1.6692641e+7*fReflXYZMuQuad[6][0])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][0]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[2][0]+2384663.0*fReflXYZMuQuad[1][0]+3.190559553141742e+22*fReflXYZMuQuad[0][0])-3.20880527843592e+30*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])))); 
  fReflXYQuad[0][10] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][1]+9.0*fReflXYZMuQuad[6][1])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][1]-1.346286087882789e+17*fReflXYZMuQuad[5][1])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][1]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])))); 
  fReflXYQuad[0][11] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][1]-4.17103711299121e+14*fReflXYZMuQuad[6][1])+5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[5][1]+2384663.0*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])+3.190559553141743e+22*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[0][12] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][2]-4188761.0*fReflXYZMuQuad[6][2])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYZMuQuad[4][2]-4.63256860547201e+14*fReflXYZMuQuad[5][2])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][2]-1.0*(8377522.0*fReflXYZMuQuad[1][2]+2.266096151179001e+23*fReflXYZMuQuad[0][2])))); 
  fReflXYQuad[0][13] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][0]-2.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0])-27.0*fReflXYZMuQuad[3][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[2][0])+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[0][14] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][2]+7.4121097687552e+14*fReflXYZMuQuad[6][2]+4.63256860547201e+14*fReflXYZMuQuad[5][2])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][2]+7.4121097687552e+14*fReflXYZMuQuad[1][2]+4.63256860547201e+14*fReflXYZMuQuad[0][2]))); 
  fReflXYQuad[0][15] = 1.719407810605221e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][0]+27.0*fReflXYZMuQuad[6][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[5][0])+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[4][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[0][16] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][1]+1.6692641e+7*fReflXYZMuQuad[6][1])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][1]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[2][1]+2384663.0*fReflXYZMuQuad[1][1]+3.190559553141743e+22*fReflXYZMuQuad[0][1])-3.20880527843592e+30*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])))); 
  fReflXYQuad[0][17] = 1.719407810605222e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][1]+fReflXYZMuQuad[5][1])-27.0*fReflXYZMuQuad[3][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[2][1])+10.0*fReflXYZMuQuad[1][1]-5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[0][18] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][2]+9.0*fReflXYZMuQuad[6][2])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][2]-1.346286087882789e+17*fReflXYZMuQuad[5][2])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][2]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])))); 
  fReflXYQuad[0][19] = 1.719407810605222e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][1]+27.0*fReflXYZMuQuad[6][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[5][1])+2.0*(5.0*fReflXYZMuQuad[3][1]-5.0*fReflXYZMuQuad[4][1])+5.0*fReflXYZMuQuad[2][1])); 
  } 

 
// node (x,y)_2 
  vcutSq_i = -(0.05*q_*(zVal*((63.63961030678928*phiWall[16]-63.63961030678928*phi[16]-47.43416490252571*phiWall[9]+47.43416490252571*phi[9])*zVal-36.74234614174767*phiWall[17]+36.74234614174767*phi[17]-21.90890230020666*phiWall[14]+21.90890230020666*phi[14]+27.38612787525831*phiWall[13]-27.38612787525831*phi[13]+32.86335345030997*phiWall[6]-32.86335345030997*phi[6]-24.49489742783179*phiWall[3]+24.49489742783179*phi[3])-21.21320343559643*phiWall[16]+21.21320343559643*phi[16]-21.21320343559643*phiWall[11]+21.21320343559643*phi[11]+15.8113883008419*phiWall[9]-15.8113883008419*phi[9]-12.64911064067352*phiWall[8]+12.64911064067352*phi[8]+15.8113883008419*phiWall[7]-15.8113883008419*phi[7]+18.97366596101028*phiWall[2]-18.97366596101028*phi[2]-14.14213562373095*phiWall[0]+14.14213562373095*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[1][0] = 0.0; 
  fReflXYQuad[1][1] = 0.0; 
  fReflXYQuad[1][2] = 0.0; 
  fReflXYQuad[1][3] = 0.0; 
  fReflXYQuad[1][4] = 0.0; 
  fReflXYQuad[1][5] = 0.0; 
  fReflXYQuad[1][6] = 0.0; 
  fReflXYQuad[1][7] = 0.0; 
  fReflXYQuad[1][8] = 0.0; 
  fReflXYQuad[1][9] = 0.0; 
  fReflXYQuad[1][10] = 0.0; 
  fReflXYQuad[1][11] = 0.0; 
  fReflXYQuad[1][12] = 0.0; 
  fReflXYQuad[1][13] = 0.0; 
  fReflXYQuad[1][14] = 0.0; 
  fReflXYQuad[1][15] = 0.0; 
  fReflXYQuad[1][16] = 0.0; 
  fReflXYQuad[1][17] = 0.0; 
  fReflXYQuad[1][18] = 0.0; 
  fReflXYQuad[1][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[1][0] = 0.05*(2.23606797749979*(6.708203932499369*f[31]+4.0*f[17]-5.0*f[16])+2.0*(5.0*f[0]-6.708203932499369*f[2])); 
  fReflXYQuad[1][1] = 0.01666666666666667*(45.0*f[56]+6.708203932499369*(4.0*f[34]-5.0*f[33])+6.0*(5.0*f[3]-6.708203932499369*f[8])); 
  fReflXYQuad[1][2] = 0.01666666666666667*(45.0*f[59]+6.708203932499369*(4.0*f[38]-5.0*f[37])+6.0*(5.0*f[4]-6.708203932499369*f[10])); 
  fReflXYQuad[1][3] = 0.01666666666666667*(45.0*f[68]+6.708203932499369*(4.0*f[44]-5.0*f[43])+6.0*(5.0*f[5]-6.708203932499369*f[13])); 
  fReflXYQuad[1][4] = 0.05*(2.23606797749979*(6.708203932499369*f[87]+4.0*f[62]-5.0*f[61])+2.0*(5.0*f[11]-6.708203932499369*f[24])); 
  fReflXYQuad[1][5] = 0.05*(2.23606797749979*(6.708203932499369*f[91]+4.0*f[71]-5.0*f[70])+2.0*(5.0*f[14]-6.708203932499369*f[27])); 
  fReflXYQuad[1][6] = 0.05*(2.23606797749979*(6.708203932499369*f[94]+4.0*f[75]-5.0*f[74])+2.0*(5.0*f[15]-6.708203932499369*f[29])); 
  fReflXYQuad[1][7] = -0.1*(6.708203932499369*f[36]-5.0*f[18]); 
  fReflXYQuad[1][8] = -0.1*(6.708203932499369*f[41]-5.0*f[19]); 
  fReflXYQuad[1][9] = -0.1*(6.708203932499369*f[48]-5.0*f[20]); 
  fReflXYQuad[1][10] = 0.01666666666666667*(45.0*f[107]+6.708203932499369*(4.0*f[97]-5.0*f[96])+6.0*(5.0*f[30]-6.708203932499369*f[55])); 
  fReflXYQuad[1][11] = -0.1*(6.708203932499369*f[64]-5.0*f[39]); 
  fReflXYQuad[1][12] = -0.1*(6.708203932499369*f[67]-5.0*f[42]); 
  fReflXYQuad[1][13] = -0.1*(6.708203932499369*f[73]-5.0*f[45]); 
  fReflXYQuad[1][14] = -0.1*(6.708203932499369*f[78]-5.0*f[46]); 
  fReflXYQuad[1][15] = -0.1*(6.708203932499369*f[82]-5.0*f[49]); 
  fReflXYQuad[1][16] = -0.1*(6.708203932499369*f[84]-5.0*f[50]); 
  fReflXYQuad[1][17] = -0.1*(6.708203932499369*f[99]-5.0*f[76]); 
  fReflXYQuad[1][18] = -0.1*(6.708203932499369*f[102]-5.0*f[79]); 
  fReflXYQuad[1][19] = -0.1*(6.708203932499369*f[106]-5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(405.0*f[107]+60.37383539249431*(4.0*(f[106]+f[99]+f[97])-5.0*(f[96]+f[94]+f[87]))+9.0*(5.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-40.24922359499622*f[55])+33.54101966249684*(4.0*(f[50]+f[39]+f[38])-5.0*f[37])+6.0*(3.0*(15.0*(f[30]+f[29]+f[24])-11.18033988749895*(f[15]+f[11]+f[10]))+25.0*f[4])))/(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]+f[73]+f[71])-5.0*(f[70]+f[68]+f[56]))+6.708203932499369*((-4.0*(f[49]+f[48]+f[45]+f[44]))+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-54.0*f[27]+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]))+2.0*(3.0*(15.0*(f[14]+f[13]+f[8])-11.18033988749895*(f[5]+f[3]+f[2]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]+f[73]+f[71])-5.0*(f[70]+f[68]+f[56]))+6.708203932499369*((-4.0*(f[49]+f[48]+f[45]+f[44]))+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-54.0*f[27]+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]))+2.0*(3.0*(15.0*(f[14]+f[13]+f[8])-11.18033988749895*(f[5]+f[3]+f[2]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]+f[73]+f[71])-5.0*(f[70]+f[68]+f[56]))+6.708203932499369*((-4.0*(f[49]+f[48]+f[45]+f[44]))+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-54.0*f[27]+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]))+2.0*(3.0*(15.0*(f[14]+f[13]+f[8])-11.18033988749895*(f[5]+f[3]+f[2]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*(f[106]+f[99]+f[97])-5.0*(f[96]+f[94]+f[87]))+9.0*(5.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-40.24922359499622*f[55])+33.54101966249684*(4.0*(f[50]+f[39]+f[38])-5.0*f[37])+6.0*(3.0*(15.0*(f[30]+f[29]+f[24])-11.18033988749895*(f[15]+f[11]+f[10]))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(60.37383539249431*f[102]+5.0*((-9.0*(f[79]+f[78]+f[67]))+6.708203932499369*(f[46]+f[42]+f[41])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]+f[73]+f[71])-5.0*(f[70]+f[68]+f[56]))+6.708203932499369*((-4.0*(f[49]+f[48]+f[45]+f[44]))+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-54.0*f[27]+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]))+2.0*(3.0*(15.0*(f[14]+f[13]+f[8])-11.18033988749895*(f[5]+f[3]+f[2]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*(f[106]+f[99]+f[97])-5.0*(f[96]+f[94]+f[87]))+9.0*(5.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-40.24922359499622*f[55])+33.54101966249684*(4.0*(f[50]+f[39]+f[38])-5.0*f[37])+6.0*(3.0*(15.0*(f[30]+f[29]+f[24])-11.18033988749895*(f[15]+f[11]+f[10]))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(60.37383539249431*f[102]+5.0*((-9.0*(f[79]+f[78]+f[67]))+6.708203932499369*(f[46]+f[42]+f[41])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[99]+f[94])+9.0*(4.0*f[84]-5.0*f[76]+4.0*f[75]-5.0*(f[74]+f[64]+f[59]))+6.708203932499369*((-4.0*f[50])+5.0*f[39]-4.0*f[38]+5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[15]+f[10])-3.0*f[29])-5.0*f[4])))/(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*(4.0*f[48]-5.0*f[45]+4.0*f[44])+5.0*(5.0*f[16]-1.0*(6.708203932499369*(f[43]+f[36]+f[31])+4.0*f[20]-5.0*f[18]+4.0*f[17])))+2.0*(3.0*(11.18033988749895*(f[5]+f[2])-15.0*f[13])-25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*(4.0*f[48]-5.0*f[45]+4.0*f[44])+5.0*(5.0*f[16]-1.0*(6.708203932499369*(f[43]+f[36]+f[31])+4.0*f[20]-5.0*f[18]+4.0*f[17])))+2.0*(3.0*(11.18033988749895*(f[5]+f[2])-15.0*f[13])-25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*(4.0*f[48]-5.0*f[45]+4.0*f[44])+5.0*(5.0*f[16]-1.0*(6.708203932499369*(f[43]+f[36]+f[31])+4.0*f[20]-5.0*f[18]+4.0*f[17])))+2.0*(3.0*(11.18033988749895*(f[5]+f[2])-15.0*f[13])-25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.008333333333333333*(60.37383539249431*(f[99]+f[94])+9.0*(4.0*f[84]-5.0*f[76]+4.0*f[75]-5.0*(f[74]+f[64]+f[59]))+6.708203932499369*((-4.0*f[50])+5.0*f[39]-4.0*f[38]+5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[15]+f[10])-3.0*f[29])-5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.05*(9.0*f[78]-6.708203932499369*(f[46]+f[41])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*(4.0*f[48]-5.0*f[45]+4.0*f[44])+5.0*(5.0*f[16]-1.0*(6.708203932499369*(f[43]+f[36]+f[31])+4.0*f[20]-5.0*f[18]+4.0*f[17])))+2.0*(3.0*(11.18033988749895*(f[5]+f[2])-15.0*f[13])-25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.008333333333333333*(60.37383539249431*(f[99]+f[94])+9.0*(4.0*f[84]-5.0*f[76]+4.0*f[75]-5.0*(f[74]+f[64]+f[59]))+6.708203932499369*((-4.0*f[50])+5.0*f[39]-4.0*f[38]+5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[15]+f[10])-3.0*f[29])-5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.05*(9.0*f[78]-6.708203932499369*(f[46]+f[41])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]+60.37383539249431*(4.0*(f[106]-1.0*f[99]+f[97])+5.0*((-1.0*f[96])+f[94]-1.0*f[87]))+9.0*(5.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*f[61])-1.0*(25.0*f[59]+40.24922359499622*f[55]))+33.54101966249684*(5.0*f[37]-4.0*(f[50]+f[39]+f[38]))+6.0*(3.0*(15.0*(f[30]-1.0*f[29]+f[24])+11.18033988749895*(f[15]-1.0*f[11]+f[10]))-25.0*f[4])))/(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]-1.0*f[73]+f[71])-5.0*f[70]+5.0*f[68])+6.708203932499369*(4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*f[33])-1.0*(33.54101966249684*f[31]+54.0*f[27]-5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17]))))+2.0*(3.0*(15.0*(f[14]-1.0*f[13]+f[8])+11.18033988749895*(f[5]-1.0*f[3]+f[2]))-25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]-1.0*f[73]+f[71])-5.0*f[70]+5.0*f[68])+6.708203932499369*(4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*f[33])-1.0*(33.54101966249684*f[31]+54.0*f[27]-5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17]))))+2.0*(3.0*(15.0*(f[14]-1.0*f[13]+f[8])+11.18033988749895*(f[5]-1.0*f[3]+f[2]))-25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]-1.0*f[73]+f[71])-5.0*f[70]+5.0*f[68])+6.708203932499369*(4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*f[33])-1.0*(33.54101966249684*f[31]+54.0*f[27]-5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17]))))+2.0*(3.0*(15.0*(f[14]-1.0*f[13]+f[8])+11.18033988749895*(f[5]-1.0*f[3]+f[2]))-25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*(f[106]-1.0*f[99]+f[97])+5.0*((-1.0*f[96])+f[94]-1.0*f[87]))+9.0*(5.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*f[61])-1.0*(25.0*f[59]+40.24922359499622*f[55]))+33.54101966249684*(5.0*f[37]-4.0*(f[50]+f[39]+f[38]))+6.0*(3.0*(15.0*(f[30]-1.0*f[29]+f[24])+11.18033988749895*(f[15]-1.0*f[11]+f[10]))-25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(60.37383539249431*f[102]+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[67])+6.708203932499369*((-1.0*f[46])+f[42]-1.0*f[41])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]-1.0*f[73]+f[71])-5.0*f[70]+5.0*f[68])+6.708203932499369*(4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*f[33])-1.0*(33.54101966249684*f[31]+54.0*f[27]-5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17]))))+2.0*(3.0*(15.0*(f[14]-1.0*f[13]+f[8])+11.18033988749895*(f[5]-1.0*f[3]+f[2]))-25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*(f[106]-1.0*f[99]+f[97])+5.0*((-1.0*f[96])+f[94]-1.0*f[87]))+9.0*(5.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*f[61])-1.0*(25.0*f[59]+40.24922359499622*f[55]))+33.54101966249684*(5.0*f[37]-4.0*(f[50]+f[39]+f[38]))+6.0*(3.0*(15.0*(f[30]-1.0*f[29]+f[24])+11.18033988749895*(f[15]-1.0*f[11]+f[10]))-25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(60.37383539249431*f[102]+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[67])+6.708203932499369*((-1.0*f[46])+f[42]-1.0*f[41])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[106]+f[87])+9.0*((-5.0*(f[85]+f[84]))+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+6.708203932499369*(5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[11]+f[10])-3.0*f[24])-5.0*f[4])))/(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*(4.0*(f[36]+f[34])-5.0*(f[49]+f[48]))+5.0*((-6.708203932499369*(f[33]+f[31]))+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[3]+f[2])-15.0*f[8])-25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*(4.0*(f[36]+f[34])-5.0*(f[49]+f[48]))+5.0*((-6.708203932499369*(f[33]+f[31]))+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[3]+f[2])-15.0*f[8])-25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*(4.0*(f[36]+f[34])-5.0*(f[49]+f[48]))+5.0*((-6.708203932499369*(f[33]+f[31]))+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[3]+f[2])-15.0*f[8])-25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.008333333333333333*(60.37383539249431*(f[106]+f[87])+9.0*((-5.0*(f[85]+f[84]))+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+6.708203932499369*(5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[11]+f[10])-3.0*f[24])-5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.05*(9.0*f[67]-6.708203932499369*(f[42]+f[41])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*(4.0*(f[36]+f[34])-5.0*(f[49]+f[48]))+5.0*((-6.708203932499369*(f[33]+f[31]))+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[3]+f[2])-15.0*f[8])-25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.008333333333333333*(60.37383539249431*(f[106]+f[87])+9.0*((-5.0*(f[85]+f[84]))+4.0*(f[64]+f[62])-5.0*(f[61]+f[59]))+6.708203932499369*(5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[11]+f[10])-3.0*f[24])-5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.05*(9.0*f[67]-6.708203932499369*(f[42]+f[41])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[106]+f[87])+9.0*(5.0*(f[84]-1.0*f[85])+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+6.708203932499369*((-5.0*f[50])+4.0*(f[39]+f[38])-5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[11]-1.0*f[10])-3.0*f[24])+5.0*f[4])))/(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*((-5.0*f[49])+5.0*f[48]+4.0*(f[34]-1.0*f[36]))+5.0*(6.708203932499369*(f[31]-1.0*f[33])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[3]-1.0*f[2])-15.0*f[8])+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*((-5.0*f[49])+5.0*f[48]+4.0*(f[34]-1.0*f[36]))+5.0*(6.708203932499369*(f[31]-1.0*f[33])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[3]-1.0*f[2])-15.0*f[8])+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*((-5.0*f[49])+5.0*f[48]+4.0*(f[34]-1.0*f[36]))+5.0*(6.708203932499369*(f[31]-1.0*f[33])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[3]-1.0*f[2])-15.0*f[8])+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(60.37383539249431*(f[106]+f[87])+9.0*(5.0*(f[84]-1.0*f[85])+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+6.708203932499369*((-5.0*f[50])+4.0*(f[39]+f[38])-5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[11]-1.0*f[10])-3.0*f[24])+5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(9.0*f[67]+6.708203932499369*(f[41]-1.0*f[42])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*((-5.0*f[49])+5.0*f[48]+4.0*(f[34]-1.0*f[36]))+5.0*(6.708203932499369*(f[31]-1.0*f[33])-5.0*f[20]+4.0*(f[18]+f[17])-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[3]-1.0*f[2])-15.0*f[8])+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(60.37383539249431*(f[106]+f[87])+9.0*(5.0*(f[84]-1.0*f[85])+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))+6.708203932499369*((-5.0*f[50])+4.0*(f[39]+f[38])-5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[11]-1.0*f[10])-3.0*f[24])+5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(9.0*f[67]+6.708203932499369*(f[41]-1.0*f[42])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]+60.37383539249431*(4.0*((-1.0*f[106])+f[99]+f[97])+5.0*(f[87]-1.0*(f[96]+f[94])))+9.0*(5.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62]))-1.0*(25.0*(f[61]+f[59])+40.24922359499622*f[55]))+33.54101966249684*(5.0*f[37]-4.0*(f[50]+f[39]+f[38]))+6.0*(3.0*(15.0*(f[30]+f[29]-1.0*f[24])+11.18033988749895*((-1.0*f[15])+f[11]+f[10]))-25.0*f[4])))/(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*((-1.0*f[82])+f[73]+f[71])+5.0*(f[56]-1.0*(f[70]+f[68])))+6.708203932499369*(4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34]))-1.0*(33.54101966249684*(f[33]+f[31])+54.0*f[27]-5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17]))))+2.0*(3.0*(15.0*(f[14]+f[13]-1.0*f[8])+11.18033988749895*((-1.0*f[5])+f[3]+f[2]))-25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*((-1.0*f[82])+f[73]+f[71])+5.0*(f[56]-1.0*(f[70]+f[68])))+6.708203932499369*(4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34]))-1.0*(33.54101966249684*(f[33]+f[31])+54.0*f[27]-5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17]))))+2.0*(3.0*(15.0*(f[14]+f[13]-1.0*f[8])+11.18033988749895*((-1.0*f[5])+f[3]+f[2]))-25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*((-1.0*f[82])+f[73]+f[71])+5.0*(f[56]-1.0*(f[70]+f[68])))+6.708203932499369*(4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34]))-1.0*(33.54101966249684*(f[33]+f[31])+54.0*f[27]-5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17]))))+2.0*(3.0*(15.0*(f[14]+f[13]-1.0*f[8])+11.18033988749895*((-1.0*f[5])+f[3]+f[2]))-25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*((-1.0*f[106])+f[99]+f[97])+5.0*(f[87]-1.0*(f[96]+f[94])))+9.0*(5.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62]))-1.0*(25.0*(f[61]+f[59])+40.24922359499622*f[55]))+33.54101966249684*(5.0*f[37]-4.0*(f[50]+f[39]+f[38]))+6.0*(3.0*(15.0*(f[30]+f[29]-1.0*f[24])+11.18033988749895*((-1.0*f[15])+f[11]+f[10]))-25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(60.37383539249431*f[102]+5.0*(9.0*(f[67]-1.0*(f[79]+f[78]))+6.708203932499369*(f[46]-1.0*(f[42]+f[41]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*((-1.0*f[82])+f[73]+f[71])+5.0*(f[56]-1.0*(f[70]+f[68])))+6.708203932499369*(4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34]))-1.0*(33.54101966249684*(f[33]+f[31])+54.0*f[27]-5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17]))))+2.0*(3.0*(15.0*(f[14]+f[13]-1.0*f[8])+11.18033988749895*((-1.0*f[5])+f[3]+f[2]))-25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*((-1.0*f[106])+f[99]+f[97])+5.0*(f[87]-1.0*(f[96]+f[94])))+9.0*(5.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62]))-1.0*(25.0*(f[61]+f[59])+40.24922359499622*f[55]))+33.54101966249684*(5.0*f[37]-4.0*(f[50]+f[39]+f[38]))+6.0*(3.0*(15.0*(f[30]+f[29]-1.0*f[24])+11.18033988749895*((-1.0*f[15])+f[11]+f[10]))-25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(60.37383539249431*f[102]+5.0*(9.0*(f[67]-1.0*(f[79]+f[78]))+6.708203932499369*(f[46]-1.0*(f[42]+f[41]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[99]+f[94])+9.0*((-1.0*(4.0*f[84]+5.0*f[76]))+4.0*f[75]+5.0*((-1.0*f[74])+f[64]+f[59]))+6.708203932499369*(4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[15]-1.0*f[10])-3.0*f[29])+5.0*f[4])))/(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*(4.0*f[44]-1.0*(4.0*f[48]+5.0*f[45]))+5.0*(6.708203932499369*((-1.0*f[43])+f[36]+f[31])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[5]-1.0*f[2])-15.0*f[13])+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*(4.0*f[44]-1.0*(4.0*f[48]+5.0*f[45]))+5.0*(6.708203932499369*((-1.0*f[43])+f[36]+f[31])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[5]-1.0*f[2])-15.0*f[13])+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*(4.0*f[44]-1.0*(4.0*f[48]+5.0*f[45]))+5.0*(6.708203932499369*((-1.0*f[43])+f[36]+f[31])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[5]-1.0*f[2])-15.0*f[13])+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.008333333333333333*(60.37383539249431*(f[99]+f[94])+9.0*((-1.0*(4.0*f[84]+5.0*f[76]))+4.0*f[75]+5.0*((-1.0*f[74])+f[64]+f[59]))+6.708203932499369*(4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[15]-1.0*f[10])-3.0*f[29])+5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.05*(9.0*f[78]+6.708203932499369*(f[41]-1.0*f[46])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*(4.0*f[44]-1.0*(4.0*f[48]+5.0*f[45]))+5.0*(6.708203932499369*((-1.0*f[43])+f[36]+f[31])+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[5]-1.0*f[2])-15.0*f[13])+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.008333333333333333*(60.37383539249431*(f[99]+f[94])+9.0*((-1.0*(4.0*f[84]+5.0*f[76]))+4.0*f[75]+5.0*((-1.0*f[74])+f[64]+f[59]))+6.708203932499369*(4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[15]-1.0*f[10])-3.0*f[29])+5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.05*(9.0*f[78]+6.708203932499369*(f[41]-1.0*f[46])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]+60.37383539249431*(4.0*(f[97]-1.0*(f[106]+f[99]))+5.0*((-1.0*f[96])+f[94]+f[87]))+9.0*(5.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))-40.24922359499622*f[55])+33.54101966249684*(4.0*(f[50]+f[39]+f[38])-5.0*f[37])+6.0*(3.0*(15.0*(f[30]-1.0*(f[29]+f[24]))+11.18033988749895*(f[15]+f[11]-1.0*f[10]))+25.0*f[4])))/(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[71]-1.0*(f[82]+f[73]))-5.0*f[70]+5.0*(f[68]+f[56]))+6.708203932499369*(4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])-5.0*f[33]+5.0*f[31])-54.0*f[27]+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]))+2.0*(3.0*(15.0*(f[14]-1.0*(f[13]+f[8]))+11.18033988749895*(f[5]+f[3]-1.0*f[2]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[71]-1.0*(f[82]+f[73]))-5.0*f[70]+5.0*(f[68]+f[56]))+6.708203932499369*(4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])-5.0*f[33]+5.0*f[31])-54.0*f[27]+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]))+2.0*(3.0*(15.0*(f[14]-1.0*(f[13]+f[8]))+11.18033988749895*(f[5]+f[3]-1.0*f[2]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[71]-1.0*(f[82]+f[73]))-5.0*f[70]+5.0*(f[68]+f[56]))+6.708203932499369*(4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])-5.0*f[33]+5.0*f[31])-54.0*f[27]+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]))+2.0*(3.0*(15.0*(f[14]-1.0*(f[13]+f[8]))+11.18033988749895*(f[5]+f[3]-1.0*f[2]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*(f[97]-1.0*(f[106]+f[99]))+5.0*((-1.0*f[96])+f[94]+f[87]))+9.0*(5.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))-40.24922359499622*f[55])+33.54101966249684*(4.0*(f[50]+f[39]+f[38])-5.0*f[37])+6.0*(3.0*(15.0*(f[30]-1.0*(f[29]+f[24]))+11.18033988749895*(f[15]+f[11]-1.0*f[10]))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(60.37383539249431*f[102]+5.0*(9.0*((-1.0*f[79])+f[78]+f[67])+6.708203932499369*(f[41]-1.0*(f[46]+f[42]))-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[71]-1.0*(f[82]+f[73]))-5.0*f[70]+5.0*(f[68]+f[56]))+6.708203932499369*(4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])-5.0*f[33]+5.0*f[31])-54.0*f[27]+5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16]))+2.0*(3.0*(15.0*(f[14]-1.0*(f[13]+f[8]))+11.18033988749895*(f[5]+f[3]-1.0*f[2]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*(f[97]-1.0*(f[106]+f[99]))+5.0*((-1.0*f[96])+f[94]+f[87]))+9.0*(5.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))-40.24922359499622*f[55])+33.54101966249684*(4.0*(f[50]+f[39]+f[38])-5.0*f[37])+6.0*(3.0*(15.0*(f[30]-1.0*(f[29]+f[24]))+11.18033988749895*(f[15]+f[11]-1.0*f[10]))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(60.37383539249431*f[102]+5.0*(9.0*((-1.0*f[79])+f[78]+f[67])+6.708203932499369*(f[41]-1.0*(f[46]+f[42]))-5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[1][0] = 0.05555555555555555*(fReflXYZMuQuad[7][0]+8.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0]+8.0*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[1][1] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][0]-4188761.0*fReflXYZMuQuad[6][0])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])-4.63256860547201e+14*fReflXYZMuQuad[5][0])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][0]-1.0*(8377522.0*fReflXYZMuQuad[1][0]+2.266096151179001e+23*fReflXYZMuQuad[0][0])))); 
  fReflXYQuad[1][2] = 0.05555555555555555*(fReflXYZMuQuad[7][1]+8.0*fReflXYZMuQuad[6][1]+fReflXYZMuQuad[5][1]+8.0*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[1][3] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][0]+7.4121097687552e+14*fReflXYZMuQuad[6][0]+4.63256860547201e+14*fReflXYZMuQuad[5][0])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][0]+7.4121097687552e+14*fReflXYZMuQuad[1][0]+4.63256860547201e+14*fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[1][4] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][1]-4188761.0*fReflXYZMuQuad[6][1])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])-4.63256860547201e+14*fReflXYZMuQuad[5][1])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][1]-1.0*(8377522.0*fReflXYZMuQuad[1][1]+2.266096151179001e+23*fReflXYZMuQuad[0][1])))); 
  fReflXYQuad[1][5] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][0]+9.0*fReflXYZMuQuad[6][0])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][0]-1.346286087882789e+17*fReflXYZMuQuad[5][0])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][0]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])))); 
  fReflXYQuad[1][6] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][1]+7.4121097687552e+14*fReflXYZMuQuad[6][1]+4.63256860547201e+14*fReflXYZMuQuad[5][1])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][1]+7.4121097687552e+14*fReflXYZMuQuad[1][1]+4.63256860547201e+14*fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[1][7] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][0]-4.17103711299121e+14*fReflXYZMuQuad[6][0])+5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[5][0]+2384663.0*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])+3.190559553141742e+22*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[1][8] = 0.05555555555555555*(fReflXYZMuQuad[7][2]+8.0*fReflXYZMuQuad[6][2]+fReflXYZMuQuad[5][2]+8.0*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[1][9] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][0]+1.6692641e+7*fReflXYZMuQuad[6][0])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][0]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[2][0]+2384663.0*fReflXYZMuQuad[1][0]+3.190559553141742e+22*fReflXYZMuQuad[0][0])-3.20880527843592e+30*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])))); 
  fReflXYQuad[1][10] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][1]+9.0*fReflXYZMuQuad[6][1])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][1]-1.346286087882789e+17*fReflXYZMuQuad[5][1])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][1]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])))); 
  fReflXYQuad[1][11] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][1]-4.17103711299121e+14*fReflXYZMuQuad[6][1])+5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[5][1]+2384663.0*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])+3.190559553141743e+22*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[1][12] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][2]-4188761.0*fReflXYZMuQuad[6][2])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYZMuQuad[4][2]-4.63256860547201e+14*fReflXYZMuQuad[5][2])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][2]-1.0*(8377522.0*fReflXYZMuQuad[1][2]+2.266096151179001e+23*fReflXYZMuQuad[0][2])))); 
  fReflXYQuad[1][13] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][0]-2.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0])-27.0*fReflXYZMuQuad[3][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[2][0])+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[1][14] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][2]+7.4121097687552e+14*fReflXYZMuQuad[6][2]+4.63256860547201e+14*fReflXYZMuQuad[5][2])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][2]+7.4121097687552e+14*fReflXYZMuQuad[1][2]+4.63256860547201e+14*fReflXYZMuQuad[0][2]))); 
  fReflXYQuad[1][15] = 1.719407810605221e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][0]+27.0*fReflXYZMuQuad[6][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[5][0])+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[4][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[1][16] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][1]+1.6692641e+7*fReflXYZMuQuad[6][1])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][1]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[2][1]+2384663.0*fReflXYZMuQuad[1][1]+3.190559553141743e+22*fReflXYZMuQuad[0][1])-3.20880527843592e+30*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])))); 
  fReflXYQuad[1][17] = 1.719407810605222e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][1]+fReflXYZMuQuad[5][1])-27.0*fReflXYZMuQuad[3][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[2][1])+10.0*fReflXYZMuQuad[1][1]-5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[1][18] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][2]+9.0*fReflXYZMuQuad[6][2])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][2]-1.346286087882789e+17*fReflXYZMuQuad[5][2])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][2]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])))); 
  fReflXYQuad[1][19] = 1.719407810605222e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][1]+27.0*fReflXYZMuQuad[6][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[5][1])+2.0*(5.0*fReflXYZMuQuad[3][1]-5.0*fReflXYZMuQuad[4][1])+5.0*fReflXYZMuQuad[2][1])); 
  } 

 
// node (x,y)_3 
  vcutSq_i = -(0.01*q_*(zVal*((426.9074841227313*phiWall[19]-426.9074841227313*phi[19]+318.1980515339465*phiWall[16]-318.1980515339465*(phi[16]+phiWall[15])+318.1980515339465*phi[15]-237.1708245126285*phiWall[9]+237.1708245126285*phi[9])*zVal-146.9693845669907*phiWall[18]+146.9693845669907*(phi[18]+phiWall[17])-146.9693845669907*phi[17]-109.5445115010333*phiWall[14]+109.5445115010333*phi[14]-109.5445115010333*phiWall[13]+109.5445115010333*phi[13]+220.454076850486*phiWall[10]-220.454076850486*phi[10]+164.3167672515499*phiWall[6]-164.3167672515499*(phi[6]+phiWall[5])+164.3167672515499*phi[5]-122.4744871391589*phiWall[3]+122.4744871391589*phi[3])-142.3024947075771*phiWall[19]+142.3024947075771*phi[19]-106.0660171779822*phiWall[16]+106.0660171779822*(phi[16]+phiWall[15])-106.0660171779822*phi[15]-84.85281374238573*phiWall[12]+84.85281374238573*(phi[12]+phiWall[11])-84.85281374238573*phi[11]+79.0569415042095*phiWall[9]-79.0569415042095*phi[9]-63.24555320336762*phiWall[8]+63.24555320336762*phi[8]-63.24555320336762*phiWall[7]+63.24555320336762*phi[7]+127.2792206135786*phiWall[4]-127.2792206135786*phi[4]+94.86832980505142*phiWall[2]-94.86832980505142*(phi[2]+phiWall[1])+94.86832980505142*phi[1]-70.71067811865477*phiWall[0]+70.71067811865477*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[2][0] = 0.0; 
  fReflXYQuad[2][1] = 0.0; 
  fReflXYQuad[2][2] = 0.0; 
  fReflXYQuad[2][3] = 0.0; 
  fReflXYQuad[2][4] = 0.0; 
  fReflXYQuad[2][5] = 0.0; 
  fReflXYQuad[2][6] = 0.0; 
  fReflXYQuad[2][7] = 0.0; 
  fReflXYQuad[2][8] = 0.0; 
  fReflXYQuad[2][9] = 0.0; 
  fReflXYQuad[2][10] = 0.0; 
  fReflXYQuad[2][11] = 0.0; 
  fReflXYQuad[2][12] = 0.0; 
  fReflXYQuad[2][13] = 0.0; 
  fReflXYQuad[2][14] = 0.0; 
  fReflXYQuad[2][15] = 0.0; 
  fReflXYQuad[2][16] = 0.0; 
  fReflXYQuad[2][17] = 0.0; 
  fReflXYQuad[2][18] = 0.0; 
  fReflXYQuad[2][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[2][0] = 0.02*(4.47213595499958*(6.708203932499369*(f[32]-1.0*f[31])+5.0*(f[17]+f[16]))+3.0*(11.18033988749895*(f[1]-1.0*f[2])-15.0*f[6])+25.0*f[0]); 
  fReflXYQuad[2][1] = 0.03333333333333333*(2.0*(9.0*(f[57]-1.0*f[56])+6.708203932499369*(f[34]+f[33]))+3.0*(3.0*(2.23606797749979*(f[7]-1.0*f[8])-3.0*f[21])+5.0*f[3])); 
  fReflXYQuad[2][2] = 0.03333333333333333*(2.0*(9.0*(f[60]-1.0*f[59])+6.708203932499369*(f[38]+f[37]))+3.0*(3.0*(2.23606797749979*(f[9]-1.0*f[10])-3.0*f[22])+5.0*f[4])); 
  fReflXYQuad[2][3] = 0.03333333333333333*(2.0*(9.0*(f[69]-1.0*f[68])+6.708203932499369*(f[44]+f[43]))+3.0*(3.0*(2.23606797749979*(f[12]-1.0*f[13])-3.0*f[25])+5.0*f[5])); 
  fReflXYQuad[2][4] = 0.02*(4.47213595499958*(6.708203932499369*(f[88]-1.0*f[87])+5.0*(f[62]+f[61]))+3.0*(11.18033988749895*(f[23]-1.0*f[24])-15.0*f[51])+25.0*f[11]); 
  fReflXYQuad[2][5] = 0.02*(4.47213595499958*(6.708203932499369*(f[92]-1.0*f[91])+5.0*(f[71]+f[70]))+3.0*(11.18033988749895*(f[26]-1.0*f[27])-15.0*f[52])+25.0*f[14]); 
  fReflXYQuad[2][6] = 0.02*(4.47213595499958*(6.708203932499369*(f[95]-1.0*f[94])+5.0*(f[75]+f[74]))+3.0*(11.18033988749895*(f[28]-1.0*f[29])-15.0*f[53])+25.0*f[15]); 
  fReflXYQuad[2][7] = -0.1*(9.0*f[58]+6.708203932499369*f[36]-1.0*(6.708203932499369*f[35]+5.0*f[18])); 
  fReflXYQuad[2][8] = -0.1*(9.0*f[65]+6.708203932499369*f[41]-1.0*(6.708203932499369*f[40]+5.0*f[19])); 
  fReflXYQuad[2][9] = -0.1*(9.0*f[80]+6.708203932499369*f[48]-1.0*(6.708203932499369*f[47]+5.0*f[20])); 
  fReflXYQuad[2][10] = 0.03333333333333333*(2.0*(9.0*(f[108]-1.0*f[107])+6.708203932499369*(f[97]+f[96]))+3.0*(3.0*(2.23606797749979*(f[54]-1.0*f[55])-3.0*f[86])+5.0*f[30])); 
  fReflXYQuad[2][11] = -0.1*(9.0*f[89]+6.708203932499369*f[64]-1.0*(6.708203932499369*f[63]+5.0*f[39])); 
  fReflXYQuad[2][12] = -0.1*(9.0*f[90]+6.708203932499369*f[67]-1.0*(6.708203932499369*f[66]+5.0*f[42])); 
  fReflXYQuad[2][13] = -0.1*(9.0*f[93]+6.708203932499369*f[73]-1.0*(6.708203932499369*f[72]+5.0*f[45])); 
  fReflXYQuad[2][14] = -0.1*(9.0*f[100]+6.708203932499369*f[78]-1.0*(6.708203932499369*f[77]+5.0*f[46])); 
  fReflXYQuad[2][15] = -0.1*(9.0*f[103]+6.708203932499369*f[82]-1.0*(6.708203932499369*f[81]+5.0*f[49])); 
  fReflXYQuad[2][16] = -0.1*(9.0*f[104]+6.708203932499369*f[84]-1.0*(6.708203932499369*f[83]+5.0*f[50])); 
  fReflXYQuad[2][17] = -0.1*(9.0*f[109]+6.708203932499369*f[99]-1.0*(6.708203932499369*f[98]+5.0*f[76])); 
  fReflXYQuad[2][18] = -0.1*(9.0*f[110]+6.708203932499369*f[102]-1.0*(6.708203932499369*f[101]+5.0*f[79])); 
  fReflXYQuad[2][19] = -0.1*(9.0*f[111]+6.708203932499369*f[106]-1.0*(6.708203932499369*f[105]+5.0*f[85])); 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]+f[108]-1.0*f[107])+60.37383539249431*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87])))+9.0*((-27.0*f[86])+10.0*((-1.0*(f[85]+f[84]))+f[83]-1.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59]))+20.12461179749811*((-1.0*f[55])+f[54]+f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]+f[29]-1.0*f[28]+f[24])-1.0*(15.0*(f[23]+f[22])+11.18033988749895*(f[15]+f[11]+f[10]-1.0*f[9])))+25.0*f[4])))/(269.9999999999999*(f[103]+f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))-45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+27.0*((-1.0*f[27])+f[26]+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]+f[13]-1.0*f[12]+f[8])-1.0*(75.0*(f[7]+f[6])+55.90169943749476*(f[5]+f[3]+f[2]-1.0*f[1])))+125.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(269.9999999999999*(f[103]+f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))-45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+27.0*((-1.0*f[27])+f[26]+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]+f[13]-1.0*f[12]+f[8])-1.0*(75.0*(f[7]+f[6])+55.90169943749476*(f[5]+f[3]+f[2]-1.0*f[1])))+125.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.002*(269.9999999999999*(f[103]+f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))-45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+27.0*((-1.0*f[27])+f[26]+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]+f[13]-1.0*f[12]+f[8])-1.0*(75.0*(f[7]+f[6])+55.90169943749476*(f[5]+f[3]+f[2]-1.0*f[1])))+125.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]-1.0*f[107])+60.37383539249431*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87])))+9.0*((-27.0*f[86])+10.0*((-1.0*(f[85]+f[84]))+f[83]-1.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59]))+20.12461179749811*((-1.0*f[55])+f[54]+f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]+f[29]-1.0*f[28]+f[24])-1.0*(15.0*(f[23]+f[22])+11.18033988749895*(f[15]+f[11]+f[10]-1.0*f[9])))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(81.0*f[110]+60.37383539249431*(f[102]-1.0*(f[101]+f[100]+f[90]))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]-1.0*f[67]+f[66]+f[65])+6.708203932499369*(f[46]+f[42]+f[41])-1.0*(6.708203932499369*f[40]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.002*(269.9999999999999*(f[103]+f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*(f[82]-1.0*(f[81]+f[80]-1.0*f[73]+f[72])+f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57]-1.0*f[56]))-45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+27.0*((-1.0*f[27])+f[26]+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]+f[13]-1.0*f[12]+f[8])-1.0*(75.0*(f[7]+f[6])+55.90169943749476*(f[5]+f[3]+f[2]-1.0*f[1])))+125.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]-1.0*f[107])+60.37383539249431*(f[106]-1.0*(f[105]+f[104]-1.0*f[99]+f[98])+f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88]-1.0*f[87])))+9.0*((-27.0*f[86])+10.0*((-1.0*(f[85]+f[84]))+f[83]-1.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59]))+20.12461179749811*((-1.0*f[55])+f[54]+f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]+f[29]-1.0*f[28]+f[24])-1.0*(15.0*(f[23]+f[22])+11.18033988749895*(f[15]+f[11]+f[10]-1.0*f[9])))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(81.0*f[110]+60.37383539249431*(f[102]-1.0*(f[101]+f[100]+f[90]))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]-1.0*f[67]+f[66]+f[65])+6.708203932499369*(f[46]+f[42]+f[41])-1.0*(6.708203932499369*f[40]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+60.37383539249431*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*f[95]-1.0*(4.0*f[94]+5.0*f[89]))+9.0*(5.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(10.0*(f[59]-1.0*f[60])-20.12461179749811*f[53]))+33.54101966249684*(5.0*f[39]-4.0*f[50])+2.0*(3.0*(3.0*(15.0*((-1.0*f[29])+f[28]+f[22])+11.18033988749895*(f[15]+f[10]))-1.0*(33.54101966249685*f[9]+25.0*f[4]))-67.08203932499369*(f[38]+f[37]))))/(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]+5.0*f[73]+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]))+6.708203932499369*(4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45]-4.0*(f[44]+f[43]))-5.0*f[36]+5.0*f[35])+2.0*(13.41640786499874*(f[31]-1.0*f[32])-27.0*f[25])+5.0*(5.0*f[18]-4.0*f[20]))+2.0*((-22.3606797749979*(f[17]+f[16]))+3.0*(15.0*((-1.0*f[13])+f[12]+f[6])+11.18033988749895*(f[5]+f[2]))-1.0*(33.54101966249685*f[1]+25.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]+5.0*f[73]+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]))+6.708203932499369*(4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45]-4.0*(f[44]+f[43]))-5.0*f[36]+5.0*f[35])+2.0*(13.41640786499874*(f[31]-1.0*f[32])-27.0*f[25])+5.0*(5.0*f[18]-4.0*f[20]))+2.0*((-22.3606797749979*(f[17]+f[16]))+3.0*(15.0*((-1.0*f[13])+f[12]+f[6])+11.18033988749895*(f[5]+f[2]))-1.0*(33.54101966249685*f[1]+25.0*f[0]))) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]+5.0*f[73]+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]))+6.708203932499369*(4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45]-4.0*(f[44]+f[43]))-5.0*f[36]+5.0*f[35])+2.0*(13.41640786499874*(f[31]-1.0*f[32])-27.0*f[25])+5.0*(5.0*f[18]-4.0*f[20]))+2.0*((-22.3606797749979*(f[17]+f[16]))+3.0*(15.0*((-1.0*f[13])+f[12]+f[6])+11.18033988749895*(f[5]+f[2]))-1.0*(33.54101966249685*f[1]+25.0*f[0]))))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.001666666666666667*(405.0*f[109]+60.37383539249431*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*f[95]-1.0*(4.0*f[94]+5.0*f[89]))+9.0*(5.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(10.0*(f[59]-1.0*f[60])-20.12461179749811*f[53]))+33.54101966249684*(5.0*f[39]-4.0*f[50])+2.0*(3.0*(3.0*(15.0*((-1.0*f[29])+f[28]+f[22])+11.18033988749895*(f[15]+f[10]))-1.0*(33.54101966249685*f[9]+25.0*f[4]))-67.08203932499369*(f[38]+f[37]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01*(60.37383539249431*f[100]+5.0*(9.0*(f[78]-1.0*(f[77]+f[65]))+6.708203932499369*(f[40]-1.0*(f[46]+f[41]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]+5.0*f[73]+4.0*f[69]-1.0*(4.0*f[68]+5.0*f[58]))+6.708203932499369*(4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45]-4.0*(f[44]+f[43]))-5.0*f[36]+5.0*f[35])+2.0*(13.41640786499874*(f[31]-1.0*f[32])-27.0*f[25])+5.0*(5.0*f[18]-4.0*f[20]))+2.0*((-22.3606797749979*(f[17]+f[16]))+3.0*(15.0*((-1.0*f[13])+f[12]+f[6])+11.18033988749895*(f[5]+f[2]))-1.0*(33.54101966249685*f[1]+25.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.001666666666666667*(405.0*f[109]+60.37383539249431*(4.0*f[104]+5.0*(f[99]-1.0*f[98])+4.0*f[95]-1.0*(4.0*f[94]+5.0*f[89]))+9.0*(5.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(10.0*(f[59]-1.0*f[60])-20.12461179749811*f[53]))+33.54101966249684*(5.0*f[39]-4.0*f[50])+2.0*(3.0*(3.0*(15.0*((-1.0*f[29])+f[28]+f[22])+11.18033988749895*(f[15]+f[10]))-1.0*(33.54101966249685*f[9]+25.0*f[4]))-67.08203932499369*(f[38]+f[37]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01*(60.37383539249431*f[100]+5.0*(9.0*(f[78]-1.0*(f[77]+f[65]))+6.708203932499369*(f[40]-1.0*(f[46]+f[41]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]-1.0*f[107])+60.37383539249431*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87]))+9.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64]-1.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+20.12461179749811*((-1.0*f[55])+f[54]-1.0*f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]-1.0*f[29]+f[28]+f[24]-1.0*f[23]+f[22])+11.18033988749895*(f[15]-1.0*f[11]+f[10]))-1.0*(33.54101966249685*f[9]+25.0*f[4]))))/(269.9999999999999*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])-45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*f[13]+f[12]+f[8]-1.0*f[7]+f[6])+55.90169943749476*(f[5]-1.0*f[3]+f[2]))-1.0*(167.7050983124843*f[1]+125.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(269.9999999999999*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])-45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*f[13]+f[12]+f[8]-1.0*f[7]+f[6])+55.90169943749476*(f[5]-1.0*f[3]+f[2]))-1.0*(167.7050983124843*f[1]+125.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.002*(269.9999999999999*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])-45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*f[13]+f[12]+f[8]-1.0*f[7]+f[6])+55.90169943749476*(f[5]-1.0*f[3]+f[2]))-1.0*(167.7050983124843*f[1]+125.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]-1.0*f[107])+60.37383539249431*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87]))+9.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64]-1.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+20.12461179749811*((-1.0*f[55])+f[54]-1.0*f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]-1.0*f[29]+f[28]+f[24]-1.0*f[23]+f[22])+11.18033988749895*(f[15]-1.0*f[11]+f[10]))-1.0*(33.54101966249685*f[9]+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(81.0*f[110]+60.37383539249431*(f[102]-1.0*f[101]+f[100]-1.0*f[90])+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*(f[77]+f[67]-1.0*f[66]+f[65]))+6.708203932499369*((-1.0*f[46])+f[42]-1.0*f[41]+f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.002*(269.9999999999999*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*(f[82]-1.0*f[81]+f[80]-1.0*f[73]+f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56])-45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*f[13]+f[12]+f[8]-1.0*f[7]+f[6])+55.90169943749476*(f[5]-1.0*f[3]+f[2]))-1.0*(167.7050983124843*f[1]+125.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]-1.0*f[107])+60.37383539249431*(f[106]-1.0*f[105]+f[104]-1.0*f[99]+f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87]))+9.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64]-1.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+20.12461179749811*((-1.0*f[55])+f[54]-1.0*f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]-1.0*f[29]+f[28]+f[24]-1.0*f[23]+f[22])+11.18033988749895*(f[15]-1.0*f[11]+f[10]))-1.0*(33.54101966249685*f[9]+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(81.0*f[110]+60.37383539249431*(f[102]-1.0*f[101]+f[100]-1.0*f[90])+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*(f[77]+f[67]-1.0*f[66]+f[65]))+6.708203932499369*((-1.0*f[46])+f[42]-1.0*f[41]+f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+60.37383539249431*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+9.0*(25.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-20.12461179749811*f[51]))+167.7050983124842*f[50]+2.0*(3.0*(3.0*(15.0*((-1.0*f[24])+f[23]+f[22])+11.18033988749895*(f[11]+f[10]))-1.0*(33.54101966249685*f[9]+25.0*f[4]))-67.08203932499369*(f[39]+f[38]+f[37]))))/(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+4.0*(f[58]+f[57]-1.0*f[56]))+6.708203932499369*(5.0*f[47]-5.0*(f[49]+f[48]))+2.0*(13.41640786499874*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])-27.0*f[21])+25.0*f[20])+2.0*((-22.3606797749979*(f[18]+f[17]+f[16]))+3.0*(15.0*((-1.0*f[8])+f[7]+f[6])+11.18033988749895*(f[3]+f[2]))-1.0*(33.54101966249685*f[1]+25.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+4.0*(f[58]+f[57]-1.0*f[56]))+6.708203932499369*(5.0*f[47]-5.0*(f[49]+f[48]))+2.0*(13.41640786499874*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])-27.0*f[21])+25.0*f[20])+2.0*((-22.3606797749979*(f[18]+f[17]+f[16]))+3.0*(15.0*((-1.0*f[8])+f[7]+f[6])+11.18033988749895*(f[3]+f[2]))-1.0*(33.54101966249685*f[1]+25.0*f[0]))) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+4.0*(f[58]+f[57]-1.0*f[56]))+6.708203932499369*(5.0*f[47]-5.0*(f[49]+f[48]))+2.0*(13.41640786499874*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])-27.0*f[21])+25.0*f[20])+2.0*((-22.3606797749979*(f[18]+f[17]+f[16]))+3.0*(15.0*((-1.0*f[8])+f[7]+f[6])+11.18033988749895*(f[3]+f[2]))-1.0*(33.54101966249685*f[1]+25.0*f[0]))))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+9.0*(25.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-20.12461179749811*f[51]))+167.7050983124842*f[50]+2.0*(3.0*(3.0*(15.0*((-1.0*f[24])+f[23]+f[22])+11.18033988749895*(f[11]+f[10]))-1.0*(33.54101966249685*f[9]+25.0*f[4]))-67.08203932499369*(f[39]+f[38]+f[37]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01*(60.37383539249431*f[90]+5.0*(9.0*(f[67]-1.0*(f[66]+f[65]))+6.708203932499369*(f[40]-1.0*(f[42]+f[41]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]-1.0*(f[81]+f[80]))+4.0*(f[58]+f[57]-1.0*f[56]))+6.708203932499369*(5.0*f[47]-5.0*(f[49]+f[48]))+2.0*(13.41640786499874*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31])-27.0*f[21])+25.0*f[20])+2.0*((-22.3606797749979*(f[18]+f[17]+f[16]))+3.0*(15.0*((-1.0*f[8])+f[7]+f[6])+11.18033988749895*(f[3]+f[2]))-1.0*(33.54101966249685*f[1]+25.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*(f[106]-1.0*(f[105]+f[104]))+4.0*(f[89]+f[88]-1.0*f[87]))+9.0*(25.0*(f[83]-1.0*(f[85]+f[84]))+2.0*(10.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59])-20.12461179749811*f[51]))+167.7050983124842*f[50]+2.0*(3.0*(3.0*(15.0*((-1.0*f[24])+f[23]+f[22])+11.18033988749895*(f[11]+f[10]))-1.0*(33.54101966249685*f[9]+25.0*f[4]))-67.08203932499369*(f[39]+f[38]+f[37]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01*(60.37383539249431*f[90]+5.0*(9.0*(f[67]-1.0*(f[66]+f[65]))+6.708203932499369*(f[40]-1.0*(f[42]+f[41]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+60.37383539249431*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+45.0*(5.0*((-1.0*f[85])+f[84]-1.0*f[83])+4.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60]))-1.0*(18.0*(10.0*f[59]+20.12461179749811*f[51])+167.7050983124842*f[50])+2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[24])+f[23]-1.0*f[22])+11.18033988749895*(f[11]-1.0*f[10]+f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]+f[80])+4.0*((-1.0*f[58])+f[57]-1.0*f[56]))+6.708203932499369*((-5.0*f[49])+5.0*f[48]+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]))-1.0*(2.0*(13.41640786499874*f[31]+27.0*f[21])+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*((-1.0*f[8])+f[7]-1.0*f[6])+11.18033988749895*(f[3]-1.0*f[2]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]+f[80])+4.0*((-1.0*f[58])+f[57]-1.0*f[56]))+6.708203932499369*((-5.0*f[49])+5.0*f[48]+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]))-1.0*(2.0*(13.41640786499874*f[31]+27.0*f[21])+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*((-1.0*f[8])+f[7]-1.0*f[6])+11.18033988749895*(f[3]-1.0*f[2]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]+f[80])+4.0*((-1.0*f[58])+f[57]-1.0*f[56]))+6.708203932499369*((-5.0*f[49])+5.0*f[48]+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]))-1.0*(2.0*(13.41640786499874*f[31]+27.0*f[21])+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*((-1.0*f[8])+f[7]-1.0*f[6])+11.18033988749895*(f[3]-1.0*f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+45.0*(5.0*((-1.0*f[85])+f[84]-1.0*f[83])+4.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60]))-1.0*(18.0*(10.0*f[59]+20.12461179749811*f[51])+167.7050983124842*f[50])+2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[24])+f[23]-1.0*f[22])+11.18033988749895*(f[11]-1.0*f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.01*(60.37383539249431*f[90]+5.0*(9.0*(f[67]-1.0*f[66]+f[65])+6.708203932499369*(f[41]-1.0*f[42])-1.0*(6.708203932499369*f[40]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]+f[80])+4.0*((-1.0*f[58])+f[57]-1.0*f[56]))+6.708203932499369*((-5.0*f[49])+5.0*f[48]+4.0*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32]))-1.0*(2.0*(13.41640786499874*f[31]+27.0*f[21])+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*((-1.0*f[8])+f[7]-1.0*f[6])+11.18033988749895*(f[3]-1.0*f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*(f[106]-1.0*f[105]+f[104])+4.0*((-1.0*f[89])+f[88]-1.0*f[87]))+45.0*(5.0*((-1.0*f[85])+f[84]-1.0*f[83])+4.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60]))-1.0*(18.0*(10.0*f[59]+20.12461179749811*f[51])+167.7050983124842*f[50])+2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[24])+f[23]-1.0*f[22])+11.18033988749895*(f[11]-1.0*f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.01*(60.37383539249431*f[90]+5.0*(9.0*(f[67]-1.0*f[66]+f[65])+6.708203932499369*(f[41]-1.0*f[42])-1.0*(6.708203932499369*f[40]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))+60.37383539249431*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87])))+9.0*(27.0*f[86]+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+20.12461179749811*(f[55]-1.0*(f[54]+f[53]-1.0*f[51])))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*(f[30]+f[29]))+f[28]+f[24]-1.0*(f[23]+f[22]))+11.18033988749895*(f[15]-1.0*(f[11]+f[10]-1.0*f[9])))+25.0*f[4])))/(269.9999999999999*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+9.0*(22.3606797749979*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31]))+27.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21]))+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*(f[14]+f[13]))+f[12]+f[8]-1.0*(f[7]+f[6]))+55.90169943749476*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))+125.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(269.9999999999999*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+9.0*(22.3606797749979*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31]))+27.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21]))+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*(f[14]+f[13]))+f[12]+f[8]-1.0*(f[7]+f[6]))+55.90169943749476*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))+125.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.002*(269.9999999999999*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+9.0*(22.3606797749979*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31]))+27.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21]))+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*(f[14]+f[13]))+f[12]+f[8]-1.0*(f[7]+f[6]))+55.90169943749476*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))+125.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))+60.37383539249431*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87])))+9.0*(27.0*f[86]+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+20.12461179749811*(f[55]-1.0*(f[54]+f[53]-1.0*f[51])))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*(f[30]+f[29]))+f[28]+f[24]-1.0*(f[23]+f[22]))+11.18033988749895*(f[15]-1.0*(f[11]+f[10]-1.0*f[9])))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(81.0*f[110]+60.37383539249431*(f[102]-1.0*(f[101]+f[100]-1.0*f[90]))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]+f[67]-1.0*(f[66]+f[65]))+6.708203932499369*(f[46]-1.0*(f[42]+f[41]-1.0*f[40]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.002*(269.9999999999999*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+9.0*(22.3606797749979*(f[82]-1.0*(f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57]-1.0*f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31]))+27.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21]))+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*(f[14]+f[13]))+f[12]+f[8]-1.0*(f[7]+f[6]))+55.90169943749476*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))+125.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))+60.37383539249431*(f[106]-1.0*(f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87])))+9.0*(27.0*f[86]+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+20.12461179749811*(f[55]-1.0*(f[54]+f[53]-1.0*f[51])))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*(f[30]+f[29]))+f[28]+f[24]-1.0*(f[23]+f[22]))+11.18033988749895*(f[15]-1.0*(f[11]+f[10]-1.0*f[9])))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(81.0*f[110]+60.37383539249431*(f[102]-1.0*(f[101]+f[100]-1.0*f[90]))+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]+f[67]-1.0*(f[66]+f[65]))+6.708203932499369*(f[46]-1.0*(f[42]+f[41]-1.0*f[40]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+60.37383539249431*((-4.0*f[104])+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94])+5.0*f[89])+9.0*(5.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(10.0*f[60]-1.0*(10.0*f[59]+20.12461179749811*f[53])))+33.54101966249684*(4.0*f[50]-5.0*f[39])+2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[29])+f[28]-1.0*f[22])+11.18033988749895*(f[15]-1.0*f[10]+f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[93]+9.0*((-4.0*f[80])+5.0*f[73]+4.0*(f[69]-1.0*f[68])+5.0*f[58])+6.708203932499369*(4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*f[36])+2.0*(13.41640786499874*f[32]-1.0*(13.41640786499874*f[31]+27.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*((-1.0*f[13])+f[12]-1.0*f[6])+11.18033988749895*(f[5]-1.0*f[2]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*((-4.0*f[80])+5.0*f[73]+4.0*(f[69]-1.0*f[68])+5.0*f[58])+6.708203932499369*(4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*f[36])+2.0*(13.41640786499874*f[32]-1.0*(13.41640786499874*f[31]+27.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*((-1.0*f[13])+f[12]-1.0*f[6])+11.18033988749895*(f[5]-1.0*f[2]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*((-4.0*f[80])+5.0*f[73]+4.0*(f[69]-1.0*f[68])+5.0*f[58])+6.708203932499369*(4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*f[36])+2.0*(13.41640786499874*f[32]-1.0*(13.41640786499874*f[31]+27.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*((-1.0*f[13])+f[12]-1.0*f[6])+11.18033988749895*(f[5]-1.0*f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.001666666666666667*(405.0*f[109]+60.37383539249431*((-4.0*f[104])+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94])+5.0*f[89])+9.0*(5.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(10.0*f[60]-1.0*(10.0*f[59]+20.12461179749811*f[53])))+33.54101966249684*(4.0*f[50]-5.0*f[39])+2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[29])+f[28]-1.0*f[22])+11.18033988749895*(f[15]-1.0*f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(60.37383539249431*f[100]+5.0*(9.0*(f[78]-1.0*f[77]+f[65])+6.708203932499369*(f[41]-1.0*f[46])-1.0*(6.708203932499369*f[40]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*((-4.0*f[80])+5.0*f[73]+4.0*(f[69]-1.0*f[68])+5.0*f[58])+6.708203932499369*(4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*f[36])+2.0*(13.41640786499874*f[32]-1.0*(13.41640786499874*f[31]+27.0*f[25]))+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*((-1.0*f[13])+f[12]-1.0*f[6])+11.18033988749895*(f[5]-1.0*f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.001666666666666667*(405.0*f[109]+60.37383539249431*((-4.0*f[104])+5.0*(f[99]-1.0*f[98])+4.0*(f[95]-1.0*f[94])+5.0*f[89])+9.0*(5.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(10.0*f[60]-1.0*(10.0*f[59]+20.12461179749811*f[53])))+33.54101966249684*(4.0*f[50]-5.0*f[39])+2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[29])+f[28]-1.0*f[22])+11.18033988749895*(f[15]-1.0*f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(60.37383539249431*f[100]+5.0*(9.0*(f[78]-1.0*f[77]+f[65])+6.708203932499369*(f[41]-1.0*f[46])-1.0*(6.708203932499369*f[40]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])+60.37383539249431*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87])))+9.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+20.12461179749811*(f[55]-1.0*f[54]+f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[30])+f[29]-1.0*f[28]+f[24]-1.0*f[23]+f[22])+11.18033988749895*(f[10]-1.0*(f[15]+f[11])))-1.0*(33.54101966249685*f[9]+25.0*f[4]))))/(269.9999999999999*(f[103]+f[93]-1.0*f[92]+f[91])+9.0*(22.3606797749979*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*(f[27]-1.0*f[26]+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6])+55.90169943749476*(f[2]-1.0*(f[5]+f[3])))-1.0*(167.7050983124843*f[1]+125.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(269.9999999999999*(f[103]+f[93]-1.0*f[92]+f[91])+9.0*(22.3606797749979*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*(f[27]-1.0*f[26]+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6])+55.90169943749476*(f[2]-1.0*(f[5]+f[3])))-1.0*(167.7050983124843*f[1]+125.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.002*(269.9999999999999*(f[103]+f[93]-1.0*f[92]+f[91])+9.0*(22.3606797749979*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*(f[27]-1.0*f[26]+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6])+55.90169943749476*(f[2]-1.0*(f[5]+f[3])))-1.0*(167.7050983124843*f[1]+125.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])+60.37383539249431*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87])))+9.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+20.12461179749811*(f[55]-1.0*f[54]+f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[30])+f[29]-1.0*f[28]+f[24]-1.0*f[23]+f[22])+11.18033988749895*(f[10]-1.0*(f[15]+f[11])))-1.0*(33.54101966249685*f[9]+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(81.0*f[110]+60.37383539249431*(f[102]-1.0*f[101]+f[100]+f[90])+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[77]+f[67]-1.0*f[66]+f[65])+6.708203932499369*(f[41]-1.0*(f[46]+f[42]))-1.0*(6.708203932499369*f[40]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.002*(269.9999999999999*(f[103]+f[93]-1.0*f[92]+f[91])+9.0*(22.3606797749979*(f[82]-1.0*f[81]+f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*(f[68]+f[58])+f[57]-1.0*f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*(f[27]-1.0*f[26]+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]-1.0*f[12]+f[8]-1.0*f[7]+f[6])+55.90169943749476*(f[2]-1.0*(f[5]+f[3])))-1.0*(167.7050983124843*f[1]+125.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])+60.37383539249431*(f[106]-1.0*f[105]+f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*(f[94]+f[89])+f[88]-1.0*f[87])))+9.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+20.12461179749811*(f[55]-1.0*f[54]+f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[30])+f[29]-1.0*f[28]+f[24]-1.0*f[23]+f[22])+11.18033988749895*(f[10]-1.0*(f[15]+f[11])))-1.0*(33.54101966249685*f[9]+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(81.0*f[110]+60.37383539249431*(f[102]-1.0*f[101]+f[100]+f[90])+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[77]+f[67]-1.0*f[66]+f[65])+6.708203932499369*(f[41]-1.0*(f[46]+f[42]))-1.0*(6.708203932499369*f[40]+5.0*f[19]))))*fac; 
   } 
  } 
  fReflXYQuad[2][0] = 0.05555555555555555*(fReflXYZMuQuad[7][0]+8.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0]+8.0*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[2][1] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][0]-4188761.0*fReflXYZMuQuad[6][0])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])-4.63256860547201e+14*fReflXYZMuQuad[5][0])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][0]-1.0*(8377522.0*fReflXYZMuQuad[1][0]+2.266096151179001e+23*fReflXYZMuQuad[0][0])))); 
  fReflXYQuad[2][2] = 0.05555555555555555*(fReflXYZMuQuad[7][1]+8.0*fReflXYZMuQuad[6][1]+fReflXYZMuQuad[5][1]+8.0*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[2][3] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][0]+7.4121097687552e+14*fReflXYZMuQuad[6][0]+4.63256860547201e+14*fReflXYZMuQuad[5][0])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][0]+7.4121097687552e+14*fReflXYZMuQuad[1][0]+4.63256860547201e+14*fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[2][4] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][1]-4188761.0*fReflXYZMuQuad[6][1])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])-4.63256860547201e+14*fReflXYZMuQuad[5][1])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][1]-1.0*(8377522.0*fReflXYZMuQuad[1][1]+2.266096151179001e+23*fReflXYZMuQuad[0][1])))); 
  fReflXYQuad[2][5] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][0]+9.0*fReflXYZMuQuad[6][0])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][0]-1.346286087882789e+17*fReflXYZMuQuad[5][0])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][0]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])))); 
  fReflXYQuad[2][6] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][1]+7.4121097687552e+14*fReflXYZMuQuad[6][1]+4.63256860547201e+14*fReflXYZMuQuad[5][1])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][1]+7.4121097687552e+14*fReflXYZMuQuad[1][1]+4.63256860547201e+14*fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[2][7] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][0]-4.17103711299121e+14*fReflXYZMuQuad[6][0])+5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[5][0]+2384663.0*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])+3.190559553141742e+22*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[2][8] = 0.05555555555555555*(fReflXYZMuQuad[7][2]+8.0*fReflXYZMuQuad[6][2]+fReflXYZMuQuad[5][2]+8.0*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[2][9] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][0]+1.6692641e+7*fReflXYZMuQuad[6][0])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][0]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[2][0]+2384663.0*fReflXYZMuQuad[1][0]+3.190559553141742e+22*fReflXYZMuQuad[0][0])-3.20880527843592e+30*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])))); 
  fReflXYQuad[2][10] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][1]+9.0*fReflXYZMuQuad[6][1])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][1]-1.346286087882789e+17*fReflXYZMuQuad[5][1])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][1]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])))); 
  fReflXYQuad[2][11] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][1]-4.17103711299121e+14*fReflXYZMuQuad[6][1])+5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[5][1]+2384663.0*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])+3.190559553141743e+22*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[2][12] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][2]-4188761.0*fReflXYZMuQuad[6][2])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYZMuQuad[4][2]-4.63256860547201e+14*fReflXYZMuQuad[5][2])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][2]-1.0*(8377522.0*fReflXYZMuQuad[1][2]+2.266096151179001e+23*fReflXYZMuQuad[0][2])))); 
  fReflXYQuad[2][13] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][0]-2.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0])-27.0*fReflXYZMuQuad[3][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[2][0])+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[2][14] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][2]+7.4121097687552e+14*fReflXYZMuQuad[6][2]+4.63256860547201e+14*fReflXYZMuQuad[5][2])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][2]+7.4121097687552e+14*fReflXYZMuQuad[1][2]+4.63256860547201e+14*fReflXYZMuQuad[0][2]))); 
  fReflXYQuad[2][15] = 1.719407810605221e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][0]+27.0*fReflXYZMuQuad[6][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[5][0])+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[4][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[2][16] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][1]+1.6692641e+7*fReflXYZMuQuad[6][1])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][1]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[2][1]+2384663.0*fReflXYZMuQuad[1][1]+3.190559553141743e+22*fReflXYZMuQuad[0][1])-3.20880527843592e+30*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])))); 
  fReflXYQuad[2][17] = 1.719407810605222e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][1]+fReflXYZMuQuad[5][1])-27.0*fReflXYZMuQuad[3][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[2][1])+10.0*fReflXYZMuQuad[1][1]-5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[2][18] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][2]+9.0*fReflXYZMuQuad[6][2])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][2]-1.346286087882789e+17*fReflXYZMuQuad[5][2])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][2]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])))); 
  fReflXYQuad[2][19] = 1.719407810605222e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][1]+27.0*fReflXYZMuQuad[6][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[5][1])+2.0*(5.0*fReflXYZMuQuad[3][1]-5.0*fReflXYZMuQuad[4][1])+5.0*fReflXYZMuQuad[2][1])); 
  } 

 
// node (x,y)_4 
  vcutSq_i = -(0.05*q_*(zVal*((63.63961030678928*phiWall[15]-63.63961030678928*phi[15]-47.43416490252571*phiWall[9]+47.43416490252571*phi[9])*zVal-36.74234614174767*phiWall[18]+36.74234614174767*phi[18]+27.38612787525831*phiWall[14]-27.38612787525831*phi[14]-21.90890230020666*phiWall[13]+21.90890230020666*phi[13]+32.86335345030997*phiWall[5]-32.86335345030997*phi[5]-24.49489742783179*phiWall[3]+24.49489742783179*phi[3])-21.21320343559643*phiWall[15]+21.21320343559643*phi[15]-21.21320343559643*phiWall[12]+21.21320343559643*phi[12]+15.8113883008419*phiWall[9]-15.8113883008419*phi[9]+15.8113883008419*phiWall[8]-15.8113883008419*phi[8]-12.64911064067352*phiWall[7]+12.64911064067352*phi[7]+18.97366596101028*phiWall[1]-18.97366596101028*phi[1]-14.14213562373095*phiWall[0]+14.14213562373095*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[3][0] = 0.0; 
  fReflXYQuad[3][1] = 0.0; 
  fReflXYQuad[3][2] = 0.0; 
  fReflXYQuad[3][3] = 0.0; 
  fReflXYQuad[3][4] = 0.0; 
  fReflXYQuad[3][5] = 0.0; 
  fReflXYQuad[3][6] = 0.0; 
  fReflXYQuad[3][7] = 0.0; 
  fReflXYQuad[3][8] = 0.0; 
  fReflXYQuad[3][9] = 0.0; 
  fReflXYQuad[3][10] = 0.0; 
  fReflXYQuad[3][11] = 0.0; 
  fReflXYQuad[3][12] = 0.0; 
  fReflXYQuad[3][13] = 0.0; 
  fReflXYQuad[3][14] = 0.0; 
  fReflXYQuad[3][15] = 0.0; 
  fReflXYQuad[3][16] = 0.0; 
  fReflXYQuad[3][17] = 0.0; 
  fReflXYQuad[3][18] = 0.0; 
  fReflXYQuad[3][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[3][0] = 0.05*(2.23606797749979*(6.708203932499369*f[32]-5.0*f[17])+2.0*(2.23606797749979*(2.0*f[16]-3.0*f[1])+5.0*f[0])); 
  fReflXYQuad[3][1] = 0.01666666666666667*(5.0*(9.0*f[57]-6.708203932499369*f[34])+2.0*(13.41640786499874*f[33]+3.0*(5.0*f[3]-6.708203932499369*f[7]))); 
  fReflXYQuad[3][2] = 0.01666666666666667*(5.0*(9.0*f[60]-6.708203932499369*f[38])+2.0*(13.41640786499874*f[37]+3.0*(5.0*f[4]-6.708203932499369*f[9]))); 
  fReflXYQuad[3][3] = 0.01666666666666667*(5.0*(9.0*f[69]-6.708203932499369*f[44])+2.0*(13.41640786499874*f[43]+3.0*(5.0*f[5]-6.708203932499369*f[12]))); 
  fReflXYQuad[3][4] = 0.05*(2.23606797749979*(6.708203932499369*f[88]-5.0*f[62])+2.0*(2.23606797749979*(2.0*f[61]-3.0*f[23])+5.0*f[11])); 
  fReflXYQuad[3][5] = 0.05*(2.23606797749979*(6.708203932499369*f[92]-5.0*f[71])+2.0*(2.23606797749979*(2.0*f[70]-3.0*f[26])+5.0*f[14])); 
  fReflXYQuad[3][6] = 0.05*(2.23606797749979*(6.708203932499369*f[95]-5.0*f[75])+2.0*(2.23606797749979*(2.0*f[74]-3.0*f[28])+5.0*f[15])); 
  fReflXYQuad[3][7] = -0.1*(6.708203932499369*f[35]-5.0*f[18]); 
  fReflXYQuad[3][8] = -0.1*(6.708203932499369*f[40]-5.0*f[19]); 
  fReflXYQuad[3][9] = -0.1*(6.708203932499369*f[47]-5.0*f[20]); 
  fReflXYQuad[3][10] = 0.01666666666666667*(5.0*(9.0*f[108]-6.708203932499369*f[97])+2.0*(13.41640786499874*f[96]+3.0*(5.0*f[30]-6.708203932499369*f[54]))); 
  fReflXYQuad[3][11] = -0.1*(6.708203932499369*f[63]-5.0*f[39]); 
  fReflXYQuad[3][12] = -0.1*(6.708203932499369*f[66]-5.0*f[42]); 
  fReflXYQuad[3][13] = -0.1*(6.708203932499369*f[72]-5.0*f[45]); 
  fReflXYQuad[3][14] = -0.1*(6.708203932499369*f[77]-5.0*f[46]); 
  fReflXYQuad[3][15] = -0.1*(6.708203932499369*f[81]-5.0*f[49]); 
  fReflXYQuad[3][16] = -0.1*(6.708203932499369*f[83]-5.0*f[50]); 
  fReflXYQuad[3][17] = -0.1*(6.708203932499369*f[98]-5.0*f[76]); 
  fReflXYQuad[3][18] = -0.1*(6.708203932499369*f[101]-5.0*f[79]); 
  fReflXYQuad[3][19] = -0.1*(6.708203932499369*f[105]-5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(405.0*f[108]+60.37383539249431*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96]-5.0*(f[95]+f[88]))+9.0*(5.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])-40.24922359499622*f[54])+33.54101966249684*(4.0*(f[50]+f[39])-5.0*f[38])+2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*(f[30]+f[28]+f[23])-11.18033988749895*(f[15]+f[11]+f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-5.0*(f[69]+f[57]))+6.708203932499369*((-4.0*(f[49]+f[47]+f[45]))+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-54.0*f[26]+5.0*(4.0*(f[20]+f[18])-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[14]+f[12]+f[7])-11.18033988749895*(f[5]+f[3]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-5.0*(f[69]+f[57]))+6.708203932499369*((-4.0*(f[49]+f[47]+f[45]))+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-54.0*f[26]+5.0*(4.0*(f[20]+f[18])-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[14]+f[12]+f[7])-11.18033988749895*(f[5]+f[3]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-5.0*(f[69]+f[57]))+6.708203932499369*((-4.0*(f[49]+f[47]+f[45]))+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-54.0*f[26]+5.0*(4.0*(f[20]+f[18])-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[14]+f[12]+f[7])-11.18033988749895*(f[5]+f[3]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.001666666666666667*(405.0*f[108]+60.37383539249431*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96]-5.0*(f[95]+f[88]))+9.0*(5.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])-40.24922359499622*f[54])+33.54101966249684*(4.0*(f[50]+f[39])-5.0*f[38])+2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*(f[30]+f[28]+f[23])-11.18033988749895*(f[15]+f[11]+f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(60.37383539249431*f[101]+5.0*((-9.0*(f[79]+f[77]+f[66]))+6.708203932499369*(f[46]+f[42]+f[40])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[81]+f[72])-5.0*f[71]+4.0*f[70]-5.0*(f[69]+f[57]))+6.708203932499369*((-4.0*(f[49]+f[47]+f[45]))+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-54.0*f[26]+5.0*(4.0*(f[20]+f[18])-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[14]+f[12]+f[7])-11.18033988749895*(f[5]+f[3]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.001666666666666667*(405.0*f[108]+60.37383539249431*(4.0*(f[105]+f[98])-5.0*f[97]+4.0*f[96]-5.0*(f[95]+f[88]))+9.0*(5.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])-40.24922359499622*f[54])+33.54101966249684*(4.0*(f[50]+f[39])-5.0*f[38])+2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*(f[30]+f[28]+f[23])-11.18033988749895*(f[15]+f[11]+f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(60.37383539249431*f[101]+5.0*((-9.0*(f[79]+f[77]+f[66]))+6.708203932499369*(f[46]+f[42]+f[40])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[98]+f[95])+9.0*(4.0*f[83]-5.0*(f[76]+f[75])+4.0*f[74]-5.0*(f[63]+f[60]))+6.708203932499369*(5.0*(f[39]+f[38])-4.0*f[50])+2.0*(3.0*(3.0*(2.23606797749979*(f[15]+f[9])-3.0*f[28])-5.0*f[4])-13.41640786499874*f[37])))/(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43])-5.0*(6.708203932499369*(f[35]+f[32])+4.0*f[20]-5.0*(f[18]+f[17])))+2.0*((-22.3606797749979*f[16])+3.0*(11.18033988749895*(f[5]+f[1])-15.0*f[12])-25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43])-5.0*(6.708203932499369*(f[35]+f[32])+4.0*f[20]-5.0*(f[18]+f[17])))+2.0*((-22.3606797749979*f[16])+3.0*(11.18033988749895*(f[5]+f[1])-15.0*f[12])-25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43])-5.0*(6.708203932499369*(f[35]+f[32])+4.0*f[20]-5.0*(f[18]+f[17])))+2.0*((-22.3606797749979*f[16])+3.0*(11.18033988749895*(f[5]+f[1])-15.0*f[12])-25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.008333333333333333*(60.37383539249431*(f[98]+f[95])+9.0*(4.0*f[83]-5.0*(f[76]+f[75])+4.0*f[74]-5.0*(f[63]+f[60]))+6.708203932499369*(5.0*(f[39]+f[38])-4.0*f[50])+2.0*(3.0*(3.0*(2.23606797749979*(f[15]+f[9])-3.0*f[28])-5.0*f[4])-13.41640786499874*f[37])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.05*(9.0*f[77]-6.708203932499369*(f[46]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*(4.0*f[47]-5.0*(f[45]+f[44])+4.0*f[43])-5.0*(6.708203932499369*(f[35]+f[32])+4.0*f[20]-5.0*(f[18]+f[17])))+2.0*((-22.3606797749979*f[16])+3.0*(11.18033988749895*(f[5]+f[1])-15.0*f[12])-25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.008333333333333333*(60.37383539249431*(f[98]+f[95])+9.0*(4.0*f[83]-5.0*(f[76]+f[75])+4.0*f[74]-5.0*(f[63]+f[60]))+6.708203932499369*(5.0*(f[39]+f[38])-4.0*f[50])+2.0*(3.0*(3.0*(2.23606797749979*(f[15]+f[9])-3.0*f[28])-5.0*f[4])-13.41640786499874*f[37])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.05*(9.0*f[77]-6.708203932499369*(f[46]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]+60.37383539249431*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96])+5.0*(f[95]-1.0*f[88]))+9.0*(5.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62])-1.0*(5.0*(4.0*f[61]+5.0*f[60])+40.24922359499622*f[54]))+33.54101966249684*(5.0*f[38]-4.0*(f[50]+f[39]))+2.0*(3.0*(3.0*(15.0*(f[30]-1.0*f[28]+f[23])+11.18033988749895*(f[15]-1.0*f[11]+f[9]))-25.0*f[4])-67.08203932499369*f[37])))/(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*f[69])+6.708203932499369*(4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34])-1.0*(6.708203932499369*(4.0*f[33]+5.0*f[32])+54.0*f[26]-5.0*(5.0*f[17]-4.0*(f[20]+f[18]))))+2.0*((-22.3606797749979*f[16])+3.0*(15.0*(f[14]-1.0*f[12]+f[7])+11.18033988749895*(f[5]-1.0*f[3]+f[1]))-25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*f[69])+6.708203932499369*(4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34])-1.0*(6.708203932499369*(4.0*f[33]+5.0*f[32])+54.0*f[26]-5.0*(5.0*f[17]-4.0*(f[20]+f[18]))))+2.0*((-22.3606797749979*f[16])+3.0*(15.0*(f[14]-1.0*f[12]+f[7])+11.18033988749895*(f[5]-1.0*f[3]+f[1]))-25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*f[69])+6.708203932499369*(4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34])-1.0*(6.708203932499369*(4.0*f[33]+5.0*f[32])+54.0*f[26]-5.0*(5.0*f[17]-4.0*(f[20]+f[18]))))+2.0*((-22.3606797749979*f[16])+3.0*(15.0*(f[14]-1.0*f[12]+f[7])+11.18033988749895*(f[5]-1.0*f[3]+f[1]))-25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.001666666666666667*(405.0*f[108]+60.37383539249431*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96])+5.0*(f[95]-1.0*f[88]))+9.0*(5.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62])-1.0*(5.0*(4.0*f[61]+5.0*f[60])+40.24922359499622*f[54]))+33.54101966249684*(5.0*f[38]-4.0*(f[50]+f[39]))+2.0*(3.0*(3.0*(15.0*(f[30]-1.0*f[28]+f[23])+11.18033988749895*(f[15]-1.0*f[11]+f[9]))-25.0*f[4])-67.08203932499369*f[37])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(60.37383539249431*f[101]+5.0*(9.0*((-1.0*f[79])+f[77]-1.0*f[66])+6.708203932499369*((-1.0*f[46])+f[42]-1.0*f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*f[81]-1.0*(4.0*f[72]+5.0*f[71]-4.0*f[70])+5.0*f[69])+6.708203932499369*(4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34])-1.0*(6.708203932499369*(4.0*f[33]+5.0*f[32])+54.0*f[26]-5.0*(5.0*f[17]-4.0*(f[20]+f[18]))))+2.0*((-22.3606797749979*f[16])+3.0*(15.0*(f[14]-1.0*f[12]+f[7])+11.18033988749895*(f[5]-1.0*f[3]+f[1]))-25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.001666666666666667*(405.0*f[108]+60.37383539249431*(4.0*f[105]-1.0*(4.0*f[98]+5.0*f[97]-4.0*f[96])+5.0*(f[95]-1.0*f[88]))+9.0*(5.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62])-1.0*(5.0*(4.0*f[61]+5.0*f[60])+40.24922359499622*f[54]))+33.54101966249684*(5.0*f[38]-4.0*(f[50]+f[39]))+2.0*(3.0*(3.0*(15.0*(f[30]-1.0*f[28]+f[23])+11.18033988749895*(f[15]-1.0*f[11]+f[9]))-25.0*f[4])-67.08203932499369*f[37])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(60.37383539249431*f[101]+5.0*(9.0*((-1.0*f[79])+f[77]-1.0*f[66])+6.708203932499369*((-1.0*f[46])+f[42]-1.0*f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[105]+f[88])+9.0*((-5.0*(f[85]+f[83]))+4.0*f[63]-5.0*f[62]+4.0*f[61]-5.0*f[60])+6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*f[38])+2.0*(3.0*(3.0*(2.23606797749979*(f[11]+f[9])-3.0*f[23])-5.0*f[4])-13.41640786499874*f[37])))/(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*((-5.0*(f[49]+f[47]))+4.0*f[35]-5.0*f[34]+4.0*f[33])+5.0*((-6.708203932499369*f[32])+5.0*f[20]-4.0*f[18]+5.0*f[17]))+2.0*((-22.3606797749979*f[16])+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[7])-25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*((-5.0*(f[49]+f[47]))+4.0*f[35]-5.0*f[34]+4.0*f[33])+5.0*((-6.708203932499369*f[32])+5.0*f[20]-4.0*f[18]+5.0*f[17]))+2.0*((-22.3606797749979*f[16])+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[7])-25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*((-5.0*(f[49]+f[47]))+4.0*f[35]-5.0*f[34]+4.0*f[33])+5.0*((-6.708203932499369*f[32])+5.0*f[20]-4.0*f[18]+5.0*f[17]))+2.0*((-22.3606797749979*f[16])+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[7])-25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.008333333333333333*(60.37383539249431*(f[105]+f[88])+9.0*((-5.0*(f[85]+f[83]))+4.0*f[63]-5.0*f[62]+4.0*f[61]-5.0*f[60])+6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*f[38])+2.0*(3.0*(3.0*(2.23606797749979*(f[11]+f[9])-3.0*f[23])-5.0*f[4])-13.41640786499874*f[37])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.05*(9.0*f[66]-6.708203932499369*(f[42]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*((-5.0*(f[49]+f[47]))+4.0*f[35]-5.0*f[34]+4.0*f[33])+5.0*((-6.708203932499369*f[32])+5.0*f[20]-4.0*f[18]+5.0*f[17]))+2.0*((-22.3606797749979*f[16])+3.0*(11.18033988749895*(f[3]+f[1])-15.0*f[7])-25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.008333333333333333*(60.37383539249431*(f[105]+f[88])+9.0*((-5.0*(f[85]+f[83]))+4.0*f[63]-5.0*f[62]+4.0*f[61]-5.0*f[60])+6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*f[38])+2.0*(3.0*(3.0*(2.23606797749979*(f[11]+f[9])-3.0*f[23])-5.0*f[4])-13.41640786499874*f[37])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.05*(9.0*f[66]-6.708203932499369*(f[42]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[105]+f[88])+9.0*(5.0*(f[83]-1.0*f[85])-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+6.708203932499369*((-5.0*f[50])+4.0*f[39]-5.0*f[38])+2.0*(13.41640786499874*f[37]+3.0*(3.0*(2.23606797749979*(f[11]-1.0*f[9])-3.0*f[23])+5.0*f[4]))))/(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*((-5.0*f[49])+5.0*f[47]-1.0*(4.0*f[35]+5.0*f[34]-4.0*f[33]))+5.0*(6.708203932499369*f[32]-5.0*f[20]+4.0*f[18]-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[3]-1.0*f[1])-15.0*f[7])+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*((-5.0*f[49])+5.0*f[47]-1.0*(4.0*f[35]+5.0*f[34]-4.0*f[33]))+5.0*(6.708203932499369*f[32]-5.0*f[20]+4.0*f[18]-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[3]-1.0*f[1])-15.0*f[7])+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*((-5.0*f[49])+5.0*f[47]-1.0*(4.0*f[35]+5.0*f[34]-4.0*f[33]))+5.0*(6.708203932499369*f[32]-5.0*f[20]+4.0*f[18]-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[3]-1.0*f[1])-15.0*f[7])+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(60.37383539249431*(f[105]+f[88])+9.0*(5.0*(f[83]-1.0*f[85])-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+6.708203932499369*((-5.0*f[50])+4.0*f[39]-5.0*f[38])+2.0*(13.41640786499874*f[37]+3.0*(3.0*(2.23606797749979*(f[11]-1.0*f[9])-3.0*f[23])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(9.0*f[66]+6.708203932499369*(f[40]-1.0*f[42])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*((-5.0*f[49])+5.0*f[47]-1.0*(4.0*f[35]+5.0*f[34]-4.0*f[33]))+5.0*(6.708203932499369*f[32]-5.0*f[20]+4.0*f[18]-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[3]-1.0*f[1])-15.0*f[7])+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.008333333333333333*(60.37383539249431*(f[105]+f[88])+9.0*(5.0*(f[83]-1.0*f[85])-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))+6.708203932499369*((-5.0*f[50])+4.0*f[39]-5.0*f[38])+2.0*(13.41640786499874*f[37]+3.0*(3.0*(2.23606797749979*(f[11]-1.0*f[9])-3.0*f[23])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.05*(9.0*f[66]+6.708203932499369*(f[40]-1.0*f[42])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]+60.37383539249431*(4.0*(f[98]-1.0*f[105])-5.0*f[97]+4.0*f[96]+5.0*(f[88]-1.0*f[95]))+9.0*(5.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61])-1.0*(25.0*f[60]+40.24922359499622*f[54]))+33.54101966249684*(5.0*f[38]-4.0*(f[50]+f[39]))+2.0*(3.0*(3.0*(15.0*(f[30]+f[28]-1.0*f[23])+11.18033988749895*((-1.0*f[15])+f[11]+f[9]))-25.0*f[4])-67.08203932499369*f[37])))/(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[72]-1.0*f[81])-5.0*f[71]+4.0*f[70]-5.0*f[69]+5.0*f[57])+6.708203932499369*(4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33])-1.0*(33.54101966249684*f[32]+54.0*f[26]-5.0*(5.0*f[17]-4.0*(f[20]+f[18]))))+2.0*((-22.3606797749979*f[16])+3.0*(15.0*(f[14]+f[12]-1.0*f[7])+11.18033988749895*((-1.0*f[5])+f[3]+f[1]))-25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[72]-1.0*f[81])-5.0*f[71]+4.0*f[70]-5.0*f[69]+5.0*f[57])+6.708203932499369*(4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33])-1.0*(33.54101966249684*f[32]+54.0*f[26]-5.0*(5.0*f[17]-4.0*(f[20]+f[18]))))+2.0*((-22.3606797749979*f[16])+3.0*(15.0*(f[14]+f[12]-1.0*f[7])+11.18033988749895*((-1.0*f[5])+f[3]+f[1]))-25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[72]-1.0*f[81])-5.0*f[71]+4.0*f[70]-5.0*f[69]+5.0*f[57])+6.708203932499369*(4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33])-1.0*(33.54101966249684*f[32]+54.0*f[26]-5.0*(5.0*f[17]-4.0*(f[20]+f[18]))))+2.0*((-22.3606797749979*f[16])+3.0*(15.0*(f[14]+f[12]-1.0*f[7])+11.18033988749895*((-1.0*f[5])+f[3]+f[1]))-25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.001666666666666667*(405.0*f[108]+60.37383539249431*(4.0*(f[98]-1.0*f[105])-5.0*f[97]+4.0*f[96]+5.0*(f[88]-1.0*f[95]))+9.0*(5.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61])-1.0*(25.0*f[60]+40.24922359499622*f[54]))+33.54101966249684*(5.0*f[38]-4.0*(f[50]+f[39]))+2.0*(3.0*(3.0*(15.0*(f[30]+f[28]-1.0*f[23])+11.18033988749895*((-1.0*f[15])+f[11]+f[9]))-25.0*f[4])-67.08203932499369*f[37])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(60.37383539249431*f[101]+5.0*(9.0*(f[66]-1.0*(f[79]+f[77]))+6.708203932499369*(f[46]-1.0*(f[42]+f[40]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[72]-1.0*f[81])-5.0*f[71]+4.0*f[70]-5.0*f[69]+5.0*f[57])+6.708203932499369*(4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33])-1.0*(33.54101966249684*f[32]+54.0*f[26]-5.0*(5.0*f[17]-4.0*(f[20]+f[18]))))+2.0*((-22.3606797749979*f[16])+3.0*(15.0*(f[14]+f[12]-1.0*f[7])+11.18033988749895*((-1.0*f[5])+f[3]+f[1]))-25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.001666666666666667*(405.0*f[108]+60.37383539249431*(4.0*(f[98]-1.0*f[105])-5.0*f[97]+4.0*f[96]+5.0*(f[88]-1.0*f[95]))+9.0*(5.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61])-1.0*(25.0*f[60]+40.24922359499622*f[54]))+33.54101966249684*(5.0*f[38]-4.0*(f[50]+f[39]))+2.0*(3.0*(3.0*(15.0*(f[30]+f[28]-1.0*f[23])+11.18033988749895*((-1.0*f[15])+f[11]+f[9]))-25.0*f[4])-67.08203932499369*f[37])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(60.37383539249431*f[101]+5.0*(9.0*(f[66]-1.0*(f[79]+f[77]))+6.708203932499369*(f[46]-1.0*(f[42]+f[40]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[98]+f[95])+9.0*((-1.0*(4.0*f[83]+5.0*(f[76]+f[75])))+4.0*f[74]+5.0*(f[63]+f[60]))+6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]))+2.0*(13.41640786499874*f[37]+3.0*(3.0*(2.23606797749979*(f[15]-1.0*f[9])-3.0*f[28])+5.0*f[4]))))/(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*(4.0*f[43]-1.0*(4.0*f[47]+5.0*(f[45]+f[44])))+5.0*(6.708203932499369*(f[35]+f[32])+4.0*f[20]-5.0*(f[18]+f[17])))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[5]-1.0*f[1])-15.0*f[12])+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*(4.0*f[43]-1.0*(4.0*f[47]+5.0*(f[45]+f[44])))+5.0*(6.708203932499369*(f[35]+f[32])+4.0*f[20]-5.0*(f[18]+f[17])))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[5]-1.0*f[1])-15.0*f[12])+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*(4.0*f[43]-1.0*(4.0*f[47]+5.0*(f[45]+f[44])))+5.0*(6.708203932499369*(f[35]+f[32])+4.0*f[20]-5.0*(f[18]+f[17])))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[5]-1.0*f[1])-15.0*f[12])+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.008333333333333333*(60.37383539249431*(f[98]+f[95])+9.0*((-1.0*(4.0*f[83]+5.0*(f[76]+f[75])))+4.0*f[74]+5.0*(f[63]+f[60]))+6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]))+2.0*(13.41640786499874*f[37]+3.0*(3.0*(2.23606797749979*(f[15]-1.0*f[9])-3.0*f[28])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.05*(9.0*f[77]+6.708203932499369*(f[40]-1.0*f[46])-5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*(4.0*f[43]-1.0*(4.0*f[47]+5.0*(f[45]+f[44])))+5.0*(6.708203932499369*(f[35]+f[32])+4.0*f[20]-5.0*(f[18]+f[17])))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[5]-1.0*f[1])-15.0*f[12])+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.008333333333333333*(60.37383539249431*(f[98]+f[95])+9.0*((-1.0*(4.0*f[83]+5.0*(f[76]+f[75])))+4.0*f[74]+5.0*(f[63]+f[60]))+6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]))+2.0*(13.41640786499874*f[37]+3.0*(3.0*(2.23606797749979*(f[15]-1.0*f[9])-3.0*f[28])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.05*(9.0*f[77]+6.708203932499369*(f[40]-1.0*f[46])-5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]+60.37383539249431*(5.0*(f[95]+f[88])-1.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96]))+9.0*(5.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))-40.24922359499622*f[54])+33.54101966249684*(4.0*(f[50]+f[39])-5.0*f[38])+2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*(f[30]-1.0*(f[28]+f[23]))+11.18033988749895*(f[15]+f[11]-1.0*f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[92]+9.0*(5.0*(f[69]+f[57])-1.0*(4.0*(f[81]+f[72])+5.0*f[71]-4.0*f[70]))+6.708203932499369*(4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))-54.0*f[26]+5.0*(4.0*(f[20]+f[18])-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[14]-1.0*(f[12]+f[7]))+11.18033988749895*(f[5]+f[3]-1.0*f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(5.0*(f[69]+f[57])-1.0*(4.0*(f[81]+f[72])+5.0*f[71]-4.0*f[70]))+6.708203932499369*(4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))-54.0*f[26]+5.0*(4.0*(f[20]+f[18])-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[14]-1.0*(f[12]+f[7]))+11.18033988749895*(f[5]+f[3]-1.0*f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(5.0*(f[69]+f[57])-1.0*(4.0*(f[81]+f[72])+5.0*f[71]-4.0*f[70]))+6.708203932499369*(4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))-54.0*f[26]+5.0*(4.0*(f[20]+f[18])-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[14]-1.0*(f[12]+f[7]))+11.18033988749895*(f[5]+f[3]-1.0*f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.001666666666666667*(405.0*f[108]+60.37383539249431*(5.0*(f[95]+f[88])-1.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96]))+9.0*(5.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))-40.24922359499622*f[54])+33.54101966249684*(4.0*(f[50]+f[39])-5.0*f[38])+2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*(f[30]-1.0*(f[28]+f[23]))+11.18033988749895*(f[15]+f[11]-1.0*f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(60.37383539249431*f[101]+5.0*(9.0*((-1.0*f[79])+f[77]+f[66])+6.708203932499369*(f[40]-1.0*(f[46]+f[42]))-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(5.0*(f[69]+f[57])-1.0*(4.0*(f[81]+f[72])+5.0*f[71]-4.0*f[70]))+6.708203932499369*(4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))-54.0*f[26]+5.0*(4.0*(f[20]+f[18])-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[14]-1.0*(f[12]+f[7]))+11.18033988749895*(f[5]+f[3]-1.0*f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.001666666666666667*(405.0*f[108]+60.37383539249431*(5.0*(f[95]+f[88])-1.0*(4.0*(f[105]+f[98])+5.0*f[97]-4.0*f[96]))+9.0*(5.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))-40.24922359499622*f[54])+33.54101966249684*(4.0*(f[50]+f[39])-5.0*f[38])+2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*(f[30]-1.0*(f[28]+f[23]))+11.18033988749895*(f[15]+f[11]-1.0*f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(60.37383539249431*f[101]+5.0*(9.0*((-1.0*f[79])+f[77]+f[66])+6.708203932499369*(f[40]-1.0*(f[46]+f[42]))-5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[3][0] = 0.05555555555555555*(fReflXYZMuQuad[7][0]+8.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0]+8.0*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[3][1] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][0]-4188761.0*fReflXYZMuQuad[6][0])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])-4.63256860547201e+14*fReflXYZMuQuad[5][0])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][0]-1.0*(8377522.0*fReflXYZMuQuad[1][0]+2.266096151179001e+23*fReflXYZMuQuad[0][0])))); 
  fReflXYQuad[3][2] = 0.05555555555555555*(fReflXYZMuQuad[7][1]+8.0*fReflXYZMuQuad[6][1]+fReflXYZMuQuad[5][1]+8.0*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[3][3] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][0]+7.4121097687552e+14*fReflXYZMuQuad[6][0]+4.63256860547201e+14*fReflXYZMuQuad[5][0])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][0]+7.4121097687552e+14*fReflXYZMuQuad[1][0]+4.63256860547201e+14*fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[3][4] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][1]-4188761.0*fReflXYZMuQuad[6][1])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])-4.63256860547201e+14*fReflXYZMuQuad[5][1])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][1]-1.0*(8377522.0*fReflXYZMuQuad[1][1]+2.266096151179001e+23*fReflXYZMuQuad[0][1])))); 
  fReflXYQuad[3][5] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][0]+9.0*fReflXYZMuQuad[6][0])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][0]-1.346286087882789e+17*fReflXYZMuQuad[5][0])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][0]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])))); 
  fReflXYQuad[3][6] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][1]+7.4121097687552e+14*fReflXYZMuQuad[6][1]+4.63256860547201e+14*fReflXYZMuQuad[5][1])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][1]+7.4121097687552e+14*fReflXYZMuQuad[1][1]+4.63256860547201e+14*fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[3][7] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][0]-4.17103711299121e+14*fReflXYZMuQuad[6][0])+5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[5][0]+2384663.0*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])+3.190559553141742e+22*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[3][8] = 0.05555555555555555*(fReflXYZMuQuad[7][2]+8.0*fReflXYZMuQuad[6][2]+fReflXYZMuQuad[5][2]+8.0*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[3][9] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][0]+1.6692641e+7*fReflXYZMuQuad[6][0])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][0]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[2][0]+2384663.0*fReflXYZMuQuad[1][0]+3.190559553141742e+22*fReflXYZMuQuad[0][0])-3.20880527843592e+30*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])))); 
  fReflXYQuad[3][10] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][1]+9.0*fReflXYZMuQuad[6][1])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][1]-1.346286087882789e+17*fReflXYZMuQuad[5][1])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][1]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])))); 
  fReflXYQuad[3][11] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][1]-4.17103711299121e+14*fReflXYZMuQuad[6][1])+5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[5][1]+2384663.0*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])+3.190559553141743e+22*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[3][12] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][2]-4188761.0*fReflXYZMuQuad[6][2])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYZMuQuad[4][2]-4.63256860547201e+14*fReflXYZMuQuad[5][2])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][2]-1.0*(8377522.0*fReflXYZMuQuad[1][2]+2.266096151179001e+23*fReflXYZMuQuad[0][2])))); 
  fReflXYQuad[3][13] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][0]-2.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0])-27.0*fReflXYZMuQuad[3][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[2][0])+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[3][14] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][2]+7.4121097687552e+14*fReflXYZMuQuad[6][2]+4.63256860547201e+14*fReflXYZMuQuad[5][2])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][2]+7.4121097687552e+14*fReflXYZMuQuad[1][2]+4.63256860547201e+14*fReflXYZMuQuad[0][2]))); 
  fReflXYQuad[3][15] = 1.719407810605221e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][0]+27.0*fReflXYZMuQuad[6][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[5][0])+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[4][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[3][16] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][1]+1.6692641e+7*fReflXYZMuQuad[6][1])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][1]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[2][1]+2384663.0*fReflXYZMuQuad[1][1]+3.190559553141743e+22*fReflXYZMuQuad[0][1])-3.20880527843592e+30*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])))); 
  fReflXYQuad[3][17] = 1.719407810605222e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][1]+fReflXYZMuQuad[5][1])-27.0*fReflXYZMuQuad[3][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[2][1])+10.0*fReflXYZMuQuad[1][1]-5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[3][18] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][2]+9.0*fReflXYZMuQuad[6][2])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][2]-1.346286087882789e+17*fReflXYZMuQuad[5][2])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][2]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])))); 
  fReflXYQuad[3][19] = 1.719407810605222e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][1]+27.0*fReflXYZMuQuad[6][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[5][1])+2.0*(5.0*fReflXYZMuQuad[3][1]-5.0*fReflXYZMuQuad[4][1])+5.0*fReflXYZMuQuad[2][1])); 
  } 

 
// node (x,y)_5 
  vcutSq_i = (0.05*q_*(zVal*((63.63961030678928*phiWall[15]-63.63961030678928*phi[15]+47.43416490252571*phiWall[9]-47.43416490252571*phi[9])*zVal-36.74234614174767*phiWall[18]+36.74234614174767*phi[18]-27.38612787525831*phiWall[14]+27.38612787525831*phi[14]+21.90890230020666*phiWall[13]-21.90890230020666*phi[13]+32.86335345030997*phiWall[5]-32.86335345030997*phi[5]+24.49489742783179*phiWall[3]-24.49489742783179*phi[3])-21.21320343559643*phiWall[15]+21.21320343559643*phi[15]-21.21320343559643*phiWall[12]+21.21320343559643*phi[12]-15.8113883008419*phiWall[9]+15.8113883008419*phi[9]-15.8113883008419*phiWall[8]+15.8113883008419*phi[8]+12.64911064067352*phiWall[7]-12.64911064067352*phi[7]+18.97366596101028*phiWall[1]-18.97366596101028*phi[1]+14.14213562373095*phiWall[0]-14.14213562373095*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[4][0] = 0.0; 
  fReflXYQuad[4][1] = 0.0; 
  fReflXYQuad[4][2] = 0.0; 
  fReflXYQuad[4][3] = 0.0; 
  fReflXYQuad[4][4] = 0.0; 
  fReflXYQuad[4][5] = 0.0; 
  fReflXYQuad[4][6] = 0.0; 
  fReflXYQuad[4][7] = 0.0; 
  fReflXYQuad[4][8] = 0.0; 
  fReflXYQuad[4][9] = 0.0; 
  fReflXYQuad[4][10] = 0.0; 
  fReflXYQuad[4][11] = 0.0; 
  fReflXYQuad[4][12] = 0.0; 
  fReflXYQuad[4][13] = 0.0; 
  fReflXYQuad[4][14] = 0.0; 
  fReflXYQuad[4][15] = 0.0; 
  fReflXYQuad[4][16] = 0.0; 
  fReflXYQuad[4][17] = 0.0; 
  fReflXYQuad[4][18] = 0.0; 
  fReflXYQuad[4][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[4][0] = -0.05*(2.23606797749979*(6.708203932499369*f[32]+5.0*f[17])-2.0*(2.23606797749979*(2.0*f[16]+3.0*f[1])+5.0*f[0])); 
  fReflXYQuad[4][1] = -0.01666666666666667*(5.0*(9.0*f[57]+6.708203932499369*f[34])-2.0*(13.41640786499874*f[33]+3.0*(6.708203932499369*f[7]+5.0*f[3]))); 
  fReflXYQuad[4][2] = -0.01666666666666667*(5.0*(9.0*f[60]+6.708203932499369*f[38])-2.0*(13.41640786499874*f[37]+3.0*(6.708203932499369*f[9]+5.0*f[4]))); 
  fReflXYQuad[4][3] = -0.01666666666666667*(5.0*(9.0*f[69]+6.708203932499369*f[44])-2.0*(13.41640786499874*f[43]+3.0*(6.708203932499369*f[12]+5.0*f[5]))); 
  fReflXYQuad[4][4] = -0.05*(2.23606797749979*(6.708203932499369*f[88]+5.0*f[62])-2.0*(2.23606797749979*(2.0*f[61]+3.0*f[23])+5.0*f[11])); 
  fReflXYQuad[4][5] = -0.05*(2.23606797749979*(6.708203932499369*f[92]+5.0*f[71])-2.0*(2.23606797749979*(2.0*f[70]+3.0*f[26])+5.0*f[14])); 
  fReflXYQuad[4][6] = -0.05*(2.23606797749979*(6.708203932499369*f[95]+5.0*f[75])-2.0*(2.23606797749979*(2.0*f[74]+3.0*f[28])+5.0*f[15])); 
  fReflXYQuad[4][7] = 0.1*(6.708203932499369*f[35]+5.0*f[18]); 
  fReflXYQuad[4][8] = 0.1*(6.708203932499369*f[40]+5.0*f[19]); 
  fReflXYQuad[4][9] = 0.1*(6.708203932499369*f[47]+5.0*f[20]); 
  fReflXYQuad[4][10] = -0.01666666666666667*(5.0*(9.0*f[108]+6.708203932499369*f[97])-2.0*(13.41640786499874*f[96]+3.0*(6.708203932499369*f[54]+5.0*f[30]))); 
  fReflXYQuad[4][11] = 0.1*(6.708203932499369*f[63]+5.0*f[39]); 
  fReflXYQuad[4][12] = 0.1*(6.708203932499369*f[66]+5.0*f[42]); 
  fReflXYQuad[4][13] = 0.1*(6.708203932499369*f[72]+5.0*f[45]); 
  fReflXYQuad[4][14] = 0.1*(6.708203932499369*f[77]+5.0*f[46]); 
  fReflXYQuad[4][15] = 0.1*(6.708203932499369*f[81]+5.0*f[49]); 
  fReflXYQuad[4][16] = 0.1*(6.708203932499369*f[83]+5.0*f[50]); 
  fReflXYQuad[4][17] = 0.1*(6.708203932499369*f[98]+5.0*f[76]); 
  fReflXYQuad[4][18] = 0.1*(6.708203932499369*f[101]+5.0*f[79]); 
  fReflXYQuad[4][19] = 0.1*(6.708203932499369*f[105]+5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(405.0*f[108]+60.37383539249431*(4.0*(f[105]+f[98])+5.0*f[97]-1.0*(4.0*f[96]+5.0*(f[95]+f[88])))+9.0*(5.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))-40.24922359499622*f[54])+33.54101966249684*(5.0*f[38]-4.0*(f[50]+f[39]))+2.0*(3.0*(3.0*(15.0*((-1.0*f[30])+f[28]+f[23])+11.18033988749895*(f[15]+f[11]))-1.0*(33.54101966249685*f[9]+25.0*f[4]))-67.08203932499369*f[37])))/(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])))+6.708203932499369*(4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))-54.0*f[26]+5.0*(5.0*f[17]-4.0*(f[20]+f[18])))+2.0*((-22.3606797749979*f[16])+3.0*(15.0*((-1.0*f[14])+f[12]+f[7])+11.18033988749895*(f[5]+f[3]))-1.0*(33.54101966249685*f[1]+25.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])))+6.708203932499369*(4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))-54.0*f[26]+5.0*(5.0*f[17]-4.0*(f[20]+f[18])))+2.0*((-22.3606797749979*f[16])+3.0*(15.0*((-1.0*f[14])+f[12]+f[7])+11.18033988749895*(f[5]+f[3]))-1.0*(33.54101966249685*f[1]+25.0*f[0]))) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])))+6.708203932499369*(4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))-54.0*f[26]+5.0*(5.0*f[17]-4.0*(f[20]+f[18])))+2.0*((-22.3606797749979*f[16])+3.0*(15.0*((-1.0*f[14])+f[12]+f[7])+11.18033988749895*(f[5]+f[3]))-1.0*(33.54101966249685*f[1]+25.0*f[0]))))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.001666666666666667*(405.0*f[108]+60.37383539249431*(4.0*(f[105]+f[98])+5.0*f[97]-1.0*(4.0*f[96]+5.0*(f[95]+f[88])))+9.0*(5.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))-40.24922359499622*f[54])+33.54101966249684*(5.0*f[38]-4.0*(f[50]+f[39]))+2.0*(3.0*(3.0*(15.0*((-1.0*f[30])+f[28]+f[23])+11.18033988749895*(f[15]+f[11]))-1.0*(33.54101966249685*f[9]+25.0*f[4]))-67.08203932499369*f[37])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(60.37383539249431*f[101]+5.0*(9.0*(f[79]-1.0*(f[77]+f[66]))+6.708203932499369*(f[40]-1.0*(f[46]+f[42]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[81]+f[72])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]+f[57])))+6.708203932499369*(4.0*(f[49]-1.0*f[47]+f[45])-5.0*f[44]+4.0*f[43]-1.0*(4.0*f[35]+5.0*f[34]-1.0*(4.0*f[33]+5.0*f[32])))-54.0*f[26]+5.0*(5.0*f[17]-4.0*(f[20]+f[18])))+2.0*((-22.3606797749979*f[16])+3.0*(15.0*((-1.0*f[14])+f[12]+f[7])+11.18033988749895*(f[5]+f[3]))-1.0*(33.54101966249685*f[1]+25.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.001666666666666667*(405.0*f[108]+60.37383539249431*(4.0*(f[105]+f[98])+5.0*f[97]-1.0*(4.0*f[96]+5.0*(f[95]+f[88])))+9.0*(5.0*(4.0*(f[85]-1.0*f[83]+f[76])-5.0*f[75]+4.0*f[74]-1.0*(4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60])))-40.24922359499622*f[54])+33.54101966249684*(5.0*f[38]-4.0*(f[50]+f[39]))+2.0*(3.0*(3.0*(15.0*((-1.0*f[30])+f[28]+f[23])+11.18033988749895*(f[15]+f[11]))-1.0*(33.54101966249685*f[9]+25.0*f[4]))-67.08203932499369*f[37])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(60.37383539249431*f[101]+5.0*(9.0*(f[79]-1.0*(f[77]+f[66]))+6.708203932499369*(f[40]-1.0*(f[46]+f[42]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[98]+f[95])+9.0*(4.0*f[83]+5.0*(f[76]+f[75])-1.0*(4.0*f[74]+5.0*(f[63]+f[60])))+6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]))+2.0*(13.41640786499874*f[37]+3.0*(3.0*(2.23606797749979*(f[9]-1.0*f[15])-3.0*f[28])+5.0*f[4]))))/(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-4.0*f[43])+5.0*((-6.708203932499369*(f[35]+f[32]))+4.0*f[20]-5.0*(f[18]+f[17])))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[1]-1.0*f[5])-15.0*f[12])+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-4.0*f[43])+5.0*((-6.708203932499369*(f[35]+f[32]))+4.0*f[20]-5.0*(f[18]+f[17])))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[1]-1.0*f[5])-15.0*f[12])+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-4.0*f[43])+5.0*((-6.708203932499369*(f[35]+f[32]))+4.0*f[20]-5.0*(f[18]+f[17])))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[1]-1.0*f[5])-15.0*f[12])+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(60.37383539249431*(f[98]+f[95])+9.0*(4.0*f[83]+5.0*(f[76]+f[75])-1.0*(4.0*f[74]+5.0*(f[63]+f[60])))+6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]))+2.0*(13.41640786499874*f[37]+3.0*(3.0*(2.23606797749979*(f[9]-1.0*f[15])-3.0*f[28])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(9.0*f[77]+6.708203932499369*f[46]-1.0*(6.708203932499369*f[40]+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*(4.0*f[47]+5.0*(f[45]+f[44])-4.0*f[43])+5.0*((-6.708203932499369*(f[35]+f[32]))+4.0*f[20]-5.0*(f[18]+f[17])))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[1]-1.0*f[5])-15.0*f[12])+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(60.37383539249431*(f[98]+f[95])+9.0*(4.0*f[83]+5.0*(f[76]+f[75])-1.0*(4.0*f[74]+5.0*(f[63]+f[60])))+6.708203932499369*(4.0*f[50]-5.0*(f[39]+f[38]))+2.0*(13.41640786499874*f[37]+3.0*(3.0*(2.23606797749979*(f[9]-1.0*f[15])-3.0*f[28])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(9.0*f[77]+6.708203932499369*f[46]-1.0*(6.708203932499369*f[40]+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]+60.37383539249431*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96]+5.0*(f[95]-1.0*f[88]))+9.0*(5.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61])-1.0*(25.0*f[60]+40.24922359499622*f[54]))+33.54101966249684*(4.0*(f[50]+f[39])-5.0*f[38])+2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*(f[23]-1.0*(f[30]+f[28]))+11.18033988749895*((-1.0*f[15])+f[11]+f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*f[69])+6.708203932499369*(4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33])-1.0*(33.54101966249684*f[32]+54.0*f[26]-5.0*(4.0*(f[20]+f[18])-5.0*f[17])))+2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[7]-1.0*(f[14]+f[12]))+11.18033988749895*((-1.0*f[5])+f[3]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*f[69])+6.708203932499369*(4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33])-1.0*(33.54101966249684*f[32]+54.0*f[26]-5.0*(4.0*(f[20]+f[18])-5.0*f[17])))+2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[7]-1.0*(f[14]+f[12]))+11.18033988749895*((-1.0*f[5])+f[3]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*f[69])+6.708203932499369*(4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33])-1.0*(33.54101966249684*f[32]+54.0*f[26]-5.0*(4.0*(f[20]+f[18])-5.0*f[17])))+2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[7]-1.0*(f[14]+f[12]))+11.18033988749895*((-1.0*f[5])+f[3]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.001666666666666667*(405.0*f[108]+60.37383539249431*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96]+5.0*(f[95]-1.0*f[88]))+9.0*(5.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61])-1.0*(25.0*f[60]+40.24922359499622*f[54]))+33.54101966249684*(4.0*(f[50]+f[39])-5.0*f[38])+2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*(f[23]-1.0*(f[30]+f[28]))+11.18033988749895*((-1.0*f[15])+f[11]+f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(60.37383539249431*f[101]+5.0*(9.0*(f[79]+f[77]-1.0*f[66])+6.708203932499369*f[46]-1.0*(6.708203932499369*(f[42]+f[40])+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[81]-1.0*f[72])+5.0*f[71]-4.0*f[70]+5.0*f[69])+6.708203932499369*(4.0*(f[49]+f[47]-1.0*f[45])+5.0*f[44]+4.0*(f[35]-1.0*f[43])-5.0*f[34]+4.0*f[33])-1.0*(33.54101966249684*f[32]+54.0*f[26]-5.0*(4.0*(f[20]+f[18])-5.0*f[17])))+2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[7]-1.0*(f[14]+f[12]))+11.18033988749895*((-1.0*f[5])+f[3]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.001666666666666667*(405.0*f[108]+60.37383539249431*(4.0*(f[105]-1.0*f[98])+5.0*f[97]-4.0*f[96]+5.0*(f[95]-1.0*f[88]))+9.0*(5.0*(4.0*(f[85]+f[83]-1.0*f[76])+5.0*f[75]+4.0*(f[63]-1.0*f[74])-5.0*f[62]+4.0*f[61])-1.0*(25.0*f[60]+40.24922359499622*f[54]))+33.54101966249684*(4.0*(f[50]+f[39])-5.0*f[38])+2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*(f[23]-1.0*(f[30]+f[28]))+11.18033988749895*((-1.0*f[15])+f[11]+f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(60.37383539249431*f[101]+5.0*(9.0*(f[79]+f[77]-1.0*f[66])+6.708203932499369*f[46]-1.0*(6.708203932499369*(f[42]+f[40])+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[105]+f[88])+9.0*(5.0*(f[85]-1.0*f[83])+4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+6.708203932499369*((-5.0*f[50])+4.0*f[39]-5.0*f[38])+2.0*(13.41640786499874*f[37]+3.0*(3.0*(2.23606797749979*(f[9]-1.0*f[11])-3.0*f[23])+5.0*f[4]))))/(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*(5.0*f[49]+4.0*f[35]+5.0*f[34]-4.0*f[33])+5.0*((-1.0*(6.708203932499369*f[32]+5.0*f[20]))+4.0*f[18]-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[7])+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*(5.0*f[49]+4.0*f[35]+5.0*f[34]-4.0*f[33])+5.0*((-1.0*(6.708203932499369*f[32]+5.0*f[20]))+4.0*f[18]-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[7])+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*(5.0*f[49]+4.0*f[35]+5.0*f[34]-4.0*f[33])+5.0*((-1.0*(6.708203932499369*f[32]+5.0*f[20]))+4.0*f[18]-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[7])+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(60.37383539249431*(f[105]+f[88])+9.0*(5.0*(f[85]-1.0*f[83])+4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+6.708203932499369*((-5.0*f[50])+4.0*f[39]-5.0*f[38])+2.0*(13.41640786499874*f[37]+3.0*(3.0*(2.23606797749979*(f[9]-1.0*f[11])-3.0*f[23])+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(9.0*f[66]+6.708203932499369*f[42]-1.0*(6.708203932499369*f[40]+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*(5.0*f[49]+4.0*f[35]+5.0*f[34]-4.0*f[33])+5.0*((-1.0*(6.708203932499369*f[32]+5.0*f[20]))+4.0*f[18]-5.0*f[17]))+2.0*(22.3606797749979*f[16]+3.0*(11.18033988749895*(f[1]-1.0*f[3])-15.0*f[7])+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(60.37383539249431*(f[105]+f[88])+9.0*(5.0*(f[85]-1.0*f[83])+4.0*f[63]+5.0*f[62]-1.0*(4.0*f[61]+5.0*f[60]))+6.708203932499369*((-5.0*f[50])+4.0*f[39]-5.0*f[38])+2.0*(13.41640786499874*f[37]+3.0*(3.0*(2.23606797749979*(f[9]-1.0*f[11])-3.0*f[23])+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(9.0*f[66]+6.708203932499369*f[42]-1.0*(6.708203932499369*f[40]+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[105]+f[88])+9.0*(5.0*(f[85]+f[83])-4.0*f[63]+5.0*f[62]-4.0*f[61]+5.0*f[60])+6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*f[38])-2.0*(13.41640786499874*f[37]+3.0*(3.0*(3.0*f[23]+2.23606797749979*(f[11]+f[9]))+5.0*f[4]))))/(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33])+5.0*(6.708203932499369*f[32]+5.0*f[20]-4.0*f[18]+5.0*f[17]))-2.0*(22.3606797749979*f[16]+3.0*(15.0*f[7]+11.18033988749895*(f[3]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33])+5.0*(6.708203932499369*f[32]+5.0*f[20]-4.0*f[18]+5.0*f[17]))-2.0*(22.3606797749979*f[16]+3.0*(15.0*f[7]+11.18033988749895*(f[3]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33])+5.0*(6.708203932499369*f[32]+5.0*f[20]-4.0*f[18]+5.0*f[17]))-2.0*(22.3606797749979*f[16]+3.0*(15.0*f[7]+11.18033988749895*(f[3]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(60.37383539249431*(f[105]+f[88])+9.0*(5.0*(f[85]+f[83])-4.0*f[63]+5.0*f[62]-4.0*f[61]+5.0*f[60])+6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*f[38])-2.0*(13.41640786499874*f[37]+3.0*(3.0*(3.0*f[23]+2.23606797749979*(f[11]+f[9]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(9.0*f[66]+6.708203932499369*(f[42]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(2.23606797749979*(45.0*(f[81]+f[57])+6.708203932499369*(5.0*(f[49]+f[47])-4.0*f[35]+5.0*f[34]-4.0*f[33])+5.0*(6.708203932499369*f[32]+5.0*f[20]-4.0*f[18]+5.0*f[17]))-2.0*(22.3606797749979*f[16]+3.0*(15.0*f[7]+11.18033988749895*(f[3]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(60.37383539249431*(f[105]+f[88])+9.0*(5.0*(f[85]+f[83])-4.0*f[63]+5.0*f[62]-4.0*f[61]+5.0*f[60])+6.708203932499369*(5.0*f[50]-4.0*f[39]+5.0*f[38])-2.0*(13.41640786499874*f[37]+3.0*(3.0*(3.0*f[23]+2.23606797749979*(f[11]+f[9]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(9.0*f[66]+6.708203932499369*(f[42]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]+60.37383539249431*(4.0*(f[98]-1.0*f[105])+5.0*f[97]-4.0*f[96]+5.0*(f[88]-1.0*f[95]))+9.0*(5.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62])-1.0*(5.0*(4.0*f[61]+5.0*f[60])+40.24922359499622*f[54]))+33.54101966249684*(4.0*(f[50]+f[39])-5.0*f[38])+2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*((-1.0*f[30])+f[28]-1.0*f[23])+11.18033988749895*(f[15]-1.0*f[11]+f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[72]-1.0*f[81])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]-1.0*f[57])))+6.708203932499369*(4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34])-1.0*(6.708203932499369*(4.0*f[33]+5.0*f[32])+54.0*f[26]-5.0*(4.0*(f[20]+f[18])-5.0*f[17])))+2.0*(22.3606797749979*f[16]+3.0*(15.0*((-1.0*f[14])+f[12]-1.0*f[7])+11.18033988749895*(f[5]-1.0*f[3]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[72]-1.0*f[81])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]-1.0*f[57])))+6.708203932499369*(4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34])-1.0*(6.708203932499369*(4.0*f[33]+5.0*f[32])+54.0*f[26]-5.0*(4.0*(f[20]+f[18])-5.0*f[17])))+2.0*(22.3606797749979*f[16]+3.0*(15.0*((-1.0*f[14])+f[12]-1.0*f[7])+11.18033988749895*(f[5]-1.0*f[3]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[72]-1.0*f[81])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]-1.0*f[57])))+6.708203932499369*(4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34])-1.0*(6.708203932499369*(4.0*f[33]+5.0*f[32])+54.0*f[26]-5.0*(4.0*(f[20]+f[18])-5.0*f[17])))+2.0*(22.3606797749979*f[16]+3.0*(15.0*((-1.0*f[14])+f[12]-1.0*f[7])+11.18033988749895*(f[5]-1.0*f[3]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.001666666666666667*(405.0*f[108]+60.37383539249431*(4.0*(f[98]-1.0*f[105])+5.0*f[97]-4.0*f[96]+5.0*(f[88]-1.0*f[95]))+9.0*(5.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62])-1.0*(5.0*(4.0*f[61]+5.0*f[60])+40.24922359499622*f[54]))+33.54101966249684*(4.0*(f[50]+f[39])-5.0*f[38])+2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*((-1.0*f[30])+f[28]-1.0*f[23])+11.18033988749895*(f[15]-1.0*f[11]+f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(60.37383539249431*f[101]+5.0*(9.0*(f[79]-1.0*f[77]+f[66])+6.708203932499369*(f[42]-1.0*f[46])-1.0*(6.708203932499369*f[40]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*(4.0*(f[72]-1.0*f[81])+5.0*f[71]-1.0*(4.0*f[70]+5.0*(f[69]-1.0*f[57])))+6.708203932499369*(4.0*((-1.0*f[49])+f[47]+f[45])-5.0*f[44]+4.0*(f[43]+f[35])+5.0*f[34])-1.0*(6.708203932499369*(4.0*f[33]+5.0*f[32])+54.0*f[26]-5.0*(4.0*(f[20]+f[18])-5.0*f[17])))+2.0*(22.3606797749979*f[16]+3.0*(15.0*((-1.0*f[14])+f[12]-1.0*f[7])+11.18033988749895*(f[5]-1.0*f[3]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.001666666666666667*(405.0*f[108]+60.37383539249431*(4.0*(f[98]-1.0*f[105])+5.0*f[97]-4.0*f[96]+5.0*(f[88]-1.0*f[95]))+9.0*(5.0*(4.0*((-1.0*f[85])+f[83]+f[76])-5.0*f[75]+4.0*(f[74]+f[63])+5.0*f[62])-1.0*(5.0*(4.0*f[61]+5.0*f[60])+40.24922359499622*f[54]))+33.54101966249684*(4.0*(f[50]+f[39])-5.0*f[38])+2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*((-1.0*f[30])+f[28]-1.0*f[23])+11.18033988749895*(f[15]-1.0*f[11]+f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(60.37383539249431*f[101]+5.0*(9.0*(f[79]-1.0*f[77]+f[66])+6.708203932499369*(f[42]-1.0*f[46])-1.0*(6.708203932499369*f[40]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[98]+f[95])+9.0*((-4.0*f[83])+5.0*(f[76]+f[75])-4.0*f[74]+5.0*(f[63]+f[60]))+6.708203932499369*(5.0*(f[39]+f[38])-4.0*f[50])-2.0*(13.41640786499874*f[37]+3.0*(3.0*(3.0*f[28]+2.23606797749979*(f[15]+f[9]))+5.0*f[4]))))/(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*((-4.0*f[47])+5.0*(f[45]+f[44])-4.0*f[43])+5.0*(6.708203932499369*(f[35]+f[32])-4.0*f[20]+5.0*(f[18]+f[17])))-2.0*(22.3606797749979*f[16]+3.0*(15.0*f[12]+11.18033988749895*(f[5]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*((-4.0*f[47])+5.0*(f[45]+f[44])-4.0*f[43])+5.0*(6.708203932499369*(f[35]+f[32])-4.0*f[20]+5.0*(f[18]+f[17])))-2.0*(22.3606797749979*f[16]+3.0*(15.0*f[12]+11.18033988749895*(f[5]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*((-4.0*f[47])+5.0*(f[45]+f[44])-4.0*f[43])+5.0*(6.708203932499369*(f[35]+f[32])-4.0*f[20]+5.0*(f[18]+f[17])))-2.0*(22.3606797749979*f[16]+3.0*(15.0*f[12]+11.18033988749895*(f[5]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.008333333333333333*(60.37383539249431*(f[98]+f[95])+9.0*((-4.0*f[83])+5.0*(f[76]+f[75])-4.0*f[74]+5.0*(f[63]+f[60]))+6.708203932499369*(5.0*(f[39]+f[38])-4.0*f[50])-2.0*(13.41640786499874*f[37]+3.0*(3.0*(3.0*f[28]+2.23606797749979*(f[15]+f[9]))+5.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.05*(9.0*f[77]+6.708203932499369*(f[46]+f[40])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(45.0*(f[72]+f[69])+6.708203932499369*((-4.0*f[47])+5.0*(f[45]+f[44])-4.0*f[43])+5.0*(6.708203932499369*(f[35]+f[32])-4.0*f[20]+5.0*(f[18]+f[17])))-2.0*(22.3606797749979*f[16]+3.0*(15.0*f[12]+11.18033988749895*(f[5]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.008333333333333333*(60.37383539249431*(f[98]+f[95])+9.0*((-4.0*f[83])+5.0*(f[76]+f[75])-4.0*f[74]+5.0*(f[63]+f[60]))+6.708203932499369*(5.0*(f[39]+f[38])-4.0*f[50])-2.0*(13.41640786499874*f[37]+3.0*(3.0*(3.0*f[28]+2.23606797749979*(f[15]+f[9]))+5.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.05*(9.0*f[77]+6.708203932499369*(f[46]+f[40])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[108]+60.37383539249431*((-4.0*(f[105]+f[98]))+5.0*f[97]-4.0*f[96]+5.0*(f[95]+f[88]))+9.0*(5.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])-40.24922359499622*f[54])+33.54101966249684*(5.0*f[38]-4.0*(f[50]+f[39]))-2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*(f[30]+f[28]+f[23])+11.18033988749895*(f[15]+f[11]+f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[92]+9.0*((-4.0*(f[81]+f[72]))+5.0*f[71]-4.0*f[70]+5.0*(f[69]+f[57]))+6.708203932499369*((-4.0*(f[49]+f[47]+f[45]))+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-54.0*f[26]+5.0*(5.0*f[17]-4.0*(f[20]+f[18])))-2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[14]+f[12]+f[7])+11.18033988749895*(f[5]+f[3]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*((-4.0*(f[81]+f[72]))+5.0*f[71]-4.0*f[70]+5.0*(f[69]+f[57]))+6.708203932499369*((-4.0*(f[49]+f[47]+f[45]))+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-54.0*f[26]+5.0*(5.0*f[17]-4.0*(f[20]+f[18])))-2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[14]+f[12]+f[7])+11.18033988749895*(f[5]+f[3]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*((-4.0*(f[81]+f[72]))+5.0*f[71]-4.0*f[70]+5.0*(f[69]+f[57]))+6.708203932499369*((-4.0*(f[49]+f[47]+f[45]))+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-54.0*f[26]+5.0*(5.0*f[17]-4.0*(f[20]+f[18])))-2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[14]+f[12]+f[7])+11.18033988749895*(f[5]+f[3]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.001666666666666667*(405.0*f[108]+60.37383539249431*((-4.0*(f[105]+f[98]))+5.0*f[97]-4.0*f[96]+5.0*(f[95]+f[88]))+9.0*(5.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])-40.24922359499622*f[54])+33.54101966249684*(5.0*f[38]-4.0*(f[50]+f[39]))-2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*(f[30]+f[28]+f[23])+11.18033988749895*(f[15]+f[11]+f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(60.37383539249431*f[101]+5.0*(9.0*(f[79]+f[77]+f[66])+6.708203932499369*(f[46]+f[42]+f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[92]+9.0*((-4.0*(f[81]+f[72]))+5.0*f[71]-4.0*f[70]+5.0*(f[69]+f[57]))+6.708203932499369*((-4.0*(f[49]+f[47]+f[45]))+5.0*f[44]-4.0*(f[43]+f[35])+5.0*f[34]-4.0*f[33]+5.0*f[32])-54.0*f[26]+5.0*(5.0*f[17]-4.0*(f[20]+f[18])))-2.0*(22.3606797749979*f[16]+3.0*(15.0*(f[14]+f[12]+f[7])+11.18033988749895*(f[5]+f[3]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.001666666666666667*(405.0*f[108]+60.37383539249431*((-4.0*(f[105]+f[98]))+5.0*f[97]-4.0*f[96]+5.0*(f[95]+f[88]))+9.0*(5.0*((-4.0*(f[85]+f[83]+f[76]))+5.0*f[75]-4.0*(f[74]+f[63])+5.0*f[62]-4.0*f[61]+5.0*f[60])-40.24922359499622*f[54])+33.54101966249684*(5.0*f[38]-4.0*(f[50]+f[39]))-2.0*(67.08203932499369*f[37]+3.0*(3.0*(15.0*(f[30]+f[28]+f[23])+11.18033988749895*(f[15]+f[11]+f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(60.37383539249431*f[101]+5.0*(9.0*(f[79]+f[77]+f[66])+6.708203932499369*(f[46]+f[42]+f[40])+5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[4][0] = 0.05555555555555555*(fReflXYZMuQuad[7][0]+8.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0]+8.0*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[4][1] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][0]-4188761.0*fReflXYZMuQuad[6][0])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])-4.63256860547201e+14*fReflXYZMuQuad[5][0])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][0]-1.0*(8377522.0*fReflXYZMuQuad[1][0]+2.266096151179001e+23*fReflXYZMuQuad[0][0])))); 
  fReflXYQuad[4][2] = 0.05555555555555555*(fReflXYZMuQuad[7][1]+8.0*fReflXYZMuQuad[6][1]+fReflXYZMuQuad[5][1]+8.0*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[4][3] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][0]+7.4121097687552e+14*fReflXYZMuQuad[6][0]+4.63256860547201e+14*fReflXYZMuQuad[5][0])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][0]+7.4121097687552e+14*fReflXYZMuQuad[1][0]+4.63256860547201e+14*fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[4][4] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][1]-4188761.0*fReflXYZMuQuad[6][1])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])-4.63256860547201e+14*fReflXYZMuQuad[5][1])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][1]-1.0*(8377522.0*fReflXYZMuQuad[1][1]+2.266096151179001e+23*fReflXYZMuQuad[0][1])))); 
  fReflXYQuad[4][5] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][0]+9.0*fReflXYZMuQuad[6][0])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][0]-1.346286087882789e+17*fReflXYZMuQuad[5][0])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][0]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])))); 
  fReflXYQuad[4][6] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][1]+7.4121097687552e+14*fReflXYZMuQuad[6][1]+4.63256860547201e+14*fReflXYZMuQuad[5][1])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][1]+7.4121097687552e+14*fReflXYZMuQuad[1][1]+4.63256860547201e+14*fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[4][7] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][0]-4.17103711299121e+14*fReflXYZMuQuad[6][0])+5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[5][0]+2384663.0*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])+3.190559553141742e+22*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[4][8] = 0.05555555555555555*(fReflXYZMuQuad[7][2]+8.0*fReflXYZMuQuad[6][2]+fReflXYZMuQuad[5][2]+8.0*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[4][9] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][0]+1.6692641e+7*fReflXYZMuQuad[6][0])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][0]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[2][0]+2384663.0*fReflXYZMuQuad[1][0]+3.190559553141742e+22*fReflXYZMuQuad[0][0])-3.20880527843592e+30*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])))); 
  fReflXYQuad[4][10] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][1]+9.0*fReflXYZMuQuad[6][1])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][1]-1.346286087882789e+17*fReflXYZMuQuad[5][1])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][1]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])))); 
  fReflXYQuad[4][11] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][1]-4.17103711299121e+14*fReflXYZMuQuad[6][1])+5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[5][1]+2384663.0*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])+3.190559553141743e+22*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[4][12] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][2]-4188761.0*fReflXYZMuQuad[6][2])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYZMuQuad[4][2]-4.63256860547201e+14*fReflXYZMuQuad[5][2])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][2]-1.0*(8377522.0*fReflXYZMuQuad[1][2]+2.266096151179001e+23*fReflXYZMuQuad[0][2])))); 
  fReflXYQuad[4][13] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][0]-2.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0])-27.0*fReflXYZMuQuad[3][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[2][0])+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[4][14] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][2]+7.4121097687552e+14*fReflXYZMuQuad[6][2]+4.63256860547201e+14*fReflXYZMuQuad[5][2])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][2]+7.4121097687552e+14*fReflXYZMuQuad[1][2]+4.63256860547201e+14*fReflXYZMuQuad[0][2]))); 
  fReflXYQuad[4][15] = 1.719407810605221e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][0]+27.0*fReflXYZMuQuad[6][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[5][0])+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[4][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[4][16] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][1]+1.6692641e+7*fReflXYZMuQuad[6][1])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][1]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[2][1]+2384663.0*fReflXYZMuQuad[1][1]+3.190559553141743e+22*fReflXYZMuQuad[0][1])-3.20880527843592e+30*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])))); 
  fReflXYQuad[4][17] = 1.719407810605222e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][1]+fReflXYZMuQuad[5][1])-27.0*fReflXYZMuQuad[3][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[2][1])+10.0*fReflXYZMuQuad[1][1]-5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[4][18] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][2]+9.0*fReflXYZMuQuad[6][2])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][2]-1.346286087882789e+17*fReflXYZMuQuad[5][2])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][2]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])))); 
  fReflXYQuad[4][19] = 1.719407810605222e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][1]+27.0*fReflXYZMuQuad[6][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[5][1])+2.0*(5.0*fReflXYZMuQuad[3][1]-5.0*fReflXYZMuQuad[4][1])+5.0*fReflXYZMuQuad[2][1])); 
  } 

 
// node (x,y)_6 
  vcutSq_i = -(0.01*q_*(zVal*((426.9074841227313*phiWall[19]-426.9074841227313*phi[19]-318.1980515339465*phiWall[16]+318.1980515339465*(phi[16]+phiWall[15])-318.1980515339465*phi[15]-237.1708245126285*phiWall[9]+237.1708245126285*phi[9])*zVal+146.9693845669907*phiWall[18]-146.9693845669907*(phi[18]+phiWall[17])+146.9693845669907*phi[17]-109.5445115010333*phiWall[14]+109.5445115010333*phi[14]-109.5445115010333*phiWall[13]+109.5445115010333*phi[13]+220.454076850486*phiWall[10]-220.454076850486*phi[10]-164.3167672515499*phiWall[6]+164.3167672515499*(phi[6]+phiWall[5])-164.3167672515499*phi[5]-122.4744871391589*phiWall[3]+122.4744871391589*phi[3])-142.3024947075771*phiWall[19]+142.3024947075771*phi[19]+106.0660171779822*phiWall[16]-106.0660171779822*(phi[16]+phiWall[15])+106.0660171779822*phi[15]+84.85281374238573*phiWall[12]-84.85281374238573*(phi[12]+phiWall[11])+84.85281374238573*phi[11]+79.0569415042095*phiWall[9]-79.0569415042095*phi[9]-63.24555320336762*phiWall[8]+63.24555320336762*phi[8]-63.24555320336762*phiWall[7]+63.24555320336762*phi[7]+127.2792206135786*phiWall[4]-127.2792206135786*phi[4]-94.86832980505142*phiWall[2]+94.86832980505142*(phi[2]+phiWall[1])-94.86832980505142*phi[1]-70.71067811865477*phiWall[0]+70.71067811865477*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[5][0] = 0.0; 
  fReflXYQuad[5][1] = 0.0; 
  fReflXYQuad[5][2] = 0.0; 
  fReflXYQuad[5][3] = 0.0; 
  fReflXYQuad[5][4] = 0.0; 
  fReflXYQuad[5][5] = 0.0; 
  fReflXYQuad[5][6] = 0.0; 
  fReflXYQuad[5][7] = 0.0; 
  fReflXYQuad[5][8] = 0.0; 
  fReflXYQuad[5][9] = 0.0; 
  fReflXYQuad[5][10] = 0.0; 
  fReflXYQuad[5][11] = 0.0; 
  fReflXYQuad[5][12] = 0.0; 
  fReflXYQuad[5][13] = 0.0; 
  fReflXYQuad[5][14] = 0.0; 
  fReflXYQuad[5][15] = 0.0; 
  fReflXYQuad[5][16] = 0.0; 
  fReflXYQuad[5][17] = 0.0; 
  fReflXYQuad[5][18] = 0.0; 
  fReflXYQuad[5][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[5][0] = -0.02*(4.47213595499958*(6.708203932499369*f[32]-1.0*(6.708203932499369*f[31]+5.0*(f[17]+f[16])))+3.0*(15.0*f[6]+11.18033988749895*(f[1]-1.0*f[2]))-25.0*f[0]); 
  fReflXYQuad[5][1] = -0.03333333333333333*(2.0*(9.0*f[57]-1.0*(9.0*f[56]+6.708203932499369*(f[34]+f[33])))+3.0*(3.0*(3.0*f[21]+2.23606797749979*(f[7]-1.0*f[8]))-5.0*f[3])); 
  fReflXYQuad[5][2] = -0.03333333333333333*(2.0*(9.0*f[60]-1.0*(9.0*f[59]+6.708203932499369*(f[38]+f[37])))+3.0*(3.0*(3.0*f[22]+2.23606797749979*(f[9]-1.0*f[10]))-5.0*f[4])); 
  fReflXYQuad[5][3] = -0.03333333333333333*(2.0*(9.0*f[69]-1.0*(9.0*f[68]+6.708203932499369*(f[44]+f[43])))+3.0*(3.0*(3.0*f[25]+2.23606797749979*(f[12]-1.0*f[13]))-5.0*f[5])); 
  fReflXYQuad[5][4] = -0.02*(4.47213595499958*(6.708203932499369*f[88]-1.0*(6.708203932499369*f[87]+5.0*(f[62]+f[61])))+3.0*(15.0*f[51]+11.18033988749895*(f[23]-1.0*f[24]))-25.0*f[11]); 
  fReflXYQuad[5][5] = -0.02*(4.47213595499958*(6.708203932499369*f[92]-1.0*(6.708203932499369*f[91]+5.0*(f[71]+f[70])))+3.0*(15.0*f[52]+11.18033988749895*(f[26]-1.0*f[27]))-25.0*f[14]); 
  fReflXYQuad[5][6] = -0.02*(4.47213595499958*(6.708203932499369*f[95]-1.0*(6.708203932499369*f[94]+5.0*(f[75]+f[74])))+3.0*(15.0*f[53]+11.18033988749895*(f[28]-1.0*f[29]))-25.0*f[15]); 
  fReflXYQuad[5][7] = -0.1*(9.0*f[58]+6.708203932499369*(f[35]-1.0*f[36])-5.0*f[18]); 
  fReflXYQuad[5][8] = -0.1*(9.0*f[65]+6.708203932499369*(f[40]-1.0*f[41])-5.0*f[19]); 
  fReflXYQuad[5][9] = -0.1*(9.0*f[80]+6.708203932499369*(f[47]-1.0*f[48])-5.0*f[20]); 
  fReflXYQuad[5][10] = -0.03333333333333333*(2.0*(9.0*f[108]-1.0*(9.0*f[107]+6.708203932499369*(f[97]+f[96])))+3.0*(3.0*(3.0*f[86]+2.23606797749979*(f[54]-1.0*f[55]))-5.0*f[30])); 
  fReflXYQuad[5][11] = -0.1*(9.0*f[89]+6.708203932499369*(f[63]-1.0*f[64])-5.0*f[39]); 
  fReflXYQuad[5][12] = -0.1*(9.0*f[90]+6.708203932499369*(f[66]-1.0*f[67])-5.0*f[42]); 
  fReflXYQuad[5][13] = -0.1*(9.0*f[93]+6.708203932499369*(f[72]-1.0*f[73])-5.0*f[45]); 
  fReflXYQuad[5][14] = -0.1*(9.0*f[100]+6.708203932499369*(f[77]-1.0*f[78])-5.0*f[46]); 
  fReflXYQuad[5][15] = -0.1*(9.0*f[103]+6.708203932499369*(f[81]-1.0*f[82])-5.0*f[49]); 
  fReflXYQuad[5][16] = -0.1*(9.0*f[104]+6.708203932499369*(f[83]-1.0*f[84])-5.0*f[50]); 
  fReflXYQuad[5][17] = -0.1*(9.0*f[109]+6.708203932499369*(f[98]-1.0*f[99])-5.0*f[76]); 
  fReflXYQuad[5][18] = -0.1*(9.0*f[110]+6.708203932499369*(f[101]-1.0*f[102])-5.0*f[79]); 
  fReflXYQuad[5][19] = -0.1*(9.0*f[111]+6.708203932499369*(f[105]-1.0*f[106])-5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])+60.37383539249431*((-1.0*f[106])+f[105]-1.0*(f[104]+f[99]-1.0*(f[98]+f[97]))+f[96]+f[95]-1.0*(f[94]+f[89]-1.0*f[88]+f[87])))+9.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+20.12461179749811*(f[55]-1.0*f[54]+f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]-1.0*f[29]+f[28]-1.0*f[24]+f[23]-1.0*f[22])+11.18033988749895*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9]))+25.0*f[4])))/(269.9999999999999*(f[103]+f[93]-1.0*f[92]+f[91])+9.0*(22.3606797749979*((-1.0*f[82])+f[81]-1.0*(f[80]+f[73]-1.0*(f[72]+f[71]))+f[70]+f[69]-1.0*(f[68]+f[58]-1.0*f[57]))-1.0*(22.3606797749979*f[56]+45.0*f[52]))+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*(f[27]-1.0*f[26]+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*f[13]+f[12]-1.0*f[8]+f[7]-1.0*f[6])+55.90169943749476*((-1.0*(f[5]+f[3]))+f[2]-1.0*f[1]))+125.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(269.9999999999999*(f[103]+f[93]-1.0*f[92]+f[91])+9.0*(22.3606797749979*((-1.0*f[82])+f[81]-1.0*(f[80]+f[73]-1.0*(f[72]+f[71]))+f[70]+f[69]-1.0*(f[68]+f[58]-1.0*f[57]))-1.0*(22.3606797749979*f[56]+45.0*f[52]))+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*(f[27]-1.0*f[26]+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*f[13]+f[12]-1.0*f[8]+f[7]-1.0*f[6])+55.90169943749476*((-1.0*(f[5]+f[3]))+f[2]-1.0*f[1]))+125.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.002*(269.9999999999999*(f[103]+f[93]-1.0*f[92]+f[91])+9.0*(22.3606797749979*((-1.0*f[82])+f[81]-1.0*(f[80]+f[73]-1.0*(f[72]+f[71]))+f[70]+f[69]-1.0*(f[68]+f[58]-1.0*f[57]))-1.0*(22.3606797749979*f[56]+45.0*f[52]))+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*(f[27]-1.0*f[26]+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*f[13]+f[12]-1.0*f[8]+f[7]-1.0*f[6])+55.90169943749476*((-1.0*(f[5]+f[3]))+f[2]-1.0*f[1]))+125.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])+60.37383539249431*((-1.0*f[106])+f[105]-1.0*(f[104]+f[99]-1.0*(f[98]+f[97]))+f[96]+f[95]-1.0*(f[94]+f[89]-1.0*f[88]+f[87])))+9.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+20.12461179749811*(f[55]-1.0*f[54]+f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]-1.0*f[29]+f[28]-1.0*f[24]+f[23]-1.0*f[22])+11.18033988749895*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9]))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(81.0*f[110]+60.37383539249431*((-1.0*f[102])+f[101]-1.0*(f[100]+f[90]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[77]+f[67]-1.0*f[66]+f[65])+6.708203932499369*(f[46]+f[42]-1.0*f[41]+f[40])-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (0.002*(269.9999999999999*(f[103]+f[93]-1.0*f[92]+f[91])+9.0*(22.3606797749979*((-1.0*f[82])+f[81]-1.0*(f[80]+f[73]-1.0*(f[72]+f[71]))+f[70]+f[69]-1.0*(f[68]+f[58]-1.0*f[57]))-1.0*(22.3606797749979*f[56]+45.0*f[52]))+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*(f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*(f[27]-1.0*f[26]+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*f[13]+f[12]-1.0*f[8]+f[7]-1.0*f[6])+55.90169943749476*((-1.0*(f[5]+f[3]))+f[2]-1.0*f[1]))+125.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*f[108]+f[107])+60.37383539249431*((-1.0*f[106])+f[105]-1.0*(f[104]+f[99]-1.0*(f[98]+f[97]))+f[96]+f[95]-1.0*(f[94]+f[89]-1.0*f[88]+f[87])))+9.0*((-27.0*f[86])+10.0*((-1.0*f[85])+f[84]-1.0*(f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+20.12461179749811*(f[55]-1.0*f[54]+f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]-1.0*f[29]+f[28]-1.0*f[24]+f[23]-1.0*f[22])+11.18033988749895*((-1.0*(f[15]+f[11]))+f[10]-1.0*f[9]))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (-0.01*(81.0*f[110]+60.37383539249431*((-1.0*f[102])+f[101]-1.0*(f[100]+f[90]))+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*f[77]+f[67]-1.0*f[66]+f[65])+6.708203932499369*(f[46]+f[42]-1.0*f[41]+f[40])-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+60.37383539249431*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95])-5.0*f[89])+9.0*(5.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(10.0*f[60]-1.0*(10.0*f[59]+20.12461179749811*f[53])))+33.54101966249684*(5.0*f[39]-4.0*f[50])+2.0*(3.0*(3.0*(15.0*(f[29]-1.0*f[28]+f[22])+11.18033988749895*(f[15]-1.0*f[10]+f[9]))-25.0*f[4])-67.08203932499369*(f[38]+f[37]))))/(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]-5.0*f[73]+5.0*f[72]+4.0*(f[68]-1.0*f[69])-5.0*f[58])+6.708203932499369*(4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*f[36])+2.0*(13.41640786499874*f[32]-1.0*(13.41640786499874*f[31]+27.0*f[25]))+5.0*(5.0*f[18]-4.0*f[20]))+2.0*((-22.3606797749979*(f[17]+f[16]))+3.0*(15.0*(f[13]-1.0*f[12]+f[6])+11.18033988749895*(f[5]-1.0*f[2]+f[1]))-25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]-5.0*f[73]+5.0*f[72]+4.0*(f[68]-1.0*f[69])-5.0*f[58])+6.708203932499369*(4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*f[36])+2.0*(13.41640786499874*f[32]-1.0*(13.41640786499874*f[31]+27.0*f[25]))+5.0*(5.0*f[18]-4.0*f[20]))+2.0*((-22.3606797749979*(f[17]+f[16]))+3.0*(15.0*(f[13]-1.0*f[12]+f[6])+11.18033988749895*(f[5]-1.0*f[2]+f[1]))-25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]-5.0*f[73]+5.0*f[72]+4.0*(f[68]-1.0*f[69])-5.0*f[58])+6.708203932499369*(4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*f[36])+2.0*(13.41640786499874*f[32]-1.0*(13.41640786499874*f[31]+27.0*f[25]))+5.0*(5.0*f[18]-4.0*f[20]))+2.0*((-22.3606797749979*(f[17]+f[16]))+3.0*(15.0*(f[13]-1.0*f[12]+f[6])+11.18033988749895*(f[5]-1.0*f[2]+f[1]))-25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.001666666666666667*(405.0*f[109]+60.37383539249431*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95])-5.0*f[89])+9.0*(5.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(10.0*f[60]-1.0*(10.0*f[59]+20.12461179749811*f[53])))+33.54101966249684*(5.0*f[39]-4.0*f[50])+2.0*(3.0*(3.0*(15.0*(f[29]-1.0*f[28]+f[22])+11.18033988749895*(f[15]-1.0*f[10]+f[9]))-25.0*f[4])-67.08203932499369*(f[38]+f[37]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01*(60.37383539249431*f[100]+5.0*(9.0*((-1.0*f[78])+f[77]-1.0*f[65])+6.708203932499369*((-1.0*f[46])+f[41]-1.0*f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]-5.0*f[73]+5.0*f[72]+4.0*(f[68]-1.0*f[69])-5.0*f[58])+6.708203932499369*(4.0*(f[47]-1.0*f[48])-5.0*f[45]+4.0*(f[44]+f[43])+5.0*f[36])+2.0*(13.41640786499874*f[32]-1.0*(13.41640786499874*f[31]+27.0*f[25]))+5.0*(5.0*f[18]-4.0*f[20]))+2.0*((-22.3606797749979*(f[17]+f[16]))+3.0*(15.0*(f[13]-1.0*f[12]+f[6])+11.18033988749895*(f[5]-1.0*f[2]+f[1]))-25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (-0.001666666666666667*(405.0*f[109]+60.37383539249431*(4.0*f[104]+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95])-5.0*f[89])+9.0*(5.0*(4.0*(f[83]-1.0*f[84])-5.0*f[76]+4.0*(f[75]+f[74])+5.0*(f[64]-1.0*f[63]))+2.0*(10.0*f[60]-1.0*(10.0*f[59]+20.12461179749811*f[53])))+33.54101966249684*(5.0*f[39]-4.0*f[50])+2.0*(3.0*(3.0*(15.0*(f[29]-1.0*f[28]+f[22])+11.18033988749895*(f[15]-1.0*f[10]+f[9]))-25.0*f[4])-67.08203932499369*(f[38]+f[37]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (0.01*(60.37383539249431*f[100]+5.0*(9.0*((-1.0*f[78])+f[77]-1.0*f[65])+6.708203932499369*((-1.0*f[46])+f[41]-1.0*f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))+60.37383539249431*((-1.0*f[106])+f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+9.0*((-27.0*f[86])+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+20.12461179749811*(f[55]-1.0*(f[54]+f[53]-1.0*f[51])))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]+f[29]-1.0*(f[28]+f[24]-1.0*(f[23]+f[22])))+11.18033988749895*(f[15]-1.0*(f[11]+f[10]-1.0*f[9])))-25.0*f[4])))/(269.9999999999999*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+9.0*(22.3606797749979*((-1.0*f[82])+f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57])-1.0*(22.3606797749979*f[56]+45.0*f[52]))+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31]))+27.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21]))-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]+f[13]-1.0*(f[12]+f[8]-1.0*(f[7]+f[6])))+55.90169943749476*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))-125.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(269.9999999999999*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+9.0*(22.3606797749979*((-1.0*f[82])+f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57])-1.0*(22.3606797749979*f[56]+45.0*f[52]))+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31]))+27.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21]))-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]+f[13]-1.0*(f[12]+f[8]-1.0*(f[7]+f[6])))+55.90169943749476*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))-125.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.002*(269.9999999999999*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+9.0*(22.3606797749979*((-1.0*f[82])+f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57])-1.0*(22.3606797749979*f[56]+45.0*f[52]))+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31]))+27.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21]))-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]+f[13]-1.0*(f[12]+f[8]-1.0*(f[7]+f[6])))+55.90169943749476*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))-125.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))+60.37383539249431*((-1.0*f[106])+f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+9.0*((-27.0*f[86])+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+20.12461179749811*(f[55]-1.0*(f[54]+f[53]-1.0*f[51])))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]+f[29]-1.0*(f[28]+f[24]-1.0*(f[23]+f[22])))+11.18033988749895*(f[15]-1.0*(f[11]+f[10]-1.0*f[9])))-25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(81.0*f[110]+60.37383539249431*((-1.0*f[102])+f[101]+f[100]-1.0*f[90])+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]+f[67]-1.0*(f[66]+f[65]))+6.708203932499369*((-1.0*f[46])+f[42]+f[41]-1.0*f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (-0.002*(269.9999999999999*(f[103]-1.0*(f[93]+f[92]-1.0*f[91]))+9.0*(22.3606797749979*((-1.0*f[82])+f[81]+f[80]+f[73]-1.0*f[72]+f[71]+f[70]-1.0*f[69]+f[68]+f[58]+f[57])-1.0*(22.3606797749979*f[56]+45.0*f[52]))+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]+f[45]+f[44]+f[43]-1.0*f[36]+f[35]-1.0*(f[34]+f[33]-1.0*f[32]+f[31]))+27.0*(f[27]-1.0*(f[26]+f[25]-1.0*f[21]))-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]+f[13]-1.0*(f[12]+f[8]-1.0*(f[7]+f[6])))+55.90169943749476*(f[5]-1.0*(f[3]+f[2]-1.0*f[1])))-125.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]-1.0*f[107]))+60.37383539249431*((-1.0*f[106])+f[105]+f[104]+f[99]-1.0*f[98]+f[97]+f[96]-1.0*f[95]+f[94]+f[89]+f[88]-1.0*f[87]))+9.0*((-27.0*f[86])+10.0*((-1.0*(f[85]+f[84]))+f[83]+f[76]+f[75]+f[74]-1.0*f[64]+f[63]-1.0*(f[62]+f[61]-1.0*f[60]+f[59]))+20.12461179749811*(f[55]-1.0*(f[54]+f[53]-1.0*f[51])))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]+f[29]-1.0*(f[28]+f[24]-1.0*(f[23]+f[22])))+11.18033988749895*(f[15]-1.0*(f[11]+f[10]-1.0*f[9])))-25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (0.01*(81.0*f[110]+60.37383539249431*((-1.0*f[102])+f[101]+f[100]-1.0*f[90])+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]+f[67]-1.0*(f[66]+f[65]))+6.708203932499369*((-1.0*f[46])+f[42]+f[41]-1.0*f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+60.37383539249431*(5.0*((-1.0*f[106])+f[105]-1.0*f[104])+4.0*(f[89]-1.0*f[88]+f[87]))+9.0*(25.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+20.12461179749811*f[51])))+167.7050983124842*f[50]+2.0*(3.0*(3.0*(15.0*(f[24]-1.0*f[23]+f[22])+11.18033988749895*(f[11]-1.0*f[10]+f[9]))-25.0*f[4])-67.08203932499369*(f[39]+f[38]+f[37]))))/(2.23606797749979*(60.37383539249431*f[103]+9.0*((-5.0*f[82])+5.0*f[81]+4.0*(f[58]-1.0*f[57]+f[56]))+6.708203932499369*(5.0*f[48]-5.0*f[49])+2.0*(13.41640786499874*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32])-1.0*(13.41640786499874*f[31]+27.0*f[21]))+25.0*f[20])+2.0*((-22.3606797749979*(f[18]+f[17]+f[16]))+3.0*(15.0*(f[8]-1.0*f[7]+f[6])+11.18033988749895*(f[3]-1.0*f[2]+f[1]))-25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*((-5.0*f[82])+5.0*f[81]+4.0*(f[58]-1.0*f[57]+f[56]))+6.708203932499369*(5.0*f[48]-5.0*f[49])+2.0*(13.41640786499874*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32])-1.0*(13.41640786499874*f[31]+27.0*f[21]))+25.0*f[20])+2.0*((-22.3606797749979*(f[18]+f[17]+f[16]))+3.0*(15.0*(f[8]-1.0*f[7]+f[6])+11.18033988749895*(f[3]-1.0*f[2]+f[1]))-25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*((-5.0*f[82])+5.0*f[81]+4.0*(f[58]-1.0*f[57]+f[56]))+6.708203932499369*(5.0*f[48]-5.0*f[49])+2.0*(13.41640786499874*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32])-1.0*(13.41640786499874*f[31]+27.0*f[21]))+25.0*f[20])+2.0*((-22.3606797749979*(f[18]+f[17]+f[16]))+3.0*(15.0*(f[8]-1.0*f[7]+f[6])+11.18033988749895*(f[3]-1.0*f[2]+f[1]))-25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*((-1.0*f[106])+f[105]-1.0*f[104])+4.0*(f[89]-1.0*f[88]+f[87]))+9.0*(25.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+20.12461179749811*f[51])))+167.7050983124842*f[50]+2.0*(3.0*(3.0*(15.0*(f[24]-1.0*f[23]+f[22])+11.18033988749895*(f[11]-1.0*f[10]+f[9]))-25.0*f[4])-67.08203932499369*(f[39]+f[38]+f[37]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01*(60.37383539249431*f[90]+5.0*(9.0*((-1.0*f[67])+f[66]-1.0*f[65])+6.708203932499369*((-1.0*f[42])+f[41]-1.0*f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*((-5.0*f[82])+5.0*f[81]+4.0*(f[58]-1.0*f[57]+f[56]))+6.708203932499369*(5.0*f[48]-5.0*f[49])+2.0*(13.41640786499874*((-1.0*f[36])+f[35]+f[34]+f[33]+f[32])-1.0*(13.41640786499874*f[31]+27.0*f[21]))+25.0*f[20])+2.0*((-22.3606797749979*(f[18]+f[17]+f[16]))+3.0*(15.0*(f[8]-1.0*f[7]+f[6])+11.18033988749895*(f[3]-1.0*f[2]+f[1]))-25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (-0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*((-1.0*f[106])+f[105]-1.0*f[104])+4.0*(f[89]-1.0*f[88]+f[87]))+9.0*(25.0*((-1.0*f[85])+f[84]-1.0*f[83])+2.0*(10.0*((-1.0*f[64])+f[63]+f[62]+f[61]+f[60])-1.0*(10.0*f[59]+20.12461179749811*f[51])))+167.7050983124842*f[50]+2.0*(3.0*(3.0*(15.0*(f[24]-1.0*f[23]+f[22])+11.18033988749895*(f[11]-1.0*f[10]+f[9]))-25.0*f[4])-67.08203932499369*(f[39]+f[38]+f[37]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (0.01*(60.37383539249431*f[90]+5.0*(9.0*((-1.0*f[67])+f[66]-1.0*f[65])+6.708203932499369*((-1.0*f[42])+f[41]-1.0*f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+60.37383539249431*(5.0*((-1.0*f[106])+f[105]+f[104])+4.0*(f[87]-1.0*(f[89]+f[88])))+45.0*(5.0*(f[83]-1.0*(f[85]+f[84]))+4.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59]))-1.0*(362.243012354966*f[51]+167.7050983124842*f[50])+2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[24]-1.0*(f[23]+f[22]))+11.18033988749895*(f[11]+f[10]-1.0*f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[103]+9.0*((-5.0*f[82])+5.0*(f[81]+f[80])+4.0*(f[56]-1.0*(f[58]+f[57])))+6.708203932499369*((-5.0*(f[49]+f[48]))+5.0*f[47]+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))-1.0*(54.0*f[21]+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[8]-1.0*(f[7]+f[6]))+11.18033988749895*(f[3]+f[2]-1.0*f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*((-5.0*f[82])+5.0*(f[81]+f[80])+4.0*(f[56]-1.0*(f[58]+f[57])))+6.708203932499369*((-5.0*(f[49]+f[48]))+5.0*f[47]+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))-1.0*(54.0*f[21]+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[8]-1.0*(f[7]+f[6]))+11.18033988749895*(f[3]+f[2]-1.0*f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*((-5.0*f[82])+5.0*(f[81]+f[80])+4.0*(f[56]-1.0*(f[58]+f[57])))+6.708203932499369*((-5.0*(f[49]+f[48]))+5.0*f[47]+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))-1.0*(54.0*f[21]+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[8]-1.0*(f[7]+f[6]))+11.18033988749895*(f[3]+f[2]-1.0*f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*((-1.0*f[106])+f[105]+f[104])+4.0*(f[87]-1.0*(f[89]+f[88])))+45.0*(5.0*(f[83]-1.0*(f[85]+f[84]))+4.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59]))-1.0*(362.243012354966*f[51]+167.7050983124842*f[50])+2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[24]-1.0*(f[23]+f[22]))+11.18033988749895*(f[11]+f[10]-1.0*f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.01*(60.37383539249431*f[90]+5.0*(9.0*((-1.0*f[67])+f[66]+f[65])+6.708203932499369*(f[40]-1.0*(f[42]+f[41]))-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*((-5.0*f[82])+5.0*(f[81]+f[80])+4.0*(f[56]-1.0*(f[58]+f[57])))+6.708203932499369*((-5.0*(f[49]+f[48]))+5.0*f[47]+4.0*(f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))-1.0*(54.0*f[21]+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[8]-1.0*(f[7]+f[6]))+11.18033988749895*(f[3]+f[2]-1.0*f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*((-1.0*f[106])+f[105]+f[104])+4.0*(f[87]-1.0*(f[89]+f[88])))+45.0*(5.0*(f[83]-1.0*(f[85]+f[84]))+4.0*(f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59]))-1.0*(362.243012354966*f[51]+167.7050983124842*f[50])+2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[24]-1.0*(f[23]+f[22]))+11.18033988749895*(f[11]+f[10]-1.0*f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (-0.01*(60.37383539249431*f[90]+5.0*(9.0*((-1.0*f[67])+f[66]+f[65])+6.708203932499369*(f[40]-1.0*(f[42]+f[41]))-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]-1.0*f[107])+60.37383539249431*((-1.0*f[106])+f[105]-1.0*f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])))+9.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64]-1.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+20.12461179749811*((-1.0*f[55])+f[54]-1.0*f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[30])+f[29]-1.0*(f[28]+f[24]-1.0*f[23]+f[22]))+11.18033988749895*(f[15]-1.0*f[11]+f[10]-1.0*f[9]))+25.0*f[4])))/(269.9999999999999*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*((-1.0*f[82])+f[81]-1.0*f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]-1.0*(f[12]+f[8]-1.0*f[7]+f[6]))+55.90169943749476*(f[5]-1.0*f[3]+f[2]-1.0*f[1]))+125.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(269.9999999999999*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*((-1.0*f[82])+f[81]-1.0*f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]-1.0*(f[12]+f[8]-1.0*f[7]+f[6]))+55.90169943749476*(f[5]-1.0*f[3]+f[2]-1.0*f[1]))+125.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.002*(269.9999999999999*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*((-1.0*f[82])+f[81]-1.0*f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]-1.0*(f[12]+f[8]-1.0*f[7]+f[6]))+55.90169943749476*(f[5]-1.0*f[3]+f[2]-1.0*f[1]))+125.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]-1.0*f[107])+60.37383539249431*((-1.0*f[106])+f[105]-1.0*f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])))+9.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64]-1.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+20.12461179749811*((-1.0*f[55])+f[54]-1.0*f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[30])+f[29]-1.0*(f[28]+f[24]-1.0*f[23]+f[22]))+11.18033988749895*(f[15]-1.0*f[11]+f[10]-1.0*f[9]))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(81.0*f[110]+60.37383539249431*((-1.0*f[102])+f[101]-1.0*f[100]+f[90])+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*(f[77]+f[67]-1.0*f[66]+f[65]))+6.708203932499369*(f[46]-1.0*f[42]+f[41]-1.0*f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.002*(269.9999999999999*(f[103]-1.0*f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*((-1.0*f[82])+f[81]-1.0*f[80]+f[73]-1.0*(f[72]+f[71]+f[70]+f[69]-1.0*f[68]+f[58]-1.0*f[57]+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*f[49])+f[48]-1.0*f[47]+f[45]+f[44]+f[43]+f[36]-1.0*(f[35]+f[34]+f[33]+f[32]-1.0*f[31]))+27.0*((-1.0*f[27])+f[26]-1.0*f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]-1.0*(f[12]+f[8]-1.0*f[7]+f[6]))+55.90169943749476*(f[5]-1.0*f[3]+f[2]-1.0*f[1]))+125.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]-1.0*f[107])+60.37383539249431*((-1.0*f[106])+f[105]-1.0*f[104]+f[99]-1.0*(f[98]+f[97]+f[96]+f[95]-1.0*f[94]+f[89]-1.0*f[88]+f[87])))+9.0*(27.0*f[86]+10.0*((-1.0*f[85])+f[84]-1.0*f[83]+f[76]+f[75]+f[74]+f[64]-1.0*(f[63]+f[62]+f[61]+f[60]-1.0*f[59]))+20.12461179749811*((-1.0*f[55])+f[54]-1.0*f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[30])+f[29]-1.0*(f[28]+f[24]-1.0*f[23]+f[22]))+11.18033988749895*(f[15]-1.0*f[11]+f[10]-1.0*f[9]))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (0.01*(81.0*f[110]+60.37383539249431*((-1.0*f[102])+f[101]-1.0*f[100]+f[90])+5.0*(9.0*((-1.0*f[79])+f[78]-1.0*(f[77]+f[67]-1.0*f[66]+f[65]))+6.708203932499369*(f[46]-1.0*f[42]+f[41]-1.0*f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+60.37383539249431*((-4.0*f[104])+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95])+5.0*f[89])+9.0*(5.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(10.0*(f[59]-1.0*f[60])-20.12461179749811*f[53]))+33.54101966249684*(4.0*f[50]-5.0*f[39])+2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*(f[29]-1.0*(f[28]+f[22]))+11.18033988749895*(f[15]+f[10]-1.0*f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[93]+9.0*((-1.0*(4.0*f[80]+5.0*f[73]))+5.0*f[72]+4.0*(f[68]-1.0*f[69])+5.0*f[58])+6.708203932499369*(4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45]-4.0*(f[44]+f[43]))-5.0*f[36]+5.0*f[35])+2.0*(13.41640786499874*(f[31]-1.0*f[32])-27.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[13]-1.0*(f[12]+f[6]))+11.18033988749895*(f[5]+f[2]-1.0*f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*((-1.0*(4.0*f[80]+5.0*f[73]))+5.0*f[72]+4.0*(f[68]-1.0*f[69])+5.0*f[58])+6.708203932499369*(4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45]-4.0*(f[44]+f[43]))-5.0*f[36]+5.0*f[35])+2.0*(13.41640786499874*(f[31]-1.0*f[32])-27.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[13]-1.0*(f[12]+f[6]))+11.18033988749895*(f[5]+f[2]-1.0*f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*((-1.0*(4.0*f[80]+5.0*f[73]))+5.0*f[72]+4.0*(f[68]-1.0*f[69])+5.0*f[58])+6.708203932499369*(4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45]-4.0*(f[44]+f[43]))-5.0*f[36]+5.0*f[35])+2.0*(13.41640786499874*(f[31]-1.0*f[32])-27.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[13]-1.0*(f[12]+f[6]))+11.18033988749895*(f[5]+f[2]-1.0*f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.001666666666666667*(405.0*f[109]+60.37383539249431*((-4.0*f[104])+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95])+5.0*f[89])+9.0*(5.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(10.0*(f[59]-1.0*f[60])-20.12461179749811*f[53]))+33.54101966249684*(4.0*f[50]-5.0*f[39])+2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*(f[29]-1.0*(f[28]+f[22]))+11.18033988749895*(f[15]+f[10]-1.0*f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(60.37383539249431*f[100]+5.0*(9.0*((-1.0*f[78])+f[77]+f[65])+6.708203932499369*(f[40]-1.0*(f[46]+f[41]))-5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*((-1.0*(4.0*f[80]+5.0*f[73]))+5.0*f[72]+4.0*(f[68]-1.0*f[69])+5.0*f[58])+6.708203932499369*(4.0*f[48]-1.0*(4.0*f[47]+5.0*f[45]-4.0*(f[44]+f[43]))-5.0*f[36]+5.0*f[35])+2.0*(13.41640786499874*(f[31]-1.0*f[32])-27.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[13]-1.0*(f[12]+f[6]))+11.18033988749895*(f[5]+f[2]-1.0*f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (0.001666666666666667*(405.0*f[109]+60.37383539249431*((-4.0*f[104])+5.0*(f[98]-1.0*f[99])+4.0*(f[94]-1.0*f[95])+5.0*f[89])+9.0*(5.0*(4.0*f[84]-1.0*(4.0*f[83]+5.0*f[76]-4.0*(f[75]+f[74]))+5.0*(f[63]-1.0*f[64]))+2.0*(10.0*(f[59]-1.0*f[60])-20.12461179749811*f[53]))+33.54101966249684*(4.0*f[50]-5.0*f[39])+2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*(f[29]-1.0*(f[28]+f[22]))+11.18033988749895*(f[15]+f[10]-1.0*f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (-0.01*(60.37383539249431*f[100]+5.0*(9.0*((-1.0*f[78])+f[77]+f[65])+6.708203932499369*(f[40]-1.0*(f[46]+f[41]))-5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]+f[108]-1.0*f[107])+60.37383539249431*((-1.0*f[106])+f[105]+f[104]-1.0*f[99]+f[98]-1.0*(f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88])+f[87])))+9.0*(27.0*f[86]+10.0*((-1.0*(f[85]+f[84]))+f[83]-1.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59]))+20.12461179749811*((-1.0*f[55])+f[54]+f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(45.0*((-1.0*(f[30]+f[29]))+f[28]-1.0*f[24]+f[23]+f[22])-1.0*(33.54101966249685*(f[15]+f[11]+f[10]-1.0*f[9])+25.0*f[4]))))/(269.9999999999999*(f[103]+f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*((-1.0*f[82])+f[81]+f[80]-1.0*f[73]+f[72]-1.0*(f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57])+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+27.0*((-1.0*f[27])+f[26]+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*(f[14]+f[13]))+f[12]-1.0*f[8]+f[7]+f[6])-55.90169943749476*(f[5]+f[3]+f[2]-1.0*f[1]))-125.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(269.9999999999999*(f[103]+f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*((-1.0*f[82])+f[81]+f[80]-1.0*f[73]+f[72]-1.0*(f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57])+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+27.0*((-1.0*f[27])+f[26]+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*(f[14]+f[13]))+f[12]-1.0*f[8]+f[7]+f[6])-55.90169943749476*(f[5]+f[3]+f[2]-1.0*f[1]))-125.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.002*(269.9999999999999*(f[103]+f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*((-1.0*f[82])+f[81]+f[80]-1.0*f[73]+f[72]-1.0*(f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57])+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+27.0*((-1.0*f[27])+f[26]+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*(f[14]+f[13]))+f[12]-1.0*f[8]+f[7]+f[6])-55.90169943749476*(f[5]+f[3]+f[2]-1.0*f[1]))-125.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]-1.0*f[107])+60.37383539249431*((-1.0*f[106])+f[105]+f[104]-1.0*f[99]+f[98]-1.0*(f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88])+f[87])))+9.0*(27.0*f[86]+10.0*((-1.0*(f[85]+f[84]))+f[83]-1.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59]))+20.12461179749811*((-1.0*f[55])+f[54]+f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(45.0*((-1.0*(f[30]+f[29]))+f[28]-1.0*f[24]+f[23]+f[22])-1.0*(33.54101966249685*(f[15]+f[11]+f[10]-1.0*f[9])+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(81.0*f[110]+60.37383539249431*((-1.0*f[102])+f[101]+f[100]+f[90])+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]-1.0*f[67]+f[66]+f[65])-1.0*(6.708203932499369*(f[46]+f[42]+f[41]-1.0*f[40])+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.002*(269.9999999999999*(f[103]+f[93]+f[92]-1.0*f[91])+9.0*(22.3606797749979*((-1.0*f[82])+f[81]+f[80]-1.0*f[73]+f[72]-1.0*(f[71]+f[70]-1.0*f[69]+f[68]-1.0*(f[58]+f[57])+f[56]))+45.0*f[52])+11.18033988749895*(13.41640786499874*((-1.0*(f[49]+f[48]))+f[47]-1.0*(f[45]+f[44]+f[43]+f[36]-1.0*f[35]+f[34]+f[33]-1.0*f[32]+f[31]))+27.0*((-1.0*f[27])+f[26]+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*(f[14]+f[13]))+f[12]-1.0*f[8]+f[7]+f[6])-55.90169943749476*(f[5]+f[3]+f[2]-1.0*f[1]))-125.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]-1.0*f[107])+60.37383539249431*((-1.0*f[106])+f[105]+f[104]-1.0*f[99]+f[98]-1.0*(f[97]+f[96]-1.0*f[95]+f[94]-1.0*(f[89]+f[88])+f[87])))+9.0*(27.0*f[86]+10.0*((-1.0*(f[85]+f[84]))+f[83]-1.0*(f[76]+f[75]+f[74]+f[64]-1.0*f[63]+f[62]+f[61]-1.0*f[60]+f[59]))+20.12461179749811*((-1.0*f[55])+f[54]+f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(45.0*((-1.0*(f[30]+f[29]))+f[28]-1.0*f[24]+f[23]+f[22])-1.0*(33.54101966249685*(f[15]+f[11]+f[10]-1.0*f[9])+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (-0.01*(81.0*f[110]+60.37383539249431*((-1.0*f[102])+f[101]+f[100]+f[90])+5.0*(9.0*((-1.0*(f[79]+f[78]))+f[77]-1.0*f[67]+f[66]+f[65])-1.0*(6.708203932499369*(f[46]+f[42]+f[41]-1.0*f[40])+5.0*f[19]))))*fac; 
   } 
  } 
  fReflXYQuad[5][0] = 0.05555555555555555*(fReflXYZMuQuad[7][0]+8.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0]+8.0*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[5][1] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][0]-4188761.0*fReflXYZMuQuad[6][0])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])-4.63256860547201e+14*fReflXYZMuQuad[5][0])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][0]-1.0*(8377522.0*fReflXYZMuQuad[1][0]+2.266096151179001e+23*fReflXYZMuQuad[0][0])))); 
  fReflXYQuad[5][2] = 0.05555555555555555*(fReflXYZMuQuad[7][1]+8.0*fReflXYZMuQuad[6][1]+fReflXYZMuQuad[5][1]+8.0*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[5][3] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][0]+7.4121097687552e+14*fReflXYZMuQuad[6][0]+4.63256860547201e+14*fReflXYZMuQuad[5][0])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][0]+7.4121097687552e+14*fReflXYZMuQuad[1][0]+4.63256860547201e+14*fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[5][4] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][1]-4188761.0*fReflXYZMuQuad[6][1])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])-4.63256860547201e+14*fReflXYZMuQuad[5][1])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][1]-1.0*(8377522.0*fReflXYZMuQuad[1][1]+2.266096151179001e+23*fReflXYZMuQuad[0][1])))); 
  fReflXYQuad[5][5] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][0]+9.0*fReflXYZMuQuad[6][0])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][0]-1.346286087882789e+17*fReflXYZMuQuad[5][0])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][0]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])))); 
  fReflXYQuad[5][6] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][1]+7.4121097687552e+14*fReflXYZMuQuad[6][1]+4.63256860547201e+14*fReflXYZMuQuad[5][1])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][1]+7.4121097687552e+14*fReflXYZMuQuad[1][1]+4.63256860547201e+14*fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[5][7] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][0]-4.17103711299121e+14*fReflXYZMuQuad[6][0])+5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[5][0]+2384663.0*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])+3.190559553141742e+22*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[5][8] = 0.05555555555555555*(fReflXYZMuQuad[7][2]+8.0*fReflXYZMuQuad[6][2]+fReflXYZMuQuad[5][2]+8.0*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[5][9] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][0]+1.6692641e+7*fReflXYZMuQuad[6][0])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][0]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[2][0]+2384663.0*fReflXYZMuQuad[1][0]+3.190559553141742e+22*fReflXYZMuQuad[0][0])-3.20880527843592e+30*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])))); 
  fReflXYQuad[5][10] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][1]+9.0*fReflXYZMuQuad[6][1])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][1]-1.346286087882789e+17*fReflXYZMuQuad[5][1])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][1]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])))); 
  fReflXYQuad[5][11] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][1]-4.17103711299121e+14*fReflXYZMuQuad[6][1])+5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[5][1]+2384663.0*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])+3.190559553141743e+22*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[5][12] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][2]-4188761.0*fReflXYZMuQuad[6][2])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYZMuQuad[4][2]-4.63256860547201e+14*fReflXYZMuQuad[5][2])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][2]-1.0*(8377522.0*fReflXYZMuQuad[1][2]+2.266096151179001e+23*fReflXYZMuQuad[0][2])))); 
  fReflXYQuad[5][13] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][0]-2.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0])-27.0*fReflXYZMuQuad[3][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[2][0])+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[5][14] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][2]+7.4121097687552e+14*fReflXYZMuQuad[6][2]+4.63256860547201e+14*fReflXYZMuQuad[5][2])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][2]+7.4121097687552e+14*fReflXYZMuQuad[1][2]+4.63256860547201e+14*fReflXYZMuQuad[0][2]))); 
  fReflXYQuad[5][15] = 1.719407810605221e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][0]+27.0*fReflXYZMuQuad[6][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[5][0])+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[4][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[5][16] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][1]+1.6692641e+7*fReflXYZMuQuad[6][1])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][1]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[2][1]+2384663.0*fReflXYZMuQuad[1][1]+3.190559553141743e+22*fReflXYZMuQuad[0][1])-3.20880527843592e+30*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])))); 
  fReflXYQuad[5][17] = 1.719407810605222e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][1]+fReflXYZMuQuad[5][1])-27.0*fReflXYZMuQuad[3][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[2][1])+10.0*fReflXYZMuQuad[1][1]-5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[5][18] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][2]+9.0*fReflXYZMuQuad[6][2])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][2]-1.346286087882789e+17*fReflXYZMuQuad[5][2])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][2]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])))); 
  fReflXYQuad[5][19] = 1.719407810605222e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][1]+27.0*fReflXYZMuQuad[6][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[5][1])+2.0*(5.0*fReflXYZMuQuad[3][1]-5.0*fReflXYZMuQuad[4][1])+5.0*fReflXYZMuQuad[2][1])); 
  } 

 
// node (x,y)_7 
  vcutSq_i = (0.05*q_*(zVal*((63.63961030678928*phiWall[16]-63.63961030678928*phi[16]+47.43416490252571*phiWall[9]-47.43416490252571*phi[9])*zVal-36.74234614174767*phiWall[17]+36.74234614174767*phi[17]+21.90890230020666*phiWall[14]-21.90890230020666*phi[14]-27.38612787525831*phiWall[13]+27.38612787525831*phi[13]+32.86335345030997*phiWall[6]-32.86335345030997*phi[6]+24.49489742783179*phiWall[3]-24.49489742783179*phi[3])-21.21320343559643*phiWall[16]+21.21320343559643*phi[16]-21.21320343559643*phiWall[11]+21.21320343559643*phi[11]-15.8113883008419*phiWall[9]+15.8113883008419*phi[9]+12.64911064067352*phiWall[8]-12.64911064067352*phi[8]-15.8113883008419*phiWall[7]+15.8113883008419*phi[7]+18.97366596101028*phiWall[2]-18.97366596101028*phi[2]+14.14213562373095*phiWall[0]-14.14213562373095*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[6][0] = 0.0; 
  fReflXYQuad[6][1] = 0.0; 
  fReflXYQuad[6][2] = 0.0; 
  fReflXYQuad[6][3] = 0.0; 
  fReflXYQuad[6][4] = 0.0; 
  fReflXYQuad[6][5] = 0.0; 
  fReflXYQuad[6][6] = 0.0; 
  fReflXYQuad[6][7] = 0.0; 
  fReflXYQuad[6][8] = 0.0; 
  fReflXYQuad[6][9] = 0.0; 
  fReflXYQuad[6][10] = 0.0; 
  fReflXYQuad[6][11] = 0.0; 
  fReflXYQuad[6][12] = 0.0; 
  fReflXYQuad[6][13] = 0.0; 
  fReflXYQuad[6][14] = 0.0; 
  fReflXYQuad[6][15] = 0.0; 
  fReflXYQuad[6][16] = 0.0; 
  fReflXYQuad[6][17] = 0.0; 
  fReflXYQuad[6][18] = 0.0; 
  fReflXYQuad[6][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[6][0] = -0.05*(2.23606797749979*(6.708203932499369*f[31]-4.0*f[17]+5.0*f[16])-2.0*(6.708203932499369*f[2]+5.0*f[0])); 
  fReflXYQuad[6][1] = -0.01666666666666667*(45.0*f[56]+6.708203932499369*(5.0*f[33]-4.0*f[34])-6.0*(6.708203932499369*f[8]+5.0*f[3])); 
  fReflXYQuad[6][2] = -0.01666666666666667*(45.0*f[59]+6.708203932499369*(5.0*f[37]-4.0*f[38])-6.0*(6.708203932499369*f[10]+5.0*f[4])); 
  fReflXYQuad[6][3] = -0.01666666666666667*(45.0*f[68]+6.708203932499369*(5.0*f[43]-4.0*f[44])-6.0*(6.708203932499369*f[13]+5.0*f[5])); 
  fReflXYQuad[6][4] = -0.05*(2.23606797749979*(6.708203932499369*f[87]-4.0*f[62]+5.0*f[61])-2.0*(6.708203932499369*f[24]+5.0*f[11])); 
  fReflXYQuad[6][5] = -0.05*(2.23606797749979*(6.708203932499369*f[91]-4.0*f[71]+5.0*f[70])-2.0*(6.708203932499369*f[27]+5.0*f[14])); 
  fReflXYQuad[6][6] = -0.05*(2.23606797749979*(6.708203932499369*f[94]-4.0*f[75]+5.0*f[74])-2.0*(6.708203932499369*f[29]+5.0*f[15])); 
  fReflXYQuad[6][7] = 0.1*(6.708203932499369*f[36]+5.0*f[18]); 
  fReflXYQuad[6][8] = 0.1*(6.708203932499369*f[41]+5.0*f[19]); 
  fReflXYQuad[6][9] = 0.1*(6.708203932499369*f[48]+5.0*f[20]); 
  fReflXYQuad[6][10] = -0.01666666666666667*(45.0*f[107]+6.708203932499369*(5.0*f[96]-4.0*f[97])-6.0*(6.708203932499369*f[55]+5.0*f[30])); 
  fReflXYQuad[6][11] = 0.1*(6.708203932499369*f[64]+5.0*f[39]); 
  fReflXYQuad[6][12] = 0.1*(6.708203932499369*f[67]+5.0*f[42]); 
  fReflXYQuad[6][13] = 0.1*(6.708203932499369*f[73]+5.0*f[45]); 
  fReflXYQuad[6][14] = 0.1*(6.708203932499369*f[78]+5.0*f[46]); 
  fReflXYQuad[6][15] = 0.1*(6.708203932499369*f[82]+5.0*f[49]); 
  fReflXYQuad[6][16] = 0.1*(6.708203932499369*f[84]+5.0*f[50]); 
  fReflXYQuad[6][17] = 0.1*(6.708203932499369*f[99]+5.0*f[76]); 
  fReflXYQuad[6][18] = 0.1*(6.708203932499369*f[102]+5.0*f[79]); 
  fReflXYQuad[6][19] = 0.1*(6.708203932499369*f[106]+5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.1924500897298753*(405.0*f[107]+60.37383539249431*(4.0*(f[106]+f[99]-1.0*f[97])+5.0*(f[96]-1.0*(f[94]+f[87])))+9.0*(5.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))-40.24922359499622*f[55])+33.54101966249684*(5.0*f[37]-4.0*(f[50]+f[39]+f[38]))+6.0*(3.0*(15.0*((-1.0*f[30])+f[29]+f[24])+11.18033988749895*(f[15]+f[11]))-1.0*(33.54101966249685*f[10]+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*(f[70]-1.0*(f[68]+f[56])))+6.708203932499369*(4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])-5.0*f[33]+5.0*f[31])-54.0*f[27]+5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17])))+2.0*(3.0*(15.0*((-1.0*f[14])+f[13]+f[8])+11.18033988749895*(f[5]+f[3]))-1.0*(33.54101966249685*f[2]+25.0*f[0]))); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*(f[70]-1.0*(f[68]+f[56])))+6.708203932499369*(4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])-5.0*f[33]+5.0*f[31])-54.0*f[27]+5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17])))+2.0*(3.0*(15.0*((-1.0*f[14])+f[13]+f[8])+11.18033988749895*(f[5]+f[3]))-1.0*(33.54101966249685*f[2]+25.0*f[0]))) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*(f[70]-1.0*(f[68]+f[56])))+6.708203932499369*(4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])-5.0*f[33]+5.0*f[31])-54.0*f[27]+5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17])))+2.0*(3.0*(15.0*((-1.0*f[14])+f[13]+f[8])+11.18033988749895*(f[5]+f[3]))-1.0*(33.54101966249685*f[2]+25.0*f[0]))))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*(f[106]+f[99]-1.0*f[97])+5.0*(f[96]-1.0*(f[94]+f[87])))+9.0*(5.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))-40.24922359499622*f[55])+33.54101966249684*(5.0*f[37]-4.0*(f[50]+f[39]+f[38]))+6.0*(3.0*(15.0*((-1.0*f[30])+f[29]+f[24])+11.18033988749895*(f[15]+f[11]))-1.0*(33.54101966249685*f[10]+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(60.37383539249431*f[102]+5.0*(9.0*(f[79]-1.0*(f[78]+f[67]))+6.708203932499369*(f[41]-1.0*(f[46]+f[42]))+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]+f[73]-1.0*f[71])+5.0*(f[70]-1.0*(f[68]+f[56])))+6.708203932499369*(4.0*(f[49]-1.0*f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[34]-1.0*f[36])-5.0*f[33]+5.0*f[31])-54.0*f[27]+5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17])))+2.0*(3.0*(15.0*((-1.0*f[14])+f[13]+f[8])+11.18033988749895*(f[5]+f[3]))-1.0*(33.54101966249685*f[2]+25.0*f[0]))))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*(f[106]+f[99]-1.0*f[97])+5.0*(f[96]-1.0*(f[94]+f[87])))+9.0*(5.0*(4.0*(f[85]-1.0*f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[62]-1.0*f[64])+5.0*(f[59]-1.0*f[61]))-40.24922359499622*f[55])+33.54101966249684*(5.0*f[37]-4.0*(f[50]+f[39]+f[38]))+6.0*(3.0*(15.0*((-1.0*f[30])+f[29]+f[24])+11.18033988749895*(f[15]+f[11]))-1.0*(33.54101966249685*f[10]+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(60.37383539249431*f[102]+5.0*(9.0*(f[79]-1.0*(f[78]+f[67]))+6.708203932499369*(f[41]-1.0*(f[46]+f[42]))+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[99]+f[94])+9.0*(4.0*f[84]+5.0*f[76]-4.0*f[75]+5.0*(f[74]-1.0*(f[64]+f[59])))+6.708203932499369*(4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[10]-1.0*f[15])-3.0*f[29])+5.0*f[4])))/(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44])+5.0*(6.708203932499369*(f[43]-1.0*(f[36]+f[31]))+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[2]-1.0*f[5])-15.0*f[13])+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44])+5.0*(6.708203932499369*(f[43]-1.0*(f[36]+f[31]))+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[2]-1.0*f[5])-15.0*f[13])+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44])+5.0*(6.708203932499369*(f[43]-1.0*(f[36]+f[31]))+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[2]-1.0*f[5])-15.0*f[13])+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(60.37383539249431*(f[99]+f[94])+9.0*(4.0*f[84]+5.0*f[76]-4.0*f[75]+5.0*(f[74]-1.0*(f[64]+f[59])))+6.708203932499369*(4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[10]-1.0*f[15])-3.0*f[29])+5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(9.0*f[78]+6.708203932499369*f[46]-1.0*(6.708203932499369*f[41]+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*(4.0*f[48]+5.0*f[45]-4.0*f[44])+5.0*(6.708203932499369*(f[43]-1.0*(f[36]+f[31]))+4.0*f[20]-5.0*f[18]+4.0*f[17]-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[2]-1.0*f[5])-15.0*f[13])+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.008333333333333333*(60.37383539249431*(f[99]+f[94])+9.0*(4.0*f[84]+5.0*f[76]-4.0*f[75]+5.0*(f[74]-1.0*(f[64]+f[59])))+6.708203932499369*(4.0*f[50]-5.0*f[39]+4.0*f[38]-5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[10]-1.0*f[15])-3.0*f[29])+5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.05*(9.0*f[78]+6.708203932499369*f[46]-1.0*(6.708203932499369*f[41]+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]+60.37383539249431*(4.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(f[96]+f[94]-1.0*f[87]))+9.0*(5.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62]))-1.0*(25.0*(f[61]+f[59])+40.24922359499622*f[55]))+33.54101966249684*(4.0*(f[50]+f[39]+f[38])-5.0*f[37])+6.0*(3.0*(15.0*(f[24]-1.0*(f[30]+f[29]))+11.18033988749895*((-1.0*f[15])+f[11]+f[10]))+25.0*f[4])))/(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68]))+6.708203932499369*(4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34]))-1.0*(33.54101966249684*(f[33]+f[31])+54.0*f[27]-5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16])))+2.0*(3.0*(15.0*(f[8]-1.0*(f[14]+f[13]))+11.18033988749895*((-1.0*f[5])+f[3]+f[2]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68]))+6.708203932499369*(4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34]))-1.0*(33.54101966249684*(f[33]+f[31])+54.0*f[27]-5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16])))+2.0*(3.0*(15.0*(f[8]-1.0*(f[14]+f[13]))+11.18033988749895*((-1.0*f[5])+f[3]+f[2]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68]))+6.708203932499369*(4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34]))-1.0*(33.54101966249684*(f[33]+f[31])+54.0*f[27]-5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16])))+2.0*(3.0*(15.0*(f[8]-1.0*(f[14]+f[13]))+11.18033988749895*((-1.0*f[5])+f[3]+f[2]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(f[96]+f[94]-1.0*f[87]))+9.0*(5.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62]))-1.0*(25.0*(f[61]+f[59])+40.24922359499622*f[55]))+33.54101966249684*(4.0*(f[50]+f[39]+f[38])-5.0*f[37])+6.0*(3.0*(15.0*(f[24]-1.0*(f[30]+f[29]))+11.18033988749895*((-1.0*f[15])+f[11]+f[10]))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(60.37383539249431*f[102]+5.0*(9.0*(f[79]+f[78]-1.0*f[67])+6.708203932499369*f[46]-1.0*(6.708203932499369*(f[42]+f[41])+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*(f[82]-1.0*(f[73]+f[71]))+5.0*(f[70]+f[68]))+6.708203932499369*(4.0*(f[49]+f[48]-1.0*(f[45]+f[44]))+5.0*f[43]+4.0*(f[36]+f[34]))-1.0*(33.54101966249684*(f[33]+f[31])+54.0*f[27]-5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16])))+2.0*(3.0*(15.0*(f[8]-1.0*(f[14]+f[13]))+11.18033988749895*((-1.0*f[5])+f[3]+f[2]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*(f[106]-1.0*(f[99]+f[97]))+5.0*(f[96]+f[94]-1.0*f[87]))+9.0*(5.0*(4.0*(f[85]+f[84]-1.0*(f[76]+f[75]))+5.0*f[74]+4.0*(f[64]+f[62]))-1.0*(25.0*(f[61]+f[59])+40.24922359499622*f[55]))+33.54101966249684*(4.0*(f[50]+f[39]+f[38])-5.0*f[37])+6.0*(3.0*(15.0*(f[24]-1.0*(f[30]+f[29]))+11.18033988749895*((-1.0*f[15])+f[11]+f[10]))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(60.37383539249431*f[102]+5.0*(9.0*(f[79]+f[78]-1.0*f[67])+6.708203932499369*f[46]-1.0*(6.708203932499369*(f[42]+f[41])+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[106]+f[87])+9.0*(5.0*(f[85]-1.0*f[84])+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+6.708203932499369*((-5.0*f[50])+4.0*(f[39]+f[38])-5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[10]-1.0*f[11])-3.0*f[24])+5.0*f[4])))/(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*(5.0*f[49]+4.0*(f[36]-1.0*f[34]))+5.0*(6.708203932499369*f[33]-1.0*(6.708203932499369*f[31]+5.0*f[20]-4.0*(f[18]+f[17]))-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[2]-1.0*f[3])-15.0*f[8])+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*(5.0*f[49]+4.0*(f[36]-1.0*f[34]))+5.0*(6.708203932499369*f[33]-1.0*(6.708203932499369*f[31]+5.0*f[20]-4.0*(f[18]+f[17]))-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[2]-1.0*f[3])-15.0*f[8])+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*(5.0*f[49]+4.0*(f[36]-1.0*f[34]))+5.0*(6.708203932499369*f[33]-1.0*(6.708203932499369*f[31]+5.0*f[20]-4.0*(f[18]+f[17]))-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[2]-1.0*f[3])-15.0*f[8])+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(60.37383539249431*(f[106]+f[87])+9.0*(5.0*(f[85]-1.0*f[84])+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+6.708203932499369*((-5.0*f[50])+4.0*(f[39]+f[38])-5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[10]-1.0*f[11])-3.0*f[24])+5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(9.0*f[67]+6.708203932499369*f[42]-1.0*(6.708203932499369*f[41]+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*(5.0*f[49]+4.0*(f[36]-1.0*f[34]))+5.0*(6.708203932499369*f[33]-1.0*(6.708203932499369*f[31]+5.0*f[20]-4.0*(f[18]+f[17]))-5.0*f[16]))+2.0*(3.0*(11.18033988749895*(f[2]-1.0*f[3])-15.0*f[8])+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.008333333333333333*(60.37383539249431*(f[106]+f[87])+9.0*(5.0*(f[85]-1.0*f[84])+4.0*(f[64]-1.0*f[62])+5.0*(f[61]-1.0*f[59]))+6.708203932499369*((-5.0*f[50])+4.0*(f[39]+f[38])-5.0*f[37])+6.0*(3.0*(2.23606797749979*(f[10]-1.0*f[11])-3.0*f[24])+5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.05*(9.0*f[67]+6.708203932499369*f[42]-1.0*(6.708203932499369*f[41]+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[106]+f[87])+9.0*(5.0*(f[85]+f[84])-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))+6.708203932499369*(5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])-6.0*(3.0*(3.0*f[24]+2.23606797749979*(f[11]+f[10]))+5.0*f[4])))/(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34]))+5.0*(6.708203932499369*(f[33]+f[31])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]))-2.0*(3.0*(15.0*f[8]+11.18033988749895*(f[3]+f[2]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34]))+5.0*(6.708203932499369*(f[33]+f[31])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]))-2.0*(3.0*(15.0*f[8]+11.18033988749895*(f[3]+f[2]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34]))+5.0*(6.708203932499369*(f[33]+f[31])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]))-2.0*(3.0*(15.0*f[8]+11.18033988749895*(f[3]+f[2]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(60.37383539249431*(f[106]+f[87])+9.0*(5.0*(f[85]+f[84])-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))+6.708203932499369*(5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])-6.0*(3.0*(3.0*f[24]+2.23606797749979*(f[11]+f[10]))+5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(9.0*f[67]+6.708203932499369*(f[42]+f[41])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(2.23606797749979*(45.0*(f[82]+f[56])+6.708203932499369*(5.0*(f[49]+f[48])-4.0*(f[36]+f[34]))+5.0*(6.708203932499369*(f[33]+f[31])+5.0*f[20]-4.0*(f[18]+f[17])+5.0*f[16]))-2.0*(3.0*(15.0*f[8]+11.18033988749895*(f[3]+f[2]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.008333333333333333*(60.37383539249431*(f[106]+f[87])+9.0*(5.0*(f[85]+f[84])-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))+6.708203932499369*(5.0*f[50]-4.0*(f[39]+f[38])+5.0*f[37])-6.0*(3.0*(3.0*f[24]+2.23606797749979*(f[11]+f[10]))+5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.05*(9.0*f[67]+6.708203932499369*(f[42]+f[41])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]+60.37383539249431*(4.0*((-1.0*f[106])+f[99]-1.0*f[97])+5.0*(f[96]-1.0*f[94]+f[87]))+9.0*(5.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*f[61])-1.0*(25.0*f[59]+40.24922359499622*f[55]))+33.54101966249684*(4.0*(f[50]+f[39]+f[38])-5.0*f[37])+6.0*(3.0*(15.0*((-1.0*f[30])+f[29]-1.0*f[24])+11.18033988749895*(f[15]-1.0*f[11]+f[10]))+25.0*f[4])))/(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*((-1.0*f[82])+f[73]-1.0*f[71])+5.0*(f[70]+f[56]))+6.708203932499369*(4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*f[33])-1.0*(33.54101966249684*f[31]+54.0*f[27]-5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16])))+2.0*(3.0*(15.0*((-1.0*f[14])+f[13]-1.0*f[8])+11.18033988749895*(f[5]-1.0*f[3]+f[2]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*((-1.0*f[82])+f[73]-1.0*f[71])+5.0*(f[70]+f[56]))+6.708203932499369*(4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*f[33])-1.0*(33.54101966249684*f[31]+54.0*f[27]-5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16])))+2.0*(3.0*(15.0*((-1.0*f[14])+f[13]-1.0*f[8])+11.18033988749895*(f[5]-1.0*f[3]+f[2]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*((-1.0*f[82])+f[73]-1.0*f[71])+5.0*(f[70]+f[56]))+6.708203932499369*(4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*f[33])-1.0*(33.54101966249684*f[31]+54.0*f[27]-5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16])))+2.0*(3.0*(15.0*((-1.0*f[14])+f[13]-1.0*f[8])+11.18033988749895*(f[5]-1.0*f[3]+f[2]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*((-1.0*f[106])+f[99]-1.0*f[97])+5.0*(f[96]-1.0*f[94]+f[87]))+9.0*(5.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*f[61])-1.0*(25.0*f[59]+40.24922359499622*f[55]))+33.54101966249684*(4.0*(f[50]+f[39]+f[38])-5.0*f[37])+6.0*(3.0*(15.0*((-1.0*f[30])+f[29]-1.0*f[24])+11.18033988749895*(f[15]-1.0*f[11]+f[10]))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(60.37383539249431*f[102]+5.0*(9.0*(f[79]-1.0*f[78]+f[67])+6.708203932499369*(f[42]-1.0*f[46])-1.0*(6.708203932499369*f[41]+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(4.0*((-1.0*f[82])+f[73]-1.0*f[71])+5.0*(f[70]+f[56]))+6.708203932499369*(4.0*((-1.0*f[49])+f[48]+f[45]+f[44])-5.0*f[43]+4.0*(f[36]-1.0*f[34])+5.0*f[33])-1.0*(33.54101966249684*f[31]+54.0*f[27]-5.0*(4.0*(f[20]+f[18]+f[17])-5.0*f[16])))+2.0*(3.0*(15.0*((-1.0*f[14])+f[13]-1.0*f[8])+11.18033988749895*(f[5]-1.0*f[3]+f[2]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (0.001666666666666667*(405.0*f[107]+60.37383539249431*(4.0*((-1.0*f[106])+f[99]-1.0*f[97])+5.0*(f[96]-1.0*f[94]+f[87]))+9.0*(5.0*(4.0*((-1.0*f[85])+f[84]+f[76]+f[75])-5.0*f[74]+4.0*(f[64]-1.0*f[62])+5.0*f[61])-1.0*(25.0*f[59]+40.24922359499622*f[55]))+33.54101966249684*(4.0*(f[50]+f[39]+f[38])-5.0*f[37])+6.0*(3.0*(15.0*((-1.0*f[30])+f[29]-1.0*f[24])+11.18033988749895*(f[15]-1.0*f[11]+f[10]))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(60.37383539249431*f[102]+5.0*(9.0*(f[79]-1.0*f[78]+f[67])+6.708203932499369*(f[42]-1.0*f[46])-1.0*(6.708203932499369*f[41]+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(60.37383539249431*(f[99]+f[94])+9.0*((-4.0*f[84])+5.0*f[76]-4.0*f[75]+5.0*(f[74]+f[64]+f[59]))+6.708203932499369*((-4.0*f[50])+5.0*f[39]-4.0*f[38]+5.0*f[37])-6.0*(3.0*(3.0*f[29]+2.23606797749979*(f[15]+f[10]))+5.0*f[4])))/(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*((-4.0*f[48])+5.0*f[45]-4.0*f[44])+5.0*(6.708203932499369*(f[43]+f[36]+f[31])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]))-2.0*(3.0*(15.0*f[13]+11.18033988749895*(f[5]+f[2]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*((-4.0*f[48])+5.0*f[45]-4.0*f[44])+5.0*(6.708203932499369*(f[43]+f[36]+f[31])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]))-2.0*(3.0*(15.0*f[13]+11.18033988749895*(f[5]+f[2]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*((-4.0*f[48])+5.0*f[45]-4.0*f[44])+5.0*(6.708203932499369*(f[43]+f[36]+f[31])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]))-2.0*(3.0*(15.0*f[13]+11.18033988749895*(f[5]+f[2]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.008333333333333333*(60.37383539249431*(f[99]+f[94])+9.0*((-4.0*f[84])+5.0*f[76]-4.0*f[75]+5.0*(f[74]+f[64]+f[59]))+6.708203932499369*((-4.0*f[50])+5.0*f[39]-4.0*f[38]+5.0*f[37])-6.0*(3.0*(3.0*f[29]+2.23606797749979*(f[15]+f[10]))+5.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.05*(9.0*f[78]+6.708203932499369*(f[46]+f[41])+5.0*f[19]))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(45.0*(f[73]+f[68])+6.708203932499369*((-4.0*f[48])+5.0*f[45]-4.0*f[44])+5.0*(6.708203932499369*(f[43]+f[36]+f[31])-4.0*f[20]+5.0*f[18]-4.0*f[17]+5.0*f[16]))-2.0*(3.0*(15.0*f[13]+11.18033988749895*(f[5]+f[2]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.008333333333333333*(60.37383539249431*(f[99]+f[94])+9.0*((-4.0*f[84])+5.0*f[76]-4.0*f[75]+5.0*(f[74]+f[64]+f[59]))+6.708203932499369*((-4.0*f[50])+5.0*f[39]-4.0*f[38]+5.0*f[37])-6.0*(3.0*(3.0*f[29]+2.23606797749979*(f[15]+f[10]))+5.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.05*(9.0*f[78]+6.708203932499369*(f[46]+f[41])+5.0*f[19]))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[107]+60.37383539249431*(5.0*(f[96]+f[94]+f[87])-4.0*(f[106]+f[99]+f[97]))+9.0*(5.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-40.24922359499622*f[55])+33.54101966249684*(5.0*f[37]-4.0*(f[50]+f[39]+f[38]))-6.0*(3.0*(15.0*(f[30]+f[29]+f[24])+11.18033988749895*(f[15]+f[11]+f[10]))+25.0*f[4])))/(2.23606797749979*(60.37383539249431*f[91]+9.0*(5.0*(f[70]+f[68]+f[56])-4.0*(f[82]+f[73]+f[71]))+6.708203932499369*((-4.0*(f[49]+f[48]+f[45]+f[44]))+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-54.0*f[27]+5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17])))-2.0*(3.0*(15.0*(f[14]+f[13]+f[8])+11.18033988749895*(f[5]+f[3]+f[2]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(5.0*(f[70]+f[68]+f[56])-4.0*(f[82]+f[73]+f[71]))+6.708203932499369*((-4.0*(f[49]+f[48]+f[45]+f[44]))+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-54.0*f[27]+5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17])))-2.0*(3.0*(15.0*(f[14]+f[13]+f[8])+11.18033988749895*(f[5]+f[3]+f[2]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(5.0*(f[70]+f[68]+f[56])-4.0*(f[82]+f[73]+f[71]))+6.708203932499369*((-4.0*(f[49]+f[48]+f[45]+f[44]))+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-54.0*f[27]+5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17])))-2.0*(3.0*(15.0*(f[14]+f[13]+f[8])+11.18033988749895*(f[5]+f[3]+f[2]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.001666666666666667*(405.0*f[107]+60.37383539249431*(5.0*(f[96]+f[94]+f[87])-4.0*(f[106]+f[99]+f[97]))+9.0*(5.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-40.24922359499622*f[55])+33.54101966249684*(5.0*f[37]-4.0*(f[50]+f[39]+f[38]))-6.0*(3.0*(15.0*(f[30]+f[29]+f[24])+11.18033988749895*(f[15]+f[11]+f[10]))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(60.37383539249431*f[102]+5.0*(9.0*(f[79]+f[78]+f[67])+6.708203932499369*(f[46]+f[42]+f[41])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[91]+9.0*(5.0*(f[70]+f[68]+f[56])-4.0*(f[82]+f[73]+f[71]))+6.708203932499369*((-4.0*(f[49]+f[48]+f[45]+f[44]))+5.0*f[43]-4.0*(f[36]+f[34])+5.0*(f[33]+f[31]))-54.0*f[27]+5.0*(5.0*f[16]-4.0*(f[20]+f[18]+f[17])))-2.0*(3.0*(15.0*(f[14]+f[13]+f[8])+11.18033988749895*(f[5]+f[3]+f[2]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (-0.001666666666666667*(405.0*f[107]+60.37383539249431*(5.0*(f[96]+f[94]+f[87])-4.0*(f[106]+f[99]+f[97]))+9.0*(5.0*((-4.0*(f[85]+f[84]+f[76]+f[75]))+5.0*f[74]-4.0*(f[64]+f[62])+5.0*(f[61]+f[59]))-40.24922359499622*f[55])+33.54101966249684*(5.0*f[37]-4.0*(f[50]+f[39]+f[38]))-6.0*(3.0*(15.0*(f[30]+f[29]+f[24])+11.18033988749895*(f[15]+f[11]+f[10]))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(60.37383539249431*f[102]+5.0*(9.0*(f[79]+f[78]+f[67])+6.708203932499369*(f[46]+f[42]+f[41])+5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[6][0] = 0.05555555555555555*(fReflXYZMuQuad[7][0]+8.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0]+8.0*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[6][1] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][0]-4188761.0*fReflXYZMuQuad[6][0])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])-4.63256860547201e+14*fReflXYZMuQuad[5][0])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][0]-1.0*(8377522.0*fReflXYZMuQuad[1][0]+2.266096151179001e+23*fReflXYZMuQuad[0][0])))); 
  fReflXYQuad[6][2] = 0.05555555555555555*(fReflXYZMuQuad[7][1]+8.0*fReflXYZMuQuad[6][1]+fReflXYZMuQuad[5][1]+8.0*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[6][3] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][0]+7.4121097687552e+14*fReflXYZMuQuad[6][0]+4.63256860547201e+14*fReflXYZMuQuad[5][0])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][0]+7.4121097687552e+14*fReflXYZMuQuad[1][0]+4.63256860547201e+14*fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[6][4] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][1]-4188761.0*fReflXYZMuQuad[6][1])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])-4.63256860547201e+14*fReflXYZMuQuad[5][1])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][1]-1.0*(8377522.0*fReflXYZMuQuad[1][1]+2.266096151179001e+23*fReflXYZMuQuad[0][1])))); 
  fReflXYQuad[6][5] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][0]+9.0*fReflXYZMuQuad[6][0])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][0]-1.346286087882789e+17*fReflXYZMuQuad[5][0])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][0]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])))); 
  fReflXYQuad[6][6] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][1]+7.4121097687552e+14*fReflXYZMuQuad[6][1]+4.63256860547201e+14*fReflXYZMuQuad[5][1])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][1]+7.4121097687552e+14*fReflXYZMuQuad[1][1]+4.63256860547201e+14*fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[6][7] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][0]-4.17103711299121e+14*fReflXYZMuQuad[6][0])+5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[5][0]+2384663.0*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])+3.190559553141742e+22*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[6][8] = 0.05555555555555555*(fReflXYZMuQuad[7][2]+8.0*fReflXYZMuQuad[6][2]+fReflXYZMuQuad[5][2]+8.0*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[6][9] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][0]+1.6692641e+7*fReflXYZMuQuad[6][0])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][0]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[2][0]+2384663.0*fReflXYZMuQuad[1][0]+3.190559553141742e+22*fReflXYZMuQuad[0][0])-3.20880527843592e+30*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])))); 
  fReflXYQuad[6][10] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][1]+9.0*fReflXYZMuQuad[6][1])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][1]-1.346286087882789e+17*fReflXYZMuQuad[5][1])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][1]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])))); 
  fReflXYQuad[6][11] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][1]-4.17103711299121e+14*fReflXYZMuQuad[6][1])+5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[5][1]+2384663.0*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])+3.190559553141743e+22*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[6][12] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][2]-4188761.0*fReflXYZMuQuad[6][2])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYZMuQuad[4][2]-4.63256860547201e+14*fReflXYZMuQuad[5][2])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][2]-1.0*(8377522.0*fReflXYZMuQuad[1][2]+2.266096151179001e+23*fReflXYZMuQuad[0][2])))); 
  fReflXYQuad[6][13] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][0]-2.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0])-27.0*fReflXYZMuQuad[3][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[2][0])+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[6][14] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][2]+7.4121097687552e+14*fReflXYZMuQuad[6][2]+4.63256860547201e+14*fReflXYZMuQuad[5][2])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][2]+7.4121097687552e+14*fReflXYZMuQuad[1][2]+4.63256860547201e+14*fReflXYZMuQuad[0][2]))); 
  fReflXYQuad[6][15] = 1.719407810605221e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][0]+27.0*fReflXYZMuQuad[6][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[5][0])+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[4][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[6][16] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][1]+1.6692641e+7*fReflXYZMuQuad[6][1])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][1]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[2][1]+2384663.0*fReflXYZMuQuad[1][1]+3.190559553141743e+22*fReflXYZMuQuad[0][1])-3.20880527843592e+30*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])))); 
  fReflXYQuad[6][17] = 1.719407810605222e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][1]+fReflXYZMuQuad[5][1])-27.0*fReflXYZMuQuad[3][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[2][1])+10.0*fReflXYZMuQuad[1][1]-5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[6][18] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][2]+9.0*fReflXYZMuQuad[6][2])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][2]-1.346286087882789e+17*fReflXYZMuQuad[5][2])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][2]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])))); 
  fReflXYQuad[6][19] = 1.719407810605222e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][1]+27.0*fReflXYZMuQuad[6][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[5][1])+2.0*(5.0*fReflXYZMuQuad[3][1]-5.0*fReflXYZMuQuad[4][1])+5.0*fReflXYZMuQuad[2][1])); 
  } 

 
// node (x,y)_8 
  vcutSq_i = (0.01*q_*(zVal*((426.9074841227313*phiWall[19]-426.9074841227313*phi[19]+318.1980515339465*phiWall[16]-318.1980515339465*phi[16]+318.1980515339465*phiWall[15]-318.1980515339465*phi[15]+237.1708245126285*phiWall[9]-237.1708245126285*phi[9])*zVal+146.9693845669907*phiWall[18]-146.9693845669907*phi[18]+146.9693845669907*phiWall[17]-146.9693845669907*phi[17]+109.5445115010333*phiWall[14]-109.5445115010333*phi[14]+109.5445115010333*phiWall[13]-109.5445115010333*phi[13]+220.454076850486*phiWall[10]-220.454076850486*phi[10]+164.3167672515499*phiWall[6]-164.3167672515499*phi[6]+164.3167672515499*phiWall[5]-164.3167672515499*phi[5]+122.4744871391589*phiWall[3]-122.4744871391589*phi[3])-142.3024947075771*phiWall[19]+142.3024947075771*phi[19]-106.0660171779822*phiWall[16]+106.0660171779822*phi[16]-106.0660171779822*phiWall[15]+106.0660171779822*phi[15]+84.85281374238573*phiWall[12]-84.85281374238573*phi[12]+84.85281374238573*phiWall[11]-84.85281374238573*phi[11]-79.0569415042095*phiWall[9]+79.0569415042095*phi[9]+63.24555320336762*phiWall[8]-63.24555320336762*phi[8]+63.24555320336762*phiWall[7]-63.24555320336762*phi[7]+127.2792206135786*phiWall[4]-127.2792206135786*phi[4]+94.86832980505142*phiWall[2]-94.86832980505142*phi[2]+94.86832980505142*phiWall[1]-94.86832980505142*phi[1]+70.71067811865477*phiWall[0]-70.71067811865477*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[7][0] = 0.0; 
  fReflXYQuad[7][1] = 0.0; 
  fReflXYQuad[7][2] = 0.0; 
  fReflXYQuad[7][3] = 0.0; 
  fReflXYQuad[7][4] = 0.0; 
  fReflXYQuad[7][5] = 0.0; 
  fReflXYQuad[7][6] = 0.0; 
  fReflXYQuad[7][7] = 0.0; 
  fReflXYQuad[7][8] = 0.0; 
  fReflXYQuad[7][9] = 0.0; 
  fReflXYQuad[7][10] = 0.0; 
  fReflXYQuad[7][11] = 0.0; 
  fReflXYQuad[7][12] = 0.0; 
  fReflXYQuad[7][13] = 0.0; 
  fReflXYQuad[7][14] = 0.0; 
  fReflXYQuad[7][15] = 0.0; 
  fReflXYQuad[7][16] = 0.0; 
  fReflXYQuad[7][17] = 0.0; 
  fReflXYQuad[7][18] = 0.0; 
  fReflXYQuad[7][19] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[7][0] = 0.02*(4.47213595499958*(6.708203932499369*(f[32]+f[31])+5.0*(f[17]+f[16]))+3.0*(15.0*f[6]+11.18033988749895*(f[2]+f[1]))+25.0*f[0]); 
  fReflXYQuad[7][1] = 0.03333333333333333*(2.0*(9.0*(f[57]+f[56])+6.708203932499369*(f[34]+f[33]))+3.0*(3.0*(3.0*f[21]+2.23606797749979*(f[8]+f[7]))+5.0*f[3])); 
  fReflXYQuad[7][2] = 0.03333333333333333*(2.0*(9.0*(f[60]+f[59])+6.708203932499369*(f[38]+f[37]))+3.0*(3.0*(3.0*f[22]+2.23606797749979*(f[10]+f[9]))+5.0*f[4])); 
  fReflXYQuad[7][3] = 0.03333333333333333*(2.0*(9.0*(f[69]+f[68])+6.708203932499369*(f[44]+f[43]))+3.0*(3.0*(3.0*f[25]+2.23606797749979*(f[13]+f[12]))+5.0*f[5])); 
  fReflXYQuad[7][4] = 0.02*(4.47213595499958*(6.708203932499369*(f[88]+f[87])+5.0*(f[62]+f[61]))+3.0*(15.0*f[51]+11.18033988749895*(f[24]+f[23]))+25.0*f[11]); 
  fReflXYQuad[7][5] = 0.02*(4.47213595499958*(6.708203932499369*(f[92]+f[91])+5.0*(f[71]+f[70]))+3.0*(15.0*f[52]+11.18033988749895*(f[27]+f[26]))+25.0*f[14]); 
  fReflXYQuad[7][6] = 0.02*(4.47213595499958*(6.708203932499369*(f[95]+f[94])+5.0*(f[75]+f[74]))+3.0*(15.0*f[53]+11.18033988749895*(f[29]+f[28]))+25.0*f[15]); 
  fReflXYQuad[7][7] = 0.1*(9.0*f[58]+6.708203932499369*(f[36]+f[35])+5.0*f[18]); 
  fReflXYQuad[7][8] = 0.1*(9.0*f[65]+6.708203932499369*(f[41]+f[40])+5.0*f[19]); 
  fReflXYQuad[7][9] = 0.1*(9.0*f[80]+6.708203932499369*(f[48]+f[47])+5.0*f[20]); 
  fReflXYQuad[7][10] = 0.03333333333333333*(2.0*(9.0*(f[108]+f[107])+6.708203932499369*(f[97]+f[96]))+3.0*(3.0*(3.0*f[86]+2.23606797749979*(f[55]+f[54]))+5.0*f[30])); 
  fReflXYQuad[7][11] = 0.1*(9.0*f[89]+6.708203932499369*(f[64]+f[63])+5.0*f[39]); 
  fReflXYQuad[7][12] = 0.1*(9.0*f[90]+6.708203932499369*(f[67]+f[66])+5.0*f[42]); 
  fReflXYQuad[7][13] = 0.1*(9.0*f[93]+6.708203932499369*(f[73]+f[72])+5.0*f[45]); 
  fReflXYQuad[7][14] = 0.1*(9.0*f[100]+6.708203932499369*(f[78]+f[77])+5.0*f[46]); 
  fReflXYQuad[7][15] = 0.1*(9.0*f[103]+6.708203932499369*(f[82]+f[81])+5.0*f[49]); 
  fReflXYQuad[7][16] = 0.1*(9.0*f[104]+6.708203932499369*(f[84]+f[83])+5.0*f[50]); 
  fReflXYQuad[7][17] = 0.1*(9.0*f[109]+6.708203932499369*(f[99]+f[98])+5.0*f[76]); 
  fReflXYQuad[7][18] = 0.1*(9.0*f[110]+6.708203932499369*(f[102]+f[101])+5.0*f[79]); 
  fReflXYQuad[7][19] = 0.1*(9.0*f[111]+6.708203932499369*(f[106]+f[105])+5.0*f[85]); 
  } else { // partial reflection 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]-1.0*(f[108]+f[107]))+60.37383539249431*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87]))))+9.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])+f[60]+f[59]))+20.12461179749811*((-1.0*(f[55]+f[54]))+f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[30])+f[29]+f[28]+f[24]+f[23]-1.0*f[22])+11.18033988749895*(f[15]+f[11]))-1.0*(33.54101966249685*(f[10]+f[9])+25.0*f[4]))))/(269.9999999999999*(f[103]+f[93]-1.0*(f[92]+f[91]))+9.0*(22.3606797749979*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))-45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31]))+27.0*((-1.0*(f[27]+f[26]))+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]+f[12]+f[8]+f[7]-1.0*f[6])+55.90169943749476*(f[5]+f[3]))-1.0*(167.7050983124843*(f[2]+f[1])+125.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(269.9999999999999*(f[103]+f[93]-1.0*(f[92]+f[91]))+9.0*(22.3606797749979*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))-45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31]))+27.0*((-1.0*(f[27]+f[26]))+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]+f[12]+f[8]+f[7]-1.0*f[6])+55.90169943749476*(f[5]+f[3]))-1.0*(167.7050983124843*(f[2]+f[1])+125.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[0][0] = 0.0; 
  fReflXYZMuQuad[0][1] = 0.0; 
  fReflXYZMuQuad[0][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.002*(269.9999999999999*(f[103]+f[93]-1.0*(f[92]+f[91]))+9.0*(22.3606797749979*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))-45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31]))+27.0*((-1.0*(f[27]+f[26]))+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]+f[12]+f[8]+f[7]-1.0*f[6])+55.90169943749476*(f[5]+f[3]))-1.0*(167.7050983124843*(f[2]+f[1])+125.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*(f[108]+f[107]))+60.37383539249431*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87]))))+9.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])+f[60]+f[59]))+20.12461179749811*((-1.0*(f[55]+f[54]))+f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[30])+f[29]+f[28]+f[24]+f[23]-1.0*f[22])+11.18033988749895*(f[15]+f[11]))-1.0*(33.54101966249685*(f[10]+f[9])+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(81.0*f[110]+60.37383539249431*(f[102]+f[101]-1.0*(f[100]+f[90]))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]+f[67]+f[66]-1.0*f[65]))+6.708203932499369*((-1.0*(f[46]+f[42]))+f[41]+f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][0] = (-0.002*(269.9999999999999*(f[103]+f[93]-1.0*(f[92]+f[91]))+9.0*(22.3606797749979*(f[82]+f[81]-1.0*f[80]+f[73]+f[72]-1.0*(f[71]+f[70]-1.0*(f[69]+f[68])+f[58]-1.0*(f[57]+f[56])))-45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]-1.0*f[45])+f[44]+f[43]-1.0*(f[36]+f[35]-1.0*(f[34]+f[33])+f[32]+f[31]))+27.0*((-1.0*(f[27]+f[26]))+f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*f[14])+f[13]+f[12]+f[8]+f[7]-1.0*f[6])+55.90169943749476*(f[5]+f[3]))-1.0*(167.7050983124843*(f[2]+f[1])+125.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]+f[109]-1.0*(f[108]+f[107]))+60.37383539249431*(f[106]+f[105]-1.0*f[104]+f[99]+f[98]-1.0*(f[97]+f[96]-1.0*(f[95]+f[94])+f[89]-1.0*(f[88]+f[87]))))+9.0*((-27.0*f[86])+10.0*(f[85]-1.0*(f[84]+f[83]-1.0*f[76])+f[75]+f[74]-1.0*(f[64]+f[63]-1.0*(f[62]+f[61])+f[60]+f[59]))+20.12461179749811*((-1.0*(f[55]+f[54]))+f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*f[30])+f[29]+f[28]+f[24]+f[23]-1.0*f[22])+11.18033988749895*(f[15]+f[11]))-1.0*(33.54101966249685*(f[10]+f[9])+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[0][2] = (0.01*(81.0*f[110]+60.37383539249431*(f[102]+f[101]-1.0*(f[100]+f[90]))+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]+f[67]+f[66]-1.0*f[65]))+6.708203932499369*((-1.0*(f[46]+f[42]))+f[41]+f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+60.37383539249431*(4.0*f[104]+5.0*(f[99]+f[98])-1.0*(4.0*(f[95]+f[94])+5.0*f[89]))+9.0*(5.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(10.0*(f[60]+f[59])-20.12461179749811*f[53]))+33.54101966249684*(4.0*f[50]-5.0*f[39])+2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*(f[22]-1.0*(f[29]+f[28]))+11.18033988749895*((-1.0*f[15])+f[10]+f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]))+6.708203932499369*(4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35])))+2.0*(13.41640786499874*(f[32]+f[31])-27.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[6]-1.0*(f[13]+f[12]))+11.18033988749895*((-1.0*f[5])+f[2]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]))+6.708203932499369*(4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35])))+2.0*(13.41640786499874*(f[32]+f[31])-27.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[6]-1.0*(f[13]+f[12]))+11.18033988749895*((-1.0*f[5])+f[2]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[1][0] = 0.0; 
  fReflXYZMuQuad[1][1] = 0.0; 
  fReflXYZMuQuad[1][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]))+6.708203932499369*(4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35])))+2.0*(13.41640786499874*(f[32]+f[31])-27.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[6]-1.0*(f[13]+f[12]))+11.18033988749895*((-1.0*f[5])+f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.001666666666666667*(405.0*f[109]+60.37383539249431*(4.0*f[104]+5.0*(f[99]+f[98])-1.0*(4.0*(f[95]+f[94])+5.0*f[89]))+9.0*(5.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(10.0*(f[60]+f[59])-20.12461179749811*f[53]))+33.54101966249684*(4.0*f[50]-5.0*f[39])+2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*(f[22]-1.0*(f[29]+f[28]))+11.18033988749895*((-1.0*f[15])+f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01*(60.37383539249431*f[100]+5.0*(9.0*(f[78]+f[77]-1.0*f[65])+6.708203932499369*f[46]-1.0*(6.708203932499369*(f[41]+f[40])+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][0] = (0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*(4.0*f[80]+5.0*(f[73]+f[72])-1.0*(4.0*(f[69]+f[68])+5.0*f[58]))+6.708203932499369*(4.0*(f[48]+f[47])+5.0*f[45]-1.0*(4.0*(f[44]+f[43])+5.0*(f[36]+f[35])))+2.0*(13.41640786499874*(f[32]+f[31])-27.0*f[25])+5.0*(4.0*f[20]-5.0*f[18]))+2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[6]-1.0*(f[13]+f[12]))+11.18033988749895*((-1.0*f[5])+f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][1] = (0.001666666666666667*(405.0*f[109]+60.37383539249431*(4.0*f[104]+5.0*(f[99]+f[98])-1.0*(4.0*(f[95]+f[94])+5.0*f[89]))+9.0*(5.0*(4.0*(f[84]+f[83])+5.0*f[76]-1.0*(4.0*(f[75]+f[74])+5.0*(f[64]+f[63])))+2.0*(10.0*(f[60]+f[59])-20.12461179749811*f[53]))+33.54101966249684*(4.0*f[50]-5.0*f[39])+2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*(f[22]-1.0*(f[29]+f[28]))+11.18033988749895*((-1.0*f[15])+f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[1][2] = (-0.01*(60.37383539249431*f[100]+5.0*(9.0*(f[78]+f[77]-1.0*f[65])+6.708203932499369*f[46]-1.0*(6.708203932499369*(f[41]+f[40])+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]+f[107]))+60.37383539249431*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87]))))+9.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-20.12461179749811*(f[55]+f[54]+f[53]-1.0*f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*(f[30]+f[29]+f[28]-1.0*f[24]))+f[23]+f[22])+11.18033988749895*((-1.0*f[15])+f[11]+f[10]+f[9]))+25.0*f[4])))/(269.9999999999999*(f[103]-1.0*(f[93]+f[92]+f[91]))+9.0*(22.3606797749979*(f[82]+f[81]+f[80])-1.0*(22.3606797749979*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+45.0*f[52]))+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31])-27.0*(f[27]+f[26]+f[25]-1.0*f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*(f[14]+f[13]+f[12]-1.0*f[8]))+f[7]+f[6])+55.90169943749476*((-1.0*f[5])+f[3]+f[2]+f[1]))+125.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(269.9999999999999*(f[103]-1.0*(f[93]+f[92]+f[91]))+9.0*(22.3606797749979*(f[82]+f[81]+f[80])-1.0*(22.3606797749979*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+45.0*f[52]))+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31])-27.0*(f[27]+f[26]+f[25]-1.0*f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*(f[14]+f[13]+f[12]-1.0*f[8]))+f[7]+f[6])+55.90169943749476*((-1.0*f[5])+f[3]+f[2]+f[1]))+125.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[2][0] = 0.0; 
  fReflXYZMuQuad[2][1] = 0.0; 
  fReflXYZMuQuad[2][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.002*(269.9999999999999*(f[103]-1.0*(f[93]+f[92]+f[91]))+9.0*(22.3606797749979*(f[82]+f[81]+f[80])-1.0*(22.3606797749979*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+45.0*f[52]))+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31])-27.0*(f[27]+f[26]+f[25]-1.0*f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*(f[14]+f[13]+f[12]-1.0*f[8]))+f[7]+f[6])+55.90169943749476*((-1.0*f[5])+f[3]+f[2]+f[1]))+125.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]+f[107]))+60.37383539249431*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87]))))+9.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-20.12461179749811*(f[55]+f[54]+f[53]-1.0*f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*(f[30]+f[29]+f[28]-1.0*f[24]))+f[23]+f[22])+11.18033988749895*((-1.0*f[15])+f[11]+f[10]+f[9]))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(81.0*f[110]+60.37383539249431*(f[102]+f[101]+f[100]-1.0*f[90])+5.0*(9.0*(f[79]+f[78]+f[77]-1.0*(f[67]+f[66]+f[65]))+6.708203932499369*f[46]-1.0*(6.708203932499369*(f[42]+f[41]+f[40])+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][0] = (0.002*(269.9999999999999*(f[103]-1.0*(f[93]+f[92]+f[91]))+9.0*(22.3606797749979*(f[82]+f[81]+f[80])-1.0*(22.3606797749979*(f[73]+f[72]+f[71]+f[70]+f[69]+f[68]-1.0*(f[58]+f[57]+f[56]))+45.0*f[52]))+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]-1.0*(f[45]+f[44]+f[43]-1.0*(f[36]+f[35]+f[34]+f[33]))+f[32]+f[31])-27.0*(f[27]+f[26]+f[25]-1.0*f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*((-1.0*(f[14]+f[13]+f[12]-1.0*f[8]))+f[7]+f[6])+55.90169943749476*((-1.0*f[5])+f[3]+f[2]+f[1]))+125.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][1] = (0.003333333333333334*(2.0*(81.0*(f[111]-1.0*(f[109]+f[108]+f[107]))+60.37383539249431*(f[106]+f[105]+f[104]-1.0*(f[99]+f[98]+f[97]+f[96]+f[95]+f[94]-1.0*(f[89]+f[88]+f[87]))))+9.0*((-27.0*f[86])+10.0*(f[85]+f[84]+f[83]-1.0*(f[76]+f[75]+f[74]-1.0*(f[64]+f[63]+f[62]+f[61]))+f[60]+f[59])-20.12461179749811*(f[55]+f[54]+f[53]-1.0*f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*((-1.0*(f[30]+f[29]+f[28]-1.0*f[24]))+f[23]+f[22])+11.18033988749895*((-1.0*f[15])+f[11]+f[10]+f[9]))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[2][2] = (-0.01*(81.0*f[110]+60.37383539249431*(f[102]+f[101]+f[100]-1.0*f[90])+5.0*(9.0*(f[79]+f[78]+f[77]-1.0*(f[67]+f[66]+f[65]))+6.708203932499369*f[46]-1.0*(6.708203932499369*(f[42]+f[41]+f[40])+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+60.37383539249431*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+45.0*(5.0*(f[85]-1.0*(f[84]+f[83]))+4.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))-1.0*(362.243012354966*f[51]+167.7050983124842*f[50])+2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[22]-1.0*(f[24]+f[23]))+11.18033988749895*((-1.0*f[11])+f[10]+f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]+f[81])+4.0*(f[58]-1.0*(f[57]+f[56])))+6.708203932499369*(5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))-1.0*(54.0*f[21]+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[6]-1.0*(f[8]+f[7]))+11.18033988749895*((-1.0*f[3])+f[2]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]+f[81])+4.0*(f[58]-1.0*(f[57]+f[56])))+6.708203932499369*(5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))-1.0*(54.0*f[21]+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[6]-1.0*(f[8]+f[7]))+11.18033988749895*((-1.0*f[3])+f[2]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[3][0] = 0.0; 
  fReflXYZMuQuad[3][1] = 0.0; 
  fReflXYZMuQuad[3][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]+f[81])+4.0*(f[58]-1.0*(f[57]+f[56])))+6.708203932499369*(5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))-1.0*(54.0*f[21]+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[6]-1.0*(f[8]+f[7]))+11.18033988749895*((-1.0*f[3])+f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+45.0*(5.0*(f[85]-1.0*(f[84]+f[83]))+4.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))-1.0*(362.243012354966*f[51]+167.7050983124842*f[50])+2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[22]-1.0*(f[24]+f[23]))+11.18033988749895*((-1.0*f[11])+f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01*(60.37383539249431*f[90]+5.0*(9.0*(f[67]+f[66]-1.0*f[65])+6.708203932499369*f[42]-1.0*(6.708203932499369*(f[41]+f[40])+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][0] = (0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]+f[81])+4.0*(f[58]-1.0*(f[57]+f[56])))+6.708203932499369*(5.0*(f[49]-1.0*(f[48]+f[47]))+4.0*(f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))-1.0*(54.0*f[21]+25.0*f[20]))+2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[6]-1.0*(f[8]+f[7]))+11.18033988749895*((-1.0*f[3])+f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][1] = (0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*(f[106]+f[105]-1.0*f[104])+4.0*(f[89]-1.0*(f[88]+f[87])))+45.0*(5.0*(f[85]-1.0*(f[84]+f[83]))+4.0*(f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))-1.0*(362.243012354966*f[51]+167.7050983124842*f[50])+2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[22]-1.0*(f[24]+f[23]))+11.18033988749895*((-1.0*f[11])+f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[3][2] = (-0.01*(60.37383539249431*f[90]+5.0*(9.0*(f[67]+f[66]-1.0*f[65])+6.708203932499369*f[42]-1.0*(6.708203932499369*(f[41]+f[40])+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[111]+60.37383539249431*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+9.0*(25.0*(f[85]+f[84]+f[83])-2.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+20.12461179749811*f[51]))+167.7050983124842*f[50]-2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[24]+f[23]+f[22])+11.18033988749895*(f[11]+f[10]+f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]+f[81]+f[80])-4.0*(f[58]+f[57]+f[56]))+33.54101966249684*(f[49]+f[48]+f[47])-2.0*(13.41640786499874*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*f[21])+25.0*f[20])-2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[8]+f[7]+f[6])+11.18033988749895*(f[3]+f[2]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]+f[81]+f[80])-4.0*(f[58]+f[57]+f[56]))+33.54101966249684*(f[49]+f[48]+f[47])-2.0*(13.41640786499874*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*f[21])+25.0*f[20])-2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[8]+f[7]+f[6])+11.18033988749895*(f[3]+f[2]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[4][0] = 0.0; 
  fReflXYZMuQuad[4][1] = 0.0; 
  fReflXYZMuQuad[4][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]+f[81]+f[80])-4.0*(f[58]+f[57]+f[56]))+33.54101966249684*(f[49]+f[48]+f[47])-2.0*(13.41640786499874*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*f[21])+25.0*f[20])-2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[8]+f[7]+f[6])+11.18033988749895*(f[3]+f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+9.0*(25.0*(f[85]+f[84]+f[83])-2.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+20.12461179749811*f[51]))+167.7050983124842*f[50]-2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[24]+f[23]+f[22])+11.18033988749895*(f[11]+f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.01*(60.37383539249431*f[90]+5.0*(9.0*(f[67]+f[66]+f[65])+6.708203932499369*(f[42]+f[41]+f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[103]+9.0*(5.0*(f[82]+f[81]+f[80])-4.0*(f[58]+f[57]+f[56]))+33.54101966249684*(f[49]+f[48]+f[47])-2.0*(13.41640786499874*(f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*f[21])+25.0*f[20])-2.0*(22.3606797749979*(f[18]+f[17]+f[16])+3.0*(15.0*(f[8]+f[7]+f[6])+11.18033988749895*(f[3]+f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][1] = (-0.001666666666666667*(405.0*f[111]+60.37383539249431*(5.0*(f[106]+f[105]+f[104])-4.0*(f[89]+f[88]+f[87]))+9.0*(25.0*(f[85]+f[84]+f[83])-2.0*(10.0*(f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+20.12461179749811*f[51]))+167.7050983124842*f[50]-2.0*(67.08203932499369*(f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[24]+f[23]+f[22])+11.18033988749895*(f[11]+f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[4][2] = (0.01*(60.37383539249431*f[90]+5.0*(9.0*(f[67]+f[66]+f[65])+6.708203932499369*(f[42]+f[41]+f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])+60.37383539249431*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87]))))+9.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+20.12461179749811*(f[55]+f[54]-1.0*f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]-1.0*(f[29]+f[28]-1.0*f[24])+f[23]-1.0*f[22])+11.18033988749895*(f[11]-1.0*f[15]))-1.0*(33.54101966249685*(f[10]+f[9])+25.0*f[4]))))/(269.9999999999999*(f[103]-1.0*f[93]+f[92]+f[91])+9.0*(22.3606797749979*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+27.0*(f[27]+f[26]-1.0*f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*(f[13]+f[12]-1.0*f[8])+f[7]-1.0*f[6])+55.90169943749476*(f[3]-1.0*f[5]))-1.0*(167.7050983124843*(f[2]+f[1])+125.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.002*(269.9999999999999*(f[103]-1.0*f[93]+f[92]+f[91])+9.0*(22.3606797749979*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+27.0*(f[27]+f[26]-1.0*f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*(f[13]+f[12]-1.0*f[8])+f[7]-1.0*f[6])+55.90169943749476*(f[3]-1.0*f[5]))-1.0*(167.7050983124843*(f[2]+f[1])+125.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[5][0] = 0.0; 
  fReflXYZMuQuad[5][1] = 0.0; 
  fReflXYZMuQuad[5][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.002*(269.9999999999999*(f[103]-1.0*f[93]+f[92]+f[91])+9.0*(22.3606797749979*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+27.0*(f[27]+f[26]-1.0*f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*(f[13]+f[12]-1.0*f[8])+f[7]-1.0*f[6])+55.90169943749476*(f[3]-1.0*f[5]))-1.0*(167.7050983124843*(f[2]+f[1])+125.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])+60.37383539249431*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87]))))+9.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+20.12461179749811*(f[55]+f[54]-1.0*f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]-1.0*(f[29]+f[28]-1.0*f[24])+f[23]-1.0*f[22])+11.18033988749895*(f[11]-1.0*f[15]))-1.0*(33.54101966249685*(f[10]+f[9])+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(81.0*f[110]+60.37383539249431*(f[102]+f[101]-1.0*f[100]+f[90])+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]-1.0*f[67])+f[66]-1.0*f[65])+6.708203932499369*(f[42]-1.0*f[46])-1.0*(6.708203932499369*(f[41]+f[40])+5.0*f[19]))))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][0] = (-0.002*(269.9999999999999*(f[103]-1.0*f[93]+f[92]+f[91])+9.0*(22.3606797749979*(f[82]+f[81]-1.0*(f[80]+f[73]+f[72]-1.0*(f[71]+f[70])+f[69]+f[68]+f[58]-1.0*(f[57]+f[56])))+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]-1.0*(f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]-1.0*(f[34]+f[33]-1.0*(f[32]+f[31]))))+27.0*(f[27]+f[26]-1.0*f[25]+f[21])-10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]-1.0*(f[13]+f[12]-1.0*f[8])+f[7]-1.0*f[6])+55.90169943749476*(f[3]-1.0*f[5]))-1.0*(167.7050983124843*(f[2]+f[1])+125.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][1] = (-0.003333333333333334*(2.0*(81.0*(f[111]-1.0*f[109]+f[108]+f[107])+60.37383539249431*(f[106]+f[105]-1.0*(f[104]+f[99]+f[98]-1.0*(f[97]+f[96])+f[95]+f[94]+f[89]-1.0*(f[88]+f[87]))))+9.0*(27.0*f[86]+10.0*(f[85]-1.0*(f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]-1.0*(f[62]+f[61]-1.0*(f[60]+f[59]))))+20.12461179749811*(f[55]+f[54]-1.0*f[53]+f[51]))-67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]-1.0*(f[29]+f[28]-1.0*f[24])+f[23]-1.0*f[22])+11.18033988749895*(f[11]-1.0*f[15]))-1.0*(33.54101966249685*(f[10]+f[9])+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[5][2] = (-0.01*(81.0*f[110]+60.37383539249431*(f[102]+f[101]-1.0*f[100]+f[90])+5.0*(9.0*(f[79]-1.0*(f[78]+f[77]-1.0*f[67])+f[66]-1.0*f[65])+6.708203932499369*(f[42]-1.0*f[46])-1.0*(6.708203932499369*(f[41]+f[40])+5.0*f[19]))))*fac; 
   } 
  } 
  xbarVal = (0.1924500897298753*(405.0*f[109]+60.37383539249431*((-4.0*f[104])+5.0*(f[99]+f[98])-4.0*(f[95]+f[94])+5.0*f[89])+9.0*(5.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))-2.0*(10.0*(f[60]+f[59])+20.12461179749811*f[53]))+33.54101966249684*(5.0*f[39]-4.0*f[50])-2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*(f[29]+f[28]+f[22])+11.18033988749895*(f[15]+f[10]+f[9]))+25.0*f[4]))))/(2.23606797749979*(60.37383539249431*f[93]+9.0*((-4.0*f[80])+5.0*(f[73]+f[72])-4.0*(f[69]+f[68])+5.0*f[58])+6.708203932499369*((-4.0*(f[48]+f[47]))+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))-2.0*(13.41640786499874*(f[32]+f[31])+27.0*f[25])+5.0*(5.0*f[18]-4.0*f[20]))-2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[13]+f[12]+f[6])+11.18033988749895*(f[5]+f[2]+f[1]))+25.0*f[0])); 
  // if f is not realizable, no reflection from this node 
  if(-0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*((-4.0*f[80])+5.0*(f[73]+f[72])-4.0*(f[69]+f[68])+5.0*f[58])+6.708203932499369*((-4.0*(f[48]+f[47]))+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))-2.0*(13.41640786499874*(f[32]+f[31])+27.0*f[25])+5.0*(5.0*f[18]-4.0*f[20]))-2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[13]+f[12]+f[6])+11.18033988749895*(f[5]+f[2]+f[1]))+25.0*f[0])) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[6][0] = 0.0; 
  fReflXYZMuQuad[6][1] = 0.0; 
  fReflXYZMuQuad[6][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*((-4.0*f[80])+5.0*(f[73]+f[72])-4.0*(f[69]+f[68])+5.0*f[58])+6.708203932499369*((-4.0*(f[48]+f[47]))+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))-2.0*(13.41640786499874*(f[32]+f[31])+27.0*f[25])+5.0*(5.0*f[18]-4.0*f[20]))-2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[13]+f[12]+f[6])+11.18033988749895*(f[5]+f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.001666666666666667*(405.0*f[109]+60.37383539249431*((-4.0*f[104])+5.0*(f[99]+f[98])-4.0*(f[95]+f[94])+5.0*f[89])+9.0*(5.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))-2.0*(10.0*(f[60]+f[59])+20.12461179749811*f[53]))+33.54101966249684*(5.0*f[39]-4.0*f[50])-2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*(f[29]+f[28]+f[22])+11.18033988749895*(f[15]+f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(60.37383539249431*f[100]+5.0*(9.0*(f[78]+f[77]+f[65])+6.708203932499369*(f[46]+f[41]+f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][0] = (-0.005*(2.23606797749979*(60.37383539249431*f[93]+9.0*((-4.0*f[80])+5.0*(f[73]+f[72])-4.0*(f[69]+f[68])+5.0*f[58])+6.708203932499369*((-4.0*(f[48]+f[47]))+5.0*f[45]-4.0*(f[44]+f[43])+5.0*(f[36]+f[35]))-2.0*(13.41640786499874*(f[32]+f[31])+27.0*f[25])+5.0*(5.0*f[18]-4.0*f[20]))-2.0*(22.3606797749979*(f[17]+f[16])+3.0*(15.0*(f[13]+f[12]+f[6])+11.18033988749895*(f[5]+f[2]+f[1]))+25.0*f[0])))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][1] = (-0.001666666666666667*(405.0*f[109]+60.37383539249431*((-4.0*f[104])+5.0*(f[99]+f[98])-4.0*(f[95]+f[94])+5.0*f[89])+9.0*(5.0*((-4.0*(f[84]+f[83]))+5.0*f[76]-4.0*(f[75]+f[74])+5.0*(f[64]+f[63]))-2.0*(10.0*(f[60]+f[59])+20.12461179749811*f[53]))+33.54101966249684*(5.0*f[39]-4.0*f[50])-2.0*(67.08203932499369*(f[38]+f[37])+3.0*(3.0*(15.0*(f[29]+f[28]+f[22])+11.18033988749895*(f[15]+f[10]+f[9]))+25.0*f[4]))))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[6][2] = (0.01*(60.37383539249431*f[100]+5.0*(9.0*(f[78]+f[77]+f[65])+6.708203932499369*(f[46]+f[41]+f[40])+5.0*f[19])))*fac; 
   } 
  } 
  xbarVal = (0.9622504486493765*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])+60.37383539249431*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87]))+9.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+20.12461179749811*(f[55]+f[54]+f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+11.18033988749895*(f[15]+f[11]+f[10]+f[9]))+25.0*f[4])))/(269.9999999999999*(f[103]+f[93]+f[92]+f[91])+9.0*(22.3606797749979*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*(f[27]+f[26]+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+55.90169943749476*(f[5]+f[3]+f[2]+f[1]))+125.0*f[0]); 
  // if f is not realizable, no reflection from this node 
  if(0.002*(269.9999999999999*(f[103]+f[93]+f[92]+f[91])+9.0*(22.3606797749979*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*(f[27]+f[26]+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+55.90169943749476*(f[5]+f[3]+f[2]+f[1]))+125.0*f[0]) <= 0. || std::abs(xbarVal)>=.95) { 
  fReflXYZMuQuad[7][0] = 0.0; 
  fReflXYZMuQuad[7][1] = 0.0; 
  fReflXYZMuQuad[7][2] = 0.0; 
  } else {
   b = invL(xbarVal); 
   if(wv > 0) {
    xc = 2.*(std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 0. : b<-500? 1. : std::abs(b)<1e-10? (1.+xc)/2. : (std::exp(b*xc)-std::exp(-b))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.002*(269.9999999999999*(f[103]+f[93]+f[92]+f[91])+9.0*(22.3606797749979*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*(f[27]+f[26]+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+55.90169943749476*(f[5]+f[3]+f[2]+f[1]))+125.0*f[0]))*fac; 
    fac = (b>500 || std::abs(b)<1e-8)? 0. : b<-500? 1. : ((b*xc-1)*std::exp(b*xc)+(b+1)*std::exp(-b))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])+60.37383539249431*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87]))+9.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+20.12461179749811*(f[55]+f[54]+f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+11.18033988749895*(f[15]+f[11]+f[10]+f[9]))+25.0*f[4])))*fac; 
    fac = (((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3-(2*(b*b+3*(b+1))*std::exp(-b))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(81.0*f[110]+60.37383539249431*(f[102]+f[101]+f[100]+f[90])+5.0*(9.0*(f[79]+f[78]+f[77]+f[67]+f[66]+f[65])+6.708203932499369*(f[46]+f[42]+f[41]+f[40])+5.0*f[19])))*fac; 
   } else { 
    xc = 2.*(-std::sqrt(vcutSq_i)-wv)/dv; 
    fac = b>500? 1. : b<-500? 0. : std::abs(b)<1e-10? (1.-xc)/2. : (std::exp(b)-std::exp(b*xc))/(2.*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][0] = (0.002*(269.9999999999999*(f[103]+f[93]+f[92]+f[91])+9.0*(22.3606797749979*(f[82]+f[81]+f[80]+f[73]+f[72]+f[71]+f[70]+f[69]+f[68]+f[58]+f[57]+f[56])+45.0*f[52])+11.18033988749895*(13.41640786499874*(f[49]+f[48]+f[47]+f[45]+f[44]+f[43]+f[36]+f[35]+f[34]+f[33]+f[32]+f[31])+27.0*(f[27]+f[26]+f[25]+f[21])+10.0*(f[20]+f[18]+f[17]+f[16]))+3.0*(75.0*(f[14]+f[13]+f[12]+f[8]+f[7]+f[6])+55.90169943749476*(f[5]+f[3]+f[2]+f[1]))+125.0*f[0]))*fac; 
    fac = b>500? 1. : (b<-500 || std::abs(b)<1e-8)? 0. : ((b-1)*std::exp(b)-(b*xc-1)*std::exp(b*xc))/2./(b*std::cosh(b)-std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][1] = (0.003333333333333334*(2.0*(81.0*(f[111]+f[109]+f[108]+f[107])+60.37383539249431*(f[106]+f[105]+f[104]+f[99]+f[98]+f[97]+f[96]+f[95]+f[94]+f[89]+f[88]+f[87]))+9.0*(27.0*f[86]+10.0*(f[85]+f[84]+f[83]+f[76]+f[75]+f[74]+f[64]+f[63]+f[62]+f[61]+f[60]+f[59])+20.12461179749811*(f[55]+f[54]+f[53]+f[51]))+67.08203932499369*(f[50]+f[39]+f[38]+f[37])+3.0*(3.0*(15.0*(f[30]+f[29]+f[28]+f[24]+f[23]+f[22])+11.18033988749895*(f[15]+f[11]+f[10]+f[9]))+25.0*f[4])))*fac; 
    fac = ((2*(b*b+3*(1-b))*std::exp(b))/3-((b*(3*b*xc*xc-(6*xc+b))+6)*std::exp(b*xc))/3)/(-4*b*std::cosh(b) + 4/3*(3 + b*b)*std::sinh(b)); 
    if(std::isnan(fac) || std::isinf(fac)) {printf("reflect fac = %LG, b=%LG, xbarVal=%LG \n", fac, b, xbarVal); fac=0.;} 
    fReflXYZMuQuad[7][2] = (0.01*(81.0*f[110]+60.37383539249431*(f[102]+f[101]+f[100]+f[90])+5.0*(9.0*(f[79]+f[78]+f[77]+f[67]+f[66]+f[65])+6.708203932499369*(f[46]+f[42]+f[41]+f[40])+5.0*f[19])))*fac; 
   } 
  } 
  fReflXYQuad[7][0] = 0.05555555555555555*(fReflXYZMuQuad[7][0]+8.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0]+8.0*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+fReflXYZMuQuad[2][0]+8.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]); 
  fReflXYQuad[7][1] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][0]-4188761.0*fReflXYZMuQuad[6][0])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])-4.63256860547201e+14*fReflXYZMuQuad[5][0])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][0]-1.0*(8377522.0*fReflXYZMuQuad[1][0]+2.266096151179001e+23*fReflXYZMuQuad[0][0])))); 
  fReflXYQuad[7][2] = 0.05555555555555555*(fReflXYZMuQuad[7][1]+8.0*fReflXYZMuQuad[6][1]+fReflXYZMuQuad[5][1]+8.0*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+fReflXYZMuQuad[2][1]+8.0*fReflXYZMuQuad[1][1]+fReflXYZMuQuad[0][1]); 
  fReflXYQuad[7][3] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][0]+7.4121097687552e+14*fReflXYZMuQuad[6][0]+4.63256860547201e+14*fReflXYZMuQuad[5][0])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][0]+7.4121097687552e+14*fReflXYZMuQuad[1][0]+4.63256860547201e+14*fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[7][4] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][1]-4188761.0*fReflXYZMuQuad[6][1])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])-4.63256860547201e+14*fReflXYZMuQuad[5][1])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][1]-1.0*(8377522.0*fReflXYZMuQuad[1][1]+2.266096151179001e+23*fReflXYZMuQuad[0][1])))); 
  fReflXYQuad[7][5] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][0]+9.0*fReflXYZMuQuad[6][0])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][0]-1.346286087882789e+17*fReflXYZMuQuad[5][0])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][0]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][0]+fReflXYZMuQuad[1][0])))); 
  fReflXYQuad[7][6] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][1]+7.4121097687552e+14*fReflXYZMuQuad[6][1]+4.63256860547201e+14*fReflXYZMuQuad[5][1])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][1]+7.4121097687552e+14*fReflXYZMuQuad[1][1]+4.63256860547201e+14*fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[7][7] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][0]-4.17103711299121e+14*fReflXYZMuQuad[6][0])+5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[5][0]+2384663.0*(fReflXYZMuQuad[4][0]-1.0*fReflXYZMuQuad[3][0])+3.190559553141742e+22*(fReflXYZMuQuad[2][0]-2.0*fReflXYZMuQuad[1][0]+fReflXYZMuQuad[0][0]))); 
  fReflXYQuad[7][8] = 0.05555555555555555*(fReflXYZMuQuad[7][2]+8.0*fReflXYZMuQuad[6][2]+fReflXYZMuQuad[5][2]+8.0*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+fReflXYZMuQuad[2][2]+8.0*fReflXYZMuQuad[1][2]+fReflXYZMuQuad[0][2]); 
  fReflXYQuad[7][9] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][0]+1.6692641e+7*fReflXYZMuQuad[6][0])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][0]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflXYZMuQuad[2][0]+2384663.0*fReflXYZMuQuad[1][0]+3.190559553141742e+22*fReflXYZMuQuad[0][0])-3.20880527843592e+30*(fReflXYZMuQuad[4][0]+fReflXYZMuQuad[3][0])))); 
  fReflXYQuad[7][10] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][1]+9.0*fReflXYZMuQuad[6][1])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][1]-1.346286087882789e+17*fReflXYZMuQuad[5][1])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][1]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][1]+fReflXYZMuQuad[1][1])))); 
  fReflXYQuad[7][11] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYZMuQuad[7][1]-4.17103711299121e+14*fReflXYZMuQuad[6][1])+5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[5][1]+2384663.0*(fReflXYZMuQuad[4][1]-1.0*fReflXYZMuQuad[3][1])+3.190559553141743e+22*(fReflXYZMuQuad[2][1]+fReflXYZMuQuad[0][1]))); 
  fReflXYQuad[7][12] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYZMuQuad[7][2]-4188761.0*fReflXYZMuQuad[6][2])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYZMuQuad[4][2]-4.63256860547201e+14*fReflXYZMuQuad[5][2])+1.6692641e+7*(2.266096151179001e+23*fReflXYZMuQuad[2][2]-1.0*(8377522.0*fReflXYZMuQuad[1][2]+2.266096151179001e+23*fReflXYZMuQuad[0][2])))); 
  fReflXYQuad[7][13] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][0]-2.0*fReflXYZMuQuad[6][0]+fReflXYZMuQuad[5][0])-27.0*fReflXYZMuQuad[3][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[2][0])+2.0*fReflXYZMuQuad[1][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[7][14] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[7][2]+7.4121097687552e+14*fReflXYZMuQuad[6][2]+4.63256860547201e+14*fReflXYZMuQuad[5][2])-1.0*(9.988783372543001e+12*(fReflXYZMuQuad[4][2]+fReflXYZMuQuad[3][2])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYZMuQuad[2][2]+7.4121097687552e+14*fReflXYZMuQuad[1][2]+4.63256860547201e+14*fReflXYZMuQuad[0][2]))); 
  fReflXYQuad[7][15] = 1.719407810605221e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][0]+27.0*fReflXYZMuQuad[6][0]+1.077028870306231e+18*((-1.0*fReflXYZMuQuad[5][0])+2.0*(fReflXYZMuQuad[3][0]-1.0*fReflXYZMuQuad[4][0])+fReflXYZMuQuad[2][0]-1.0*fReflXYZMuQuad[0][0])); 
  fReflXYQuad[7][16] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYZMuQuad[7][1]+1.6692641e+7*fReflXYZMuQuad[6][1])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYZMuQuad[5][1]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflXYZMuQuad[2][1]+2384663.0*fReflXYZMuQuad[1][1]+3.190559553141743e+22*fReflXYZMuQuad[0][1])-3.20880527843592e+30*(fReflXYZMuQuad[4][1]+fReflXYZMuQuad[3][1])))); 
  fReflXYQuad[7][17] = 1.719407810605222e-19*(1.077028870306231e+18*(fReflXYZMuQuad[7][1]+fReflXYZMuQuad[5][1])-27.0*fReflXYZMuQuad[3][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[2][1])+10.0*fReflXYZMuQuad[1][1]-5.0*fReflXYZMuQuad[0][1])); 
  fReflXYQuad[7][18] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYZMuQuad[7][2]+9.0*fReflXYZMuQuad[6][2])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYZMuQuad[4][2]-1.346286087882789e+17*fReflXYZMuQuad[5][2])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYZMuQuad[0][2]-1.0*(9.93730136036331e+15*fReflXYZMuQuad[2][2]+fReflXYZMuQuad[1][2])))); 
  fReflXYQuad[7][19] = 1.719407810605222e-19*(1.077028870306231e+18*fReflXYZMuQuad[7][1]+27.0*fReflXYZMuQuad[6][1]+2.154057740612463e+17*((-5.0*fReflXYZMuQuad[5][1])+2.0*(5.0*fReflXYZMuQuad[3][1]-5.0*fReflXYZMuQuad[4][1])+5.0*fReflXYZMuQuad[2][1])); 
  } 

 
  fRefl[0] = 0.05555555555555555*(fReflXYQuad[7][0]+8.0*fReflXYQuad[6][0]+fReflXYQuad[5][0]+8.0*(fReflXYQuad[4][0]+fReflXYQuad[3][0])+fReflXYQuad[2][0]+8.0*fReflXYQuad[1][0]+fReflXYQuad[0][0]); 
  fRefl[1] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][0]-4188761.0*fReflXYQuad[6][0])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYQuad[4][0]-1.0*fReflXYQuad[3][0])-4.63256860547201e+14*fReflXYQuad[5][0])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][0]-1.0*(8377522.0*fReflXYQuad[1][0]+2.266096151179001e+23*fReflXYQuad[0][0])))); 
  fRefl[2] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][0]+7.4121097687552e+14*fReflXYQuad[6][0]+4.63256860547201e+14*fReflXYQuad[5][0])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][0]+fReflXYQuad[3][0])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][0]+7.4121097687552e+14*fReflXYQuad[1][0]+4.63256860547201e+14*fReflXYQuad[0][0]))); 
  fRefl[3] = 0.05555555555555555*(fReflXYQuad[7][1]+8.0*fReflXYQuad[6][1]+fReflXYQuad[5][1]+8.0*(fReflXYQuad[4][1]+fReflXYQuad[3][1])+fReflXYQuad[2][1]+8.0*fReflXYQuad[1][1]+fReflXYQuad[0][1]); 
  fRefl[4] = 0.05555555555555555*(fReflXYQuad[7][2]+8.0*fReflXYQuad[6][2]+fReflXYQuad[5][2]+8.0*(fReflXYQuad[4][2]+fReflXYQuad[3][2])+fReflXYQuad[2][2]+8.0*fReflXYQuad[1][2]+fReflXYQuad[0][2]); 
  fRefl[5] = 0.05555555555555555*(fReflXYQuad[7][3]+8.0*fReflXYQuad[6][3]+fReflXYQuad[5][3]+8.0*(fReflXYQuad[4][3]+fReflXYQuad[3][3])+fReflXYQuad[2][3]+8.0*fReflXYQuad[1][3]+fReflXYQuad[0][3]); 
  fRefl[6] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][0]+9.0*fReflXYQuad[6][0])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][0]-1.346286087882789e+17*fReflXYQuad[5][0])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][0]-1.0*(9.93730136036331e+15*fReflXYQuad[2][0]+fReflXYQuad[1][0])))); 
  fRefl[7] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][1]-4188761.0*fReflXYQuad[6][1])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYQuad[4][1]-1.0*fReflXYQuad[3][1])-4.63256860547201e+14*fReflXYQuad[5][1])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][1]-1.0*(8377522.0*fReflXYQuad[1][1]+2.266096151179001e+23*fReflXYQuad[0][1])))); 
  fRefl[8] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][1]+7.4121097687552e+14*fReflXYQuad[6][1]+4.63256860547201e+14*fReflXYQuad[5][1])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][1]+fReflXYQuad[3][1])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][1]+7.4121097687552e+14*fReflXYQuad[1][1]+4.63256860547201e+14*fReflXYQuad[0][1]))); 
  fRefl[9] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][2]-4188761.0*fReflXYQuad[6][2])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYQuad[4][2]-1.0*fReflXYQuad[3][2])-4.63256860547201e+14*fReflXYQuad[5][2])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][2]-1.0*(8377522.0*fReflXYQuad[1][2]+2.266096151179001e+23*fReflXYQuad[0][2])))); 
  fRefl[10] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][2]+7.4121097687552e+14*fReflXYQuad[6][2]+4.63256860547201e+14*fReflXYQuad[5][2])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][2]+fReflXYQuad[3][2])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][2]+7.4121097687552e+14*fReflXYQuad[1][2]+4.63256860547201e+14*fReflXYQuad[0][2]))); 
  fRefl[11] = 0.05555555555555555*(fReflXYQuad[7][4]+8.0*fReflXYQuad[6][4]+fReflXYQuad[5][4]+8.0*(fReflXYQuad[4][4]+fReflXYQuad[3][4])+fReflXYQuad[2][4]+8.0*fReflXYQuad[1][4]+fReflXYQuad[0][4]); 
  fRefl[12] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][3]-4188761.0*fReflXYQuad[6][3])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYQuad[4][3]-1.0*fReflXYQuad[3][3])-4.63256860547201e+14*fReflXYQuad[5][3])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][3]-1.0*(8377522.0*fReflXYQuad[1][3]+2.266096151179001e+23*fReflXYQuad[0][3])))); 
  fRefl[13] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][3]+7.4121097687552e+14*fReflXYQuad[6][3]+4.63256860547201e+14*fReflXYQuad[5][3])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][3]+fReflXYQuad[3][3])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][3]+7.4121097687552e+14*fReflXYQuad[1][3]+4.63256860547201e+14*fReflXYQuad[0][3]))); 
  fRefl[14] = 0.05555555555555555*(fReflXYQuad[7][5]+8.0*fReflXYQuad[6][5]+fReflXYQuad[5][5]+8.0*(fReflXYQuad[4][5]+fReflXYQuad[3][5])+fReflXYQuad[2][5]+8.0*fReflXYQuad[1][5]+fReflXYQuad[0][5]); 
  fRefl[15] = 0.05555555555555555*(fReflXYQuad[7][6]+8.0*fReflXYQuad[6][6]+fReflXYQuad[5][6]+8.0*(fReflXYQuad[4][6]+fReflXYQuad[3][6])+fReflXYQuad[2][6]+8.0*fReflXYQuad[1][6]+fReflXYQuad[0][6]); 
  fRefl[16] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYQuad[7][0]-4.17103711299121e+14*fReflXYQuad[6][0])+5.028593299999999e+7*(3.190559553141742e+22*fReflXYQuad[5][0]+2384663.0*(fReflXYQuad[4][0]-1.0*fReflXYQuad[3][0])+3.190559553141742e+22*(fReflXYQuad[2][0]-2.0*fReflXYQuad[1][0]+fReflXYQuad[0][0]))); 
  fRefl[17] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYQuad[7][0]+1.6692641e+7*fReflXYQuad[6][0])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYQuad[5][0]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflXYQuad[2][0]+2384663.0*fReflXYQuad[1][0]+3.190559553141742e+22*fReflXYQuad[0][0])-3.20880527843592e+30*(fReflXYQuad[4][0]+fReflXYQuad[3][0])))); 
  fRefl[18] = 0.05555555555555555*(fReflXYQuad[7][7]+8.0*fReflXYQuad[6][7]+fReflXYQuad[5][7]+8.0*(fReflXYQuad[4][7]+fReflXYQuad[3][7])+fReflXYQuad[2][7]+8.0*fReflXYQuad[1][7]+fReflXYQuad[0][7]); 
  fRefl[19] = 0.05555555555555555*(fReflXYQuad[7][8]+8.0*fReflXYQuad[6][8]+fReflXYQuad[5][8]+8.0*(fReflXYQuad[4][8]+fReflXYQuad[3][8])+fReflXYQuad[2][8]+8.0*fReflXYQuad[1][8]+fReflXYQuad[0][8]); 
  fRefl[20] = 0.05555555555555555*(fReflXYQuad[7][9]+8.0*fReflXYQuad[6][9]+fReflXYQuad[5][9]+8.0*(fReflXYQuad[4][9]+fReflXYQuad[3][9])+fReflXYQuad[2][9]+8.0*fReflXYQuad[1][9]+fReflXYQuad[0][9]); 
  fRefl[21] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][1]+9.0*fReflXYQuad[6][1])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][1]-1.346286087882789e+17*fReflXYQuad[5][1])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][1]-1.0*(9.93730136036331e+15*fReflXYQuad[2][1]+fReflXYQuad[1][1])))); 
  fRefl[22] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][2]+9.0*fReflXYQuad[6][2])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][2]-1.346286087882789e+17*fReflXYQuad[5][2])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][2]-1.0*(9.93730136036331e+15*fReflXYQuad[2][2]+fReflXYQuad[1][2])))); 
  fRefl[23] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][4]-4188761.0*fReflXYQuad[6][4])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYQuad[4][4]-1.0*fReflXYQuad[3][4])-4.63256860547201e+14*fReflXYQuad[5][4])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][4]-1.0*(8377522.0*fReflXYQuad[1][4]+2.266096151179001e+23*fReflXYQuad[0][4])))); 
  fRefl[24] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][4]+7.4121097687552e+14*fReflXYQuad[6][4]+4.63256860547201e+14*fReflXYQuad[5][4])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][4]+fReflXYQuad[3][4])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][4]+7.4121097687552e+14*fReflXYQuad[1][4]+4.63256860547201e+14*fReflXYQuad[0][4]))); 
  fRefl[25] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][3]+9.0*fReflXYQuad[6][3])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][3]-1.346286087882789e+17*fReflXYQuad[5][3])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][3]-1.0*(9.93730136036331e+15*fReflXYQuad[2][3]+fReflXYQuad[1][3])))); 
  fRefl[26] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][5]-4188761.0*fReflXYQuad[6][5])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYQuad[4][5]-1.0*fReflXYQuad[3][5])-4.63256860547201e+14*fReflXYQuad[5][5])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][5]-1.0*(8377522.0*fReflXYQuad[1][5]+2.266096151179001e+23*fReflXYQuad[0][5])))); 
  fRefl[27] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][5]+7.4121097687552e+14*fReflXYQuad[6][5]+4.63256860547201e+14*fReflXYQuad[5][5])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][5]+fReflXYQuad[3][5])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][5]+7.4121097687552e+14*fReflXYQuad[1][5]+4.63256860547201e+14*fReflXYQuad[0][5]))); 
  fRefl[28] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][6]-4188761.0*fReflXYQuad[6][6])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYQuad[4][6]-1.0*fReflXYQuad[3][6])-4.63256860547201e+14*fReflXYQuad[5][6])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][6]-1.0*(8377522.0*fReflXYQuad[1][6]+2.266096151179001e+23*fReflXYQuad[0][6])))); 
  fRefl[29] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][6]+7.4121097687552e+14*fReflXYQuad[6][6]+4.63256860547201e+14*fReflXYQuad[5][6])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][6]+fReflXYQuad[3][6])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][6]+7.4121097687552e+14*fReflXYQuad[1][6]+4.63256860547201e+14*fReflXYQuad[0][6]))); 
  fRefl[30] = 0.05555555555555555*(fReflXYQuad[7][10]+8.0*fReflXYQuad[6][10]+fReflXYQuad[5][10]+8.0*(fReflXYQuad[4][10]+fReflXYQuad[3][10])+fReflXYQuad[2][10]+8.0*fReflXYQuad[1][10]+fReflXYQuad[0][10]); 
  fRefl[31] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflXYQuad[7][0]-2.0*fReflXYQuad[6][0]+fReflXYQuad[5][0])-27.0*fReflXYQuad[3][0]+1.077028870306231e+18*((-1.0*fReflXYQuad[2][0])+2.0*fReflXYQuad[1][0]-1.0*fReflXYQuad[0][0])); 
  fRefl[32] = 1.719407810605221e-19*(1.077028870306231e+18*fReflXYQuad[7][0]+27.0*fReflXYQuad[6][0]+1.077028870306231e+18*((-1.0*fReflXYQuad[5][0])+2.0*(fReflXYQuad[3][0]-1.0*fReflXYQuad[4][0])+fReflXYQuad[2][0]-1.0*fReflXYQuad[0][0])); 
  fRefl[33] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYQuad[7][1]-4.17103711299121e+14*fReflXYQuad[6][1])+5.028593299999999e+7*(3.190559553141743e+22*fReflXYQuad[5][1]+2384663.0*(fReflXYQuad[4][1]-1.0*fReflXYQuad[3][1])+2.127039702094495e+21*(15.0*fReflXYQuad[2][1]+15.0*fReflXYQuad[0][1]))); 
  fRefl[34] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYQuad[7][1]+1.6692641e+7*fReflXYQuad[6][1])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYQuad[5][1]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflXYQuad[2][1]+2384663.0*fReflXYQuad[1][1]+3.190559553141743e+22*fReflXYQuad[0][1])-3.20880527843592e+30*(fReflXYQuad[4][1]+fReflXYQuad[3][1])))); 
  fRefl[35] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][7]-4188761.0*fReflXYQuad[6][7])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYQuad[4][7]-4.63256860547201e+14*fReflXYQuad[5][7])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][7]-1.0*(8377522.0*fReflXYQuad[1][7]+2.266096151179001e+23*fReflXYQuad[0][7])))); 
  fRefl[36] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][7]+7.4121097687552e+14*fReflXYQuad[6][7]+4.63256860547201e+14*fReflXYQuad[5][7])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][7]+fReflXYQuad[3][7])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][7]+7.4121097687552e+14*fReflXYQuad[1][7]+4.63256860547201e+14*fReflXYQuad[0][7]))); 
  fRefl[37] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYQuad[7][2]-4.17103711299121e+14*fReflXYQuad[6][2])+5.028593299999999e+7*(3.190559553141743e+22*fReflXYQuad[5][2]+2384663.0*(fReflXYQuad[4][2]-1.0*fReflXYQuad[3][2])+2.127039702094495e+21*(15.0*fReflXYQuad[2][2]+15.0*fReflXYQuad[0][2]))); 
  fRefl[38] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYQuad[7][2]+1.6692641e+7*fReflXYQuad[6][2])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYQuad[5][2]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflXYQuad[2][2]+2384663.0*fReflXYQuad[1][2]+3.190559553141743e+22*fReflXYQuad[0][2])-3.20880527843592e+30*(fReflXYQuad[4][2]+fReflXYQuad[3][2])))); 
  fRefl[39] = 0.05555555555555555*(fReflXYQuad[7][11]+8.0*fReflXYQuad[6][11]+fReflXYQuad[5][11]+8.0*(fReflXYQuad[4][11]+fReflXYQuad[3][11])+fReflXYQuad[2][11]+8.0*fReflXYQuad[1][11]+fReflXYQuad[0][11]); 
  fRefl[40] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][8]-4188761.0*fReflXYQuad[6][8])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYQuad[4][8]-4.63256860547201e+14*fReflXYQuad[5][8])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][8]-1.0*(8377522.0*fReflXYQuad[1][8]+2.266096151179001e+23*fReflXYQuad[0][8])))); 
  fRefl[41] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][8]+7.4121097687552e+14*fReflXYQuad[6][8]+4.63256860547201e+14*fReflXYQuad[5][8])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][8]+fReflXYQuad[3][8])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][8]+7.4121097687552e+14*fReflXYQuad[1][8]+4.63256860547201e+14*fReflXYQuad[0][8]))); 
  fRefl[42] = 0.05555555555555555*(fReflXYQuad[7][12]+8.0*fReflXYQuad[6][12]+fReflXYQuad[5][12]+8.0*(fReflXYQuad[4][12]+fReflXYQuad[3][12])+fReflXYQuad[2][12]+8.0*fReflXYQuad[1][12]+fReflXYQuad[0][12]); 
  fRefl[43] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYQuad[7][3]-4.17103711299121e+14*fReflXYQuad[6][3])+5.028593299999999e+7*(3.190559553141743e+22*fReflXYQuad[5][3]+2384663.0*(fReflXYQuad[4][3]-1.0*fReflXYQuad[3][3])+2.127039702094495e+21*(15.0*fReflXYQuad[2][3]+15.0*fReflXYQuad[0][3]))); 
  fRefl[44] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYQuad[7][3]+1.6692641e+7*fReflXYQuad[6][3])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYQuad[5][3]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflXYQuad[2][3]+2384663.0*fReflXYQuad[1][3]+3.190559553141743e+22*fReflXYQuad[0][3])-3.20880527843592e+30*(fReflXYQuad[4][3]+fReflXYQuad[3][3])))); 
  fRefl[45] = 0.05555555555555555*(fReflXYQuad[7][13]+8.0*fReflXYQuad[6][13]+fReflXYQuad[5][13]+8.0*(fReflXYQuad[4][13]+fReflXYQuad[3][13])+fReflXYQuad[2][13]+8.0*fReflXYQuad[1][13]+fReflXYQuad[0][13]); 
  fRefl[46] = 0.05555555555555555*(fReflXYQuad[7][14]+8.0*fReflXYQuad[6][14]+fReflXYQuad[5][14]+8.0*(fReflXYQuad[4][14]+fReflXYQuad[3][14])+fReflXYQuad[2][14]+8.0*fReflXYQuad[1][14]+fReflXYQuad[0][14]); 
  fRefl[47] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][9]-4188761.0*fReflXYQuad[6][9])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYQuad[4][9]-4.63256860547201e+14*fReflXYQuad[5][9])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][9]-1.0*(8377522.0*fReflXYQuad[1][9]+2.266096151179001e+23*fReflXYQuad[0][9])))); 
  fRefl[48] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][9]+7.4121097687552e+14*fReflXYQuad[6][9]+4.63256860547201e+14*fReflXYQuad[5][9])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][9]+fReflXYQuad[3][9])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][9]+7.4121097687552e+14*fReflXYQuad[1][9]+4.63256860547201e+14*fReflXYQuad[0][9]))); 
  fRefl[49] = 0.05555555555555555*(fReflXYQuad[7][15]+8.0*fReflXYQuad[6][15]+fReflXYQuad[5][15]+8.0*(fReflXYQuad[4][15]+fReflXYQuad[3][15])+fReflXYQuad[2][15]+8.0*fReflXYQuad[1][15]+fReflXYQuad[0][15]); 
  fRefl[50] = 0.05555555555555555*(fReflXYQuad[7][16]+8.0*fReflXYQuad[6][16]+fReflXYQuad[5][16]+8.0*(fReflXYQuad[4][16]+fReflXYQuad[3][16])+fReflXYQuad[2][16]+8.0*fReflXYQuad[1][16]+fReflXYQuad[0][16]); 
  fRefl[51] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][4]+9.0*fReflXYQuad[6][4])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][4]-1.346286087882789e+17*fReflXYQuad[5][4])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][4]-1.0*(9.93730136036331e+15*fReflXYQuad[2][4]+fReflXYQuad[1][4])))); 
  fRefl[52] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][5]+9.0*fReflXYQuad[6][5])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][5]-1.346286087882789e+17*fReflXYQuad[5][5])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][5]-1.0*(9.93730136036331e+15*fReflXYQuad[2][5]+fReflXYQuad[1][5])))); 
  fRefl[53] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][6]+9.0*fReflXYQuad[6][6])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][6]-1.346286087882789e+17*fReflXYQuad[5][6])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][6]-1.0*(9.93730136036331e+15*fReflXYQuad[2][6]+fReflXYQuad[1][6])))); 
  fRefl[54] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][10]-4188761.0*fReflXYQuad[6][10])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*(fReflXYQuad[4][10]-1.0*fReflXYQuad[3][10])-4.63256860547201e+14*fReflXYQuad[5][10])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][10]-1.0*(8377522.0*fReflXYQuad[1][10]+2.266096151179001e+23*fReflXYQuad[0][10])))); 
  fRefl[55] = 9.295228288552239e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][10]+7.4121097687552e+14*fReflXYQuad[6][10]+4.63256860547201e+14*fReflXYQuad[5][10])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][10]+fReflXYQuad[3][10])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][10]+7.4121097687552e+14*fReflXYQuad[1][10]+4.63256860547201e+14*fReflXYQuad[0][10]))); 
  fRefl[56] = 1.719407810605222e-19*(2.154057740612463e+17*(5.0*fReflXYQuad[7][1]+5.0*fReflXYQuad[5][1])-27.0*fReflXYQuad[3][1]+2.154057740612463e+17*((-5.0*fReflXYQuad[2][1])+10.0*fReflXYQuad[1][1]-5.0*fReflXYQuad[0][1])); 
  fRefl[57] = 1.719407810605222e-19*(1.077028870306231e+18*fReflXYQuad[7][1]+27.0*fReflXYQuad[6][1]+2.154057740612463e+17*((-5.0*fReflXYQuad[5][1])+2.0*(5.0*fReflXYQuad[3][1]-5.0*fReflXYQuad[4][1])+5.0*fReflXYQuad[2][1])); 
  fRefl[58] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][7]+9.0*fReflXYQuad[6][7])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][7]-1.346286087882789e+17*fReflXYQuad[5][7])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][7]-1.0*(9.93730136036331e+15*fReflXYQuad[2][7]+fReflXYQuad[1][7])))); 
  fRefl[59] = 1.719407810605222e-19*(2.154057740612463e+17*(5.0*fReflXYQuad[7][2]+5.0*fReflXYQuad[5][2])-27.0*fReflXYQuad[3][2]+2.154057740612463e+17*((-5.0*fReflXYQuad[2][2])+10.0*fReflXYQuad[1][2]-5.0*fReflXYQuad[0][2])); 
  fRefl[60] = 1.719407810605222e-19*(1.077028870306231e+18*fReflXYQuad[7][2]+27.0*fReflXYQuad[6][2]+2.154057740612463e+17*((-5.0*fReflXYQuad[5][2])+2.0*(5.0*fReflXYQuad[3][2]-5.0*fReflXYQuad[4][2])+5.0*fReflXYQuad[2][2])); 
  fRefl[61] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYQuad[7][4]-4.17103711299121e+14*fReflXYQuad[6][4])+5.028593299999999e+7*(3.190559553141742e+22*fReflXYQuad[5][4]+2384663.0*(fReflXYQuad[4][4]-1.0*fReflXYQuad[3][4])+3.190559553141742e+22*(fReflXYQuad[2][4]-2.0*fReflXYQuad[1][4]+fReflXYQuad[0][4]))); 
  fRefl[62] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYQuad[7][4]+1.6692641e+7*fReflXYQuad[6][4])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYQuad[5][4]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflXYQuad[2][4]+2384663.0*fReflXYQuad[1][4]+3.190559553141742e+22*fReflXYQuad[0][4])-3.20880527843592e+30*(fReflXYQuad[4][4]+fReflXYQuad[3][4])))); 
  fRefl[63] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][11]-4188761.0*fReflXYQuad[6][11])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYQuad[4][11]-4.63256860547201e+14*fReflXYQuad[5][11])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][11]-1.0*(8377522.0*fReflXYQuad[1][11]+2.266096151179001e+23*fReflXYQuad[0][11])))); 
  fRefl[64] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][11]+7.4121097687552e+14*fReflXYQuad[6][11]+4.63256860547201e+14*fReflXYQuad[5][11])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][11]+fReflXYQuad[3][11])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][11]+7.4121097687552e+14*fReflXYQuad[1][11]+4.63256860547201e+14*fReflXYQuad[0][11]))); 
  fRefl[65] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][8]+9.0*fReflXYQuad[6][8])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][8]-1.346286087882789e+17*fReflXYQuad[5][8])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][8]-1.0*(9.93730136036331e+15*fReflXYQuad[2][8]+fReflXYQuad[1][8])))); 
  fRefl[66] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][12]-4188761.0*fReflXYQuad[6][12])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYQuad[4][12]-4.63256860547201e+14*fReflXYQuad[5][12])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][12]-1.0*(8377522.0*fReflXYQuad[1][12]+2.266096151179001e+23*fReflXYQuad[0][12])))); 
  fRefl[67] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][12]+7.4121097687552e+14*fReflXYQuad[6][12]+4.63256860547201e+14*fReflXYQuad[5][12])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][12]+fReflXYQuad[3][12])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][12]+7.4121097687552e+14*fReflXYQuad[1][12]+4.63256860547201e+14*fReflXYQuad[0][12]))); 
  fRefl[68] = 1.719407810605222e-19*(2.154057740612463e+17*(5.0*fReflXYQuad[7][3]+5.0*fReflXYQuad[5][3])-27.0*fReflXYQuad[3][3]+2.154057740612463e+17*((-5.0*fReflXYQuad[2][3])+10.0*fReflXYQuad[1][3]-5.0*fReflXYQuad[0][3])); 
  fRefl[69] = 1.719407810605222e-19*(1.077028870306231e+18*fReflXYQuad[7][3]+27.0*fReflXYQuad[6][3]+2.154057740612463e+17*((-5.0*fReflXYQuad[5][3])+2.0*(5.0*fReflXYQuad[3][3]-5.0*fReflXYQuad[4][3])+5.0*fReflXYQuad[2][3])); 
  fRefl[70] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYQuad[7][5]-4.17103711299121e+14*fReflXYQuad[6][5])+5.028593299999999e+7*(3.190559553141742e+22*fReflXYQuad[5][5]+2384663.0*(fReflXYQuad[4][5]-1.0*fReflXYQuad[3][5])+3.190559553141742e+22*(fReflXYQuad[2][5]-2.0*fReflXYQuad[1][5]+fReflXYQuad[0][5]))); 
  fRefl[71] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYQuad[7][5]+1.6692641e+7*fReflXYQuad[6][5])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYQuad[5][5]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflXYQuad[2][5]+2384663.0*fReflXYQuad[1][5]+3.190559553141742e+22*fReflXYQuad[0][5])-3.20880527843592e+30*(fReflXYQuad[4][5]+fReflXYQuad[3][5])))); 
  fRefl[72] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][13]-4188761.0*fReflXYQuad[6][13])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYQuad[4][13]-4.63256860547201e+14*fReflXYQuad[5][13])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][13]-1.0*(8377522.0*fReflXYQuad[1][13]+2.266096151179001e+23*fReflXYQuad[0][13])))); 
  fRefl[73] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][13]+7.4121097687552e+14*fReflXYQuad[6][13]+4.63256860547201e+14*fReflXYQuad[5][13])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][13]+fReflXYQuad[3][13])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][13]+7.4121097687552e+14*fReflXYQuad[1][13]+4.63256860547201e+14*fReflXYQuad[0][13]))); 
  fRefl[74] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYQuad[7][6]-4.17103711299121e+14*fReflXYQuad[6][6])+5.028593299999999e+7*(3.190559553141742e+22*fReflXYQuad[5][6]+2384663.0*(fReflXYQuad[4][6]-1.0*fReflXYQuad[3][6])+3.190559553141742e+22*(fReflXYQuad[2][6]-2.0*fReflXYQuad[1][6]+fReflXYQuad[0][6]))); 
  fRefl[75] = 6.326811562591036e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYQuad[7][6]+1.6692641e+7*fReflXYQuad[6][6])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYQuad[5][6]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141742e+22*fReflXYQuad[2][6]+2384663.0*fReflXYQuad[1][6]+3.190559553141742e+22*fReflXYQuad[0][6])-3.20880527843592e+30*(fReflXYQuad[4][6]+fReflXYQuad[3][6])))); 
  fRefl[76] = 0.05555555555555555*(fReflXYQuad[7][17]+8.0*fReflXYQuad[6][17]+fReflXYQuad[5][17]+8.0*(fReflXYQuad[4][17]+fReflXYQuad[3][17])+fReflXYQuad[2][17]+8.0*fReflXYQuad[1][17]+fReflXYQuad[0][17]); 
  fRefl[77] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][14]-4188761.0*fReflXYQuad[6][14])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYQuad[4][14]-4.63256860547201e+14*fReflXYQuad[5][14])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][14]-1.0*(8377522.0*fReflXYQuad[1][14]+2.266096151179001e+23*fReflXYQuad[0][14])))); 
  fRefl[78] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][14]+7.4121097687552e+14*fReflXYQuad[6][14]+4.63256860547201e+14*fReflXYQuad[5][14])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][14]+fReflXYQuad[3][14])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][14]+7.4121097687552e+14*fReflXYQuad[1][14]+4.63256860547201e+14*fReflXYQuad[0][14]))); 
  fRefl[79] = 0.05555555555555555*(fReflXYQuad[7][18]+8.0*fReflXYQuad[6][18]+fReflXYQuad[5][18]+8.0*(fReflXYQuad[4][18]+fReflXYQuad[3][18])+fReflXYQuad[2][18]+8.0*fReflXYQuad[1][18]+fReflXYQuad[0][18]); 
  fRefl[80] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][9]+9.0*fReflXYQuad[6][9])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][9]-1.346286087882789e+17*fReflXYQuad[5][9])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][9]-1.0*(9.93730136036331e+15*fReflXYQuad[2][9]+fReflXYQuad[1][9])))); 
  fRefl[81] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][15]-4188761.0*fReflXYQuad[6][15])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYQuad[4][15]-4.63256860547201e+14*fReflXYQuad[5][15])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][15]-1.0*(8377522.0*fReflXYQuad[1][15]+2.266096151179001e+23*fReflXYQuad[0][15])))); 
  fRefl[82] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][15]+7.4121097687552e+14*fReflXYQuad[6][15]+4.63256860547201e+14*fReflXYQuad[5][15])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][15]+fReflXYQuad[3][15])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][15]+7.4121097687552e+14*fReflXYQuad[1][15]+4.63256860547201e+14*fReflXYQuad[0][15]))); 
  fRefl[83] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][16]-4188761.0*fReflXYQuad[6][16])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYQuad[4][16]-4.63256860547201e+14*fReflXYQuad[5][16])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][16]-1.0*(8377522.0*fReflXYQuad[1][16]+2.266096151179001e+23*fReflXYQuad[0][16])))); 
  fRefl[84] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][16]+7.4121097687552e+14*fReflXYQuad[6][16]+4.63256860547201e+14*fReflXYQuad[5][16])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][16]+fReflXYQuad[3][16])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][16]+7.4121097687552e+14*fReflXYQuad[1][16]+4.63256860547201e+14*fReflXYQuad[0][16]))); 
  fRefl[85] = 0.05555555555555555*(fReflXYQuad[7][19]+8.0*fReflXYQuad[6][19]+fReflXYQuad[5][19]+8.0*(fReflXYQuad[4][19]+fReflXYQuad[3][19])+fReflXYQuad[2][19]+8.0*fReflXYQuad[1][19]+fReflXYQuad[0][19]); 
  fRefl[86] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][10]+9.0*fReflXYQuad[6][10])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][10]-1.346286087882789e+17*fReflXYQuad[5][10])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][10]-1.0*(9.93730136036331e+15*fReflXYQuad[2][10]+fReflXYQuad[1][10])))); 
  fRefl[87] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflXYQuad[7][4]-2.0*fReflXYQuad[6][4]+fReflXYQuad[5][4])-27.0*fReflXYQuad[3][4]+1.077028870306231e+18*((-1.0*fReflXYQuad[2][4])+2.0*fReflXYQuad[1][4]-1.0*fReflXYQuad[0][4])); 
  fRefl[88] = 1.719407810605221e-19*(1.077028870306231e+18*fReflXYQuad[7][4]+27.0*fReflXYQuad[6][4]+1.077028870306231e+18*((-1.0*fReflXYQuad[5][4])+2.0*(fReflXYQuad[3][4]-1.0*fReflXYQuad[4][4])+fReflXYQuad[2][4]-1.0*fReflXYQuad[0][4])); 
  fRefl[89] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][11]+9.0*fReflXYQuad[6][11])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][11]-1.346286087882789e+17*fReflXYQuad[5][11])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][11]-1.0*(9.93730136036331e+15*fReflXYQuad[2][11]+fReflXYQuad[1][11])))); 
  fRefl[90] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][12]+9.0*fReflXYQuad[6][12])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][12]-1.346286087882789e+17*fReflXYQuad[5][12])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][12]-1.0*(9.93730136036331e+15*fReflXYQuad[2][12]+fReflXYQuad[1][12])))); 
  fRefl[91] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflXYQuad[7][5]-2.0*fReflXYQuad[6][5]+fReflXYQuad[5][5])-27.0*fReflXYQuad[3][5]+1.077028870306231e+18*((-1.0*fReflXYQuad[2][5])+2.0*fReflXYQuad[1][5]-1.0*fReflXYQuad[0][5])); 
  fRefl[92] = 1.719407810605221e-19*(1.077028870306231e+18*fReflXYQuad[7][5]+27.0*fReflXYQuad[6][5]+1.077028870306231e+18*((-1.0*fReflXYQuad[5][5])+2.0*(fReflXYQuad[3][5]-1.0*fReflXYQuad[4][5])+fReflXYQuad[2][5]-1.0*fReflXYQuad[0][5])); 
  fRefl[93] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][13]+9.0*fReflXYQuad[6][13])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][13]-1.346286087882789e+17*fReflXYQuad[5][13])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][13]-1.0*(9.93730136036331e+15*fReflXYQuad[2][13]+fReflXYQuad[1][13])))); 
  fRefl[94] = 1.719407810605221e-19*(1.077028870306231e+18*(fReflXYQuad[7][6]-2.0*fReflXYQuad[6][6]+fReflXYQuad[5][6])-27.0*fReflXYQuad[3][6]+1.077028870306231e+18*((-1.0*fReflXYQuad[2][6])+2.0*fReflXYQuad[1][6]-1.0*fReflXYQuad[0][6])); 
  fRefl[95] = 1.719407810605221e-19*(1.077028870306231e+18*fReflXYQuad[7][6]+27.0*fReflXYQuad[6][6]+1.077028870306231e+18*((-1.0*fReflXYQuad[5][6])+2.0*(fReflXYQuad[3][6]-1.0*fReflXYQuad[4][6])+fReflXYQuad[2][6]-1.0*fReflXYQuad[0][6])); 
  fRefl[96] = 1.548563879333187e-31*(7.693063359330223e+15*(2.0855185564956e+14*fReflXYQuad[7][10]-4.17103711299121e+14*fReflXYQuad[6][10])+5.028593299999999e+7*(3.190559553141743e+22*fReflXYQuad[5][10]+2384663.0*(fReflXYQuad[4][10]-1.0*fReflXYQuad[3][10])+2.127039702094495e+21*(15.0*fReflXYQuad[2][10]+15.0*fReflXYQuad[0][10]))); 
  fRefl[97] = 6.32681156259104e-55*(9.450856442503457e+29*(4.155147325021804e+23*fReflXYQuad[7][10]+1.6692641e+7*fReflXYQuad[6][10])+2.087265252531456e+15*(1.881394845173929e+38*fReflXYQuad[5][10]+1.17264507e+8*(5.028593299999999e+7*(3.190559553141743e+22*fReflXYQuad[2][10]+2384663.0*fReflXYQuad[1][10]+3.190559553141743e+22*fReflXYQuad[0][10])-3.20880527843592e+30*(fReflXYQuad[4][10]+fReflXYQuad[3][10])))); 
  fRefl[98] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][17]-4188761.0*fReflXYQuad[6][17])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYQuad[4][17]-4.63256860547201e+14*fReflXYQuad[5][17])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][17]-1.0*(8377522.0*fReflXYQuad[1][17]+2.266096151179001e+23*fReflXYQuad[0][17])))); 
  fRefl[99] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][17]+7.4121097687552e+14*fReflXYQuad[6][17]+4.63256860547201e+14*fReflXYQuad[5][17])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][17]+fReflXYQuad[3][17])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][17]+7.4121097687552e+14*fReflXYQuad[1][17]+4.63256860547201e+14*fReflXYQuad[0][17]))); 
  fRefl[100] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][14]+9.0*fReflXYQuad[6][14])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][14]-1.346286087882789e+17*fReflXYQuad[5][14])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][14]-1.0*(9.93730136036331e+15*fReflXYQuad[2][14]+fReflXYQuad[1][14])))); 
  fRefl[101] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][18]-4188761.0*fReflXYQuad[6][18])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYQuad[4][18]-4.63256860547201e+14*fReflXYQuad[5][18])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][18]-1.0*(8377522.0*fReflXYQuad[1][18]+2.266096151179001e+23*fReflXYQuad[0][18])))); 
  fRefl[102] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][18]+7.4121097687552e+14*fReflXYQuad[6][18]+4.63256860547201e+14*fReflXYQuad[5][18])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][18]+fReflXYQuad[3][18])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][18]+7.4121097687552e+14*fReflXYQuad[1][18]+4.63256860547201e+14*fReflXYQuad[0][18]))); 
  fRefl[103] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][15]+9.0*fReflXYQuad[6][15])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][15]-1.346286087882789e+17*fReflXYQuad[5][15])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][15]-1.0*(9.93730136036331e+15*fReflXYQuad[2][15]+fReflXYQuad[1][15])))); 
  fRefl[104] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][16]+9.0*fReflXYQuad[6][16])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][16]-1.346286087882789e+17*fReflXYQuad[5][16])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][16]-1.0*(9.93730136036331e+15*fReflXYQuad[2][16]+fReflXYQuad[1][16])))); 
  fRefl[105] = 1.138357133969994e-46*(1.947190939890761e+22*(9.340587065745829e+22*fReflXYQuad[7][19]-4188761.0*fReflXYQuad[6][19])+4.80816459958139e+14*(8.165476379223229e+15*(7.4121097687552e+14*fReflXYQuad[4][19]-4.63256860547201e+14*fReflXYQuad[5][19])+1.6692641e+7*(2.266096151179001e+23*fReflXYQuad[2][19]-1.0*(8377522.0*fReflXYQuad[1][19]+2.266096151179001e+23*fReflXYQuad[0][19])))); 
  fRefl[106] = 9.295228288552242e-31*(4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[7][19]+7.4121097687552e+14*fReflXYQuad[6][19]+4.63256860547201e+14*fReflXYQuad[5][19])-1.0*(9.988783372543001e+12*(fReflXYQuad[4][19]+fReflXYQuad[3][19])+4.80816459958139e+14*(4.63256860547201e+14*fReflXYQuad[2][19]+7.4121097687552e+14*fReflXYQuad[1][19]+4.63256860547201e+14*fReflXYQuad[0][19]))); 
  fRefl[107] = 1.719407810605222e-19*(2.154057740612463e+17*(5.0*fReflXYQuad[7][10]+5.0*fReflXYQuad[5][10])-27.0*fReflXYQuad[3][10]+2.154057740612463e+17*((-5.0*fReflXYQuad[2][10])+10.0*fReflXYQuad[1][10]-5.0*fReflXYQuad[0][10])); 
  fRefl[108] = 1.719407810605222e-19*(1.077028870306231e+18*fReflXYQuad[7][10]+27.0*fReflXYQuad[6][10]+2.154057740612463e+17*((-5.0*fReflXYQuad[5][10])+2.0*(5.0*fReflXYQuad[3][10]-5.0*fReflXYQuad[4][10])+5.0*fReflXYQuad[2][10])); 
  fRefl[109] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][17]+9.0*fReflXYQuad[6][17])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][17]-1.346286087882789e+17*fReflXYQuad[5][17])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][17]-1.0*(9.93730136036331e+15*fReflXYQuad[2][17]+fReflXYQuad[1][17])))); 
  fRefl[110] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][18]+9.0*fReflXYQuad[6][18])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][18]-1.346286087882789e+17*fReflXYQuad[5][18])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][18]-1.0*(9.93730136036331e+15*fReflXYQuad[2][18]+fReflXYQuad[1][18])))); 
  fRefl[111] = 6.325146169955158e-49*(3.34461264313896e+30*(1.31304952186069e+17*fReflXYQuad[7][19]+9.0*fReflXYQuad[6][19])+3.282623804651726e+15*(9.93730136036331e+14*(9.0*fReflXYQuad[4][19]-1.346286087882789e+17*fReflXYQuad[5][19])+1.346286087882789e+16*(9.93730136036331e+15*fReflXYQuad[0][19]-1.0*(9.93730136036331e+15*fReflXYQuad[2][19]+fReflXYQuad[1][19])))); 
}
