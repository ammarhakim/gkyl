#include <cmath> 
#include <GyrokineticModDecl.h> 
// approximation for inverse Langevin function 
double invL(double x) { 
  // from Kroger 
  return (3.*x-x*x*x*(6. + x*x - 2.*x*x*x*x)/5.)/(1.-x*x); 
}

void calcSheathReflection3x2vSer_P1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl) 
{ 
  double vcutSq_i, xc, b; 
  double fReflXYQuad[4][4]; 
  double fReflXYMuQuad[2][2]; 
  

// quadrature node (x,y)_i=1 
  vcutSq_i = (0.5*q_*((2.449489742783178*phiWall[7]-2.449489742783178*(phi[7]+phiWall[6])+2.449489742783178*phi[6]-2.449489742783178*phiWall[5]+2.449489742783178*(phi[5]+phiWall[3])-2.449489742783178*phi[3])*zVal+1.414213562373095*phiWall[4]-1.414213562373095*(phi[4]+phiWall[2])+1.414213562373095*phi[2]-1.414213562373095*phiWall[1]+1.414213562373095*(phi[1]+phiWall[0])-1.414213562373095*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[0][0] = 0.0; 
  fReflXYQuad[0][1] = 0.0; 
  fReflXYQuad[0][2] = 0.0; 
  fReflXYQuad[0][3] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[0][0] = 0.3535533905932737*(1.732050807568877*(f[16]-1.0*(f[8]+f[7])+f[3])*zVal+f[6]-1.0*(f[2]+f[1])+f[0]); 
  fReflXYQuad[0][1] = 0.3535533905932737*(1.732050807568877*(f[26]-1.0*(f[19]+f[18])+f[11])*zVal+f[17]-1.0*(f[10]+f[9])+f[4]); 
  fReflXYQuad[0][2] = 0.3535533905932737*(1.732050807568877*(f[27]-1.0*(f[22]+f[21])+f[14])*zVal+f[20]-1.0*(f[13]+f[12])+f[5]); 
  fReflXYQuad[0][3] = 0.3535533905932737*(1.732050807568877*(f[31]-1.0*(f[30]+f[29])+f[25])*zVal+f[28]-1.0*(f[24]+f[23])+f[15]); 
  } else { // partial reflection 
  b = invL((0.5773502691896258*(1.732050807568877*(f[31]-1.0*(f[30]+f[29]+f[26])+f[25]+f[19]+f[18]-1.0*f[11])*zVal+f[28]-1.0*(f[24]+f[23]+f[17])+f[15]+f[10]+f[9]-1.0*f[4]))/(1.732050807568877*(f[27]-1.0*(f[22]+f[21]+f[16])+f[14]+f[8]+f[7]-1.0*f[3])*zVal+f[20]-1.0*(f[13]+f[12]+f[6])+f[5]+f[2]+f[1]-1.0*f[0])); 
  if(wv > 0) {
  xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
  fReflXYMuQuad[0][0] = -0.25*(1.732050807568877*(f[27]-1.0*(f[22]+f[21]+f[16])+f[14]+f[8]+f[7]-1.0*f[3])*zVal+f[20]-1.0*(f[13]+f[12]+f[6])+f[5]+f[2]+f[1]-1.0*f[0])*(exp(b*xc)-exp(-b))/sinh(b); 
  fReflXYMuQuad[0][1] = -0.25*(1.732050807568877*(f[31]-1.0*(f[30]+f[29]+f[26])+f[25]+f[19]+f[18]-1.0*f[11])*zVal+f[28]-1.0*(f[24]+f[23]+f[17])+f[15]+f[10]+f[9]-1.0*f[4])*((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/(b*cosh(b)-sinh(b)); 
  } else { 
  xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
  fReflXYMuQuad[0][0] = -0.25*(1.732050807568877*(f[27]-1.0*(f[22]+f[21]+f[16])+f[14]+f[8]+f[7]-1.0*f[3])*zVal+f[20]-1.0*(f[13]+f[12]+f[6])+f[5]+f[2]+f[1]-1.0*f[0])*(exp(b)-exp(b*xc))/sinh(b); 
  fReflXYMuQuad[0][1] = -0.25*(1.732050807568877*(f[31]-1.0*(f[30]+f[29]+f[26])+f[25]+f[19]+f[18]-1.0*f[11])*zVal+f[28]-1.0*(f[24]+f[23]+f[17])+f[15]+f[10]+f[9]-1.0*f[4])*((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/(b*cosh(b)-sinh(b)); 
  } 
  b = invL((0.5773502691896258*(1.732050807568877*(f[31]-1.0*(f[30]+f[29])+f[26]+f[25]-1.0*(f[19]+f[18])+f[11])*zVal+f[28]-1.0*(f[24]+f[23])+f[17]+f[15]-1.0*(f[10]+f[9])+f[4]))/(1.732050807568877*(f[27]-1.0*(f[22]+f[21])+f[16]+f[14]-1.0*(f[8]+f[7])+f[3])*zVal+f[20]-1.0*(f[13]+f[12])+f[6]+f[5]-1.0*(f[2]+f[1])+f[0])); 
  if(wv > 0) {
  xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
  fReflXYMuQuad[1][0] = 0.25*(1.732050807568877*(f[27]-1.0*(f[22]+f[21])+f[16]+f[14]-1.0*(f[8]+f[7])+f[3])*zVal+f[20]-1.0*(f[13]+f[12])+f[6]+f[5]-1.0*(f[2]+f[1])+f[0])*(exp(b*xc)-exp(-b))/sinh(b); 
  fReflXYMuQuad[1][1] = 0.25*(1.732050807568877*(f[31]-1.0*(f[30]+f[29])+f[26]+f[25]-1.0*(f[19]+f[18])+f[11])*zVal+f[28]-1.0*(f[24]+f[23])+f[17]+f[15]-1.0*(f[10]+f[9])+f[4])*((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/(b*cosh(b)-sinh(b)); 
  } else { 
  xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
  fReflXYMuQuad[1][0] = 0.25*(1.732050807568877*(f[27]-1.0*(f[22]+f[21])+f[16]+f[14]-1.0*(f[8]+f[7])+f[3])*zVal+f[20]-1.0*(f[13]+f[12])+f[6]+f[5]-1.0*(f[2]+f[1])+f[0])*(exp(b)-exp(b*xc))/sinh(b); 
  fReflXYMuQuad[1][1] = 0.25*(1.732050807568877*(f[31]-1.0*(f[30]+f[29])+f[26]+f[25]-1.0*(f[19]+f[18])+f[11])*zVal+f[28]-1.0*(f[24]+f[23])+f[17]+f[15]-1.0*(f[10]+f[9])+f[4])*((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/(b*cosh(b)-sinh(b)); 
  } 
  fReflXYQuad[0][0] = 0.7071067811865468*(fReflXYMuQuad[1][0]+fReflXYMuQuad[0][0]); 
  fReflXYQuad[0][1] = 0.7071067811865468*(fReflXYMuQuad[1][1]+fReflXYMuQuad[0][1]); 
  fReflXYQuad[0][2] = 0.7071067811865468*(fReflXYMuQuad[1][0]-1.0*fReflXYMuQuad[0][0]); 
  fReflXYQuad[0][3] = 0.7071067811865468*(fReflXYMuQuad[1][1]-1.0*fReflXYMuQuad[0][1]); 
  } 

 
// quadrature node (x,y)_i=2 
  vcutSq_i = -(0.5*q_*((2.449489742783178*phiWall[7]-2.449489742783178*phi[7]+2.449489742783178*phiWall[6]-2.449489742783178*(phi[6]+phiWall[5])+2.449489742783178*phi[5]-2.449489742783178*phiWall[3]+2.449489742783178*phi[3])*zVal+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]+1.414213562373095*phiWall[2]-1.414213562373095*(phi[2]+phiWall[1])+1.414213562373095*phi[1]-1.414213562373095*phiWall[0]+1.414213562373095*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[1][0] = 0.0; 
  fReflXYQuad[1][1] = 0.0; 
  fReflXYQuad[1][2] = 0.0; 
  fReflXYQuad[1][3] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[1][0] = -0.3535533905932737*(1.732050807568877*(f[16]+f[8]-1.0*(f[7]+f[3]))*zVal+f[6]+f[2]-1.0*(f[1]+f[0])); 
  fReflXYQuad[1][1] = -0.3535533905932737*(1.732050807568877*(f[26]+f[19]-1.0*(f[18]+f[11]))*zVal+f[17]+f[10]-1.0*(f[9]+f[4])); 
  fReflXYQuad[1][2] = -0.3535533905932737*(1.732050807568877*(f[27]+f[22]-1.0*(f[21]+f[14]))*zVal+f[20]+f[13]-1.0*(f[12]+f[5])); 
  fReflXYQuad[1][3] = -0.3535533905932737*(1.732050807568877*(f[31]+f[30]-1.0*(f[29]+f[25]))*zVal+f[28]+f[24]-1.0*(f[23]+f[15])); 
  } else { // partial reflection 
  b = invL((0.5773502691896258*(1.732050807568877*(f[31]+f[30]-1.0*(f[29]+f[26]+f[25]+f[19])+f[18]+f[11])*zVal+f[28]+f[24]-1.0*(f[23]+f[17]+f[15]+f[10])+f[9]+f[4]))/(1.732050807568877*(f[27]+f[22]-1.0*(f[21]+f[16]+f[14]+f[8])+f[7]+f[3])*zVal+f[20]+f[13]-1.0*(f[12]+f[6]+f[5]+f[2])+f[1]+f[0])); 
  if(wv > 0) {
  xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
  fReflXYMuQuad[0][0] = 0.25*(1.732050807568877*(f[27]+f[22]-1.0*(f[21]+f[16]+f[14]+f[8])+f[7]+f[3])*zVal+f[20]+f[13]-1.0*(f[12]+f[6]+f[5]+f[2])+f[1]+f[0])*(exp(b*xc)-exp(-b))/sinh(b); 
  fReflXYMuQuad[0][1] = 0.25*(1.732050807568877*(f[31]+f[30]-1.0*(f[29]+f[26]+f[25]+f[19])+f[18]+f[11])*zVal+f[28]+f[24]-1.0*(f[23]+f[17]+f[15]+f[10])+f[9]+f[4])*((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/(b*cosh(b)-sinh(b)); 
  } else { 
  xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
  fReflXYMuQuad[0][0] = 0.25*(1.732050807568877*(f[27]+f[22]-1.0*(f[21]+f[16]+f[14]+f[8])+f[7]+f[3])*zVal+f[20]+f[13]-1.0*(f[12]+f[6]+f[5]+f[2])+f[1]+f[0])*(exp(b)-exp(b*xc))/sinh(b); 
  fReflXYMuQuad[0][1] = 0.25*(1.732050807568877*(f[31]+f[30]-1.0*(f[29]+f[26]+f[25]+f[19])+f[18]+f[11])*zVal+f[28]+f[24]-1.0*(f[23]+f[17]+f[15]+f[10])+f[9]+f[4])*((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/(b*cosh(b)-sinh(b)); 
  } 
  b = invL((0.5773502691896258*(1.732050807568877*(f[31]+f[30]-1.0*f[29]+f[26]-1.0*f[25]+f[19]-1.0*(f[18]+f[11]))*zVal+f[28]+f[24]-1.0*f[23]+f[17]-1.0*f[15]+f[10]-1.0*(f[9]+f[4])))/(1.732050807568877*(f[27]+f[22]-1.0*f[21]+f[16]-1.0*f[14]+f[8]-1.0*(f[7]+f[3]))*zVal+f[20]+f[13]-1.0*f[12]+f[6]-1.0*f[5]+f[2]-1.0*(f[1]+f[0]))); 
  if(wv > 0) {
  xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
  fReflXYMuQuad[1][0] = -0.25*(1.732050807568877*(f[27]+f[22]-1.0*f[21]+f[16]-1.0*f[14]+f[8]-1.0*(f[7]+f[3]))*zVal+f[20]+f[13]-1.0*f[12]+f[6]-1.0*f[5]+f[2]-1.0*(f[1]+f[0]))*(exp(b*xc)-exp(-b))/sinh(b); 
  fReflXYMuQuad[1][1] = -0.25*(1.732050807568877*(f[31]+f[30]-1.0*f[29]+f[26]-1.0*f[25]+f[19]-1.0*(f[18]+f[11]))*zVal+f[28]+f[24]-1.0*f[23]+f[17]-1.0*f[15]+f[10]-1.0*(f[9]+f[4]))*((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/(b*cosh(b)-sinh(b)); 
  } else { 
  xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
  fReflXYMuQuad[1][0] = -0.25*(1.732050807568877*(f[27]+f[22]-1.0*f[21]+f[16]-1.0*f[14]+f[8]-1.0*(f[7]+f[3]))*zVal+f[20]+f[13]-1.0*f[12]+f[6]-1.0*f[5]+f[2]-1.0*(f[1]+f[0]))*(exp(b)-exp(b*xc))/sinh(b); 
  fReflXYMuQuad[1][1] = -0.25*(1.732050807568877*(f[31]+f[30]-1.0*f[29]+f[26]-1.0*f[25]+f[19]-1.0*(f[18]+f[11]))*zVal+f[28]+f[24]-1.0*f[23]+f[17]-1.0*f[15]+f[10]-1.0*(f[9]+f[4]))*((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/(b*cosh(b)-sinh(b)); 
  } 
  fReflXYQuad[1][0] = 0.7071067811865468*(fReflXYMuQuad[1][0]+fReflXYMuQuad[0][0]); 
  fReflXYQuad[1][1] = 0.7071067811865468*(fReflXYMuQuad[1][1]+fReflXYMuQuad[0][1]); 
  fReflXYQuad[1][2] = 0.7071067811865468*(fReflXYMuQuad[1][0]-1.0*fReflXYMuQuad[0][0]); 
  fReflXYQuad[1][3] = 0.7071067811865468*(fReflXYMuQuad[1][1]-1.0*fReflXYMuQuad[0][1]); 
  } 

 
// quadrature node (x,y)_i=3 
  vcutSq_i = -(0.5*q_*((2.449489742783178*phiWall[7]-2.449489742783178*(phi[7]+phiWall[6])+2.449489742783178*(phi[6]+phiWall[5])-2.449489742783178*(phi[5]+phiWall[3])+2.449489742783178*phi[3])*zVal+1.414213562373095*phiWall[4]-1.414213562373095*(phi[4]+phiWall[2])+1.414213562373095*(phi[2]+phiWall[1])-1.414213562373095*(phi[1]+phiWall[0])+1.414213562373095*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[2][0] = 0.0; 
  fReflXYQuad[2][1] = 0.0; 
  fReflXYQuad[2][2] = 0.0; 
  fReflXYQuad[2][3] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[2][0] = -0.3535533905932737*(1.732050807568877*(f[16]-1.0*f[8]+f[7]-1.0*f[3])*zVal+f[6]-1.0*f[2]+f[1]-1.0*f[0]); 
  fReflXYQuad[2][1] = -0.3535533905932737*(1.732050807568877*(f[26]-1.0*f[19]+f[18]-1.0*f[11])*zVal+f[17]-1.0*f[10]+f[9]-1.0*f[4]); 
  fReflXYQuad[2][2] = -0.3535533905932737*(1.732050807568877*(f[27]-1.0*f[22]+f[21]-1.0*f[14])*zVal+f[20]-1.0*f[13]+f[12]-1.0*f[5]); 
  fReflXYQuad[2][3] = -0.3535533905932737*(1.732050807568877*(f[31]-1.0*f[30]+f[29]-1.0*f[25])*zVal+f[28]-1.0*f[24]+f[23]-1.0*f[15]); 
  } else { // partial reflection 
  b = invL((0.5773502691896258*(1.732050807568877*(f[31]-1.0*f[30]+f[29]-1.0*(f[26]+f[25])+f[19]-1.0*f[18]+f[11])*zVal+f[28]-1.0*f[24]+f[23]-1.0*(f[17]+f[15])+f[10]-1.0*f[9]+f[4]))/(1.732050807568877*(f[27]-1.0*f[22]+f[21]-1.0*(f[16]+f[14])+f[8]-1.0*f[7]+f[3])*zVal+f[20]-1.0*f[13]+f[12]-1.0*(f[6]+f[5])+f[2]-1.0*f[1]+f[0])); 
  if(wv > 0) {
  xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
  fReflXYMuQuad[0][0] = 0.25*(1.732050807568877*(f[27]-1.0*f[22]+f[21]-1.0*(f[16]+f[14])+f[8]-1.0*f[7]+f[3])*zVal+f[20]-1.0*f[13]+f[12]-1.0*(f[6]+f[5])+f[2]-1.0*f[1]+f[0])*(exp(b*xc)-exp(-b))/sinh(b); 
  fReflXYMuQuad[0][1] = 0.25*(1.732050807568877*(f[31]-1.0*f[30]+f[29]-1.0*(f[26]+f[25])+f[19]-1.0*f[18]+f[11])*zVal+f[28]-1.0*f[24]+f[23]-1.0*(f[17]+f[15])+f[10]-1.0*f[9]+f[4])*((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/(b*cosh(b)-sinh(b)); 
  } else { 
  xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
  fReflXYMuQuad[0][0] = 0.25*(1.732050807568877*(f[27]-1.0*f[22]+f[21]-1.0*(f[16]+f[14])+f[8]-1.0*f[7]+f[3])*zVal+f[20]-1.0*f[13]+f[12]-1.0*(f[6]+f[5])+f[2]-1.0*f[1]+f[0])*(exp(b)-exp(b*xc))/sinh(b); 
  fReflXYMuQuad[0][1] = 0.25*(1.732050807568877*(f[31]-1.0*f[30]+f[29]-1.0*(f[26]+f[25])+f[19]-1.0*f[18]+f[11])*zVal+f[28]-1.0*f[24]+f[23]-1.0*(f[17]+f[15])+f[10]-1.0*f[9]+f[4])*((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/(b*cosh(b)-sinh(b)); 
  } 
  b = invL((0.5773502691896258*(1.732050807568877*(f[31]-1.0*f[30]+f[29]+f[26]-1.0*(f[25]+f[19])+f[18]-1.0*f[11])*zVal+f[28]-1.0*f[24]+f[23]+f[17]-1.0*(f[15]+f[10])+f[9]-1.0*f[4]))/(1.732050807568877*(f[27]-1.0*f[22]+f[21]+f[16]-1.0*(f[14]+f[8])+f[7]-1.0*f[3])*zVal+f[20]-1.0*f[13]+f[12]+f[6]-1.0*(f[5]+f[2])+f[1]-1.0*f[0])); 
  if(wv > 0) {
  xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
  fReflXYMuQuad[1][0] = -0.25*(1.732050807568877*(f[27]-1.0*f[22]+f[21]+f[16]-1.0*(f[14]+f[8])+f[7]-1.0*f[3])*zVal+f[20]-1.0*f[13]+f[12]+f[6]-1.0*(f[5]+f[2])+f[1]-1.0*f[0])*(exp(b*xc)-exp(-b))/sinh(b); 
  fReflXYMuQuad[1][1] = -0.25*(1.732050807568877*(f[31]-1.0*f[30]+f[29]+f[26]-1.0*(f[25]+f[19])+f[18]-1.0*f[11])*zVal+f[28]-1.0*f[24]+f[23]+f[17]-1.0*(f[15]+f[10])+f[9]-1.0*f[4])*((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/(b*cosh(b)-sinh(b)); 
  } else { 
  xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
  fReflXYMuQuad[1][0] = -0.25*(1.732050807568877*(f[27]-1.0*f[22]+f[21]+f[16]-1.0*(f[14]+f[8])+f[7]-1.0*f[3])*zVal+f[20]-1.0*f[13]+f[12]+f[6]-1.0*(f[5]+f[2])+f[1]-1.0*f[0])*(exp(b)-exp(b*xc))/sinh(b); 
  fReflXYMuQuad[1][1] = -0.25*(1.732050807568877*(f[31]-1.0*f[30]+f[29]+f[26]-1.0*(f[25]+f[19])+f[18]-1.0*f[11])*zVal+f[28]-1.0*f[24]+f[23]+f[17]-1.0*(f[15]+f[10])+f[9]-1.0*f[4])*((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/(b*cosh(b)-sinh(b)); 
  } 
  fReflXYQuad[2][0] = 0.7071067811865468*(fReflXYMuQuad[1][0]+fReflXYMuQuad[0][0]); 
  fReflXYQuad[2][1] = 0.7071067811865468*(fReflXYMuQuad[1][1]+fReflXYMuQuad[0][1]); 
  fReflXYQuad[2][2] = 0.7071067811865468*(fReflXYMuQuad[1][0]-1.0*fReflXYMuQuad[0][0]); 
  fReflXYQuad[2][3] = 0.7071067811865468*(fReflXYMuQuad[1][1]-1.0*fReflXYMuQuad[0][1]); 
  } 

 
// quadrature node (x,y)_i=4 
  vcutSq_i = (0.5*q_*((2.449489742783178*phiWall[7]-2.449489742783178*phi[7]+2.449489742783178*phiWall[6]-2.449489742783178*phi[6]+2.449489742783178*phiWall[5]-2.449489742783178*phi[5]+2.449489742783178*phiWall[3]-2.449489742783178*phi[3])*zVal+1.414213562373095*phiWall[4]-1.414213562373095*phi[4]+1.414213562373095*phiWall[2]-1.414213562373095*phi[2]+1.414213562373095*phiWall[1]-1.414213562373095*phi[1]+1.414213562373095*phiWall[0]-1.414213562373095*phi[0]))/m_; 
  if(vcutSq_i <= vlowerSq) { // absorb (no reflection) 
  fReflXYQuad[3][0] = 0.0; 
  fReflXYQuad[3][1] = 0.0; 
  fReflXYQuad[3][2] = 0.0; 
  fReflXYQuad[3][3] = 0.0; 
  } else if(vcutSq_i > vupperSq) { // full reflection 
  fReflXYQuad[3][0] = 0.3535533905932737*(1.732050807568877*(f[16]+f[8]+f[7]+f[3])*zVal+f[6]+f[2]+f[1]+f[0]); 
  fReflXYQuad[3][1] = 0.3535533905932737*(1.732050807568877*(f[26]+f[19]+f[18]+f[11])*zVal+f[17]+f[10]+f[9]+f[4]); 
  fReflXYQuad[3][2] = 0.3535533905932737*(1.732050807568877*(f[27]+f[22]+f[21]+f[14])*zVal+f[20]+f[13]+f[12]+f[5]); 
  fReflXYQuad[3][3] = 0.3535533905932737*(1.732050807568877*(f[31]+f[30]+f[29]+f[25])*zVal+f[28]+f[24]+f[23]+f[15]); 
  } else { // partial reflection 
  b = invL((0.5773502691896258*(1.732050807568877*(f[31]+f[30]+f[29]-1.0*f[26]+f[25]-1.0*(f[19]+f[18]+f[11]))*zVal+f[28]+f[24]+f[23]-1.0*f[17]+f[15]-1.0*(f[10]+f[9]+f[4])))/(1.732050807568877*(f[27]+f[22]+f[21]-1.0*f[16]+f[14]-1.0*(f[8]+f[7]+f[3]))*zVal+f[20]+f[13]+f[12]-1.0*f[6]+f[5]-1.0*(f[2]+f[1]+f[0]))); 
  if(wv > 0) {
  xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
  fReflXYMuQuad[0][0] = -0.25*(1.732050807568877*(f[27]+f[22]+f[21]-1.0*f[16]+f[14]-1.0*(f[8]+f[7]+f[3]))*zVal+f[20]+f[13]+f[12]-1.0*f[6]+f[5]-1.0*(f[2]+f[1]+f[0]))*(exp(b*xc)-exp(-b))/sinh(b); 
  fReflXYMuQuad[0][1] = -0.25*(1.732050807568877*(f[31]+f[30]+f[29]-1.0*f[26]+f[25]-1.0*(f[19]+f[18]+f[11]))*zVal+f[28]+f[24]+f[23]-1.0*f[17]+f[15]-1.0*(f[10]+f[9]+f[4]))*((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/(b*cosh(b)-sinh(b)); 
  } else { 
  xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
  fReflXYMuQuad[0][0] = -0.25*(1.732050807568877*(f[27]+f[22]+f[21]-1.0*f[16]+f[14]-1.0*(f[8]+f[7]+f[3]))*zVal+f[20]+f[13]+f[12]-1.0*f[6]+f[5]-1.0*(f[2]+f[1]+f[0]))*(exp(b)-exp(b*xc))/sinh(b); 
  fReflXYMuQuad[0][1] = -0.25*(1.732050807568877*(f[31]+f[30]+f[29]-1.0*f[26]+f[25]-1.0*(f[19]+f[18]+f[11]))*zVal+f[28]+f[24]+f[23]-1.0*f[17]+f[15]-1.0*(f[10]+f[9]+f[4]))*((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/(b*cosh(b)-sinh(b)); 
  } 
  b = invL((0.5773502691896258*(1.732050807568877*(f[31]+f[30]+f[29]+f[26]+f[25]+f[19]+f[18]+f[11])*zVal+f[28]+f[24]+f[23]+f[17]+f[15]+f[10]+f[9]+f[4]))/(1.732050807568877*(f[27]+f[22]+f[21]+f[16]+f[14]+f[8]+f[7]+f[3])*zVal+f[20]+f[13]+f[12]+f[6]+f[5]+f[2]+f[1]+f[0])); 
  if(wv > 0) {
  xc = (2*(sqrt(vcutSq_i)-wv))/dv; 
  fReflXYMuQuad[1][0] = 0.25*(1.732050807568877*(f[27]+f[22]+f[21]+f[16]+f[14]+f[8]+f[7]+f[3])*zVal+f[20]+f[13]+f[12]+f[6]+f[5]+f[2]+f[1]+f[0])*(exp(b*xc)-exp(-b))/sinh(b); 
  fReflXYMuQuad[1][1] = 0.25*(1.732050807568877*(f[31]+f[30]+f[29]+f[26]+f[25]+f[19]+f[18]+f[11])*zVal+f[28]+f[24]+f[23]+f[17]+f[15]+f[10]+f[9]+f[4])*((b*xc-1)*exp(b*xc)+(b+1)*exp(-b))/(b*cosh(b)-sinh(b)); 
  } else { 
  xc = (2*((-wv)-sqrt(vcutSq_i)))/dv; 
  fReflXYMuQuad[1][0] = 0.25*(1.732050807568877*(f[27]+f[22]+f[21]+f[16]+f[14]+f[8]+f[7]+f[3])*zVal+f[20]+f[13]+f[12]+f[6]+f[5]+f[2]+f[1]+f[0])*(exp(b)-exp(b*xc))/sinh(b); 
  fReflXYMuQuad[1][1] = 0.25*(1.732050807568877*(f[31]+f[30]+f[29]+f[26]+f[25]+f[19]+f[18]+f[11])*zVal+f[28]+f[24]+f[23]+f[17]+f[15]+f[10]+f[9]+f[4])*((b-1)*exp(b)-(b*xc-1)*exp(b*xc))/(b*cosh(b)-sinh(b)); 
  } 
  fReflXYQuad[3][0] = 0.7071067811865468*(fReflXYMuQuad[1][0]+fReflXYMuQuad[0][0]); 
  fReflXYQuad[3][1] = 0.7071067811865468*(fReflXYMuQuad[1][1]+fReflXYMuQuad[0][1]); 
  fReflXYQuad[3][2] = 0.7071067811865468*(fReflXYMuQuad[1][0]-1.0*fReflXYMuQuad[0][0]); 
  fReflXYQuad[3][3] = 0.7071067811865468*(fReflXYMuQuad[1][1]-1.0*fReflXYMuQuad[0][1]); 
  } 

 
  fRefl[0] = 0.7071067811865475*(fReflXYQuad[3][0]+fReflXYQuad[2][0]+fReflXYQuad[1][0]+fReflXYQuad[0][0]); 
  fRefl[1] = 0.7071067811865475*(fReflXYQuad[3][0]-1.0*fReflXYQuad[2][0]+fReflXYQuad[1][0]-1.0*fReflXYQuad[0][0]); 
  fRefl[2] = 0.7071067811865475*(fReflXYQuad[3][0]+fReflXYQuad[2][0]-1.0*fReflXYQuad[1][0]-1.0*fReflXYQuad[0][0]); 
  fRefl[3] = 0.0; 
  fRefl[4] = 0.7071067811865475*(fReflXYQuad[3][1]+fReflXYQuad[2][1]+fReflXYQuad[1][1]+fReflXYQuad[0][1]); 
  fRefl[5] = 0.7071067811865475*(fReflXYQuad[3][2]+fReflXYQuad[2][2]+fReflXYQuad[1][2]+fReflXYQuad[0][2]); 
  fRefl[6] = 0.7071067811865475*(fReflXYQuad[3][0]-1.0*fReflXYQuad[2][0]-1.0*fReflXYQuad[1][0]+fReflXYQuad[0][0]); 
  fRefl[7] = 0.0; 
  fRefl[8] = 0.0; 
  fRefl[9] = 0.7071067811865475*(fReflXYQuad[3][1]-1.0*fReflXYQuad[2][1]+fReflXYQuad[1][1]-1.0*fReflXYQuad[0][1]); 
  fRefl[10] = 0.7071067811865475*(fReflXYQuad[3][1]+fReflXYQuad[2][1]-1.0*fReflXYQuad[1][1]-1.0*fReflXYQuad[0][1]); 
  fRefl[11] = 0.0; 
  fRefl[12] = 0.7071067811865475*(fReflXYQuad[3][2]-1.0*fReflXYQuad[2][2]+fReflXYQuad[1][2]-1.0*fReflXYQuad[0][2]); 
  fRefl[13] = 0.7071067811865475*(fReflXYQuad[3][2]+fReflXYQuad[2][2]-1.0*fReflXYQuad[1][2]-1.0*fReflXYQuad[0][2]); 
  fRefl[14] = 0.0; 
  fRefl[15] = 0.7071067811865475*(fReflXYQuad[3][3]+fReflXYQuad[2][3]+fReflXYQuad[1][3]+fReflXYQuad[0][3]); 
  fRefl[16] = 0.0; 
  fRefl[17] = 0.7071067811865475*(fReflXYQuad[3][1]-1.0*fReflXYQuad[2][1]-1.0*fReflXYQuad[1][1]+fReflXYQuad[0][1]); 
  fRefl[18] = 0.0; 
  fRefl[19] = 0.0; 
  fRefl[20] = 0.7071067811865475*(fReflXYQuad[3][2]-1.0*fReflXYQuad[2][2]-1.0*fReflXYQuad[1][2]+fReflXYQuad[0][2]); 
  fRefl[21] = 0.0; 
  fRefl[22] = 0.0; 
  fRefl[23] = 0.7071067811865475*(fReflXYQuad[3][3]-1.0*fReflXYQuad[2][3]+fReflXYQuad[1][3]-1.0*fReflXYQuad[0][3]); 
  fRefl[24] = 0.7071067811865475*(fReflXYQuad[3][3]+fReflXYQuad[2][3]-1.0*fReflXYQuad[1][3]-1.0*fReflXYQuad[0][3]); 
  fRefl[25] = 0.0; 
  fRefl[26] = 0.0; 
  fRefl[27] = 0.0; 
  fRefl[28] = 0.7071067811865475*(fReflXYQuad[3][3]-1.0*fReflXYQuad[2][3]-1.0*fReflXYQuad[1][3]+fReflXYQuad[0][3]); 
  fRefl[29] = 0.0; 
  fRefl[30] = 0.0; 
  fRefl[31] = 0.0; 
}
