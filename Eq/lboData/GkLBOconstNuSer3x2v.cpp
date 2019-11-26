#include <GkLBOModDecl.h> 
double GkLBOconstNuVol3x2vSerP1(const double m_, const double *w, const double *dxv, double *cflRateByDir, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[5]:            Cell-center coordinates. 
  // dxv[5]:          Cell spacing. 
  // cflRatebyDir[5]: CFL rate in each direction. 
  // nuSum:           collisionalities added (self and cross species collisionalities). 
  // nuUSum:          sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum:       sum of thermal speeds squared time their respective collisionalities. 
  // f[32]:            Input distribution function. 
  // out[32]:          Incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0]   = 2.0/dxv[3]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[0]   = rdv2[0]; 
  rdv2[1]   = 2.0/dxv[4]; 
  rdvSq4[1] = rdv2[1]*rdv2[1]; 
  rdv2[1]   = rdv2[1]; 

  double alphaVpar[32]; 
  alphaVpar[0] = rdv2[0]*(2.0*nuUSum[0]-5.656854249492382*w[3]*nuSum); 
  alphaVpar[1] = 2.0*rdv2[0]*nuUSum[1]; 
  alphaVpar[2] = 2.0*rdv2[0]*nuUSum[2]; 
  alphaVpar[3] = 2.0*rdv2[0]*nuUSum[3]; 
  alphaVpar[4] = -3.265986323710906*nuSum; 
  alphaVpar[6] = 2.0*rdv2[0]*nuUSum[4]; 
  alphaVpar[7] = 2.0*rdv2[0]*nuUSum[5]; 
  alphaVpar[8] = 2.0*rdv2[0]*nuUSum[6]; 
  alphaVpar[16] = 2.0*rdv2[0]*nuUSum[7]; 

  double alphaMu[32]; 
  alphaMu[0] = (-11.31370849898477*rdv2[1]*w[4]*nuSum)+1.414213562373095*rdv2[1]*BmagInv[7]*nuVtSqSum[7]*m_+1.414213562373095*rdv2[1]*BmagInv[6]*nuVtSqSum[6]*m_+1.414213562373095*rdv2[1]*BmagInv[5]*nuVtSqSum[5]*m_+1.414213562373095*rdv2[1]*BmagInv[4]*nuVtSqSum[4]*m_+1.414213562373095*rdv2[1]*BmagInv[3]*nuVtSqSum[3]*m_+1.414213562373095*rdv2[1]*BmagInv[2]*nuVtSqSum[2]*m_+1.414213562373095*BmagInv[1]*nuVtSqSum[1]*rdv2[1]*m_+1.414213562373095*BmagInv[0]*nuVtSqSum[0]*rdv2[1]*m_; 
  alphaMu[1] = 1.414213562373095*rdv2[1]*BmagInv[6]*nuVtSqSum[7]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[6]*BmagInv[7]*m_+1.414213562373095*rdv2[1]*BmagInv[3]*nuVtSqSum[5]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[3]*BmagInv[5]*m_+1.414213562373095*rdv2[1]*BmagInv[2]*nuVtSqSum[4]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[2]*BmagInv[4]*m_+1.414213562373095*BmagInv[0]*nuVtSqSum[1]*rdv2[1]*m_+1.414213562373095*nuVtSqSum[0]*BmagInv[1]*rdv2[1]*m_; 
  alphaMu[2] = 1.414213562373095*rdv2[1]*BmagInv[5]*nuVtSqSum[7]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[5]*BmagInv[7]*m_+1.414213562373095*rdv2[1]*BmagInv[3]*nuVtSqSum[6]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[3]*BmagInv[6]*m_+1.414213562373095*BmagInv[1]*rdv2[1]*nuVtSqSum[4]*m_+1.414213562373095*nuVtSqSum[1]*rdv2[1]*BmagInv[4]*m_+1.414213562373095*BmagInv[0]*rdv2[1]*nuVtSqSum[2]*m_+1.414213562373095*nuVtSqSum[0]*rdv2[1]*BmagInv[2]*m_; 
  alphaMu[3] = 1.414213562373095*rdv2[1]*BmagInv[4]*nuVtSqSum[7]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[4]*BmagInv[7]*m_+1.414213562373095*rdv2[1]*BmagInv[2]*nuVtSqSum[6]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[2]*BmagInv[6]*m_+1.414213562373095*BmagInv[1]*rdv2[1]*nuVtSqSum[5]*m_+1.414213562373095*nuVtSqSum[1]*rdv2[1]*BmagInv[5]*m_+1.414213562373095*BmagInv[0]*rdv2[1]*nuVtSqSum[3]*m_+1.414213562373095*nuVtSqSum[0]*rdv2[1]*BmagInv[3]*m_; 
  alphaMu[5] = -6.531972647421813*nuSum; 
  alphaMu[6] = 1.414213562373095*rdv2[1]*BmagInv[3]*nuVtSqSum[7]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[3]*BmagInv[7]*m_+1.414213562373095*rdv2[1]*BmagInv[5]*nuVtSqSum[6]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[5]*BmagInv[6]*m_+1.414213562373095*BmagInv[0]*rdv2[1]*nuVtSqSum[4]*m_+1.414213562373095*nuVtSqSum[0]*rdv2[1]*BmagInv[4]*m_+1.414213562373095*BmagInv[1]*rdv2[1]*nuVtSqSum[2]*m_+1.414213562373095*nuVtSqSum[1]*rdv2[1]*BmagInv[2]*m_; 
  alphaMu[7] = 1.414213562373095*rdv2[1]*BmagInv[2]*nuVtSqSum[7]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[2]*BmagInv[7]*m_+1.414213562373095*rdv2[1]*BmagInv[4]*nuVtSqSum[6]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[4]*BmagInv[6]*m_+1.414213562373095*BmagInv[0]*rdv2[1]*nuVtSqSum[5]*m_+1.414213562373095*nuVtSqSum[0]*rdv2[1]*BmagInv[5]*m_+1.414213562373095*BmagInv[1]*rdv2[1]*nuVtSqSum[3]*m_+1.414213562373095*nuVtSqSum[1]*rdv2[1]*BmagInv[3]*m_; 
  alphaMu[8] = 1.414213562373095*BmagInv[1]*rdv2[1]*nuVtSqSum[7]*m_+1.414213562373095*nuVtSqSum[1]*rdv2[1]*BmagInv[7]*m_+1.414213562373095*BmagInv[0]*rdv2[1]*nuVtSqSum[6]*m_+1.414213562373095*nuVtSqSum[0]*rdv2[1]*BmagInv[6]*m_+1.414213562373095*rdv2[1]*BmagInv[4]*nuVtSqSum[5]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[4]*BmagInv[5]*m_+1.414213562373095*rdv2[1]*BmagInv[2]*nuVtSqSum[3]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[2]*BmagInv[3]*m_; 
  alphaMu[16] = 1.414213562373095*BmagInv[0]*rdv2[1]*nuVtSqSum[7]*m_+1.414213562373095*nuVtSqSum[0]*rdv2[1]*BmagInv[7]*m_+1.414213562373095*BmagInv[1]*rdv2[1]*nuVtSqSum[6]*m_+1.414213562373095*nuVtSqSum[1]*rdv2[1]*BmagInv[6]*m_+1.414213562373095*rdv2[1]*BmagInv[2]*nuVtSqSum[5]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[2]*BmagInv[5]*m_+1.414213562373095*rdv2[1]*BmagInv[3]*nuVtSqSum[4]*m_+1.414213562373095*rdv2[1]*nuVtSqSum[3]*BmagInv[4]*m_; 

  cflRateByDir[0] = 0.; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points. 
  alphaL = 0.03125*((-0.1767766952966368*alphaVpar[16])+0.1767766952966368*(alphaVpar[8]+alphaVpar[7]+alphaVpar[6])-0.3061862178478971*alphaVpar[4]-0.1767766952966368*(alphaVpar[3]+alphaVpar[2]+alphaVpar[1])+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*(alphaVpar[16]+alphaVpar[8])-0.1767766952966368*(alphaVpar[7]+alphaVpar[6])-0.3061862178478971*alphaVpar[4]-0.1767766952966368*(alphaVpar[3]+alphaVpar[2])+0.1767766952966368*(alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*alphaVpar[16]-0.1767766952966368*alphaVpar[8]+0.1767766952966368*alphaVpar[7]-0.1767766952966368*alphaVpar[6]-0.3061862178478971*alphaVpar[4]-0.1767766952966368*alphaVpar[3]+0.1767766952966368*alphaVpar[2]-0.1767766952966368*alphaVpar[1]+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*(alphaVpar[16]+alphaVpar[8]+alphaVpar[7]))+0.1767766952966368*alphaVpar[6]-0.3061862178478971*alphaVpar[4]-0.1767766952966368*alphaVpar[3]+0.1767766952966368*(alphaVpar[2]+alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*alphaVpar[16]-0.1767766952966368*(alphaVpar[8]+alphaVpar[7])+0.1767766952966368*alphaVpar[6]-0.3061862178478971*alphaVpar[4]+0.1767766952966368*alphaVpar[3]-0.1767766952966368*(alphaVpar[2]+alphaVpar[1])+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*(alphaVpar[16]+alphaVpar[8]))+0.1767766952966368*alphaVpar[7]-0.1767766952966368*alphaVpar[6]-0.3061862178478971*alphaVpar[4]+0.1767766952966368*alphaVpar[3]-0.1767766952966368*alphaVpar[2]+0.1767766952966368*(alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*alphaVpar[16])+0.1767766952966368*alphaVpar[8]-0.1767766952966368*(alphaVpar[7]+alphaVpar[6])-0.3061862178478971*alphaVpar[4]+0.1767766952966368*(alphaVpar[3]+alphaVpar[2])-0.1767766952966368*alphaVpar[1]+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*(alphaVpar[16]+alphaVpar[8]+alphaVpar[7]+alphaVpar[6])-0.3061862178478971*alphaVpar[4]+0.1767766952966368*(alphaVpar[3]+alphaVpar[2]+alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*alphaVpar[16])+0.1767766952966368*(alphaVpar[8]+alphaVpar[7]+alphaVpar[6])-0.3061862178478971*alphaVpar[4]-0.1767766952966368*(alphaVpar[3]+alphaVpar[2]+alphaVpar[1])+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*(alphaVpar[16]+alphaVpar[8])-0.1767766952966368*(alphaVpar[7]+alphaVpar[6])-0.3061862178478971*alphaVpar[4]-0.1767766952966368*(alphaVpar[3]+alphaVpar[2])+0.1767766952966368*(alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*alphaVpar[16]-0.1767766952966368*alphaVpar[8]+0.1767766952966368*alphaVpar[7]-0.1767766952966368*alphaVpar[6]-0.3061862178478971*alphaVpar[4]-0.1767766952966368*alphaVpar[3]+0.1767766952966368*alphaVpar[2]-0.1767766952966368*alphaVpar[1]+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*(alphaVpar[16]+alphaVpar[8]+alphaVpar[7]))+0.1767766952966368*alphaVpar[6]-0.3061862178478971*alphaVpar[4]-0.1767766952966368*alphaVpar[3]+0.1767766952966368*(alphaVpar[2]+alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*alphaVpar[16]-0.1767766952966368*(alphaVpar[8]+alphaVpar[7])+0.1767766952966368*alphaVpar[6]-0.3061862178478971*alphaVpar[4]+0.1767766952966368*alphaVpar[3]-0.1767766952966368*(alphaVpar[2]+alphaVpar[1])+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*(alphaVpar[16]+alphaVpar[8]))+0.1767766952966368*alphaVpar[7]-0.1767766952966368*alphaVpar[6]-0.3061862178478971*alphaVpar[4]+0.1767766952966368*alphaVpar[3]-0.1767766952966368*alphaVpar[2]+0.1767766952966368*(alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*alphaVpar[16])+0.1767766952966368*alphaVpar[8]-0.1767766952966368*(alphaVpar[7]+alphaVpar[6])-0.3061862178478971*alphaVpar[4]+0.1767766952966368*(alphaVpar[3]+alphaVpar[2])-0.1767766952966368*alphaVpar[1]+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*(alphaVpar[16]+alphaVpar[8]+alphaVpar[7]+alphaVpar[6])-0.3061862178478971*alphaVpar[4]+0.1767766952966368*(alphaVpar[3]+alphaVpar[2]+alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaL); 
  // Evaluate alpha at right surface quadrature points. 
  alphaR = 0.03125*((-0.1767766952966368*alphaVpar[16])+0.1767766952966368*(alphaVpar[8]+alphaVpar[7]+alphaVpar[6])+0.3061862178478971*alphaVpar[4]-0.1767766952966368*(alphaVpar[3]+alphaVpar[2]+alphaVpar[1])+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*(alphaVpar[16]+alphaVpar[8])-0.1767766952966368*(alphaVpar[7]+alphaVpar[6])+0.3061862178478971*alphaVpar[4]-0.1767766952966368*(alphaVpar[3]+alphaVpar[2])+0.1767766952966368*(alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*alphaVpar[16]-0.1767766952966368*alphaVpar[8]+0.1767766952966368*alphaVpar[7]-0.1767766952966368*alphaVpar[6]+0.3061862178478971*alphaVpar[4]-0.1767766952966368*alphaVpar[3]+0.1767766952966368*alphaVpar[2]-0.1767766952966368*alphaVpar[1]+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*(alphaVpar[16]+alphaVpar[8]+alphaVpar[7]))+0.1767766952966368*alphaVpar[6]+0.3061862178478971*alphaVpar[4]-0.1767766952966368*alphaVpar[3]+0.1767766952966368*(alphaVpar[2]+alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*alphaVpar[16]-0.1767766952966368*(alphaVpar[8]+alphaVpar[7])+0.1767766952966368*alphaVpar[6]+0.3061862178478971*alphaVpar[4]+0.1767766952966368*alphaVpar[3]-0.1767766952966368*(alphaVpar[2]+alphaVpar[1])+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*(alphaVpar[16]+alphaVpar[8]))+0.1767766952966368*alphaVpar[7]-0.1767766952966368*alphaVpar[6]+0.3061862178478971*alphaVpar[4]+0.1767766952966368*alphaVpar[3]-0.1767766952966368*alphaVpar[2]+0.1767766952966368*(alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*alphaVpar[16])+0.1767766952966368*alphaVpar[8]-0.1767766952966368*(alphaVpar[7]+alphaVpar[6])+0.3061862178478971*alphaVpar[4]+0.1767766952966368*(alphaVpar[3]+alphaVpar[2])-0.1767766952966368*alphaVpar[1]+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*(alphaVpar[16]+alphaVpar[8]+alphaVpar[7]+alphaVpar[6])+0.3061862178478971*alphaVpar[4]+0.1767766952966368*(alphaVpar[3]+alphaVpar[2]+alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*alphaVpar[16])+0.1767766952966368*(alphaVpar[8]+alphaVpar[7]+alphaVpar[6])+0.3061862178478971*alphaVpar[4]-0.1767766952966368*(alphaVpar[3]+alphaVpar[2]+alphaVpar[1])+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*(alphaVpar[16]+alphaVpar[8])-0.1767766952966368*(alphaVpar[7]+alphaVpar[6])+0.3061862178478971*alphaVpar[4]-0.1767766952966368*(alphaVpar[3]+alphaVpar[2])+0.1767766952966368*(alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*alphaVpar[16]-0.1767766952966368*alphaVpar[8]+0.1767766952966368*alphaVpar[7]-0.1767766952966368*alphaVpar[6]+0.3061862178478971*alphaVpar[4]-0.1767766952966368*alphaVpar[3]+0.1767766952966368*alphaVpar[2]-0.1767766952966368*alphaVpar[1]+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*(alphaVpar[16]+alphaVpar[8]+alphaVpar[7]))+0.1767766952966368*alphaVpar[6]+0.3061862178478971*alphaVpar[4]-0.1767766952966368*alphaVpar[3]+0.1767766952966368*(alphaVpar[2]+alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*alphaVpar[16]-0.1767766952966368*(alphaVpar[8]+alphaVpar[7])+0.1767766952966368*alphaVpar[6]+0.3061862178478971*alphaVpar[4]+0.1767766952966368*alphaVpar[3]-0.1767766952966368*(alphaVpar[2]+alphaVpar[1])+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*(alphaVpar[16]+alphaVpar[8]))+0.1767766952966368*alphaVpar[7]-0.1767766952966368*alphaVpar[6]+0.3061862178478971*alphaVpar[4]+0.1767766952966368*alphaVpar[3]-0.1767766952966368*alphaVpar[2]+0.1767766952966368*(alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*alphaVpar[16])+0.1767766952966368*alphaVpar[8]-0.1767766952966368*(alphaVpar[7]+alphaVpar[6])+0.3061862178478971*alphaVpar[4]+0.1767766952966368*(alphaVpar[3]+alphaVpar[2])-0.1767766952966368*alphaVpar[1]+0.1767766952966368*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*(alphaVpar[16]+alphaVpar[8]+alphaVpar[7]+alphaVpar[6])+0.3061862178478971*alphaVpar[4]+0.1767766952966368*(alphaVpar[3]+alphaVpar[2]+alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaR); 
  // Evaluate alpha at left surface quadrature points. 
  alphaL = 0.03125*((-0.1767766952966368*alphaMu[16])+0.1767766952966368*(alphaMu[8]+alphaMu[7]+alphaMu[6])-0.3061862178478971*alphaMu[5]-0.1767766952966368*(alphaMu[3]+alphaMu[2]+alphaMu[1])+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*(alphaMu[16]+alphaMu[8])-0.1767766952966368*(alphaMu[7]+alphaMu[6])-0.3061862178478971*alphaMu[5]-0.1767766952966368*(alphaMu[3]+alphaMu[2])+0.1767766952966368*(alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*alphaMu[16]-0.1767766952966368*alphaMu[8]+0.1767766952966368*alphaMu[7]-0.1767766952966368*alphaMu[6]-0.3061862178478971*alphaMu[5]-0.1767766952966368*alphaMu[3]+0.1767766952966368*alphaMu[2]-0.1767766952966368*alphaMu[1]+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*(alphaMu[16]+alphaMu[8]+alphaMu[7]))+0.1767766952966368*alphaMu[6]-0.3061862178478971*alphaMu[5]-0.1767766952966368*alphaMu[3]+0.1767766952966368*(alphaMu[2]+alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*alphaMu[16]-0.1767766952966368*(alphaMu[8]+alphaMu[7])+0.1767766952966368*alphaMu[6]-0.3061862178478971*alphaMu[5]+0.1767766952966368*alphaMu[3]-0.1767766952966368*(alphaMu[2]+alphaMu[1])+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*(alphaMu[16]+alphaMu[8]))+0.1767766952966368*alphaMu[7]-0.1767766952966368*alphaMu[6]-0.3061862178478971*alphaMu[5]+0.1767766952966368*alphaMu[3]-0.1767766952966368*alphaMu[2]+0.1767766952966368*(alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*alphaMu[16])+0.1767766952966368*alphaMu[8]-0.1767766952966368*(alphaMu[7]+alphaMu[6])-0.3061862178478971*alphaMu[5]+0.1767766952966368*(alphaMu[3]+alphaMu[2])-0.1767766952966368*alphaMu[1]+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*(alphaMu[16]+alphaMu[8]+alphaMu[7]+alphaMu[6])-0.3061862178478971*alphaMu[5]+0.1767766952966368*(alphaMu[3]+alphaMu[2]+alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*alphaMu[16])+0.1767766952966368*(alphaMu[8]+alphaMu[7]+alphaMu[6])-0.3061862178478971*alphaMu[5]-0.1767766952966368*(alphaMu[3]+alphaMu[2]+alphaMu[1])+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*(alphaMu[16]+alphaMu[8])-0.1767766952966368*(alphaMu[7]+alphaMu[6])-0.3061862178478971*alphaMu[5]-0.1767766952966368*(alphaMu[3]+alphaMu[2])+0.1767766952966368*(alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*alphaMu[16]-0.1767766952966368*alphaMu[8]+0.1767766952966368*alphaMu[7]-0.1767766952966368*alphaMu[6]-0.3061862178478971*alphaMu[5]-0.1767766952966368*alphaMu[3]+0.1767766952966368*alphaMu[2]-0.1767766952966368*alphaMu[1]+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*(alphaMu[16]+alphaMu[8]+alphaMu[7]))+0.1767766952966368*alphaMu[6]-0.3061862178478971*alphaMu[5]-0.1767766952966368*alphaMu[3]+0.1767766952966368*(alphaMu[2]+alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*alphaMu[16]-0.1767766952966368*(alphaMu[8]+alphaMu[7])+0.1767766952966368*alphaMu[6]-0.3061862178478971*alphaMu[5]+0.1767766952966368*alphaMu[3]-0.1767766952966368*(alphaMu[2]+alphaMu[1])+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*(alphaMu[16]+alphaMu[8]))+0.1767766952966368*alphaMu[7]-0.1767766952966368*alphaMu[6]-0.3061862178478971*alphaMu[5]+0.1767766952966368*alphaMu[3]-0.1767766952966368*alphaMu[2]+0.1767766952966368*(alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*((-0.1767766952966368*alphaMu[16])+0.1767766952966368*alphaMu[8]-0.1767766952966368*(alphaMu[7]+alphaMu[6])-0.3061862178478971*alphaMu[5]+0.1767766952966368*(alphaMu[3]+alphaMu[2])-0.1767766952966368*alphaMu[1]+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.03125*(0.1767766952966368*(alphaMu[16]+alphaMu[8]+alphaMu[7]+alphaMu[6])-0.3061862178478971*alphaMu[5]+0.1767766952966368*(alphaMu[3]+alphaMu[2]+alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  // Evaluate alpha at right surface quadrature points. 
  alphaR = 0.03125*((-0.1767766952966368*alphaMu[16])+0.1767766952966368*(alphaMu[8]+alphaMu[7]+alphaMu[6])+0.3061862178478971*alphaMu[5]-0.1767766952966368*(alphaMu[3]+alphaMu[2]+alphaMu[1])+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*(alphaMu[16]+alphaMu[8])-0.1767766952966368*(alphaMu[7]+alphaMu[6])+0.3061862178478971*alphaMu[5]-0.1767766952966368*(alphaMu[3]+alphaMu[2])+0.1767766952966368*(alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*alphaMu[16]-0.1767766952966368*alphaMu[8]+0.1767766952966368*alphaMu[7]-0.1767766952966368*alphaMu[6]+0.3061862178478971*alphaMu[5]-0.1767766952966368*alphaMu[3]+0.1767766952966368*alphaMu[2]-0.1767766952966368*alphaMu[1]+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*(alphaMu[16]+alphaMu[8]+alphaMu[7]))+0.1767766952966368*alphaMu[6]+0.3061862178478971*alphaMu[5]-0.1767766952966368*alphaMu[3]+0.1767766952966368*(alphaMu[2]+alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*alphaMu[16]-0.1767766952966368*(alphaMu[8]+alphaMu[7])+0.1767766952966368*alphaMu[6]+0.3061862178478971*alphaMu[5]+0.1767766952966368*alphaMu[3]-0.1767766952966368*(alphaMu[2]+alphaMu[1])+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*(alphaMu[16]+alphaMu[8]))+0.1767766952966368*alphaMu[7]-0.1767766952966368*alphaMu[6]+0.3061862178478971*alphaMu[5]+0.1767766952966368*alphaMu[3]-0.1767766952966368*alphaMu[2]+0.1767766952966368*(alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*alphaMu[16])+0.1767766952966368*alphaMu[8]-0.1767766952966368*(alphaMu[7]+alphaMu[6])+0.3061862178478971*alphaMu[5]+0.1767766952966368*(alphaMu[3]+alphaMu[2])-0.1767766952966368*alphaMu[1]+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*(alphaMu[16]+alphaMu[8]+alphaMu[7]+alphaMu[6])+0.3061862178478971*alphaMu[5]+0.1767766952966368*(alphaMu[3]+alphaMu[2]+alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*alphaMu[16])+0.1767766952966368*(alphaMu[8]+alphaMu[7]+alphaMu[6])+0.3061862178478971*alphaMu[5]-0.1767766952966368*(alphaMu[3]+alphaMu[2]+alphaMu[1])+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*(alphaMu[16]+alphaMu[8])-0.1767766952966368*(alphaMu[7]+alphaMu[6])+0.3061862178478971*alphaMu[5]-0.1767766952966368*(alphaMu[3]+alphaMu[2])+0.1767766952966368*(alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*alphaMu[16]-0.1767766952966368*alphaMu[8]+0.1767766952966368*alphaMu[7]-0.1767766952966368*alphaMu[6]+0.3061862178478971*alphaMu[5]-0.1767766952966368*alphaMu[3]+0.1767766952966368*alphaMu[2]-0.1767766952966368*alphaMu[1]+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*(alphaMu[16]+alphaMu[8]+alphaMu[7]))+0.1767766952966368*alphaMu[6]+0.3061862178478971*alphaMu[5]-0.1767766952966368*alphaMu[3]+0.1767766952966368*(alphaMu[2]+alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*alphaMu[16]-0.1767766952966368*(alphaMu[8]+alphaMu[7])+0.1767766952966368*alphaMu[6]+0.3061862178478971*alphaMu[5]+0.1767766952966368*alphaMu[3]-0.1767766952966368*(alphaMu[2]+alphaMu[1])+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*(alphaMu[16]+alphaMu[8]))+0.1767766952966368*alphaMu[7]-0.1767766952966368*alphaMu[6]+0.3061862178478971*alphaMu[5]+0.1767766952966368*alphaMu[3]-0.1767766952966368*alphaMu[2]+0.1767766952966368*(alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*((-0.1767766952966368*alphaMu[16])+0.1767766952966368*alphaMu[8]-0.1767766952966368*(alphaMu[7]+alphaMu[6])+0.3061862178478971*alphaMu[5]+0.1767766952966368*(alphaMu[3]+alphaMu[2])-0.1767766952966368*alphaMu[1]+0.1767766952966368*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.03125*(0.1767766952966368*(alphaMu[16]+alphaMu[8]+alphaMu[7]+alphaMu[6])+0.3061862178478971*alphaMu[5]+0.1767766952966368*(alphaMu[3]+alphaMu[2]+alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaR); 

  out[4] += 0.3061862178478971*(alphaVpar[16]*f[16]+alphaVpar[8]*f[8]+alphaVpar[7]*f[7]+alphaVpar[6]*f[6]+alphaVpar[4]*f[4]+alphaVpar[3]*f[3]+alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[5] += 0.3061862178478971*(alphaMu[16]*f[16]+alphaMu[8]*f[8]+alphaMu[7]*f[7]+alphaMu[6]*f[6]+alphaMu[5]*f[5]+alphaMu[3]*f[3]+alphaMu[2]*f[2]+alphaMu[1]*f[1]+alphaMu[0]*f[0]); 
  out[9] += 0.3061862178478971*(alphaVpar[8]*f[16]+f[8]*alphaVpar[16]+alphaVpar[4]*f[9]+alphaVpar[3]*f[7]+f[3]*alphaVpar[7]+alphaVpar[2]*f[6]+f[2]*alphaVpar[6]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 
  out[10] += 0.3061862178478971*(alphaVpar[7]*f[16]+f[7]*alphaVpar[16]+alphaVpar[4]*f[10]+alphaVpar[3]*f[8]+f[3]*alphaVpar[8]+alphaVpar[1]*f[6]+f[1]*alphaVpar[6]+alphaVpar[0]*f[2]+f[0]*alphaVpar[2]); 
  out[11] += 0.3061862178478971*(alphaVpar[6]*f[16]+f[6]*alphaVpar[16]+alphaVpar[4]*f[11]+alphaVpar[2]*f[8]+f[2]*alphaVpar[8]+alphaVpar[1]*f[7]+f[1]*alphaVpar[7]+alphaVpar[0]*f[3]+f[0]*alphaVpar[3]); 
  out[12] += 0.3061862178478971*(alphaMu[8]*f[16]+f[8]*alphaMu[16]+alphaMu[5]*f[12]+alphaMu[3]*f[7]+f[3]*alphaMu[7]+alphaMu[2]*f[6]+f[2]*alphaMu[6]+alphaMu[0]*f[1]+f[0]*alphaMu[1]); 
  out[13] += 0.3061862178478971*(alphaMu[7]*f[16]+f[7]*alphaMu[16]+alphaMu[5]*f[13]+alphaMu[3]*f[8]+f[3]*alphaMu[8]+alphaMu[1]*f[6]+f[1]*alphaMu[6]+alphaMu[0]*f[2]+f[0]*alphaMu[2]); 
  out[14] += 0.3061862178478971*(alphaMu[6]*f[16]+f[6]*alphaMu[16]+alphaMu[5]*f[14]+alphaMu[2]*f[8]+f[2]*alphaMu[8]+alphaMu[1]*f[7]+f[1]*alphaMu[7]+alphaMu[0]*f[3]+f[0]*alphaMu[3]); 
  out[15] += 0.3061862178478971*(alphaVpar[16]*f[27]+alphaMu[16]*f[26]+alphaVpar[8]*f[22]+alphaVpar[7]*f[21]+alphaVpar[6]*f[20]+alphaMu[8]*f[19]+alphaMu[7]*f[18]+alphaMu[6]*f[17]+(alphaMu[5]+alphaVpar[4])*f[15]+alphaVpar[3]*f[14]+alphaVpar[2]*f[13]+alphaVpar[1]*f[12]+alphaMu[3]*f[11]+alphaMu[2]*f[10]+alphaMu[1]*f[9]+alphaVpar[0]*f[5]+alphaMu[0]*f[4]); 
  out[17] += 0.3061862178478971*(alphaVpar[4]*f[17]+alphaVpar[3]*f[16]+f[3]*alphaVpar[16]+alphaVpar[7]*f[8]+f[7]*alphaVpar[8]+alphaVpar[0]*f[6]+f[0]*alphaVpar[6]+alphaVpar[1]*f[2]+f[1]*alphaVpar[2]); 
  out[18] += 0.3061862178478971*(alphaVpar[4]*f[18]+alphaVpar[2]*f[16]+f[2]*alphaVpar[16]+alphaVpar[6]*f[8]+f[6]*alphaVpar[8]+alphaVpar[0]*f[7]+f[0]*alphaVpar[7]+alphaVpar[1]*f[3]+f[1]*alphaVpar[3]); 
  out[19] += 0.3061862178478971*(alphaVpar[4]*f[19]+alphaVpar[1]*f[16]+f[1]*alphaVpar[16]+alphaVpar[0]*f[8]+f[0]*alphaVpar[8]+alphaVpar[6]*f[7]+f[6]*alphaVpar[7]+alphaVpar[2]*f[3]+f[2]*alphaVpar[3]); 
  out[20] += 0.3061862178478971*(alphaMu[5]*f[20]+alphaMu[3]*f[16]+f[3]*alphaMu[16]+alphaMu[7]*f[8]+f[7]*alphaMu[8]+alphaMu[0]*f[6]+f[0]*alphaMu[6]+alphaMu[1]*f[2]+f[1]*alphaMu[2]); 
  out[21] += 0.3061862178478971*(alphaMu[5]*f[21]+alphaMu[2]*f[16]+f[2]*alphaMu[16]+alphaMu[6]*f[8]+f[6]*alphaMu[8]+alphaMu[0]*f[7]+f[0]*alphaMu[7]+alphaMu[1]*f[3]+f[1]*alphaMu[3]); 
  out[22] += 0.3061862178478971*(alphaMu[5]*f[22]+alphaMu[1]*f[16]+f[1]*alphaMu[16]+alphaMu[0]*f[8]+f[0]*alphaMu[8]+alphaMu[6]*f[7]+f[6]*alphaMu[7]+alphaMu[2]*f[3]+f[2]*alphaMu[3]); 
  out[23] += 0.3061862178478971*(alphaVpar[8]*f[27]+alphaMu[8]*f[26]+(alphaMu[5]+alphaVpar[4])*f[23]+alphaVpar[16]*f[22]+alphaVpar[3]*f[21]+alphaVpar[2]*f[20]+alphaMu[16]*f[19]+alphaMu[3]*f[18]+alphaMu[2]*f[17]+alphaVpar[7]*f[14]+alphaVpar[6]*f[13]+alphaVpar[0]*f[12]+alphaMu[7]*f[11]+alphaMu[6]*f[10]+alphaMu[0]*f[9]+alphaVpar[1]*f[5]+alphaMu[1]*f[4]); 
  out[24] += 0.3061862178478971*(alphaVpar[7]*f[27]+alphaMu[7]*f[26]+(alphaMu[5]+alphaVpar[4])*f[24]+alphaVpar[3]*f[22]+alphaVpar[16]*f[21]+alphaVpar[1]*f[20]+alphaMu[3]*f[19]+alphaMu[16]*f[18]+alphaMu[1]*f[17]+alphaVpar[8]*f[14]+alphaVpar[0]*f[13]+alphaVpar[6]*f[12]+alphaMu[8]*f[11]+alphaMu[0]*f[10]+alphaMu[6]*f[9]+alphaVpar[2]*f[5]+alphaMu[2]*f[4]); 
  out[25] += 0.3061862178478971*(alphaVpar[6]*f[27]+alphaMu[6]*f[26]+(alphaMu[5]+alphaVpar[4])*f[25]+alphaVpar[2]*f[22]+alphaVpar[1]*f[21]+alphaVpar[16]*f[20]+alphaMu[2]*f[19]+alphaMu[1]*f[18]+alphaMu[16]*f[17]+alphaVpar[0]*f[14]+alphaVpar[8]*f[13]+alphaVpar[7]*f[12]+alphaMu[0]*f[11]+alphaMu[8]*f[10]+alphaMu[7]*f[9]+alphaVpar[3]*f[5]+alphaMu[3]*f[4]); 
  out[26] += 0.3061862178478971*(alphaVpar[4]*f[26]+alphaVpar[0]*f[16]+f[0]*alphaVpar[16]+alphaVpar[1]*f[8]+f[1]*alphaVpar[8]+alphaVpar[2]*f[7]+f[2]*alphaVpar[7]+alphaVpar[3]*f[6]+f[3]*alphaVpar[6]); 
  out[27] += 0.3061862178478971*(alphaMu[5]*f[27]+alphaMu[0]*f[16]+f[0]*alphaMu[16]+alphaMu[1]*f[8]+f[1]*alphaMu[8]+alphaMu[2]*f[7]+f[2]*alphaMu[7]+alphaMu[3]*f[6]+f[3]*alphaMu[6]); 
  out[28] += 0.3061862178478971*((alphaMu[5]+alphaVpar[4])*f[28]+alphaVpar[3]*f[27]+alphaMu[3]*f[26]+alphaVpar[7]*f[22]+alphaVpar[8]*f[21]+alphaVpar[0]*f[20]+alphaMu[7]*f[19]+alphaMu[8]*f[18]+alphaMu[0]*f[17]+f[14]*alphaVpar[16]+f[11]*alphaMu[16]+alphaVpar[1]*f[13]+alphaVpar[2]*f[12]+alphaMu[1]*f[10]+alphaMu[2]*f[9]+f[5]*alphaVpar[6]+f[4]*alphaMu[6]); 
  out[29] += 0.3061862178478971*((alphaMu[5]+alphaVpar[4])*f[29]+alphaVpar[2]*f[27]+alphaMu[2]*f[26]+alphaVpar[6]*f[22]+alphaVpar[0]*f[21]+alphaVpar[8]*f[20]+alphaMu[6]*f[19]+alphaMu[0]*f[18]+alphaMu[8]*f[17]+f[13]*alphaVpar[16]+f[10]*alphaMu[16]+alphaVpar[1]*f[14]+alphaVpar[3]*f[12]+alphaMu[1]*f[11]+alphaMu[3]*f[9]+f[5]*alphaVpar[7]+f[4]*alphaMu[7]); 
  out[30] += 0.3061862178478971*((alphaMu[5]+alphaVpar[4])*f[30]+alphaVpar[1]*f[27]+alphaMu[1]*f[26]+alphaVpar[0]*f[22]+alphaVpar[6]*f[21]+alphaVpar[7]*f[20]+alphaMu[0]*f[19]+alphaMu[6]*f[18]+alphaMu[7]*f[17]+f[12]*alphaVpar[16]+f[9]*alphaMu[16]+alphaVpar[2]*f[14]+alphaVpar[3]*f[13]+alphaMu[2]*f[11]+alphaMu[3]*f[10]+f[5]*alphaVpar[8]+f[4]*alphaMu[8]); 
  out[31] += 0.3061862178478971*((alphaMu[5]+alphaVpar[4])*f[31]+alphaVpar[0]*f[27]+alphaMu[0]*f[26]+alphaVpar[1]*f[22]+alphaVpar[2]*f[21]+alphaVpar[3]*f[20]+alphaMu[1]*f[19]+alphaMu[2]*f[18]+alphaMu[3]*f[17]+f[5]*alphaVpar[16]+f[4]*alphaMu[16]+alphaVpar[6]*f[14]+alphaVpar[7]*f[13]+alphaVpar[8]*f[12]+alphaMu[6]*f[11]+alphaMu[7]*f[10]+alphaMu[8]*f[9]); 

  return std::abs(0.1767766952966368*alphaVpar[0]) + 2.0*rdv2[1]*w[4]*nuSum+rdvSq4[1]*w[4]*(0.3333333333333333*BmagInv[7]*nuVtSqSum[7]+0.3333333333333333*BmagInv[6]*nuVtSqSum[6]+0.3333333333333333*BmagInv[5]*nuVtSqSum[5]+0.3333333333333333*BmagInv[4]*nuVtSqSum[4]+0.3333333333333333*BmagInv[3]*nuVtSqSum[3]+0.3333333333333333*BmagInv[2]*nuVtSqSum[2]+0.3333333333333333*BmagInv[1]*nuVtSqSum[1]+0.3333333333333333*BmagInv[0]*nuVtSqSum[0])*m_+0.4714045207910312*nuVtSqSum[0]*rdvSq4[0]; 

} 
