#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmCrossPrimMomentsGreene3x3vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[24]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
  uSMuO[2] = uSelf[2]-1.0*uOther[2]; 
  uSMuO[3] = uSelf[3]-1.0*uOther[3]; 
  uSMuO[4] = uSelf[4]-1.0*uOther[4]; 
  uSMuO[5] = uSelf[5]-1.0*uOther[5]; 
  uSMuO[6] = uSelf[6]-1.0*uOther[6]; 
  uSMuO[7] = uSelf[7]-1.0*uOther[7]; 
  uSMuO[8] = uSelf[8]-1.0*uOther[8]; 
  uSMuO[9] = uSelf[9]-1.0*uOther[9]; 
  uSMuO[10] = uSelf[10]-1.0*uOther[10]; 
  uSMuO[11] = uSelf[11]-1.0*uOther[11]; 
  uSMuO[12] = uSelf[12]-1.0*uOther[12]; 
  uSMuO[13] = uSelf[13]-1.0*uOther[13]; 
  uSMuO[14] = uSelf[14]-1.0*uOther[14]; 
  uSMuO[15] = uSelf[15]-1.0*uOther[15]; 
  uSMuO[16] = uSelf[16]-1.0*uOther[16]; 
  uSMuO[17] = uSelf[17]-1.0*uOther[17]; 
  uSMuO[18] = uSelf[18]-1.0*uOther[18]; 
  uSMuO[19] = uSelf[19]-1.0*uOther[19]; 
  uSMuO[20] = uSelf[20]-1.0*uOther[20]; 
  uSMuO[21] = uSelf[21]-1.0*uOther[21]; 
  uSMuO[22] = uSelf[22]-1.0*uOther[22]; 
  uSMuO[23] = uSelf[23]-1.0*uOther[23]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[8]; 
  for (unsigned short int k=0; k<8; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    uSMuOSq[0] += 0.3535533905932737*uSMuO[a0+7]*uSMuO[a0+7]+0.3535533905932737*uSMuO[a0+6]*uSMuO[a0+6]+0.3535533905932737*uSMuO[a0+5]*uSMuO[a0+5]+0.3535533905932737*uSMuO[a0+4]*uSMuO[a0+4]+0.3535533905932737*uSMuO[a0+3]*uSMuO[a0+3]+0.3535533905932737*uSMuO[a0+2]*uSMuO[a0+2]+0.3535533905932737*uSMuO[a0+1]*uSMuO[a0+1]+0.3535533905932737*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += 0.7071067811865475*uSMuO[a0+6]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+4]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+1]; 
    uSMuOSq[2] += 0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+4]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+2]; 
    uSMuOSq[3] += 0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+3]; 
    uSMuOSq[4] += 0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+4]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+2]; 
    uSMuOSq[5] += 0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+3]; 
    uSMuOSq[6] += 0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+3]; 
    uSMuOSq[7] += 0.7071067811865475*uSMuO[a0]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+4]; 
  } 
 
  // ..... Get the cross flow velocity uSelf2 ..... // 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    uCross[a0] = 0.5*(uSelf[a0]+uOther[a0])-0.5*uSMuO[a0]*beta; 
    uCross[a0+1] = 0.5*(uSelf[a0+1]+uOther[a0+1])-0.5*uSMuO[a0+1]*beta; 
    uCross[a0+2] = 0.5*(uSelf[a0+2]+uOther[a0+2])-0.5*uSMuO[a0+2]*beta; 
    uCross[a0+3] = 0.5*(uSelf[a0+3]+uOther[a0+3])-0.5*uSMuO[a0+3]*beta; 
    uCross[a0+4] = 0.5*(uSelf[a0+4]+uOther[a0+4])-0.5*uSMuO[a0+4]*beta; 
    uCross[a0+5] = 0.5*(uSelf[a0+5]+uOther[a0+5])-0.5*uSMuO[a0+5]*beta; 
    uCross[a0+6] = 0.5*(uSelf[a0+6]+uOther[a0+6])-0.5*uSMuO[a0+6]*beta; 
    uCross[a0+7] = 0.5*(uSelf[a0+7]+uOther[a0+7])-0.5*uSMuO[a0+7]*beta; 
 
  } 
 
  double mBetaFrac = (0.5*(beta+1.0))/(mRat+1.0); 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = (vtSqOther[0]+0.3333333333333333*uSMuOSq[0])*mBetaFrac*mRat-1.0*vtSqSelf[0]*mBetaFrac+vtSqSelf[0]; 
  vtSqCross[1] = (vtSqOther[1]+0.3333333333333333*uSMuOSq[1])*mBetaFrac*mRat-1.0*vtSqSelf[1]*mBetaFrac+vtSqSelf[1]; 
  vtSqCross[2] = (vtSqOther[2]+0.3333333333333333*uSMuOSq[2])*mBetaFrac*mRat-1.0*vtSqSelf[2]*mBetaFrac+vtSqSelf[2]; 
  vtSqCross[3] = (vtSqOther[3]+0.3333333333333333*uSMuOSq[3])*mBetaFrac*mRat-1.0*vtSqSelf[3]*mBetaFrac+vtSqSelf[3]; 
  vtSqCross[4] = (vtSqOther[4]+0.3333333333333333*uSMuOSq[4])*mBetaFrac*mRat-1.0*vtSqSelf[4]*mBetaFrac+vtSqSelf[4]; 
  vtSqCross[5] = (vtSqOther[5]+0.3333333333333333*uSMuOSq[5])*mBetaFrac*mRat-1.0*vtSqSelf[5]*mBetaFrac+vtSqSelf[5]; 
  vtSqCross[6] = (vtSqOther[6]+0.3333333333333333*uSMuOSq[6])*mBetaFrac*mRat-1.0*vtSqSelf[6]*mBetaFrac+vtSqSelf[6]; 
  vtSqCross[7] = (vtSqOther[7]+0.3333333333333333*uSMuOSq[7])*mBetaFrac*mRat-1.0*vtSqSelf[7]*mBetaFrac+vtSqSelf[7]; 
 
} 
 
void VmCrossPrimMomentsGreene3x3vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[60]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
  uSMuO[2] = uSelf[2]-1.0*uOther[2]; 
  uSMuO[3] = uSelf[3]-1.0*uOther[3]; 
  uSMuO[4] = uSelf[4]-1.0*uOther[4]; 
  uSMuO[5] = uSelf[5]-1.0*uOther[5]; 
  uSMuO[6] = uSelf[6]-1.0*uOther[6]; 
  uSMuO[7] = uSelf[7]-1.0*uOther[7]; 
  uSMuO[8] = uSelf[8]-1.0*uOther[8]; 
  uSMuO[9] = uSelf[9]-1.0*uOther[9]; 
  uSMuO[10] = uSelf[10]-1.0*uOther[10]; 
  uSMuO[11] = uSelf[11]-1.0*uOther[11]; 
  uSMuO[12] = uSelf[12]-1.0*uOther[12]; 
  uSMuO[13] = uSelf[13]-1.0*uOther[13]; 
  uSMuO[14] = uSelf[14]-1.0*uOther[14]; 
  uSMuO[15] = uSelf[15]-1.0*uOther[15]; 
  uSMuO[16] = uSelf[16]-1.0*uOther[16]; 
  uSMuO[17] = uSelf[17]-1.0*uOther[17]; 
  uSMuO[18] = uSelf[18]-1.0*uOther[18]; 
  uSMuO[19] = uSelf[19]-1.0*uOther[19]; 
  uSMuO[20] = uSelf[20]-1.0*uOther[20]; 
  uSMuO[21] = uSelf[21]-1.0*uOther[21]; 
  uSMuO[22] = uSelf[22]-1.0*uOther[22]; 
  uSMuO[23] = uSelf[23]-1.0*uOther[23]; 
  uSMuO[24] = uSelf[24]-1.0*uOther[24]; 
  uSMuO[25] = uSelf[25]-1.0*uOther[25]; 
  uSMuO[26] = uSelf[26]-1.0*uOther[26]; 
  uSMuO[27] = uSelf[27]-1.0*uOther[27]; 
  uSMuO[28] = uSelf[28]-1.0*uOther[28]; 
  uSMuO[29] = uSelf[29]-1.0*uOther[29]; 
  uSMuO[30] = uSelf[30]-1.0*uOther[30]; 
  uSMuO[31] = uSelf[31]-1.0*uOther[31]; 
  uSMuO[32] = uSelf[32]-1.0*uOther[32]; 
  uSMuO[33] = uSelf[33]-1.0*uOther[33]; 
  uSMuO[34] = uSelf[34]-1.0*uOther[34]; 
  uSMuO[35] = uSelf[35]-1.0*uOther[35]; 
  uSMuO[36] = uSelf[36]-1.0*uOther[36]; 
  uSMuO[37] = uSelf[37]-1.0*uOther[37]; 
  uSMuO[38] = uSelf[38]-1.0*uOther[38]; 
  uSMuO[39] = uSelf[39]-1.0*uOther[39]; 
  uSMuO[40] = uSelf[40]-1.0*uOther[40]; 
  uSMuO[41] = uSelf[41]-1.0*uOther[41]; 
  uSMuO[42] = uSelf[42]-1.0*uOther[42]; 
  uSMuO[43] = uSelf[43]-1.0*uOther[43]; 
  uSMuO[44] = uSelf[44]-1.0*uOther[44]; 
  uSMuO[45] = uSelf[45]-1.0*uOther[45]; 
  uSMuO[46] = uSelf[46]-1.0*uOther[46]; 
  uSMuO[47] = uSelf[47]-1.0*uOther[47]; 
  uSMuO[48] = uSelf[48]-1.0*uOther[48]; 
  uSMuO[49] = uSelf[49]-1.0*uOther[49]; 
  uSMuO[50] = uSelf[50]-1.0*uOther[50]; 
  uSMuO[51] = uSelf[51]-1.0*uOther[51]; 
  uSMuO[52] = uSelf[52]-1.0*uOther[52]; 
  uSMuO[53] = uSelf[53]-1.0*uOther[53]; 
  uSMuO[54] = uSelf[54]-1.0*uOther[54]; 
  uSMuO[55] = uSelf[55]-1.0*uOther[55]; 
  uSMuO[56] = uSelf[56]-1.0*uOther[56]; 
  uSMuO[57] = uSelf[57]-1.0*uOther[57]; 
  uSMuO[58] = uSelf[58]-1.0*uOther[58]; 
  uSMuO[59] = uSelf[59]-1.0*uOther[59]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[20]; 
  for (unsigned short int k=0; k<20; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 20*vd; 
    uSMuOSq[0] += 0.3535533905932737*uSMuO[a0+19]*uSMuO[a0+19]+0.3535533905932737*uSMuO[a0+18]*uSMuO[a0+18]+0.3535533905932737*uSMuO[a0+17]*uSMuO[a0+17]+0.3535533905932737*uSMuO[a0+16]*uSMuO[a0+16]+0.3535533905932737*uSMuO[a0+15]*uSMuO[a0+15]+0.3535533905932737*uSMuO[a0+14]*uSMuO[a0+14]+0.3535533905932737*uSMuO[a0+13]*uSMuO[a0+13]+0.3535533905932737*uSMuO[a0+12]*uSMuO[a0+12]+0.3535533905932737*uSMuO[a0+11]*uSMuO[a0+11]+0.3535533905932737*uSMuO[a0+10]*uSMuO[a0+10]+0.3535533905932737*uSMuO[a0+9]*uSMuO[a0+9]+0.3535533905932737*uSMuO[a0+8]*uSMuO[a0+8]+0.3535533905932737*uSMuO[a0+7]*uSMuO[a0+7]+0.3535533905932737*uSMuO[a0+6]*uSMuO[a0+6]+0.3535533905932737*uSMuO[a0+5]*uSMuO[a0+5]+0.3535533905932737*uSMuO[a0+4]*uSMuO[a0+4]+0.3535533905932737*uSMuO[a0+3]*uSMuO[a0+3]+0.3535533905932737*uSMuO[a0+2]*uSMuO[a0+2]+0.3535533905932737*uSMuO[a0+1]*uSMuO[a0+1]+0.3535533905932737*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += 0.7071067811865475*uSMuO[a0+16]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0+14]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+10]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+9]*uSMuO[a0+15]+0.632455532033676*uSMuO[a0+5]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0+8]*uSMuO[a0+12]+0.632455532033676*uSMuO[a0+4]*uSMuO[a0+11]+0.7071067811865475*uSMuO[a0+6]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+1]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+4]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+1]; 
    uSMuOSq[2] += 0.7071067811865475*uSMuO[a0+15]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+10]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0+13]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+9]*uSMuO[a0+16]+0.632455532033676*uSMuO[a0+6]*uSMuO[a0+14]+0.632455532033676*uSMuO[a0+4]*uSMuO[a0+12]+0.7071067811865475*uSMuO[a0+7]*uSMuO[a0+11]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+2]*uSMuO[a0+8]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+4]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+2]; 
    uSMuOSq[3] += 0.6324555320336759*uSMuO[a0+10]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0+12]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0+11]*uSMuO[a0+17]+0.632455532033676*uSMuO[a0+6]*uSMuO[a0+16]+0.632455532033676*uSMuO[a0+5]*uSMuO[a0+15]+0.7071067811865475*uSMuO[a0+8]*uSMuO[a0+14]+0.7071067811865475*uSMuO[a0+7]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+3]*uSMuO[a0+9]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+3]; 
    uSMuOSq[4] += 0.7071067811865475*uSMuO[a0+9]*uSMuO[a0+19]+0.5656854249492381*uSMuO[a0+17]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+6]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+5]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+15]*uSMuO[a0+16]+0.632455532033676*uSMuO[a0+10]*uSMuO[a0+14]+0.632455532033676*uSMuO[a0+10]*uSMuO[a0+13]+0.5656854249492381*uSMuO[a0+11]*uSMuO[a0+12]+0.632455532033676*uSMuO[a0+2]*uSMuO[a0+12]+0.632455532033676*uSMuO[a0+1]*uSMuO[a0+11]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+4]*uSMuO[a0+8]+0.6324555320336759*uSMuO[a0+4]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+4]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+2]; 
    uSMuOSq[5] += 0.5656854249492381*uSMuO[a0+17]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+6]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0+8]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+4]*uSMuO[a0+17]+0.632455532033676*uSMuO[a0+10]*uSMuO[a0+16]+0.5656854249492381*uSMuO[a0+13]*uSMuO[a0+15]+0.632455532033676*uSMuO[a0+3]*uSMuO[a0+15]+0.7071067811865475*uSMuO[a0+12]*uSMuO[a0+14]+0.632455532033676*uSMuO[a0+1]*uSMuO[a0+13]+0.632455532033676*uSMuO[a0+10]*uSMuO[a0+11]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+5]*uSMuO[a0+9]+0.6324555320336759*uSMuO[a0+5]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+3]; 
    uSMuOSq[6] += 0.5656854249492381*uSMuO[a0+18]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+5]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+4]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0+7]*uSMuO[a0+17]+0.5656854249492381*uSMuO[a0+14]*uSMuO[a0+16]+0.632455532033676*uSMuO[a0+3]*uSMuO[a0+16]+0.632455532033676*uSMuO[a0+10]*uSMuO[a0+15]+0.632455532033676*uSMuO[a0+2]*uSMuO[a0+14]+0.7071067811865475*uSMuO[a0+11]*uSMuO[a0+13]+0.632455532033676*uSMuO[a0+10]*uSMuO[a0+12]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+6]*uSMuO[a0+9]+0.6324555320336759*uSMuO[a0+6]*uSMuO[a0+8]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+3]; 
    uSMuOSq[7] += 0.3162277660168379*uSMuO[a0+19]*uSMuO[a0+19]+0.3162277660168379*uSMuO[a0+18]*uSMuO[a0+18]+0.2258769757263128*uSMuO[a0+17]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+6]*uSMuO[a0+17]+0.3162277660168379*uSMuO[a0+15]*uSMuO[a0+15]+0.2258769757263128*uSMuO[a0+13]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+13]+0.3162277660168379*uSMuO[a0+12]*uSMuO[a0+12]+0.2258769757263128*uSMuO[a0+11]*uSMuO[a0+11]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+11]+0.3162277660168379*uSMuO[a0+10]*uSMuO[a0+10]+0.2258769757263128*uSMuO[a0+7]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+7]+0.3162277660168379*uSMuO[a0+5]*uSMuO[a0+5]+0.3162277660168379*uSMuO[a0+4]*uSMuO[a0+4]+0.3162277660168379*uSMuO[a0+1]*uSMuO[a0+1]; 
    uSMuOSq[8] += 0.3162277660168379*uSMuO[a0+19]*uSMuO[a0+19]+0.2258769757263128*uSMuO[a0+18]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+18]+0.3162277660168379*uSMuO[a0+17]*uSMuO[a0+17]+0.3162277660168379*uSMuO[a0+16]*uSMuO[a0+16]+0.2258769757263128*uSMuO[a0+14]*uSMuO[a0+14]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+14]+0.2258769757263128*uSMuO[a0+12]*uSMuO[a0+12]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+12]+0.3162277660168379*uSMuO[a0+11]*uSMuO[a0+11]+0.3162277660168379*uSMuO[a0+10]*uSMuO[a0+10]+0.2258769757263128*uSMuO[a0+8]*uSMuO[a0+8]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+8]+0.3162277660168379*uSMuO[a0+6]*uSMuO[a0+6]+0.3162277660168379*uSMuO[a0+4]*uSMuO[a0+4]+0.3162277660168379*uSMuO[a0+2]*uSMuO[a0+2]; 
    uSMuOSq[9] += 0.2258769757263128*uSMuO[a0+19]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+19]+0.3162277660168379*uSMuO[a0+18]*uSMuO[a0+18]+0.3162277660168379*uSMuO[a0+17]*uSMuO[a0+17]+0.2258769757263128*uSMuO[a0+16]*uSMuO[a0+16]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+16]+0.2258769757263128*uSMuO[a0+15]*uSMuO[a0+15]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+15]+0.3162277660168379*uSMuO[a0+14]*uSMuO[a0+14]+0.3162277660168379*uSMuO[a0+13]*uSMuO[a0+13]+0.3162277660168379*uSMuO[a0+10]*uSMuO[a0+10]+0.2258769757263128*uSMuO[a0+9]*uSMuO[a0+9]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+9]+0.3162277660168379*uSMuO[a0+6]*uSMuO[a0+6]+0.3162277660168379*uSMuO[a0+5]*uSMuO[a0+5]+0.3162277660168379*uSMuO[a0+3]*uSMuO[a0+3]; 
    uSMuOSq[10] += 0.5656854249492382*uSMuO[a0+14]*uSMuO[a0+19]+0.5656854249492382*uSMuO[a0+13]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+3]*uSMuO[a0+19]+0.5656854249492382*uSMuO[a0+16]*uSMuO[a0+18]+0.5656854249492382*uSMuO[a0+11]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+2]*uSMuO[a0+18]+0.5656854249492382*uSMuO[a0+15]*uSMuO[a0+17]+0.5656854249492382*uSMuO[a0+12]*uSMuO[a0+17]+0.6324555320336759*uSMuO[a0+1]*uSMuO[a0+17]+0.632455532033676*uSMuO[a0+5]*uSMuO[a0+16]+0.632455532033676*uSMuO[a0+6]*uSMuO[a0+15]+0.632455532033676*uSMuO[a0+4]*uSMuO[a0+14]+0.632455532033676*uSMuO[a0+4]*uSMuO[a0+13]+0.632455532033676*uSMuO[a0+6]*uSMuO[a0+12]+0.632455532033676*uSMuO[a0+5]*uSMuO[a0+11]+0.6324555320336759*uSMuO[a0+9]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+8]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+7]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+4]; 
    uSMuOSq[11] += 0.6324555320336759*uSMuO[a0+15]*uSMuO[a0+19]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+14]*uSMuO[a0+17]+0.4517539514526256*uSMuO[a0+13]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+6]*uSMuO[a0+13]+0.5656854249492381*uSMuO[a0+4]*uSMuO[a0+12]+0.6324555320336759*uSMuO[a0+8]*uSMuO[a0+11]+0.4517539514526256*uSMuO[a0+7]*uSMuO[a0+11]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+11]+0.632455532033676*uSMuO[a0+5]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+7]+0.632455532033676*uSMuO[a0+1]*uSMuO[a0+4]; 
    uSMuOSq[12] += 0.6324555320336759*uSMuO[a0+16]*uSMuO[a0+19]+0.4517539514526256*uSMuO[a0+14]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+13]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+18]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+14]+0.4517539514526256*uSMuO[a0+8]*uSMuO[a0+12]+0.6324555320336759*uSMuO[a0+7]*uSMuO[a0+12]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+12]+0.5656854249492381*uSMuO[a0+4]*uSMuO[a0+11]+0.632455532033676*uSMuO[a0+6]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+8]+0.632455532033676*uSMuO[a0+2]*uSMuO[a0+4]; 
    uSMuOSq[13] += 0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+12]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+16]*uSMuO[a0+17]+0.4517539514526256*uSMuO[a0+11]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+17]+0.5656854249492381*uSMuO[a0+5]*uSMuO[a0+15]+0.6324555320336759*uSMuO[a0+9]*uSMuO[a0+13]+0.4517539514526256*uSMuO[a0+7]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0+6]*uSMuO[a0+11]+0.632455532033676*uSMuO[a0+4]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+7]+0.632455532033676*uSMuO[a0+1]*uSMuO[a0+5]; 
    uSMuOSq[14] += 0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+15]*uSMuO[a0+18]+0.4517539514526256*uSMuO[a0+12]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+11]*uSMuO[a0+17]+0.5656854249492381*uSMuO[a0+6]*uSMuO[a0+16]+0.6324555320336759*uSMuO[a0+9]*uSMuO[a0+14]+0.4517539514526256*uSMuO[a0+8]*uSMuO[a0+14]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+14]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+12]+0.632455532033676*uSMuO[a0+4]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+8]+0.632455532033676*uSMuO[a0+2]*uSMuO[a0+6]; 
    uSMuOSq[15] += 0.4517539514526256*uSMuO[a0+16]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+11]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+14]*uSMuO[a0+18]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+16]+0.4517539514526256*uSMuO[a0+9]*uSMuO[a0+15]+0.6324555320336759*uSMuO[a0+7]*uSMuO[a0+15]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+15]+0.5656854249492381*uSMuO[a0+5]*uSMuO[a0+13]+0.632455532033676*uSMuO[a0+6]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+9]+0.632455532033676*uSMuO[a0+3]*uSMuO[a0+5]; 
    uSMuOSq[16] += 0.4517539514526256*uSMuO[a0+15]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+12]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+19]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+13]*uSMuO[a0+17]+0.4517539514526256*uSMuO[a0+9]*uSMuO[a0+16]+0.6324555320336759*uSMuO[a0+8]*uSMuO[a0+16]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+16]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+15]+0.5656854249492381*uSMuO[a0+6]*uSMuO[a0+14]+0.632455532033676*uSMuO[a0+5]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+9]+0.632455532033676*uSMuO[a0+3]*uSMuO[a0+6]; 
    uSMuOSq[17] += 0.5059644256269408*uSMuO[a0+18]*uSMuO[a0+19]+0.5656854249492381*uSMuO[a0+5]*uSMuO[a0+19]+0.5656854249492381*uSMuO[a0+4]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+9]*uSMuO[a0+17]+0.6324555320336759*uSMuO[a0+8]*uSMuO[a0+17]+0.4517539514526256*uSMuO[a0+7]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+17]+0.6324555320336759*uSMuO[a0+13]*uSMuO[a0+16]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+15]+0.6324555320336759*uSMuO[a0+11]*uSMuO[a0+14]+0.4517539514526256*uSMuO[a0+11]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+13]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+12]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+11]+0.6324555320336759*uSMuO[a0+1]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+6]*uSMuO[a0+7]+0.6324555320336759*uSMuO[a0+4]*uSMuO[a0+5]; 
    uSMuOSq[18] += 0.5059644256269408*uSMuO[a0+17]*uSMuO[a0+19]+0.5656854249492381*uSMuO[a0+6]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+9]*uSMuO[a0+18]+0.4517539514526256*uSMuO[a0+8]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+7]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+18]+0.5656854249492381*uSMuO[a0+4]*uSMuO[a0+17]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+16]+0.6324555320336759*uSMuO[a0+14]*uSMuO[a0+15]+0.4517539514526256*uSMuO[a0+12]*uSMuO[a0+14]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+14]+0.6324555320336759*uSMuO[a0+12]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+12]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+11]+0.6324555320336759*uSMuO[a0+2]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+8]+0.6324555320336759*uSMuO[a0+4]*uSMuO[a0+6]; 
    uSMuOSq[19] += 0.4517539514526256*uSMuO[a0+9]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+8]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+7]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+19]+0.5059644256269408*uSMuO[a0+17]*uSMuO[a0+18]+0.5656854249492381*uSMuO[a0+6]*uSMuO[a0+18]+0.5656854249492381*uSMuO[a0+5]*uSMuO[a0+17]+0.4517539514526256*uSMuO[a0+15]*uSMuO[a0+16]+0.6324555320336759*uSMuO[a0+12]*uSMuO[a0+16]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+16]+0.6324555320336759*uSMuO[a0+11]*uSMuO[a0+15]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+15]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+14]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+13]+0.6324555320336759*uSMuO[a0+3]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+9]+0.6324555320336759*uSMuO[a0+5]*uSMuO[a0+6]; 
  } 
 
  // ..... Get the cross flow velocity uSelf2 ..... // 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 20*vd; 
    uCross[a0] = 0.5*(uSelf[a0]+uOther[a0])-0.5*uSMuO[a0]*beta; 
    uCross[a0+1] = 0.5*(uSelf[a0+1]+uOther[a0+1])-0.5*uSMuO[a0+1]*beta; 
    uCross[a0+2] = 0.5*(uSelf[a0+2]+uOther[a0+2])-0.5*uSMuO[a0+2]*beta; 
    uCross[a0+3] = 0.5*(uSelf[a0+3]+uOther[a0+3])-0.5*uSMuO[a0+3]*beta; 
    uCross[a0+4] = 0.5*(uSelf[a0+4]+uOther[a0+4])-0.5*uSMuO[a0+4]*beta; 
    uCross[a0+5] = 0.5*(uSelf[a0+5]+uOther[a0+5])-0.5*uSMuO[a0+5]*beta; 
    uCross[a0+6] = 0.5*(uSelf[a0+6]+uOther[a0+6])-0.5*uSMuO[a0+6]*beta; 
    uCross[a0+7] = 0.5*(uSelf[a0+7]+uOther[a0+7])-0.5*uSMuO[a0+7]*beta; 
    uCross[a0+8] = 0.5*(uSelf[a0+8]+uOther[a0+8])-0.5*uSMuO[a0+8]*beta; 
    uCross[a0+9] = 0.5*(uSelf[a0+9]+uOther[a0+9])-0.5*uSMuO[a0+9]*beta; 
    uCross[a0+10] = 0.5*(uSelf[a0+10]+uOther[a0+10])-0.5*uSMuO[a0+10]*beta; 
    uCross[a0+11] = 0.5*(uSelf[a0+11]+uOther[a0+11])-0.5*uSMuO[a0+11]*beta; 
    uCross[a0+12] = 0.5*(uSelf[a0+12]+uOther[a0+12])-0.5*uSMuO[a0+12]*beta; 
    uCross[a0+13] = 0.5*(uSelf[a0+13]+uOther[a0+13])-0.5*uSMuO[a0+13]*beta; 
    uCross[a0+14] = 0.5*(uSelf[a0+14]+uOther[a0+14])-0.5*uSMuO[a0+14]*beta; 
    uCross[a0+15] = 0.5*(uSelf[a0+15]+uOther[a0+15])-0.5*uSMuO[a0+15]*beta; 
    uCross[a0+16] = 0.5*(uSelf[a0+16]+uOther[a0+16])-0.5*uSMuO[a0+16]*beta; 
    uCross[a0+17] = 0.5*(uSelf[a0+17]+uOther[a0+17])-0.5*uSMuO[a0+17]*beta; 
    uCross[a0+18] = 0.5*(uSelf[a0+18]+uOther[a0+18])-0.5*uSMuO[a0+18]*beta; 
    uCross[a0+19] = 0.5*(uSelf[a0+19]+uOther[a0+19])-0.5*uSMuO[a0+19]*beta; 
 
  } 
 
  double mBetaFrac = (0.5*(beta+1.0))/(mRat+1.0); 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = (vtSqOther[0]+0.3333333333333333*uSMuOSq[0])*mBetaFrac*mRat-1.0*vtSqSelf[0]*mBetaFrac+vtSqSelf[0]; 
  vtSqCross[1] = (vtSqOther[1]+0.3333333333333333*uSMuOSq[1])*mBetaFrac*mRat-1.0*vtSqSelf[1]*mBetaFrac+vtSqSelf[1]; 
  vtSqCross[2] = (vtSqOther[2]+0.3333333333333333*uSMuOSq[2])*mBetaFrac*mRat-1.0*vtSqSelf[2]*mBetaFrac+vtSqSelf[2]; 
  vtSqCross[3] = (vtSqOther[3]+0.3333333333333333*uSMuOSq[3])*mBetaFrac*mRat-1.0*vtSqSelf[3]*mBetaFrac+vtSqSelf[3]; 
  vtSqCross[4] = (vtSqOther[4]+0.3333333333333333*uSMuOSq[4])*mBetaFrac*mRat-1.0*vtSqSelf[4]*mBetaFrac+vtSqSelf[4]; 
  vtSqCross[5] = (vtSqOther[5]+0.3333333333333333*uSMuOSq[5])*mBetaFrac*mRat-1.0*vtSqSelf[5]*mBetaFrac+vtSqSelf[5]; 
  vtSqCross[6] = (vtSqOther[6]+0.3333333333333333*uSMuOSq[6])*mBetaFrac*mRat-1.0*vtSqSelf[6]*mBetaFrac+vtSqSelf[6]; 
  vtSqCross[7] = (vtSqOther[7]+0.3333333333333333*uSMuOSq[7])*mBetaFrac*mRat-1.0*vtSqSelf[7]*mBetaFrac+vtSqSelf[7]; 
  vtSqCross[8] = (vtSqOther[8]+0.3333333333333333*uSMuOSq[8])*mBetaFrac*mRat-1.0*vtSqSelf[8]*mBetaFrac+vtSqSelf[8]; 
  vtSqCross[9] = (vtSqOther[9]+0.3333333333333333*uSMuOSq[9])*mBetaFrac*mRat-1.0*vtSqSelf[9]*mBetaFrac+vtSqSelf[9]; 
  vtSqCross[10] = (vtSqOther[10]+0.3333333333333333*uSMuOSq[10])*mBetaFrac*mRat-1.0*vtSqSelf[10]*mBetaFrac+vtSqSelf[10]; 
  vtSqCross[11] = (vtSqOther[11]+0.3333333333333333*uSMuOSq[11])*mBetaFrac*mRat-1.0*vtSqSelf[11]*mBetaFrac+vtSqSelf[11]; 
  vtSqCross[12] = (vtSqOther[12]+0.3333333333333333*uSMuOSq[12])*mBetaFrac*mRat-1.0*vtSqSelf[12]*mBetaFrac+vtSqSelf[12]; 
  vtSqCross[13] = (vtSqOther[13]+0.3333333333333333*uSMuOSq[13])*mBetaFrac*mRat-1.0*vtSqSelf[13]*mBetaFrac+vtSqSelf[13]; 
  vtSqCross[14] = (vtSqOther[14]+0.3333333333333333*uSMuOSq[14])*mBetaFrac*mRat-1.0*vtSqSelf[14]*mBetaFrac+vtSqSelf[14]; 
  vtSqCross[15] = (vtSqOther[15]+0.3333333333333333*uSMuOSq[15])*mBetaFrac*mRat-1.0*vtSqSelf[15]*mBetaFrac+vtSqSelf[15]; 
  vtSqCross[16] = (vtSqOther[16]+0.3333333333333333*uSMuOSq[16])*mBetaFrac*mRat-1.0*vtSqSelf[16]*mBetaFrac+vtSqSelf[16]; 
  vtSqCross[17] = (vtSqOther[17]+0.3333333333333333*uSMuOSq[17])*mBetaFrac*mRat-1.0*vtSqSelf[17]*mBetaFrac+vtSqSelf[17]; 
  vtSqCross[18] = (vtSqOther[18]+0.3333333333333333*uSMuOSq[18])*mBetaFrac*mRat-1.0*vtSqSelf[18]*mBetaFrac+vtSqSelf[18]; 
  vtSqCross[19] = (vtSqOther[19]+0.3333333333333333*uSMuOSq[19])*mBetaFrac*mRat-1.0*vtSqSelf[19]*mBetaFrac+vtSqSelf[19]; 
 
} 
 
void VmCrossPrimMomentsHeavyIon3x3vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[24]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
  uSMuO[2] = uSelf[2]-1.0*uOther[2]; 
  uSMuO[3] = uSelf[3]-1.0*uOther[3]; 
  uSMuO[4] = uSelf[4]-1.0*uOther[4]; 
  uSMuO[5] = uSelf[5]-1.0*uOther[5]; 
  uSMuO[6] = uSelf[6]-1.0*uOther[6]; 
  uSMuO[7] = uSelf[7]-1.0*uOther[7]; 
  uSMuO[8] = uSelf[8]-1.0*uOther[8]; 
  uSMuO[9] = uSelf[9]-1.0*uOther[9]; 
  uSMuO[10] = uSelf[10]-1.0*uOther[10]; 
  uSMuO[11] = uSelf[11]-1.0*uOther[11]; 
  uSMuO[12] = uSelf[12]-1.0*uOther[12]; 
  uSMuO[13] = uSelf[13]-1.0*uOther[13]; 
  uSMuO[14] = uSelf[14]-1.0*uOther[14]; 
  uSMuO[15] = uSelf[15]-1.0*uOther[15]; 
  uSMuO[16] = uSelf[16]-1.0*uOther[16]; 
  uSMuO[17] = uSelf[17]-1.0*uOther[17]; 
  uSMuO[18] = uSelf[18]-1.0*uOther[18]; 
  uSMuO[19] = uSelf[19]-1.0*uOther[19]; 
  uSMuO[20] = uSelf[20]-1.0*uOther[20]; 
  uSMuO[21] = uSelf[21]-1.0*uOther[21]; 
  uSMuO[22] = uSelf[22]-1.0*uOther[22]; 
  uSMuO[23] = uSelf[23]-1.0*uOther[23]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[8]; 
  for (unsigned short int k=0; k<8; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    uSMuOSq[0] += 0.3535533905932737*uSMuO[a0+7]*uSMuO[a0+7]+0.3535533905932737*uSMuO[a0+6]*uSMuO[a0+6]+0.3535533905932737*uSMuO[a0+5]*uSMuO[a0+5]+0.3535533905932737*uSMuO[a0+4]*uSMuO[a0+4]+0.3535533905932737*uSMuO[a0+3]*uSMuO[a0+3]+0.3535533905932737*uSMuO[a0+2]*uSMuO[a0+2]+0.3535533905932737*uSMuO[a0+1]*uSMuO[a0+1]+0.3535533905932737*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += 0.7071067811865475*uSMuO[a0+6]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+4]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+1]; 
    uSMuOSq[2] += 0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+4]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+2]; 
    uSMuOSq[3] += 0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+3]; 
    uSMuOSq[4] += 0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+4]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+2]; 
    uSMuOSq[5] += 0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+3]; 
    uSMuOSq[6] += 0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+3]; 
    uSMuOSq[7] += 0.7071067811865475*uSMuO[a0]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+4]; 
  } 
 
  // ..... Get the cross flow velocity uCross ..... // 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    uCross[a0] = uOther[a0]; 
    uCross[a0+1] = uOther[a0+1]; 
    uCross[a0+2] = uOther[a0+2]; 
    uCross[a0+3] = uOther[a0+3]; 
    uCross[a0+4] = uOther[a0+4]; 
    uCross[a0+5] = uOther[a0+5]; 
    uCross[a0+6] = uOther[a0+6]; 
    uCross[a0+7] = uOther[a0+7]; 
 
  } 
 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = 0.3333333333333333*uSMuOSq[0]*mRat+vtSqSelf[0]; 
  vtSqCross[1] = 0.3333333333333333*uSMuOSq[1]*mRat+vtSqSelf[1]; 
  vtSqCross[2] = 0.3333333333333333*uSMuOSq[2]*mRat+vtSqSelf[2]; 
  vtSqCross[3] = 0.3333333333333333*uSMuOSq[3]*mRat+vtSqSelf[3]; 
  vtSqCross[4] = 0.3333333333333333*uSMuOSq[4]*mRat+vtSqSelf[4]; 
  vtSqCross[5] = 0.3333333333333333*uSMuOSq[5]*mRat+vtSqSelf[5]; 
  vtSqCross[6] = 0.3333333333333333*uSMuOSq[6]*mRat+vtSqSelf[6]; 
  vtSqCross[7] = 0.3333333333333333*uSMuOSq[7]*mRat+vtSqSelf[7]; 
 
} 
 
void VmCrossPrimMomentsHeavyIon3x3vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[60]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
  uSMuO[2] = uSelf[2]-1.0*uOther[2]; 
  uSMuO[3] = uSelf[3]-1.0*uOther[3]; 
  uSMuO[4] = uSelf[4]-1.0*uOther[4]; 
  uSMuO[5] = uSelf[5]-1.0*uOther[5]; 
  uSMuO[6] = uSelf[6]-1.0*uOther[6]; 
  uSMuO[7] = uSelf[7]-1.0*uOther[7]; 
  uSMuO[8] = uSelf[8]-1.0*uOther[8]; 
  uSMuO[9] = uSelf[9]-1.0*uOther[9]; 
  uSMuO[10] = uSelf[10]-1.0*uOther[10]; 
  uSMuO[11] = uSelf[11]-1.0*uOther[11]; 
  uSMuO[12] = uSelf[12]-1.0*uOther[12]; 
  uSMuO[13] = uSelf[13]-1.0*uOther[13]; 
  uSMuO[14] = uSelf[14]-1.0*uOther[14]; 
  uSMuO[15] = uSelf[15]-1.0*uOther[15]; 
  uSMuO[16] = uSelf[16]-1.0*uOther[16]; 
  uSMuO[17] = uSelf[17]-1.0*uOther[17]; 
  uSMuO[18] = uSelf[18]-1.0*uOther[18]; 
  uSMuO[19] = uSelf[19]-1.0*uOther[19]; 
  uSMuO[20] = uSelf[20]-1.0*uOther[20]; 
  uSMuO[21] = uSelf[21]-1.0*uOther[21]; 
  uSMuO[22] = uSelf[22]-1.0*uOther[22]; 
  uSMuO[23] = uSelf[23]-1.0*uOther[23]; 
  uSMuO[24] = uSelf[24]-1.0*uOther[24]; 
  uSMuO[25] = uSelf[25]-1.0*uOther[25]; 
  uSMuO[26] = uSelf[26]-1.0*uOther[26]; 
  uSMuO[27] = uSelf[27]-1.0*uOther[27]; 
  uSMuO[28] = uSelf[28]-1.0*uOther[28]; 
  uSMuO[29] = uSelf[29]-1.0*uOther[29]; 
  uSMuO[30] = uSelf[30]-1.0*uOther[30]; 
  uSMuO[31] = uSelf[31]-1.0*uOther[31]; 
  uSMuO[32] = uSelf[32]-1.0*uOther[32]; 
  uSMuO[33] = uSelf[33]-1.0*uOther[33]; 
  uSMuO[34] = uSelf[34]-1.0*uOther[34]; 
  uSMuO[35] = uSelf[35]-1.0*uOther[35]; 
  uSMuO[36] = uSelf[36]-1.0*uOther[36]; 
  uSMuO[37] = uSelf[37]-1.0*uOther[37]; 
  uSMuO[38] = uSelf[38]-1.0*uOther[38]; 
  uSMuO[39] = uSelf[39]-1.0*uOther[39]; 
  uSMuO[40] = uSelf[40]-1.0*uOther[40]; 
  uSMuO[41] = uSelf[41]-1.0*uOther[41]; 
  uSMuO[42] = uSelf[42]-1.0*uOther[42]; 
  uSMuO[43] = uSelf[43]-1.0*uOther[43]; 
  uSMuO[44] = uSelf[44]-1.0*uOther[44]; 
  uSMuO[45] = uSelf[45]-1.0*uOther[45]; 
  uSMuO[46] = uSelf[46]-1.0*uOther[46]; 
  uSMuO[47] = uSelf[47]-1.0*uOther[47]; 
  uSMuO[48] = uSelf[48]-1.0*uOther[48]; 
  uSMuO[49] = uSelf[49]-1.0*uOther[49]; 
  uSMuO[50] = uSelf[50]-1.0*uOther[50]; 
  uSMuO[51] = uSelf[51]-1.0*uOther[51]; 
  uSMuO[52] = uSelf[52]-1.0*uOther[52]; 
  uSMuO[53] = uSelf[53]-1.0*uOther[53]; 
  uSMuO[54] = uSelf[54]-1.0*uOther[54]; 
  uSMuO[55] = uSelf[55]-1.0*uOther[55]; 
  uSMuO[56] = uSelf[56]-1.0*uOther[56]; 
  uSMuO[57] = uSelf[57]-1.0*uOther[57]; 
  uSMuO[58] = uSelf[58]-1.0*uOther[58]; 
  uSMuO[59] = uSelf[59]-1.0*uOther[59]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[20]; 
  for (unsigned short int k=0; k<20; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 20*vd; 
    uSMuOSq[0] += 0.3535533905932737*uSMuO[a0+19]*uSMuO[a0+19]+0.3535533905932737*uSMuO[a0+18]*uSMuO[a0+18]+0.3535533905932737*uSMuO[a0+17]*uSMuO[a0+17]+0.3535533905932737*uSMuO[a0+16]*uSMuO[a0+16]+0.3535533905932737*uSMuO[a0+15]*uSMuO[a0+15]+0.3535533905932737*uSMuO[a0+14]*uSMuO[a0+14]+0.3535533905932737*uSMuO[a0+13]*uSMuO[a0+13]+0.3535533905932737*uSMuO[a0+12]*uSMuO[a0+12]+0.3535533905932737*uSMuO[a0+11]*uSMuO[a0+11]+0.3535533905932737*uSMuO[a0+10]*uSMuO[a0+10]+0.3535533905932737*uSMuO[a0+9]*uSMuO[a0+9]+0.3535533905932737*uSMuO[a0+8]*uSMuO[a0+8]+0.3535533905932737*uSMuO[a0+7]*uSMuO[a0+7]+0.3535533905932737*uSMuO[a0+6]*uSMuO[a0+6]+0.3535533905932737*uSMuO[a0+5]*uSMuO[a0+5]+0.3535533905932737*uSMuO[a0+4]*uSMuO[a0+4]+0.3535533905932737*uSMuO[a0+3]*uSMuO[a0+3]+0.3535533905932737*uSMuO[a0+2]*uSMuO[a0+2]+0.3535533905932737*uSMuO[a0+1]*uSMuO[a0+1]+0.3535533905932737*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += 0.7071067811865475*uSMuO[a0+16]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0+14]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+10]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+9]*uSMuO[a0+15]+0.632455532033676*uSMuO[a0+5]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0+8]*uSMuO[a0+12]+0.632455532033676*uSMuO[a0+4]*uSMuO[a0+11]+0.7071067811865475*uSMuO[a0+6]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+1]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+4]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+1]; 
    uSMuOSq[2] += 0.7071067811865475*uSMuO[a0+15]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+10]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0+13]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+9]*uSMuO[a0+16]+0.632455532033676*uSMuO[a0+6]*uSMuO[a0+14]+0.632455532033676*uSMuO[a0+4]*uSMuO[a0+12]+0.7071067811865475*uSMuO[a0+7]*uSMuO[a0+11]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+2]*uSMuO[a0+8]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+4]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+2]; 
    uSMuOSq[3] += 0.6324555320336759*uSMuO[a0+10]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0+12]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0+11]*uSMuO[a0+17]+0.632455532033676*uSMuO[a0+6]*uSMuO[a0+16]+0.632455532033676*uSMuO[a0+5]*uSMuO[a0+15]+0.7071067811865475*uSMuO[a0+8]*uSMuO[a0+14]+0.7071067811865475*uSMuO[a0+7]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+3]*uSMuO[a0+9]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+3]; 
    uSMuOSq[4] += 0.7071067811865475*uSMuO[a0+9]*uSMuO[a0+19]+0.5656854249492381*uSMuO[a0+17]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+6]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+5]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+15]*uSMuO[a0+16]+0.632455532033676*uSMuO[a0+10]*uSMuO[a0+14]+0.632455532033676*uSMuO[a0+10]*uSMuO[a0+13]+0.5656854249492381*uSMuO[a0+11]*uSMuO[a0+12]+0.632455532033676*uSMuO[a0+2]*uSMuO[a0+12]+0.632455532033676*uSMuO[a0+1]*uSMuO[a0+11]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+4]*uSMuO[a0+8]+0.6324555320336759*uSMuO[a0+4]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+4]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+2]; 
    uSMuOSq[5] += 0.5656854249492381*uSMuO[a0+17]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+6]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0+8]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+4]*uSMuO[a0+17]+0.632455532033676*uSMuO[a0+10]*uSMuO[a0+16]+0.5656854249492381*uSMuO[a0+13]*uSMuO[a0+15]+0.632455532033676*uSMuO[a0+3]*uSMuO[a0+15]+0.7071067811865475*uSMuO[a0+12]*uSMuO[a0+14]+0.632455532033676*uSMuO[a0+1]*uSMuO[a0+13]+0.632455532033676*uSMuO[a0+10]*uSMuO[a0+11]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+5]*uSMuO[a0+9]+0.6324555320336759*uSMuO[a0+5]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+3]; 
    uSMuOSq[6] += 0.5656854249492381*uSMuO[a0+18]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+5]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+4]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0+7]*uSMuO[a0+17]+0.5656854249492381*uSMuO[a0+14]*uSMuO[a0+16]+0.632455532033676*uSMuO[a0+3]*uSMuO[a0+16]+0.632455532033676*uSMuO[a0+10]*uSMuO[a0+15]+0.632455532033676*uSMuO[a0+2]*uSMuO[a0+14]+0.7071067811865475*uSMuO[a0+11]*uSMuO[a0+13]+0.632455532033676*uSMuO[a0+10]*uSMuO[a0+12]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+6]*uSMuO[a0+9]+0.6324555320336759*uSMuO[a0+6]*uSMuO[a0+8]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+3]; 
    uSMuOSq[7] += 0.3162277660168379*uSMuO[a0+19]*uSMuO[a0+19]+0.3162277660168379*uSMuO[a0+18]*uSMuO[a0+18]+0.2258769757263128*uSMuO[a0+17]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+6]*uSMuO[a0+17]+0.3162277660168379*uSMuO[a0+15]*uSMuO[a0+15]+0.2258769757263128*uSMuO[a0+13]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+13]+0.3162277660168379*uSMuO[a0+12]*uSMuO[a0+12]+0.2258769757263128*uSMuO[a0+11]*uSMuO[a0+11]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+11]+0.3162277660168379*uSMuO[a0+10]*uSMuO[a0+10]+0.2258769757263128*uSMuO[a0+7]*uSMuO[a0+7]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+7]+0.3162277660168379*uSMuO[a0+5]*uSMuO[a0+5]+0.3162277660168379*uSMuO[a0+4]*uSMuO[a0+4]+0.3162277660168379*uSMuO[a0+1]*uSMuO[a0+1]; 
    uSMuOSq[8] += 0.3162277660168379*uSMuO[a0+19]*uSMuO[a0+19]+0.2258769757263128*uSMuO[a0+18]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+18]+0.3162277660168379*uSMuO[a0+17]*uSMuO[a0+17]+0.3162277660168379*uSMuO[a0+16]*uSMuO[a0+16]+0.2258769757263128*uSMuO[a0+14]*uSMuO[a0+14]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+14]+0.2258769757263128*uSMuO[a0+12]*uSMuO[a0+12]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+12]+0.3162277660168379*uSMuO[a0+11]*uSMuO[a0+11]+0.3162277660168379*uSMuO[a0+10]*uSMuO[a0+10]+0.2258769757263128*uSMuO[a0+8]*uSMuO[a0+8]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+8]+0.3162277660168379*uSMuO[a0+6]*uSMuO[a0+6]+0.3162277660168379*uSMuO[a0+4]*uSMuO[a0+4]+0.3162277660168379*uSMuO[a0+2]*uSMuO[a0+2]; 
    uSMuOSq[9] += 0.2258769757263128*uSMuO[a0+19]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+19]+0.3162277660168379*uSMuO[a0+18]*uSMuO[a0+18]+0.3162277660168379*uSMuO[a0+17]*uSMuO[a0+17]+0.2258769757263128*uSMuO[a0+16]*uSMuO[a0+16]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+16]+0.2258769757263128*uSMuO[a0+15]*uSMuO[a0+15]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+15]+0.3162277660168379*uSMuO[a0+14]*uSMuO[a0+14]+0.3162277660168379*uSMuO[a0+13]*uSMuO[a0+13]+0.3162277660168379*uSMuO[a0+10]*uSMuO[a0+10]+0.2258769757263128*uSMuO[a0+9]*uSMuO[a0+9]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+9]+0.3162277660168379*uSMuO[a0+6]*uSMuO[a0+6]+0.3162277660168379*uSMuO[a0+5]*uSMuO[a0+5]+0.3162277660168379*uSMuO[a0+3]*uSMuO[a0+3]; 
    uSMuOSq[10] += 0.5656854249492382*uSMuO[a0+14]*uSMuO[a0+19]+0.5656854249492382*uSMuO[a0+13]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+3]*uSMuO[a0+19]+0.5656854249492382*uSMuO[a0+16]*uSMuO[a0+18]+0.5656854249492382*uSMuO[a0+11]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+2]*uSMuO[a0+18]+0.5656854249492382*uSMuO[a0+15]*uSMuO[a0+17]+0.5656854249492382*uSMuO[a0+12]*uSMuO[a0+17]+0.6324555320336759*uSMuO[a0+1]*uSMuO[a0+17]+0.632455532033676*uSMuO[a0+5]*uSMuO[a0+16]+0.632455532033676*uSMuO[a0+6]*uSMuO[a0+15]+0.632455532033676*uSMuO[a0+4]*uSMuO[a0+14]+0.632455532033676*uSMuO[a0+4]*uSMuO[a0+13]+0.632455532033676*uSMuO[a0+6]*uSMuO[a0+12]+0.632455532033676*uSMuO[a0+5]*uSMuO[a0+11]+0.6324555320336759*uSMuO[a0+9]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+8]*uSMuO[a0+10]+0.6324555320336759*uSMuO[a0+7]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+6]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+5]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+4]; 
    uSMuOSq[11] += 0.6324555320336759*uSMuO[a0+15]*uSMuO[a0+19]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+14]*uSMuO[a0+17]+0.4517539514526256*uSMuO[a0+13]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+6]*uSMuO[a0+13]+0.5656854249492381*uSMuO[a0+4]*uSMuO[a0+12]+0.6324555320336759*uSMuO[a0+8]*uSMuO[a0+11]+0.4517539514526256*uSMuO[a0+7]*uSMuO[a0+11]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+11]+0.632455532033676*uSMuO[a0+5]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+7]+0.632455532033676*uSMuO[a0+1]*uSMuO[a0+4]; 
    uSMuOSq[12] += 0.6324555320336759*uSMuO[a0+16]*uSMuO[a0+19]+0.4517539514526256*uSMuO[a0+14]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+13]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+18]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+14]+0.4517539514526256*uSMuO[a0+8]*uSMuO[a0+12]+0.6324555320336759*uSMuO[a0+7]*uSMuO[a0+12]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+12]+0.5656854249492381*uSMuO[a0+4]*uSMuO[a0+11]+0.632455532033676*uSMuO[a0+6]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+8]+0.632455532033676*uSMuO[a0+2]*uSMuO[a0+4]; 
    uSMuOSq[13] += 0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+12]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+16]*uSMuO[a0+17]+0.4517539514526256*uSMuO[a0+11]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+17]+0.5656854249492381*uSMuO[a0+5]*uSMuO[a0+15]+0.6324555320336759*uSMuO[a0+9]*uSMuO[a0+13]+0.4517539514526256*uSMuO[a0+7]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0+6]*uSMuO[a0+11]+0.632455532033676*uSMuO[a0+4]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+7]+0.632455532033676*uSMuO[a0+1]*uSMuO[a0+5]; 
    uSMuOSq[14] += 0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+15]*uSMuO[a0+18]+0.4517539514526256*uSMuO[a0+12]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+11]*uSMuO[a0+17]+0.5656854249492381*uSMuO[a0+6]*uSMuO[a0+16]+0.6324555320336759*uSMuO[a0+9]*uSMuO[a0+14]+0.4517539514526256*uSMuO[a0+8]*uSMuO[a0+14]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+14]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+12]+0.632455532033676*uSMuO[a0+4]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+8]+0.632455532033676*uSMuO[a0+2]*uSMuO[a0+6]; 
    uSMuOSq[15] += 0.4517539514526256*uSMuO[a0+16]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+11]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+14]*uSMuO[a0+18]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+16]+0.4517539514526256*uSMuO[a0+9]*uSMuO[a0+15]+0.6324555320336759*uSMuO[a0+7]*uSMuO[a0+15]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+15]+0.5656854249492381*uSMuO[a0+5]*uSMuO[a0+13]+0.632455532033676*uSMuO[a0+6]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+9]+0.632455532033676*uSMuO[a0+3]*uSMuO[a0+5]; 
    uSMuOSq[16] += 0.4517539514526256*uSMuO[a0+15]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+12]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+19]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+13]*uSMuO[a0+17]+0.4517539514526256*uSMuO[a0+9]*uSMuO[a0+16]+0.6324555320336759*uSMuO[a0+8]*uSMuO[a0+16]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+16]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+15]+0.5656854249492381*uSMuO[a0+6]*uSMuO[a0+14]+0.632455532033676*uSMuO[a0+5]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+9]+0.632455532033676*uSMuO[a0+3]*uSMuO[a0+6]; 
    uSMuOSq[17] += 0.5059644256269408*uSMuO[a0+18]*uSMuO[a0+19]+0.5656854249492381*uSMuO[a0+5]*uSMuO[a0+19]+0.5656854249492381*uSMuO[a0+4]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+9]*uSMuO[a0+17]+0.6324555320336759*uSMuO[a0+8]*uSMuO[a0+17]+0.4517539514526256*uSMuO[a0+7]*uSMuO[a0+17]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+17]+0.6324555320336759*uSMuO[a0+13]*uSMuO[a0+16]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+15]+0.6324555320336759*uSMuO[a0+11]*uSMuO[a0+14]+0.4517539514526256*uSMuO[a0+11]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+13]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+12]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+11]+0.6324555320336759*uSMuO[a0+1]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+6]*uSMuO[a0+7]+0.6324555320336759*uSMuO[a0+4]*uSMuO[a0+5]; 
    uSMuOSq[18] += 0.5059644256269408*uSMuO[a0+17]*uSMuO[a0+19]+0.5656854249492381*uSMuO[a0+6]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+9]*uSMuO[a0+18]+0.4517539514526256*uSMuO[a0+8]*uSMuO[a0+18]+0.6324555320336759*uSMuO[a0+7]*uSMuO[a0+18]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+18]+0.5656854249492381*uSMuO[a0+4]*uSMuO[a0+17]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+16]+0.6324555320336759*uSMuO[a0+14]*uSMuO[a0+15]+0.4517539514526256*uSMuO[a0+12]*uSMuO[a0+14]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+14]+0.6324555320336759*uSMuO[a0+12]*uSMuO[a0+13]+0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+12]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+11]+0.6324555320336759*uSMuO[a0+2]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+5]*uSMuO[a0+8]+0.6324555320336759*uSMuO[a0+4]*uSMuO[a0+6]; 
    uSMuOSq[19] += 0.4517539514526256*uSMuO[a0+9]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+8]*uSMuO[a0+19]+0.6324555320336759*uSMuO[a0+7]*uSMuO[a0+19]+0.7071067811865475*uSMuO[a0]*uSMuO[a0+19]+0.5059644256269408*uSMuO[a0+17]*uSMuO[a0+18]+0.5656854249492381*uSMuO[a0+6]*uSMuO[a0+18]+0.5656854249492381*uSMuO[a0+5]*uSMuO[a0+17]+0.4517539514526256*uSMuO[a0+15]*uSMuO[a0+16]+0.6324555320336759*uSMuO[a0+12]*uSMuO[a0+16]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+16]+0.6324555320336759*uSMuO[a0+11]*uSMuO[a0+15]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+15]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+14]+0.5656854249492382*uSMuO[a0+10]*uSMuO[a0+13]+0.6324555320336759*uSMuO[a0+3]*uSMuO[a0+10]+0.7071067811865475*uSMuO[a0+4]*uSMuO[a0+9]+0.6324555320336759*uSMuO[a0+5]*uSMuO[a0+6]; 
  } 
 
  // ..... Get the cross flow velocity uCross ..... // 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 20*vd; 
    uCross[a0] = uOther[a0]; 
    uCross[a0+1] = uOther[a0+1]; 
    uCross[a0+2] = uOther[a0+2]; 
    uCross[a0+3] = uOther[a0+3]; 
    uCross[a0+4] = uOther[a0+4]; 
    uCross[a0+5] = uOther[a0+5]; 
    uCross[a0+6] = uOther[a0+6]; 
    uCross[a0+7] = uOther[a0+7]; 
    uCross[a0+8] = uOther[a0+8]; 
    uCross[a0+9] = uOther[a0+9]; 
    uCross[a0+10] = uOther[a0+10]; 
    uCross[a0+11] = uOther[a0+11]; 
    uCross[a0+12] = uOther[a0+12]; 
    uCross[a0+13] = uOther[a0+13]; 
    uCross[a0+14] = uOther[a0+14]; 
    uCross[a0+15] = uOther[a0+15]; 
    uCross[a0+16] = uOther[a0+16]; 
    uCross[a0+17] = uOther[a0+17]; 
    uCross[a0+18] = uOther[a0+18]; 
    uCross[a0+19] = uOther[a0+19]; 
 
  } 
 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = 0.3333333333333333*uSMuOSq[0]*mRat+vtSqSelf[0]; 
  vtSqCross[1] = 0.3333333333333333*uSMuOSq[1]*mRat+vtSqSelf[1]; 
  vtSqCross[2] = 0.3333333333333333*uSMuOSq[2]*mRat+vtSqSelf[2]; 
  vtSqCross[3] = 0.3333333333333333*uSMuOSq[3]*mRat+vtSqSelf[3]; 
  vtSqCross[4] = 0.3333333333333333*uSMuOSq[4]*mRat+vtSqSelf[4]; 
  vtSqCross[5] = 0.3333333333333333*uSMuOSq[5]*mRat+vtSqSelf[5]; 
  vtSqCross[6] = 0.3333333333333333*uSMuOSq[6]*mRat+vtSqSelf[6]; 
  vtSqCross[7] = 0.3333333333333333*uSMuOSq[7]*mRat+vtSqSelf[7]; 
  vtSqCross[8] = 0.3333333333333333*uSMuOSq[8]*mRat+vtSqSelf[8]; 
  vtSqCross[9] = 0.3333333333333333*uSMuOSq[9]*mRat+vtSqSelf[9]; 
  vtSqCross[10] = 0.3333333333333333*uSMuOSq[10]*mRat+vtSqSelf[10]; 
  vtSqCross[11] = 0.3333333333333333*uSMuOSq[11]*mRat+vtSqSelf[11]; 
  vtSqCross[12] = 0.3333333333333333*uSMuOSq[12]*mRat+vtSqSelf[12]; 
  vtSqCross[13] = 0.3333333333333333*uSMuOSq[13]*mRat+vtSqSelf[13]; 
  vtSqCross[14] = 0.3333333333333333*uSMuOSq[14]*mRat+vtSqSelf[14]; 
  vtSqCross[15] = 0.3333333333333333*uSMuOSq[15]*mRat+vtSqSelf[15]; 
  vtSqCross[16] = 0.3333333333333333*uSMuOSq[16]*mRat+vtSqSelf[16]; 
  vtSqCross[17] = 0.3333333333333333*uSMuOSq[17]*mRat+vtSqSelf[17]; 
  vtSqCross[18] = 0.3333333333333333*uSMuOSq[18]*mRat+vtSqSelf[18]; 
  vtSqCross[19] = 0.3333333333333333*uSMuOSq[19]*mRat+vtSqSelf[19]; 
 
} 
 
