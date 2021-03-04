#include <ProjectFluxOnGhosts.h> 
#include <math.h> 

void ProjectFluxOnGhosts1x1vDir1Ser_P1(const double zVal, const double *fIn, double *fHat) 
{ 
  // fIn[4]:    input distf. 
  // fHat[4]:   projection of flux on ghost cell. 
 
  fHat[0] = 0.3333333333333333*(3.0*fIn[4]*zVal+1.732050807568877*fIn[3]); 
  fHat[1] = zVal*(1.732050807568877*fIn[4]*zVal+fIn[3]); 
  fHat[2] = 0.3333333333333333*(3.0*fIn[2]*zVal+1.732050807568877*fIn[1]); 
  fHat[3] = zVal*(1.732050807568877*fIn[2]*zVal+fIn[1]); 
 
}

void ProjectFluxOnGhosts1x2vDir1Ser_P1(const double zVal, const double *fIn, double *fHat) 
{ 
  // fIn[8]:    input distf. 
  // fHat[8]:   projection of flux on ghost cell. 
 
  fHat[0] = 0.3333333333333333*(3.0*fIn[5]*zVal+1.732050807568877*fIn[3]); 
  fHat[1] = zVal*(1.732050807568877*fIn[5]*zVal+fIn[3]); 
  fHat[2] = 0.3333333333333333*(3.0*fIn[2]*zVal+1.732050807568877*fIn[1]); 
  fHat[3] = 0.3333333333333333*(3.0*fIn[8]*zVal+1.732050807568877*fIn[7]); 
  fHat[4] = zVal*(1.732050807568877*fIn[2]*zVal+fIn[1]); 
  fHat[5] = zVal*(1.732050807568877*fIn[8]*zVal+fIn[7]); 
  fHat[6] = 0.3333333333333333*(3.0*fIn[6]*zVal+1.732050807568877*fIn[4]); 
  fHat[7] = zVal*(1.732050807568877*fIn[6]*zVal+fIn[4]); 
 
}

void ProjectFluxOnGhosts1x3vDir1Ser_P1(const double zVal, const double *fIn, double *fHat) 
{ 
  // fIn[16]:    input distf. 
  // fHat[16]:   projection of flux on ghost cell. 
 
  fHat[0] = 0.3333333333333333*(3.0*fIn[6]*zVal+1.732050807568877*fIn[3]); 
  fHat[1] = zVal*(1.732050807568877*fIn[6]*zVal+fIn[3]); 
  fHat[2] = 0.3333333333333333*(3.0*fIn[2]*zVal+1.732050807568877*fIn[1]); 
  fHat[3] = 0.3333333333333333*(3.0*fIn[12]*zVal+1.732050807568877*fIn[8]); 
  fHat[4] = 0.3333333333333333*(3.0*fIn[13]*zVal+1.732050807568877*fIn[10]); 
  fHat[5] = zVal*(1.732050807568877*fIn[2]*zVal+fIn[1]); 
  fHat[6] = zVal*(1.732050807568877*fIn[12]*zVal+fIn[8]); 
  fHat[7] = 0.3333333333333333*(3.0*fIn[7]*zVal+1.732050807568877*fIn[4]); 
  fHat[8] = zVal*(1.732050807568877*fIn[13]*zVal+fIn[10]); 
  fHat[9] = 0.3333333333333333*(3.0*fIn[9]*zVal+1.732050807568877*fIn[5]); 
  fHat[10] = 0.3333333333333333*(3.0*fIn[16]*zVal+1.732050807568877*fIn[15]); 
  fHat[11] = zVal*(1.732050807568877*fIn[7]*zVal+fIn[4]); 
  fHat[12] = zVal*(1.732050807568877*fIn[9]*zVal+fIn[5]); 
  fHat[13] = zVal*(1.732050807568877*fIn[16]*zVal+fIn[15]); 
  fHat[14] = 0.3333333333333333*(3.0*fIn[14]*zVal+1.732050807568877*fIn[11]); 
  fHat[15] = zVal*(1.732050807568877*fIn[14]*zVal+fIn[11]); 
 
}

void ProjectFluxOnGhosts3x3vDir3Ser_P1(const double zVal, const double *fIn, double *fHat) 
{ 
  // fIn[64]:    input distf. 
  // fHat[64]:   projection of flux on ghost cell. 
 
  fHat[0] = 0.3333333333333333*(3.0*fIn[20]*zVal+1.732050807568877*fIn[7]); 
  fHat[1] = 0.3333333333333333*(3.0*fIn[34]*zVal+1.732050807568877*fIn[18]); 
  fHat[2] = 0.3333333333333333*(3.0*fIn[35]*zVal+1.732050807568877*fIn[19]); 
  fHat[3] = zVal*(1.732050807568877*fIn[20]*zVal+fIn[7]); 
  fHat[4] = 0.3333333333333333*(3.0*fIn[38]*zVal+1.732050807568877*fIn[21]); 
  fHat[5] = 0.3333333333333333*(3.0*fIn[41]*zVal+1.732050807568877*fIn[22]); 
  fHat[6] = 0.3333333333333333*(3.0*fIn[4]*zVal+1.732050807568877*fIn[1]); 
  fHat[7] = 0.3333333333333333*(3.0*fIn[48]*zVal+1.732050807568877*fIn[33]); 
  fHat[8] = zVal*(1.732050807568877*fIn[34]*zVal+fIn[18]); 
  fHat[9] = zVal*(1.732050807568877*fIn[35]*zVal+fIn[19]); 
  fHat[10] = 0.3333333333333333*(3.0*fIn[50]*zVal+1.732050807568877*fIn[36]); 
  fHat[11] = 0.3333333333333333*(3.0*fIn[51]*zVal+1.732050807568877*fIn[37]); 
  fHat[12] = zVal*(1.732050807568877*fIn[38]*zVal+fIn[21]); 
  fHat[13] = 0.3333333333333333*(3.0*fIn[53]*zVal+1.732050807568877*fIn[39]); 
  fHat[14] = 0.3333333333333333*(3.0*fIn[54]*zVal+1.732050807568877*fIn[40]); 
  fHat[15] = zVal*(1.732050807568877*fIn[41]*zVal+fIn[22]); 
  fHat[16] = 0.3333333333333333*(3.0*fIn[57]*zVal+1.732050807568877*fIn[42]); 
  fHat[17] = 0.3333333333333333*(3.0*fIn[9]*zVal+1.732050807568877*fIn[2]); 
  fHat[18] = 0.3333333333333333*(3.0*fIn[10]*zVal+1.732050807568877*fIn[3]); 
  fHat[19] = zVal*(1.732050807568877*fIn[4]*zVal+fIn[1]); 
  fHat[20] = 0.3333333333333333*(3.0*fIn[13]*zVal+1.732050807568877*fIn[5]); 
  fHat[21] = 0.3333333333333333*(3.0*fIn[16]*zVal+1.732050807568877*fIn[6]); 
  fHat[22] = zVal*(1.732050807568877*fIn[48]*zVal+fIn[33]); 
  fHat[23] = 0.3333333333333333*(3.0*fIn[59]*zVal+1.732050807568877*fIn[49]); 
  fHat[24] = zVal*(1.732050807568877*fIn[50]*zVal+fIn[36]); 
  fHat[25] = zVal*(1.732050807568877*fIn[51]*zVal+fIn[37]); 
  fHat[26] = 0.3333333333333333*(3.0*fIn[60]*zVal+1.732050807568877*fIn[52]); 
  fHat[27] = zVal*(1.732050807568877*fIn[53]*zVal+fIn[39]); 
  fHat[28] = zVal*(1.732050807568877*fIn[54]*zVal+fIn[40]); 
  fHat[29] = 0.3333333333333333*(3.0*fIn[62]*zVal+1.732050807568877*fIn[55]); 
  fHat[30] = 0.3333333333333333*(3.0*fIn[63]*zVal+1.732050807568877*fIn[56]); 
  fHat[31] = zVal*(1.732050807568877*fIn[57]*zVal+fIn[42]); 
  fHat[32] = 0.3333333333333333*(3.0*fIn[23]*zVal+1.732050807568877*fIn[8]); 
  fHat[33] = zVal*(1.732050807568877*fIn[9]*zVal+fIn[2]); 
  fHat[34] = zVal*(1.732050807568877*fIn[10]*zVal+fIn[3]); 
  fHat[35] = 0.3333333333333333*(3.0*fIn[25]*zVal+1.732050807568877*fIn[11]); 
  fHat[36] = 0.3333333333333333*(3.0*fIn[26]*zVal+1.732050807568877*fIn[12]); 
  fHat[37] = zVal*(1.732050807568877*fIn[13]*zVal+fIn[5]); 
  fHat[38] = 0.3333333333333333*(3.0*fIn[28]*zVal+1.732050807568877*fIn[14]); 
  fHat[39] = 0.3333333333333333*(3.0*fIn[29]*zVal+1.732050807568877*fIn[15]); 
  fHat[40] = zVal*(1.732050807568877*fIn[16]*zVal+fIn[6]); 
  fHat[41] = 0.3333333333333333*(3.0*fIn[32]*zVal+1.732050807568877*fIn[17]); 
  fHat[42] = zVal*(1.732050807568877*fIn[59]*zVal+fIn[49]); 
  fHat[43] = zVal*(1.732050807568877*fIn[60]*zVal+fIn[52]); 
  fHat[44] = 0.3333333333333333*(3.0*fIn[64]*zVal+1.732050807568877*fIn[61]); 
  fHat[45] = zVal*(1.732050807568877*fIn[62]*zVal+fIn[55]); 
  fHat[46] = zVal*(1.732050807568877*fIn[63]*zVal+fIn[56]); 
  fHat[47] = zVal*(1.732050807568877*fIn[23]*zVal+fIn[8]); 
  fHat[48] = 0.3333333333333333*(3.0*fIn[43]*zVal+1.732050807568877*fIn[24]); 
  fHat[49] = zVal*(1.732050807568877*fIn[25]*zVal+fIn[11]); 
  fHat[50] = zVal*(1.732050807568877*fIn[26]*zVal+fIn[12]); 
  fHat[51] = 0.3333333333333333*(3.0*fIn[44]*zVal+1.732050807568877*fIn[27]); 
  fHat[52] = zVal*(1.732050807568877*fIn[28]*zVal+fIn[14]); 
  fHat[53] = zVal*(1.732050807568877*fIn[29]*zVal+fIn[15]); 
  fHat[54] = 0.3333333333333333*(3.0*fIn[46]*zVal+1.732050807568877*fIn[30]); 
  fHat[55] = 0.3333333333333333*(3.0*fIn[47]*zVal+1.732050807568877*fIn[31]); 
  fHat[56] = zVal*(1.732050807568877*fIn[32]*zVal+fIn[17]); 
  fHat[57] = zVal*(1.732050807568877*fIn[64]*zVal+fIn[61]); 
  fHat[58] = zVal*(1.732050807568877*fIn[43]*zVal+fIn[24]); 
  fHat[59] = zVal*(1.732050807568877*fIn[44]*zVal+fIn[27]); 
  fHat[60] = 0.3333333333333333*(3.0*fIn[58]*zVal+1.732050807568877*fIn[45]); 
  fHat[61] = zVal*(1.732050807568877*fIn[46]*zVal+fIn[30]); 
  fHat[62] = zVal*(1.732050807568877*fIn[47]*zVal+fIn[31]); 
  fHat[63] = zVal*(1.732050807568877*fIn[58]*zVal+fIn[45]); 
 
}
