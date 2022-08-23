#include <CartFieldIntegratedQuantCalcImpl.h>
#include <cmath>

void
gkylCartFieldIntQuantV(int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out)
{
  double vol = 1.0;
  if (nb>1)
  {
    for (unsigned d=0; d<ndim; ++d)
      vol *= dxv[d]/(std::sqrt(2.0));
  }
  else
  { // for polyOrder = 0 no normalization is applied
    for (unsigned d=0; d<ndim; ++d)
      vol *= dxv[d];
  }

  for (unsigned c=0; c<nc; ++c)
    out[c] += fIn[c*nb]*vol;
}

void
gkylCartFieldIntQuantAbsV(int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out)
{
  double vol = 1.0;
  if (nb>1)
  {
    for (unsigned d=0; d<ndim; ++d)
      vol *= dxv[d]/(std::sqrt(2.0));
  }
  else
  { // for polyOrder = 0 no normalization is applied
    for (unsigned d=0; d<ndim; ++d)
      vol *= dxv[d];
  }

  for (unsigned c=0; c<nc; ++c)
    out[c] += std::abs(fIn[c*nb])*vol;
}

void
gkylCartFieldIntQuantV2(int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out)
{
  double vol = 1.0;
  if (nb>1)
  {
    for (unsigned d=0; d<ndim; ++d)
      vol *= dxv[d]/2.0;
  }
  else
  { // for polyOrder = 0 no normalization is applied
    for (unsigned d=0; d<ndim; ++d)
      vol *= dxv[d];
  }

  for (unsigned c=0; c<nc; ++c)
  {
    double v2 = 0.0;
    for (unsigned b=0; b<nb; ++b)
      v2 += fIn[c*nb+b]*fIn[c*nb+b];
    out[c] += v2*vol;
  }
}

void
gkylCartFieldIntQuantGradPerpV2_2x_p1(int ndim, unsigned nc, unsigned nb, const double *dxv, const double *inw, const double *fIn, double *out)
{
  double vol = 1.0;
  // assume polyOrder >= 1
  for (unsigned d=0; d<2; ++d)
    vol *= dxv[d]/2.0;

  double dfac2_x = 4.0/dxv[0]/dxv[0];
  double dfac2_y = 4.0/dxv[1]/dxv[1];

  double phixSq[4] = {0.};
  phixSq[0] = 0.5*(3.0*fIn[3]*fIn[3]+3.0*fIn[1]*fIn[1]);
  phixSq[2] = 3.0*fIn[1]*fIn[3];

  double phiySq[4] = {0.};
  phiySq[0] = 0.5*(3.0*fIn[3]*fIn[3]+3.0*fIn[2]*fIn[2]);
  phiySq[1] = 3.0*fIn[2]*fIn[3];

  // assume 1 component
  out[0] += vol*((inw[1]*phiySq[1]+inw[0]*phiySq[0])*dfac2_y+(inw[2]*phixSq[2]+inw[0]*phixSq[0])*dfac2_x);
}

void
gkylCartFieldIntQuantGradPerpV2_2x_p2(int ndim, unsigned nc, unsigned nb, const double *dxv, const double *inw, const double *fIn, double *out)
{
  double vol = 1.0;
  // assume polyOrder >= 1
  for (unsigned d=0; d<2; ++d)
    vol *= dxv[d]/2.0;

  double dfac2_x = 4.0/dxv[0]/dxv[0];
  double dfac2_y = 4.0/dxv[1]/dxv[1];

  double phixSq[8] = {0.};
  phixSq[0] = 0.5*(3.0*fIn[7]*fIn[7]+15.0*fIn[6]*fIn[6]+15.0*fIn[4]*fIn[4]+3.0*fIn[3]*fIn[3]+3.0*fIn[1]*fIn[1]);
  phixSq[1] = 6.708203932499369*fIn[3]*fIn[6]+6.708203932499369*fIn[1]*fIn[4];
  phixSq[2] = 0.2*(13.41640786499874*fIn[3]*fIn[7]+75.00000000000001*fIn[4]*fIn[6]+15.0*fIn[1]*fIn[3]);
  phixSq[3] = 6.0*fIn[6]*fIn[7]+6.708203932499369*fIn[1]*fIn[6]+6.708203932499369*fIn[3]*fIn[4];
  phixSq[4] = 6.708203932499369*fIn[6]*fIn[6]+6.708203932499369*fIn[4]*fIn[4];
  phixSq[5] = 0.02857142857142857*(33.54101966249685*fIn[7]*fIn[7]+105.0*fIn[1]*fIn[7]+234.787137637478*fIn[6]*fIn[6]+46.95742752749558*fIn[3]*fIn[3]);
  phixSq[6] = 13.41640786499874*fIn[4]*fIn[6];
  phixSq[7] = 6.708203932499369*fIn[4]*fIn[7]+6.0*fIn[3]*fIn[6];

  double phiySq[8] = {0.};
  phiySq[0] = 0.5*(15.0*fIn[7]*fIn[7]+3.0*fIn[6]*fIn[6]+15.0*fIn[5]*fIn[5]+3.0*fIn[3]*fIn[3]+3.0*fIn[2]*fIn[2]);
  phiySq[1] = 0.2*(75.00000000000001*fIn[5]*fIn[7]+13.41640786499874*fIn[3]*fIn[6]+15.0*fIn[2]*fIn[3]);
  phiySq[2] = 6.708203932499369*fIn[3]*fIn[7]+6.708203932499369*fIn[2]*fIn[5];
  phiySq[3] = (6.0*fIn[6]+6.708203932499369*fIn[2])*fIn[7]+6.708203932499369*fIn[3]*fIn[5];
  phiySq[4] = 0.02857142857142857*(234.787137637478*fIn[7]*fIn[7]+33.54101966249685*fIn[6]*fIn[6]+105.0*fIn[2]*fIn[6]+46.95742752749558*fIn[3]*fIn[3]);
  phiySq[5] = 6.708203932499369*fIn[7]*fIn[7]+6.708203932499369*fIn[5]*fIn[5];
  phiySq[6] = 6.0*fIn[3]*fIn[7]+6.708203932499369*fIn[5]*fIn[6];
  phiySq[7] = 13.41640786499874*fIn[5]*fIn[7];

  // assume 1 component
  out[0] += vol*((inw[7]*phiySq[7]+inw[6]*phiySq[6]+inw[5]*phiySq[5]+inw[4]*phiySq[4]+inw[3]*phiySq[3]+inw[2]*phiySq[2]+inw[1]*phiySq[1]+inw[0]*phiySq[0])*dfac2_y+(inw[7]*phixSq[7]+inw[6]*phixSq[6]+inw[5]*phixSq[5]+inw[4]*phixSq[4]+inw[3]*phixSq[3]+inw[2]*phixSq[2]+inw[1]*phixSq[1]+inw[0]*phixSq[0])*dfac2_x);
}

void
gkylCartFieldIntQuantGradPerpV2_3x_p1(int ndim, unsigned nc, unsigned nb, const double *dxv, const double *inw, const double *fIn, double *out)
{
  double vol = 1.0;
  // assume polyOrder >= 1
  for (unsigned d=0; d<3; ++d)
    vol *= dxv[d]/2.0;

  double dfac2_x = 4.0/dxv[0]/dxv[0];
  double dfac2_y = 4.0/dxv[1]/dxv[1];

  double phixSq[8] = {0.};
  phixSq[0] = 0.3535533905932737*(3.0*fIn[7]*fIn[7]+3.0*fIn[5]*fIn[5]+3.0*fIn[4]*fIn[4]+3.0*fIn[1]*fIn[1]);
  phixSq[2] = 0.7071067811865475*(3.0*fIn[5]*fIn[7]+3.0*fIn[1]*fIn[4]);
  phixSq[3] = 0.7071067811865475*(3.0*fIn[4]*fIn[7]+3.0*fIn[1]*fIn[5]);
  phixSq[6] = 0.7071067811865475*(3.0*fIn[1]*fIn[7]+3.0*fIn[4]*fIn[5]);

  double phiySq[8] = {0.};
  phiySq[0] = 0.3535533905932737*(3.0*fIn[7]*fIn[7]+3.0*fIn[6]*fIn[6]+3.0*fIn[4]*fIn[4]+3.0*fIn[2]*fIn[2]);
  phiySq[1] = 0.7071067811865475*(3.0*fIn[6]*fIn[7]+3.0*fIn[2]*fIn[4]);
  phiySq[3] = 0.7071067811865475*(3.0*fIn[4]*fIn[7]+3.0*fIn[2]*fIn[6]);
  phiySq[5] = 0.7071067811865475*(3.0*fIn[2]*fIn[7]+3.0*fIn[4]*fIn[6]);

  // assume 1 component
  out[0] += vol*((inw[5]*phiySq[5]+inw[3]*phiySq[3]+inw[1]*phiySq[1]+inw[0]*phiySq[0])*dfac2_y+(inw[6]*phixSq[6]+inw[3]*phixSq[3]+inw[2]*phixSq[2]+inw[0]*phixSq[0])*dfac2_x);
}

void
gkylCartFieldIntQuantGradPerpV2_3x_p2(int ndim, unsigned nc, unsigned nb, const double *dxv, const double *inw, const double *fIn, double *out)
{
  double vol = 1.0;
  // assume polyOrder >= 1
  for (unsigned d=0; d<3; ++d)
    vol *= dxv[d]/2.0;

  double dfac2_x = 4.0/dxv[0]/dxv[0];
  double dfac2_y = 4.0/dxv[1]/dxv[1];

  double phixSq[20] = {0.};
  phixSq[0] = 0.3535533905932737*(3.0*fIn[19]*fIn[19]+3.0*fIn[18]*fIn[18]+15.0*fIn[17]*fIn[17]+3.0*fIn[15]*fIn[15]+15.0*fIn[13]*fIn[13]+3.0*fIn[12]*fIn[12]+15.0*fIn[11]*fIn[11]+3.0*fIn[10]*fIn[10]+15.0*fIn[7]*fIn[7]+3.0*fIn[5]*fIn[5]+3.0*fIn[4]*fIn[4]+3.0*fIn[1]*fIn[1]);
  phixSq[1] = 0.7071067811865475*(6.708203932499369*fIn[10]*fIn[17]+6.708203932499369*fIn[5]*fIn[13]+6.708203932499369*fIn[4]*fIn[11]+6.708203932499369*fIn[1]*fIn[7]);
  phixSq[2] = 0.1414213562373095*(15.0*fIn[15]*fIn[19]+13.41640786499874*fIn[10]*fIn[18]+75.00000000000001*fIn[13]*fIn[17]+13.41640786499874*fIn[4]*fIn[12]+75.00000000000001*fIn[7]*fIn[11]+15.0*fIn[5]*fIn[10]+15.0*fIn[1]*fIn[4]);
  phixSq[3] = 0.1414213562373095*(13.41640786499874*fIn[10]*fIn[19]+15.0*fIn[12]*fIn[18]+75.00000000000001*fIn[11]*fIn[17]+13.41640786499874*fIn[5]*fIn[15]+75.00000000000001*fIn[7]*fIn[13]+15.0*fIn[4]*fIn[10]+15.0*fIn[1]*fIn[5]);
  phixSq[4] = 0.7071067811865475*(6.0*fIn[17]*fIn[18]+6.708203932499369*fIn[5]*fIn[17]+6.708203932499369*fIn[10]*fIn[13]+6.0*fIn[11]*fIn[12]+6.708203932499369*fIn[1]*fIn[11]+6.708203932499369*fIn[4]*fIn[7]);
  phixSq[5] = 0.7071067811865475*(6.0*fIn[17]*fIn[19]+6.708203932499369*fIn[4]*fIn[17]+6.0*fIn[13]*fIn[15]+6.708203932499369*fIn[1]*fIn[13]+6.708203932499369*fIn[10]*fIn[11]+6.708203932499369*fIn[5]*fIn[7]);
  phixSq[6] = 0.1414213562373095*((12.0*fIn[18]+13.41640786499874*fIn[5])*fIn[19]+13.41640786499874*fIn[4]*fIn[18]+75.0*fIn[7]*fIn[17]+13.41640786499874*fIn[10]*fIn[15]+75.0*fIn[11]*fIn[13]+13.41640786499874*fIn[10]*fIn[12]+15.0*fIn[1]*fIn[10]+15.0*fIn[4]*fIn[5]);
  phixSq[7] = 0.7071067811865475*(6.708203932499369*fIn[17]*fIn[17]+6.708203932499369*fIn[13]*fIn[13]+6.708203932499369*fIn[11]*fIn[11]+6.708203932499369*fIn[7]*fIn[7]);
  phixSq[8] = 0.02020305089104421*(46.95742752749558*fIn[19]*fIn[19]+33.54101966249685*fIn[18]*fIn[18]+105.0*fIn[5]*fIn[18]+234.787137637478*fIn[17]*fIn[17]+33.54101966249685*fIn[12]*fIn[12]+105.0*fIn[1]*fIn[12]+234.787137637478*fIn[11]*fIn[11]+46.95742752749558*fIn[10]*fIn[10]+46.95742752749558*fIn[4]*fIn[4]);
  phixSq[9] = 0.02020305089104421*(33.54101966249685*fIn[19]*fIn[19]+105.0*fIn[4]*fIn[19]+46.95742752749558*fIn[18]*fIn[18]+234.787137637478*fIn[17]*fIn[17]+33.54101966249685*fIn[15]*fIn[15]+105.0*fIn[1]*fIn[15]+234.787137637478*fIn[13]*fIn[13]+46.95742752749558*fIn[10]*fIn[10]+46.95742752749558*fIn[5]*fIn[5]);
  phixSq[10] = 0.1414213562373095*(30.0*fIn[13]*fIn[19]+30.0*fIn[11]*fIn[18]+(30.0*fIn[15]+30.0*fIn[12]+33.54101966249685*fIn[1])*fIn[17]+33.54101966249684*fIn[4]*fIn[13]+33.54101966249684*fIn[5]*fIn[11]+33.54101966249685*fIn[7]*fIn[10]);
  phixSq[11] = 0.7071067811865475*(13.41640786499874*fIn[13]*fIn[17]+13.41640786499874*fIn[7]*fIn[11]);
  phixSq[12] = 0.1414213562373095*(33.54101966249685*fIn[13]*fIn[18]+30.0*fIn[10]*fIn[17]+33.54101966249685*fIn[7]*fIn[12]+30.0*fIn[4]*fIn[11]);
  phixSq[13] = 0.7071067811865475*(13.41640786499874*fIn[11]*fIn[17]+13.41640786499874*fIn[7]*fIn[13]);
  phixSq[14] = 0.004040610178208843*(420.0000000000001*fIn[10]*fIn[19]+(469.5742752749559*fIn[15]+335.4101966249685*fIn[12]+525.0000000000001*fIn[1])*fIn[18]+2347.87137637478*fIn[11]*fIn[17]+525.0*fIn[5]*fIn[12]+469.5742752749558*fIn[4]*fIn[10]);
  phixSq[15] = 0.1414213562373095*(33.54101966249685*fIn[11]*fIn[19]+30.0*fIn[10]*fIn[17]+33.54101966249685*fIn[7]*fIn[15]+30.0*fIn[5]*fIn[13]);
  phixSq[16] = 0.004040610178208843*((335.4101966249685*fIn[15]+469.5742752749559*fIn[12]+525.0000000000001*fIn[1])*fIn[19]+420.0000000000001*fIn[10]*fIn[18]+2347.87137637478*fIn[13]*fIn[17]+525.0*fIn[4]*fIn[15]+469.5742752749558*fIn[5]*fIn[10]);
  phixSq[17] = 0.7071067811865475*(13.41640786499874*fIn[7]*fIn[17]+13.41640786499874*fIn[11]*fIn[13]);
  phixSq[18] = 0.1414213562373095*(26.83281572999748*fIn[17]*fIn[19]+33.54101966249685*fIn[7]*fIn[18]+30.0*fIn[4]*fIn[17]+33.54101966249685*fIn[12]*fIn[13]+30.0*fIn[10]*fIn[11]);
  phixSq[19] = 0.1414213562373095*(33.54101966249685*fIn[7]*fIn[19]+26.83281572999748*fIn[17]*fIn[18]+30.0*fIn[5]*fIn[17]+33.54101966249685*fIn[11]*fIn[15]+30.0*fIn[10]*fIn[13]);

  double phiySq[20] = {0.};
  phiySq[0] = 0.3535533905932737*(3.0*fIn[19]*fIn[19]+15.0*fIn[18]*fIn[18]+3.0*fIn[17]*fIn[17]+3.0*fIn[16]*fIn[16]+15.0*fIn[14]*fIn[14]+15.0*fIn[12]*fIn[12]+3.0*fIn[11]*fIn[11]+3.0*fIn[10]*fIn[10]+15.0*fIn[8]*fIn[8]+3.0*fIn[6]*fIn[6]+3.0*fIn[4]*fIn[4]+3.0*fIn[2]*fIn[2]);
  phiySq[1] = 0.1414213562373095*(15.0*fIn[16]*fIn[19]+75.00000000000001*fIn[14]*fIn[18]+13.41640786499874*fIn[10]*fIn[17]+75.00000000000001*fIn[8]*fIn[12]+13.41640786499874*fIn[4]*fIn[11]+15.0*fIn[6]*fIn[10]+15.0*fIn[2]*fIn[4]);
  phiySq[2] = 0.7071067811865475*(6.708203932499369*fIn[10]*fIn[18]+6.708203932499369*fIn[6]*fIn[14]+6.708203932499369*fIn[4]*fIn[12]+6.708203932499369*fIn[2]*fIn[8]);
  phiySq[3] = 0.1414213562373095*(13.41640786499874*fIn[10]*fIn[19]+75.00000000000001*fIn[12]*fIn[18]+15.0*fIn[11]*fIn[17]+13.41640786499874*fIn[6]*fIn[16]+75.00000000000001*fIn[8]*fIn[14]+15.0*fIn[4]*fIn[10]+15.0*fIn[2]*fIn[6]);
  phiySq[4] = 0.7071067811865475*((6.0*fIn[17]+6.708203932499369*fIn[6])*fIn[18]+6.708203932499369*fIn[10]*fIn[14]+(6.0*fIn[11]+6.708203932499369*fIn[2])*fIn[12]+6.708203932499369*fIn[4]*fIn[8]);
  phiySq[5] = 0.1414213562373095*((12.0*fIn[17]+13.41640786499874*fIn[6])*fIn[19]+75.0*fIn[8]*fIn[18]+13.41640786499874*fIn[4]*fIn[17]+13.41640786499874*fIn[10]*fIn[16]+75.0*fIn[12]*fIn[14]+13.41640786499874*fIn[10]*fIn[11]+15.0*fIn[2]*fIn[10]+15.0*fIn[4]*fIn[6]);
  phiySq[6] = 0.7071067811865475*(6.0*fIn[18]*fIn[19]+6.708203932499369*fIn[4]*fIn[18]+6.0*fIn[14]*fIn[16]+6.708203932499369*fIn[2]*fIn[14]+6.708203932499369*fIn[10]*fIn[12]+6.708203932499369*fIn[6]*fIn[8]);
  phiySq[7] = 0.02020305089104421*(46.95742752749558*fIn[19]*fIn[19]+234.787137637478*fIn[18]*fIn[18]+33.54101966249685*fIn[17]*fIn[17]+105.0*fIn[6]*fIn[17]+234.787137637478*fIn[12]*fIn[12]+33.54101966249685*fIn[11]*fIn[11]+105.0*fIn[2]*fIn[11]+46.95742752749558*fIn[10]*fIn[10]+46.95742752749558*fIn[4]*fIn[4]);
  phiySq[8] = 0.7071067811865475*(6.708203932499369*fIn[18]*fIn[18]+6.708203932499369*fIn[14]*fIn[14]+6.708203932499369*fIn[12]*fIn[12]+6.708203932499369*fIn[8]*fIn[8]);
  phiySq[9] = 0.02020305089104421*(33.54101966249685*fIn[19]*fIn[19]+105.0*fIn[4]*fIn[19]+234.787137637478*fIn[18]*fIn[18]+46.95742752749558*fIn[17]*fIn[17]+33.54101966249685*fIn[16]*fIn[16]+105.0*fIn[2]*fIn[16]+234.787137637478*fIn[14]*fIn[14]+46.95742752749558*fIn[10]*fIn[10]+46.95742752749558*fIn[6]*fIn[6]);
  phiySq[10] = 0.1414213562373095*(30.0*fIn[14]*fIn[19]+(30.0*fIn[16]+30.0*fIn[11]+33.54101966249685*fIn[2])*fIn[18]+30.0*fIn[12]*fIn[17]+33.54101966249684*fIn[4]*fIn[14]+33.54101966249684*fIn[6]*fIn[12]+33.54101966249685*fIn[8]*fIn[10]);
  phiySq[11] = 0.1414213562373095*(30.0*fIn[10]*fIn[18]+33.54101966249685*fIn[14]*fIn[17]+30.0*fIn[4]*fIn[12]+33.54101966249685*fIn[8]*fIn[11]);
  phiySq[12] = 0.7071067811865475*(13.41640786499874*fIn[14]*fIn[18]+13.41640786499874*fIn[8]*fIn[12]);
  phiySq[13] = 0.004040610178208843*(420.0000000000001*fIn[10]*fIn[19]+2347.87137637478*fIn[12]*fIn[18]+(469.5742752749559*fIn[16]+335.4101966249685*fIn[11]+525.0000000000001*fIn[2])*fIn[17]+525.0*fIn[6]*fIn[11]+469.5742752749558*fIn[4]*fIn[10]);
  phiySq[14] = 0.7071067811865475*(13.41640786499874*fIn[12]*fIn[18]+13.41640786499874*fIn[8]*fIn[14]);
  phiySq[15] = 0.004040610178208843*((335.4101966249685*fIn[16]+469.5742752749559*fIn[11]+525.0000000000001*fIn[2])*fIn[19]+2347.87137637478*fIn[14]*fIn[18]+420.0000000000001*fIn[10]*fIn[17]+525.0*fIn[4]*fIn[16]+469.5742752749558*fIn[6]*fIn[10]);
  phiySq[16] = 0.1414213562373095*(33.54101966249685*fIn[12]*fIn[19]+30.0*fIn[10]*fIn[18]+33.54101966249685*fIn[8]*fIn[16]+30.0*fIn[6]*fIn[14]);
  phiySq[17] = 0.1414213562373095*(26.83281572999748*fIn[18]*fIn[19]+30.0*fIn[4]*fIn[18]+33.54101966249685*fIn[8]*fIn[17]+33.54101966249685*fIn[11]*fIn[14]+30.0*fIn[10]*fIn[12]);
  phiySq[18] = 0.7071067811865475*(13.41640786499874*fIn[8]*fIn[18]+13.41640786499874*fIn[12]*fIn[14]);
  phiySq[19] = 0.1414213562373095*(33.54101966249685*fIn[8]*fIn[19]+(26.83281572999748*fIn[17]+30.0*fIn[6])*fIn[18]+33.54101966249685*fIn[12]*fIn[16]+30.0*fIn[10]*fIn[14]);

  // assume 1 component
  out[0] += vol*((inw[19]*phiySq[19]+inw[18]*phiySq[18]+inw[17]*phiySq[17]+inw[16]*phiySq[16]+inw[15]*phiySq[15]+inw[14]*phiySq[14]+inw[13]*phiySq[13]+inw[12]*phiySq[12]+inw[11]*phiySq[11]+inw[10]*phiySq[10]+inw[9]*phiySq[9]+inw[8]*phiySq[8]+inw[7]*phiySq[7]+inw[6]*phiySq[6]+inw[5]*phiySq[5]+inw[4]*phiySq[4]+inw[3]*phiySq[3]+inw[2]*phiySq[2]+inw[1]*phiySq[1]+inw[0]*phiySq[0])*dfac2_y+(inw[19]*phixSq[19]+inw[18]*phixSq[18]+inw[17]*phixSq[17]+inw[16]*phixSq[16]+inw[15]*phixSq[15]+inw[14]*phixSq[14]+inw[13]*phixSq[13]+inw[12]*phixSq[12]+inw[11]*phixSq[11]+inw[10]*phixSq[10]+inw[9]*phixSq[9]+inw[8]*phixSq[8]+inw[7]*phixSq[7]+inw[6]*phixSq[6]+inw[5]*phixSq[5]+inw[4]*phixSq[4]+inw[3]*phixSq[3]+inw[2]*phixSq[2]+inw[1]*phixSq[1]+inw[0]*phixSq[0])*dfac2_x);
}
