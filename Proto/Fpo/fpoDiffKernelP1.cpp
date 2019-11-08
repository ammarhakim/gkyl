#include <math.h>
#include <fpoKernelsDecl.h>

void fpoDiffKernelP1(const double dt, const double *dv, const double *fTL, const double *fT, const double *fTR, const double *fL, const double *fC, const double *fR, const double *fBL, const double *fB, const double *fBR, const double *gTL, const double *gT, const double *gTR, const double *gL, const double *gC, const double *gR, const double *gBL, const double *gB, const double *gBR, const int isTopEdge, const int isBotEdge, const int isLeftEdge, const int isRightEdge, double *fOut) {
  double Jxx = 16/(dv[0]*dv[0]*dv[0]*dv[0]);
  double Jyy = 16/(dv[1]*dv[1]*dv[1]*dv[1]);
  double Jxy = 16/(dv[0]*dv[0]*dv[1]*dv[1]);
  double Jyx = Jxy;

  fOut[0] += ((-0.29296875*fT[3]*gT[3]*Jyy)+0.05859375*fC[3]*gT[3]*Jyy+0.3247595264191644*fT[1]*gT[3]*Jyy+0.0811898816047911*fC[1]*gT[3]*Jyy-0.05859375*fT[3]*gC[3]*Jyy+0.5859375*fC[3]*gC[3]*Jyy-0.05859375*fB[3]*gC[3]*Jyy+0.0811898816047911*fT[1]*gC[3]*Jyy-0.0811898816047911*fB[1]*gC[3]*Jyy+0.05859375*fC[3]*gB[3]*Jyy-0.29296875*fB[3]*gB[3]*Jyy-0.0811898816047911*fC[1]*gB[3]*Jyy-0.3247595264191644*fB[1]*gB[3]*Jyy+0.1014873520059889*gT[1]*fT[3]*Jyy-0.1014873520059889*gC[1]*fT[3]*Jyy-0.1014873520059889*gT[1]*fC[3]*Jyy+0.1014873520059889*gB[1]*fC[3]*Jyy+0.1014873520059889*gC[1]*fB[3]*Jyy-0.1014873520059889*gB[1]*fB[3]*Jyy-0.29296875*fT[2]*gT[2]*Jyy+0.05859375*fC[2]*gT[2]*Jyy+0.3247595264191644*fT[0]*gT[2]*Jyy+0.0811898816047911*fC[0]*gT[2]*Jyy-0.05859375*fT[2]*gC[2]*Jyy+0.5859375*fC[2]*gC[2]*Jyy-0.05859375*fB[2]*gC[2]*Jyy+0.0811898816047911*fT[0]*gC[2]*Jyy-0.0811898816047911*fB[0]*gC[2]*Jyy+0.05859375*fC[2]*gB[2]*Jyy-0.29296875*fB[2]*gB[2]*Jyy-0.0811898816047911*fC[0]*gB[2]*Jyy-0.3247595264191644*fB[0]*gB[2]*Jyy+0.1014873520059889*gT[0]*fT[2]*Jyy-0.1014873520059889*gC[0]*fT[2]*Jyy-0.1014873520059889*gT[0]*fC[2]*Jyy+0.1014873520059889*gB[0]*fC[2]*Jyy+0.1014873520059889*gC[0]*fB[2]*Jyy-0.1014873520059889*gB[0]*fB[2]*Jyy-0.1171875*fT[1]*gT[1]*Jyy-0.1171875*fC[1]*gT[1]*Jyy+0.1171875*fT[1]*gC[1]*Jyy+0.234375*fC[1]*gC[1]*Jyy+0.1171875*fB[1]*gC[1]*Jyy-0.1171875*fC[1]*gB[1]*Jyy-0.1171875*fB[1]*gB[1]*Jyy-0.1171875*fT[0]*gT[0]*Jyy-0.1171875*fC[0]*gT[0]*Jyy+0.1171875*fT[0]*gC[0]*Jyy+0.234375*fC[0]*gC[0]*Jyy+0.1171875*fB[0]*gC[0]*Jyy-0.1171875*fC[0]*gB[0]*Jyy-0.1171875*fB[0]*gB[0]*Jyy+0.01171875*fTR[3]*gTR[3]*Jyx+0.01953125*fR[3]*gTR[3]*Jyx-0.01014873520059889*fTR[2]*gTR[3]*Jyx-0.01691455866766482*fR[2]*gTR[3]*Jyx+0.02480801937924173*fTR[1]*gTR[3]*Jyx+0.04285021529141753*fR[1]*gTR[3]*Jyx-0.021484375*fTR[0]*gTR[3]*Jyx-0.037109375*fR[0]*gTR[3]*Jyx+0.01171875*fTL[3]*gTL[3]*Jyx+0.01953125*fL[3]*gTL[3]*Jyx+0.01014873520059889*fTL[2]*gTL[3]*Jyx+0.01691455866766482*fL[2]*gTL[3]*Jyx+0.02480801937924173*fTL[1]*gTL[3]*Jyx+0.04285021529141753*fL[1]*gTL[3]*Jyx+0.021484375*fTL[0]*gTL[3]*Jyx+0.037109375*fL[0]*gTL[3]*Jyx-0.0234375*fT[3]*gT[3]*Jyx-0.0390625*fC[3]*gT[3]*Jyx-0.04961603875848346*fT[1]*gT[3]*Jyx-0.08570043058283507*fC[1]*gT[3]*Jyx-0.01953125*fTR[3]*gR[3]*Jyx-0.0234375*fR[3]*gR[3]*Jyx-0.01953125*fBR[3]*gR[3]*Jyx+0.01691455866766482*fTR[2]*gR[3]*Jyx+0.02029747040119777*fR[2]*gR[3]*Jyx+0.01691455866766482*fBR[2]*gR[3]*Jyx+0.04285021529141753*fTR[1]*gR[3]*Jyx-0.04285021529141753*fBR[1]*gR[3]*Jyx-0.037109375*fTR[0]*gR[3]*Jyx+0.037109375*fBR[0]*gR[3]*Jyx-0.01953125*fTL[3]*gL[3]*Jyx-0.01171875*fL[3]*gL[3]*Jyx-0.01691455866766482*fTL[2]*gL[3]*Jyx-0.01014873520059889*fL[2]*gL[3]*Jyx+0.04285021529141753*fTL[1]*gL[3]*Jyx+0.02480801937924173*fL[1]*gL[3]*Jyx+0.037109375*fTL[0]*gL[3]*Jyx+0.021484375*fL[0]*gL[3]*Jyx+0.0390625*fT[3]*gC[3]*Jyx+0.046875*fC[3]*gC[3]*Jyx+0.0390625*fB[3]*gC[3]*Jyx-0.08570043058283507*fT[1]*gC[3]*Jyx+0.08570043058283507*fB[1]*gC[3]*Jyx+0.01953125*fR[3]*gBR[3]*Jyx+0.01171875*fBR[3]*gBR[3]*Jyx-0.01691455866766482*fR[2]*gBR[3]*Jyx-0.01014873520059889*fBR[2]*gBR[3]*Jyx-0.04285021529141753*fR[1]*gBR[3]*Jyx-0.02480801937924173*fBR[1]*gBR[3]*Jyx+0.037109375*fR[0]*gBR[3]*Jyx+0.021484375*fBR[0]*gBR[3]*Jyx-0.1353164693413185*fBL[1]*gBL[3]*Jyx-0.1171875*fBL[0]*gBL[3]*Jyx-0.0390625*fC[3]*gB[3]*Jyx-0.0234375*fB[3]*gB[3]*Jyx+0.08570043058283507*fC[1]*gB[3]*Jyx+0.04961603875848346*fB[1]*gB[3]*Jyx+0.02029747040119777*gTR[1]*fTR[3]*Jyx-0.02029747040119777*gR[1]*fTR[3]*Jyx+0.02029747040119777*gTL[1]*fTL[3]*Jyx-0.02029747040119777*gL[1]*fTL[3]*Jyx-0.04059494080239555*gT[1]*fT[3]*Jyx+0.04059494080239555*gC[1]*fT[3]*Jyx-0.02029747040119777*gTR[1]*fR[3]*Jyx+0.02029747040119777*gBR[1]*fR[3]*Jyx-0.02029747040119777*gTL[1]*fL[3]*Jyx+0.02029747040119777*gL[1]*fL[3]*Jyx+0.04059494080239555*gT[1]*fC[3]*Jyx-0.04059494080239555*gB[1]*fC[3]*Jyx+0.02029747040119777*gR[1]*fBR[3]*Jyx-0.02029747040119777*gBR[1]*fBR[3]*Jyx-0.04059494080239555*gC[1]*fB[3]*Jyx+0.04059494080239555*gB[1]*fB[3]*Jyx-0.017578125*gTR[1]*fTR[2]*Jyx+0.017578125*gR[1]*fTR[2]*Jyx+0.017578125*gTL[1]*fTL[2]*Jyx-0.017578125*gL[1]*fTL[2]*Jyx+0.017578125*gTR[1]*fR[2]*Jyx-0.017578125*gBR[1]*fR[2]*Jyx-0.017578125*gTL[1]*fL[2]*Jyx+0.017578125*gL[1]*fL[2]*Jyx-0.017578125*gR[1]*fBR[2]*Jyx+0.017578125*gBR[1]*fBR[2]*Jyx-0.03515625*fTR[1]*gTR[1]*Jyx-0.03515625*fR[1]*gTR[1]*Jyx+0.03044620560179666*fTR[0]*gTR[1]*Jyx+0.03044620560179666*fR[0]*gTR[1]*Jyx-0.03515625*fTL[1]*gTL[1]*Jyx-0.03515625*fL[1]*gTL[1]*Jyx-0.03044620560179666*fTL[0]*gTL[1]*Jyx-0.03044620560179666*fL[0]*gTL[1]*Jyx+0.0703125*fT[1]*gT[1]*Jyx+0.0703125*fC[1]*gT[1]*Jyx+0.03515625*fTR[1]*gR[1]*Jyx+0.0703125*fR[1]*gR[1]*Jyx+0.03515625*fBR[1]*gR[1]*Jyx-0.03044620560179666*fTR[0]*gR[1]*Jyx-0.06089241120359332*fR[0]*gR[1]*Jyx-0.03044620560179666*fBR[0]*gR[1]*Jyx+0.03515625*fTL[1]*gL[1]*Jyx+0.03515625*fL[1]*gL[1]*Jyx+0.03044620560179666*fTL[0]*gL[1]*Jyx+0.03044620560179666*fL[0]*gL[1]*Jyx-0.0703125*fT[1]*gC[1]*Jyx-0.140625*fC[1]*gC[1]*Jyx-0.0703125*fB[1]*gC[1]*Jyx-0.03515625*fR[1]*gBR[1]*Jyx-0.03515625*fBR[1]*gBR[1]*Jyx+0.03044620560179666*fR[0]*gBR[1]*Jyx+0.03044620560179666*fBR[0]*gBR[1]*Jyx+0.0703125*fC[1]*gB[1]*Jyx+0.0703125*fB[1]*gB[1]*Jyx+0.01171875*fTR[3]*gTR[3]*Jxy+0.01953125*fT[3]*gTR[3]*Jxy+0.02480801937924173*fTR[2]*gTR[3]*Jxy+0.04285021529141753*fT[2]*gTR[3]*Jxy-0.01014873520059889*fTR[1]*gTR[3]*Jxy-0.01691455866766482*fT[1]*gTR[3]*Jxy-0.021484375*fTR[0]*gTR[3]*Jxy-0.037109375*fT[0]*gTR[3]*Jxy+0.01171875*fTL[3]*gTL[3]*Jxy+0.01953125*fT[3]*gTL[3]*Jxy-0.02480801937924173*fTL[2]*gTL[3]*Jxy-0.04285021529141753*fT[2]*gTL[3]*Jxy-0.01014873520059889*fTL[1]*gTL[3]*Jxy-0.01691455866766482*fT[1]*gTL[3]*Jxy+0.021484375*fTL[0]*gTL[3]*Jxy+0.037109375*fT[0]*gTL[3]*Jxy-0.01953125*fTR[3]*gT[3]*Jxy-0.01953125*fTL[3]*gT[3]*Jxy-0.0234375*fT[3]*gT[3]*Jxy+0.04285021529141753*fTR[2]*gT[3]*Jxy-0.04285021529141753*fTL[2]*gT[3]*Jxy+0.01691455866766482*fTR[1]*gT[3]*Jxy+0.01691455866766482*fTL[1]*gT[3]*Jxy+0.02029747040119777*fT[1]*gT[3]*Jxy-0.037109375*fTR[0]*gT[3]*Jxy+0.037109375*fTL[0]*gT[3]*Jxy-0.0234375*fR[3]*gR[3]*Jxy-0.0390625*fC[3]*gR[3]*Jxy-0.04961603875848346*fR[2]*gR[3]*Jxy-0.08570043058283507*fC[2]*gR[3]*Jxy-0.0234375*fL[3]*gL[3]*Jxy-0.0390625*fC[3]*gL[3]*Jxy+0.04961603875848346*fL[2]*gL[3]*Jxy+0.08570043058283507*fC[2]*gL[3]*Jxy+0.0390625*fR[3]*gC[3]*Jxy+0.0390625*fL[3]*gC[3]*Jxy+0.046875*fC[3]*gC[3]*Jxy-0.08570043058283507*fR[2]*gC[3]*Jxy+0.08570043058283507*fL[2]*gC[3]*Jxy+0.01171875*fBR[3]*gBR[3]*Jxy+0.01953125*fB[3]*gBR[3]*Jxy+0.02480801937924173*fBR[2]*gBR[3]*Jxy+0.04285021529141753*fB[2]*gBR[3]*Jxy+0.01014873520059889*fBR[1]*gBR[3]*Jxy+0.01691455866766482*fB[1]*gBR[3]*Jxy+0.021484375*fBR[0]*gBR[3]*Jxy+0.037109375*fB[0]*gBR[3]*Jxy+0.01171875*fBL[3]*gBL[3]*Jxy+0.01953125*fB[3]*gBL[3]*Jxy-0.02480801937924173*fBL[2]*gBL[3]*Jxy-0.04285021529141753*fB[2]*gBL[3]*Jxy+0.01014873520059889*fBL[1]*gBL[3]*Jxy+0.01691455866766482*fB[1]*gBL[3]*Jxy-0.021484375*fBL[0]*gBL[3]*Jxy-0.037109375*fB[0]*gBL[3]*Jxy-0.01953125*fBR[3]*gB[3]*Jxy-0.01953125*fBL[3]*gB[3]*Jxy-0.0234375*fB[3]*gB[3]*Jxy+0.04285021529141753*fBR[2]*gB[3]*Jxy-0.04285021529141753*fBL[2]*gB[3]*Jxy-0.01691455866766482*fBR[1]*gB[3]*Jxy-0.01691455866766482*fBL[1]*gB[3]*Jxy-0.02029747040119777*fB[1]*gB[3]*Jxy+0.037109375*fBR[0]*gB[3]*Jxy-0.037109375*fBL[0]*gB[3]*Jxy+0.02029747040119777*gTR[2]*fTR[3]*Jxy-0.02029747040119777*gT[2]*fTR[3]*Jxy-0.02029747040119777*gTL[2]*fTL[3]*Jxy+0.02029747040119777*gT[2]*fTL[3]*Jxy-0.02029747040119777*gTR[2]*fT[3]*Jxy+0.02029747040119777*gTL[2]*fT[3]*Jxy-0.04059494080239555*gR[2]*fR[3]*Jxy+0.04059494080239555*gC[2]*fR[3]*Jxy+0.04059494080239555*gL[2]*fL[3]*Jxy-0.04059494080239555*gC[2]*fL[3]*Jxy+0.04059494080239555*gR[2]*fC[3]*Jxy-0.04059494080239555*gL[2]*fC[3]*Jxy+0.02029747040119777*gBR[2]*fBR[3]*Jxy-0.02029747040119777*gB[2]*fBR[3]*Jxy-0.02029747040119777*gBL[2]*fBL[3]*Jxy+0.02029747040119777*gB[2]*fBL[3]*Jxy-0.02029747040119777*gBR[2]*fB[3]*Jxy+0.02029747040119777*gBL[2]*fB[3]*Jxy-0.03515625*fTR[2]*gTR[2]*Jxy-0.03515625*fT[2]*gTR[2]*Jxy-0.017578125*fTR[1]*gTR[2]*Jxy+0.017578125*fT[1]*gTR[2]*Jxy+0.03044620560179666*fTR[0]*gTR[2]*Jxy+0.03044620560179666*fT[0]*gTR[2]*Jxy-0.03515625*fTL[2]*gTL[2]*Jxy-0.03515625*fT[2]*gTL[2]*Jxy+0.017578125*fTL[1]*gTL[2]*Jxy-0.017578125*fT[1]*gTL[2]*Jxy+0.03044620560179666*fTL[0]*gTL[2]*Jxy+0.03044620560179666*fT[0]*gTL[2]*Jxy+0.03515625*fTR[2]*gT[2]*Jxy+0.03515625*fTL[2]*gT[2]*Jxy+0.0703125*fT[2]*gT[2]*Jxy+0.017578125*fTR[1]*gT[2]*Jxy-0.017578125*fTL[1]*gT[2]*Jxy-0.03044620560179666*fTR[0]*gT[2]*Jxy-0.03044620560179666*fTL[0]*gT[2]*Jxy-0.06089241120359332*fT[0]*gT[2]*Jxy+0.0703125*fR[2]*gR[2]*Jxy+0.0703125*fC[2]*gR[2]*Jxy+0.0703125*fL[2]*gL[2]*Jxy+0.0703125*fC[2]*gL[2]*Jxy-0.0703125*fR[2]*gC[2]*Jxy-0.0703125*fL[2]*gC[2]*Jxy-0.140625*fC[2]*gC[2]*Jxy-0.03515625*fBR[2]*gBR[2]*Jxy-0.03515625*fB[2]*gBR[2]*Jxy+0.017578125*fBR[1]*gBR[2]*Jxy-0.017578125*fB[1]*gBR[2]*Jxy-0.03044620560179666*fBR[0]*gBR[2]*Jxy-0.03044620560179666*fB[0]*gBR[2]*Jxy-0.03515625*fBL[2]*gBL[2]*Jxy-0.03515625*fB[2]*gBL[2]*Jxy-0.017578125*fBL[1]*gBL[2]*Jxy+0.017578125*fB[1]*gBL[2]*Jxy-0.03044620560179666*fBL[0]*gBL[2]*Jxy-0.03044620560179666*fB[0]*gBL[2]*Jxy+0.03515625*fBR[2]*gB[2]*Jxy+0.03515625*fBL[2]*gB[2]*Jxy+0.0703125*fB[2]*gB[2]*Jxy-0.017578125*fBR[1]*gB[2]*Jxy+0.017578125*fBL[1]*gB[2]*Jxy+0.03044620560179666*fBR[0]*gB[2]*Jxy+0.03044620560179666*fBL[0]*gB[2]*Jxy+0.06089241120359332*fB[0]*gB[2]*Jxy-0.29296875*fR[3]*gR[3]*Jxx+0.05859375*fC[3]*gR[3]*Jxx+0.3247595264191644*fR[2]*gR[3]*Jxx+0.0811898816047911*fC[2]*gR[3]*Jxx-0.29296875*fL[3]*gL[3]*Jxx+0.05859375*fC[3]*gL[3]*Jxx-0.3247595264191644*fL[2]*gL[3]*Jxx-0.0811898816047911*fC[2]*gL[3]*Jxx-0.05859375*fR[3]*gC[3]*Jxx-0.05859375*fL[3]*gC[3]*Jxx+0.5859375*fC[3]*gC[3]*Jxx+0.0811898816047911*fR[2]*gC[3]*Jxx-0.0811898816047911*fL[2]*gC[3]*Jxx+0.1014873520059889*gR[2]*fR[3]*Jxx-0.1014873520059889*gC[2]*fR[3]*Jxx-0.1014873520059889*gL[2]*fL[3]*Jxx+0.1014873520059889*gC[2]*fL[3]*Jxx-0.1014873520059889*gR[2]*fC[3]*Jxx+0.1014873520059889*gL[2]*fC[3]*Jxx-0.1171875*fR[2]*gR[2]*Jxx-0.1171875*fC[2]*gR[2]*Jxx-0.1171875*fL[2]*gL[2]*Jxx-0.1171875*fC[2]*gL[2]*Jxx+0.1171875*fR[2]*gC[2]*Jxx+0.1171875*fL[2]*gC[2]*Jxx+0.234375*fC[2]*gC[2]*Jxx-0.29296875*fR[1]*gR[1]*Jxx+0.05859375*fC[1]*gR[1]*Jxx+0.3247595264191644*fR[0]*gR[1]*Jxx+0.0811898816047911*fC[0]*gR[1]*Jxx-0.29296875*fL[1]*gL[1]*Jxx+0.05859375*fC[1]*gL[1]*Jxx-0.3247595264191644*fL[0]*gL[1]*Jxx-0.0811898816047911*fC[0]*gL[1]*Jxx-0.05859375*fR[1]*gC[1]*Jxx-0.05859375*fL[1]*gC[1]*Jxx+0.5859375*fC[1]*gC[1]*Jxx+0.0811898816047911*fR[0]*gC[1]*Jxx-0.0811898816047911*fL[0]*gC[1]*Jxx+0.1014873520059889*gR[0]*fR[1]*Jxx-0.1014873520059889*gC[0]*fR[1]*Jxx-0.1014873520059889*gL[0]*fL[1]*Jxx+0.1014873520059889*gC[0]*fL[1]*Jxx-0.1014873520059889*gR[0]*fC[1]*Jxx+0.1014873520059889*gL[0]*fC[1]*Jxx-0.1171875*fR[0]*gR[0]*Jxx-0.1171875*fC[0]*gR[0]*Jxx-0.1171875*fL[0]*gL[0]*Jxx-0.1171875*fC[0]*gL[0]*Jxx+0.1171875*fR[0]*gC[0]*Jxx+0.1171875*fL[0]*gC[0]*Jxx+0.234375*fC[0]*gC[0]*Jxx)*dt;
  fOut[1] += ((-0.29296875*fT[2]*gT[3]*Jyy)+0.05859375*fC[2]*gT[3]*Jyy+0.3247595264191644*fT[0]*gT[3]*Jyy+0.08118988160479111*fC[0]*gT[3]*Jyy-0.05859375*fT[2]*gC[3]*Jyy+0.5859375*fC[2]*gC[3]*Jyy-0.05859375*fB[2]*gC[3]*Jyy+0.08118988160479111*fT[0]*gC[3]*Jyy-0.08118988160479111*fB[0]*gC[3]*Jyy+0.05859375*fC[2]*gB[3]*Jyy-0.29296875*fB[2]*gB[3]*Jyy-0.08118988160479111*fC[0]*gB[3]*Jyy-0.3247595264191644*fB[0]*gB[3]*Jyy-0.29296875*gT[2]*fT[3]*Jyy-0.05859375*gC[2]*fT[3]*Jyy+0.1014873520059889*gT[0]*fT[3]*Jyy-0.1014873520059889*gC[0]*fT[3]*Jyy+0.05859375*gT[2]*fC[3]*Jyy+0.5859375*gC[2]*fC[3]*Jyy+0.05859375*gB[2]*fC[3]*Jyy-0.1014873520059889*gT[0]*fC[3]*Jyy+0.1014873520059889*gB[0]*fC[3]*Jyy-0.05859375*gC[2]*fB[3]*Jyy-0.29296875*gB[2]*fB[3]*Jyy+0.1014873520059889*gC[0]*fB[3]*Jyy-0.1014873520059889*gB[0]*fB[3]*Jyy+0.3247595264191644*fT[1]*gT[2]*Jyy+0.08118988160479111*fC[1]*gT[2]*Jyy+0.08118988160479111*fT[1]*gC[2]*Jyy-0.08118988160479111*fB[1]*gC[2]*Jyy-0.08118988160479111*fC[1]*gB[2]*Jyy-0.3247595264191644*fB[1]*gB[2]*Jyy+0.1014873520059889*gT[1]*fT[2]*Jyy-0.1014873520059889*gC[1]*fT[2]*Jyy-0.1014873520059889*gT[1]*fC[2]*Jyy+0.1014873520059889*gB[1]*fC[2]*Jyy+0.1014873520059889*gC[1]*fB[2]*Jyy-0.1014873520059889*gB[1]*fB[2]*Jyy-0.1171875*fT[0]*gT[1]*Jyy-0.1171875*fC[0]*gT[1]*Jyy+0.1171875*fT[0]*gC[1]*Jyy+0.234375*fC[0]*gC[1]*Jyy+0.1171875*fB[0]*gC[1]*Jyy-0.1171875*fC[0]*gB[1]*Jyy-0.1171875*fB[0]*gB[1]*Jyy-0.1171875*gT[0]*fT[1]*Jyy+0.1171875*gC[0]*fT[1]*Jyy-0.1171875*gT[0]*fC[1]*Jyy+0.234375*gC[0]*fC[1]*Jyy-0.1171875*gB[0]*fC[1]*Jyy+0.1171875*gC[0]*fB[1]*Jyy-0.1171875*gB[0]*fB[1]*Jyy+0.02029747040119778*fTR[3]*gTR[3]*Jyx+0.03382911733532963*fR[3]*gTR[3]*Jyx-0.017578125*fTR[2]*gTR[3]*Jyx-0.029296875*fR[2]*gTR[3]*Jyx+0.04296875*fTR[1]*gTR[3]*Jyx+0.07421875*fR[1]*gTR[3]*Jyx-0.0372120290688626*fTR[0]*gTR[3]*Jyx-0.0642753229371263*fR[0]*gTR[3]*Jyx-0.02029747040119778*fTL[3]*gTL[3]*Jyx-0.03382911733532963*fL[3]*gTL[3]*Jyx-0.017578125*fTL[2]*gTL[3]*Jyx-0.029296875*fL[2]*gTL[3]*Jyx-0.04296875*fTL[1]*gTL[3]*Jyx-0.07421875*fL[1]*gTL[3]*Jyx-0.0372120290688626*fTL[0]*gTL[3]*Jyx-0.0642753229371263*fL[0]*gTL[3]*Jyx+0.10546875*fT[2]*gT[3]*Jyx+0.17578125*fC[2]*gT[3]*Jyx+0.2232721744131756*fT[0]*gT[3]*Jyx+0.3856519376227578*fC[0]*gT[3]*Jyx-0.03382911733532963*fTR[3]*gR[3]*Jyx-0.04059494080239556*fR[3]*gR[3]*Jyx-0.03382911733532963*fBR[3]*gR[3]*Jyx+0.029296875*fTR[2]*gR[3]*Jyx+0.03515625*fR[2]*gR[3]*Jyx+0.029296875*fBR[2]*gR[3]*Jyx+0.07421875*fTR[1]*gR[3]*Jyx-0.07421875*fBR[1]*gR[3]*Jyx-0.0642753229371263*fTR[0]*gR[3]*Jyx+0.0642753229371263*fBR[0]*gR[3]*Jyx+0.03382911733532963*fTL[3]*gL[3]*Jyx+0.02029747040119778*fL[3]*gL[3]*Jyx+0.029296875*fTL[2]*gL[3]*Jyx+0.017578125*fL[2]*gL[3]*Jyx-0.07421875*fTL[1]*gL[3]*Jyx-0.04296875*fL[1]*gL[3]*Jyx-0.0642753229371263*fTL[0]*gL[3]*Jyx-0.0372120290688626*fL[0]*gL[3]*Jyx-0.17578125*fT[2]*gC[3]*Jyx-0.2109375*fC[2]*gC[3]*Jyx-0.17578125*fB[2]*gC[3]*Jyx+0.3856519376227578*fT[0]*gC[3]*Jyx-0.3856519376227578*fB[0]*gC[3]*Jyx+0.03382911733532963*fR[3]*gBR[3]*Jyx+0.02029747040119778*fBR[3]*gBR[3]*Jyx-0.029296875*fR[2]*gBR[3]*Jyx-0.017578125*fBR[2]*gBR[3]*Jyx-0.07421875*fR[1]*gBR[3]*Jyx-0.04296875*fBR[1]*gBR[3]*Jyx+0.0642753229371263*fR[0]*gBR[3]*Jyx+0.0372120290688626*fBR[0]*gBR[3]*Jyx+0.234375*fBL[1]*gBL[3]*Jyx+0.2029747040119778*fBL[0]*gBL[3]*Jyx+0.17578125*fC[2]*gB[3]*Jyx+0.10546875*fB[2]*gB[3]*Jyx-0.3856519376227578*fC[0]*gB[3]*Jyx-0.2232721744131756*fB[0]*gB[3]*Jyx+0.03515625*gTR[1]*fTR[3]*Jyx-0.03515625*gR[1]*fTR[3]*Jyx-0.03515625*gTL[1]*fTL[3]*Jyx+0.03515625*gL[1]*fTL[3]*Jyx-0.03515625*gTR[1]*fR[3]*Jyx+0.03515625*gBR[1]*fR[3]*Jyx+0.03515625*gTL[1]*fL[3]*Jyx-0.03515625*gL[1]*fL[3]*Jyx+0.03515625*gR[1]*fBR[3]*Jyx-0.03515625*gBR[1]*fBR[3]*Jyx-0.03044620560179666*gTR[1]*fTR[2]*Jyx+0.03044620560179666*gR[1]*fTR[2]*Jyx-0.03044620560179666*gTL[1]*fTL[2]*Jyx+0.03044620560179666*gL[1]*fTL[2]*Jyx+0.1826772336107799*gT[1]*fT[2]*Jyx-0.1826772336107799*gC[1]*fT[2]*Jyx+0.03044620560179666*gTR[1]*fR[2]*Jyx-0.03044620560179666*gBR[1]*fR[2]*Jyx+0.03044620560179666*gTL[1]*fL[2]*Jyx-0.03044620560179666*gL[1]*fL[2]*Jyx-0.1826772336107799*gT[1]*fC[2]*Jyx+0.1826772336107799*gB[1]*fC[2]*Jyx-0.03044620560179666*gR[1]*fBR[2]*Jyx+0.03044620560179666*gBR[1]*fBR[2]*Jyx+0.1826772336107799*gC[1]*fB[2]*Jyx-0.1826772336107799*gB[1]*fB[2]*Jyx-0.06089241120359332*fTR[1]*gTR[1]*Jyx-0.06089241120359332*fR[1]*gTR[1]*Jyx+0.052734375*fTR[0]*gTR[1]*Jyx+0.052734375*fR[0]*gTR[1]*Jyx+0.06089241120359332*fTL[1]*gTL[1]*Jyx+0.06089241120359332*fL[1]*gTL[1]*Jyx+0.052734375*fTL[0]*gTL[1]*Jyx+0.052734375*fL[0]*gTL[1]*Jyx-0.31640625*fT[0]*gT[1]*Jyx-0.31640625*fC[0]*gT[1]*Jyx+0.06089241120359332*fTR[1]*gR[1]*Jyx+0.1217848224071866*fR[1]*gR[1]*Jyx+0.06089241120359332*fBR[1]*gR[1]*Jyx-0.052734375*fTR[0]*gR[1]*Jyx-0.10546875*fR[0]*gR[1]*Jyx-0.052734375*fBR[0]*gR[1]*Jyx-0.06089241120359332*fTL[1]*gL[1]*Jyx-0.06089241120359332*fL[1]*gL[1]*Jyx-0.052734375*fTL[0]*gL[1]*Jyx-0.052734375*fL[0]*gL[1]*Jyx+0.31640625*fT[0]*gC[1]*Jyx+0.6328125*fC[0]*gC[1]*Jyx+0.31640625*fB[0]*gC[1]*Jyx-0.06089241120359332*fR[1]*gBR[1]*Jyx-0.06089241120359332*fBR[1]*gBR[1]*Jyx+0.052734375*fR[0]*gBR[1]*Jyx+0.052734375*fBR[0]*gBR[1]*Jyx-0.31640625*fC[0]*gB[1]*Jyx-0.31640625*fB[0]*gB[1]*Jyx+0.02029747040119778*fTR[3]*gTR[3]*Jxy+0.03382911733532963*fT[3]*gTR[3]*Jxy+0.04296875*fTR[2]*gTR[3]*Jxy+0.07421875*fT[2]*gTR[3]*Jxy-0.017578125*fTR[1]*gTR[3]*Jxy-0.029296875*fT[1]*gTR[3]*Jxy-0.0372120290688626*fTR[0]*gTR[3]*Jxy-0.0642753229371263*fT[0]*gTR[3]*Jxy-0.02029747040119778*fTL[3]*gTL[3]*Jxy-0.03382911733532963*fT[3]*gTL[3]*Jxy+0.04296875*fTL[2]*gTL[3]*Jxy+0.07421875*fT[2]*gTL[3]*Jxy+0.017578125*fTL[1]*gTL[3]*Jxy+0.029296875*fT[1]*gTL[3]*Jxy-0.0372120290688626*fTL[0]*gTL[3]*Jxy-0.0642753229371263*fT[0]*gTL[3]*Jxy-0.03382911733532963*fTR[3]*gT[3]*Jxy+0.03382911733532963*fTL[3]*gT[3]*Jxy+0.07421875*fTR[2]*gT[3]*Jxy+0.07421875*fTL[2]*gT[3]*Jxy+0.0859375*fT[2]*gT[3]*Jxy+0.029296875*fTR[1]*gT[3]*Jxy-0.029296875*fTL[1]*gT[3]*Jxy-0.0642753229371263*fTR[0]*gT[3]*Jxy-0.0642753229371263*fTL[0]*gT[3]*Jxy-0.0744240581377252*fT[0]*gT[3]*Jxy-0.04059494080239556*fR[3]*gR[3]*Jxy-0.06765823467065926*fC[3]*gR[3]*Jxy-0.0859375*fR[2]*gR[3]*Jxy-0.1484375*fC[2]*gR[3]*Jxy+0.04059494080239556*fL[3]*gL[3]*Jxy+0.06765823467065926*fC[3]*gL[3]*Jxy-0.0859375*fL[2]*gL[3]*Jxy-0.1484375*fC[2]*gL[3]*Jxy+0.06765823467065926*fR[3]*gC[3]*Jxy-0.06765823467065926*fL[3]*gC[3]*Jxy-0.1484375*fR[2]*gC[3]*Jxy-0.1484375*fL[2]*gC[3]*Jxy-0.171875*fC[2]*gC[3]*Jxy+0.02029747040119778*fBR[3]*gBR[3]*Jxy+0.03382911733532963*fB[3]*gBR[3]*Jxy+0.04296875*fBR[2]*gBR[3]*Jxy+0.07421875*fB[2]*gBR[3]*Jxy+0.017578125*fBR[1]*gBR[3]*Jxy+0.029296875*fB[1]*gBR[3]*Jxy+0.0372120290688626*fBR[0]*gBR[3]*Jxy+0.0642753229371263*fB[0]*gBR[3]*Jxy-0.02029747040119778*fBL[3]*gBL[3]*Jxy-0.03382911733532963*fB[3]*gBL[3]*Jxy+0.04296875*fBL[2]*gBL[3]*Jxy+0.07421875*fB[2]*gBL[3]*Jxy-0.017578125*fBL[1]*gBL[3]*Jxy-0.029296875*fB[1]*gBL[3]*Jxy+0.0372120290688626*fBL[0]*gBL[3]*Jxy+0.0642753229371263*fB[0]*gBL[3]*Jxy-0.03382911733532963*fBR[3]*gB[3]*Jxy+0.03382911733532963*fBL[3]*gB[3]*Jxy+0.07421875*fBR[2]*gB[3]*Jxy+0.07421875*fBL[2]*gB[3]*Jxy+0.0859375*fB[2]*gB[3]*Jxy-0.029296875*fBR[1]*gB[3]*Jxy+0.029296875*fBL[1]*gB[3]*Jxy+0.0642753229371263*fBR[0]*gB[3]*Jxy+0.0642753229371263*fBL[0]*gB[3]*Jxy+0.0744240581377252*fB[0]*gB[3]*Jxy+0.03515625*gTR[2]*fTR[3]*Jxy-0.03515625*gT[2]*fTR[3]*Jxy+0.03515625*gTL[2]*fTL[3]*Jxy-0.03515625*gT[2]*fTL[3]*Jxy-0.03515625*gTR[2]*fT[3]*Jxy-0.03515625*gTL[2]*fT[3]*Jxy+0.0703125*gT[2]*fT[3]*Jxy-0.0703125*gR[2]*fR[3]*Jxy+0.0703125*gC[2]*fR[3]*Jxy-0.0703125*gL[2]*fL[3]*Jxy+0.0703125*gC[2]*fL[3]*Jxy+0.0703125*gR[2]*fC[3]*Jxy+0.0703125*gL[2]*fC[3]*Jxy-0.140625*gC[2]*fC[3]*Jxy+0.03515625*gBR[2]*fBR[3]*Jxy-0.03515625*gB[2]*fBR[3]*Jxy+0.03515625*gBL[2]*fBL[3]*Jxy-0.03515625*gB[2]*fBL[3]*Jxy-0.03515625*gBR[2]*fB[3]*Jxy-0.03515625*gBL[2]*fB[3]*Jxy+0.0703125*gB[2]*fB[3]*Jxy-0.06089241120359332*fTR[2]*gTR[2]*Jxy-0.06089241120359332*fT[2]*gTR[2]*Jxy-0.03044620560179666*fTR[1]*gTR[2]*Jxy+0.03044620560179666*fT[1]*gTR[2]*Jxy+0.052734375*fTR[0]*gTR[2]*Jxy+0.052734375*fT[0]*gTR[2]*Jxy+0.06089241120359332*fTL[2]*gTL[2]*Jxy+0.06089241120359332*fT[2]*gTL[2]*Jxy-0.03044620560179666*fTL[1]*gTL[2]*Jxy+0.03044620560179666*fT[1]*gTL[2]*Jxy-0.052734375*fTL[0]*gTL[2]*Jxy-0.052734375*fT[0]*gTL[2]*Jxy+0.06089241120359332*fTR[2]*gT[2]*Jxy-0.06089241120359332*fTL[2]*gT[2]*Jxy+0.03044620560179666*fTR[1]*gT[2]*Jxy+0.03044620560179666*fTL[1]*gT[2]*Jxy-0.06089241120359332*fT[1]*gT[2]*Jxy-0.052734375*fTR[0]*gT[2]*Jxy+0.052734375*fTL[0]*gT[2]*Jxy+0.1217848224071866*fR[2]*gR[2]*Jxy+0.1217848224071866*fC[2]*gR[2]*Jxy-0.1217848224071866*fL[2]*gL[2]*Jxy-0.1217848224071866*fC[2]*gL[2]*Jxy-0.1217848224071866*fR[2]*gC[2]*Jxy+0.1217848224071866*fL[2]*gC[2]*Jxy-0.06089241120359332*fBR[2]*gBR[2]*Jxy-0.06089241120359332*fB[2]*gBR[2]*Jxy+0.03044620560179666*fBR[1]*gBR[2]*Jxy-0.03044620560179666*fB[1]*gBR[2]*Jxy-0.052734375*fBR[0]*gBR[2]*Jxy-0.052734375*fB[0]*gBR[2]*Jxy+0.06089241120359332*fBL[2]*gBL[2]*Jxy+0.06089241120359332*fB[2]*gBL[2]*Jxy+0.03044620560179666*fBL[1]*gBL[2]*Jxy-0.03044620560179666*fB[1]*gBL[2]*Jxy+0.052734375*fBL[0]*gBL[2]*Jxy+0.052734375*fB[0]*gBL[2]*Jxy+0.06089241120359332*fBR[2]*gB[2]*Jxy-0.06089241120359332*fBL[2]*gB[2]*Jxy-0.03044620560179666*fBR[1]*gB[2]*Jxy-0.03044620560179666*fBL[1]*gB[2]*Jxy+0.06089241120359332*fB[1]*gB[2]*Jxy+0.052734375*fBR[0]*gB[2]*Jxy-0.052734375*fBL[0]*gB[2]*Jxy-0.1962088805449119*fR[3]*gR[3]*Jxx+0.1962088805449119*fC[3]*gR[3]*Jxx+0.3515625*fR[2]*gR[3]*Jxx+0.1640625*fC[2]*gR[3]*Jxx+0.1962088805449119*fL[3]*gL[3]*Jxx-0.1962088805449119*fC[3]*gL[3]*Jxx+0.3515625*fL[2]*gL[3]*Jxx+0.1640625*fC[2]*gL[3]*Jxx-0.006765823467065927*fR[3]*gC[3]*Jxx+0.006765823467065927*fL[3]*gC[3]*Jxx+0.1171875*fR[2]*gC[3]*Jxx+0.1171875*fL[2]*gC[3]*Jxx+1.546875*fC[2]*gC[3]*Jxx+0.05859375*gR[2]*fR[3]*Jxx-0.05859375*gC[2]*fR[3]*Jxx+0.05859375*gL[2]*fL[3]*Jxx-0.05859375*gC[2]*fL[3]*Jxx-0.29296875*gR[2]*fC[3]*Jxx-0.29296875*gL[2]*fC[3]*Jxx+0.5859375*gC[2]*fC[3]*Jxx-0.1353164693413185*fR[2]*gR[2]*Jxx-0.270632938682637*fC[2]*gR[2]*Jxx+0.1353164693413185*fL[2]*gL[2]*Jxx+0.270632938682637*fC[2]*gL[2]*Jxx+0.1353164693413185*fR[2]*gC[2]*Jxx-0.1353164693413185*fL[2]*gC[2]*Jxx-0.1962088805449119*fR[1]*gR[1]*Jxx+0.1962088805449119*fC[1]*gR[1]*Jxx+0.3515625*fR[0]*gR[1]*Jxx+0.1640625*fC[0]*gR[1]*Jxx+0.1962088805449119*fL[1]*gL[1]*Jxx-0.1962088805449119*fC[1]*gL[1]*Jxx+0.3515625*fL[0]*gL[1]*Jxx+0.1640625*fC[0]*gL[1]*Jxx-0.006765823467065927*fR[1]*gC[1]*Jxx+0.006765823467065927*fL[1]*gC[1]*Jxx+0.1171875*fR[0]*gC[1]*Jxx+0.1171875*fL[0]*gC[1]*Jxx+1.546875*fC[0]*gC[1]*Jxx+0.05859375*gR[0]*fR[1]*Jxx-0.05859375*gC[0]*fR[1]*Jxx+0.05859375*gL[0]*fL[1]*Jxx-0.05859375*gC[0]*fL[1]*Jxx-0.29296875*gR[0]*fC[1]*Jxx-0.29296875*gL[0]*fC[1]*Jxx+0.5859375*gC[0]*fC[1]*Jxx-0.1353164693413185*fR[0]*gR[0]*Jxx-0.270632938682637*fC[0]*gR[0]*Jxx+0.1353164693413185*fL[0]*gL[0]*Jxx+0.270632938682637*fC[0]*gL[0]*Jxx+0.1353164693413185*fR[0]*gC[0]*Jxx-0.1353164693413185*fL[0]*gC[0]*Jxx)*dt;
  fOut[0] += 0.0;
  fOut[1] += 0.0;
}