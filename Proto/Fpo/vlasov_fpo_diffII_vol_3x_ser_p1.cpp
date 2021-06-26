#include <math.h>
#include <vlasov_fpo.h>

double vlasov_fpo_diffII_vol_3x_ser_p1(const double dt, const double* dv,
                                      const double* fIn,
                                      const double* gIn,
                                      double *fOut) {
  const double Jxx = 16/(dv[0]*dv[0]*dv[0]*dv[0]);
  const double Jyy = 16/(dv[1]*dv[1]*dv[1]*dv[1]);
  const double Jzz = 16/(dv[2]*dv[2]*dv[2]*dv[2]);
  const double Jxy = 16/(dv[0]*dv[0]*dv[1]*dv[1]);
  const double Jyx = Jxy;
  const double Jxz = 16/(dv[0]*dv[0]*dv[2]*dv[2]);
  const double Jzx = Jxz;
  const double Jyz = 16/(dv[1]*dv[1]*dv[2]*dv[2]);
  const double Jzy = Jyz;

  fOut[0] += 0.0;
  fOut[1] += ((-1.590990257669731*fIn[6]*gIn[7]*Jzx)-1.590990257669731*fIn[3]*gIn[5]*Jzx-1.590990257669731*fIn[6]*gIn[7]*Jyx-1.590990257669731*fIn[2]*gIn[4]*Jyx)*dt;
  fOut[2] += ((-1.590990257669731*fIn[5]*gIn[7]*Jzy)-1.590990257669731*fIn[3]*gIn[6]*Jzy-1.590990257669731*fIn[5]*gIn[7]*Jxy-1.590990257669731*fIn[1]*gIn[4]*Jxy)*dt;
  fOut[3] += ((-1.590990257669731*fIn[4]*gIn[7]*Jyz)-1.590990257669731*fIn[2]*gIn[6]*Jyz-1.590990257669731*fIn[4]*gIn[7]*Jxz-1.590990257669731*fIn[1]*gIn[5]*Jxz)*dt;
  fOut[4] += ((-1.590990257669731*fIn[3]*gIn[7]*Jzy)-1.590990257669731*fIn[5]*gIn[6]*Jzy-1.590990257669731*fIn[3]*gIn[7]*Jzx-1.590990257669731*gIn[5]*fIn[6]*Jzx)*dt;
  fOut[5] += ((-1.590990257669731*fIn[2]*gIn[7]*Jyz)-1.590990257669731*fIn[4]*gIn[6]*Jyz-1.590990257669731*fIn[2]*gIn[7]*Jyx-1.590990257669731*gIn[4]*fIn[6]*Jyx)*dt;
  fOut[6] += ((-1.590990257669731*fIn[1]*gIn[7]*Jxz)-1.590990257669731*fIn[4]*gIn[5]*Jxz-1.590990257669731*fIn[1]*gIn[7]*Jxy-1.590990257669731*gIn[4]*fIn[5]*Jxy)*dt;
  fOut[7] += 0.0;

  return 0.7954951288348656*gIn[6]*Jyz+0.7954951288348656*gIn[5]*Jxz+0.7954951288348656*gIn[4]*Jxy;
}