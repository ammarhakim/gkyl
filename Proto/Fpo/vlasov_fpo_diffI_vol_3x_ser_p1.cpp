#include <math.h>
#include <vlasov_fpo.h>

double vlasov_fpo_diffI_vol_3x_ser_p1(const double dt, const double* dv,
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
  fOut[1] += 0.0;
  fOut[2] += 0.0;
  fOut[3] += 0.0;
  fOut[4] += (3.181980515339463*fIn[3]*gIn[7]*Jxy+3.181980515339463*fIn[0]*gIn[4]*Jxy)*dt;
  fOut[5] += (3.181980515339463*fIn[2]*gIn[7]*Jxz+3.181980515339463*fIn[0]*gIn[5]*Jxz)*dt;
  fOut[6] += (3.181980515339463*fIn[1]*gIn[7]*Jyz+3.181980515339463*fIn[0]*gIn[6]*Jyz)*dt;
  fOut[7] += (3.181980515339463*fIn[0]*gIn[7]*Jyz+3.181980515339463*fIn[1]*gIn[6]*Jyz+3.181980515339463*fIn[0]*gIn[7]*Jxz+3.181980515339463*fIn[2]*gIn[5]*Jxz+3.181980515339463*fIn[0]*gIn[7]*Jxy+3.181980515339463*fIn[3]*gIn[4]*Jxy)*dt;

  return 0.7954951288348656*gIn[6]*Jyz+0.7954951288348656*gIn[5]*Jxz+0.7954951288348656*gIn[4]*Jxy;
}