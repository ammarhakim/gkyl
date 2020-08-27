#include <math.h>
#include <fpoKernelsDecl.h>

double fpoDiffVolSer3xP1(const double dt, const double* dv,
                         const double* fCCC,
                         const double* gCCC,
                         double *fOut) {
  double Jxx = 16/(dv[0]*dv[0]*dv[0]*dv[0]);
  double Jyy = 16/(dv[1]*dv[1]*dv[1]*dv[1]);
  double Jzz = 16/(dv[2]*dv[2]*dv[2]*dv[2]);
  double Jxy = 16/(dv[0]*dv[0]*dv[1]*dv[1]);
  double Jyx = Jxy;
  double Jxz = 16/(dv[0]*dv[0]*dv[2]*dv[2]);
  double Jzx = Jxz;
  double Jyz = 16/(dv[1]*dv[1]*dv[2]*dv[2]);
  double Jzy = Jyz;

  fOut[0] += 0.0;
  fOut[1] += 0.0;
  fOut[2] += 0.0;
  fOut[3] += 0.0;
  fOut[4] += (3.181980515339463*fCCC[3]*gCCC[7]*Jxy+3.181980515339463*fCCC[0]*gCCC[4]*Jxy)*dt;
  fOut[5] += (3.181980515339463*fCCC[2]*gCCC[7]*Jxz+3.181980515339463*fCCC[0]*gCCC[5]*Jxz)*dt;
  fOut[6] += (3.181980515339463*fCCC[1]*gCCC[7]*Jyz+3.181980515339463*fCCC[0]*gCCC[6]*Jyz)*dt;
  fOut[7] += (3.181980515339463*fCCC[0]*gCCC[7]*Jyz+3.181980515339463*fCCC[1]*gCCC[6]*Jyz+3.181980515339463*fCCC[0]*gCCC[7]*Jxz+3.181980515339463*fCCC[2]*gCCC[5]*Jxz+3.181980515339463*fCCC[0]*gCCC[7]*Jxy+3.181980515339463*fCCC[3]*gCCC[4]*Jxy)*dt;

  return 0.7954951288348656*gCCC[6]*Jyz+0.7954951288348656*gCCC[5]*Jxz+0.7954951288348656*gCCC[4]*Jxy;
}