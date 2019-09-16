#include <math.h>
#include <fpoKernelsDecl.h>

void fpoDiffKernelP1(const double dt, const double *dv, const double *f, const double *fL, const double *fR, const double *fT, const double *fB, const double *g, const double *gL, const double *gR, const double *gT, const double *gB, const int isTopEdge, const int isBotEdge, const int isLeftEdge, const int isRightEdge, double *fOut) {
  if (isLeftEdge) {
    fOut[0] += -0.125*((45.0*fR[3]-15.0*f[3]-41.56921938165305*fR[2]-10.39230484541326*f[2])*gR[3]+(15.0*fR[3]-45.0*f[3]-10.39230484541326*fR[2]-41.56921938165305*f[2])*g[3]+(17.32050807568877*g[2]-17.32050807568877*gR[2])*fR[3]+(17.32050807568877*gR[2]-17.32050807568877*g[2])*f[3]+(15.0*fR[2]+15.0*f[2])*gR[2]+((-15.0*fR[2])-15.0*f[2])*g[2]+(45.0*fR[1]-15.0*f[1]-41.56921938165305*fR[0]-10.39230484541326*f[0])*gR[1]+(15.0*fR[1]-45.0*f[1]-10.39230484541326*fR[0]-41.56921938165305*f[0])*g[1]+(17.32050807568877*g[0]-17.32050807568877*gR[0])*fR[1]+(17.32050807568877*gR[0]-17.32050807568877*g[0])*f[1]+(15.0*fR[0]+15.0*f[0])*gR[0]+((-15.0*fR[0])-15.0*f[0])*g[0])*gkyl_ipow(dv[0],-4)*dt;
    fOut[1] += -0.07216878364870323*((111.0*fR[3]-21.0*f[3]-103.9230484541326*fR[2]-10.39230484541326*f[2])*gR[3]+(69.0*fR[3]-159.0*f[3]-51.96152422706631*fR[2]-145.4922678357857*f[2])*g[3]+(51.96152422706631*g[2]-51.96152422706631*gR[2])*fR[3]+(51.96152422706631*gR[2]-51.96152422706631*g[2])*f[3]+(45.0*fR[2]+45.0*f[2])*gR[2]+((-45.0*fR[2])-45.0*f[2])*g[2]+(111.0*fR[1]-21.0*f[1]-103.9230484541326*fR[0]-10.39230484541326*f[0])*gR[1]+(69.0*fR[1]-159.0*f[1]-51.96152422706631*fR[0]-145.4922678357857*f[0])*g[1]+(51.96152422706631*g[0]-51.96152422706631*gR[0])*fR[1]+(51.96152422706631*gR[0]-51.96152422706631*g[0])*f[1]+(45.0*fR[0]+45.0*f[0])*gR[0]+((-45.0*fR[0])-45.0*f[0])*g[0])*gkyl_ipow(dv[0],-4)*dt;
    fOut[2] += -0.125*((45.0*fR[1]-15.0*f[1]-41.56921938165305*fR[0]-10.39230484541326*f[0])*gR[3]+(15.0*fR[1]-45.0*f[1]-10.39230484541326*fR[0]-41.56921938165305*f[0])*g[3]+(45.0*gR[1]+15.0*g[1]-17.32050807568877*gR[0]+17.32050807568877*g[0])*fR[3]+((-15.0*gR[1])-45.0*g[1]+17.32050807568877*gR[0]-17.32050807568877*g[0])*f[3]+((-17.32050807568877*fR[1])+17.32050807568877*f[1]+15.0*fR[0]+15.0*f[0])*gR[2]+(17.32050807568877*fR[1]-17.32050807568877*f[1]-15.0*fR[0]-15.0*f[0])*g[2]+((-41.56921938165305*gR[1])-10.39230484541326*g[1]+15.0*gR[0]-15.0*g[0])*fR[2]+((-10.39230484541326*gR[1])-41.56921938165305*g[1]+15.0*gR[0]-15.0*g[0])*f[2])*gkyl_ipow(dv[0],-4)*dt;
    fOut[3] += -0.125*((64.08587988004845*fR[1]-12.12435565298214*f[1]-60.0*fR[0]-6.0*f[0])*gR[3]+(39.83716857408417*fR[1]-91.7986928011505*f[1]-30.0*fR[0]-84.0*f[0])*g[3]+(64.08587988004845*gR[1]+39.83716857408417*g[1]-30.0*gR[0]+30.0*g[0])*fR[3]+((-12.12435565298214*gR[1])-91.7986928011505*g[1]+30.0*gR[0]-30.0*g[0])*f[3]+((-30.0*fR[1])+30.0*f[1]+25.98076211353316*fR[0]+25.98076211353316*f[0])*gR[2]+(30.0*fR[1]-30.0*f[1]-25.98076211353316*fR[0]-25.98076211353316*f[0])*g[2]+((-60.0*gR[1])-30.0*g[1]+25.98076211353316*gR[0]-25.98076211353316*g[0])*fR[2]+((-6.0*gR[1])-84.0*g[1]+25.98076211353316*gR[0]-25.98076211353316*g[0])*f[2])*gkyl_ipow(dv[0],-4)*dt;
  } else if (isRightEdge) {
    fOut[0] += -0.125*((45.0*fL[3]-15.0*f[3]+41.56921938165305*fL[2]+10.39230484541326*f[2])*gL[3]+(15.0*fL[3]-45.0*f[3]+10.39230484541326*fL[2]+41.56921938165305*f[2])*g[3]+(17.32050807568877*gL[2]-17.32050807568877*g[2])*fL[3]+(17.32050807568877*g[2]-17.32050807568877*gL[2])*f[3]+(15.0*fL[2]+15.0*f[2])*gL[2]+((-15.0*fL[2])-15.0*f[2])*g[2]+(45.0*fL[1]-15.0*f[1]+41.56921938165305*fL[0]+10.39230484541326*f[0])*gL[1]+(15.0*fL[1]-45.0*f[1]+10.39230484541326*fL[0]+41.56921938165305*f[0])*g[1]+(17.32050807568877*gL[0]-17.32050807568877*g[0])*fL[1]+(17.32050807568877*g[0]-17.32050807568877*gL[0])*f[1]+(15.0*fL[0]+15.0*f[0])*gL[0]+((-15.0*fL[0])-15.0*f[0])*g[0])*gkyl_ipow(dv[0],-4)*dt;
    fOut[1] += 0.07216878364870323*((111.0*fL[3]-21.0*f[3]+103.9230484541326*fL[2]+10.39230484541326*f[2])*gL[3]+(69.0*fL[3]-159.0*f[3]+51.96152422706631*fL[2]+145.4922678357857*f[2])*g[3]+(51.96152422706631*gL[2]-51.96152422706631*g[2])*fL[3]+(51.96152422706631*g[2]-51.96152422706631*gL[2])*f[3]+(45.0*fL[2]+45.0*f[2])*gL[2]+((-45.0*fL[2])-45.0*f[2])*g[2]+(111.0*fL[1]-21.0*f[1]+103.9230484541326*fL[0]+10.39230484541326*f[0])*gL[1]+(69.0*fL[1]-159.0*f[1]+51.96152422706631*fL[0]+145.4922678357857*f[0])*g[1]+(51.96152422706631*gL[0]-51.96152422706631*g[0])*fL[1]+(51.96152422706631*g[0]-51.96152422706631*gL[0])*f[1]+(45.0*fL[0]+45.0*f[0])*gL[0]+((-45.0*fL[0])-45.0*f[0])*g[0])*gkyl_ipow(dv[0],-4)*dt;
    fOut[2] += -0.125*((45.0*fL[1]-15.0*f[1]+41.56921938165305*fL[0]+10.39230484541326*f[0])*gL[3]+(15.0*fL[1]-45.0*f[1]+10.39230484541326*fL[0]+41.56921938165305*f[0])*g[3]+(45.0*gL[1]+15.0*g[1]+17.32050807568877*gL[0]-17.32050807568877*g[0])*fL[3]+((-15.0*gL[1])-45.0*g[1]-17.32050807568877*gL[0]+17.32050807568877*g[0])*f[3]+(17.32050807568877*fL[1]-17.32050807568877*f[1]+15.0*fL[0]+15.0*f[0])*gL[2]+((-17.32050807568877*fL[1])+17.32050807568877*f[1]-15.0*fL[0]-15.0*f[0])*g[2]+(41.56921938165305*gL[1]+10.39230484541326*g[1]+15.0*gL[0]-15.0*g[0])*fL[2]+(10.39230484541326*gL[1]+41.56921938165305*g[1]+15.0*gL[0]-15.0*g[0])*f[2])*gkyl_ipow(dv[0],-4)*dt;
    fOut[3] += 0.125*((64.08587988004845*fL[1]-12.12435565298214*f[1]+60.0*fL[0]+6.0*f[0])*gL[3]+(39.83716857408417*fL[1]-91.7986928011505*f[1]+30.0*fL[0]+84.0*f[0])*g[3]+(64.08587988004845*gL[1]+39.83716857408417*g[1]+30.0*gL[0]-30.0*g[0])*fL[3]+((-12.12435565298214*gL[1])-91.7986928011505*g[1]-30.0*gL[0]+30.0*g[0])*f[3]+(30.0*fL[1]-30.0*f[1]+25.98076211353316*fL[0]+25.98076211353316*f[0])*gL[2]+((-30.0*fL[1])+30.0*f[1]-25.98076211353316*fL[0]-25.98076211353316*f[0])*g[2]+(60.0*gL[1]+30.0*g[1]+25.98076211353316*gL[0]-25.98076211353316*g[0])*fL[2]+(6.0*gL[1]+84.0*g[1]+25.98076211353316*gL[0]-25.98076211353316*g[0])*f[2])*gkyl_ipow(dv[0],-4)*dt;
  } else {
    fOut[0] += -0.125*((45.0*fR[3]-15.0*f[3]-41.56921938165305*fR[2]-10.39230484541326*f[2])*gR[3]+(45.0*fL[3]-15.0*f[3]+41.56921938165305*fL[2]+10.39230484541326*f[2])*gL[3]+(15.0*fR[3]+15.0*fL[3]-90.0*f[3]-10.39230484541326*fR[2]+10.39230484541326*fL[2])*g[3]+(17.32050807568877*g[2]-17.32050807568877*gR[2])*fR[3]+(17.32050807568877*gL[2]-17.32050807568877*g[2])*fL[3]+(17.32050807568877*gR[2]-17.32050807568877*gL[2])*f[3]+(15.0*fR[2]+15.0*f[2])*gR[2]+(15.0*fL[2]+15.0*f[2])*gL[2]+((-15.0*fR[2])-15.0*fL[2]-30.0*f[2])*g[2]+(45.0*fR[1]-15.0*f[1]-41.56921938165305*fR[0]-10.39230484541326*f[0])*gR[1]+(45.0*fL[1]-15.0*f[1]+41.56921938165305*fL[0]+10.39230484541326*f[0])*gL[1]+(15.0*fR[1]+15.0*fL[1]-90.0*f[1]-10.39230484541326*fR[0]+10.39230484541326*fL[0])*g[1]+(17.32050807568877*g[0]-17.32050807568877*gR[0])*fR[1]+(17.32050807568877*gL[0]-17.32050807568877*g[0])*fL[1]+(17.32050807568877*gR[0]-17.32050807568877*gL[0])*f[1]+(15.0*fR[0]+15.0*f[0])*gR[0]+(15.0*fL[0]+15.0*f[0])*gL[0]+((-15.0*fR[0])-15.0*fL[0]-30.0*f[0])*g[0])*gkyl_ipow(dv[0],-4)*dt;
    fOut[1] += -0.125*((64.08587988004845*fR[3]-12.12435565298214*f[3]-60.0*fR[2]-6.0*f[2])*gR[3]+((-64.08587988004845*fL[3])+12.12435565298214*f[3]-60.0*fL[2]-6.0*f[2])*gL[3]+(39.83716857408417*fR[3]-39.83716857408417*fL[3]-30.0*fR[2]-30.0*fL[2]-168.0*f[2])*g[3]+(30.0*g[2]-30.0*gR[2])*fR[3]+(30.0*g[2]-30.0*gL[2])*fL[3]+(30.0*gR[2]+30.0*gL[2]-60.0*g[2])*f[3]+(25.98076211353316*fR[2]+25.98076211353316*f[2])*gR[2]+((-25.98076211353316*fL[2])-25.98076211353316*f[2])*gL[2]+(25.98076211353316*fL[2]-25.98076211353316*fR[2])*g[2]+(64.08587988004845*fR[1]-12.12435565298214*f[1]-60.0*fR[0]-6.0*f[0])*gR[1]+((-64.08587988004845*fL[1])+12.12435565298214*f[1]-60.0*fL[0]-6.0*f[0])*gL[1]+(39.83716857408417*fR[1]-39.83716857408417*fL[1]-30.0*fR[0]-30.0*fL[0]-168.0*f[0])*g[1]+(30.0*g[0]-30.0*gR[0])*fR[1]+(30.0*g[0]-30.0*gL[0])*fL[1]+(30.0*gR[0]+30.0*gL[0]-60.0*g[0])*f[1]+(25.98076211353316*fR[0]+25.98076211353316*f[0])*gR[0]+((-25.98076211353316*fL[0])-25.98076211353316*f[0])*gL[0]+(25.98076211353316*fL[0]-25.98076211353316*fR[0])*g[0])*gkyl_ipow(dv[0],-4)*dt;
    fOut[2] += -0.125*((45.0*fR[1]-15.0*f[1]-41.56921938165305*fR[0]-10.39230484541326*f[0])*gR[3]+(45.0*fL[1]-15.0*f[1]+41.56921938165305*fL[0]+10.39230484541326*f[0])*gL[3]+(15.0*fR[1]+15.0*fL[1]-90.0*f[1]-10.39230484541326*fR[0]+10.39230484541326*fL[0])*g[3]+(45.0*gR[1]+15.0*g[1]-17.32050807568877*gR[0]+17.32050807568877*g[0])*fR[3]+(45.0*gL[1]+15.0*g[1]+17.32050807568877*gL[0]-17.32050807568877*g[0])*fL[3]+((-15.0*gR[1])-15.0*gL[1]-90.0*g[1]+17.32050807568877*gR[0]-17.32050807568877*gL[0])*f[3]+((-17.32050807568877*fR[1])+17.32050807568877*f[1]+15.0*fR[0]+15.0*f[0])*gR[2]+(17.32050807568877*fL[1]-17.32050807568877*f[1]+15.0*fL[0]+15.0*f[0])*gL[2]+(17.32050807568877*fR[1]-17.32050807568877*fL[1]-15.0*fR[0]-15.0*fL[0]-30.0*f[0])*g[2]+((-41.56921938165305*gR[1])-10.39230484541326*g[1]+15.0*gR[0]-15.0*g[0])*fR[2]+(41.56921938165305*gL[1]+10.39230484541326*g[1]+15.0*gL[0]-15.0*g[0])*fL[2]+((-10.39230484541326*gR[1])+10.39230484541326*gL[1]+15.0*gR[0]+15.0*gL[0]-30.0*g[0])*f[2])*gkyl_ipow(dv[0],-4)*dt;
    fOut[3] += -0.125*((64.08587988004845*fR[1]-12.12435565298214*f[1]-60.0*fR[0]-6.0*f[0])*gR[3]+((-64.08587988004845*fL[1])+12.12435565298214*f[1]-60.0*fL[0]-6.0*f[0])*gL[3]+(39.83716857408417*fR[1]-39.83716857408417*fL[1]-30.0*fR[0]-30.0*fL[0]-168.0*f[0])*g[3]+(64.08587988004845*gR[1]+39.83716857408417*g[1]-30.0*gR[0]+30.0*g[0])*fR[3]+((-64.08587988004845*gL[1])-39.83716857408417*g[1]-30.0*gL[0]+30.0*g[0])*fL[3]+((-12.12435565298214*gR[1])+12.12435565298214*gL[1]+30.0*gR[0]+30.0*gL[0]-60.0*g[0])*f[3]+((-30.0*fR[1])+30.0*f[1]+25.98076211353316*fR[0]+25.98076211353316*f[0])*gR[2]+((-30.0*fL[1])+30.0*f[1]-25.98076211353316*fL[0]-25.98076211353316*f[0])*gL[2]+(30.0*fR[1]+30.0*fL[1]-60.0*f[1]-25.98076211353316*fR[0]+25.98076211353316*fL[0])*g[2]+((-60.0*gR[1])-30.0*g[1]+25.98076211353316*gR[0]-25.98076211353316*g[0])*fR[2]+((-60.0*gL[1])-30.0*g[1]-25.98076211353316*gL[0]+25.98076211353316*g[0])*fL[2]+((-6.0*gR[1])-6.0*gL[1]-168.0*g[1]+25.98076211353316*gR[0]-25.98076211353316*gL[0])*f[2])*gkyl_ipow(dv[0],-4)*dt;
  }

  fOut[0] += 0.0;
  fOut[1] += 0.0;
  fOut[2] += 0.0;
  fOut[3] += 0.0;

  if (isLeftEdge) {
    fOut[0] += 0.125*((30.0*fR[3]-30.0*f[3]-25.98076211353316*fR[2]-25.98076211353316*f[2])*gR[3]+(30.0*fR[3]-30.0*f[3]-25.98076211353316*fR[2]-25.98076211353316*f[2])*g[3]+(31.17691453623978*g[2]-31.17691453623978*gR[2])*fR[3]+(31.17691453623978*gR[2]-31.17691453623978*g[2])*f[3]+(27.0*fR[2]+27.0*f[2])*gR[2]+((-27.0*fR[2])-27.0*f[2])*g[2])*gkyl_ipow(dv[0],-4)*dt;
    fOut[1] += 0.07216878364870323*((90.0*fR[3]-90.0*f[3]-77.94228634059945*fR[2]-77.94228634059945*f[2])*gR[3]+(90.0*fR[3]-90.0*f[3]-77.94228634059945*fR[2]-77.94228634059945*f[2])*g[3]+(93.53074360871933*g[2]-93.53074360871933*gR[2])*fR[3]+(93.53074360871933*gR[2]-93.53074360871933*g[2])*f[3]+(81.0*fR[2]+81.0*f[2])*gR[2]+((-81.0*fR[2])-81.0*f[2])*g[2])*gkyl_ipow(dv[0],-4)*dt;
    fOut[2] += -0.07216878364870323*((51.96152422706631*fR[1]-51.96152422706631*f[1]-45.0*fR[0]-45.0*f[0])*gR[3]+(51.96152422706631*fR[1]-51.96152422706631*f[1]-45.0*fR[0]-45.0*f[0])*g[3]+((-54.0*fR[1])+54.0*f[1]+46.76537180435967*fR[0]+46.76537180435967*f[0])*gR[2]+(54.0*fR[1]-54.0*f[1]-46.76537180435967*fR[0]-46.76537180435967*f[0])*g[2])*gkyl_ipow(dv[0],-4)*dt;
    fOut[3] += -0.125*((51.96152422706631*fR[1]-51.96152422706631*f[1]-45.0*fR[0]-45.0*f[0])*gR[3]+(51.96152422706631*fR[1]-51.96152422706631*f[1]-45.0*fR[0]-45.0*f[0])*g[3]+((-54.0*fR[1])+54.0*f[1]+46.76537180435967*fR[0]+46.76537180435967*f[0])*gR[2]+(54.0*fR[1]-54.0*f[1]-46.76537180435967*fR[0]-46.76537180435967*f[0])*g[2])*gkyl_ipow(dv[0],-4)*dt;
  } else if (isRightEdge) {
    fOut[0] += 0.125*((30.0*fL[3]-30.0*f[3]+25.98076211353316*fL[2]+25.98076211353316*f[2])*gL[3]+(30.0*fL[3]-30.0*f[3]+25.98076211353316*fL[2]+25.98076211353316*f[2])*g[3]+(31.17691453623978*gL[2]-31.17691453623978*g[2])*fL[3]+(31.17691453623978*g[2]-31.17691453623978*gL[2])*f[3]+(27.0*fL[2]+27.0*f[2])*gL[2]+((-27.0*fL[2])-27.0*f[2])*g[2])*gkyl_ipow(dv[0],-4)*dt;
    fOut[1] += -0.07216878364870323*((90.0*fL[3]-90.0*f[3]+77.94228634059945*fL[2]+77.94228634059945*f[2])*gL[3]+(90.0*fL[3]-90.0*f[3]+77.94228634059945*fL[2]+77.94228634059945*f[2])*g[3]+(93.53074360871933*gL[2]-93.53074360871933*g[2])*fL[3]+(93.53074360871933*g[2]-93.53074360871933*gL[2])*f[3]+(81.0*fL[2]+81.0*f[2])*gL[2]+((-81.0*fL[2])-81.0*f[2])*g[2])*gkyl_ipow(dv[0],-4)*dt;
    fOut[2] += -0.07216878364870323*((51.96152422706631*fL[1]-51.96152422706631*f[1]+45.0*fL[0]+45.0*f[0])*gL[3]+(51.96152422706631*fL[1]-51.96152422706631*f[1]+45.0*fL[0]+45.0*f[0])*g[3]+(54.0*fL[1]-54.0*f[1]+46.76537180435967*fL[0]+46.76537180435967*f[0])*gL[2]+((-54.0*fL[1])+54.0*f[1]-46.76537180435967*fL[0]-46.76537180435967*f[0])*g[2])*gkyl_ipow(dv[0],-4)*dt;
    fOut[3] += 0.125*((51.96152422706631*fL[1]-51.96152422706631*f[1]+45.0*fL[0]+45.0*f[0])*gL[3]+(51.96152422706631*fL[1]-51.96152422706631*f[1]+45.0*fL[0]+45.0*f[0])*g[3]+(54.0*fL[1]-54.0*f[1]+46.76537180435967*fL[0]+46.76537180435967*f[0])*gL[2]+((-54.0*fL[1])+54.0*f[1]-46.76537180435967*fL[0]-46.76537180435967*f[0])*g[2])*gkyl_ipow(dv[0],-4)*dt;
  } else {
    fOut[0] += 0.125*((30.0*fR[3]-30.0*f[3]-25.98076211353316*fR[2]-25.98076211353316*f[2])*gR[3]+(30.0*fL[3]-30.0*f[3]+25.98076211353316*fL[2]+25.98076211353316*f[2])*gL[3]+(30.0*fR[3]+30.0*fL[3]-60.0*f[3]-25.98076211353316*fR[2]+25.98076211353316*fL[2])*g[3]+(31.17691453623978*g[2]-31.17691453623978*gR[2])*fR[3]+(31.17691453623978*gL[2]-31.17691453623978*g[2])*fL[3]+(31.17691453623978*gR[2]-31.17691453623978*gL[2])*f[3]+(27.0*fR[2]+27.0*f[2])*gR[2]+(27.0*fL[2]+27.0*f[2])*gL[2]+((-27.0*fR[2])-27.0*fL[2]-54.0*f[2])*g[2])*gkyl_ipow(dv[0],-4)*dt;
    fOut[1] += 0.125*((51.96152422706631*fR[3]-51.96152422706631*f[3]-45.0*fR[2]-45.0*f[2])*gR[3]+((-51.96152422706631*fL[3])+51.96152422706631*f[3]-45.0*fL[2]-45.0*f[2])*gL[3]+(51.96152422706631*fR[3]-51.96152422706631*fL[3]-45.0*fR[2]-45.0*fL[2]-90.0*f[2])*g[3]+(54.0*g[2]-54.0*gR[2])*fR[3]+(54.0*g[2]-54.0*gL[2])*fL[3]+(54.0*gR[2]+54.0*gL[2]-108.0*g[2])*f[3]+(46.76537180435967*fR[2]+46.76537180435967*f[2])*gR[2]+((-46.76537180435967*fL[2])-46.76537180435967*f[2])*gL[2]+(46.76537180435967*fL[2]-46.76537180435967*fR[2])*g[2])*gkyl_ipow(dv[0],-4)*dt;
    fOut[2] += -0.125*((30.0*fR[1]-30.0*f[1]-25.98076211353316*fR[0]-25.98076211353316*f[0])*gR[3]+(30.0*fL[1]-30.0*f[1]+25.98076211353316*fL[0]+25.98076211353316*f[0])*gL[3]+(30.0*fR[1]+30.0*fL[1]-60.0*f[1]-25.98076211353316*fR[0]+25.98076211353316*fL[0])*g[3]+((-31.17691453623978*fR[1])+31.17691453623978*f[1]+27.0*fR[0]+27.0*f[0])*gR[2]+(31.17691453623978*fL[1]-31.17691453623978*f[1]+27.0*fL[0]+27.0*f[0])*gL[2]+(31.17691453623978*fR[1]-31.17691453623978*fL[1]-27.0*fR[0]-27.0*fL[0]-54.0*f[0])*g[2])*gkyl_ipow(dv[0],-4)*dt;
    fOut[3] += -0.125*((51.96152422706631*fR[1]-51.96152422706631*f[1]-45.0*fR[0]-45.0*f[0])*gR[3]+((-51.96152422706631*fL[1])+51.96152422706631*f[1]-45.0*fL[0]-45.0*f[0])*gL[3]+(51.96152422706631*fR[1]-51.96152422706631*fL[1]-45.0*fR[0]-45.0*fL[0]-90.0*f[0])*g[3]+((-54.0*fR[1])+54.0*f[1]+46.76537180435967*fR[0]+46.76537180435967*f[0])*gR[2]+((-54.0*fL[1])+54.0*f[1]-46.76537180435967*fL[0]-46.76537180435967*f[0])*gL[2]+(54.0*fR[1]+54.0*fL[1]-108.0*f[1]-46.76537180435967*fR[0]+46.76537180435967*fL[0])*g[2])*gkyl_ipow(dv[0],-4)*dt;
  }

  fOut[0] += 0.0;
  fOut[1] += 0.0;
  fOut[2] += 0.0;
  fOut[3] += -1.0*((10.39230484541326*f[1]+6.0*f[0])*gR[3]+(6.0*f[0]-10.39230484541326*f[1])*gL[3]-12.0*f[0]*g[3]+((-9.0*f[1])-5.196152422706631*f[0])*gR[2]+(5.196152422706631*f[0]-9.0*f[1])*gL[2]+18.0*f[1]*g[2])*gkyl_ipow(dv[0],-4)*dt;

  if (isBotEdge) {
    fOut[0] += -0.125*((45.0*fT[3]-15.0*f[3]-41.56921938165305*fT[1]-10.39230484541326*f[1])*gT[3]+(15.0*fT[3]-45.0*f[3]-10.39230484541326*fT[1]-41.56921938165305*f[1])*g[3]+(17.32050807568877*g[1]-17.32050807568877*gT[1])*fT[3]+(17.32050807568877*gT[1]-17.32050807568877*g[1])*f[3]+(45.0*fT[2]-15.0*f[2]-41.56921938165305*fT[0]-10.39230484541326*f[0])*gT[2]+(15.0*fT[2]-45.0*f[2]-10.39230484541326*fT[0]-41.56921938165305*f[0])*g[2]+(17.32050807568877*g[0]-17.32050807568877*gT[0])*fT[2]+(17.32050807568877*gT[0]-17.32050807568877*g[0])*f[2]+(15.0*fT[1]+15.0*f[1])*gT[1]+((-15.0*fT[1])-15.0*f[1])*g[1]+(15.0*fT[0]+15.0*f[0])*gT[0]+((-15.0*fT[0])-15.0*f[0])*g[0])*gkyl_ipow(dv[1],-4)*dt;
    fOut[1] += -0.125*((45.0*fT[2]-15.0*f[2]-41.56921938165305*fT[0]-10.39230484541326*f[0])*gT[3]+(15.0*fT[2]-45.0*f[2]-10.39230484541326*fT[0]-41.56921938165305*f[0])*g[3]+(45.0*gT[2]+15.0*g[2]-17.32050807568877*gT[0]+17.32050807568877*g[0])*fT[3]+((-15.0*gT[2])-45.0*g[2]+17.32050807568877*gT[0]-17.32050807568877*g[0])*f[3]+((-41.56921938165305*fT[1])-10.39230484541326*f[1])*gT[2]+((-10.39230484541326*fT[1])-41.56921938165305*f[1])*g[2]+(17.32050807568877*g[1]-17.32050807568877*gT[1])*fT[2]+(17.32050807568877*gT[1]-17.32050807568877*g[1])*f[2]+(15.0*fT[0]+15.0*f[0])*gT[1]+((-15.0*fT[0])-15.0*f[0])*g[1]+(15.0*gT[0]-15.0*g[0])*fT[1]+(15.0*gT[0]-15.0*g[0])*f[1])*gkyl_ipow(dv[1],-4)*dt;
    fOut[2] += -0.07216878364870323*((111.0*fT[3]-21.0*f[3]-103.9230484541326*fT[1]-10.39230484541326*f[1])*gT[3]+(69.0*fT[3]-159.0*f[3]-51.96152422706631*fT[1]-145.4922678357857*f[1])*g[3]+(51.96152422706631*g[1]-51.96152422706631*gT[1])*fT[3]+(51.96152422706631*gT[1]-51.96152422706631*g[1])*f[3]+(111.0*fT[2]-21.0*f[2]-103.9230484541326*fT[0]-10.39230484541326*f[0])*gT[2]+(69.0*fT[2]-159.0*f[2]-51.96152422706631*fT[0]-145.4922678357857*f[0])*g[2]+(51.96152422706631*g[0]-51.96152422706631*gT[0])*fT[2]+(51.96152422706631*gT[0]-51.96152422706631*g[0])*f[2]+(45.0*fT[1]+45.0*f[1])*gT[1]+((-45.0*fT[1])-45.0*f[1])*g[1]+(45.0*fT[0]+45.0*f[0])*gT[0]+((-45.0*fT[0])-45.0*f[0])*g[0])*gkyl_ipow(dv[1],-4)*dt;
    fOut[3] += -0.125*((64.08587988004845*fT[2]-12.12435565298214*f[2]-60.0*fT[0]-6.0*f[0])*gT[3]+(39.83716857408417*fT[2]-91.7986928011505*f[2]-30.0*fT[0]-84.0*f[0])*g[3]+(64.08587988004845*gT[2]+39.83716857408417*g[2]-30.0*gT[0]+30.0*g[0])*fT[3]+((-12.12435565298214*gT[2])-91.7986928011505*g[2]+30.0*gT[0]-30.0*g[0])*f[3]+((-60.0*fT[1])-6.0*f[1])*gT[2]+((-30.0*fT[1])-84.0*f[1])*g[2]+(30.0*g[1]-30.0*gT[1])*fT[2]+(30.0*gT[1]-30.0*g[1])*f[2]+(25.98076211353316*fT[0]+25.98076211353316*f[0])*gT[1]+((-25.98076211353316*fT[0])-25.98076211353316*f[0])*g[1]+(25.98076211353316*gT[0]-25.98076211353316*g[0])*fT[1]+(25.98076211353316*gT[0]-25.98076211353316*g[0])*f[1])*gkyl_ipow(dv[1],-4)*dt;
  } else if (isTopEdge) {
    fOut[0] += -0.125*((45.0*fB[3]-15.0*f[3]+41.56921938165305*fB[1]+10.39230484541326*f[1])*gB[3]+(15.0*fB[3]-45.0*f[3]+10.39230484541326*fB[1]+41.56921938165305*f[1])*g[3]+(17.32050807568877*gB[1]-17.32050807568877*g[1])*fB[3]+(17.32050807568877*g[1]-17.32050807568877*gB[1])*f[3]+(45.0*fB[2]-15.0*f[2]+41.56921938165305*fB[0]+10.39230484541326*f[0])*gB[2]+(15.0*fB[2]-45.0*f[2]+10.39230484541326*fB[0]+41.56921938165305*f[0])*g[2]+(17.32050807568877*gB[0]-17.32050807568877*g[0])*fB[2]+(17.32050807568877*g[0]-17.32050807568877*gB[0])*f[2]+(15.0*fB[1]+15.0*f[1])*gB[1]+((-15.0*fB[1])-15.0*f[1])*g[1]+(15.0*fB[0]+15.0*f[0])*gB[0]+((-15.0*fB[0])-15.0*f[0])*g[0])*gkyl_ipow(dv[1],-4)*dt;
    fOut[1] += -0.125*((45.0*fB[2]-15.0*f[2]+41.56921938165305*fB[0]+10.39230484541326*f[0])*gB[3]+(15.0*fB[2]-45.0*f[2]+10.39230484541326*fB[0]+41.56921938165305*f[0])*g[3]+(45.0*gB[2]+15.0*g[2]+17.32050807568877*gB[0]-17.32050807568877*g[0])*fB[3]+((-15.0*gB[2])-45.0*g[2]-17.32050807568877*gB[0]+17.32050807568877*g[0])*f[3]+(41.56921938165305*fB[1]+10.39230484541326*f[1])*gB[2]+(10.39230484541326*fB[1]+41.56921938165305*f[1])*g[2]+(17.32050807568877*gB[1]-17.32050807568877*g[1])*fB[2]+(17.32050807568877*g[1]-17.32050807568877*gB[1])*f[2]+(15.0*fB[0]+15.0*f[0])*gB[1]+((-15.0*fB[0])-15.0*f[0])*g[1]+(15.0*gB[0]-15.0*g[0])*fB[1]+(15.0*gB[0]-15.0*g[0])*f[1])*gkyl_ipow(dv[1],-4)*dt;
    fOut[2] += 0.07216878364870323*((111.0*fB[3]-21.0*f[3]+103.9230484541326*fB[1]+10.39230484541326*f[1])*gB[3]+(69.0*fB[3]-159.0*f[3]+51.96152422706631*fB[1]+145.4922678357857*f[1])*g[3]+(51.96152422706631*gB[1]-51.96152422706631*g[1])*fB[3]+(51.96152422706631*g[1]-51.96152422706631*gB[1])*f[3]+(111.0*fB[2]-21.0*f[2]+103.9230484541326*fB[0]+10.39230484541326*f[0])*gB[2]+(69.0*fB[2]-159.0*f[2]+51.96152422706631*fB[0]+145.4922678357857*f[0])*g[2]+(51.96152422706631*gB[0]-51.96152422706631*g[0])*fB[2]+(51.96152422706631*g[0]-51.96152422706631*gB[0])*f[2]+(45.0*fB[1]+45.0*f[1])*gB[1]+((-45.0*fB[1])-45.0*f[1])*g[1]+(45.0*fB[0]+45.0*f[0])*gB[0]+((-45.0*fB[0])-45.0*f[0])*g[0])*gkyl_ipow(dv[1],-4)*dt;
    fOut[3] += 0.125*((64.08587988004845*fB[2]-12.12435565298214*f[2]+60.0*fB[0]+6.0*f[0])*gB[3]+(39.83716857408417*fB[2]-91.7986928011505*f[2]+30.0*fB[0]+84.0*f[0])*g[3]+(64.08587988004845*gB[2]+39.83716857408417*g[2]+30.0*gB[0]-30.0*g[0])*fB[3]+((-12.12435565298214*gB[2])-91.7986928011505*g[2]-30.0*gB[0]+30.0*g[0])*f[3]+(60.0*fB[1]+6.0*f[1])*gB[2]+(30.0*fB[1]+84.0*f[1])*g[2]+(30.0*gB[1]-30.0*g[1])*fB[2]+(30.0*g[1]-30.0*gB[1])*f[2]+(25.98076211353316*fB[0]+25.98076211353316*f[0])*gB[1]+((-25.98076211353316*fB[0])-25.98076211353316*f[0])*g[1]+(25.98076211353316*gB[0]-25.98076211353316*g[0])*fB[1]+(25.98076211353316*gB[0]-25.98076211353316*g[0])*f[1])*gkyl_ipow(dv[1],-4)*dt;
  } else {
    fOut[0] += -0.125*((45.0*fT[3]-15.0*f[3]-41.56921938165305*fT[1]-10.39230484541326*f[1])*gT[3]+(45.0*fB[3]-15.0*f[3]+41.56921938165305*fB[1]+10.39230484541326*f[1])*gB[3]+(15.0*fT[3]+15.0*fB[3]-90.0*f[3]-10.39230484541326*fT[1]+10.39230484541326*fB[1])*g[3]+(17.32050807568877*g[1]-17.32050807568877*gT[1])*fT[3]+(17.32050807568877*gB[1]-17.32050807568877*g[1])*fB[3]+(17.32050807568877*gT[1]-17.32050807568877*gB[1])*f[3]+(45.0*fT[2]-15.0*f[2]-41.56921938165305*fT[0]-10.39230484541326*f[0])*gT[2]+(45.0*fB[2]-15.0*f[2]+41.56921938165305*fB[0]+10.39230484541326*f[0])*gB[2]+(15.0*fT[2]+15.0*fB[2]-90.0*f[2]-10.39230484541326*fT[0]+10.39230484541326*fB[0])*g[2]+(17.32050807568877*g[0]-17.32050807568877*gT[0])*fT[2]+(17.32050807568877*gB[0]-17.32050807568877*g[0])*fB[2]+(17.32050807568877*gT[0]-17.32050807568877*gB[0])*f[2]+(15.0*fT[1]+15.0*f[1])*gT[1]+(15.0*fB[1]+15.0*f[1])*gB[1]+((-15.0*fT[1])-15.0*fB[1]-30.0*f[1])*g[1]+(15.0*fT[0]+15.0*f[0])*gT[0]+(15.0*fB[0]+15.0*f[0])*gB[0]+((-15.0*fT[0])-15.0*fB[0]-30.0*f[0])*g[0])*gkyl_ipow(dv[1],-4)*dt;
    fOut[1] += -0.125*((45.0*fT[2]-15.0*f[2]-41.56921938165305*fT[0]-10.39230484541326*f[0])*gT[3]+(45.0*fB[2]-15.0*f[2]+41.56921938165305*fB[0]+10.39230484541326*f[0])*gB[3]+(15.0*fT[2]+15.0*fB[2]-90.0*f[2]-10.39230484541326*fT[0]+10.39230484541326*fB[0])*g[3]+(45.0*gT[2]+15.0*g[2]-17.32050807568877*gT[0]+17.32050807568877*g[0])*fT[3]+(45.0*gB[2]+15.0*g[2]+17.32050807568877*gB[0]-17.32050807568877*g[0])*fB[3]+((-15.0*gT[2])-15.0*gB[2]-90.0*g[2]+17.32050807568877*gT[0]-17.32050807568877*gB[0])*f[3]+((-41.56921938165305*fT[1])-10.39230484541326*f[1])*gT[2]+(41.56921938165305*fB[1]+10.39230484541326*f[1])*gB[2]+(10.39230484541326*fB[1]-10.39230484541326*fT[1])*g[2]+(17.32050807568877*g[1]-17.32050807568877*gT[1])*fT[2]+(17.32050807568877*gB[1]-17.32050807568877*g[1])*fB[2]+(17.32050807568877*gT[1]-17.32050807568877*gB[1])*f[2]+(15.0*fT[0]+15.0*f[0])*gT[1]+(15.0*fB[0]+15.0*f[0])*gB[1]+((-15.0*fT[0])-15.0*fB[0]-30.0*f[0])*g[1]+(15.0*gT[0]-15.0*g[0])*fT[1]+(15.0*gB[0]-15.0*g[0])*fB[1]+(15.0*gT[0]+15.0*gB[0]-30.0*g[0])*f[1])*gkyl_ipow(dv[1],-4)*dt;
    fOut[2] += -0.125*((64.08587988004845*fT[3]-12.12435565298214*f[3]-60.0*fT[1]-6.0*f[1])*gT[3]+((-64.08587988004845*fB[3])+12.12435565298214*f[3]-60.0*fB[1]-6.0*f[1])*gB[3]+(39.83716857408417*fT[3]-39.83716857408417*fB[3]-30.0*fT[1]-30.0*fB[1]-168.0*f[1])*g[3]+(30.0*g[1]-30.0*gT[1])*fT[3]+(30.0*g[1]-30.0*gB[1])*fB[3]+(30.0*gT[1]+30.0*gB[1]-60.0*g[1])*f[3]+(64.08587988004845*fT[2]-12.12435565298214*f[2]-60.0*fT[0]-6.0*f[0])*gT[2]+((-64.08587988004845*fB[2])+12.12435565298214*f[2]-60.0*fB[0]-6.0*f[0])*gB[2]+(39.83716857408417*fT[2]-39.83716857408417*fB[2]-30.0*fT[0]-30.0*fB[0]-168.0*f[0])*g[2]+(30.0*g[0]-30.0*gT[0])*fT[2]+(30.0*g[0]-30.0*gB[0])*fB[2]+(30.0*gT[0]+30.0*gB[0]-60.0*g[0])*f[2]+(25.98076211353316*fT[1]+25.98076211353316*f[1])*gT[1]+((-25.98076211353316*fB[1])-25.98076211353316*f[1])*gB[1]+(25.98076211353316*fB[1]-25.98076211353316*fT[1])*g[1]+(25.98076211353316*fT[0]+25.98076211353316*f[0])*gT[0]+((-25.98076211353316*fB[0])-25.98076211353316*f[0])*gB[0]+(25.98076211353316*fB[0]-25.98076211353316*fT[0])*g[0])*gkyl_ipow(dv[1],-4)*dt;
    fOut[3] += -0.125*((64.08587988004845*fT[2]-12.12435565298214*f[2]-60.0*fT[0]-6.0*f[0])*gT[3]+((-64.08587988004845*fB[2])+12.12435565298214*f[2]-60.0*fB[0]-6.0*f[0])*gB[3]+(39.83716857408417*fT[2]-39.83716857408417*fB[2]-30.0*fT[0]-30.0*fB[0]-168.0*f[0])*g[3]+(64.08587988004845*gT[2]+39.83716857408417*g[2]-30.0*gT[0]+30.0*g[0])*fT[3]+((-64.08587988004845*gB[2])-39.83716857408417*g[2]-30.0*gB[0]+30.0*g[0])*fB[3]+((-12.12435565298214*gT[2])+12.12435565298214*gB[2]+30.0*gT[0]+30.0*gB[0]-60.0*g[0])*f[3]+((-60.0*fT[1])-6.0*f[1])*gT[2]+((-60.0*fB[1])-6.0*f[1])*gB[2]+((-30.0*fT[1])-30.0*fB[1]-168.0*f[1])*g[2]+(30.0*g[1]-30.0*gT[1])*fT[2]+(30.0*g[1]-30.0*gB[1])*fB[2]+(30.0*gT[1]+30.0*gB[1]-60.0*g[1])*f[2]+(25.98076211353316*fT[0]+25.98076211353316*f[0])*gT[1]+((-25.98076211353316*fB[0])-25.98076211353316*f[0])*gB[1]+(25.98076211353316*fB[0]-25.98076211353316*fT[0])*g[1]+(25.98076211353316*gT[0]-25.98076211353316*g[0])*fT[1]+(25.98076211353316*g[0]-25.98076211353316*gB[0])*fB[1]+(25.98076211353316*gT[0]-25.98076211353316*gB[0])*f[1])*gkyl_ipow(dv[1],-4)*dt;
  }

  fOut[0] += 0.0;
  fOut[1] += 0.0;
  fOut[2] += 0.0;
  fOut[3] += 0.0;
  if (isBotEdge) {
    fOut[0] += 0.125*((30.0*fT[3]-30.0*f[3]-25.98076211353316*fT[1]-25.98076211353316*f[1])*gT[3]+(30.0*fT[3]-30.0*f[3]-25.98076211353316*fT[1]-25.98076211353316*f[1])*g[3]+(31.17691453623978*g[1]-31.17691453623978*gT[1])*fT[3]+(31.17691453623978*gT[1]-31.17691453623978*g[1])*f[3]+(27.0*fT[1]+27.0*f[1])*gT[1]+((-27.0*fT[1])-27.0*f[1])*g[1])*gkyl_ipow(dv[1],-4)*dt;
    fOut[1] += -0.07216878364870323*((51.96152422706631*fT[2]-51.96152422706631*f[2]-45.0*fT[0]-45.0*f[0])*gT[3]+(51.96152422706631*fT[2]-51.96152422706631*f[2]-45.0*fT[0]-45.0*f[0])*g[3]+(54.0*g[1]-54.0*gT[1])*fT[2]+(54.0*gT[1]-54.0*g[1])*f[2]+(46.76537180435967*fT[0]+46.76537180435967*f[0])*gT[1]+((-46.76537180435967*fT[0])-46.76537180435967*f[0])*g[1])*gkyl_ipow(dv[1],-4)*dt;
    fOut[2] += 0.07216878364870323*((90.0*fT[3]-90.0*f[3]-77.94228634059945*fT[1]-77.94228634059945*f[1])*gT[3]+(90.0*fT[3]-90.0*f[3]-77.94228634059945*fT[1]-77.94228634059945*f[1])*g[3]+(93.53074360871933*g[1]-93.53074360871933*gT[1])*fT[3]+(93.53074360871933*gT[1]-93.53074360871933*g[1])*f[3]+(81.0*fT[1]+81.0*f[1])*gT[1]+((-81.0*fT[1])-81.0*f[1])*g[1])*gkyl_ipow(dv[1],-4)*dt;
    fOut[3] += -0.125*((51.96152422706631*fT[2]-51.96152422706631*f[2]-45.0*fT[0]-45.0*f[0])*gT[3]+(51.96152422706631*fT[2]-51.96152422706631*f[2]-45.0*fT[0]-45.0*f[0])*g[3]+(54.0*g[1]-54.0*gT[1])*fT[2]+(54.0*gT[1]-54.0*g[1])*f[2]+(46.76537180435967*fT[0]+46.76537180435967*f[0])*gT[1]+((-46.76537180435967*fT[0])-46.76537180435967*f[0])*g[1])*gkyl_ipow(dv[1],-4)*dt;
  } else if (isTopEdge) {
    fOut[0] += 0.125*((30.0*fB[3]-30.0*f[3]+25.98076211353316*fB[1]+25.98076211353316*f[1])*gB[3]+(30.0*fB[3]-30.0*f[3]+25.98076211353316*fB[1]+25.98076211353316*f[1])*g[3]+(31.17691453623978*gB[1]-31.17691453623978*g[1])*fB[3]+(31.17691453623978*g[1]-31.17691453623978*gB[1])*f[3]+(27.0*fB[1]+27.0*f[1])*gB[1]+((-27.0*fB[1])-27.0*f[1])*g[1])*gkyl_ipow(dv[1],-4)*dt;
    fOut[1] += -0.07216878364870323*((51.96152422706631*fB[2]-51.96152422706631*f[2]+45.0*fB[0]+45.0*f[0])*gB[3]+(51.96152422706631*fB[2]-51.96152422706631*f[2]+45.0*fB[0]+45.0*f[0])*g[3]+(54.0*gB[1]-54.0*g[1])*fB[2]+(54.0*g[1]-54.0*gB[1])*f[2]+(46.76537180435967*fB[0]+46.76537180435967*f[0])*gB[1]+((-46.76537180435967*fB[0])-46.76537180435967*f[0])*g[1])*gkyl_ipow(dv[1],-4)*dt;
    fOut[2] += -0.07216878364870323*((90.0*fB[3]-90.0*f[3]+77.94228634059945*fB[1]+77.94228634059945*f[1])*gB[3]+(90.0*fB[3]-90.0*f[3]+77.94228634059945*fB[1]+77.94228634059945*f[1])*g[3]+(93.53074360871933*gB[1]-93.53074360871933*g[1])*fB[3]+(93.53074360871933*g[1]-93.53074360871933*gB[1])*f[3]+(81.0*fB[1]+81.0*f[1])*gB[1]+((-81.0*fB[1])-81.0*f[1])*g[1])*gkyl_ipow(dv[1],-4)*dt;
    fOut[3] += 0.125*((51.96152422706631*fB[2]-51.96152422706631*f[2]+45.0*fB[0]+45.0*f[0])*gB[3]+(51.96152422706631*fB[2]-51.96152422706631*f[2]+45.0*fB[0]+45.0*f[0])*g[3]+(54.0*gB[1]-54.0*g[1])*fB[2]+(54.0*g[1]-54.0*gB[1])*f[2]+(46.76537180435967*fB[0]+46.76537180435967*f[0])*gB[1]+((-46.76537180435967*fB[0])-46.76537180435967*f[0])*g[1])*gkyl_ipow(dv[1],-4)*dt;
  } else {
    fOut[0] += 0.125*((30.0*fT[3]-30.0*f[3]-25.98076211353316*fT[1]-25.98076211353316*f[1])*gT[3]+(30.0*fB[3]-30.0*f[3]+25.98076211353316*fB[1]+25.98076211353316*f[1])*gB[3]+(30.0*fT[3]+30.0*fB[3]-60.0*f[3]-25.98076211353316*fT[1]+25.98076211353316*fB[1])*g[3]+(31.17691453623978*g[1]-31.17691453623978*gT[1])*fT[3]+(31.17691453623978*gB[1]-31.17691453623978*g[1])*fB[3]+(31.17691453623978*gT[1]-31.17691453623978*gB[1])*f[3]+(27.0*fT[1]+27.0*f[1])*gT[1]+(27.0*fB[1]+27.0*f[1])*gB[1]+((-27.0*fT[1])-27.0*fB[1]-54.0*f[1])*g[1])*gkyl_ipow(dv[1],-4)*dt;
    fOut[1] += -0.125*((30.0*fT[2]-30.0*f[2]-25.98076211353316*fT[0]-25.98076211353316*f[0])*gT[3]+(30.0*fB[2]-30.0*f[2]+25.98076211353316*fB[0]+25.98076211353316*f[0])*gB[3]+(30.0*fT[2]+30.0*fB[2]-60.0*f[2]-25.98076211353316*fT[0]+25.98076211353316*fB[0])*g[3]+(31.17691453623978*g[1]-31.17691453623978*gT[1])*fT[2]+(31.17691453623978*gB[1]-31.17691453623978*g[1])*fB[2]+(31.17691453623978*gT[1]-31.17691453623978*gB[1])*f[2]+(27.0*fT[0]+27.0*f[0])*gT[1]+(27.0*fB[0]+27.0*f[0])*gB[1]+((-27.0*fT[0])-27.0*fB[0]-54.0*f[0])*g[1])*gkyl_ipow(dv[1],-4)*dt;
    fOut[2] += 0.125*((51.96152422706631*fT[3]-51.96152422706631*f[3]-45.0*fT[1]-45.0*f[1])*gT[3]+((-51.96152422706631*fB[3])+51.96152422706631*f[3]-45.0*fB[1]-45.0*f[1])*gB[3]+(51.96152422706631*fT[3]-51.96152422706631*fB[3]-45.0*fT[1]-45.0*fB[1]-90.0*f[1])*g[3]+(54.0*g[1]-54.0*gT[1])*fT[3]+(54.0*g[1]-54.0*gB[1])*fB[3]+(54.0*gT[1]+54.0*gB[1]-108.0*g[1])*f[3]+(46.76537180435967*fT[1]+46.76537180435967*f[1])*gT[1]+((-46.76537180435967*fB[1])-46.76537180435967*f[1])*gB[1]+(46.76537180435967*fB[1]-46.76537180435967*fT[1])*g[1])*gkyl_ipow(dv[1],-4)*dt;
    fOut[3] += -0.125*((51.96152422706631*fT[2]-51.96152422706631*f[2]-45.0*fT[0]-45.0*f[0])*gT[3]+((-51.96152422706631*fB[2])+51.96152422706631*f[2]-45.0*fB[0]-45.0*f[0])*gB[3]+(51.96152422706631*fT[2]-51.96152422706631*fB[2]-45.0*fT[0]-45.0*fB[0]-90.0*f[0])*g[3]+(54.0*g[1]-54.0*gT[1])*fT[2]+(54.0*g[1]-54.0*gB[1])*fB[2]+(54.0*gT[1]+54.0*gB[1]-108.0*g[1])*f[2]+(46.76537180435967*fT[0]+46.76537180435967*f[0])*gT[1]+((-46.76537180435967*fB[0])-46.76537180435967*f[0])*gB[1]+(46.76537180435967*fB[0]-46.76537180435967*fT[0])*g[1])*gkyl_ipow(dv[1],-4)*dt;
  }

  fOut[0] += 0.0;
  fOut[1] += 0.0;
  fOut[2] += 0.0;
  fOut[3] += -1.0*((10.39230484541326*f[2]+6.0*f[0])*gT[3]+(6.0*f[0]-10.39230484541326*f[2])*gB[3]-12.0*f[0]*g[3]+((-9.0*gT[1])-9.0*gB[1]+18.0*g[1])*f[2]-5.196152422706631*f[0]*gT[1]+5.196152422706631*f[0]*gB[1])*gkyl_ipow(dv[1],-4)*dt;
}