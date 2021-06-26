#include <math.h>
#include <vlasov_fpo.h>

void vlasov_fpo_moms_3x_ser_p1(const double* dv, const double* vc, const double* f, double* out) {
  out[0] += 0.3535533905932737*dv[0]*f[0]*dv[1]*dv[2];
  out[1] += 0.1020620726159657*dv[1]*f[1]*dv[2]*pow(dv[0],2)+0.3535533905932737*dv[0]*f[0]*vc[0]*dv[1]*dv[2];
  out[2] += 0.1020620726159657*dv[0]*dv[2]*f[2]*pow(dv[1],2)+0.3535533905932737*dv[0]*f[0]*dv[1]*vc[1]*dv[2];
  out[3] += 0.1020620726159657*dv[0]*dv[1]*f[3]*pow(dv[2],2)+0.3535533905932737*dv[0]*f[0]*dv[1]*dv[2]*vc[2];
  out[4] += 0.3535533905932737*dv[0]*f[0]*dv[1]*dv[2]*pow(vc[2],2)+0.02946278254943947*dv[0]*f[0]*dv[1]*pow(dv[2],3)+0.2041241452319315*dv[0]*dv[1]*vc[2]*f[3]*pow(dv[2],2)+0.3535533905932737*dv[0]*f[0]*dv[1]*dv[2]*pow(vc[1],2)+0.02946278254943947*dv[0]*f[0]*dv[2]*pow(dv[1],3)+0.2041241452319315*dv[0]*vc[1]*dv[2]*f[2]*pow(dv[1],2)+0.3535533905932737*dv[0]*f[0]*dv[1]*dv[2]*pow(vc[0],2)+0.02946278254943947*f[0]*dv[1]*dv[2]*pow(dv[0],3)+0.2041241452319315*vc[0]*dv[1]*f[1]*dv[2]*pow(dv[0],2);
}

double vlasov_fpo_momM0mod_3x_ser_p1(const double* dv, const double* vc, const double* f,
                                     int isXloEdge, int isXupEdge,
                                     int isYloEdge, int isYupEdge,
                                     int isZloEdge, int isZupEdge) {
  double out = 3 * 0.3535533905932737*dv[0]*f[0]*dv[1]*dv[2];
  if (isXloEdge) {
    out += (-0.6123724356957944*vc[0]*dv[1]*f[1]*dv[2])+0.3061862178478971*dv[0]*dv[1]*f[1]*dv[2]+0.3535533905932737*f[0]*vc[0]*dv[1]*dv[2]-0.1767766952966368*dv[0]*f[0]*dv[1]*dv[2];
  } else if (isXupEdge) {
    out += (-0.6123724356957944*vc[0]*dv[1]*f[1]*dv[2])-0.3061862178478971*dv[0]*dv[1]*f[1]*dv[2]-0.3535533905932737*f[0]*vc[0]*dv[1]*dv[2]-0.1767766952966368*dv[0]*f[0]*dv[1]*dv[2];
  } else if (isYloEdge) {
    out += (-0.6123724356957944*dv[0]*vc[1]*dv[2]*f[2])+0.3061862178478971*dv[0]*dv[1]*dv[2]*f[2]+0.3535533905932737*dv[0]*f[0]*vc[1]*dv[2]-0.1767766952966368*dv[0]*f[0]*dv[1]*dv[2];
  } else if (isYupEdge) {
    out += (-0.6123724356957944*dv[0]*vc[1]*dv[2]*f[2])-0.3061862178478971*dv[0]*dv[1]*dv[2]*f[2]-0.3535533905932737*dv[0]*f[0]*vc[1]*dv[2]-0.1767766952966368*dv[0]*f[0]*dv[1]*dv[2];
  } else if (isZloEdge) {
    out += (-0.6123724356957944*dv[0]*dv[1]*vc[2]*f[3])+0.3061862178478971*dv[0]*dv[1]*dv[2]*f[3]+0.3535533905932737*dv[0]*f[0]*dv[1]*vc[2]-0.1767766952966368*dv[0]*f[0]*dv[1]*dv[2];
  } else {
    out += (-0.6123724356957944*dv[0]*dv[1]*vc[2]*f[3])-0.3061862178478971*dv[0]*dv[1]*dv[2]*f[3]-0.3535533905932737*dv[0]*f[0]*dv[1]*vc[2]-0.1767766952966368*dv[0]*f[0]*dv[1]*dv[2];
  }
  return out;
}

double vlasov_fpo_momM2_3x_ser_p1(const double* dv, const double* vc, const double* f) {
  double out = 0.3535533905932737*dv[0]*f[0]*dv[1]*dv[2]*pow(vc[2],2)+0.02946278254943947*dv[0]*f[0]*dv[1]*pow(dv[2],3)+0.2041241452319315*dv[0]*dv[1]*vc[2]*f[3]*pow(dv[2],2)+0.3535533905932737*dv[0]*f[0]*dv[1]*dv[2]*pow(vc[1],2)+0.02946278254943947*dv[0]*f[0]*dv[2]*pow(dv[1],3)+0.2041241452319315*dv[0]*vc[1]*dv[2]*f[2]*pow(dv[1],2)+0.3535533905932737*dv[0]*f[0]*dv[1]*dv[2]*pow(vc[0],2)+0.02946278254943947*f[0]*dv[1]*dv[2]*pow(dv[0],3)+0.2041241452319315*vc[0]*dv[1]*f[1]*dv[2]*pow(dv[0],2);
  return out;
}
