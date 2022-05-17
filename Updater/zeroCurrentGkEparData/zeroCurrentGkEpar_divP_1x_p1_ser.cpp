#include <zeroCurrentGkEpar_mod_decl.h>

void zeroCurrentGkEpar_divP_1x_p1_ser(const double dz, const double charge, const double *delparFac, const double *dlnbmagdz, const double *m2par, const double *m2perp, double *divP) 
{ 
  // dz: cell length along the field line.
  // charge: electric charge.
  // delparFac: coefficient multiplying parallel gradient (cmag/(J*B)).
  // dlnbmagdz: d(ln bmag)/dz.
  // m2par: vpar^2 moment of the distribution function.
  // m2perp: mu*B/m moment of the distribution function.

  double divPpar[2];
  divPpar[0] = (2.449489742783178*delparFac[0]*m2par[1])/dz+(2.449489742783178*m2par[0]*delparFac[1])/dz; 

  double divPperp[2];
  divPperp[0] = 0.5*delparFac[0]*dlnbmagdz[1]*m2perp[1]+0.5*dlnbmagdz[0]*delparFac[1]*m2perp[1]+0.5*m2perp[0]*delparFac[1]*dlnbmagdz[1]+0.5*delparFac[0]*dlnbmagdz[0]*m2perp[0]; 
  divPperp[1] = 0.9*delparFac[1]*dlnbmagdz[1]*m2perp[1]+0.5*delparFac[0]*dlnbmagdz[0]*m2perp[1]+0.5*delparFac[0]*m2perp[0]*dlnbmagdz[1]+0.5*dlnbmagdz[0]*m2perp[0]*delparFac[1]; 

  divP[0] += divPperp[0]*charge+divPpar[0]*charge; 
  divP[1] += divPperp[1]*charge; 

}

