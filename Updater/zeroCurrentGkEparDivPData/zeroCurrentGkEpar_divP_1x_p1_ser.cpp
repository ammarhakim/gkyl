#include <zeroCurrentGkEpar_mod_decl.h>

void zeroCurrentGkEpar_divP_1x_p1_ser(const double dz, const double charge, const double *jacobGeoInv, const double *delparFac, const double *dlnbmagdz, const double *m2par, const double *m2perp, const double *divP) 
{ 
  // dz: cell length along the field line.
  // charge: electric charge.
  // jacobGeoInv: reciprocal of the conf-space jacobian.
  // delparFac: coefficient multiplying parallel gradient (cmag/(J*B)).
  // dlnbmagdz: d(ln bmag)/dz.
  // m2par: vpar^2 moment of the distribution function.
  // m2perp: mu*B/m moment of the distribution function.

  double divPpar[2];
  divPpar[0] = 0.6123724356957944*delparFac[0]*m2par[1]*dz+0.6123724356957944*m2par[0]*delparFac[1]*dz; 

  double divPperp[2];
  divPperp[0] = 0.5*delparFac[0]*dlnbmagdz[1]*m2perp[1]+0.5*dlnbmagdz[0]*delparFac[1]*m2perp[1]+0.5*m2perp[0]*delparFac[1]*dlnbmagdz[1]+0.5*delparFac[0]*dlnbmagdz[0]*m2perp[0]; 
  divPperp[1] = 0.9*delparFac[1]*dlnbmagdz[1]*m2perp[1]+0.5*delparFac[0]*dlnbmagdz[0]*m2perp[1]+0.5*delparFac[0]*m2perp[0]*dlnbmagdz[1]+0.5*dlnbmagdz[0]*m2perp[0]*delparFac[1]; 

  divP[0] += 0.7071067811865475*divPperp[1]*jacobGeoInv[1]*charge+0.7071067811865475*divPperp[0]*jacobGeoInv[0]*charge+0.7071067811865475*divPpar[0]*jacobGeoInv[0]*charge; 
  divP[1] += 0.7071067811865475*divPperp[0]*jacobGeoInv[1]*charge+0.7071067811865475*divPpar[0]*jacobGeoInv[1]*charge+0.7071067811865475*jacobGeoInv[0]*divPperp[1]*charge; 

}

