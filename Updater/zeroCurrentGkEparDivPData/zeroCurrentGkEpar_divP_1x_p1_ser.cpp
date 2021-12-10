#include <zeroCurrentGkEpar_mod_decl.h>

void zeroCurrentGkEpar_divP_1x_p1_ser(const double dz, const double charge, const double *cmag, const double *jacobTotInv, const double *jacobGeoInv, const double *bmagInv, const double *dbmagdz, const double *m2par, const double *m2perp, const double *divP) 
{ 
  // dz: cell length along the field line.
  // charge: electric charge.
  // cmag: coefficient multiplying parallel gradient.
  // jacobTotInv: reciprocal of the conf-space jacobian times the guiding center coordinate Jacobian.
  // jacobGeoInv: reciprocal of the conf-space jacobian.
  // bmagInv: reciprocal of magnetic field magnitude.
  // dbmagdz: d(bmag)/dz.
  // m2par: vpar^2 moment of the distribution function.
  // m2perp: mu*B/m moment of the distribution function.

  double divPpar[2];
  divPpar[0] = 0.7794228634059945*cmag[1]*jacobTotInv[1]*m2par[1]*dz+0.4330127018922193*cmag[0]*jacobTotInv[0]*m2par[1]*dz+0.4330127018922193*cmag[0]*m2par[0]*jacobTotInv[1]*dz+0.4330127018922193*jacobTotInv[0]*m2par[0]*cmag[1]*dz; 

  double divPperp[2];
  divPperp[0] = 0.45*bmagInv[0]*cmag[1]*dbmagdz[1]*jacobTotInv[1]*m2perp[1]+0.45*cmag[0]*bmagInv[1]*dbmagdz[1]*jacobTotInv[1]*m2perp[1]+0.45*dbmagdz[0]*bmagInv[1]*cmag[1]*jacobTotInv[1]*m2perp[1]+0.25*bmagInv[0]*cmag[0]*dbmagdz[0]*jacobTotInv[1]*m2perp[1]+0.45*jacobTotInv[0]*bmagInv[1]*cmag[1]*dbmagdz[1]*m2perp[1]+0.25*bmagInv[0]*cmag[0]*jacobTotInv[0]*dbmagdz[1]*m2perp[1]+0.25*bmagInv[0]*dbmagdz[0]*jacobTotInv[0]*cmag[1]*m2perp[1]+0.25*cmag[0]*dbmagdz[0]*jacobTotInv[0]*bmagInv[1]*m2perp[1]+0.45*m2perp[0]*bmagInv[1]*cmag[1]*dbmagdz[1]*jacobTotInv[1]+0.25*bmagInv[0]*cmag[0]*m2perp[0]*dbmagdz[1]*jacobTotInv[1]+0.25*bmagInv[0]*dbmagdz[0]*m2perp[0]*cmag[1]*jacobTotInv[1]+0.25*cmag[0]*dbmagdz[0]*m2perp[0]*bmagInv[1]*jacobTotInv[1]+0.25*bmagInv[0]*jacobTotInv[0]*m2perp[0]*cmag[1]*dbmagdz[1]+0.25*cmag[0]*jacobTotInv[0]*m2perp[0]*bmagInv[1]*dbmagdz[1]+0.25*dbmagdz[0]*jacobTotInv[0]*m2perp[0]*bmagInv[1]*cmag[1]+0.25*bmagInv[0]*cmag[0]*dbmagdz[0]*jacobTotInv[0]*m2perp[0]; 
  divPperp[1] = 0.9642857142857143*bmagInv[1]*cmag[1]*dbmagdz[1]*jacobTotInv[1]*m2perp[1]+0.45*bmagInv[0]*cmag[0]*dbmagdz[1]*jacobTotInv[1]*m2perp[1]+0.45*bmagInv[0]*dbmagdz[0]*cmag[1]*jacobTotInv[1]*m2perp[1]+0.45*cmag[0]*dbmagdz[0]*bmagInv[1]*jacobTotInv[1]*m2perp[1]+0.45*bmagInv[0]*jacobTotInv[0]*cmag[1]*dbmagdz[1]*m2perp[1]+0.45*cmag[0]*jacobTotInv[0]*bmagInv[1]*dbmagdz[1]*m2perp[1]+0.45*dbmagdz[0]*jacobTotInv[0]*bmagInv[1]*cmag[1]*m2perp[1]+0.25*bmagInv[0]*cmag[0]*dbmagdz[0]*jacobTotInv[0]*m2perp[1]+0.45*bmagInv[0]*m2perp[0]*cmag[1]*dbmagdz[1]*jacobTotInv[1]+0.45*cmag[0]*m2perp[0]*bmagInv[1]*dbmagdz[1]*jacobTotInv[1]+0.45*dbmagdz[0]*m2perp[0]*bmagInv[1]*cmag[1]*jacobTotInv[1]+0.25*bmagInv[0]*cmag[0]*dbmagdz[0]*m2perp[0]*jacobTotInv[1]+0.45*jacobTotInv[0]*m2perp[0]*bmagInv[1]*cmag[1]*dbmagdz[1]+0.25*bmagInv[0]*cmag[0]*jacobTotInv[0]*m2perp[0]*dbmagdz[1]+0.25*bmagInv[0]*dbmagdz[0]*jacobTotInv[0]*m2perp[0]*cmag[1]+0.25*cmag[0]*dbmagdz[0]*jacobTotInv[0]*m2perp[0]*bmagInv[1]; 

  divP[0] += 0.7071067811865475*divPperp[1]*jacobGeoInv[1]*charge+0.7071067811865475*divPperp[0]*jacobGeoInv[0]*charge+0.7071067811865475*divPpar[0]*jacobGeoInv[0]*charge; 
  divP[1] += 0.7071067811865475*divPperp[0]*jacobGeoInv[1]*charge+0.7071067811865475*divPpar[0]*jacobGeoInv[1]*charge+0.7071067811865475*jacobGeoInv[0]*divPperp[1]*charge; 

}

