#include <zeroCurrentGkEpar_divP_mod_decl.h>

void zeroCurrentGkEpar_divP_1x_p2_ser(const double dz, const double charge, const double *delparFac, const double *dlnbmagdz, const double *m2par, const double *m2perp, double *divP) 
{ 
  // dz: cell length along the field line.
  // charge: electric charge.
  // delparFac: coefficient multiplying parallel gradient (cmag/(J*B)).
  // dlnbmagdz: d(ln bmag)/dz.
  // m2par: vpar^2 moment of the distribution function.
  // m2perp: mu*B/m moment of the distribution function.

  double divPpar[3];
  divPpar[0] = (2.190890230020665*delparFac[1]*m2par[2])/dz+(2.190890230020665*m2par[1]*delparFac[2])/dz+(2.449489742783178*delparFac[0]*m2par[1])/dz+(2.449489742783178*m2par[0]*delparFac[1])/dz;
  divPpar[1] = (3.499271061118827*delparFac[2]*m2par[2])/dz+(5.477225575051662*delparFac[0]*m2par[2])/dz+(5.477225575051662*m2par[0]*delparFac[2])/dz+(4.898979485566357*delparFac[1]*m2par[1])/dz;

  double divPperp[3];
  divPperp[0] = 0.31943828249997*delparFac[2]*dlnbmagdz[2]*m2perp[2]+0.5*delparFac[0]*dlnbmagdz[2]*m2perp[2]+0.5*dlnbmagdz[0]*delparFac[2]*m2perp[2]+0.4472135954999579*delparFac[1]*dlnbmagdz[1]*m2perp[2]+0.5*m2perp[0]*delparFac[2]*dlnbmagdz[2]+0.4472135954999579*delparFac[1]*m2perp[1]*dlnbmagdz[2]+0.4472135954999579*dlnbmagdz[1]*m2perp[1]*delparFac[2]+0.5*delparFac[0]*dlnbmagdz[1]*m2perp[1]+0.5*dlnbmagdz[0]*delparFac[1]*m2perp[1]+0.5*m2perp[0]*delparFac[1]*dlnbmagdz[1]+0.5*delparFac[0]*dlnbmagdz[0]*m2perp[0]; 
  divPperp[1] = 0.7857142857142857*delparFac[1]*dlnbmagdz[2]*m2perp[2]+0.7857142857142857*dlnbmagdz[1]*delparFac[2]*m2perp[2]+0.4472135954999579*delparFac[0]*dlnbmagdz[1]*m2perp[2]+0.4472135954999579*dlnbmagdz[0]*delparFac[1]*m2perp[2]+0.7857142857142857*m2perp[1]*delparFac[2]*dlnbmagdz[2]+0.4472135954999579*delparFac[0]*m2perp[1]*dlnbmagdz[2]+0.4472135954999579*m2perp[0]*delparFac[1]*dlnbmagdz[2]+0.4472135954999579*dlnbmagdz[0]*m2perp[1]*delparFac[2]+0.4472135954999579*m2perp[0]*dlnbmagdz[1]*delparFac[2]+0.9*delparFac[1]*dlnbmagdz[1]*m2perp[1]+0.5*delparFac[0]*dlnbmagdz[0]*m2perp[1]+0.5*delparFac[0]*m2perp[0]*dlnbmagdz[1]+0.5*dlnbmagdz[0]*m2perp[0]*delparFac[1]; 
  divPperp[2] = 1.071428571428571*delparFac[2]*dlnbmagdz[2]*m2perp[2]+0.31943828249997*delparFac[0]*dlnbmagdz[2]*m2perp[2]+0.31943828249997*dlnbmagdz[0]*delparFac[2]*m2perp[2]+0.7857142857142857*delparFac[1]*dlnbmagdz[1]*m2perp[2]+0.5*delparFac[0]*dlnbmagdz[0]*m2perp[2]+0.31943828249997*m2perp[0]*delparFac[2]*dlnbmagdz[2]+0.7857142857142857*delparFac[1]*m2perp[1]*dlnbmagdz[2]+0.5*delparFac[0]*m2perp[0]*dlnbmagdz[2]+0.7857142857142857*dlnbmagdz[1]*m2perp[1]*delparFac[2]+0.5*dlnbmagdz[0]*m2perp[0]*delparFac[2]+0.4472135954999579*delparFac[0]*dlnbmagdz[1]*m2perp[1]+0.4472135954999579*dlnbmagdz[0]*delparFac[1]*m2perp[1]+0.4472135954999579*m2perp[0]*delparFac[1]*dlnbmagdz[1]; 

  divP[0] += divPperp[0]*charge+divPpar[0]*charge; 
  divP[1] += divPperp[1]*charge+divPpar[1]*charge; 
  divP[2] += divPperp[2]*charge; 

}

