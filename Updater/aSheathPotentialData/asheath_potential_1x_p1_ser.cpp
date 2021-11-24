#include <asheath_potential_mod_decl.h>

void asheath_potential_gauss_phi_sheath_lower_1x_p1_ser(const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS) 
{ 
  // q_e:        electron change.
  // m_e:        electron mass.
  // T_e:        electron temperature.
  // jacInv:     reciprocal of the geometry Jacobian (1/J).
  // GammaJac_i: ion particle flux (times the Jacobian) through sheath entrance.
  // m0JacIon:   ion density (times the geometry Jacobian).
  // m0IonS:     ion density at the sheath entrance.
  // phiS:       electrostatic sheath potential.

  double GammaJacIonB[1];
  GammaJacIonB[0] = 1.224744871391589*GammaJac_i[1]+0.7071067811865475*GammaJac_i[0]; 

  double m0JacIonB[1];
  m0JacIonB[0] = 0.7071067811865475*m0JacIon[0]-1.224744871391589*m0JacIon[1]; 

  double phiS_qp[1];
  if ((std::isfinite(GammaJacIonB[0]/2)) && (GammaJacIonB[0]/2>0.) && (m0JacIonB[0]/2>0.)) {
    phiS_qp[0] = (T_e*log((2.506628274631001*GammaJacIonB[0])/(m0JacIonB[0]*sqrt(T_e/m_e))))/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }


  phiS[0] = 1.414213562373095*phiS_qp[0]; 

  m0IonS[0] = 0.5*((1.414213562373095*jacInv[1]-2.449489742783178*jacInv[0])*m0JacIon[1]+m0JacIon[0]*(1.414213562373095*jacInv[0]-2.449489742783178*jacInv[1])); 
}

void asheath_potential_gauss_phi_sheath_upper_1x_p1_ser(const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS) 
{ 
  // q_e:        electron change.
  // m_e:        electron mass.
  // T_e:        electron temperature.
  // jacInv:     reciprocal of the geometry Jacobian (1/J).
  // GammaJac_i: ion particle flux (times the Jacobian) through sheath entrance.
  // m0JacIon:   ion density (times the geometry Jacobian).
  // m0IonS:     ion density at the sheath entrance.
  // phiS:       electrostatic sheath potential.

  double GammaJacIonB[1];
  GammaJacIonB[0] = 0.7071067811865475*GammaJac_i[0]-1.224744871391589*GammaJac_i[1]; 

  double m0JacIonB[1];
  m0JacIonB[0] = 1.224744871391589*m0JacIon[1]+0.7071067811865475*m0JacIon[0]; 

  double phiS_qp[1];
  if ((std::isfinite(GammaJacIonB[0]/2)) && (GammaJacIonB[0]/2>0.) && (m0JacIonB[0]/2>0.)) {
    phiS_qp[0] = (T_e*log((2.506628274631001*GammaJacIonB[0])/(m0JacIonB[0]*sqrt(T_e/m_e))))/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }


  phiS[0] = 1.414213562373095*phiS_qp[0]; 

  m0IonS[0] = 0.5*((1.414213562373095*jacInv[1]+2.449489742783178*jacInv[0])*m0JacIon[1]+m0JacIon[0]*(2.449489742783178*jacInv[1]+1.414213562373095*jacInv[0])); 
}

void asheath_potential_gauss_phi_1x_p1_ser(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi) 
{ 
  // q_e:      electron change.
  // T_e:      electron temperature.
  // jacInv:   reciprocal of the geometry Jacobian (1/J).
  // m0JacIon: ion density.
  // m0IonS:   ion density at the sheath entrance.
  // phiS:     electrostatic sheath potential.
  // phi:      electrostatic potential in domain volume.

  double m0Ion[2];
  m0Ion[0] = 0.7071067811865475*jacInv[1]*m0JacIon[1]+0.7071067811865475*jacInv[0]*m0JacIon[0]; 
  m0Ion[1] = 0.7071067811865475*jacInv[0]*m0JacIon[1]+0.7071067811865475*m0JacIon[0]*jacInv[1]; 

  double phi_qp[2];
  phi_qp[0] = (-(1.0*log(m0Ion[0]/(m0IonS[0]-0.9999999999999999*m0IonS[1])-(0.9999999999999999*m0Ion[1])/(m0IonS[0]-0.9999999999999999*m0IonS[1]))*T_e)/q_e)-0.7071067811865474*phiS[1]+0.7071067811865475*phiS[0]; 
  phi_qp[1] = (-(1.0*log((0.9999999999999999*m0Ion[1])/(0.9999999999999999*m0IonS[1]+m0IonS[0])+m0Ion[0]/(0.9999999999999999*m0IonS[1]+m0IonS[0]))*T_e)/q_e)+0.7071067811865474*phiS[1]+0.7071067811865475*phiS[0]; 

  phi[0] = 0.7071067811865475*(phi_qp[1]+phi_qp[0]); 
  phi[1] = 0.7071067811865474*phi_qp[1]-0.7071067811865474*phi_qp[0]; 
}

