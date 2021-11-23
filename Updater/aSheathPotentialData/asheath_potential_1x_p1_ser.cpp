#include <asheath_potential_mod_decl.h>

void asheath_potential_gauss_phi_sheath_lower_1x_p1_ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS) 
{ 
  // q_e:     electron change.
  // m_e:     electron mass.
  // T_e:     electron temperature.
  // Gamma_i: ion particle flux through sheath entrance.
  // m0Ion:   ion density.
  // m0IonS:  ion density at the sheath entrance.
  // phiS:    electrostatic sheath potential.

  double GammaIonB[1];
  GammaIonB[0] = 1.224744871391589*Gamma_i[1]+0.7071067811865475*Gamma_i[0]; 

  double m0IonB[1];
  m0IonB[0] = 0.7071067811865475*m0Ion[0]-1.224744871391589*m0Ion[1]; 


  phiS[0] = (1.414213562373095*T_e*log((2.506628274631001*GammaIonB[0])/(m0IonB[0]*sqrt(T_e/m_e))))/q_e; 

  m0IonS[0] = m0Ion[0]-1.732050807568877*m0Ion[1]; 
}

void asheath_potential_gauss_phi_sheath_upper_1x_p1_ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS) 
{ 
  // q_e:     electron change.
  // m_e:     electron mass.
  // T_e:     electron temperature.
  // Gamma_i: ion particle flux through sheath entrance.
  // m0Ion:   ion density.
  // m0IonS:  ion density at the sheath entrance.
  // phiS:    electrostatic sheath potential.

  double GammaIonB[1];
  GammaIonB[0] = 0.7071067811865475*Gamma_i[0]-1.224744871391589*Gamma_i[1]; 

  double m0IonB[1];
  m0IonB[0] = 1.224744871391589*m0Ion[1]+0.7071067811865475*m0Ion[0]; 


  phiS[0] = (1.414213562373095*T_e*log((2.506628274631001*GammaIonB[0])/(m0IonB[0]*sqrt(T_e/m_e))))/q_e; 

  m0IonS[0] = 1.732050807568877*m0Ion[1]+m0Ion[0]; 
}

void asheath_potential_gauss_phi_1x_p1_ser(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi) 
{ 
  // q_e:    electron change.
  // T_e:    electron temperature.
  // m0Ion:  ion density.
  // m0IonS: ion density at the sheath entrance.
  // phiS:   electrostatic sheath potential.
  // phi:    electrostatic potential in domain volume.

  double phi_qp[2];
  phi_qp[0] = (-(1.0*log(m0Ion[0]/(m0IonS[0]-0.9999999999999999*m0IonS[1])-(0.9999999999999999*m0Ion[1])/(m0IonS[0]-0.9999999999999999*m0IonS[1]))*T_e)/q_e)-0.7071067811865474*phiS[1]+0.7071067811865475*phiS[0]; 
  phi_qp[1] = (-(1.0*log((0.9999999999999999*m0Ion[1])/(0.9999999999999999*m0IonS[1]+m0IonS[0])+m0Ion[0]/(0.9999999999999999*m0IonS[1]+m0IonS[0]))*T_e)/q_e)+0.7071067811865474*phiS[1]+0.7071067811865475*phiS[0]; 

  phi[0] = 0.7071067811865475*(phi_qp[1]+phi_qp[0]); 
  phi[1] = 0.7071067811865474*phi_qp[1]-0.7071067811865474*phi_qp[0]; 
}

