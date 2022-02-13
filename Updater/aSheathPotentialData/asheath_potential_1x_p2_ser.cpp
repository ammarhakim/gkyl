#include <asheath_potential_mod_decl.h>

void asheath_potential_gauss_phi_sheath_lower_1x_p2_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS) 
{ 
  // sheathDirDx: cell length in direction of the sheath.
  // q_e:         electron change.
  // m_e:         electron mass.
  // T_e:         electron temperature.
  // jacInv:      reciprocal of the geometry Jacobian (1/J).
  // GammaJac_i:  ion particle flux (times the Jacobian) through sheath entrance.
  // m0JacIon:    ion density (times the geometry Jacobian).
  // m0IonS:      ion density at the sheath entrance.
  // phiS:        electrostatic sheath potential.

  double GammaJacIonB[1];
  GammaJacIonB[0] = 0.7905694150420947*GammaJac_i[2]*sheathDirDx+0.6123724356957944*GammaJac_i[1]*sheathDirDx+0.3535533905932737*GammaJac_i[0]*sheathDirDx; 

  double m0JacIonB[1];
  m0JacIonB[0] = 1.58113883008419*m0JacIon[2]-1.224744871391589*m0JacIon[1]+0.7071067811865475*m0JacIon[0]; 

  double phiS_qp[1];
  if ((std::isfinite(GammaJacIonB[0]/2)) && (GammaJacIonB[0]/2>0.) && (m0JacIonB[0]/2>0.)) {
    phiS_qp[0] = (T_e*log((2.506628274631001*GammaJacIonB[0])/(m0JacIonB[0]*sqrt(T_e/m_e))))/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }


  phiS[0] = 1.414213562373095*phiS_qp[0]; 

  m0IonS[0] = 0.031943828249997*((53.75872022286246*jacInv[2]-34.2928563989645*jacInv[1]+49.49747468305833*jacInv[0])*m0JacIon[2]+(49.49747468305833*m0JacIon[0]-34.2928563989645*m0JacIon[1])*jacInv[2]+(66.40783086353598*jacInv[1]-38.34057902536163*jacInv[0])*m0JacIon[1]+m0JacIon[0]*(22.13594362117866*jacInv[0]-38.34057902536163*jacInv[1])); 
}

void asheath_potential_gauss_phi_sheath_upper_1x_p2_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS) 
{ 
  // sheathDirDx: cell length in direction of the sheath.
  // q_e:         electron change.
  // m_e:         electron mass.
  // T_e:         electron temperature.
  // jacInv:      reciprocal of the geometry Jacobian (1/J).
  // GammaJac_i:  ion particle flux (times the Jacobian) through sheath entrance.
  // m0JacIon:    ion density (times the geometry Jacobian).
  // m0IonS:      ion density at the sheath entrance.
  // phiS:        electrostatic sheath potential.

  double GammaJacIonB[1];
  GammaJacIonB[0] = 0.7905694150420947*GammaJac_i[2]*sheathDirDx-0.6123724356957944*GammaJac_i[1]*sheathDirDx+0.3535533905932737*GammaJac_i[0]*sheathDirDx; 

  double m0JacIonB[1];
  m0JacIonB[0] = 1.58113883008419*m0JacIon[2]+1.224744871391589*m0JacIon[1]+0.7071067811865475*m0JacIon[0]; 

  double phiS_qp[1];
  if ((std::isfinite(GammaJacIonB[0]/2)) && (GammaJacIonB[0]/2>0.) && (m0JacIonB[0]/2>0.)) {
    phiS_qp[0] = (T_e*log((2.506628274631001*GammaJacIonB[0])/(m0JacIonB[0]*sqrt(T_e/m_e))))/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }


  phiS[0] = 1.414213562373095*phiS_qp[0]; 

  m0IonS[0] = 0.031943828249997*((53.75872022286246*jacInv[2]+34.2928563989645*jacInv[1]+49.49747468305833*jacInv[0])*m0JacIon[2]+(34.2928563989645*m0JacIon[1]+49.49747468305833*m0JacIon[0])*jacInv[2]+(66.40783086353598*jacInv[1]+38.34057902536163*jacInv[0])*m0JacIon[1]+m0JacIon[0]*(38.34057902536163*jacInv[1]+22.13594362117866*jacInv[0])); 
}

void asheath_potential_gauss_phi_1x_p2_ser(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi) 
{ 
  // q_e:      electron change.
  // T_e:      electron temperature.
  // jacInv:   reciprocal of the geometry Jacobian (1/J).
  // m0JacIon: ion density.
  // m0IonS:   ion density at the sheath entrance.
  // phiS:     electrostatic sheath potential.
  // phi:      electrostatic potential in domain volume.

  double m0Ion[3];
  m0Ion[0] = 0.7071067811865475*jacInv[2]*m0JacIon[2]+0.7071067811865475*jacInv[1]*m0JacIon[1]+0.7071067811865475*jacInv[0]*m0JacIon[0]; 
  m0Ion[1] = 0.6324555320336759*jacInv[1]*m0JacIon[2]+0.6324555320336759*m0JacIon[1]*jacInv[2]+0.7071067811865475*jacInv[0]*m0JacIon[1]+0.7071067811865475*m0JacIon[0]*jacInv[1]; 
  m0Ion[2] = 0.4517539514526256*jacInv[2]*m0JacIon[2]+0.7071067811865475*jacInv[0]*m0JacIon[2]+0.7071067811865475*m0JacIon[0]*jacInv[2]+0.6324555320336759*jacInv[1]*m0JacIon[1]; 

  double phi_qp[3];
  phi_qp[0] = (-(1.0*log(m0Ion[0]/(m0IonS[0]-1.118033988749895*m0IonS[2])-(2.23606797749979*m0Ion[2])/(2.0*m0IonS[0]-2.23606797749979*m0IonS[2]))*T_e)/q_e)-0.7905694150420947*phiS[2]+0.7071067811865475*phiS[0]; 
  phi_qp[1] = (-(1.0*log((1.788854381999832*m0Ion[2])/(1.788854381999832*m0IonS[2]-2.683281572999748*m0IonS[1]+2.0*m0IonS[0])-(1.341640786499874*m0Ion[1])/(0.8944271909999162*m0IonS[2]-1.341640786499874*m0IonS[1]+m0IonS[0])+m0Ion[0]/(0.8944271909999162*m0IonS[2]-1.341640786499874*m0IonS[1]+m0IonS[0]))*T_e)/q_e)+0.6324555320336759*phiS[2]-0.9486832980505137*phiS[1]+0.7071067811865475*phiS[0]; 
  phi_qp[2] = (-(1.0*log((1.788854381999832*m0Ion[2])/(1.788854381999832*m0IonS[2]+2.683281572999748*m0IonS[1]+2.0*m0IonS[0])+(1.341640786499874*m0Ion[1])/(0.8944271909999162*m0IonS[2]+1.341640786499874*m0IonS[1]+m0IonS[0])+m0Ion[0]/(0.8944271909999162*m0IonS[2]+1.341640786499874*m0IonS[1]+m0IonS[0]))*T_e)/q_e)+0.6324555320336759*phiS[2]+0.9486832980505137*phiS[1]+0.7071067811865475*phiS[0]; 

  phi[0] = 0.392837100659193*(phi_qp[2]+phi_qp[1])+0.6285393610547089*phi_qp[0]; 
  phi[1] = 0.5270462766947299*phi_qp[2]-0.5270462766947299*phi_qp[1]; 
  phi[2] = 0.3513641844631533*(phi_qp[2]+phi_qp[1])-0.7027283689263064*phi_qp[0]; 
}

