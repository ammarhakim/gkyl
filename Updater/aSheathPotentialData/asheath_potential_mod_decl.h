#ifndef ASHEATHPOTENTIAL_MOD_DELC_H
#define ASHEATHPOTENTIAL_MOD_DELC_H

#include <cmath>

extern "C" { 

  void asheath_potential_gauss_phi_sheath_lower_1x_p1_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_1x_p1_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_1x_p1_ser(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_2x_p1_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_2x_p1_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_2x_p1_ser(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_3x_p1_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_3x_p1_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_3x_p1_ser(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi); 


  void asheath_potential_gauss_phi_sheath_lower_1x_p2_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_1x_p2_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_1x_p2_ser(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_2x_p2_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_2x_p2_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_2x_p2_ser(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_3x_p2_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_3x_p2_ser(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_3x_p2_ser(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi); 


  void asheath_potential_gauss_phi_sheath_lower_1x_p1_tensor(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_1x_p1_tensor(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_1x_p1_tensor(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_2x_p1_tensor(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_2x_p1_tensor(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_2x_p1_tensor(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_3x_p1_tensor(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_3x_p1_tensor(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_3x_p1_tensor(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi); 


  void asheath_potential_gauss_phi_sheath_lower_1x_p2_tensor(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_1x_p2_tensor(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_1x_p2_tensor(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_2x_p2_tensor(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_2x_p2_tensor(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_2x_p2_tensor(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_3x_p2_tensor(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_3x_p2_tensor(const double sheathDirDx, const double q_e, const double m_e, const double T_e, const int dirSign, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_3x_p2_tensor(const double q_e, const double T_e, const double *jacInv, const double *m0JacIon, const double *m0IonS, const double *phiS, double *phi); 



 
}

#endif 
