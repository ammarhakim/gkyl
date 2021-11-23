#ifndef ASHEATHPOTENTIAL_MOD_DELC_H 
#define ASHEATHPOTENTIAL_MOD_DELC_H 

#include <cmath>

extern "C" { 

  void asheath_potential_gauss_phi_sheath_lower_1x_p1_Ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_1x_p1_Ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_1x_p1_Ser(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_2x_p1_Ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_2x_p1_Ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_2x_p1_Ser(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_3x_p1_Ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_3x_p1_Ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_3x_p1_Ser(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi); 


  void asheath_potential_gauss_phi_sheath_lower_1x_p2_Ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_1x_p2_Ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_1x_p2_Ser(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_2x_p2_Ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_2x_p2_Ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_2x_p2_Ser(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_3x_p2_Ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_3x_p2_Ser(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_3x_p2_Ser(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi); 


  void asheath_potential_gauss_phi_sheath_lower_1x_p1_Tensor(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_1x_p1_Tensor(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_1x_p1_Tensor(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_2x_p1_Tensor(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_2x_p1_Tensor(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_2x_p1_Tensor(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_3x_p1_Tensor(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_3x_p1_Tensor(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_3x_p1_Tensor(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi); 


  void asheath_potential_gauss_phi_sheath_lower_1x_p2_Tensor(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_1x_p2_Tensor(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_1x_p2_Tensor(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_2x_p2_Tensor(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_2x_p2_Tensor(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_2x_p2_Tensor(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi); 

  void asheath_potential_gauss_phi_sheath_lower_3x_p2_Tensor(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_sheath_upper_3x_p2_Tensor(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS); 
  void asheath_potential_gauss_phi_3x_p2_Tensor(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi); 



 
}

#endif 
