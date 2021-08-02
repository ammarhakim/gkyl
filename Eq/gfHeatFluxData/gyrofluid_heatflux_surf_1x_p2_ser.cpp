#include <gyrofluid_heatflux_mod_decl.h>

void gyrofluid_heatflux_surf_1x_p2_ser_x(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *wL1, const double *wR1, const double *dxL1, const double *dxR1, const double cMaxIn, const double *rBmagL1, const double *rBmagR1, const double *rBmagSqL1, const double *rBmagSqR1, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, double *outL, double *outR) 
{ 
  // q_,m_:              species charge and mass.
  // kappaPar,kappaPerp: heat conductivity coefficients.
  // kperpSq:            k_perp^2.
  // wL,wR:              cell-center in left and right cells.
  // dxL,dxR:            cell length in left and right cells.
  // cMaxIn:             maximum sound speed (or some factor like it).
  // rBmag:              reciprocal of magnetic field magnitude (1/B).
  // rBmagSq:            rBmag^2.
  // sMomL,sMomR:        stepped moments (times Jacobian) in left and right cells.
  // phiL,phiR:          electrostatic potential in left and right cells.
  // primMomL,primMomR:  primitive moments (upar, Tpar, Tperp) in left and right cells.
  // outL/outR:          output increment in left and right cells.

  double rdx2L = 2.0/dxL1[0];
  double rdxSq4L = rdx2L*rdx2L;
  double rdx2R = 2.0/dxR1[0];
  double rdxSq4R = rdx2R*rdx2R;

  double GheatF3[3];

  GheatF3[0] = -0.00625*(((169.7056274847715*phiR1[2]-169.7056274847715*phiL1[2]-301.2474066278414*phiR1[1]-301.2474066278414*phiL1[1]+237.1708245126285*phiR1[0]-237.1708245126285*phiL1[0])*rBmagL1[2]+(131.4534138012399*rBmagL1[1]+75.89466384404115*rBmagL1[0])*phiR1[2]+((-131.4534138012399*rBmagL1[1])-75.89466384404115*rBmagL1[0])*phiL1[2]+((-233.3452377915607*phiR1[1])-233.3452377915607*phiL1[1]+183.7117307087383*phiR1[0]-183.7117307087383*phiL1[0])*rBmagL1[1]-134.7219358530748*rBmagL1[0]*phiR1[1]-134.7219358530748*rBmagL1[0]*phiL1[1]+(106.0660171779821*phiR1[0]-106.0660171779821*phiL1[0])*rBmagL1[0])*kappaPerp*kperpSq*q_+(((-339.411254969543*rBmagL1[2])-262.9068276024798*rBmagL1[1]-151.7893276880823*rBmagL1[0])*primMomR1[8]+(339.411254969543*rBmagL1[2]+262.9068276024798*rBmagL1[1]+151.7893276880823*rBmagL1[0])*primMomL1[8]+(602.494813255683*rBmagL1[2]+466.6904755831215*rBmagL1[1]+269.4438717061497*rBmagL1[0])*primMomR1[7]+(602.494813255683*rBmagL1[2]+466.6904755831215*rBmagL1[1]+269.4438717061497*rBmagL1[0])*primMomL1[7]+((-474.3416490252571*rBmagL1[2])-367.4234614174767*rBmagL1[1]-212.1320343559643*rBmagL1[0])*primMomR1[6]+(474.3416490252571*rBmagL1[2]+367.4234614174767*rBmagL1[1]+212.1320343559643*rBmagL1[0])*primMomL1[6])*kappaPerp+(((-169.7056274847715*rBmagL1[2])-131.4534138012399*rBmagL1[1]-75.89466384404115*rBmagL1[0])*primMomR1[5]+(169.7056274847715*rBmagL1[2]+131.4534138012399*rBmagL1[1]+75.89466384404115*rBmagL1[0])*primMomL1[5]+(301.2474066278414*rBmagL1[2]+233.3452377915607*rBmagL1[1]+134.7219358530748*rBmagL1[0])*primMomR1[4]+(301.2474066278414*rBmagL1[2]+233.3452377915607*rBmagL1[1]+134.7219358530748*rBmagL1[0])*primMomL1[4]+((-237.1708245126285*rBmagL1[2])-183.7117307087383*rBmagL1[1]-106.0660171779821*rBmagL1[0])*primMomR1[3]+(237.1708245126285*rBmagL1[2]+183.7117307087383*rBmagL1[1]+106.0660171779821*rBmagL1[0])*primMomL1[3])*kappaPar)*rdx2L; 

  double GheatF4[3];

  GheatF4[0] = -0.00625*(((169.7056274847715*phiR1[2]-169.7056274847715*phiL1[2]-301.2474066278414*phiR1[1]-301.2474066278414*phiL1[1]+237.1708245126285*phiR1[0]-237.1708245126285*phiL1[0])*rBmagSqL1[2]+(131.4534138012399*rBmagSqL1[1]+75.89466384404115*rBmagSqL1[0])*phiR1[2]+((-131.4534138012399*rBmagSqL1[1])-75.89466384404115*rBmagSqL1[0])*phiL1[2]+((-233.3452377915607*phiR1[1])-233.3452377915607*phiL1[1]+183.7117307087383*phiR1[0]-183.7117307087383*phiL1[0])*rBmagSqL1[1]-134.7219358530748*rBmagSqL1[0]*phiR1[1]-134.7219358530748*rBmagSqL1[0]*phiL1[1]+(106.0660171779821*phiR1[0]-106.0660171779821*phiL1[0])*rBmagSqL1[0])*kappaPerp*kperpSq*q_+(((-339.411254969543*rBmagSqL1[2])-262.9068276024798*rBmagSqL1[1]-151.7893276880823*rBmagSqL1[0])*primMomR1[8]+(339.411254969543*rBmagSqL1[2]+262.9068276024798*rBmagSqL1[1]+151.7893276880823*rBmagSqL1[0])*primMomL1[8]+(602.494813255683*rBmagSqL1[2]+466.6904755831215*rBmagSqL1[1]+269.4438717061497*rBmagSqL1[0])*primMomR1[7]+(602.494813255683*rBmagSqL1[2]+466.6904755831215*rBmagSqL1[1]+269.4438717061497*rBmagSqL1[0])*primMomL1[7]+((-474.3416490252571*rBmagSqL1[2])-367.4234614174767*rBmagSqL1[1]-212.1320343559643*rBmagSqL1[0])*primMomR1[6]+(474.3416490252571*rBmagSqL1[2]+367.4234614174767*rBmagSqL1[1]+212.1320343559643*rBmagSqL1[0])*primMomL1[6])*kappaPerp)*rdx2L; 

  double incr3[3];
  incr3[0] = -0.5*GheatF3[0]; 
  incr3[1] = 0.8660254037844386*GheatF3[0]; 
  incr3[2] = -1.118033988749895*GheatF3[0]; 

  double incrNonFlux3[3];
  incrNonFlux3[1] = -0.00390625*(kappaPerp*(((85.73214099741124*(phiR1[2]+phiL1[2])-123.3288287465668*phiR1[1]+123.3288287465668*phiL1[1]+87.63560920082662*(phiR1[0]+phiL1[0]))*rBmagL1[2]+(66.40783086353598*rBmagL1[1]+38.34057902536163*rBmagL1[0])*(phiR1[2]+phiL1[2])+((-95.53009996854395*phiR1[1])+95.53009996854395*phiL1[1]+67.8822509939086*(phiR1[0]+phiL1[0]))*rBmagL1[1]+rBmagL1[0]*((-55.15432893255071*phiR1[1])+55.15432893255071*phiL1[1]+39.19183588453087*(phiR1[0]+phiL1[0])))*kperpSq*q_+((-171.4642819948225*rBmagL1[2])-132.815661727072*rBmagL1[1]-76.68115805072327*rBmagL1[0])*(primMomR1[8]+primMomL1[8])+(246.6576574931337*rBmagL1[2]+191.0601999370879*rBmagL1[1]+110.3086578651014*rBmagL1[0])*primMomR1[7]+((-246.6576574931337*rBmagL1[2])-191.0601999370879*rBmagL1[1]-110.3086578651014*rBmagL1[0])*primMomL1[7]+((-175.2712184016533*rBmagL1[2])-135.7645019878172*rBmagL1[1]-78.38367176906175*rBmagL1[0])*(primMomR1[6]+primMomL1[6]))+(((-85.73214099741124*rBmagL1[2])-66.40783086353598*rBmagL1[1]-38.34057902536163*rBmagL1[0])*(primMomR1[5]+primMomL1[5])+(123.3288287465668*rBmagL1[2]+95.53009996854395*rBmagL1[1]+55.15432893255071*rBmagL1[0])*primMomR1[4]+((-123.3288287465668*rBmagL1[2])-95.53009996854395*rBmagL1[1]-55.15432893255071*rBmagL1[0])*primMomL1[4]+((-87.63560920082662*rBmagL1[2])-67.8822509939086*rBmagL1[1]-39.19183588453087*rBmagL1[0])*(primMomR1[3]+primMomL1[3]))*kappaPar); 
  incrNonFlux3[2] = 0.00390625*(kappaPerp*(((332.0391543176799*(phiR1[2]+phiL1[2])-477.6504998427197*phiR1[1]+477.6504998427197*phiL1[1]+339.411254969543*(phiR1[0]+phiL1[0]))*rBmagL1[2]+(257.1964229922337*rBmagL1[1]+148.492424049175*rBmagL1[0])*(phiR1[2]+phiL1[2])+((-369.9864862397004*phiR1[1])+369.9864862397004*phiL1[1]+262.9068276024798*(phiR1[0]+phiL1[0]))*rBmagL1[1]+rBmagL1[0]*((-213.6117974270148*phiR1[1])+213.6117974270148*phiL1[1]+151.7893276880823*(phiR1[0]+phiL1[0])))*kperpSq*q_+((-664.0783086353599*rBmagL1[2])-514.3928459844674*rBmagL1[1]-296.98484809835*rBmagL1[0])*(primMomR1[8]+primMomL1[8])+(955.3009996854395*rBmagL1[2]+739.9729724794009*rBmagL1[1]+427.2235948540296*rBmagL1[0])*primMomR1[7]+((-955.3009996854395*rBmagL1[2])-739.9729724794009*rBmagL1[1]-427.2235948540296*rBmagL1[0])*primMomL1[7]+((-678.8225099390861*rBmagL1[2])-525.8136552049598*rBmagL1[1]-303.5786553761646*rBmagL1[0])*(primMomR1[6]+primMomL1[6]))+(((-332.0391543176799*rBmagL1[2])-257.1964229922337*rBmagL1[1]-148.492424049175*rBmagL1[0])*(primMomR1[5]+primMomL1[5])+(477.6504998427197*rBmagL1[2]+369.9864862397004*rBmagL1[1]+213.6117974270148*rBmagL1[0])*primMomR1[4]+((-477.6504998427197*rBmagL1[2])-369.9864862397004*rBmagL1[1]-213.6117974270148*rBmagL1[0])*primMomL1[4]+((-339.411254969543*rBmagL1[2])-262.9068276024798*rBmagL1[1]-151.7893276880823*rBmagL1[0])*(primMomR1[3]+primMomL1[3]))*kappaPar); 

  double incr4[3];
  incr4[0] = -0.5*GheatF4[0]; 
  incr4[1] = 0.8660254037844386*GheatF4[0]; 
  incr4[2] = -1.118033988749895*GheatF4[0]; 

  double incrNonFlux4[3];
  incrNonFlux4[1] = -0.00390625*kappaPerp*(((85.73214099741124*(phiR1[2]+phiL1[2])-123.3288287465668*phiR1[1]+123.3288287465668*phiL1[1]+87.63560920082662*(phiR1[0]+phiL1[0]))*rBmagSqL1[2]+(66.40783086353598*rBmagSqL1[1]+38.34057902536163*rBmagSqL1[0])*(phiR1[2]+phiL1[2])+((-95.53009996854395*phiR1[1])+95.53009996854395*phiL1[1]+67.8822509939086*(phiR1[0]+phiL1[0]))*rBmagSqL1[1]+rBmagSqL1[0]*((-55.15432893255071*phiR1[1])+55.15432893255071*phiL1[1]+39.19183588453087*(phiR1[0]+phiL1[0])))*kperpSq*q_+((-171.4642819948225*rBmagSqL1[2])-132.815661727072*rBmagSqL1[1]-76.68115805072327*rBmagSqL1[0])*(primMomR1[8]+primMomL1[8])+(246.6576574931337*rBmagSqL1[2]+191.0601999370879*rBmagSqL1[1]+110.3086578651014*rBmagSqL1[0])*primMomR1[7]+((-246.6576574931337*rBmagSqL1[2])-191.0601999370879*rBmagSqL1[1]-110.3086578651014*rBmagSqL1[0])*primMomL1[7]+((-175.2712184016533*rBmagSqL1[2])-135.7645019878172*rBmagSqL1[1]-78.38367176906175*rBmagSqL1[0])*(primMomR1[6]+primMomL1[6])); 
  incrNonFlux4[2] = 0.00390625*kappaPerp*(((332.0391543176799*(phiR1[2]+phiL1[2])-477.6504998427197*phiR1[1]+477.6504998427197*phiL1[1]+339.411254969543*(phiR1[0]+phiL1[0]))*rBmagSqL1[2]+(257.1964229922337*rBmagSqL1[1]+148.492424049175*rBmagSqL1[0])*(phiR1[2]+phiL1[2])+((-369.9864862397004*phiR1[1])+369.9864862397004*phiL1[1]+262.9068276024798*(phiR1[0]+phiL1[0]))*rBmagSqL1[1]+rBmagSqL1[0]*((-213.6117974270148*phiR1[1])+213.6117974270148*phiL1[1]+151.7893276880823*(phiR1[0]+phiL1[0])))*kperpSq*q_+((-664.0783086353599*rBmagSqL1[2])-514.3928459844674*rBmagSqL1[1]-296.98484809835*rBmagSqL1[0])*(primMomR1[8]+primMomL1[8])+(955.3009996854395*rBmagSqL1[2]+739.9729724794009*rBmagSqL1[1]+427.2235948540296*rBmagSqL1[0])*primMomR1[7]+((-955.3009996854395*rBmagSqL1[2])-739.9729724794009*rBmagSqL1[1]-427.2235948540296*rBmagSqL1[0])*primMomL1[7]+((-678.8225099390861*rBmagSqL1[2])-525.8136552049598*rBmagSqL1[1]-303.5786553761646*rBmagSqL1[0])*(primMomR1[6]+primMomL1[6])); 

  outR[6] += incr3[0]*rdx2R; 
  outR[7] += incrNonFlux3[1]*rdxSq4R+incr3[1]*rdx2R; 
  outR[8] += incrNonFlux3[2]*rdxSq4R+incr3[2]*rdx2R; 

  outL[6] += -1.0*incr3[0]*rdx2L; 
  outL[7] += incr3[1]*rdx2L-1.0*incrNonFlux3[1]*rdxSq4L; 
  outL[8] += incrNonFlux3[2]*rdxSq4L-1.0*incr3[2]*rdx2L; 

  outR[9] += incr4[0]*rdx2R; 
  outR[10] += incrNonFlux4[1]*rdxSq4R+incr4[1]*rdx2R; 
  outR[11] += incrNonFlux4[2]*rdxSq4R+incr4[2]*rdx2R; 

  outL[9] += -1.0*incr4[0]*rdx2L; 
  outL[10] += incr4[1]*rdx2L-1.0*incrNonFlux4[1]*rdxSq4L; 
  outL[11] += incrNonFlux4[2]*rdxSq4L-1.0*incr4[2]*rdx2L; 

}
