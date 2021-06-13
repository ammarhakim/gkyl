#include <gyrofluid_mod_decl.h>

double gyrofluid_surf_1x_p2_ser_x(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *wL1, const double *dxL1, const double *wR1, const double *dxR1, const double cMaxIn, const double *jacL, const double *rBmagL, const double *jacDbmagL, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, double *outL, double *outR) 
{ 
  // q_,m_:              species charge and mass.
  // kappaPar,kappaPerp: heat conductivity coefficients.
  // kperpSq:            k_perp^2.
  // wL,wR:              cell-center in left and right cells.
  // dxL,dxR:            cell length in left and right cells.
  // cMaxIn:             maximum sound speed (or some factor like it).
  // jac:                jacobian.
  // rBmag:              reciprocal of magnetic field magnitude (1/B).
  // jacDbmag:           jacobian divided by B (J/B).
  // sMomL,sMomR:        stepped moments (times Jacobian) in left and right cells.
  // phiL,phiR:          electrostatic potential in left and right cells.
  // primMomL,primMomR:  primitive moments (upar, Tpar, Tperp) in left and right cells.
  // outL/outR:          output increment in left and right cells.

  double wxL = wL1[0];
  double rdx2L = 2.0/dxL1[0];
  double rdxSq4L = rdx2L*rdx2L;
  double wxR = wR1[0];
  double rdx2R = 2.0/dxR1[0];
  double rdxSq4R = rdx2R*rdx2R;

  double uparL[1]; 
  uparL[0] = 0.7071067811865475*(2.23606797749979*primMomL1[2]+1.732050807568877*primMomL1[1]+primMomL1[0]); 

  double uparR[1]; 
  uparR[0] = 0.7071067811865475*(2.23606797749979*primMomR1[2]-1.732050807568877*primMomR1[1]+primMomR1[0]); 

  double cMax = cMaxIn + fabs(fmax(1.58113883008419*primMomL1[2]+1.224744871391589*primMomL1[1]+0.7071067811865475*primMomL1[0],1.58113883008419*primMomR1[2]-1.224744871391589*primMomR1[1]+0.7071067811865475*primMomR1[0]));
  double sMom1Favg[1];
  sMom1Favg[0] = 0.25*(3.16227766016838*(sMomR1[5]+sMomL1[5])-2.449489742783178*sMomR1[4]+2.449489742783178*sMomL1[4]+1.414213562373095*(sMomR1[3]+sMomL1[3])); 

  double momHat1[1];
  momHat1[0] = -0.3535533905932737*((2.23606797749979*sMomR1[2]-2.23606797749979*sMomL1[2]-1.732050807568877*(sMomR1[1]+sMomL1[1])+sMomR1[0]-1.0*sMomL1[0])*cMax-2.828427124746191*sMom1Favg[0]); 

  double sMom2Favg[1];
  sMom2Favg[0] = (0.01428571428571429*((110.6797181058933*(sMomR1[8]+sMomL1[8])-85.73214099741124*sMomR1[7]+85.73214099741124*sMomL1[7]+49.49747468305833*(sMomR1[6]+sMomL1[6]))*m_+((-85.0*sMomR1[2])+54.22176684690384*sMomR1[1]-78.26237921249266*sMomR1[0])*primMomR1[8]+((-85.0*sMomL1[2])-54.22176684690384*sMomL1[1]-78.26237921249266*sMomL1[0])*primMomL1[8]+(54.22176684690384*sMomR1[2]-105.0*sMomR1[1]+60.6217782649107*sMomR1[0])*primMomR1[7]+((-54.22176684690384*sMomL1[2])-105.0*sMomL1[1]-60.6217782649107*sMomL1[0])*primMomL1[7]+((-78.26237921249266*sMomR1[2])+60.6217782649107*sMomR1[1]-35.0*sMomR1[0])*primMomR1[6]+((-78.26237921249266*sMomL1[2])-60.6217782649107*sMomL1[1]-35.0*sMomL1[0])*primMomL1[6]))/m_; 

  double momHat2[1];
  momHat2[0] = -0.3535533905932737*((2.23606797749979*sMomR1[5]-2.23606797749979*sMomL1[5]-1.732050807568877*(sMomR1[4]+sMomL1[4])+sMomR1[3]-1.0*sMomL1[3])*cMax-2.828427124746191*sMom2Favg[0]); 

  double sMom3Favg[1];
  sMom3Favg[0] = (5.102040816326531e-4*(((1190.0*primMomR1[2]-759.1047358566536*primMomR1[1]+1095.673308974897*primMomR1[0])*sMomR1[8]+(1190.0*primMomL1[2]+759.1047358566536*primMomL1[1]+1095.673308974897*primMomL1[0])*sMomL1[8]+((-759.1047358566536*primMomR1[2])+1470.0*primMomR1[1]-848.7048957087499*primMomR1[0])*sMomR1[7]+(759.1047358566536*primMomL1[2]+1470.0*primMomL1[1]+848.7048957087499*primMomL1[0])*sMomL1[7]+(1095.673308974897*primMomR1[2]-848.7048957087499*primMomR1[1]+490.0*primMomR1[0])*sMomR1[6]+(1095.673308974897*primMomL1[2]+848.7048957087499*primMomL1[1]+490.0*primMomL1[0])*sMomL1[6])*m_+((1312.345228969878*primMomR1[2]-943.0535509715237*primMomR1[1]+841.4570696119916*primMomR1[0])*sMomR1[2]+(841.4570696119916*sMomR1[0]-480.099989585503*sMomR1[1])*primMomR1[2]+(929.7096320895039*primMomR1[1]-536.768106355063*primMomR1[0])*sMomR1[1]+sMomR1[0]*(774.7580267412532*primMomR1[0]-536.768106355063*primMomR1[1]))*primMomR1[5]+((1312.345228969878*primMomL1[2]+943.0535509715237*primMomL1[1]+841.4570696119916*primMomL1[0])*sMomL1[2]+(480.099989585503*sMomL1[1]+841.4570696119916*sMomL1[0])*primMomL1[2]+(929.7096320895039*primMomL1[1]+536.768106355063*primMomL1[0])*sMomL1[1]+sMomL1[0]*(536.768106355063*primMomL1[1]+774.7580267412532*primMomL1[0]))*primMomL1[5]+(((-480.099989585503*primMomR1[2])+929.7096320895039*primMomR1[1]-536.768106355063*primMomR1[0])*sMomR1[2]+(1527.380109861327*sMomR1[1]-536.768106355063*sMomR1[0])*primMomR1[2]+(1039.446968344225*primMomR1[0]-1080.224976567381*primMomR1[1])*sMomR1[1]+sMomR1[0]*(1039.446968344225*primMomR1[1]-600.1249869818786*primMomR1[0]))*primMomR1[4]+((480.099989585503*primMomL1[2]+929.7096320895039*primMomL1[1]+536.768106355063*primMomL1[0])*sMomL1[2]+(1527.380109861327*sMomL1[1]+536.768106355063*sMomL1[0])*primMomL1[2]+(1080.224976567381*primMomL1[1]+1039.446968344225*primMomL1[0])*sMomL1[1]+sMomL1[0]*(1039.446968344225*primMomL1[1]+600.1249869818786*primMomL1[0]))*primMomL1[4]+((841.4570696119916*primMomR1[2]-536.768106355063*primMomR1[1]+774.7580267412532*primMomR1[0])*sMomR1[2]+(774.7580267412532*sMomR1[0]-536.768106355063*sMomR1[1])*primMomR1[2]+(1039.446968344225*primMomR1[1]-600.1249869818786*primMomR1[0])*sMomR1[1]+sMomR1[0]*(346.4823227814083*primMomR1[0]-600.1249869818786*primMomR1[1]))*primMomR1[3]+((841.4570696119916*primMomL1[2]+536.768106355063*primMomL1[1]+774.7580267412532*primMomL1[0])*sMomL1[2]+(536.768106355063*sMomL1[1]+774.7580267412532*sMomL1[0])*primMomL1[2]+(1039.446968344225*primMomL1[1]+600.1249869818786*primMomL1[0])*sMomL1[1]+sMomL1[0]*(600.1249869818786*primMomL1[1]+346.4823227814083*primMomL1[0]))*primMomL1[3]))/m_; 

  double momHat3[1];
  momHat3[0] = -0.3535533905932737*((2.23606797749979*sMomR1[8]-2.23606797749979*sMomL1[8]-1.732050807568877*(sMomR1[7]+sMomL1[7])+sMomR1[6]-1.0*sMomL1[6])*cMax-2.828427124746191*sMom3Favg[0]); 

  double sMom4Favg[1];
  sMom4Favg[0] = 0.007142857142857143*((85.0*primMomR1[2]-54.22176684690384*primMomR1[1]+78.26237921249266*primMomR1[0])*sMomR1[11]+(85.0*primMomL1[2]+54.22176684690384*primMomL1[1]+78.26237921249266*primMomL1[0])*sMomL1[11]+((-54.22176684690384*primMomR1[2])+105.0*primMomR1[1]-60.6217782649107*primMomR1[0])*sMomR1[10]+(54.22176684690384*primMomL1[2]+105.0*primMomL1[1]+60.6217782649107*primMomL1[0])*sMomL1[10]+(78.26237921249266*primMomR1[2]-60.6217782649107*primMomR1[1]+35.0*primMomR1[0])*sMomR1[9]+(78.26237921249266*primMomL1[2]+60.6217782649107*primMomL1[1]+35.0*primMomL1[0])*sMomL1[9]); 

  double momHat4[1];
  momHat4[0] = -0.3535533905932737*((2.23606797749979*sMomR1[11]-2.23606797749979*sMomL1[11]-1.732050807568877*(sMomR1[10]+sMomL1[10])+sMomR1[9]-1.0*sMomL1[9])*cMax-2.828427124746191*sMom4Favg[0]); 

  double GheatF3[3];

  GheatF3[0] = -0.00625*(((169.7056274847715*jacL[2]+131.4534138012399*jacL[1]+75.89466384404115*jacL[0])*phiR1[2]+((-169.7056274847715*jacL[2])-131.4534138012399*jacL[1]-75.89466384404115*jacL[0])*phiL1[2]+((-301.2474066278414*phiR1[1])-301.2474066278414*phiL1[1]+237.1708245126285*phiR1[0]-237.1708245126285*phiL1[0])*jacL[2]+((-233.3452377915607*jacL[1])-134.7219358530748*jacL[0])*phiR1[1]+((-233.3452377915607*jacL[1])-134.7219358530748*jacL[0])*phiL1[1]+(183.7117307087383*phiR1[0]-183.7117307087383*phiL1[0])*jacL[1]+106.0660171779821*jacL[0]*phiR1[0]-106.0660171779821*jacL[0]*phiL1[0])*kappaPerp*kperpSq*q_+(((-339.411254969543*jacL[2])-262.9068276024798*jacL[1]-151.7893276880823*jacL[0])*primMomR1[8]+(339.411254969543*jacL[2]+262.9068276024798*jacL[1]+151.7893276880823*jacL[0])*primMomL1[8]+(602.494813255683*jacL[2]+466.6904755831215*jacL[1]+269.4438717061497*jacL[0])*primMomR1[7]+(602.494813255683*jacL[2]+466.6904755831215*jacL[1]+269.4438717061497*jacL[0])*primMomL1[7]+((-474.3416490252571*jacL[2])-367.4234614174767*jacL[1]-212.1320343559643*jacL[0])*primMomR1[6]+(474.3416490252571*jacL[2]+367.4234614174767*jacL[1]+212.1320343559643*jacL[0])*primMomL1[6])*kappaPerp+(((-169.7056274847715*jacL[2])-131.4534138012399*jacL[1]-75.89466384404115*jacL[0])*primMomR1[5]+(169.7056274847715*jacL[2]+131.4534138012399*jacL[1]+75.89466384404115*jacL[0])*primMomL1[5]+(301.2474066278414*jacL[2]+233.3452377915607*jacL[1]+134.7219358530748*jacL[0])*primMomR1[4]+(301.2474066278414*jacL[2]+233.3452377915607*jacL[1]+134.7219358530748*jacL[0])*primMomL1[4]+((-237.1708245126285*jacL[2])-183.7117307087383*jacL[1]-106.0660171779821*jacL[0])*primMomR1[3]+(237.1708245126285*jacL[2]+183.7117307087383*jacL[1]+106.0660171779821*jacL[0])*primMomL1[3])*kappaPar)*rdx2L; 

  double GheatF4[3];

  GheatF4[0] = -0.00625*(((169.7056274847715*jacDbmagL[2]+131.4534138012399*jacDbmagL[1]+75.89466384404115*jacDbmagL[0])*phiR1[2]+((-169.7056274847715*jacDbmagL[2])-131.4534138012399*jacDbmagL[1]-75.89466384404115*jacDbmagL[0])*phiL1[2]+((-301.2474066278414*phiR1[1])-301.2474066278414*phiL1[1]+237.1708245126285*phiR1[0]-237.1708245126285*phiL1[0])*jacDbmagL[2]+((-233.3452377915607*jacDbmagL[1])-134.7219358530748*jacDbmagL[0])*phiR1[1]+((-233.3452377915607*jacDbmagL[1])-134.7219358530748*jacDbmagL[0])*phiL1[1]+(183.7117307087383*phiR1[0]-183.7117307087383*phiL1[0])*jacDbmagL[1]+106.0660171779821*jacDbmagL[0]*phiR1[0]-106.0660171779821*jacDbmagL[0]*phiL1[0])*kappaPerp*kperpSq*q_+(((-339.411254969543*jacDbmagL[2])-262.9068276024798*jacDbmagL[1]-151.7893276880823*jacDbmagL[0])*primMomR1[8]+(339.411254969543*jacDbmagL[2]+262.9068276024798*jacDbmagL[1]+151.7893276880823*jacDbmagL[0])*primMomL1[8]+(602.494813255683*jacDbmagL[2]+466.6904755831215*jacDbmagL[1]+269.4438717061497*jacDbmagL[0])*primMomR1[7]+(602.494813255683*jacDbmagL[2]+466.6904755831215*jacDbmagL[1]+269.4438717061497*jacDbmagL[0])*primMomL1[7]+((-474.3416490252571*jacDbmagL[2])-367.4234614174767*jacDbmagL[1]-212.1320343559643*jacDbmagL[0])*primMomR1[6]+(474.3416490252571*jacDbmagL[2]+367.4234614174767*jacDbmagL[1]+212.1320343559643*jacDbmagL[0])*primMomL1[6])*kappaPerp)*rdx2L; 

  double incr1[3];
  incr1[0] = 0.7071067811865475*momHat1[0]; 
  incr1[1] = -1.224744871391589*momHat1[0]; 
  incr1[2] = 1.58113883008419*momHat1[0]; 

  double incrNonFlux1[3];

  double incr2[3];
  incr2[0] = 0.7071067811865475*momHat2[0]; 
  incr2[1] = -1.224744871391589*momHat2[0]; 
  incr2[2] = 1.58113883008419*momHat2[0]; 

  double incrNonFlux2[3];

  double incr3[3];
  incr3[0] = 0.7071067811865475*momHat3[0]-0.5*GheatF3[0]; 
  incr3[1] = 0.8660254037844386*GheatF3[0]-1.224744871391589*momHat3[0]; 
  incr3[2] = 1.58113883008419*momHat3[0]-1.118033988749895*GheatF3[0]; 

  double incrNonFlux3[3];
  incrNonFlux3[1] = -0.00390625*(kappaPerp*(((85.73214099741124*jacL[2]+66.40783086353598*jacL[1]+38.34057902536163*jacL[0])*(phiR1[2]+phiL1[2])+((-123.3288287465668*phiR1[1])+123.3288287465668*phiL1[1]+87.63560920082662*(phiR1[0]+phiL1[0]))*jacL[2]+((-95.53009996854395*jacL[1])-55.15432893255071*jacL[0])*phiR1[1]+(95.53009996854395*jacL[1]+55.15432893255071*jacL[0])*phiL1[1]+(phiR1[0]+phiL1[0])*(67.8822509939086*jacL[1]+39.19183588453087*jacL[0]))*kperpSq*q_+((-171.4642819948225*jacL[2])-132.815661727072*jacL[1]-76.68115805072327*jacL[0])*(primMomR1[8]+primMomL1[8])+(246.6576574931337*jacL[2]+191.0601999370879*jacL[1]+110.3086578651014*jacL[0])*primMomR1[7]+((-246.6576574931337*jacL[2])-191.0601999370879*jacL[1]-110.3086578651014*jacL[0])*primMomL1[7]+((-175.2712184016533*jacL[2])-135.7645019878172*jacL[1]-78.38367176906175*jacL[0])*(primMomR1[6]+primMomL1[6]))+(((-85.73214099741124*jacL[2])-66.40783086353598*jacL[1]-38.34057902536163*jacL[0])*(primMomR1[5]+primMomL1[5])+(123.3288287465668*jacL[2]+95.53009996854395*jacL[1]+55.15432893255071*jacL[0])*primMomR1[4]+((-123.3288287465668*jacL[2])-95.53009996854395*jacL[1]-55.15432893255071*jacL[0])*primMomL1[4]+((-87.63560920082662*jacL[2])-67.8822509939086*jacL[1]-39.19183588453087*jacL[0])*(primMomR1[3]+primMomL1[3]))*kappaPar); 
  incrNonFlux3[2] = 0.00390625*(kappaPerp*(((332.0391543176799*jacL[2]+257.1964229922337*jacL[1]+148.492424049175*jacL[0])*(phiR1[2]+phiL1[2])+((-477.6504998427197*phiR1[1])+477.6504998427197*phiL1[1]+339.411254969543*(phiR1[0]+phiL1[0]))*jacL[2]+((-369.9864862397004*jacL[1])-213.6117974270148*jacL[0])*phiR1[1]+(369.9864862397004*jacL[1]+213.6117974270148*jacL[0])*phiL1[1]+(phiR1[0]+phiL1[0])*(262.9068276024798*jacL[1]+151.7893276880823*jacL[0]))*kperpSq*q_+((-664.0783086353599*jacL[2])-514.3928459844674*jacL[1]-296.98484809835*jacL[0])*(primMomR1[8]+primMomL1[8])+(955.3009996854395*jacL[2]+739.9729724794009*jacL[1]+427.2235948540296*jacL[0])*primMomR1[7]+((-955.3009996854395*jacL[2])-739.9729724794009*jacL[1]-427.2235948540296*jacL[0])*primMomL1[7]+((-678.8225099390861*jacL[2])-525.8136552049598*jacL[1]-303.5786553761646*jacL[0])*(primMomR1[6]+primMomL1[6]))+(((-332.0391543176799*jacL[2])-257.1964229922337*jacL[1]-148.492424049175*jacL[0])*(primMomR1[5]+primMomL1[5])+(477.6504998427197*jacL[2]+369.9864862397004*jacL[1]+213.6117974270148*jacL[0])*primMomR1[4]+((-477.6504998427197*jacL[2])-369.9864862397004*jacL[1]-213.6117974270148*jacL[0])*primMomL1[4]+((-339.411254969543*jacL[2])-262.9068276024798*jacL[1]-151.7893276880823*jacL[0])*(primMomR1[3]+primMomL1[3]))*kappaPar); 

  double incr4[3];
  incr4[0] = 0.7071067811865475*momHat4[0]-0.5*GheatF4[0]; 
  incr4[1] = 0.8660254037844386*GheatF4[0]-1.224744871391589*momHat4[0]; 
  incr4[2] = 1.58113883008419*momHat4[0]-1.118033988749895*GheatF4[0]; 

  double incrNonFlux4[3];
  incrNonFlux4[1] = -0.00390625*kappaPerp*(((85.73214099741124*jacDbmagL[2]+66.40783086353598*jacDbmagL[1]+38.34057902536163*jacDbmagL[0])*(phiR1[2]+phiL1[2])+((-123.3288287465668*phiR1[1])+123.3288287465668*phiL1[1]+87.63560920082662*(phiR1[0]+phiL1[0]))*jacDbmagL[2]+((-95.53009996854395*jacDbmagL[1])-55.15432893255071*jacDbmagL[0])*phiR1[1]+(95.53009996854395*jacDbmagL[1]+55.15432893255071*jacDbmagL[0])*phiL1[1]+(phiR1[0]+phiL1[0])*(67.8822509939086*jacDbmagL[1]+39.19183588453087*jacDbmagL[0]))*kperpSq*q_+((-171.4642819948225*jacDbmagL[2])-132.815661727072*jacDbmagL[1]-76.68115805072327*jacDbmagL[0])*(primMomR1[8]+primMomL1[8])+(246.6576574931337*jacDbmagL[2]+191.0601999370879*jacDbmagL[1]+110.3086578651014*jacDbmagL[0])*primMomR1[7]+((-246.6576574931337*jacDbmagL[2])-191.0601999370879*jacDbmagL[1]-110.3086578651014*jacDbmagL[0])*primMomL1[7]+((-175.2712184016533*jacDbmagL[2])-135.7645019878172*jacDbmagL[1]-78.38367176906175*jacDbmagL[0])*(primMomR1[6]+primMomL1[6])); 
  incrNonFlux4[2] = 0.00390625*kappaPerp*(((332.0391543176799*jacDbmagL[2]+257.1964229922337*jacDbmagL[1]+148.492424049175*jacDbmagL[0])*(phiR1[2]+phiL1[2])+((-477.6504998427197*phiR1[1])+477.6504998427197*phiL1[1]+339.411254969543*(phiR1[0]+phiL1[0]))*jacDbmagL[2]+((-369.9864862397004*jacDbmagL[1])-213.6117974270148*jacDbmagL[0])*phiR1[1]+(369.9864862397004*jacDbmagL[1]+213.6117974270148*jacDbmagL[0])*phiL1[1]+(phiR1[0]+phiL1[0])*(262.9068276024798*jacDbmagL[1]+151.7893276880823*jacDbmagL[0]))*kperpSq*q_+((-664.0783086353599*jacDbmagL[2])-514.3928459844674*jacDbmagL[1]-296.98484809835*jacDbmagL[0])*(primMomR1[8]+primMomL1[8])+(955.3009996854395*jacDbmagL[2]+739.9729724794009*jacDbmagL[1]+427.2235948540296*jacDbmagL[0])*primMomR1[7]+((-955.3009996854395*jacDbmagL[2])-739.9729724794009*jacDbmagL[1]-427.2235948540296*jacDbmagL[0])*primMomL1[7]+((-678.8225099390861*jacDbmagL[2])-525.8136552049598*jacDbmagL[1]-303.5786553761646*jacDbmagL[0])*(primMomR1[6]+primMomL1[6])); 

  outR[0] += incr1[0]*rdx2R; 
  outR[1] += incr1[1]*rdx2R; 
  outR[2] += incr1[2]*rdx2R; 

  outL[0] += -1.0*incr1[0]*rdx2L; 
  outL[1] += incr1[1]*rdx2L; 
  outL[2] += -1.0*incr1[2]*rdx2L; 

  outR[3] += incr2[0]*rdx2R; 
  outR[4] += incr2[1]*rdx2R; 
  outR[5] += incr2[2]*rdx2R; 

  outL[3] += -1.0*incr2[0]*rdx2L; 
  outL[4] += incr2[1]*rdx2L; 
  outL[5] += -1.0*incr2[2]*rdx2L; 

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

  return 0.5*sqrt((-(6.324555320336761*primMomL1[8])/m_)+(5.656854249492382*primMomL1[6])/m_-(9.48683298050514*primMomL1[5])/m_+(8.485281374238571*primMomL1[3])/m_); 
}
