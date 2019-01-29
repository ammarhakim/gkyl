#include <GkLBOModDecl.h> 
double GkLBOconstNuVol1x1vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[2]:    Cell-center coordinates. 
  // dxv[2]:  Cell spacing. 
  // nu:      diffusion coefficient (collisionality). 
  // u[2]:    bulk velocity. 
  // vtSq[2]: thermal speed squared. 
  // f[4]:    Input distribution function. 
  // out[4]:  Incremented output 
  double rdv2nu[1]; 
  double rdvSq4nu[1]; 
  rdv2nu[0] = 2.0/dxv[1]; 
  rdvSq4nu[0] = nu*rdv2nu[0]*rdv2nu[0]; 
  rdv2nu[0] = nu*rdv2nu[0]; 

  double alphaVpar[4]; 
  alphaVpar[0] = rdv2nu[0]*(1.414213562373095*u[0]-2.0*w[1]); 
  alphaVpar[1] = 1.414213562373095*rdv2nu[0]*u[1]; 
  alphaVpar[2] = -1.154700538379252*nu; 

  out[2] += 0.8660254037844386*(alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphaVpar[2]*f[3]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 

  return std::abs(0.5*alphaVpar[0]) + 0.9428090415820625*rdvSq4nu[0]*vtSq[0]; 

} 
double GkLBOconstNuVol1x1vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[2]:    Cell-center coordinates. 
  // dxv[2]:  Cell spacing. 
  // nu:      diffusion coefficient (collisionality). 
  // u[3]:    bulk velocity. 
  // vtSq[3]: thermal speed squared. 
  // f[8]:    Input distribution function. 
  // out[8]:  Incremented output 
  double rdv2nu[1]; 
  double rdvSq4nu[1]; 
  rdv2nu[0] = 2.0/dxv[1]; 
  rdvSq4nu[0] = nu*rdv2nu[0]*rdv2nu[0]; 
  rdv2nu[0] = nu*rdv2nu[0]; 

  double alphaVpar[8]; 
  alphaVpar[0] = rdv2nu[0]*(1.414213562373095*u[0]-2.0*w[1]); 
  alphaVpar[1] = 1.414213562373095*rdv2nu[0]*u[1]; 
  alphaVpar[2] = -1.154700538379252*nu; 
  alphaVpar[4] = 1.414213562373095*rdv2nu[0]*u[2]; 

  double facVpar[3]; 
  facVpar[0] = rdvSq4nu[0]*vtSq[0]; 
  facVpar[1] = rdvSq4nu[0]*vtSq[1]; 
  facVpar[2] = rdvSq4nu[0]*vtSq[2]; 

  out[2] += 0.8660254037844386*(alphaVpar[4]*f[4]+alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[3] += 0.7745966692414833*(alphaVpar[1]*f[4]+f[1]*alphaVpar[4])+0.8660254037844386*(alphaVpar[2]*f[3]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 
  out[5] += 1.936491673103709*alphaVpar[4]*f[6]+1.732050807568877*alphaVpar[2]*f[5]+4.743416490252569*facVpar[2]*f[4]+1.936491673103709*(alphaVpar[1]*f[3]+alphaVpar[0]*f[2]+f[0]*alphaVpar[2])+4.743416490252569*(f[1]*facVpar[1]+f[0]*facVpar[0]); 
  out[6] += 0.8660254037844386*alphaVpar[2]*f[6]+0.5532833351724881*alphaVpar[4]*f[4]+0.8660254037844386*(alphaVpar[0]*f[4]+f[0]*alphaVpar[4])+0.7745966692414833*alphaVpar[1]*f[1]; 
  out[7] += 1.732050807568877*(alphaVpar[2]*f[7]+alphaVpar[1]*f[6])+4.242640687119286*facVpar[1]*f[4]+f[3]*(1.732050807568877*alphaVpar[4]+1.936491673103709*alphaVpar[0])+4.242640687119286*f[1]*facVpar[2]+1.936491673103709*(alphaVpar[1]*f[2]+f[1]*alphaVpar[2])+4.743416490252569*(f[0]*facVpar[1]+facVpar[0]*f[1]); 

  return std::abs(0.5*alphaVpar[0]-0.5590169943749475*alphaVpar[4]) + 1.272792206135784*facVpar[0]-1.42302494707577*facVpar[2]; 

} 
double GkLBOconstNuVol1x1vSerP3(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[2]:    Cell-center coordinates. 
  // dxv[2]:  Cell spacing. 
  // nu:      diffusion coefficient (collisionality). 
  // u[4]:    bulk velocity. 
  // vtSq[4]: thermal speed squared. 
  // f[12]:    Input distribution function. 
  // out[12]:  Incremented output 
  double rdv2nu[1]; 
  double rdvSq4nu[1]; 
  rdv2nu[0] = 2.0/dxv[1]; 
  rdvSq4nu[0] = nu*rdv2nu[0]*rdv2nu[0]; 
  rdv2nu[0] = nu*rdv2nu[0]; 

  double alphaVpar[12]; 
  alphaVpar[0] = rdv2nu[0]*(1.414213562373095*u[0]-2.0*w[1]); 
  alphaVpar[1] = 1.414213562373095*rdv2nu[0]*u[1]; 
  alphaVpar[2] = -1.154700538379252*nu; 
  alphaVpar[4] = 1.414213562373095*rdv2nu[0]*u[2]; 
  alphaVpar[8] = 1.414213562373095*rdv2nu[0]*u[3]; 

  double facVpar[4]; 
  facVpar[0] = rdvSq4nu[0]*vtSq[0]; 
  facVpar[1] = rdvSq4nu[0]*vtSq[1]; 
  facVpar[2] = rdvSq4nu[0]*vtSq[2]; 
  facVpar[3] = rdvSq4nu[0]*vtSq[3]; 

  out[2] += 0.8660254037844386*(alphaVpar[8]*f[8]+alphaVpar[4]*f[4]+alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[3] += 0.7606388292556648*(alphaVpar[4]*f[8]+f[4]*alphaVpar[8])+0.7745966692414833*(alphaVpar[1]*f[4]+f[1]*alphaVpar[4])+0.8660254037844386*(alphaVpar[2]*f[3]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 
  out[5] += 1.936491673103709*alphaVpar[8]*f[10]+4.743416490252569*facVpar[3]*f[8]+1.936491673103709*alphaVpar[4]*f[6]+1.732050807568877*alphaVpar[2]*f[5]+4.743416490252569*facVpar[2]*f[4]+1.936491673103709*(alphaVpar[1]*f[3]+alphaVpar[0]*f[2]+f[0]*alphaVpar[2])+4.743416490252569*(f[1]*facVpar[1]+f[0]*facVpar[0]); 
  out[6] += 0.5163977794943223*alphaVpar[8]*f[8]+0.7606388292556648*(alphaVpar[1]*f[8]+f[1]*alphaVpar[8])+0.8660254037844386*alphaVpar[2]*f[6]+0.5532833351724881*alphaVpar[4]*f[4]+0.8660254037844386*(alphaVpar[0]*f[4]+f[0]*alphaVpar[4])+0.7745966692414833*alphaVpar[1]*f[1]; 
  out[7] += 1.700840128541522*alphaVpar[4]*f[10]+4.166190448976479*facVpar[2]*f[8]+1.700840128541522*f[6]*alphaVpar[8]+1.732050807568877*(alphaVpar[2]*f[7]+alphaVpar[1]*f[6])+(4.166190448976479*facVpar[3]+4.242640687119286*facVpar[1])*f[4]+f[3]*(1.732050807568877*alphaVpar[4]+1.936491673103709*alphaVpar[0])+4.242640687119286*f[1]*facVpar[2]+1.936491673103709*(alphaVpar[1]*f[2]+f[1]*alphaVpar[2])+4.743416490252569*(f[0]*facVpar[1]+facVpar[0]*f[1]); 
  out[9] += 16.20185174601965*facVpar[3]*f[10]+2.598076211353316*alphaVpar[2]*f[9]+1.322875655532295*alphaVpar[8]*f[8]+2.958039891549809*alphaVpar[1]*f[7]+16.20185174601965*facVpar[2]*f[6]+2.958039891549809*alphaVpar[0]*f[5]+1.322875655532295*alphaVpar[4]*f[4]+16.20185174601965*facVpar[1]*f[3]+(3.968626966596886*alphaVpar[2]+16.20185174601965*facVpar[0])*f[2]+1.322875655532295*(alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[10] += 0.8660254037844386*alphaVpar[2]*f[10]+(0.5163977794943223*alphaVpar[4]+0.8660254037844386*alphaVpar[0])*f[8]+(0.5163977794943223*f[4]+0.8660254037844386*f[0])*alphaVpar[8]+0.7606388292556648*(alphaVpar[1]*f[4]+f[1]*alphaVpar[4]); 
  out[11] += 2.598076211353316*alphaVpar[2]*f[11]+14.23024947075771*facVpar[2]*f[10]+1.161895003862225*(alphaVpar[4]*f[8]+f[4]*alphaVpar[8])+(2.645751311064591*alphaVpar[4]+2.958039891549809*alphaVpar[0])*f[7]+(14.23024947075771*facVpar[3]+14.49137674618944*facVpar[1])*f[6]+2.958039891549809*alphaVpar[1]*f[5]+1.183215956619923*(alphaVpar[1]*f[4]+f[1]*alphaVpar[4])+(14.49137674618944*facVpar[2]+3.968626966596886*alphaVpar[2])*f[3]+16.20185174601965*(facVpar[0]*f[3]+facVpar[1]*f[2])+1.322875655532295*(alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 

  return std::abs(0.5*alphaVpar[0]-0.5590169943749475*alphaVpar[4]) + 1.616244071283536*facVpar[0]-1.807015805810502*facVpar[2]; 

} 
