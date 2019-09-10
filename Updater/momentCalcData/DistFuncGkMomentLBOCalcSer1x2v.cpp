#include <DistFuncMomentCalcModDecl.h> 
#include <cmath> 
void GkM0Star1x2vSer_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[3]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[2]*intFac*(wr[1]-wl[1]); 
 
  out[0] += ((-0.5773502691896258*fr[2])+0.5773502691896258*fl[2]+0.5*fr[0]+0.5*fl[0])*dS; 
  out[1] += ((-0.5773502691896258*fr[4])+0.5773502691896258*fl[4]+0.5*fr[1]+0.5*fl[1])*dS; 
 
} 
 
void GkM0StarPositivity1x2vSer_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[3]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[2]*intFac*(wr[1]-wl[1]); 
 
  if ( ((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]+0.05892556509887893*fr[5]+0.05892556509887893*fl[5]+0.1178511301977579*fr[4]-0.1178511301977579*fl[4]-0.1020620726159657*fr[3]-0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && (0.06804138174397717*fr[7]-0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]-0.05892556509887893*fr[5]-0.05892556509887893*fl[5]+0.1178511301977579*fr[4]-0.1178511301977579*fl[4]+0.1020620726159657*fr[3]+0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && (0.06804138174397717*fr[7]-0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]-0.05892556509887893*fr[5]-0.05892556509887893*fl[5]-0.1178511301977579*fr[4]+0.1178511301977579*fl[4]-0.1020620726159657*fr[3]-0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && ((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]+0.05892556509887893*fr[5]+0.05892556509887893*fl[5]-0.1178511301977579*fr[4]+0.1178511301977579*fl[4]+0.1020620726159657*fr[3]+0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) ) {
    out[0] += ((-0.5773502691896258*fr[2])+0.5773502691896258*fl[2]+0.5*fr[0]+0.5*fl[0])*dS; 
    out[1] += ((-0.5773502691896258*fr[4])+0.5773502691896258*fl[4]+0.5*fr[1]+0.5*fl[1])*dS; 
  } else {

    double xBar[4];
    xBar[0] = (0.05103103630798284*fr[7]+0.05103103630798284*fl[7]-0.08838834764831849*fr[6]-0.08838834764831849*fl[6]+0.08838834764831843*fr[5]-0.08838834764831843*fl[5]-0.08838834764831849*fr[4]-0.08838834764831849*fl[4]-0.1530931089239485*fr[3]+0.1530931089239485*fl[3]+0.1530931089239485*fr[2]+0.1530931089239485*fl[2]-0.1530931089239485*fr[1]+0.1530931089239485*fl[1]+0.2651650429449551*fr[0]-0.2651650429449551*fl[0])/(0.5*(1.224744871391589*((-0.3333333333333333*(1.732050807568877*fr[4]-1.0*fr[7]))-0.5773502691896258*fr[6]+fr[2])-1.224744871391589*((-0.3333333333333333*(1.732050807568877*fl[4]-1.0*fl[7]))-0.5773502691896258*fl[6]+fl[2]))-0.25*(2.449489742783178*((-0.3333333333333333*(1.732050807568877*fr[4]-1.0*fr[7]))-0.5773502691896258*fr[6]+fr[2])-2.449489742783178*((-0.3333333333333333*(1.732050807568877*fl[4]-1.0*fl[7]))-0.5773502691896258*fl[6]+fl[2])-2.121320343559642*((-0.3333333333333333*(1.732050807568877*fr[1]-1.0*fr[5]))-0.5773502691896258*fr[3]+fr[0])-2.121320343559642*((-0.3333333333333333*(1.732050807568877*fl[1]-1.0*fl[5]))-0.5773502691896258*fl[3]+fl[0]))); 
    xBar[1] = ((-0.05103103630798284*fr[7])-0.05103103630798284*fl[7]+0.08838834764831849*fr[6]+0.08838834764831849*fl[6]-0.08838834764831843*fr[5]+0.08838834764831843*fl[5]-0.08838834764831849*fr[4]-0.08838834764831849*fl[4]+0.1530931089239485*fr[3]-0.1530931089239485*fl[3]+0.1530931089239485*fr[2]+0.1530931089239485*fl[2]-0.1530931089239485*fr[1]+0.1530931089239485*fl[1]+0.2651650429449551*fr[0]-0.2651650429449551*fl[0])/(0.5*(1.224744871391589*((-0.3333333333333333*(fr[7]+1.732050807568877*fr[4]))+0.5773502691896258*fr[6]+fr[2])-1.224744871391589*((-0.3333333333333333*(fl[7]+1.732050807568877*fl[4]))+0.5773502691896258*fl[6]+fl[2]))-0.25*(2.449489742783178*((-0.3333333333333333*(fr[7]+1.732050807568877*fr[4]))+0.5773502691896258*fr[6]+fr[2])-2.449489742783178*((-0.3333333333333333*(fl[7]+1.732050807568877*fl[4]))+0.5773502691896258*fl[6]+fl[2])-2.121320343559642*((-0.3333333333333333*(fr[5]+1.732050807568877*fr[1]))+0.5773502691896258*fr[3]+fr[0])-2.121320343559642*((-0.3333333333333333*(fl[5]+1.732050807568877*fl[1]))+0.5773502691896258*fl[3]+fl[0]))); 
    xBar[2] = ((-0.05103103630798284*fr[7])-0.05103103630798284*fl[7]-0.08838834764831849*fr[6]-0.08838834764831849*fl[6]-0.08838834764831843*fr[5]+0.08838834764831843*fl[5]+0.08838834764831849*fr[4]+0.08838834764831849*fl[4]-0.1530931089239485*fr[3]+0.1530931089239485*fl[3]+0.1530931089239485*fr[2]+0.1530931089239485*fl[2]+0.1530931089239485*fr[1]-0.1530931089239485*fl[1]+0.2651650429449551*fr[0]-0.2651650429449551*fl[0])/(0.5*(1.224744871391589*(0.3333333333333333*(1.732050807568877*fr[4]-1.0*fr[7])-0.5773502691896258*fr[6]+fr[2])-1.224744871391589*(0.3333333333333333*(1.732050807568877*fl[4]-1.0*fl[7])-0.5773502691896258*fl[6]+fl[2]))-0.25*(2.449489742783178*(0.3333333333333333*(1.732050807568877*fr[4]-1.0*fr[7])-0.5773502691896258*fr[6]+fr[2])-2.449489742783178*(0.3333333333333333*(1.732050807568877*fl[4]-1.0*fl[7])-0.5773502691896258*fl[6]+fl[2])-2.121320343559642*(0.3333333333333333*(1.732050807568877*fr[1]-1.0*fr[5])-0.5773502691896258*fr[3]+fr[0])-2.121320343559642*(0.3333333333333333*(1.732050807568877*fl[1]-1.0*fl[5])-0.5773502691896258*fl[3]+fl[0]))); 
    xBar[3] = (0.05103103630798284*fr[7]+0.05103103630798284*fl[7]+0.08838834764831849*fr[6]+0.08838834764831849*fl[6]+0.08838834764831843*fr[5]-0.08838834764831843*fl[5]+0.08838834764831849*fr[4]+0.08838834764831849*fl[4]+0.1530931089239485*fr[3]-0.1530931089239485*fl[3]+0.1530931089239485*fr[2]+0.1530931089239485*fl[2]+0.1530931089239485*fr[1]-0.1530931089239485*fl[1]+0.2651650429449551*fr[0]-0.2651650429449551*fl[0])/(0.5*(1.224744871391589*(0.3333333333333333*(fr[7]+1.732050807568877*fr[4])+0.5773502691896258*fr[6]+fr[2])-1.224744871391589*(0.3333333333333333*(fl[7]+1.732050807568877*fl[4])+0.5773502691896258*fl[6]+fl[2]))-0.25*(2.449489742783178*(0.3333333333333333*(fr[7]+1.732050807568877*fr[4])+0.5773502691896258*fr[6]+fr[2])-2.449489742783178*(0.3333333333333333*(fl[7]+1.732050807568877*fl[4])+0.5773502691896258*fl[6]+fl[2])-2.121320343559642*(0.3333333333333333*(fr[5]+1.732050807568877*fr[1])+0.5773502691896258*fr[3]+fr[0])-2.121320343559642*(0.3333333333333333*(fl[5]+1.732050807568877*fl[1])+0.5773502691896258*fl[3]+fl[0]))); 

    double xBarSq[4];
    xBarSq[0] = xBar[0]*xBar[0]; 
    xBarSq[1] = xBar[1]*xBar[1]; 
    xBarSq[2] = xBar[2]*xBar[2]; 
    xBarSq[3] = xBar[3]*xBar[3]; 

    double g1[4];
    g1[0] = (3.0*xBar[0])/(1.0-1.0*xBarSq[0])-(1.0*xBar[0]*xBarSq[0])/(1.0-1.0*xBarSq[0]); 
    g1[1] = (3.0*xBar[1])/(1.0-1.0*xBarSq[1])-(1.0*xBar[1]*xBarSq[1])/(1.0-1.0*xBarSq[1]); 
    g1[2] = (3.0*xBar[2])/(1.0-1.0*xBarSq[2])-(1.0*xBar[2]*xBarSq[2])/(1.0-1.0*xBarSq[2]); 
    g1[3] = (3.0*xBar[3])/(1.0-1.0*xBarSq[3])-(1.0*xBar[3]*xBarSq[3])/(1.0-1.0*xBarSq[3]); 

    double gBound[4];

    if (std::abs(g1[0]) > 1.0e-15) {
      double g1Sq = g1[0]*g1[0];
      gBound[0] = (-(1.387778780781446e-17*g1[0]*fr[7])/std::sinh(g1[0]))+(1.387778780781446e-17*g1[0]*fl[7])/std::sinh(g1[0])+(2.775557561562891e-17*g1[0]*fr[6])/std::sinh(g1[0])-(2.775557561562891e-17*g1[0]*fl[6])/std::sinh(g1[0])+(0.05892556509887895*g1[0]*fr[5])/std::sinh(g1[0])+(0.05892556509887895*g1[0]*fl[5])/std::sinh(g1[0])+(1.387778780781446e-17*g1[0]*fr[4])/std::sinh(g1[0])-(1.387778780781446e-17*g1[0]*fl[4])/std::sinh(g1[0])-(0.1020620726159658*g1[0]*fr[3])/std::sinh(g1[0])-(0.1020620726159658*g1[0]*fl[3])/std::sinh(g1[0])-(2.775557561562891e-17*g1[0]*fr[2])/std::sinh(g1[0])+(2.775557561562891e-17*g1[0]*fl[2])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fr[1])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fl[1])/std::sinh(g1[0])+(0.1767766952966369*fr[0]*g1[0])/std::sinh(g1[0])+(0.1767766952966369*fl[0]*g1[0])/std::sinh(g1[0]); 
    } else {
      gBound[0] = (-1.387778780781446e-17*fr[7])+1.387778780781446e-17*fl[7]+2.775557561562891e-17*fr[6]-2.775557561562891e-17*fl[6]+0.05892556509887895*fr[5]+0.05892556509887895*fl[5]+1.387778780781446e-17*fr[4]-1.387778780781446e-17*fl[4]-0.1020620726159658*fr[3]-0.1020620726159658*fl[3]-2.775557561562891e-17*fr[2]+2.775557561562891e-17*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966369*fr[0]+0.1767766952966369*fl[0]; 
    };

    if (std::abs(g1[1]) > 1.0e-15) {
      double g1Sq = g1[1]*g1[1];
      gBound[1] = (1.387778780781446e-17*g1[1]*fr[7])/std::sinh(g1[1])-(1.387778780781446e-17*g1[1]*fl[7])/std::sinh(g1[1])-(2.775557561562891e-17*g1[1]*fr[6])/std::sinh(g1[1])+(2.775557561562891e-17*g1[1]*fl[6])/std::sinh(g1[1])-(0.05892556509887895*g1[1]*fr[5])/std::sinh(g1[1])-(0.05892556509887895*g1[1]*fl[5])/std::sinh(g1[1])+(1.387778780781446e-17*g1[1]*fr[4])/std::sinh(g1[1])-(1.387778780781446e-17*g1[1]*fl[4])/std::sinh(g1[1])+(0.1020620726159658*g1[1]*fr[3])/std::sinh(g1[1])+(0.1020620726159658*g1[1]*fl[3])/std::sinh(g1[1])-(2.775557561562891e-17*g1[1]*fr[2])/std::sinh(g1[1])+(2.775557561562891e-17*g1[1]*fl[2])/std::sinh(g1[1])-(0.1020620726159657*fr[1]*g1[1])/std::sinh(g1[1])-(0.1020620726159657*fl[1]*g1[1])/std::sinh(g1[1])+(0.1767766952966369*fr[0]*g1[1])/std::sinh(g1[1])+(0.1767766952966369*fl[0]*g1[1])/std::sinh(g1[1]); 
    } else {
      gBound[1] = 1.387778780781446e-17*fr[7]-1.387778780781446e-17*fl[7]-2.775557561562891e-17*fr[6]+2.775557561562891e-17*fl[6]-0.05892556509887895*fr[5]-0.05892556509887895*fl[5]+1.387778780781446e-17*fr[4]-1.387778780781446e-17*fl[4]+0.1020620726159658*fr[3]+0.1020620726159658*fl[3]-2.775557561562891e-17*fr[2]+2.775557561562891e-17*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966369*fr[0]+0.1767766952966369*fl[0]; 
    };

    if (std::abs(g1[2]) > 1.0e-15) {
      double g1Sq = g1[2]*g1[2];
      gBound[2] = (1.387778780781446e-17*g1[2]*fr[7])/std::sinh(g1[2])-(1.387778780781446e-17*g1[2]*fl[7])/std::sinh(g1[2])+(2.775557561562891e-17*g1[2]*fr[6])/std::sinh(g1[2])-(2.775557561562891e-17*g1[2]*fl[6])/std::sinh(g1[2])-(0.05892556509887895*g1[2]*fr[5])/std::sinh(g1[2])-(0.05892556509887895*g1[2]*fl[5])/std::sinh(g1[2])-(1.387778780781446e-17*g1[2]*fr[4])/std::sinh(g1[2])+(1.387778780781446e-17*g1[2]*fl[4])/std::sinh(g1[2])-(0.1020620726159658*g1[2]*fr[3])/std::sinh(g1[2])-(0.1020620726159658*g1[2]*fl[3])/std::sinh(g1[2])-(2.775557561562891e-17*fr[2]*g1[2])/std::sinh(g1[2])+(2.775557561562891e-17*fl[2]*g1[2])/std::sinh(g1[2])+(0.1020620726159657*fr[1]*g1[2])/std::sinh(g1[2])+(0.1020620726159657*fl[1]*g1[2])/std::sinh(g1[2])+(0.1767766952966369*fr[0]*g1[2])/std::sinh(g1[2])+(0.1767766952966369*fl[0]*g1[2])/std::sinh(g1[2]); 
    } else {
      gBound[2] = 1.387778780781446e-17*fr[7]-1.387778780781446e-17*fl[7]+2.775557561562891e-17*fr[6]-2.775557561562891e-17*fl[6]-0.05892556509887895*fr[5]-0.05892556509887895*fl[5]-1.387778780781446e-17*fr[4]+1.387778780781446e-17*fl[4]-0.1020620726159658*fr[3]-0.1020620726159658*fl[3]-2.775557561562891e-17*fr[2]+2.775557561562891e-17*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966369*fr[0]+0.1767766952966369*fl[0]; 
    };

    if (std::abs(g1[3]) > 1.0e-15) {
      double g1Sq = g1[3]*g1[3];
      gBound[3] = (-(1.387778780781446e-17*g1[3]*fr[7])/std::sinh(g1[3]))+(1.387778780781446e-17*g1[3]*fl[7])/std::sinh(g1[3])-(2.775557561562891e-17*g1[3]*fr[6])/std::sinh(g1[3])+(2.775557561562891e-17*g1[3]*fl[6])/std::sinh(g1[3])+(0.05892556509887895*g1[3]*fr[5])/std::sinh(g1[3])+(0.05892556509887895*g1[3]*fl[5])/std::sinh(g1[3])-(1.387778780781446e-17*g1[3]*fr[4])/std::sinh(g1[3])+(1.387778780781446e-17*g1[3]*fl[4])/std::sinh(g1[3])+(0.1020620726159658*fr[3]*g1[3])/std::sinh(g1[3])+(0.1020620726159658*fl[3]*g1[3])/std::sinh(g1[3])-(2.775557561562891e-17*fr[2]*g1[3])/std::sinh(g1[3])+(2.775557561562891e-17*fl[2]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fr[1]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fl[1]*g1[3])/std::sinh(g1[3])+(0.1767766952966369*fr[0]*g1[3])/std::sinh(g1[3])+(0.1767766952966369*fl[0]*g1[3])/std::sinh(g1[3]); 
    } else {
      gBound[3] = (-1.387778780781446e-17*fr[7])+1.387778780781446e-17*fl[7]-2.775557561562891e-17*fr[6]+2.775557561562891e-17*fl[6]+0.05892556509887895*fr[5]+0.05892556509887895*fl[5]-1.387778780781446e-17*fr[4]+1.387778780781446e-17*fl[4]+0.1020620726159658*fr[3]+0.1020620726159658*fl[3]-2.775557561562891e-17*fr[2]+2.775557561562891e-17*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966369*fr[0]+0.1767766952966369*fl[0]; 
    };

    out[0] += (0.7071067811865475*gBound[3]+0.7071067811865475*gBound[2]+0.7071067811865475*gBound[1]+0.7071067811865475*gBound[0])*dS; 
    out[1] += (1.224744871391589*gBound[3]+1.224744871391589*gBound[2]-1.224744871391589*gBound[1]-1.224744871391589*gBound[0])*dS; 
  };
 
} 
 
void GkM1iM2Star1x2vSer(const double *w, const double *dxv, const double intFac, const double m_, const double *Bmag, const double *f, double *outM1i, double *outM2) 
{ 
  // w[3]:    Cell-center coordinates. 
  // dxv[3]:  Cell length in each direciton. 
  // intFac:  =2pi/m for gyrokinetics. 
  // m_:      mass. 
  // Bmag[2]: Magnetic field magnitude. 
  // f:       Distribution function. 
  // outM1i:  Contribution to M_1^star from this cell. 
  // outM2:   Contribution to M_2^star from this cell. 
 
  const double volFact = intFac*0.25*dxv[1]*dxv[2]; 
  double wvSq[2]; 
  wvSq[0]  = w[1]*w[1]; 
  wvSq[1]  = w[2]*w[2]; 
  double dvSq[2]; 
  dvSq[0] = dxv[1]*dxv[1]; 
  dvSq[1] = dxv[2]*dxv[2]; 
 
  outM1i[0] += 2.0*f[0]*w[1]*volFact; 
  outM1i[1] += 2.0*f[1]*w[1]*volFact; 
 
  double tmp[2]; 
  tmp[0] = 0.5773502691896258*dxv[2]*f[3]+2.0*f[0]*w[2]; 
  tmp[1] = 0.5773502691896258*dxv[2]*f[5]+2.0*f[1]*w[2]; 
 
  outM2[0] += ((1.414213562373095*Bmag[1]*tmp[1]+1.414213562373095*Bmag[0]*tmp[0])/m_+0.5773502691896258*dxv[1]*w[1]*f[2]+2.0*f[0]*wvSq[0])*volFact; 
  outM2[1] += ((1.414213562373095*Bmag[0]*tmp[1]+1.414213562373095*tmp[0]*Bmag[1])/m_+0.5773502691896258*dxv[1]*w[1]*f[4]+2.0*f[1]*wvSq[0])*volFact; 
 
} 
void GkBoundaryIntegral1x2vSer_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[8]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[2]+fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[4]+fIn[1])*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vSer_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[20]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[8])+1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
    out[1] += ((-2.23606797749979*fIn[12])+1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[11]-1.0*fIn[7])*dS; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS; 
    out[1] += (2.23606797749979*fIn[12]+1.732050807568877*fIn[4]+fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[11]+fIn[7])*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vSer_F_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[18]-2.23606797749979*fIn[8]+1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
    out[1] += (2.645751311064591*fIn[24]-2.23606797749979*fIn[12]+1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[11]-1.0*fIn[7])*dS; 
    out[3] += (1.732050807568877*fIn[23]-1.0*fIn[17])*dS; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[18]+2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS; 
    out[1] += (2.645751311064591*fIn[24]+2.23606797749979*fIn[12]+1.732050807568877*fIn[4]+fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[11]+fIn[7])*dS; 
    out[3] += (1.732050807568877*fIn[23]+fIn[17])*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vSer_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[8]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[2]-1.0*fIn[0])*dS*vBoundary+(0.8660254037844386*dxv[1]*fIn[2]-0.5*fIn[0]*dxv[1])*dS; 
    out[1] += (1.732050807568877*fIn[4]-1.0*fIn[1])*dS*vBoundary+(0.8660254037844386*dxv[1]*fIn[4]-0.5*dxv[1]*fIn[1])*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[2]+fIn[0])*dS*vBoundary+((-0.8660254037844386*dxv[1]*fIn[2])-0.5*fIn[0]*dxv[1])*dS; 
    out[1] += (1.732050807568877*fIn[4]+fIn[1])*dS*vBoundary+((-0.8660254037844386*dxv[1]*fIn[4])-0.5*dxv[1]*fIn[1])*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vSer_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[20]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[8])+1.732050807568877*fIn[2]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += ((-2.23606797749979*fIn[12])+1.732050807568877*fIn[4]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[11]-1.0*fIn[7])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS*vBoundary; 
    out[1] += (2.23606797749979*fIn[12]+1.732050807568877*fIn[4]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[11]+fIn[7])*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vSer_vF_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[18]-2.23606797749979*fIn[8]+1.732050807568877*fIn[2]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (2.645751311064591*fIn[24]-2.23606797749979*fIn[12]+1.732050807568877*fIn[4]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[11]-1.0*fIn[7])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[23]-1.0*fIn[17])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[18]+2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS*vBoundary; 
    out[1] += (2.645751311064591*fIn[24]+2.23606797749979*fIn[12]+1.732050807568877*fIn[4]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[11]+fIn[7])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[23]+fIn[17])*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vSer_vF_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[8]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]*intFac; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[5]-1.0*fIn[1])*dS*vBoundary; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[5]+fIn[1])*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vSer_vF_VY_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[20]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[9])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += ((-2.23606797749979*fIn[15])+1.732050807568877*fIn[5]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]-1.0*fIn[7])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[9]+1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (2.23606797749979*fIn[15]+1.732050807568877*fIn[5]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]+fIn[7])*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vSer_vF_VY_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]*intFac; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[19]-2.23606797749979*fIn[9]+1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (2.645751311064591*fIn[27]-2.23606797749979*fIn[15]+1.732050807568877*fIn[5]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]-1.0*fIn[7])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[25]-1.0*fIn[17])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[19]+2.23606797749979*fIn[9]+1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (2.645751311064591*fIn[27]+2.23606797749979*fIn[15]+1.732050807568877*fIn[5]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]+fIn[7])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[25]+fIn[17])*dS*vBoundary; 
 
  }
 
} 
 
