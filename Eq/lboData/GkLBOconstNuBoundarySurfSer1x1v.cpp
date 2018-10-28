#include <GkLBOModDecl.h> 
double GkLBOconstNuBoundarySurf1x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dxv[2]:    Cell spacing.
  // idx[2]:    current grid index.
  // nu:        constant collisionality. 
  // vMuMidMax: maximum midpoint value of v-u. 
  // u[1*2]:    bulk velocity (in 1 directions). 
  // vtSq[2]:   thermal speed squared. 
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  double alpha[4]; 
  alpha[0] = 1.414213562373095*vtSq[0]; 
  alpha[1] = 1.414213562373095*vtSq[1]; 

  if (idxr[1] == 1) {

    outr[2] += (0.4330127018922193*(alpha[1]*fr[1]+alpha[0]*fr[0])-0.75*(alpha[1]*fr[3]+alpha[0]*fr[2]))*rdvSq4nuR; 
    outr[3] += (0.4330127018922193*(alpha[0]*fr[1]+fr[0]*alpha[1])-0.75*(alpha[0]*fr[3]+alpha[1]*fr[2]))*rdvSq4nuR; 

  } else {

    outl[2] += ((-0.75*(alpha[1]*fl[3]+alpha[0]*fl[2]))-0.4330127018922193*(alpha[1]*fl[1]+alpha[0]*fl[0]))*rdvSq4nuL; 
    outl[3] += ((-0.75*(alpha[0]*fl[3]+alpha[1]*fl[2]))-0.4330127018922193*(alpha[0]*fl[1]+fl[0]*alpha[1]))*rdvSq4nuL; 

  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurf1x1vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dxv[2]:    Cell spacing.
  // idx[2]:    current grid index.
  // nu:        constant collisionality. 
  // vMuMidMax: maximum midpoint value of v-u. 
  // u[1*3]:    bulk velocity (in 1 directions). 
  // vtSq[3]:   thermal speed squared. 
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  double alpha[8]; 
  alpha[0] = 1.414213562373095*vtSq[0]; 
  alpha[1] = 1.414213562373095*vtSq[1]; 
  alpha[4] = 1.414213562373095*vtSq[2]; 

  if (idxr[1] == 1) {

    outr[2] += (0.9682458365518543*alpha[1]*fr[7]-0.75*alpha[4]*fr[6]+0.9682458365518543*alpha[0]*fr[5]+0.4330127018922193*alpha[4]*fr[4]-0.75*(alpha[1]*fr[3]+alpha[0]*fr[2])+0.4330127018922193*(alpha[1]*fr[1]+alpha[0]*fr[0]))*rdvSq4nuR; 
    outr[3] += ((0.8660254037844387*alpha[4]+0.9682458365518543*alpha[0])*fr[7]+alpha[1]*((-0.6708203932499369*fr[6])+0.9682458365518543*fr[5]+0.3872983346207416*fr[4])+(0.3872983346207416*fr[1]-0.6708203932499369*fr[3])*alpha[4]-0.75*(alpha[0]*fr[3]+alpha[1]*fr[2])+0.4330127018922193*(alpha[0]*fr[1]+fr[0]*alpha[1]))*rdvSq4nuR; 
    outr[5] += ((-3.75*alpha[1]*fr[7])+2.904737509655563*alpha[4]*fr[6]-3.75*alpha[0]*fr[5]-1.677050983124842*alpha[4]*fr[4]+2.904737509655563*(alpha[1]*fr[3]+alpha[0]*fr[2])-1.677050983124842*(alpha[1]*fr[1]+alpha[0]*fr[0]))*rdvSq4nuR; 
    outr[6] += (0.8660254037844386*alpha[1]*fr[7]+((-0.479157423749955*alpha[4])-0.75*alpha[0])*fr[6]+0.9682458365518543*alpha[4]*fr[5]+(0.276641667586244*alpha[4]+0.4330127018922194*alpha[0])*fr[4]+(0.4330127018922194*fr[0]-0.75*fr[2])*alpha[4]+alpha[1]*(0.3872983346207417*fr[1]-0.6708203932499369*fr[3]))*rdvSq4nuR; 
    outr[7] += (((-3.354101966249685*alpha[4])-3.75*alpha[0])*fr[7]+alpha[1]*(2.598076211353316*fr[6]-3.75*fr[5]-1.5*fr[4])+(2.598076211353316*fr[3]-1.5*fr[1])*alpha[4]+2.904737509655563*(alpha[0]*fr[3]+alpha[1]*fr[2])-1.677050983124842*(alpha[0]*fr[1]+fr[0]*alpha[1]))*rdvSq4nuR; 

  } else {

    outl[2] += ((-0.9682458365518543*alpha[1]*fl[7])-0.75*alpha[4]*fl[6]-0.9682458365518543*alpha[0]*fl[5]-0.4330127018922193*alpha[4]*fl[4]-0.75*(alpha[1]*fl[3]+alpha[0]*fl[2])-0.4330127018922193*(alpha[1]*fl[1]+alpha[0]*fl[0]))*rdvSq4nuL; 
    outl[3] += (((-0.8660254037844387*alpha[4])-0.9682458365518543*alpha[0])*fl[7]+alpha[1]*((-0.6708203932499369*fl[6])-0.9682458365518543*fl[5]-0.3872983346207416*fl[4])+((-0.6708203932499369*fl[3])-0.3872983346207416*fl[1])*alpha[4]-0.75*(alpha[0]*fl[3]+alpha[1]*fl[2])-0.4330127018922193*(alpha[0]*fl[1]+fl[0]*alpha[1]))*rdvSq4nuL; 
    outl[5] += ((-3.75*alpha[1]*fl[7])-2.904737509655563*alpha[4]*fl[6]-3.75*alpha[0]*fl[5]-1.677050983124842*alpha[4]*fl[4]-2.904737509655563*(alpha[1]*fl[3]+alpha[0]*fl[2])-1.677050983124842*(alpha[1]*fl[1]+alpha[0]*fl[0]))*rdvSq4nuL; 
    outl[6] += ((-0.8660254037844386*alpha[1]*fl[7])+((-0.479157423749955*alpha[4])-0.75*alpha[0])*fl[6]-0.9682458365518543*alpha[4]*fl[5]+((-0.276641667586244*alpha[4])-0.4330127018922194*alpha[0])*fl[4]+((-0.75*fl[2])-0.4330127018922194*fl[0])*alpha[4]+alpha[1]*((-0.6708203932499369*fl[3])-0.3872983346207417*fl[1]))*rdvSq4nuL; 
    outl[7] += (((-3.354101966249685*alpha[4])-3.75*alpha[0])*fl[7]+alpha[1]*((-2.598076211353316*fl[6])-3.75*fl[5]-1.5*fl[4])+((-2.598076211353316*fl[3])-1.5*fl[1])*alpha[4]-2.904737509655563*(alpha[0]*fl[3]+alpha[1]*fl[2])-1.677050983124842*(alpha[0]*fl[1]+fl[0]*alpha[1]))*rdvSq4nuL; 

  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurf1x1vSer_Vpar_P3(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dxv[2]:    Cell spacing.
  // idx[2]:    current grid index.
  // nu:        constant collisionality. 
  // vMuMidMax: maximum midpoint value of v-u. 
  // u[1*4]:    bulk velocity (in 1 directions). 
  // vtSq[4]:   thermal speed squared. 
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  double alpha[12]; 
  alpha[0] = 1.414213562373095*vtSq[0]; 
  alpha[1] = 1.414213562373095*vtSq[1]; 
  alpha[4] = 1.414213562373095*vtSq[2]; 
  alpha[8] = 1.414213562373095*vtSq[3]; 

  if (idxr[1] == 1) {

    outr[2] += ((-1.14564392373896*alpha[1]*fr[11])-0.7499999999999999*alpha[8]*fr[10]-1.14564392373896*alpha[0]*fr[9]+0.4330127018922193*alpha[8]*fr[8]+0.9682458365518543*alpha[1]*fr[7]-0.75*alpha[4]*fr[6]+0.9682458365518543*alpha[0]*fr[5]+0.4330127018922193*alpha[4]*fr[4]-0.75*(alpha[1]*fr[3]+alpha[0]*fr[2])+0.4330127018922193*(alpha[1]*fr[1]+alpha[0]*fr[0]))*rdvSq4nuR; 
    outr[3] += (((-1.02469507659596*alpha[4])-1.14564392373896*alpha[0])*fr[11]-0.6587325492402599*alpha[4]*fr[10]-1.14564392373896*alpha[1]*fr[9]+0.3803194146278324*alpha[4]*fr[8]+(0.3803194146278324*fr[4]-0.6587325492402599*fr[6])*alpha[8]+(0.8660254037844387*alpha[4]+0.9682458365518543*alpha[0])*fr[7]+alpha[1]*((-0.6708203932499369*fr[6])+0.9682458365518543*fr[5]+0.3872983346207416*fr[4])+(0.3872983346207416*fr[1]-0.6708203932499369*fr[3])*alpha[4]-0.75*(alpha[0]*fr[3]+alpha[1]*fr[2])+0.4330127018922193*(alpha[0]*fr[1]+fr[0]*alpha[1]))*rdvSq4nuR; 
    outr[5] += (4.437059837324712*alpha[1]*fr[11]+2.904737509655563*alpha[8]*fr[10]+4.437059837324712*alpha[0]*fr[9]-1.677050983124842*alpha[8]*fr[8]-3.75*alpha[1]*fr[7]+2.904737509655563*alpha[4]*fr[6]-3.75*alpha[0]*fr[5]-1.677050983124842*alpha[4]*fr[4]+2.904737509655563*(alpha[1]*fr[3]+alpha[0]*fr[2])-1.677050983124842*(alpha[1]*fr[1]+alpha[0]*fr[0]))*rdvSq4nuR; 
    outr[6] += (((-1.006230589874905*alpha[8])-1.02469507659596*alpha[1])*fr[11]+((-0.4472135954999578*alpha[8])-0.6587325492402597*alpha[1])*fr[10]-1.14564392373896*alpha[4]*fr[9]+(0.2581988897471611*alpha[8]+0.3803194146278324*alpha[1])*fr[8]+(0.8504200642707612*fr[7]-0.6587325492402599*fr[3]+0.3803194146278324*fr[1])*alpha[8]+0.8660254037844386*alpha[1]*fr[7]+((-0.479157423749955*alpha[4])-0.75*alpha[0])*fr[6]+0.9682458365518543*alpha[4]*fr[5]+(0.276641667586244*alpha[4]+0.4330127018922194*alpha[0])*fr[4]+(0.4330127018922194*fr[0]-0.75*fr[2])*alpha[4]+alpha[1]*(0.3872983346207417*fr[1]-0.6708203932499369*fr[3]))*rdvSq4nuR; 
    outr[7] += ((3.968626966596886*alpha[4]+4.437059837324712*alpha[0])*fr[11]+2.551260192812284*alpha[4]*fr[10]+4.437059837324712*alpha[1]*fr[9]-1.472970759092948*alpha[4]*fr[8]+(2.551260192812284*fr[6]-1.472970759092948*fr[4])*alpha[8]+((-3.354101966249685*alpha[4])-3.75*alpha[0])*fr[7]+alpha[1]*(2.598076211353316*fr[6]-3.75*fr[5]-1.5*fr[4])+(2.598076211353316*fr[3]-1.5*fr[1])*alpha[4]+2.904737509655563*(alpha[0]*fr[3]+alpha[1]*fr[2])-1.677050983124842*(alpha[0]*fr[1]+fr[0]*alpha[1]))*rdvSq4nuR; 
    outr[9] += ((-10.5*alpha[1]*fr[11])-6.87386354243376*alpha[8]*fr[10]-10.5*alpha[0]*fr[9]+3.968626966596886*alpha[8]*fr[8]+8.874119674649425*alpha[1]*fr[7]-6.87386354243376*alpha[4]*fr[6]+8.874119674649425*alpha[0]*fr[5]+3.968626966596886*alpha[4]*fr[4]-6.873863542433759*(alpha[1]*fr[3]+alpha[0]*fr[2])+3.968626966596886*(alpha[1]*fr[1]+alpha[0]*fr[0]))*rdvSq4nuR; 
    outr[10] += ((-1.006230589874905*alpha[4]*fr[11])+((-0.4472135954999579*alpha[4])-0.75*alpha[0])*fr[10]-1.14564392373896*alpha[8]*fr[9]+(0.258198889747161*alpha[4]+0.4330127018922193*alpha[0])*fr[8]+((-0.4472135954999578*fr[6])+0.9682458365518541*fr[5]+0.258198889747161*fr[4]-0.7499999999999999*fr[2]+0.4330127018922193*fr[0])*alpha[8]+0.8504200642707612*alpha[4]*fr[7]+alpha[1]*(0.3803194146278324*fr[4]-0.6587325492402597*fr[6])+(0.3803194146278324*fr[1]-0.6587325492402599*fr[3])*alpha[4])*rdvSq4nuR; 
    outr[11] += (((-9.391485505499116*alpha[4])-10.5*alpha[0])*fr[11]-6.037383539249432*alpha[4]*fr[10]-10.5*alpha[1]*fr[9]+3.485685011586675*alpha[4]*fr[8]+(3.485685011586675*fr[4]-6.037383539249432*fr[6])*alpha[8]+(7.937253933193772*alpha[4]+8.874119674649425*alpha[0])*fr[7]+alpha[1]*((-6.148170459575759*fr[6])+8.874119674649425*fr[5]+3.549647869859769*fr[4])+(3.549647869859769*fr[1]-6.148170459575758*fr[3])*alpha[4]-6.873863542433758*(alpha[0]*fr[3]+alpha[1]*fr[2])+3.968626966596886*(alpha[0]*fr[1]+fr[0]*alpha[1]))*rdvSq4nuR; 

  } else {

    outl[2] += ((-1.14564392373896*alpha[1]*fl[11])-0.7499999999999999*alpha[8]*fl[10]-1.14564392373896*alpha[0]*fl[9]-0.4330127018922193*alpha[8]*fl[8]-0.9682458365518543*alpha[1]*fl[7]-0.75*alpha[4]*fl[6]-0.9682458365518543*alpha[0]*fl[5]-0.4330127018922193*alpha[4]*fl[4]-0.75*(alpha[1]*fl[3]+alpha[0]*fl[2])-0.4330127018922193*(alpha[1]*fl[1]+alpha[0]*fl[0]))*rdvSq4nuL; 
    outl[3] += (((-1.02469507659596*alpha[4])-1.14564392373896*alpha[0])*fl[11]-0.6587325492402599*alpha[4]*fl[10]-1.14564392373896*alpha[1]*fl[9]-0.3803194146278324*alpha[4]*fl[8]+((-0.6587325492402599*fl[6])-0.3803194146278324*fl[4])*alpha[8]+((-0.8660254037844387*alpha[4])-0.9682458365518543*alpha[0])*fl[7]+alpha[1]*((-0.6708203932499369*fl[6])-0.9682458365518543*fl[5]-0.3872983346207416*fl[4])+((-0.6708203932499369*fl[3])-0.3872983346207416*fl[1])*alpha[4]-0.75*(alpha[0]*fl[3]+alpha[1]*fl[2])-0.4330127018922193*(alpha[0]*fl[1]+fl[0]*alpha[1]))*rdvSq4nuL; 
    outl[5] += ((-4.437059837324712*alpha[1]*fl[11])-2.904737509655563*alpha[8]*fl[10]-4.437059837324712*alpha[0]*fl[9]-1.677050983124842*alpha[8]*fl[8]-3.75*alpha[1]*fl[7]-2.904737509655563*alpha[4]*fl[6]-3.75*alpha[0]*fl[5]-1.677050983124842*alpha[4]*fl[4]-2.904737509655563*(alpha[1]*fl[3]+alpha[0]*fl[2])-1.677050983124842*(alpha[1]*fl[1]+alpha[0]*fl[0]))*rdvSq4nuL; 
    outl[6] += (((-1.006230589874905*alpha[8])-1.02469507659596*alpha[1])*fl[11]+((-0.4472135954999578*alpha[8])-0.6587325492402597*alpha[1])*fl[10]-1.14564392373896*alpha[4]*fl[9]+((-0.2581988897471611*alpha[8])-0.3803194146278324*alpha[1])*fl[8]+((-0.8504200642707612*fl[7])-0.6587325492402599*fl[3]-0.3803194146278324*fl[1])*alpha[8]-0.8660254037844386*alpha[1]*fl[7]+((-0.479157423749955*alpha[4])-0.75*alpha[0])*fl[6]-0.9682458365518543*alpha[4]*fl[5]+((-0.276641667586244*alpha[4])-0.4330127018922194*alpha[0])*fl[4]+((-0.75*fl[2])-0.4330127018922194*fl[0])*alpha[4]+alpha[1]*((-0.6708203932499369*fl[3])-0.3872983346207417*fl[1]))*rdvSq4nuL; 
    outl[7] += (((-3.968626966596886*alpha[4])-4.437059837324712*alpha[0])*fl[11]-2.551260192812284*alpha[4]*fl[10]-4.437059837324712*alpha[1]*fl[9]-1.472970759092948*alpha[4]*fl[8]+((-2.551260192812284*fl[6])-1.472970759092948*fl[4])*alpha[8]+((-3.354101966249685*alpha[4])-3.75*alpha[0])*fl[7]+alpha[1]*((-2.598076211353316*fl[6])-3.75*fl[5]-1.5*fl[4])+((-2.598076211353316*fl[3])-1.5*fl[1])*alpha[4]-2.904737509655563*(alpha[0]*fl[3]+alpha[1]*fl[2])-1.677050983124842*(alpha[0]*fl[1]+fl[0]*alpha[1]))*rdvSq4nuL; 
    outl[9] += ((-10.5*alpha[1]*fl[11])-6.87386354243376*alpha[8]*fl[10]-10.5*alpha[0]*fl[9]-3.968626966596886*alpha[8]*fl[8]-8.874119674649425*alpha[1]*fl[7]-6.87386354243376*alpha[4]*fl[6]-8.874119674649425*alpha[0]*fl[5]-3.968626966596886*alpha[4]*fl[4]-6.873863542433759*(alpha[1]*fl[3]+alpha[0]*fl[2])-3.968626966596886*(alpha[1]*fl[1]+alpha[0]*fl[0]))*rdvSq4nuL; 
    outl[10] += ((-1.006230589874905*alpha[4]*fl[11])+((-0.4472135954999579*alpha[4])-0.75*alpha[0])*fl[10]-1.14564392373896*alpha[8]*fl[9]+((-0.258198889747161*alpha[4])-0.4330127018922193*alpha[0])*fl[8]+((-0.4472135954999578*fl[6])-0.9682458365518541*fl[5]-0.258198889747161*fl[4]-0.7499999999999999*fl[2]-0.4330127018922193*fl[0])*alpha[8]-0.8504200642707612*alpha[4]*fl[7]+alpha[1]*((-0.6587325492402597*fl[6])-0.3803194146278324*fl[4])+((-0.6587325492402599*fl[3])-0.3803194146278324*fl[1])*alpha[4])*rdvSq4nuL; 
    outl[11] += (((-9.391485505499116*alpha[4])-10.5*alpha[0])*fl[11]-6.037383539249432*alpha[4]*fl[10]-10.5*alpha[1]*fl[9]-3.485685011586675*alpha[4]*fl[8]+((-6.037383539249432*fl[6])-3.485685011586675*fl[4])*alpha[8]+((-7.937253933193772*alpha[4])-8.874119674649425*alpha[0])*fl[7]+alpha[1]*((-6.148170459575759*fl[6])-8.874119674649425*fl[5]-3.549647869859769*fl[4])+((-6.148170459575758*fl[3])-3.549647869859769*fl[1])*alpha[4]-6.873863542433758*(alpha[0]*fl[3]+alpha[1]*fl[2])-3.968626966596886*(alpha[0]*fl[1]+fl[0]*alpha[1]))*rdvSq4nuL; 

  }
  return 0.0; 
} 
