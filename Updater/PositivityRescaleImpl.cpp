#include <PositivityRescaleImpl.h> 
double findMinNodalValue(const double *fIn, int ndim) { 
  double fmin = 0.0; 
  if (ndim == 1) { 
  double fVal[2]; // fVal = array of vales of fIn evaluated at each node 
  fVal[0] = -0.2886751345948129*(1.414213562373095*fIn[1]-2.449489742783178*fIn[0]); 
  fVal[1] = 0.2886751345948129*(1.414213562373095*fIn[1]+2.449489742783178*fIn[0]); 
  fmin = *std::min_element(fVal, fVal+2); 
  } 
  else if(ndim == 2) { 
  double fVal[4]; // fVal = array of vales of fIn evaluated at each node 
  fVal[0] = 0.1666666666666667*(fIn[3]-1.732050807568877*(fIn[2]+fIn[1])+3.0*fIn[0]); 
  fVal[1] = -0.1666666666666667*(fIn[3]+1.732050807568877*fIn[2]-1.732050807568877*fIn[1]-3.0*fIn[0]); 
  fVal[2] = -0.1666666666666667*(fIn[3]-1.732050807568877*fIn[2]+1.732050807568877*fIn[1]-3.0*fIn[0]); 
  fVal[3] = 0.1666666666666667*(fIn[3]+1.732050807568877*(fIn[2]+fIn[1])+3.0*fIn[0]); 
  fmin = *std::min_element(fVal, fVal+4); 
  } 
  else if(ndim == 3) { 
  double fVal[8]; // fVal = array of vales of fIn evaluated at each node 
  fVal[0] = -0.02777777777777778*(2.449489742783178*fIn[7]-4.242640687119286*(fIn[6]+fIn[5]+fIn[4])+7.348469228349534*(fIn[3]+fIn[2]+fIn[1])-12.72792206135786*fIn[0]); 
  fVal[1] = 0.02777777777777778*(2.449489742783178*fIn[7]+4.242640687119286*fIn[6]-4.242640687119286*(fIn[5]+fIn[4])-7.348469228349534*(fIn[3]+fIn[2])+7.348469228349534*fIn[1]+12.72792206135786*fIn[0]); 
  fVal[2] = 0.02777777777777778*(2.449489742783178*fIn[7]-4.242640687119286*fIn[6]+4.242640687119286*fIn[5]-4.242640687119286*fIn[4]-7.348469228349534*fIn[3]+7.348469228349534*fIn[2]-7.348469228349534*fIn[1]+12.72792206135786*fIn[0]); 
  fVal[3] = -0.02777777777777778*(2.449489742783178*fIn[7]+4.242640687119286*(fIn[6]+fIn[5])-4.242640687119286*fIn[4]+7.348469228349534*fIn[3]-7.348469228349534*(fIn[2]+fIn[1])-12.72792206135786*fIn[0]); 
  fVal[4] = 0.02777777777777778*(2.449489742783178*fIn[7]-4.242640687119286*(fIn[6]+fIn[5])+4.242640687119286*fIn[4]+7.348469228349534*fIn[3]-7.348469228349534*(fIn[2]+fIn[1])+12.72792206135786*fIn[0]); 
  fVal[5] = -0.02777777777777778*(2.449489742783178*fIn[7]+4.242640687119286*fIn[6]-4.242640687119286*fIn[5]+4.242640687119286*fIn[4]-7.348469228349534*fIn[3]+7.348469228349534*fIn[2]-7.348469228349534*fIn[1]-12.72792206135786*fIn[0]); 
  fVal[6] = -0.02777777777777778*(2.449489742783178*fIn[7]-4.242640687119286*fIn[6]+4.242640687119286*(fIn[5]+fIn[4])-7.348469228349534*(fIn[3]+fIn[2])+7.348469228349534*fIn[1]-12.72792206135786*fIn[0]); 
  fVal[7] = 0.02777777777777778*(2.449489742783178*fIn[7]+4.242640687119286*(fIn[6]+fIn[5]+fIn[4])+7.348469228349534*(fIn[3]+fIn[2]+fIn[1])+12.72792206135786*fIn[0]); 
  fmin = *std::min_element(fVal, fVal+8); 
  } 
  else if(ndim == 4) { 
  double fVal[16]; // fVal = array of vales of fIn evaluated at each node 
  fVal[0] = 0.02777777777777778*(fIn[15]-1.732050807568877*(fIn[14]+fIn[13]+fIn[12]+fIn[11])+3.0*(fIn[10]+fIn[9]+fIn[8]+fIn[7]+fIn[6]+fIn[5])-5.196152422706631*(fIn[4]+fIn[3]+fIn[2]+fIn[1])+9.0*fIn[0]); 
  fVal[1] = -0.02777777777777778*(fIn[15]+1.732050807568877*fIn[14]-1.732050807568877*(fIn[13]+fIn[12]+fIn[11])-3.0*(fIn[10]+fIn[9])+3.0*fIn[8]-3.0*fIn[7]+3.0*(fIn[6]+fIn[5])+5.196152422706631*(fIn[4]+fIn[3]+fIn[2])-5.196152422706631*fIn[1]-9.0*fIn[0]); 
  fVal[2] = -0.02777777777777778*(fIn[15]-1.732050807568877*fIn[14]+1.732050807568877*fIn[13]-1.732050807568877*(fIn[12]+fIn[11])-3.0*fIn[10]+3.0*fIn[9]-3.0*fIn[8]+3.0*fIn[7]-3.0*fIn[6]+3.0*fIn[5]+5.196152422706631*(fIn[4]+fIn[3])-5.196152422706631*fIn[2]+5.196152422706631*fIn[1]-9.0*fIn[0]); 
  fVal[3] = 0.02777777777777778*(fIn[15]+1.732050807568877*(fIn[14]+fIn[13])-1.732050807568877*(fIn[12]+fIn[11])+3.0*fIn[10]-3.0*(fIn[9]+fIn[8]+fIn[7]+fIn[6])+3.0*fIn[5]-5.196152422706631*(fIn[4]+fIn[3])+5.196152422706631*(fIn[2]+fIn[1])+9.0*fIn[0]); 
  fVal[4] = -0.02777777777777778*(fIn[15]-1.732050807568877*(fIn[14]+fIn[13])+1.732050807568877*fIn[12]-1.732050807568877*fIn[11]+3.0*fIn[10]-3.0*(fIn[9]+fIn[8])+3.0*(fIn[7]+fIn[6])-3.0*fIn[5]+5.196152422706631*fIn[4]-5.196152422706631*fIn[3]+5.196152422706631*(fIn[2]+fIn[1])-9.0*fIn[0]); 
  fVal[5] = 0.02777777777777778*(fIn[15]+1.732050807568877*fIn[14]-1.732050807568877*fIn[13]+1.732050807568877*fIn[12]-1.732050807568877*fIn[11]-3.0*fIn[10]+3.0*fIn[9]-3.0*(fIn[8]+fIn[7])+3.0*fIn[6]-3.0*fIn[5]-5.196152422706631*fIn[4]+5.196152422706631*fIn[3]-5.196152422706631*fIn[2]+5.196152422706631*fIn[1]+9.0*fIn[0]); 
  fVal[6] = 0.02777777777777778*(fIn[15]-1.732050807568877*fIn[14]+1.732050807568877*(fIn[13]+fIn[12])-1.732050807568877*fIn[11]-3.0*(fIn[10]+fIn[9])+3.0*(fIn[8]+fIn[7])-3.0*(fIn[6]+fIn[5])-5.196152422706631*fIn[4]+5.196152422706631*(fIn[3]+fIn[2])-5.196152422706631*fIn[1]+9.0*fIn[0]); 
  fVal[7] = -0.02777777777777778*(fIn[15]+1.732050807568877*(fIn[14]+fIn[13]+fIn[12])-1.732050807568877*fIn[11]+3.0*(fIn[10]+fIn[9]+fIn[8])-3.0*(fIn[7]+fIn[6]+fIn[5])+5.196152422706631*fIn[4]-5.196152422706631*(fIn[3]+fIn[2]+fIn[1])-9.0*fIn[0]); 
  fVal[8] = -0.02777777777777778*(fIn[15]-1.732050807568877*(fIn[14]+fIn[13]+fIn[12])+1.732050807568877*fIn[11]+3.0*(fIn[10]+fIn[9]+fIn[8])-3.0*(fIn[7]+fIn[6]+fIn[5])-5.196152422706631*fIn[4]+5.196152422706631*(fIn[3]+fIn[2]+fIn[1])-9.0*fIn[0]); 
  fVal[9] = 0.02777777777777778*(fIn[15]+1.732050807568877*fIn[14]-1.732050807568877*(fIn[13]+fIn[12])+1.732050807568877*fIn[11]-3.0*(fIn[10]+fIn[9])+3.0*(fIn[8]+fIn[7])-3.0*(fIn[6]+fIn[5])+5.196152422706631*fIn[4]-5.196152422706631*(fIn[3]+fIn[2])+5.196152422706631*fIn[1]+9.0*fIn[0]); 
  fVal[10] = 0.02777777777777778*(fIn[15]-1.732050807568877*fIn[14]+1.732050807568877*fIn[13]-1.732050807568877*fIn[12]+1.732050807568877*fIn[11]-3.0*fIn[10]+3.0*fIn[9]-3.0*(fIn[8]+fIn[7])+3.0*fIn[6]-3.0*fIn[5]+5.196152422706631*fIn[4]-5.196152422706631*fIn[3]+5.196152422706631*fIn[2]-5.196152422706631*fIn[1]+9.0*fIn[0]); 
  fVal[11] = -0.02777777777777778*(fIn[15]+1.732050807568877*(fIn[14]+fIn[13])-1.732050807568877*fIn[12]+1.732050807568877*fIn[11]+3.0*fIn[10]-3.0*(fIn[9]+fIn[8])+3.0*(fIn[7]+fIn[6])-3.0*fIn[5]-5.196152422706631*fIn[4]+5.196152422706631*fIn[3]-5.196152422706631*(fIn[2]+fIn[1])-9.0*fIn[0]); 
  fVal[12] = 0.02777777777777778*(fIn[15]-1.732050807568877*(fIn[14]+fIn[13])+1.732050807568877*(fIn[12]+fIn[11])+3.0*fIn[10]-3.0*(fIn[9]+fIn[8]+fIn[7]+fIn[6])+3.0*fIn[5]+5.196152422706631*(fIn[4]+fIn[3])-5.196152422706631*(fIn[2]+fIn[1])+9.0*fIn[0]); 
  fVal[13] = -0.02777777777777778*(fIn[15]+1.732050807568877*fIn[14]-1.732050807568877*fIn[13]+1.732050807568877*(fIn[12]+fIn[11])-3.0*fIn[10]+3.0*fIn[9]-3.0*fIn[8]+3.0*fIn[7]-3.0*fIn[6]+3.0*fIn[5]-5.196152422706631*(fIn[4]+fIn[3])+5.196152422706631*fIn[2]-5.196152422706631*fIn[1]-9.0*fIn[0]); 
  fVal[14] = -0.02777777777777778*(fIn[15]-1.732050807568877*fIn[14]+1.732050807568877*(fIn[13]+fIn[12]+fIn[11])-3.0*(fIn[10]+fIn[9])+3.0*fIn[8]-3.0*fIn[7]+3.0*(fIn[6]+fIn[5])-5.196152422706631*(fIn[4]+fIn[3]+fIn[2])+5.196152422706631*fIn[1]-9.0*fIn[0]); 
  fVal[15] = 0.02777777777777778*(fIn[15]+1.732050807568877*(fIn[14]+fIn[13]+fIn[12]+fIn[11])+3.0*(fIn[10]+fIn[9]+fIn[8]+fIn[7]+fIn[6]+fIn[5])+5.196152422706631*(fIn[4]+fIn[3]+fIn[2]+fIn[1])+9.0*fIn[0]); 
  fmin = *std::min_element(fVal, fVal+16); 
  } 
  else if(ndim == 5) { 
  double fVal[32]; // fVal = array of vales of fIn evaluated at each node 
  fVal[0] = -0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*(fIn[30]+fIn[29]+fIn[28]+fIn[27]+fIn[26])+7.348469228349534*(fIn[25]+fIn[24]+fIn[23]+fIn[22]+fIn[21]+fIn[20]+fIn[19]+fIn[18]+fIn[17]+fIn[16])-12.72792206135786*(fIn[15]+fIn[14]+fIn[13]+fIn[12]+fIn[11]+fIn[10]+fIn[9]+fIn[8]+fIn[7]+fIn[6])+22.0454076850486*(fIn[5]+fIn[4]+fIn[3]+fIn[2]+fIn[1])-38.18376618407357*fIn[0]); 
  fVal[1] = 0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*fIn[30]-4.242640687119286*(fIn[29]+fIn[28]+fIn[27]+fIn[26])-7.348469228349534*(fIn[25]+fIn[24])+7.348469228349534*fIn[23]-7.348469228349534*fIn[22]+7.348469228349534*(fIn[21]+fIn[20])-7.348469228349534*fIn[19]+7.348469228349534*(fIn[18]+fIn[17]+fIn[16])+12.72792206135786*(fIn[15]+fIn[14]+fIn[13])-12.72792206135786*fIn[12]+12.72792206135786*(fIn[11]+fIn[10])-12.72792206135786*fIn[9]+12.72792206135786*fIn[8]-12.72792206135786*(fIn[7]+fIn[6])-22.0454076850486*(fIn[5]+fIn[4]+fIn[3]+fIn[2])+22.0454076850486*fIn[1]+38.18376618407357*fIn[0]); 
  fVal[2] = 0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*fIn[30]+4.242640687119286*fIn[29]-4.242640687119286*(fIn[28]+fIn[27]+fIn[26])-7.348469228349534*fIn[25]+7.348469228349534*fIn[24]-7.348469228349534*fIn[23]+7.348469228349534*fIn[22]-7.348469228349534*fIn[21]+7.348469228349534*(fIn[20]+fIn[19])-7.348469228349534*fIn[18]+7.348469228349534*(fIn[17]+fIn[16])+12.72792206135786*(fIn[15]+fIn[14])-12.72792206135786*fIn[13]+12.72792206135786*(fIn[12]+fIn[11])-12.72792206135786*fIn[10]+12.72792206135786*fIn[9]-12.72792206135786*fIn[8]+12.72792206135786*fIn[7]-12.72792206135786*fIn[6]-22.0454076850486*(fIn[5]+fIn[4]+fIn[3])+22.0454076850486*fIn[2]-22.0454076850486*fIn[1]+38.18376618407357*fIn[0]); 
  fVal[3] = -0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*(fIn[30]+fIn[29])-4.242640687119286*(fIn[28]+fIn[27]+fIn[26])+7.348469228349534*fIn[25]-7.348469228349534*(fIn[24]+fIn[23]+fIn[22]+fIn[21])+7.348469228349534*fIn[20]-7.348469228349534*(fIn[19]+fIn[18])+7.348469228349534*(fIn[17]+fIn[16])-12.72792206135786*(fIn[15]+fIn[14])+12.72792206135786*(fIn[13]+fIn[12])-12.72792206135786*fIn[11]+12.72792206135786*(fIn[10]+fIn[9]+fIn[8]+fIn[7])-12.72792206135786*fIn[6]+22.0454076850486*(fIn[5]+fIn[4]+fIn[3])-22.0454076850486*(fIn[2]+fIn[1])-38.18376618407357*fIn[0]); 
  fVal[4] = 0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*(fIn[30]+fIn[29])+4.242640687119286*fIn[28]-4.242640687119286*(fIn[27]+fIn[26])+7.348469228349534*fIn[25]-7.348469228349534*(fIn[24]+fIn[23])+7.348469228349534*(fIn[22]+fIn[21])-7.348469228349534*fIn[20]+7.348469228349534*(fIn[19]+fIn[18])-7.348469228349534*fIn[17]+7.348469228349534*fIn[16]+12.72792206135786*fIn[15]-12.72792206135786*fIn[14]+12.72792206135786*(fIn[13]+fIn[12])-12.72792206135786*fIn[11]+12.72792206135786*(fIn[10]+fIn[9])-12.72792206135786*(fIn[8]+fIn[7])+12.72792206135786*fIn[6]-22.0454076850486*(fIn[5]+fIn[4])+22.0454076850486*fIn[3]-22.0454076850486*(fIn[2]+fIn[1])+38.18376618407357*fIn[0]); 
  fVal[5] = -0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*fIn[30]-4.242640687119286*fIn[29]+4.242640687119286*fIn[28]-4.242640687119286*(fIn[27]+fIn[26])-7.348469228349534*fIn[25]+7.348469228349534*fIn[24]-7.348469228349534*(fIn[23]+fIn[22])+7.348469228349534*fIn[21]-7.348469228349534*(fIn[20]+fIn[19])+7.348469228349534*fIn[18]-7.348469228349534*fIn[17]+7.348469228349534*fIn[16]-12.72792206135786*fIn[15]+12.72792206135786*fIn[14]-12.72792206135786*fIn[13]+12.72792206135786*(fIn[12]+fIn[11])-12.72792206135786*fIn[10]+12.72792206135786*(fIn[9]+fIn[8])-12.72792206135786*fIn[7]+12.72792206135786*fIn[6]+22.0454076850486*(fIn[5]+fIn[4])-22.0454076850486*fIn[3]+22.0454076850486*fIn[2]-22.0454076850486*fIn[1]-38.18376618407357*fIn[0]); 
  fVal[6] = -0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*fIn[30]+4.242640687119286*(fIn[29]+fIn[28])-4.242640687119286*(fIn[27]+fIn[26])-7.348469228349534*(fIn[25]+fIn[24])+7.348469228349534*(fIn[23]+fIn[22])-7.348469228349534*(fIn[21]+fIn[20])+7.348469228349534*fIn[19]-7.348469228349534*(fIn[18]+fIn[17])+7.348469228349534*fIn[16]-12.72792206135786*fIn[15]+12.72792206135786*(fIn[14]+fIn[13])-12.72792206135786*fIn[12]+12.72792206135786*(fIn[11]+fIn[10])-12.72792206135786*(fIn[9]+fIn[8])+12.72792206135786*(fIn[7]+fIn[6])+22.0454076850486*(fIn[5]+fIn[4])-22.0454076850486*(fIn[3]+fIn[2])+22.0454076850486*fIn[1]-38.18376618407357*fIn[0]); 
  fVal[7] = 0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*(fIn[30]+fIn[29]+fIn[28])-4.242640687119286*(fIn[27]+fIn[26])+7.348469228349534*(fIn[25]+fIn[24]+fIn[23])-7.348469228349534*(fIn[22]+fIn[21]+fIn[20]+fIn[19]+fIn[18]+fIn[17])+7.348469228349534*fIn[16]+12.72792206135786*fIn[15]-12.72792206135786*(fIn[14]+fIn[13]+fIn[12]+fIn[11]+fIn[10]+fIn[9])+12.72792206135786*(fIn[8]+fIn[7]+fIn[6])-22.0454076850486*(fIn[5]+fIn[4])+22.0454076850486*(fIn[3]+fIn[2]+fIn[1])+38.18376618407357*fIn[0]); 
  fVal[8] = 0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*(fIn[30]+fIn[29]+fIn[28])+4.242640687119286*fIn[27]-4.242640687119286*fIn[26]+7.348469228349534*(fIn[25]+fIn[24]+fIn[23])-7.348469228349534*(fIn[22]+fIn[21]+fIn[20])+7.348469228349534*(fIn[19]+fIn[18]+fIn[17])-7.348469228349534*fIn[16]-12.72792206135786*fIn[15]+12.72792206135786*(fIn[14]+fIn[13]+fIn[12])-12.72792206135786*(fIn[11]+fIn[10]+fIn[9])+12.72792206135786*(fIn[8]+fIn[7]+fIn[6])-22.0454076850486*fIn[5]+22.0454076850486*fIn[4]-22.0454076850486*(fIn[3]+fIn[2]+fIn[1])+38.18376618407357*fIn[0]); 
  fVal[9] = -0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*fIn[30]-4.242640687119286*(fIn[29]+fIn[28])+4.242640687119286*fIn[27]-4.242640687119286*fIn[26]-7.348469228349534*(fIn[25]+fIn[24])+7.348469228349534*(fIn[23]+fIn[22])-7.348469228349534*(fIn[21]+fIn[20]+fIn[19])+7.348469228349534*(fIn[18]+fIn[17])-7.348469228349534*fIn[16]+12.72792206135786*fIn[15]-12.72792206135786*(fIn[14]+fIn[13])+12.72792206135786*(fIn[12]+fIn[11]+fIn[10])-12.72792206135786*(fIn[9]+fIn[8])+12.72792206135786*(fIn[7]+fIn[6])+22.0454076850486*fIn[5]-22.0454076850486*fIn[4]+22.0454076850486*(fIn[3]+fIn[2])-22.0454076850486*fIn[1]-38.18376618407357*fIn[0]); 
  fVal[10] = -0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*fIn[30]+4.242640687119286*fIn[29]-4.242640687119286*fIn[28]+4.242640687119286*fIn[27]-4.242640687119286*fIn[26]-7.348469228349534*fIn[25]+7.348469228349534*fIn[24]-7.348469228349534*(fIn[23]+fIn[22])+7.348469228349534*fIn[21]-7.348469228349534*fIn[20]+7.348469228349534*fIn[19]-7.348469228349534*fIn[18]+7.348469228349534*fIn[17]-7.348469228349534*fIn[16]+12.72792206135786*fIn[15]-12.72792206135786*fIn[14]+12.72792206135786*fIn[13]-12.72792206135786*fIn[12]+12.72792206135786*fIn[11]-12.72792206135786*fIn[10]+12.72792206135786*(fIn[9]+fIn[8])-12.72792206135786*fIn[7]+12.72792206135786*fIn[6]+22.0454076850486*fIn[5]-22.0454076850486*fIn[4]+22.0454076850486*fIn[3]-22.0454076850486*fIn[2]+22.0454076850486*fIn[1]-38.18376618407357*fIn[0]); 
  fVal[11] = 0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*(fIn[30]+fIn[29])-4.242640687119286*fIn[28]+4.242640687119286*fIn[27]-4.242640687119286*fIn[26]+7.348469228349534*fIn[25]-7.348469228349534*(fIn[24]+fIn[23])+7.348469228349534*(fIn[22]+fIn[21])-7.348469228349534*(fIn[20]+fIn[19]+fIn[18])+7.348469228349534*fIn[17]-7.348469228349534*fIn[16]-12.72792206135786*fIn[15]+12.72792206135786*fIn[14]-12.72792206135786*(fIn[13]+fIn[12]+fIn[11])+12.72792206135786*(fIn[10]+fIn[9])-12.72792206135786*(fIn[8]+fIn[7])+12.72792206135786*fIn[6]-22.0454076850486*fIn[5]+22.0454076850486*fIn[4]-22.0454076850486*fIn[3]+22.0454076850486*(fIn[2]+fIn[1])+38.18376618407357*fIn[0]); 
  fVal[12] = -0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*(fIn[30]+fIn[29])+4.242640687119286*(fIn[28]+fIn[27])-4.242640687119286*fIn[26]+7.348469228349534*fIn[25]-7.348469228349534*(fIn[24]+fIn[23]+fIn[22]+fIn[21])+7.348469228349534*(fIn[20]+fIn[19]+fIn[18])-7.348469228349534*(fIn[17]+fIn[16])+12.72792206135786*(fIn[15]+fIn[14])-12.72792206135786*(fIn[13]+fIn[12]+fIn[11])+12.72792206135786*(fIn[10]+fIn[9]+fIn[8]+fIn[7])-12.72792206135786*fIn[6]+22.0454076850486*fIn[5]-22.0454076850486*(fIn[4]+fIn[3])+22.0454076850486*(fIn[2]+fIn[1])-38.18376618407357*fIn[0]); 
  fVal[13] = 0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*fIn[30]-4.242640687119286*fIn[29]+4.242640687119286*(fIn[28]+fIn[27])-4.242640687119286*fIn[26]-7.348469228349534*fIn[25]+7.348469228349534*fIn[24]-7.348469228349534*fIn[23]+7.348469228349534*fIn[22]-7.348469228349534*fIn[21]+7.348469228349534*fIn[20]-7.348469228349534*fIn[19]+7.348469228349534*fIn[18]-7.348469228349534*(fIn[17]+fIn[16])-12.72792206135786*(fIn[15]+fIn[14])+12.72792206135786*fIn[13]-12.72792206135786*fIn[12]+12.72792206135786*fIn[11]-12.72792206135786*fIn[10]+12.72792206135786*fIn[9]-12.72792206135786*fIn[8]+12.72792206135786*fIn[7]-12.72792206135786*fIn[6]-22.0454076850486*fIn[5]+22.0454076850486*(fIn[4]+fIn[3])-22.0454076850486*fIn[2]+22.0454076850486*fIn[1]+38.18376618407357*fIn[0]); 
  fVal[14] = 0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*fIn[30]+4.242640687119286*(fIn[29]+fIn[28]+fIn[27])-4.242640687119286*fIn[26]-7.348469228349534*(fIn[25]+fIn[24])+7.348469228349534*fIn[23]-7.348469228349534*fIn[22]+7.348469228349534*(fIn[21]+fIn[20]+fIn[19])-7.348469228349534*(fIn[18]+fIn[17]+fIn[16])-12.72792206135786*(fIn[15]+fIn[14]+fIn[13])+12.72792206135786*(fIn[12]+fIn[11]+fIn[10])-12.72792206135786*fIn[9]+12.72792206135786*fIn[8]-12.72792206135786*(fIn[7]+fIn[6])-22.0454076850486*fIn[5]+22.0454076850486*(fIn[4]+fIn[3]+fIn[2])-22.0454076850486*fIn[1]+38.18376618407357*fIn[0]); 
  fVal[15] = -0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*(fIn[30]+fIn[29]+fIn[28]+fIn[27])-4.242640687119286*fIn[26]+7.348469228349534*(fIn[25]+fIn[24]+fIn[23]+fIn[22]+fIn[21]+fIn[20])-7.348469228349534*(fIn[19]+fIn[18]+fIn[17]+fIn[16])+12.72792206135786*(fIn[15]+fIn[14]+fIn[13]+fIn[12])-12.72792206135786*(fIn[11]+fIn[10]+fIn[9]+fIn[8]+fIn[7]+fIn[6])+22.0454076850486*fIn[5]-22.0454076850486*(fIn[4]+fIn[3]+fIn[2]+fIn[1])-38.18376618407357*fIn[0]); 
  fVal[16] = 0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*(fIn[30]+fIn[29]+fIn[28]+fIn[27])+4.242640687119286*fIn[26]+7.348469228349534*(fIn[25]+fIn[24]+fIn[23]+fIn[22]+fIn[21]+fIn[20])-7.348469228349534*(fIn[19]+fIn[18]+fIn[17]+fIn[16])-12.72792206135786*(fIn[15]+fIn[14]+fIn[13]+fIn[12])+12.72792206135786*(fIn[11]+fIn[10]+fIn[9]+fIn[8]+fIn[7]+fIn[6])+22.0454076850486*fIn[5]-22.0454076850486*(fIn[4]+fIn[3]+fIn[2]+fIn[1])+38.18376618407357*fIn[0]); 
  fVal[17] = -0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*fIn[30]-4.242640687119286*(fIn[29]+fIn[28]+fIn[27])+4.242640687119286*fIn[26]-7.348469228349534*(fIn[25]+fIn[24])+7.348469228349534*fIn[23]-7.348469228349534*fIn[22]+7.348469228349534*(fIn[21]+fIn[20]+fIn[19])-7.348469228349534*(fIn[18]+fIn[17]+fIn[16])+12.72792206135786*(fIn[15]+fIn[14]+fIn[13])-12.72792206135786*(fIn[12]+fIn[11]+fIn[10])+12.72792206135786*fIn[9]-12.72792206135786*fIn[8]+12.72792206135786*(fIn[7]+fIn[6])-22.0454076850486*fIn[5]+22.0454076850486*(fIn[4]+fIn[3]+fIn[2])-22.0454076850486*fIn[1]-38.18376618407357*fIn[0]); 
  fVal[18] = -0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*fIn[30]+4.242640687119286*fIn[29]-4.242640687119286*(fIn[28]+fIn[27])+4.242640687119286*fIn[26]-7.348469228349534*fIn[25]+7.348469228349534*fIn[24]-7.348469228349534*fIn[23]+7.348469228349534*fIn[22]-7.348469228349534*fIn[21]+7.348469228349534*fIn[20]-7.348469228349534*fIn[19]+7.348469228349534*fIn[18]-7.348469228349534*(fIn[17]+fIn[16])+12.72792206135786*(fIn[15]+fIn[14])-12.72792206135786*fIn[13]+12.72792206135786*fIn[12]-12.72792206135786*fIn[11]+12.72792206135786*fIn[10]-12.72792206135786*fIn[9]+12.72792206135786*fIn[8]-12.72792206135786*fIn[7]+12.72792206135786*fIn[6]-22.0454076850486*fIn[5]+22.0454076850486*(fIn[4]+fIn[3])-22.0454076850486*fIn[2]+22.0454076850486*fIn[1]-38.18376618407357*fIn[0]); 
  fVal[19] = 0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*(fIn[30]+fIn[29])-4.242640687119286*(fIn[28]+fIn[27])+4.242640687119286*fIn[26]+7.348469228349534*fIn[25]-7.348469228349534*(fIn[24]+fIn[23]+fIn[22]+fIn[21])+7.348469228349534*(fIn[20]+fIn[19]+fIn[18])-7.348469228349534*(fIn[17]+fIn[16])-12.72792206135786*(fIn[15]+fIn[14])+12.72792206135786*(fIn[13]+fIn[12]+fIn[11])-12.72792206135786*(fIn[10]+fIn[9]+fIn[8]+fIn[7])+12.72792206135786*fIn[6]+22.0454076850486*fIn[5]-22.0454076850486*(fIn[4]+fIn[3])+22.0454076850486*(fIn[2]+fIn[1])+38.18376618407357*fIn[0]); 
  fVal[20] = -0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*(fIn[30]+fIn[29])+4.242640687119286*fIn[28]-4.242640687119286*fIn[27]+4.242640687119286*fIn[26]+7.348469228349534*fIn[25]-7.348469228349534*(fIn[24]+fIn[23])+7.348469228349534*(fIn[22]+fIn[21])-7.348469228349534*(fIn[20]+fIn[19]+fIn[18])+7.348469228349534*fIn[17]-7.348469228349534*fIn[16]+12.72792206135786*fIn[15]-12.72792206135786*fIn[14]+12.72792206135786*(fIn[13]+fIn[12]+fIn[11])-12.72792206135786*(fIn[10]+fIn[9])+12.72792206135786*(fIn[8]+fIn[7])-12.72792206135786*fIn[6]-22.0454076850486*fIn[5]+22.0454076850486*fIn[4]-22.0454076850486*fIn[3]+22.0454076850486*(fIn[2]+fIn[1])-38.18376618407357*fIn[0]); 
  fVal[21] = 0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*fIn[30]-4.242640687119286*fIn[29]+4.242640687119286*fIn[28]-4.242640687119286*fIn[27]+4.242640687119286*fIn[26]-7.348469228349534*fIn[25]+7.348469228349534*fIn[24]-7.348469228349534*(fIn[23]+fIn[22])+7.348469228349534*fIn[21]-7.348469228349534*fIn[20]+7.348469228349534*fIn[19]-7.348469228349534*fIn[18]+7.348469228349534*fIn[17]-7.348469228349534*fIn[16]-12.72792206135786*fIn[15]+12.72792206135786*fIn[14]-12.72792206135786*fIn[13]+12.72792206135786*fIn[12]-12.72792206135786*fIn[11]+12.72792206135786*fIn[10]-12.72792206135786*(fIn[9]+fIn[8])+12.72792206135786*fIn[7]-12.72792206135786*fIn[6]+22.0454076850486*fIn[5]-22.0454076850486*fIn[4]+22.0454076850486*fIn[3]-22.0454076850486*fIn[2]+22.0454076850486*fIn[1]+38.18376618407357*fIn[0]); 
  fVal[22] = 0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*fIn[30]+4.242640687119286*(fIn[29]+fIn[28])-4.242640687119286*fIn[27]+4.242640687119286*fIn[26]-7.348469228349534*(fIn[25]+fIn[24])+7.348469228349534*(fIn[23]+fIn[22])-7.348469228349534*(fIn[21]+fIn[20]+fIn[19])+7.348469228349534*(fIn[18]+fIn[17])-7.348469228349534*fIn[16]-12.72792206135786*fIn[15]+12.72792206135786*(fIn[14]+fIn[13])-12.72792206135786*(fIn[12]+fIn[11]+fIn[10])+12.72792206135786*(fIn[9]+fIn[8])-12.72792206135786*(fIn[7]+fIn[6])+22.0454076850486*fIn[5]-22.0454076850486*fIn[4]+22.0454076850486*(fIn[3]+fIn[2])-22.0454076850486*fIn[1]+38.18376618407357*fIn[0]); 
  fVal[23] = -0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*(fIn[30]+fIn[29]+fIn[28])-4.242640687119286*fIn[27]+4.242640687119286*fIn[26]+7.348469228349534*(fIn[25]+fIn[24]+fIn[23])-7.348469228349534*(fIn[22]+fIn[21]+fIn[20])+7.348469228349534*(fIn[19]+fIn[18]+fIn[17])-7.348469228349534*fIn[16]+12.72792206135786*fIn[15]-12.72792206135786*(fIn[14]+fIn[13]+fIn[12])+12.72792206135786*(fIn[11]+fIn[10]+fIn[9])-12.72792206135786*(fIn[8]+fIn[7]+fIn[6])-22.0454076850486*fIn[5]+22.0454076850486*fIn[4]-22.0454076850486*(fIn[3]+fIn[2]+fIn[1])-38.18376618407357*fIn[0]); 
  fVal[24] = -0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*(fIn[30]+fIn[29]+fIn[28])+4.242640687119286*(fIn[27]+fIn[26])+7.348469228349534*(fIn[25]+fIn[24]+fIn[23])-7.348469228349534*(fIn[22]+fIn[21]+fIn[20]+fIn[19]+fIn[18]+fIn[17])+7.348469228349534*fIn[16]-12.72792206135786*fIn[15]+12.72792206135786*(fIn[14]+fIn[13]+fIn[12]+fIn[11]+fIn[10]+fIn[9])-12.72792206135786*(fIn[8]+fIn[7]+fIn[6])-22.0454076850486*(fIn[5]+fIn[4])+22.0454076850486*(fIn[3]+fIn[2]+fIn[1])-38.18376618407357*fIn[0]); 
  fVal[25] = 0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*fIn[30]-4.242640687119286*(fIn[29]+fIn[28])+4.242640687119286*(fIn[27]+fIn[26])-7.348469228349534*(fIn[25]+fIn[24])+7.348469228349534*(fIn[23]+fIn[22])-7.348469228349534*(fIn[21]+fIn[20])+7.348469228349534*fIn[19]-7.348469228349534*(fIn[18]+fIn[17])+7.348469228349534*fIn[16]+12.72792206135786*fIn[15]-12.72792206135786*(fIn[14]+fIn[13])+12.72792206135786*fIn[12]-12.72792206135786*(fIn[11]+fIn[10])+12.72792206135786*(fIn[9]+fIn[8])-12.72792206135786*(fIn[7]+fIn[6])+22.0454076850486*(fIn[5]+fIn[4])-22.0454076850486*(fIn[3]+fIn[2])+22.0454076850486*fIn[1]+38.18376618407357*fIn[0]); 
  fVal[26] = 0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*fIn[30]+4.242640687119286*fIn[29]-4.242640687119286*fIn[28]+4.242640687119286*(fIn[27]+fIn[26])-7.348469228349534*fIn[25]+7.348469228349534*fIn[24]-7.348469228349534*(fIn[23]+fIn[22])+7.348469228349534*fIn[21]-7.348469228349534*(fIn[20]+fIn[19])+7.348469228349534*fIn[18]-7.348469228349534*fIn[17]+7.348469228349534*fIn[16]+12.72792206135786*fIn[15]-12.72792206135786*fIn[14]+12.72792206135786*fIn[13]-12.72792206135786*(fIn[12]+fIn[11])+12.72792206135786*fIn[10]-12.72792206135786*(fIn[9]+fIn[8])+12.72792206135786*fIn[7]-12.72792206135786*fIn[6]+22.0454076850486*(fIn[5]+fIn[4])-22.0454076850486*fIn[3]+22.0454076850486*fIn[2]-22.0454076850486*fIn[1]+38.18376618407357*fIn[0]); 
  fVal[27] = -0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*(fIn[30]+fIn[29])-4.242640687119286*fIn[28]+4.242640687119286*(fIn[27]+fIn[26])+7.348469228349534*fIn[25]-7.348469228349534*(fIn[24]+fIn[23])+7.348469228349534*(fIn[22]+fIn[21])-7.348469228349534*fIn[20]+7.348469228349534*(fIn[19]+fIn[18])-7.348469228349534*fIn[17]+7.348469228349534*fIn[16]-12.72792206135786*fIn[15]+12.72792206135786*fIn[14]-12.72792206135786*(fIn[13]+fIn[12])+12.72792206135786*fIn[11]-12.72792206135786*(fIn[10]+fIn[9])+12.72792206135786*(fIn[8]+fIn[7])-12.72792206135786*fIn[6]-22.0454076850486*(fIn[5]+fIn[4])+22.0454076850486*fIn[3]-22.0454076850486*(fIn[2]+fIn[1])-38.18376618407357*fIn[0]); 
  fVal[28] = 0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*(fIn[30]+fIn[29])+4.242640687119286*(fIn[28]+fIn[27]+fIn[26])+7.348469228349534*fIn[25]-7.348469228349534*(fIn[24]+fIn[23]+fIn[22]+fIn[21])+7.348469228349534*fIn[20]-7.348469228349534*(fIn[19]+fIn[18])+7.348469228349534*(fIn[17]+fIn[16])+12.72792206135786*(fIn[15]+fIn[14])-12.72792206135786*(fIn[13]+fIn[12])+12.72792206135786*fIn[11]-12.72792206135786*(fIn[10]+fIn[9]+fIn[8]+fIn[7])+12.72792206135786*fIn[6]+22.0454076850486*(fIn[5]+fIn[4]+fIn[3])-22.0454076850486*(fIn[2]+fIn[1])+38.18376618407357*fIn[0]); 
  fVal[29] = -0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*fIn[30]-4.242640687119286*fIn[29]+4.242640687119286*(fIn[28]+fIn[27]+fIn[26])-7.348469228349534*fIn[25]+7.348469228349534*fIn[24]-7.348469228349534*fIn[23]+7.348469228349534*fIn[22]-7.348469228349534*fIn[21]+7.348469228349534*(fIn[20]+fIn[19])-7.348469228349534*fIn[18]+7.348469228349534*(fIn[17]+fIn[16])-12.72792206135786*(fIn[15]+fIn[14])+12.72792206135786*fIn[13]-12.72792206135786*(fIn[12]+fIn[11])+12.72792206135786*fIn[10]-12.72792206135786*fIn[9]+12.72792206135786*fIn[8]-12.72792206135786*fIn[7]+12.72792206135786*fIn[6]-22.0454076850486*(fIn[5]+fIn[4]+fIn[3])+22.0454076850486*fIn[2]-22.0454076850486*fIn[1]-38.18376618407357*fIn[0]); 
  fVal[30] = -0.004629629629629629*(2.449489742783178*fIn[31]-4.242640687119286*fIn[30]+4.242640687119286*(fIn[29]+fIn[28]+fIn[27]+fIn[26])-7.348469228349534*(fIn[25]+fIn[24])+7.348469228349534*fIn[23]-7.348469228349534*fIn[22]+7.348469228349534*(fIn[21]+fIn[20])-7.348469228349534*fIn[19]+7.348469228349534*(fIn[18]+fIn[17]+fIn[16])-12.72792206135786*(fIn[15]+fIn[14]+fIn[13])+12.72792206135786*fIn[12]-12.72792206135786*(fIn[11]+fIn[10])+12.72792206135786*fIn[9]-12.72792206135786*fIn[8]+12.72792206135786*(fIn[7]+fIn[6])-22.0454076850486*(fIn[5]+fIn[4]+fIn[3]+fIn[2])+22.0454076850486*fIn[1]-38.18376618407357*fIn[0]); 
  fVal[31] = 0.004629629629629629*(2.449489742783178*fIn[31]+4.242640687119286*(fIn[30]+fIn[29]+fIn[28]+fIn[27]+fIn[26])+7.348469228349534*(fIn[25]+fIn[24]+fIn[23]+fIn[22]+fIn[21]+fIn[20]+fIn[19]+fIn[18]+fIn[17]+fIn[16])+12.72792206135786*(fIn[15]+fIn[14]+fIn[13]+fIn[12]+fIn[11]+fIn[10]+fIn[9]+fIn[8]+fIn[7]+fIn[6])+22.0454076850486*(fIn[5]+fIn[4]+fIn[3]+fIn[2]+fIn[1])+38.18376618407357*fIn[0]); 
  fmin = *std::min_element(fVal, fVal+32); 
  } 
  return fmin; 
}

// check positivity of cell average and control nodes
bool check(const double *fIn, int ndim, int numBasis, int *idx, double tCurr, int rkIdx)
{
  double f0 = fIn[0]*std::pow(0.7071067811865475,ndim);
  bool status = true;

  if(f0 < 0.) {
     if(ndim == 1) {
       printf("WARNING: negative cell avg %e in cell %2d, tCurr = %e\n", f0, idx[0], tCurr);
     } else if( ndim == 2) {
       printf("WARNING: negative cell avg %e in cell %2d %2d, tCurr = %e, rkIdx = %d\n", f0, idx[0], idx[1], tCurr, rkIdx);
     } else if( ndim == 3) {
       printf("WARNING: negative cell avg %e in cell %2d %2d %2d, tCurr = %e\n", f0, idx[0], idx[1], idx[2], tCurr);
     } else if( ndim == 4) {
       printf("WARNING: negative cell avg %e in cell %2d %2d %2d %2d, tCurr = %e\n", f0, idx[0], idx[1], idx[2], idx[3], tCurr);
     } else if( ndim == 5) {
       printf("WARNING: negative cell avg %e in cell %2d %2d %2d %2d %2d, tCurr = %e, rkIdx = %d\n", f0, idx[0], idx[1], idx[2], idx[3], idx[4], tCurr, rkIdx);
     }
     status = false;
  }

  double fmin = findMinNodalValue(fIn, ndim);
  if (fmin < 0. && status) {
     if(ndim == 1) {
       printf("warning: negative control node %e in cell %2d, tCurr = %e \n", fmin, idx[0], tCurr);
     } else if(ndim == 2) {
       printf("warning: negative control node %e in cell %2d %2d, tCurr = %e\n", fmin, idx[0], idx[1], tCurr);
     } else if(ndim == 3) {
       printf("warning: negative control node %e in cell %2d %2d %2d, tCurr = %e\n", fmin, idx[0], idx[1], idx[2], tCurr);
     } else if(ndim == 4) {
       printf("warning: negative control node %e in cell %2d %2d %2d %2d, tCurr = %e\n", fmin, idx[0], idx[1], idx[2], idx[3], tCurr);
     } else if(ndim == 5) {
       printf("warning: negative control node %e in cell %2d %2d %2d %2d %2d, tCurr = %e\n", fmin, idx[0], idx[1], idx[2], idx[3], idx[4], tCurr);
     }
  }
  return status;
}


double rescale(const double *fIn, double *fOut, int ndim, int numBasis, int *idx, double tCurr)
{
  double f0 = fIn[0]*std::pow(0.7071067811865475,ndim);
  if (f0 < 0.) return 0.;

  //if(f0 < 0.) {
  //   if(ndim == 1) {
  //     printf("WARNING: negative cell avg %e in cell %2d, tCurr = %e\n", f0, idx[0], tCurr);
  //   } else if( ndim == 2) {
  //     printf("WARNING: negative cell avg %e in cell %2d %2d, tCurr = %e\n", f0, idx[0], idx[1], tCurr);
  //   } else if( ndim == 3) {
  //     printf("WARNING: negative cell avg %e in cell %2d %2d %2d, tCurr = %e\n", f0, idx[0], idx[1], idx[2], tCurr);
  //   } else if( ndim == 4) {
  //     printf("WARNING: negative cell avg %e in cell %2d %2d %2d %2d, tCurr = %e\n", f0, idx[0], idx[1], idx[2], idx[3], tCurr);
  //   } else if( ndim == 5) {
  //     printf("WARNING: negative cell avg %e in cell %2d %2d %2d %2d %2d, tCurr = %e\n", f0, idx[0], idx[1], idx[2], idx[3], idx[4], tCurr);
  //   }
  //}
  double fmin = findMinNodalValue(fIn, ndim);
  double del2Change = 0.;
  int j = 0;

  while (fmin < 0.) {
     double theta = std::min(1.0, f0/(f0 - fmin + EPSILON));

     // modify moments. note no change to cell average
     fOut[0] = fIn[0]; 
     double del2 = 0.0;
     for(int i=1; i<numBasis; i++) {
       if(theta < 1 && j==0) {
          del2 += fIn[i]*fIn[i];
          fOut[i] = theta*fIn[i];
       } else {
          fOut[i] = theta*fOut[i];
       }
     }
     del2Change += del2*(1-theta)*(1-theta);
     
     fmin = findMinNodalValue(fOut, ndim);
     j++;
  }
  
  return del2Change;
}
