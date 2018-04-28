#include <GyrokineticModDecl.h> 
double calcSheathDeltaPhi1xSer_P1(const double *phi, const double *phiWall, const double zVal) 
{ 
  return 1.224744871391589*phi[1]*zVal-0.7071067811865475*phiWall[0]+0.7071067811865475*phi[0]; 
}

void calcSheathPartialReflection1x1vSer_P1(const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat) 
{ 
  double zVal2 = zVal*zVal; 
  double wv2 = wv*wv; 
  double wv3 = wv2*wv; 
  double wv4 = wv3*wv; 
  double wv5 = wv4*wv; 
  double dv2 = dv*dv; 
  double dv3 = dv2*dv; 
  double dv4 = dv3*dv; 
  double dv5 = dv4*dv; 
  double vcut2 = vcut*vcut; 
  double vcut3 = vcut2*vcut; 
  double vcut4 = vcut3*vcut; 
  double vcut5 = vcut4*vcut; 
  double denom = 8*wv3+((-24*vcut)-12*dv)*wv2+(24*vcut2+24*dv*vcut+6*dv2)*wv-8*vcut3-12*dv*vcut2-6*dv2*vcut-dv3; 
  fhat[0] = -(1.0*((55.42562584220407*f[1]*dv*wv2+((24.0*f[3]+27.71281292110204*f[1])*dv2-110.8512516844081*f[1]*dv*vcut)*wv+55.42562584220407*f[1]*dv*vcut2+((-24.0*f[3])-27.71281292110204*f[1])*dv2*vcut+(12.0*f[3]+13.85640646055102*f[1])*dv3)*zVal+32.0*f[0]*dv*wv2+((13.85640646055102*f[2]+16.0*f[0])*dv2-64.0*f[0]*dv*vcut)*wv+32.0*f[0]*dv*vcut2+((-13.85640646055102*f[2])-16.0*f[0])*dv2*vcut+(6.928203230275509*f[2]+8.0*f[0])*dv3))/denom; 
  fhat[1] = 0.0; 
  fhat[2] = -(1.0*((24.0*f[1]*dv2*wv-24.0*f[1]*dv2*vcut+(13.85640646055102*f[3]+12.0*f[1])*dv3)*zVal+13.85640646055102*f[0]*dv2*wv-13.85640646055102*f[0]*dv2*vcut+(8.0*f[2]+6.928203230275509*f[0])*dv3))/denom; 
  fhat[3] = 0.0; 
} 

double calcSheathDeltaPhi1xSer_P2(const double *phi, const double *phiWall, const double zVal) 
{ 
  double zVal2 = zVal*zVal; 
  return 2.371708245126284*phi[2]*(zVal2-0.3333333333333333)+1.224744871391589*phi[1]*zVal+0.7905694150420947*phiWall[2]-0.7071067811865475*phiWall[0]+0.7071067811865475*phi[0]; 
}

void calcSheathPartialReflection1x1vSer_P2(const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat) 
{ 
  double zVal2 = zVal*zVal; 
  double wv2 = wv*wv; 
  double wv3 = wv2*wv; 
  double wv4 = wv3*wv; 
  double wv5 = wv4*wv; 
  double dv2 = dv*dv; 
  double dv3 = dv2*dv; 
  double dv4 = dv3*dv; 
  double dv5 = dv4*dv; 
  double vcut2 = vcut*vcut; 
  double vcut3 = vcut2*vcut; 
  double vcut4 = vcut3*vcut; 
  double vcut5 = vcut4*vcut; 
  double denom = 64*wv5+((-320*vcut)-160*dv)*wv4+(640*vcut2+640*dv*vcut+160*dv2)*wv3+((-640*vcut3)-960*dv*vcut2-480*dv2*vcut-80*dv3)*wv2+(320*vcut4+640*dv*vcut3+480*dv2*vcut2+160*dv3*vcut+20*dv4)*wv-64*vcut5-160*dv*vcut4-160*dv2*vcut3-80*dv3*vcut2-20*dv4*vcut-2*dv5; 
  fhat[0] = -(1.0*((1931.962732559819*f[4]*dv*wv4+((2230.838407415472*f[6]+3863.925465119637*f[4])*dv2-7727.850930239274*f[4]*dv*vcut)*wv3+(11591.77639535891*f[4]*dv*vcut2+((-6692.515222246417*f[6])-11591.77639535891*f[4])*dv2*vcut+(4461.676814830945*f[6]+5903.219460599446*f[4])*dv3)*wv2+((-7727.850930239274*f[4]*dv*vcut3)+(6692.515222246417*f[6]+11591.77639535891*f[4])*dv2*vcut2+((-8923.35362966189*f[6])-11806.43892119889*f[4])*dv3*vcut+(3160.354410505252*f[6]+3112.606624679707*f[4])*dv4)*wv+1931.962732559819*f[4]*dv*vcut4+((-2230.838407415472*f[6])-3863.925465119637*f[4])*dv2*vcut3+(4461.676814830945*f[6]+5903.219460599446*f[4])*dv3*vcut2+((-3160.354410505252*f[6])-3112.606624679707*f[4])*dv4*vcut+(743.6128024718241*f[6]+657.4039853849382*f[4])*dv5)*zVal2+(997.661265159673*f[1]*dv*wv4+((1152.0*f[3]+1995.322530319346*f[1])*dv2-3990.645060638692*f[1]*dv*vcut)*wv3+(5985.967590958037*f[1]*dv*vcut2+((-3456.0*f[3])-5985.967590958037*f[1])*dv2*vcut+(247.8709341572747*f[7]+2304.0*f[3]+3048.409421321224*f[1])*dv3)*wv2+((-3990.645060638692*f[1]*dv*vcut3)+(3456.0*f[3]+5985.967590958037*f[1])*dv2*vcut2+((-495.7418683145494*f[7])-4608.0*f[3]-6096.818842642448*f[1])*dv3*vcut+(495.7418683145494*f[7]+1632.0*f[3]+1607.343149423918*f[1])*dv4)*wv+997.661265159673*f[1]*dv*vcut4+((-1152.0*f[3])-1995.322530319346*f[1])*dv2*vcut3+(247.8709341572747*f[7]+2304.0*f[3]+3048.409421321224*f[1])*dv3*vcut2+((-495.7418683145494*f[7])-1632.0*f[3]-1607.343149423918*f[1])*dv4*vcut+(185.903200617956*f[7]+384.0*f[3]+339.4819582834999*f[1])*dv5)*zVal+(576.0*f[0]-643.9875775199395*f[4])*dv*wv4+((2575.950310079758*f[4]-2304.0*f[0])*dv*vcut+((-743.6128024718241*f[6])-1287.975155039879*f[4]+665.1075101064488*f[2]+1152.0*f[0])*dv2)*wv3+((3456.0*f[0]-3863.925465119637*f[4])*dv*vcut2+(2230.838407415472*f[6]+3863.925465119637*f[4]-1995.322530319346*f[2]-3456.0*f[0])*dv2*vcut+((-1487.225604943648*f[6])+143.1083505599865*f[5]-1967.739820199815*f[4]+1330.215020212898*f[2]+1760.0*f[0])*dv3)*wv2+((2575.950310079758*f[4]-2304.0*f[0])*dv*vcut3+((-2230.838407415472*f[6])-3863.925465119637*f[4]+1995.322530319346*f[2]+3456.0*f[0])*dv2*vcut2+(2974.451209887296*f[6]-286.2167011199731*f[5]+3935.479640399631*f[4]-2660.430040425795*f[2]-3520.0*f[0])*dv3*vcut+((-1053.451470168417*f[6])+286.2167011199731*f[5]-1037.535541559902*f[4]+942.2356393174692*f[2]+928.0*f[0])*dv4)*wv+(576.0*f[0]-643.9875775199395*f[4])*dv*vcut4+(743.6128024718241*f[6]+1287.975155039879*f[4]-665.1075101064488*f[2]-1152.0*f[0])*dv2*vcut3+((-1487.225604943648*f[6])+143.1083505599865*f[5]-1967.739820199815*f[4]+1330.215020212898*f[2]+1760.0*f[0])*dv3*vcut2+(1053.451470168417*f[6]-286.2167011199731*f[5]+1037.535541559902*f[4]-942.2356393174692*f[2]-928.0*f[0])*dv4*vcut+((-247.8709341572747*f[6])+107.3312629199899*f[5]-219.1346617949794*f[4]+221.7025033688163*f[2]+196.0*f[0])*dv5))/denom; 
  fhat[1] = 0.0; 
  fhat[2] = -(0.3333333333333333*((6692.515222246414*f[4]*dv2*wv3+((10303.80124031903*f[6]+13385.03044449283*f[4])*dv3-20077.54566673924*f[4]*dv2*vcut)*wv2+(20077.54566673924*f[4]*dv2*vcut2+((-20607.60248063806*f[6])-26770.06088898566*f[4])*dv3*vcut+(9015.82608527915*f[6]+9481.063231515753*f[4])*dv4)*wv-6692.515222246414*f[4]*dv2*vcut3+(10303.80124031903*f[6]+13385.03044449283*f[4])*dv3*vcut2+((-9015.82608527915*f[6])-9481.063231515753*f[4])*dv4*vcut+(2575.950310079757*f[6]+2230.838407415471*f[4])*dv5)*zVal2+(3456.0*f[1]*dv2*wv3+((5320.86008085159*f[3]+6912.0*f[1])*dv3-10368.0*f[1]*dv2*vcut)*wv2+(10368.0*f[1]*dv2*vcut2+((-10641.72016170318*f[3])-13824.0*f[1])*dv3*vcut+(1287.975155039879*f[7]+4655.752570745141*f[3]+4896.0*f[1])*dv4)*wv-3456.0*f[1]*dv2*vcut3+(5320.86008085159*f[3]+6912.0*f[1])*dv3*vcut2+((-1287.975155039879*f[7])-4655.752570745141*f[3]-4896.0*f[1])*dv4*vcut+(643.9875775199394*f[7]+1330.215020212898*f[3]+1152.0*f[1])*dv5)*zVal+(1995.322530319346*f[0]-2230.838407415471*f[4])*dv2*wv3+((6692.515222246414*f[4]-5985.967590958037*f[0])*dv2*vcut+((-3434.600413439677*f[6])-4461.676814830943*f[4]+3072.0*f[2]+3990.645060638692*f[0])*dv3)*wv2+((5985.967590958037*f[0]-6692.515222246414*f[4])*dv2*vcut2+(6869.200826879353*f[6]+8923.353629661886*f[4]-6144.0*f[2]-7981.290121277384*f[0])*dv3*vcut+((-3005.275361759717*f[6])+743.612802471824*f[5]-3160.354410505252*f[4]+2688.0*f[2]+2826.706917952407*f[0])*dv4)*wv+(2230.838407415471*f[4]-1995.322530319346*f[0])*dv2*vcut3+((-3434.600413439677*f[6])-4461.676814830943*f[4]+3072.0*f[2]+3990.645060638692*f[0])*dv3*vcut2+(3005.275361759717*f[6]-743.612802471824*f[5]+3160.354410505252*f[4]-2688.0*f[2]-2826.706917952407*f[0])*dv4*vcut+((-858.6501033599192*f[6])+371.806401235912*f[5]-743.612802471824*f[4]+768.0*f[2]+665.1075101064488*f[0])*dv5))/denom; 
  fhat[3] = 0.0; 
  fhat[4] = 0.0; 
  fhat[5] = -(0.2*((2400.0*f[4]*dv3*wv2+((4156.921938165307*f[6]+4800.0*f[4])*dv4-4800.0*f[4]*dv3*vcut)*wv+2400.0*f[4]*dv3*vcut2+((-4156.921938165307*f[6])-4800.0*f[4])*dv4*vcut+(2078.460969082653*f[6]+1800.0*f[4])*dv5)*zVal2+(1239.354670786374*f[1]*dv3*wv2+((2146.625258399798*f[3]+2478.709341572747*f[1])*dv4-2478.709341572747*f[1]*dv3*vcut)*wv+1239.354670786374*f[1]*dv3*vcut2+((-2146.625258399798*f[3])-2478.709341572747*f[1])*dv4*vcut+(554.2562584220408*f[7]+1073.312629199899*f[3]+929.5160030897802*f[1])*dv5)*zVal+(715.5417527999329*f[0]-800.0*f[4])*dv3*wv2+((1600.0*f[4]-1431.083505599866*f[0])*dv3*vcut+((-1385.640646055102*f[6])-1600.0*f[4]+1239.354670786374*f[2]+1431.083505599866*f[0])*dv4)*wv+(715.5417527999329*f[0]-800.0*f[4])*dv3*vcut2+(1385.640646055102*f[6]+1600.0*f[4]-1239.354670786374*f[2]-1431.083505599866*f[0])*dv4*vcut+((-692.8203230275511*f[6])+320.0*f[5]-600.0*f[4]+619.6773353931868*f[2]+536.6563145999496*f[0])*dv5))/denom; 
  fhat[6] = 0.0; 
  fhat[7] = 0.0; 
} 

double calcSheathDeltaPhi2xSer_P1(const double *phi, const double *phiWall, const double zVal) 
{ 
  return 0.8660254037844386*phi[2]*zVal-0.5*phiWall[0]+0.5*phi[0]; 
}

void calcSheathPartialReflection2x2vSer_P1(const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat) 
{ 
  double zVal2 = zVal*zVal; 
  double wv2 = wv*wv; 
  double wv3 = wv2*wv; 
  double wv4 = wv3*wv; 
  double wv5 = wv4*wv; 
  double dv2 = dv*dv; 
  double dv3 = dv2*dv; 
  double dv4 = dv3*dv; 
  double dv5 = dv4*dv; 
  double vcut2 = vcut*vcut; 
  double vcut3 = vcut2*vcut; 
  double vcut4 = vcut3*vcut; 
  double vcut5 = vcut4*vcut; 
  double denom = 8*wv3+((-24*vcut)-12*dv)*wv2+(24*vcut2+24*dv*vcut+6*dv2)*wv-8*vcut3-12*dv*vcut2-6*dv2*vcut-dv3; 
  fhat[0] = -(1.0*((55.42562584220407*f[2]*dv*wv2+((24.0*f[7]+27.71281292110204*f[2])*dv2-110.8512516844081*f[2]*dv*vcut)*wv+55.42562584220407*f[2]*dv*vcut2+((-24.0*f[7])-27.71281292110204*f[2])*dv2*vcut+(12.0*f[7]+13.85640646055102*f[2])*dv3)*zVal+32.0*f[0]*dv*wv2+((13.85640646055102*f[3]+16.0*f[0])*dv2-64.0*f[0]*dv*vcut)*wv+32.0*f[0]*dv*vcut2+((-13.85640646055102*f[3])-16.0*f[0])*dv2*vcut+(6.928203230275509*f[3]+8.0*f[0])*dv3))/denom; 
  fhat[1] = -(1.0*((55.42562584220407*f[5]*dv*wv2+((24.0*f[11]+27.71281292110204*f[5])*dv2-110.8512516844081*f[5]*dv*vcut)*wv+55.42562584220407*f[5]*dv*vcut2+((-24.0*f[11])-27.71281292110204*f[5])*dv2*vcut+(12.0*f[11]+13.85640646055102*f[5])*dv3)*zVal+32.0*f[1]*dv*wv2+((13.85640646055102*f[6]+16.0*f[1])*dv2-64.0*f[1]*dv*vcut)*wv+32.0*f[1]*dv*vcut2+((-13.85640646055102*f[6])-16.0*f[1])*dv2*vcut+(6.928203230275509*f[6]+8.0*f[1])*dv3))/denom; 
  fhat[2] = 0.0; 
  fhat[3] = -(1.0*((24.0*f[2]*dv2*wv-24.0*f[2]*dv2*vcut+(13.85640646055102*f[7]+12.0*f[2])*dv3)*zVal+13.85640646055102*f[0]*dv2*wv-13.85640646055102*f[0]*dv2*vcut+(8.0*f[3]+6.928203230275509*f[0])*dv3))/denom; 
  fhat[4] = -(1.0*((55.42562584220407*f[9]*dv*wv2+((24.0*f[14]+27.71281292110204*f[9])*dv2-110.8512516844081*f[9]*dv*vcut)*wv+55.42562584220407*f[9]*dv*vcut2+((-24.0*f[14])-27.71281292110204*f[9])*dv2*vcut+(12.0*f[14]+13.85640646055102*f[9])*dv3)*zVal+32.0*f[4]*dv*wv2+((13.85640646055102*f[10]+16.0*f[4])*dv2-64.0*f[4]*dv*vcut)*wv+32.0*f[4]*dv*vcut2+((-13.85640646055102*f[10])-16.0*f[4])*dv2*vcut+(6.928203230275509*f[10]+8.0*f[4])*dv3))/denom; 
  fhat[5] = 0.0; 
  fhat[6] = -(1.0*((24.0*f[5]*dv2*wv-24.0*f[5]*dv2*vcut+(13.85640646055102*f[11]+12.0*f[5])*dv3)*zVal+13.85640646055102*f[1]*dv2*wv-13.85640646055102*f[1]*dv2*vcut+(8.0*f[6]+6.928203230275509*f[1])*dv3))/denom; 
  fhat[7] = 0.0; 
  fhat[8] = -(1.0*((55.42562584220407*f[12]*dv*wv2+((24.0*f[15]+27.71281292110204*f[12])*dv2-110.8512516844081*f[12]*dv*vcut)*wv+55.42562584220407*f[12]*dv*vcut2+((-24.0*f[15])-27.71281292110204*f[12])*dv2*vcut+(12.0*f[15]+13.85640646055102*f[12])*dv3)*zVal+32.0*f[8]*dv*wv2+((13.85640646055102*f[13]+16.0*f[8])*dv2-64.0*f[8]*dv*vcut)*wv+32.0*f[8]*dv*vcut2+((-13.85640646055102*f[13])-16.0*f[8])*dv2*vcut+(6.928203230275509*f[13]+8.0*f[8])*dv3))/denom; 
  fhat[9] = 0.0; 
  fhat[10] = -(1.0*((24.0*f[9]*dv2*wv-24.0*f[9]*dv2*vcut+(13.85640646055102*f[14]+12.0*f[9])*dv3)*zVal+13.85640646055102*f[4]*dv2*wv-13.85640646055102*f[4]*dv2*vcut+(8.0*f[10]+6.928203230275509*f[4])*dv3))/denom; 
  fhat[11] = 0.0; 
  fhat[12] = 0.0; 
  fhat[13] = -(1.0*((24.0*f[12]*dv2*wv-24.0*f[12]*dv2*vcut+(13.85640646055102*f[15]+12.0*f[12])*dv3)*zVal+13.85640646055102*f[8]*dv2*wv-13.85640646055102*f[8]*dv2*vcut+(8.0*f[13]+6.928203230275509*f[8])*dv3))/denom; 
  fhat[14] = 0.0; 
  fhat[15] = 0.0; 
} 

double calcSheathDeltaPhi3xSer_P1(const double *phi, const double *phiWall, const double zVal) 
{ 
  return 0.6123724356957944*phi[3]*zVal-0.3535533905932737*phiWall[0]+0.3535533905932737*phi[0]; 
}

void calcSheathPartialReflection3x2vSer_P1(const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat) 
{ 
  double zVal2 = zVal*zVal; 
  double wv2 = wv*wv; 
  double wv3 = wv2*wv; 
  double wv4 = wv3*wv; 
  double wv5 = wv4*wv; 
  double dv2 = dv*dv; 
  double dv3 = dv2*dv; 
  double dv4 = dv3*dv; 
  double dv5 = dv4*dv; 
  double vcut2 = vcut*vcut; 
  double vcut3 = vcut2*vcut; 
  double vcut4 = vcut3*vcut; 
  double vcut5 = vcut4*vcut; 
  double denom = 16*wv3+((-48*vcut)-24*dv)*wv2+(48*vcut2+48*dv*vcut+12*dv2)*wv-16*vcut3-24*dv*vcut2-12*dv2*vcut-2*dv3; 
  fhat[0] = -(1.0*((110.8512516844081*f[3]*dv*wv2+((48.0*f[11]+55.42562584220407*f[3])*dv2-221.7025033688163*f[3]*dv*vcut)*wv+110.8512516844081*f[3]*dv*vcut2+((-48.0*f[11])-55.42562584220407*f[3])*dv2*vcut+(24.0*f[11]+27.71281292110204*f[3])*dv3)*zVal+64.0*f[0]*dv*wv2+((27.71281292110204*f[4]+32.0*f[0])*dv2-128.0*f[0]*dv*vcut)*wv+64.0*f[0]*dv*vcut2+((-27.71281292110204*f[4])-32.0*f[0])*dv2*vcut+(13.85640646055102*f[4]+16.0*f[0])*dv3))/denom; 
  fhat[1] = -(1.0*((110.8512516844081*f[7]*dv*wv2+((48.0*f[18]+55.42562584220407*f[7])*dv2-221.7025033688163*f[7]*dv*vcut)*wv+110.8512516844081*f[7]*dv*vcut2+((-48.0*f[18])-55.42562584220407*f[7])*dv2*vcut+(24.0*f[18]+27.71281292110204*f[7])*dv3)*zVal+64.0*f[1]*dv*wv2+((27.71281292110204*f[9]+32.0*f[1])*dv2-128.0*f[1]*dv*vcut)*wv+64.0*f[1]*dv*vcut2+((-27.71281292110204*f[9])-32.0*f[1])*dv2*vcut+(13.85640646055102*f[9]+16.0*f[1])*dv3))/denom; 
  fhat[2] = -(1.0*((110.8512516844081*f[8]*dv*wv2+((48.0*f[19]+55.42562584220407*f[8])*dv2-221.7025033688163*f[8]*dv*vcut)*wv+110.8512516844081*f[8]*dv*vcut2+((-48.0*f[19])-55.42562584220407*f[8])*dv2*vcut+(24.0*f[19]+27.71281292110204*f[8])*dv3)*zVal+64.0*f[2]*dv*wv2+((27.71281292110204*f[10]+32.0*f[2])*dv2-128.0*f[2]*dv*vcut)*wv+64.0*f[2]*dv*vcut2+((-27.71281292110204*f[10])-32.0*f[2])*dv2*vcut+(13.85640646055102*f[10]+16.0*f[2])*dv3))/denom; 
  fhat[3] = 0.0; 
  fhat[4] = -(1.0*((48.0*f[3]*dv2*wv-48.0*f[3]*dv2*vcut+(27.71281292110204*f[11]+24.0*f[3])*dv3)*zVal+27.71281292110204*f[0]*dv2*wv-27.71281292110204*f[0]*dv2*vcut+(16.0*f[4]+13.85640646055102*f[0])*dv3))/denom; 
  fhat[5] = -(1.0*((110.8512516844081*f[14]*dv*wv2+((48.0*f[25]+55.42562584220407*f[14])*dv2-221.7025033688163*f[14]*dv*vcut)*wv+110.8512516844081*f[14]*dv*vcut2+((-48.0*f[25])-55.42562584220407*f[14])*dv2*vcut+(24.0*f[25]+27.71281292110204*f[14])*dv3)*zVal+64.0*f[5]*dv*wv2+((27.71281292110204*f[15]+32.0*f[5])*dv2-128.0*f[5]*dv*vcut)*wv+64.0*f[5]*dv*vcut2+((-27.71281292110204*f[15])-32.0*f[5])*dv2*vcut+(13.85640646055102*f[15]+16.0*f[5])*dv3))/denom; 
  fhat[6] = -(1.0*((110.8512516844081*f[16]*dv*wv2+((48.0*f[26]+55.42562584220407*f[16])*dv2-221.7025033688163*f[16]*dv*vcut)*wv+110.8512516844081*f[16]*dv*vcut2+((-48.0*f[26])-55.42562584220407*f[16])*dv2*vcut+(24.0*f[26]+27.71281292110204*f[16])*dv3)*zVal+64.0*f[6]*dv*wv2+((27.71281292110204*f[17]+32.0*f[6])*dv2-128.0*f[6]*dv*vcut)*wv+64.0*f[6]*dv*vcut2+((-27.71281292110204*f[17])-32.0*f[6])*dv2*vcut+(13.85640646055102*f[17]+16.0*f[6])*dv3))/denom; 
  fhat[7] = 0.0; 
  fhat[8] = 0.0; 
  fhat[9] = -(1.0*((48.0*f[7]*dv2*wv-48.0*f[7]*dv2*vcut+(27.71281292110204*f[18]+24.0*f[7])*dv3)*zVal+27.71281292110204*f[1]*dv2*wv-27.71281292110204*f[1]*dv2*vcut+(16.0*f[9]+13.85640646055102*f[1])*dv3))/denom; 
  fhat[10] = -(1.0*((48.0*f[8]*dv2*wv-48.0*f[8]*dv2*vcut+(27.71281292110204*f[19]+24.0*f[8])*dv3)*zVal+27.71281292110204*f[2]*dv2*wv-27.71281292110204*f[2]*dv2*vcut+(16.0*f[10]+13.85640646055102*f[2])*dv3))/denom; 
  fhat[11] = 0.0; 
  fhat[12] = -(1.0*((110.8512516844081*f[21]*dv*wv2+((48.0*f[29]+55.42562584220407*f[21])*dv2-221.7025033688163*f[21]*dv*vcut)*wv+110.8512516844081*f[21]*dv*vcut2+((-48.0*f[29])-55.42562584220407*f[21])*dv2*vcut+(24.0*f[29]+27.71281292110204*f[21])*dv3)*zVal+64.0*f[12]*dv*wv2+((27.71281292110204*f[23]+32.0*f[12])*dv2-128.0*f[12]*dv*vcut)*wv+64.0*f[12]*dv*vcut2+((-27.71281292110204*f[23])-32.0*f[12])*dv2*vcut+(13.85640646055102*f[23]+16.0*f[12])*dv3))/denom; 
  fhat[13] = -(1.0*((110.8512516844081*f[22]*dv*wv2+((48.0*f[30]+55.42562584220407*f[22])*dv2-221.7025033688163*f[22]*dv*vcut)*wv+110.8512516844081*f[22]*dv*vcut2+((-48.0*f[30])-55.42562584220407*f[22])*dv2*vcut+(24.0*f[30]+27.71281292110204*f[22])*dv3)*zVal+64.0*f[13]*dv*wv2+((27.71281292110204*f[24]+32.0*f[13])*dv2-128.0*f[13]*dv*vcut)*wv+64.0*f[13]*dv*vcut2+((-27.71281292110204*f[24])-32.0*f[13])*dv2*vcut+(13.85640646055102*f[24]+16.0*f[13])*dv3))/denom; 
  fhat[14] = 0.0; 
  fhat[15] = -(1.0*((48.0*f[14]*dv2*wv-48.0*f[14]*dv2*vcut+(27.71281292110204*f[25]+24.0*f[14])*dv3)*zVal+27.71281292110204*f[5]*dv2*wv-27.71281292110204*f[5]*dv2*vcut+(16.0*f[15]+13.85640646055102*f[5])*dv3))/denom; 
  fhat[16] = 0.0; 
  fhat[17] = -(1.0*((48.0*f[16]*dv2*wv-48.0*f[16]*dv2*vcut+(27.71281292110204*f[26]+24.0*f[16])*dv3)*zVal+27.71281292110204*f[6]*dv2*wv-27.71281292110204*f[6]*dv2*vcut+(16.0*f[17]+13.85640646055102*f[6])*dv3))/denom; 
  fhat[18] = 0.0; 
  fhat[19] = 0.0; 
  fhat[20] = -(1.0*((110.8512516844081*f[27]*dv*wv2+((48.0*f[31]+55.42562584220407*f[27])*dv2-221.7025033688163*f[27]*dv*vcut)*wv+110.8512516844081*f[27]*dv*vcut2+((-48.0*f[31])-55.42562584220407*f[27])*dv2*vcut+(24.0*f[31]+27.71281292110204*f[27])*dv3)*zVal+64.0*f[20]*dv*wv2+((27.71281292110204*f[28]+32.0*f[20])*dv2-128.0*f[20]*dv*vcut)*wv+64.0*f[20]*dv*vcut2+((-27.71281292110204*f[28])-32.0*f[20])*dv2*vcut+(13.85640646055102*f[28]+16.0*f[20])*dv3))/denom; 
  fhat[21] = 0.0; 
  fhat[22] = 0.0; 
  fhat[23] = -(1.0*((48.0*f[21]*dv2*wv-48.0*f[21]*dv2*vcut+(27.71281292110204*f[29]+24.0*f[21])*dv3)*zVal+27.71281292110204*f[12]*dv2*wv-27.71281292110204*f[12]*dv2*vcut+(16.0*f[23]+13.85640646055102*f[12])*dv3))/denom; 
  fhat[24] = -(1.0*((48.0*f[22]*dv2*wv-48.0*f[22]*dv2*vcut+(27.71281292110204*f[30]+24.0*f[22])*dv3)*zVal+27.71281292110204*f[13]*dv2*wv-27.71281292110204*f[13]*dv2*vcut+(16.0*f[24]+13.85640646055102*f[13])*dv3))/denom; 
  fhat[25] = 0.0; 
  fhat[26] = 0.0; 
  fhat[27] = 0.0; 
  fhat[28] = -(1.0*((48.0*f[27]*dv2*wv-48.0*f[27]*dv2*vcut+(27.71281292110204*f[31]+24.0*f[27])*dv3)*zVal+27.71281292110204*f[20]*dv2*wv-27.71281292110204*f[20]*dv2*vcut+(16.0*f[28]+13.85640646055102*f[20])*dv3))/denom; 
  fhat[29] = 0.0; 
  fhat[30] = 0.0; 
  fhat[31] = 0.0; 
} 

