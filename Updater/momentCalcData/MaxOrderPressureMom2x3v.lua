local _M = { pressureTensor = {}, energy = {} } 
 
_M.pressureTensor[1] = function (f, out, dv, w) 
   local dv1, dv2, dv3 = dv[1], dv[2], dv[3] 
   local w1, w2, w3 = w[1], w[2], w[3] 
   out[1] = 2.828427124746191*f[1]*w1^2+1.632993161855453*f[4]*dv1*w1+0.2357022603955158*f[1]*dv1^2 
   out[2] = 2.828427124746191*f[2]*w1^2+0.2357022603955158*f[2]*dv1^2 
   out[3] = 2.828427124746191*f[3]*w1^2+0.2357022603955158*f[3]*dv1^2 
   out[4] = 2.828427124746191*f[1]*w1*w2+0.8164965809277261*f[4]*dv1*w2+0.8164965809277261*f[5]*dv2*w1 
   out[5] = 2.828427124746191*f[2]*w1*w2 
   out[6] = 2.828427124746191*f[3]*w1*w2 
   out[7] = 2.828427124746191*f[1]*w1*w3+0.8164965809277261*f[4]*dv1*w3+0.8164965809277261*f[6]*dv3*w1 
   out[8] = 2.828427124746191*f[2]*w1*w3 
   out[9] = 2.828427124746191*f[3]*w1*w3 
   out[10] = 2.828427124746191*f[1]*w2^2+1.632993161855453*f[5]*dv2*w2+0.2357022603955158*f[1]*dv2^2 
   out[11] = 2.828427124746191*f[2]*w2^2+0.2357022603955158*f[2]*dv2^2 
   out[12] = 2.828427124746191*f[3]*w2^2+0.2357022603955158*f[3]*dv2^2 
   out[13] = 2.828427124746191*f[1]*w2*w3+0.8164965809277261*f[5]*dv2*w3+0.8164965809277261*f[6]*dv3*w2 
   out[14] = 2.828427124746191*f[2]*w2*w3 
   out[15] = 2.828427124746191*f[3]*w2*w3 
   out[16] = 2.828427124746191*f[1]*w3^2+1.632993161855453*f[6]*dv3*w3+0.2357022603955158*f[1]*dv3^2 
   out[17] = 2.828427124746191*f[2]*w3^2+0.2357022603955158*f[2]*dv3^2 
   out[18] = 2.828427124746191*f[3]*w3^2+0.2357022603955158*f[3]*dv3^2 
end 
_M.energy[1] = function (f, out, dv, w) 
   local dv1, dv2, dv3 = dv[1], dv[2], dv[3] 
   local w1, w2, w3 = w[1], w[2], w[3] 
   out[1] = 2.828427124746191*f[1]*w3^2+1.632993161855453*f[6]*dv3*w3+2.828427124746191*f[1]*w2^2+1.632993161855453*f[5]*dv2*w2+2.828427124746191*f[1]*w1^2+1.632993161855453*f[4]*dv1*w1+0.2357022603955158*f[1]*dv3^2+0.2357022603955158*f[1]*dv2^2+0.2357022603955158*f[1]*dv1^2 
   out[2] = 2.828427124746191*f[2]*w3^2+2.828427124746191*f[2]*w2^2+2.828427124746191*f[2]*w1^2+0.2357022603955158*f[2]*dv3^2+0.2357022603955158*f[2]*dv2^2+0.2357022603955158*f[2]*dv1^2 
   out[3] = 2.828427124746191*f[3]*w3^2+2.828427124746191*f[3]*w2^2+2.828427124746191*f[3]*w1^2+0.2357022603955158*f[3]*dv3^2+0.2357022603955158*f[3]*dv2^2+0.2357022603955158*f[3]*dv1^2 
end 
_M.pressureTensor[2] = function (f, out, dv, w) 
   local dv1, dv2, dv3 = dv[1], dv[2], dv[3] 
   local w1, w2, w3 = w[1], w[2], w[3] 
   out[1] = 2.828427124746191*f[1]*w1^2+1.632993161855453*f[4]*dv1*w1+0.210818510677892*f[19]*dv1^2+0.2357022603955158*f[1]*dv1^2 
   out[2] = 2.828427124746191*f[2]*w1^2+1.632993161855453*f[8]*dv1*w1+0.2357022603955158*f[2]*dv1^2 
   out[3] = 2.828427124746191*f[3]*w1^2+1.632993161855453*f[9]*dv1*w1+0.2357022603955158*f[3]*dv1^2 
   out[4] = 2.828427124746191*f[7]*w1^2+0.2357022603955158*f[7]*dv1^2 
   out[5] = 2.828427124746191*f[17]*w1^2+0.2357022603955158*f[17]*dv1^2 
   out[6] = 2.828427124746191*f[18]*w1^2+0.2357022603955158*f[18]*dv1^2 
   out[7] = 2.828427124746191*f[1]*w1*w2+0.8164965809277261*f[4]*dv1*w2+0.8164965809277261*f[5]*dv2*w1+0.2357022603955158*f[12]*dv1*dv2 
   out[8] = 2.828427124746191*f[2]*w1*w2+0.8164965809277261*f[8]*dv1*w2+0.8164965809277261*f[10]*dv2*w1 
   out[9] = 2.828427124746191*f[3]*w1*w2+0.8164965809277261*f[9]*dv1*w2+0.8164965809277261*f[11]*dv2*w1 
   out[10] = 2.828427124746191*f[7]*w1*w2 
   out[11] = 2.828427124746191*f[17]*w1*w2 
   out[12] = 2.828427124746191*f[18]*w1*w2 
   out[13] = 2.828427124746191*f[1]*w1*w3+0.8164965809277261*f[4]*dv1*w3+0.8164965809277261*f[6]*dv3*w1+0.2357022603955158*f[15]*dv1*dv3 
   out[14] = 2.828427124746191*f[2]*w1*w3+0.8164965809277261*f[8]*dv1*w3+0.8164965809277261*f[13]*dv3*w1 
   out[15] = 2.828427124746191*f[3]*w1*w3+0.8164965809277261*f[9]*dv1*w3+0.8164965809277261*f[14]*dv3*w1 
   out[16] = 2.828427124746191*f[7]*w1*w3 
   out[17] = 2.828427124746191*f[17]*w1*w3 
   out[18] = 2.828427124746191*f[18]*w1*w3 
   out[19] = 2.828427124746191*f[1]*w2^2+1.632993161855453*f[5]*dv2*w2+0.210818510677892*f[20]*dv2^2+0.2357022603955158*f[1]*dv2^2 
   out[20] = 2.828427124746191*f[2]*w2^2+1.632993161855453*f[10]*dv2*w2+0.2357022603955158*f[2]*dv2^2 
   out[21] = 2.828427124746191*f[3]*w2^2+1.632993161855453*f[11]*dv2*w2+0.2357022603955158*f[3]*dv2^2 
   out[22] = 2.828427124746191*f[7]*w2^2+0.2357022603955158*f[7]*dv2^2 
   out[23] = 2.828427124746191*f[17]*w2^2+0.2357022603955158*f[17]*dv2^2 
   out[24] = 2.828427124746191*f[18]*w2^2+0.2357022603955158*f[18]*dv2^2 
   out[25] = 2.828427124746191*f[1]*w2*w3+0.8164965809277261*f[5]*dv2*w3+0.8164965809277261*f[6]*dv3*w2+0.2357022603955158*f[16]*dv2*dv3 
   out[26] = 2.828427124746191*f[2]*w2*w3+0.8164965809277261*f[10]*dv2*w3+0.8164965809277261*f[13]*dv3*w2 
   out[27] = 2.828427124746191*f[3]*w2*w3+0.8164965809277261*f[11]*dv2*w3+0.8164965809277261*f[14]*dv3*w2 
   out[28] = 2.828427124746191*f[7]*w2*w3 
   out[29] = 2.828427124746191*f[17]*w2*w3 
   out[30] = 2.828427124746191*f[18]*w2*w3 
   out[31] = 2.828427124746191*f[1]*w3^2+1.632993161855453*f[6]*dv3*w3+0.210818510677892*f[21]*dv3^2+0.2357022603955158*f[1]*dv3^2 
   out[32] = 2.828427124746191*f[2]*w3^2+1.632993161855453*f[13]*dv3*w3+0.2357022603955158*f[2]*dv3^2 
   out[33] = 2.828427124746191*f[3]*w3^2+1.632993161855453*f[14]*dv3*w3+0.2357022603955158*f[3]*dv3^2 
   out[34] = 2.828427124746191*f[7]*w3^2+0.2357022603955158*f[7]*dv3^2 
   out[35] = 2.828427124746191*f[17]*w3^2+0.2357022603955158*f[17]*dv3^2 
   out[36] = 2.828427124746191*f[18]*w3^2+0.2357022603955158*f[18]*dv3^2 
end 
_M.energy[2] = function (f, out, dv, w) 
   local dv1, dv2, dv3 = dv[1], dv[2], dv[3] 
   local w1, w2, w3 = w[1], w[2], w[3] 
   out[1] = 2.828427124746191*f[1]*w3^2+1.632993161855453*f[6]*dv3*w3+2.828427124746191*f[1]*w2^2+1.632993161855453*f[5]*dv2*w2+2.828427124746191*f[1]*w1^2+1.632993161855453*f[4]*dv1*w1+0.210818510677892*f[21]*dv3^2+0.2357022603955158*f[1]*dv3^2+0.210818510677892*f[20]*dv2^2+0.2357022603955158*f[1]*dv2^2+0.210818510677892*f[19]*dv1^2+0.2357022603955158*f[1]*dv1^2 
   out[2] = 2.828427124746191*f[2]*w3^2+1.632993161855453*f[13]*dv3*w3+2.828427124746191*f[2]*w2^2+1.632993161855453*f[10]*dv2*w2+2.828427124746191*f[2]*w1^2+1.632993161855453*f[8]*dv1*w1+0.2357022603955158*f[2]*dv3^2+0.2357022603955158*f[2]*dv2^2+0.2357022603955158*f[2]*dv1^2 
   out[3] = 2.828427124746191*f[3]*w3^2+1.632993161855453*f[14]*dv3*w3+2.828427124746191*f[3]*w2^2+1.632993161855453*f[11]*dv2*w2+2.828427124746191*f[3]*w1^2+1.632993161855453*f[9]*dv1*w1+0.2357022603955158*f[3]*dv3^2+0.2357022603955158*f[3]*dv2^2+0.2357022603955158*f[3]*dv1^2 
   out[4] = 2.828427124746191*f[7]*w3^2+2.828427124746191*f[7]*w2^2+2.828427124746191*f[7]*w1^2+0.2357022603955158*f[7]*dv3^2+0.2357022603955158*f[7]*dv2^2+0.2357022603955158*f[7]*dv1^2 
   out[5] = 2.828427124746191*f[17]*w3^2+2.828427124746191*f[17]*w2^2+2.828427124746191*f[17]*w1^2+0.2357022603955158*f[17]*dv3^2+0.2357022603955158*f[17]*dv2^2+0.2357022603955158*f[17]*dv1^2 
   out[6] = 2.828427124746191*f[18]*w3^2+2.828427124746191*f[18]*w2^2+2.828427124746191*f[18]*w1^2+0.2357022603955158*f[18]*dv3^2+0.2357022603955158*f[18]*dv2^2+0.2357022603955158*f[18]*dv1^2 
end 
_M.pressureTensor[3] = function (f, out, dv, w) 
   local dv1, dv2, dv3 = dv[1], dv[2], dv[3] 
   local w1, w2, w3 = w[1], w[2], w[3] 
   out[1] = 2.828427124746191*f[1]*w1^2+1.632993161855453*f[4]*dv1*w1+0.210818510677892*f[19]*dv1^2+0.2357022603955158*f[1]*dv1^2 
   out[2] = 2.828427124746191*f[2]*w1^2+1.632993161855453*f[8]*dv1*w1+0.2108185106778921*f[36]*dv1^2+0.2357022603955158*f[2]*dv1^2 
   out[3] = 2.828427124746191*f[3]*w1^2+1.632993161855453*f[9]*dv1*w1+0.2108185106778921*f[37]*dv1^2+0.2357022603955158*f[3]*dv1^2 
   out[4] = 2.828427124746191*f[7]*w1^2+1.632993161855453*f[22]*dv1*w1+0.2357022603955158*f[7]*dv1^2 
   out[5] = 2.828427124746191*f[17]*w1^2+1.632993161855453*f[34]*dv1*w1+0.2357022603955158*f[17]*dv1^2 
   out[6] = 2.828427124746191*f[18]*w1^2+1.632993161855453*f[35]*dv1*w1+0.2357022603955158*f[18]*dv1^2 
   out[7] = 2.828427124746191*f[32]*w1^2+0.2357022603955158*f[32]*dv1^2 
   out[8] = 2.828427124746191*f[33]*w1^2+0.2357022603955158*f[33]*dv1^2 
   out[9] = 2.828427124746191*f[52]*w1^2+0.2357022603955158*f[52]*dv1^2 
   out[10] = 2.828427124746191*f[53]*w1^2+0.2357022603955158*f[53]*dv1^2 
   out[11] = 2.828427124746191*f[1]*w1*w2+0.8164965809277261*f[4]*dv1*w2+0.8164965809277261*f[5]*dv2*w1+0.2357022603955158*f[12]*dv1*dv2 
   out[12] = 2.828427124746191*f[2]*w1*w2+0.8164965809277261*f[8]*dv1*w2+0.8164965809277261*f[10]*dv2*w1+0.2357022603955158*f[24]*dv1*dv2 
   out[13] = 2.828427124746191*f[3]*w1*w2+0.8164965809277261*f[9]*dv1*w2+0.8164965809277261*f[11]*dv2*w1+0.2357022603955158*f[25]*dv1*dv2 
   out[14] = 2.828427124746191*f[7]*w1*w2+0.8164965809277261*f[22]*dv1*w2+0.8164965809277261*f[23]*dv2*w1 
   out[15] = 2.828427124746191*f[17]*w1*w2+0.816496580927726*f[34]*dv1*w2+0.816496580927726*f[38]*dv2*w1 
   out[16] = 2.828427124746191*f[18]*w1*w2+0.816496580927726*f[35]*dv1*w2+0.816496580927726*f[39]*dv2*w1 
   out[17] = 2.828427124746191*f[32]*w1*w2 
   out[18] = 2.828427124746191*f[33]*w1*w2 
   out[19] = 2.828427124746191*f[52]*w1*w2 
   out[20] = 2.828427124746191*f[53]*w1*w2 
   out[21] = 2.828427124746191*f[1]*w1*w3+0.8164965809277261*f[4]*dv1*w3+0.8164965809277261*f[6]*dv3*w1+0.2357022603955158*f[15]*dv1*dv3 
   out[22] = 2.828427124746191*f[2]*w1*w3+0.8164965809277261*f[8]*dv1*w3+0.8164965809277261*f[13]*dv3*w1+0.2357022603955158*f[27]*dv1*dv3 
   out[23] = 2.828427124746191*f[3]*w1*w3+0.8164965809277261*f[9]*dv1*w3+0.8164965809277261*f[14]*dv3*w1+0.2357022603955158*f[28]*dv1*dv3 
   out[24] = 2.828427124746191*f[7]*w1*w3+0.8164965809277261*f[22]*dv1*w3+0.8164965809277261*f[26]*dv3*w1 
   out[25] = 2.828427124746191*f[17]*w1*w3+0.816496580927726*f[34]*dv1*w3+0.816496580927726*f[44]*dv3*w1 
   out[26] = 2.828427124746191*f[18]*w1*w3+0.816496580927726*f[35]*dv1*w3+0.816496580927726*f[45]*dv3*w1 
   out[27] = 2.828427124746191*f[32]*w1*w3 
   out[28] = 2.828427124746191*f[33]*w1*w3 
   out[29] = 2.828427124746191*f[52]*w1*w3 
   out[30] = 2.828427124746191*f[53]*w1*w3 
   out[31] = 2.828427124746191*f[1]*w2^2+1.632993161855453*f[5]*dv2*w2+0.210818510677892*f[20]*dv2^2+0.2357022603955158*f[1]*dv2^2 
   out[32] = 2.828427124746191*f[2]*w2^2+1.632993161855453*f[10]*dv2*w2+0.2108185106778921*f[41]*dv2^2+0.2357022603955158*f[2]*dv2^2 
   out[33] = 2.828427124746191*f[3]*w2^2+1.632993161855453*f[11]*dv2*w2+0.2108185106778921*f[42]*dv2^2+0.2357022603955158*f[3]*dv2^2 
   out[34] = 2.828427124746191*f[7]*w2^2+1.632993161855453*f[23]*dv2*w2+0.2357022603955158*f[7]*dv2^2 
   out[35] = 2.828427124746191*f[17]*w2^2+1.632993161855453*f[38]*dv2*w2+0.2357022603955158*f[17]*dv2^2 
   out[36] = 2.828427124746191*f[18]*w2^2+1.632993161855453*f[39]*dv2*w2+0.2357022603955158*f[18]*dv2^2 
   out[37] = 2.828427124746191*f[32]*w2^2+0.2357022603955158*f[32]*dv2^2 
   out[38] = 2.828427124746191*f[33]*w2^2+0.2357022603955158*f[33]*dv2^2 
   out[39] = 2.828427124746191*f[52]*w2^2+0.2357022603955158*f[52]*dv2^2 
   out[40] = 2.828427124746191*f[53]*w2^2+0.2357022603955158*f[53]*dv2^2 
   out[41] = 2.828427124746191*f[1]*w2*w3+0.8164965809277261*f[5]*dv2*w3+0.8164965809277261*f[6]*dv3*w2+0.2357022603955158*f[16]*dv2*dv3 
   out[42] = 2.828427124746191*f[2]*w2*w3+0.8164965809277261*f[10]*dv2*w3+0.8164965809277261*f[13]*dv3*w2+0.2357022603955158*f[29]*dv2*dv3 
   out[43] = 2.828427124746191*f[3]*w2*w3+0.8164965809277261*f[11]*dv2*w3+0.8164965809277261*f[14]*dv3*w2+0.2357022603955158*f[30]*dv2*dv3 
   out[44] = 2.828427124746191*f[7]*w2*w3+0.8164965809277261*f[23]*dv2*w3+0.8164965809277261*f[26]*dv3*w2 
   out[45] = 2.828427124746191*f[17]*w2*w3+0.816496580927726*f[38]*dv2*w3+0.816496580927726*f[44]*dv3*w2 
   out[46] = 2.828427124746191*f[18]*w2*w3+0.816496580927726*f[39]*dv2*w3+0.816496580927726*f[45]*dv3*w2 
   out[47] = 2.828427124746191*f[32]*w2*w3 
   out[48] = 2.828427124746191*f[33]*w2*w3 
   out[49] = 2.828427124746191*f[52]*w2*w3 
   out[50] = 2.828427124746191*f[53]*w2*w3 
   out[51] = 2.828427124746191*f[1]*w3^2+1.632993161855453*f[6]*dv3*w3+0.210818510677892*f[21]*dv3^2+0.2357022603955158*f[1]*dv3^2 
   out[52] = 2.828427124746191*f[2]*w3^2+1.632993161855453*f[13]*dv3*w3+0.2108185106778921*f[48]*dv3^2+0.2357022603955158*f[2]*dv3^2 
   out[53] = 2.828427124746191*f[3]*w3^2+1.632993161855453*f[14]*dv3*w3+0.2108185106778921*f[49]*dv3^2+0.2357022603955158*f[3]*dv3^2 
   out[54] = 2.828427124746191*f[7]*w3^2+1.632993161855453*f[26]*dv3*w3+0.2357022603955158*f[7]*dv3^2 
   out[55] = 2.828427124746191*f[17]*w3^2+1.632993161855453*f[44]*dv3*w3+0.2357022603955158*f[17]*dv3^2 
   out[56] = 2.828427124746191*f[18]*w3^2+1.632993161855453*f[45]*dv3*w3+0.2357022603955158*f[18]*dv3^2 
   out[57] = 2.828427124746191*f[32]*w3^2+0.2357022603955158*f[32]*dv3^2 
   out[58] = 2.828427124746191*f[33]*w3^2+0.2357022603955158*f[33]*dv3^2 
   out[59] = 2.828427124746191*f[52]*w3^2+0.2357022603955158*f[52]*dv3^2 
   out[60] = 2.828427124746191*f[53]*w3^2+0.2357022603955158*f[53]*dv3^2 
end 
_M.energy[3] = function (f, out, dv, w) 
   local dv1, dv2, dv3 = dv[1], dv[2], dv[3] 
   local w1, w2, w3 = w[1], w[2], w[3] 
   out[1] = 2.828427124746191*f[1]*w3^2+1.632993161855453*f[6]*dv3*w3+2.828427124746191*f[1]*w2^2+1.632993161855453*f[5]*dv2*w2+2.828427124746191*f[1]*w1^2+1.632993161855453*f[4]*dv1*w1+0.210818510677892*f[21]*dv3^2+0.2357022603955158*f[1]*dv3^2+0.210818510677892*f[20]*dv2^2+0.2357022603955158*f[1]*dv2^2+0.210818510677892*f[19]*dv1^2+0.2357022603955158*f[1]*dv1^2 
   out[2] = 2.828427124746191*f[2]*w3^2+1.632993161855453*f[13]*dv3*w3+2.828427124746191*f[2]*w2^2+1.632993161855453*f[10]*dv2*w2+2.828427124746191*f[2]*w1^2+1.632993161855453*f[8]*dv1*w1+0.2108185106778921*f[48]*dv3^2+0.2357022603955158*f[2]*dv3^2+0.2108185106778921*f[41]*dv2^2+0.2357022603955158*f[2]*dv2^2+0.2108185106778921*f[36]*dv1^2+0.2357022603955158*f[2]*dv1^2 
   out[3] = 2.828427124746191*f[3]*w3^2+1.632993161855453*f[14]*dv3*w3+2.828427124746191*f[3]*w2^2+1.632993161855453*f[11]*dv2*w2+2.828427124746191*f[3]*w1^2+1.632993161855453*f[9]*dv1*w1+0.2108185106778921*f[49]*dv3^2+0.2357022603955158*f[3]*dv3^2+0.2108185106778921*f[42]*dv2^2+0.2357022603955158*f[3]*dv2^2+0.2108185106778921*f[37]*dv1^2+0.2357022603955158*f[3]*dv1^2 
   out[4] = 2.828427124746191*f[7]*w3^2+1.632993161855453*f[26]*dv3*w3+2.828427124746191*f[7]*w2^2+1.632993161855453*f[23]*dv2*w2+2.828427124746191*f[7]*w1^2+1.632993161855453*f[22]*dv1*w1+0.2357022603955158*f[7]*dv3^2+0.2357022603955158*f[7]*dv2^2+0.2357022603955158*f[7]*dv1^2 
   out[5] = 2.828427124746191*f[17]*w3^2+1.632993161855453*f[44]*dv3*w3+2.828427124746191*f[17]*w2^2+1.632993161855453*f[38]*dv2*w2+2.828427124746191*f[17]*w1^2+1.632993161855453*f[34]*dv1*w1+0.2357022603955158*f[17]*dv3^2+0.2357022603955158*f[17]*dv2^2+0.2357022603955158*f[17]*dv1^2 
   out[6] = 2.828427124746191*f[18]*w3^2+1.632993161855453*f[45]*dv3*w3+2.828427124746191*f[18]*w2^2+1.632993161855453*f[39]*dv2*w2+2.828427124746191*f[18]*w1^2+1.632993161855453*f[35]*dv1*w1+0.2357022603955158*f[18]*dv3^2+0.2357022603955158*f[18]*dv2^2+0.2357022603955158*f[18]*dv1^2 
   out[7] = 2.828427124746191*f[32]*w3^2+2.828427124746191*f[32]*w2^2+2.828427124746191*f[32]*w1^2+0.2357022603955158*f[32]*dv3^2+0.2357022603955158*f[32]*dv2^2+0.2357022603955158*f[32]*dv1^2 
   out[8] = 2.828427124746191*f[33]*w3^2+2.828427124746191*f[33]*w2^2+2.828427124746191*f[33]*w1^2+0.2357022603955158*f[33]*dv3^2+0.2357022603955158*f[33]*dv2^2+0.2357022603955158*f[33]*dv1^2 
   out[9] = 2.828427124746191*f[52]*w3^2+2.828427124746191*f[52]*w2^2+2.828427124746191*f[52]*w1^2+0.2357022603955158*f[52]*dv3^2+0.2357022603955158*f[52]*dv2^2+0.2357022603955158*f[52]*dv1^2 
   out[10] = 2.828427124746191*f[53]*w3^2+2.828427124746191*f[53]*w2^2+2.828427124746191*f[53]*w1^2+0.2357022603955158*f[53]*dv3^2+0.2357022603955158*f[53]*dv2^2+0.2357022603955158*f[53]*dv1^2 
end 
_M.pressureTensor[4] = function (f, out, dv, w) 
   local dv1, dv2, dv3 = dv[1], dv[2], dv[3] 
   local w1, w2, w3 = w[1], w[2], w[3] 
   out[1] = 2.828427124746191*f[1]*w1^2+1.632993161855453*f[4]*dv1*w1+0.210818510677892*f[19]*dv1^2+0.2357022603955158*f[1]*dv1^2 
   out[2] = 2.828427124746191*f[2]*w1^2+1.632993161855453*f[8]*dv1*w1+0.2108185106778921*f[36]*dv1^2+0.2357022603955158*f[2]*dv1^2 
   out[3] = 2.828427124746191*f[3]*w1^2+1.632993161855453*f[9]*dv1*w1+0.2108185106778921*f[37]*dv1^2+0.2357022603955158*f[3]*dv1^2 
   out[4] = 2.828427124746191*f[7]*w1^2+1.632993161855453*f[22]*dv1*w1+0.210818510677892*f[64]*dv1^2+0.2357022603955158*f[7]*dv1^2 
   out[5] = 2.828427124746191*f[17]*w1^2+1.632993161855453*f[34]*dv1*w1+0.210818510677892*f[93]*dv1^2+0.2357022603955158*f[17]*dv1^2 
   out[6] = 2.828427124746191*f[18]*w1^2+1.632993161855453*f[35]*dv1*w1+0.210818510677892*f[94]*dv1^2+0.2357022603955158*f[18]*dv1^2 
   out[7] = 2.828427124746191*f[32]*w1^2+1.632993161855453*f[62]*dv1*w1+0.2357022603955158*f[32]*dv1^2 
   out[8] = 2.828427124746191*f[33]*w1^2+1.632993161855453*f[63]*dv1*w1+0.2357022603955158*f[33]*dv1^2 
   out[9] = 2.828427124746191*f[52]*w1^2+1.632993161855452*f[104]*dv1*w1+0.2357022603955158*f[52]*dv1^2 
   out[10] = 2.828427124746191*f[53]*w1^2+1.632993161855452*f[105]*dv1*w1+0.2357022603955158*f[53]*dv1^2 
   out[11] = 2.828427124746191*f[92]*w1^2+0.2357022603955158*f[92]*dv1^2 
   out[12] = 2.828427124746191*f[102]*w1^2+0.2357022603955158*f[102]*dv1^2 
   out[13] = 2.828427124746191*f[103]*w1^2+0.2357022603955158*f[103]*dv1^2 
   out[14] = 2.828427124746191*f[122]*w1^2+0.2357022603955158*f[122]*dv1^2 
   out[15] = 2.828427124746191*f[123]*w1^2+0.2357022603955158*f[123]*dv1^2 
   out[16] = 2.828427124746191*f[1]*w1*w2+0.8164965809277261*f[4]*dv1*w2+0.8164965809277261*f[5]*dv2*w1+0.2357022603955158*f[12]*dv1*dv2 
   out[17] = 2.828427124746191*f[2]*w1*w2+0.8164965809277261*f[8]*dv1*w2+0.8164965809277261*f[10]*dv2*w1+0.2357022603955158*f[24]*dv1*dv2 
   out[18] = 2.828427124746191*f[3]*w1*w2+0.8164965809277261*f[9]*dv1*w2+0.8164965809277261*f[11]*dv2*w1+0.2357022603955158*f[25]*dv1*dv2 
   out[19] = 2.828427124746191*f[7]*w1*w2+0.8164965809277261*f[22]*dv1*w2+0.8164965809277261*f[23]*dv2*w1+0.2357022603955158*f[57]*dv1*dv2 
   out[20] = 2.828427124746191*f[17]*w1*w2+0.816496580927726*f[34]*dv1*w2+0.816496580927726*f[38]*dv2*w1+0.2357022603955158*f[67]*dv1*dv2 
   out[21] = 2.828427124746191*f[18]*w1*w2+0.816496580927726*f[35]*dv1*w2+0.816496580927726*f[39]*dv2*w1+0.2357022603955158*f[68]*dv1*dv2 
   out[22] = 2.828427124746191*f[32]*w1*w2+0.816496580927726*f[62]*dv1*w2+0.816496580927726*f[65]*dv2*w1 
   out[23] = 2.828427124746191*f[33]*w1*w2+0.816496580927726*f[63]*dv1*w2+0.816496580927726*f[66]*dv2*w1 
   out[24] = 2.828427124746191*f[52]*w1*w2+0.8164965809277258*f[104]*dv1*w2+0.8164965809277258*f[108]*dv2*w1 
   out[25] = 2.828427124746191*f[53]*w1*w2+0.8164965809277258*f[105]*dv1*w2+0.8164965809277258*f[109]*dv2*w1 
   out[26] = 2.828427124746191*f[92]*w1*w2 
   out[27] = 2.828427124746191*f[102]*w1*w2 
   out[28] = 2.828427124746191*f[103]*w1*w2 
   out[29] = 2.828427124746191*f[122]*w1*w2 
   out[30] = 2.828427124746191*f[123]*w1*w2 
   out[31] = 2.828427124746191*f[1]*w1*w3+0.8164965809277261*f[4]*dv1*w3+0.8164965809277261*f[6]*dv3*w1+0.2357022603955158*f[15]*dv1*dv3 
   out[32] = 2.828427124746191*f[2]*w1*w3+0.8164965809277261*f[8]*dv1*w3+0.8164965809277261*f[13]*dv3*w1+0.2357022603955158*f[27]*dv1*dv3 
   out[33] = 2.828427124746191*f[3]*w1*w3+0.8164965809277261*f[9]*dv1*w3+0.8164965809277261*f[14]*dv3*w1+0.2357022603955158*f[28]*dv1*dv3 
   out[34] = 2.828427124746191*f[7]*w1*w3+0.8164965809277261*f[22]*dv1*w3+0.8164965809277261*f[26]*dv3*w1+0.2357022603955158*f[58]*dv1*dv3 
   out[35] = 2.828427124746191*f[17]*w1*w3+0.816496580927726*f[34]*dv1*w3+0.816496580927726*f[44]*dv3*w1+0.2357022603955158*f[76]*dv1*dv3 
   out[36] = 2.828427124746191*f[18]*w1*w3+0.816496580927726*f[35]*dv1*w3+0.816496580927726*f[45]*dv3*w1+0.2357022603955158*f[77]*dv1*dv3 
   out[37] = 2.828427124746191*f[32]*w1*w3+0.816496580927726*f[62]*dv1*w3+0.816496580927726*f[74]*dv3*w1 
   out[38] = 2.828427124746191*f[33]*w1*w3+0.816496580927726*f[63]*dv1*w3+0.816496580927726*f[75]*dv3*w1 
   out[39] = 2.828427124746191*f[52]*w1*w3+0.8164965809277258*f[104]*dv1*w3+0.8164965809277258*f[114]*dv3*w1 
   out[40] = 2.828427124746191*f[53]*w1*w3+0.8164965809277258*f[105]*dv1*w3+0.8164965809277258*f[115]*dv3*w1 
   out[41] = 2.828427124746191*f[92]*w1*w3 
   out[42] = 2.828427124746191*f[102]*w1*w3 
   out[43] = 2.828427124746191*f[103]*w1*w3 
   out[44] = 2.828427124746191*f[122]*w1*w3 
   out[45] = 2.828427124746191*f[123]*w1*w3 
   out[46] = 2.828427124746191*f[1]*w2^2+1.632993161855453*f[5]*dv2*w2+0.210818510677892*f[20]*dv2^2+0.2357022603955158*f[1]*dv2^2 
   out[47] = 2.828427124746191*f[2]*w2^2+1.632993161855453*f[10]*dv2*w2+0.2108185106778921*f[41]*dv2^2+0.2357022603955158*f[2]*dv2^2 
   out[48] = 2.828427124746191*f[3]*w2^2+1.632993161855453*f[11]*dv2*w2+0.2108185106778921*f[42]*dv2^2+0.2357022603955158*f[3]*dv2^2 
   out[49] = 2.828427124746191*f[7]*w2^2+1.632993161855453*f[23]*dv2*w2+0.210818510677892*f[71]*dv2^2+0.2357022603955158*f[7]*dv2^2 
   out[50] = 2.828427124746191*f[17]*w2^2+1.632993161855453*f[38]*dv2*w2+0.210818510677892*f[95]*dv2^2+0.2357022603955158*f[17]*dv2^2 
   out[51] = 2.828427124746191*f[18]*w2^2+1.632993161855453*f[39]*dv2*w2+0.210818510677892*f[96]*dv2^2+0.2357022603955158*f[18]*dv2^2 
   out[52] = 2.828427124746191*f[32]*w2^2+1.632993161855453*f[65]*dv2*w2+0.2357022603955158*f[32]*dv2^2 
   out[53] = 2.828427124746191*f[33]*w2^2+1.632993161855453*f[66]*dv2*w2+0.2357022603955158*f[33]*dv2^2 
   out[54] = 2.828427124746191*f[52]*w2^2+1.632993161855452*f[108]*dv2*w2+0.2357022603955158*f[52]*dv2^2 
   out[55] = 2.828427124746191*f[53]*w2^2+1.632993161855452*f[109]*dv2*w2+0.2357022603955158*f[53]*dv2^2 
   out[56] = 2.828427124746191*f[92]*w2^2+0.2357022603955158*f[92]*dv2^2 
   out[57] = 2.828427124746191*f[102]*w2^2+0.2357022603955158*f[102]*dv2^2 
   out[58] = 2.828427124746191*f[103]*w2^2+0.2357022603955158*f[103]*dv2^2 
   out[59] = 2.828427124746191*f[122]*w2^2+0.2357022603955158*f[122]*dv2^2 
   out[60] = 2.828427124746191*f[123]*w2^2+0.2357022603955158*f[123]*dv2^2 
   out[61] = 2.828427124746191*f[1]*w2*w3+0.8164965809277261*f[5]*dv2*w3+0.8164965809277261*f[6]*dv3*w2+0.2357022603955158*f[16]*dv2*dv3 
   out[62] = 2.828427124746191*f[2]*w2*w3+0.8164965809277261*f[10]*dv2*w3+0.8164965809277261*f[13]*dv3*w2+0.2357022603955158*f[29]*dv2*dv3 
   out[63] = 2.828427124746191*f[3]*w2*w3+0.8164965809277261*f[11]*dv2*w3+0.8164965809277261*f[14]*dv3*w2+0.2357022603955158*f[30]*dv2*dv3 
   out[64] = 2.828427124746191*f[7]*w2*w3+0.8164965809277261*f[23]*dv2*w3+0.8164965809277261*f[26]*dv3*w2+0.2357022603955158*f[59]*dv2*dv3 
   out[65] = 2.828427124746191*f[17]*w2*w3+0.816496580927726*f[38]*dv2*w3+0.816496580927726*f[44]*dv3*w2+0.2357022603955158*f[80]*dv2*dv3 
   out[66] = 2.828427124746191*f[18]*w2*w3+0.816496580927726*f[39]*dv2*w3+0.816496580927726*f[45]*dv3*w2+0.2357022603955158*f[81]*dv2*dv3 
   out[67] = 2.828427124746191*f[32]*w2*w3+0.816496580927726*f[65]*dv2*w3+0.816496580927726*f[74]*dv3*w2 
   out[68] = 2.828427124746191*f[33]*w2*w3+0.816496580927726*f[66]*dv2*w3+0.816496580927726*f[75]*dv3*w2 
   out[69] = 2.828427124746191*f[52]*w2*w3+0.8164965809277258*f[108]*dv2*w3+0.8164965809277258*f[114]*dv3*w2 
   out[70] = 2.828427124746191*f[53]*w2*w3+0.8164965809277258*f[109]*dv2*w3+0.8164965809277258*f[115]*dv3*w2 
   out[71] = 2.828427124746191*f[92]*w2*w3 
   out[72] = 2.828427124746191*f[102]*w2*w3 
   out[73] = 2.828427124746191*f[103]*w2*w3 
   out[74] = 2.828427124746191*f[122]*w2*w3 
   out[75] = 2.828427124746191*f[123]*w2*w3 
   out[76] = 2.828427124746191*f[1]*w3^2+1.632993161855453*f[6]*dv3*w3+0.210818510677892*f[21]*dv3^2+0.2357022603955158*f[1]*dv3^2 
   out[77] = 2.828427124746191*f[2]*w3^2+1.632993161855453*f[13]*dv3*w3+0.2108185106778921*f[48]*dv3^2+0.2357022603955158*f[2]*dv3^2 
   out[78] = 2.828427124746191*f[3]*w3^2+1.632993161855453*f[14]*dv3*w3+0.2108185106778921*f[49]*dv3^2+0.2357022603955158*f[3]*dv3^2 
   out[79] = 2.828427124746191*f[7]*w3^2+1.632993161855453*f[26]*dv3*w3+0.210818510677892*f[86]*dv3^2+0.2357022603955158*f[7]*dv3^2 
   out[80] = 2.828427124746191*f[17]*w3^2+1.632993161855453*f[44]*dv3*w3+0.210818510677892*f[98]*dv3^2+0.2357022603955158*f[17]*dv3^2 
   out[81] = 2.828427124746191*f[18]*w3^2+1.632993161855453*f[45]*dv3*w3+0.210818510677892*f[99]*dv3^2+0.2357022603955158*f[18]*dv3^2 
   out[82] = 2.828427124746191*f[32]*w3^2+1.632993161855453*f[74]*dv3*w3+0.2357022603955158*f[32]*dv3^2 
   out[83] = 2.828427124746191*f[33]*w3^2+1.632993161855453*f[75]*dv3*w3+0.2357022603955158*f[33]*dv3^2 
   out[84] = 2.828427124746191*f[52]*w3^2+1.632993161855452*f[114]*dv3*w3+0.2357022603955158*f[52]*dv3^2 
   out[85] = 2.828427124746191*f[53]*w3^2+1.632993161855452*f[115]*dv3*w3+0.2357022603955158*f[53]*dv3^2 
   out[86] = 2.828427124746191*f[92]*w3^2+0.2357022603955158*f[92]*dv3^2 
   out[87] = 2.828427124746191*f[102]*w3^2+0.2357022603955158*f[102]*dv3^2 
   out[88] = 2.828427124746191*f[103]*w3^2+0.2357022603955158*f[103]*dv3^2 
   out[89] = 2.828427124746191*f[122]*w3^2+0.2357022603955158*f[122]*dv3^2 
   out[90] = 2.828427124746191*f[123]*w3^2+0.2357022603955158*f[123]*dv3^2 
end 
_M.energy[4] = function (f, out, dv, w) 
   local dv1, dv2, dv3 = dv[1], dv[2], dv[3] 
   local w1, w2, w3 = w[1], w[2], w[3] 
   out[1] = 2.828427124746191*f[1]*w3^2+1.632993161855453*f[6]*dv3*w3+2.828427124746191*f[1]*w2^2+1.632993161855453*f[5]*dv2*w2+2.828427124746191*f[1]*w1^2+1.632993161855453*f[4]*dv1*w1+0.210818510677892*f[21]*dv3^2+0.2357022603955158*f[1]*dv3^2+0.210818510677892*f[20]*dv2^2+0.2357022603955158*f[1]*dv2^2+0.210818510677892*f[19]*dv1^2+0.2357022603955158*f[1]*dv1^2 
   out[2] = 2.828427124746191*f[2]*w3^2+1.632993161855453*f[13]*dv3*w3+2.828427124746191*f[2]*w2^2+1.632993161855453*f[10]*dv2*w2+2.828427124746191*f[2]*w1^2+1.632993161855453*f[8]*dv1*w1+0.2108185106778921*f[48]*dv3^2+0.2357022603955158*f[2]*dv3^2+0.2108185106778921*f[41]*dv2^2+0.2357022603955158*f[2]*dv2^2+0.2108185106778921*f[36]*dv1^2+0.2357022603955158*f[2]*dv1^2 
   out[3] = 2.828427124746191*f[3]*w3^2+1.632993161855453*f[14]*dv3*w3+2.828427124746191*f[3]*w2^2+1.632993161855453*f[11]*dv2*w2+2.828427124746191*f[3]*w1^2+1.632993161855453*f[9]*dv1*w1+0.2108185106778921*f[49]*dv3^2+0.2357022603955158*f[3]*dv3^2+0.2108185106778921*f[42]*dv2^2+0.2357022603955158*f[3]*dv2^2+0.2108185106778921*f[37]*dv1^2+0.2357022603955158*f[3]*dv1^2 
   out[4] = 2.828427124746191*f[7]*w3^2+1.632993161855453*f[26]*dv3*w3+2.828427124746191*f[7]*w2^2+1.632993161855453*f[23]*dv2*w2+2.828427124746191*f[7]*w1^2+1.632993161855453*f[22]*dv1*w1+0.210818510677892*f[86]*dv3^2+0.2357022603955158*f[7]*dv3^2+0.210818510677892*f[71]*dv2^2+0.2357022603955158*f[7]*dv2^2+0.210818510677892*f[64]*dv1^2+0.2357022603955158*f[7]*dv1^2 
   out[5] = 2.828427124746191*f[17]*w3^2+1.632993161855453*f[44]*dv3*w3+2.828427124746191*f[17]*w2^2+1.632993161855453*f[38]*dv2*w2+2.828427124746191*f[17]*w1^2+1.632993161855453*f[34]*dv1*w1+0.210818510677892*f[98]*dv3^2+0.2357022603955158*f[17]*dv3^2+0.210818510677892*f[95]*dv2^2+0.2357022603955158*f[17]*dv2^2+0.210818510677892*f[93]*dv1^2+0.2357022603955158*f[17]*dv1^2 
   out[6] = 2.828427124746191*f[18]*w3^2+1.632993161855453*f[45]*dv3*w3+2.828427124746191*f[18]*w2^2+1.632993161855453*f[39]*dv2*w2+2.828427124746191*f[18]*w1^2+1.632993161855453*f[35]*dv1*w1+0.210818510677892*f[99]*dv3^2+0.2357022603955158*f[18]*dv3^2+0.210818510677892*f[96]*dv2^2+0.2357022603955158*f[18]*dv2^2+0.210818510677892*f[94]*dv1^2+0.2357022603955158*f[18]*dv1^2 
   out[7] = 2.828427124746191*f[32]*w3^2+1.632993161855453*f[74]*dv3*w3+2.828427124746191*f[32]*w2^2+1.632993161855453*f[65]*dv2*w2+2.828427124746191*f[32]*w1^2+1.632993161855453*f[62]*dv1*w1+0.2357022603955158*f[32]*dv3^2+0.2357022603955158*f[32]*dv2^2+0.2357022603955158*f[32]*dv1^2 
   out[8] = 2.828427124746191*f[33]*w3^2+1.632993161855453*f[75]*dv3*w3+2.828427124746191*f[33]*w2^2+1.632993161855453*f[66]*dv2*w2+2.828427124746191*f[33]*w1^2+1.632993161855453*f[63]*dv1*w1+0.2357022603955158*f[33]*dv3^2+0.2357022603955158*f[33]*dv2^2+0.2357022603955158*f[33]*dv1^2 
   out[9] = 2.828427124746191*f[52]*w3^2+1.632993161855452*f[114]*dv3*w3+2.828427124746191*f[52]*w2^2+1.632993161855452*f[108]*dv2*w2+2.828427124746191*f[52]*w1^2+1.632993161855452*f[104]*dv1*w1+0.2357022603955158*f[52]*dv3^2+0.2357022603955158*f[52]*dv2^2+0.2357022603955158*f[52]*dv1^2 
   out[10] = 2.828427124746191*f[53]*w3^2+1.632993161855452*f[115]*dv3*w3+2.828427124746191*f[53]*w2^2+1.632993161855452*f[109]*dv2*w2+2.828427124746191*f[53]*w1^2+1.632993161855452*f[105]*dv1*w1+0.2357022603955158*f[53]*dv3^2+0.2357022603955158*f[53]*dv2^2+0.2357022603955158*f[53]*dv1^2 
   out[11] = 2.828427124746191*f[92]*w3^2+2.828427124746191*f[92]*w2^2+2.828427124746191*f[92]*w1^2+0.2357022603955158*f[92]*dv3^2+0.2357022603955158*f[92]*dv2^2+0.2357022603955158*f[92]*dv1^2 
   out[12] = 2.828427124746191*f[102]*w3^2+2.828427124746191*f[102]*w2^2+2.828427124746191*f[102]*w1^2+0.2357022603955158*f[102]*dv3^2+0.2357022603955158*f[102]*dv2^2+0.2357022603955158*f[102]*dv1^2 
   out[13] = 2.828427124746191*f[103]*w3^2+2.828427124746191*f[103]*w2^2+2.828427124746191*f[103]*w1^2+0.2357022603955158*f[103]*dv3^2+0.2357022603955158*f[103]*dv2^2+0.2357022603955158*f[103]*dv1^2 
   out[14] = 2.828427124746191*f[122]*w3^2+2.828427124746191*f[122]*w2^2+2.828427124746191*f[122]*w1^2+0.2357022603955158*f[122]*dv3^2+0.2357022603955158*f[122]*dv2^2+0.2357022603955158*f[122]*dv1^2 
   out[15] = 2.828427124746191*f[123]*w3^2+2.828427124746191*f[123]*w2^2+2.828427124746191*f[123]*w1^2+0.2357022603955158*f[123]*dv3^2+0.2357022603955158*f[123]*dv2^2+0.2357022603955158*f[123]*dv1^2 
end 
return _M 
