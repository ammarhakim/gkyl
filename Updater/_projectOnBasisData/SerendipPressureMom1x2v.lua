local _M = { pressureTensor = {}, energy = {} } 
 
_M.pressureTensor[1] = function (f, out, dv, w) 
   local dv1, dv2 = dv[1], dv[2] 
   local w1, w2 = w[1], w[2] 
   out[1] = 2.0*f[1]*w1^2+1.154700538379252*f[3]*dv1*w1+0.1666666666666667*f[1]*dv1^2 
   out[2] = 2.0*f[2]*w1^2+1.154700538379252*f[5]*dv1*w1+0.1666666666666667*f[2]*dv1^2 
   out[3] = 2.0*f[1]*w1*w2+0.5773502691896258*f[3]*dv1*w2+0.5773502691896258*f[4]*dv2*w1+0.1666666666666667*f[7]*dv1*dv2 
   out[4] = 2.0*f[2]*w1*w2+0.5773502691896258*f[5]*dv1*w2+0.5773502691896258*f[6]*dv2*w1+0.1666666666666667*f[8]*dv1*dv2 
   out[5] = 2.0*f[1]*w2^2+1.154700538379252*f[4]*dv2*w2+0.1666666666666667*f[1]*dv2^2 
   out[6] = 2.0*f[2]*w2^2+1.154700538379252*f[6]*dv2*w2+0.1666666666666667*f[2]*dv2^2 
end 
_M.energy[1] = function (f, out, dv, w) 
   local dv1, dv2 = dv[1], dv[2] 
   local w1, w2 = w[1], w[2] 
   out[1] = f[1]*w2^2+0.5773502691896258*f[4]*dv2*w2+f[1]*w1^2+0.5773502691896258*f[3]*dv1*w1+0.08333333333333333*f[1]*dv2^2+0.08333333333333333*f[1]*dv1^2 
   out[2] = f[2]*w2^2+0.5773502691896258*f[6]*dv2*w2+f[2]*w1^2+0.5773502691896258*f[5]*dv1*w1+0.08333333333333333*f[2]*dv2^2+0.08333333333333333*f[2]*dv1^2 
end 
_M.pressureTensor[2] = function (f, out, dv, w) 
   local dv1, dv2 = dv[1], dv[2] 
   local w1, w2 = w[1], w[2] 
   out[1] = 2.0*f[1]*w1^2+1.154700538379252*f[3]*dv1*w1+0.149071198499986*f[9]*dv1^2+0.1666666666666667*f[1]*dv1^2 
   out[2] = 2.0*f[2]*w1^2+1.154700538379252*f[5]*dv1*w1+0.149071198499986*f[13]*dv1^2+0.1666666666666667*f[2]*dv1^2 
   out[3] = 2.0*f[8]*w1^2+1.154700538379252*f[12]*dv1*w1+0.1666666666666667*f[8]*dv1^2 
   out[4] = 2.0*f[1]*w1*w2+0.5773502691896258*f[3]*dv1*w2+0.5773502691896258*f[4]*dv2*w1+0.1666666666666667*f[7]*dv1*dv2 
   out[5] = 2.0*f[2]*w1*w2+0.5773502691896258*f[5]*dv1*w2+0.5773502691896258*f[6]*dv2*w1+0.1666666666666667*f[11]*dv1*dv2 
   out[6] = 2.0*f[8]*w1*w2+0.5773502691896258*f[12]*dv1*w2+0.5773502691896258*f[14]*dv2*w1+0.1666666666666667*f[18]*dv1*dv2 
   out[7] = 2.0*f[1]*w2^2+1.154700538379252*f[4]*dv2*w2+0.149071198499986*f[10]*dv2^2+0.1666666666666667*f[1]*dv2^2 
   out[8] = 2.0*f[2]*w2^2+1.154700538379252*f[6]*dv2*w2+0.149071198499986*f[16]*dv2^2+0.1666666666666667*f[2]*dv2^2 
   out[9] = 2.0*f[8]*w2^2+1.154700538379252*f[14]*dv2*w2+0.1666666666666667*f[8]*dv2^2 
end 
_M.energy[2] = function (f, out, dv, w) 
   local dv1, dv2 = dv[1], dv[2] 
   local w1, w2 = w[1], w[2] 
   out[1] = f[1]*w2^2+0.5773502691896258*f[4]*dv2*w2+f[1]*w1^2+0.5773502691896258*f[3]*dv1*w1+0.07453559924999298*f[10]*dv2^2+0.08333333333333333*f[1]*dv2^2+0.07453559924999298*f[9]*dv1^2+0.08333333333333333*f[1]*dv1^2 
   out[2] = f[2]*w2^2+0.5773502691896258*f[6]*dv2*w2+f[2]*w1^2+0.5773502691896258*f[5]*dv1*w1+0.07453559924999298*f[16]*dv2^2+0.08333333333333333*f[2]*dv2^2+0.07453559924999298*f[13]*dv1^2+0.08333333333333333*f[2]*dv1^2 
   out[3] = f[8]*w2^2+0.5773502691896258*f[14]*dv2*w2+f[8]*w1^2+0.5773502691896258*f[12]*dv1*w1+0.08333333333333333*f[8]*dv2^2+0.08333333333333333*f[8]*dv1^2 
end 
_M.pressureTensor[3] = function (f, out, dv, w) 
   local dv1, dv2 = dv[1], dv[2] 
   local w1, w2 = w[1], w[2] 
   out[1] = 2.0*f[1]*w1^2+1.154700538379252*f[3]*dv1*w1+0.149071198499986*f[9]*dv1^2+0.1666666666666667*f[1]*dv1^2 
   out[2] = 2.0*f[2]*w1^2+1.154700538379252*f[5]*dv1*w1+0.149071198499986*f[13]*dv1^2+0.1666666666666667*f[2]*dv1^2 
   out[3] = 2.0*f[8]*w1^2+1.154700538379252*f[12]*dv1*w1+0.1666666666666667*f[8]*dv1^2 
   out[4] = 2.0*f[18]*w1^2+1.154700538379252*f[24]*dv1*w1+0.1666666666666667*f[18]*dv1^2 
   out[5] = 2.0*f[1]*w1*w2+0.5773502691896258*f[3]*dv1*w2+0.5773502691896258*f[4]*dv2*w1+0.1666666666666667*f[7]*dv1*dv2 
   out[6] = 2.0*f[2]*w1*w2+0.5773502691896258*f[5]*dv1*w2+0.5773502691896258*f[6]*dv2*w1+0.1666666666666667*f[11]*dv1*dv2 
   out[7] = 2.0*f[8]*w1*w2+0.5773502691896258*f[12]*dv1*w2+0.5773502691896258*f[14]*dv2*w1+0.1666666666666667*f[21]*dv1*dv2 
   out[8] = 2.0*f[18]*w1*w2+0.5773502691896258*f[24]*dv1*w2+0.5773502691896258*f[26]*dv2*w1+0.1666666666666667*f[30]*dv1*dv2 
   out[9] = 2.0*f[1]*w2^2+1.154700538379252*f[4]*dv2*w2+0.149071198499986*f[10]*dv2^2+0.1666666666666667*f[1]*dv2^2 
   out[10] = 2.0*f[2]*w2^2+1.154700538379252*f[6]*dv2*w2+0.149071198499986*f[16]*dv2^2+0.1666666666666667*f[2]*dv2^2 
   out[11] = 2.0*f[8]*w2^2+1.154700538379252*f[14]*dv2*w2+0.1666666666666667*f[8]*dv2^2 
   out[12] = 2.0*f[18]*w2^2+1.154700538379252*f[26]*dv2*w2+0.1666666666666667*f[18]*dv2^2 
end 
_M.energy[3] = function (f, out, dv, w) 
   local dv1, dv2 = dv[1], dv[2] 
   local w1, w2 = w[1], w[2] 
   out[1] = f[1]*w2^2+0.5773502691896258*f[4]*dv2*w2+f[1]*w1^2+0.5773502691896258*f[3]*dv1*w1+0.07453559924999298*f[10]*dv2^2+0.08333333333333333*f[1]*dv2^2+0.07453559924999298*f[9]*dv1^2+0.08333333333333333*f[1]*dv1^2 
   out[2] = f[2]*w2^2+0.5773502691896258*f[6]*dv2*w2+f[2]*w1^2+0.5773502691896258*f[5]*dv1*w1+0.07453559924999298*f[16]*dv2^2+0.08333333333333333*f[2]*dv2^2+0.07453559924999298*f[13]*dv1^2+0.08333333333333333*f[2]*dv1^2 
   out[3] = f[8]*w2^2+0.5773502691896258*f[14]*dv2*w2+f[8]*w1^2+0.5773502691896258*f[12]*dv1*w1+0.08333333333333333*f[8]*dv2^2+0.08333333333333333*f[8]*dv1^2 
   out[4] = f[18]*w2^2+0.5773502691896258*f[26]*dv2*w2+f[18]*w1^2+0.5773502691896258*f[24]*dv1*w1+0.08333333333333333*f[18]*dv2^2+0.08333333333333333*f[18]*dv1^2 
end 
_M.pressureTensor[4] = function (f, out, dv, w) 
   local dv1, dv2 = dv[1], dv[2] 
   local w1, w2 = w[1], w[2] 
   out[1] = 2.0*f[1]*w1^2+1.154700538379252*f[3]*dv1*w1+0.149071198499986*f[9]*dv1^2+0.1666666666666667*f[1]*dv1^2 
   out[2] = 2.0*f[2]*w1^2+1.154700538379252*f[5]*dv1*w1+0.149071198499986*f[13]*dv1^2+0.1666666666666667*f[2]*dv1^2 
   out[3] = 2.0*f[8]*w1^2+1.154700538379252*f[12]*dv1*w1+0.149071198499986*f[24]*dv1^2+0.1666666666666667*f[8]*dv1^2 
   out[4] = 2.0*f[18]*w1^2+1.154700538379252*f[27]*dv1*w1+0.1666666666666667*f[18]*dv1^2 
   out[5] = 2.0*f[33]*w1^2+1.154700538379252*f[42]*dv1*w1+0.1666666666666667*f[33]*dv1^2 
   out[6] = 2.0*f[1]*w1*w2+0.5773502691896258*f[3]*dv1*w2+0.5773502691896258*f[4]*dv2*w1+0.1666666666666667*f[7]*dv1*dv2 
   out[7] = 2.0*f[2]*w1*w2+0.5773502691896258*f[5]*dv1*w2+0.5773502691896258*f[6]*dv2*w1+0.1666666666666667*f[11]*dv1*dv2 
   out[8] = 2.0*f[8]*w1*w2+0.5773502691896258*f[12]*dv1*w2+0.5773502691896258*f[14]*dv2*w1+0.1666666666666667*f[21]*dv1*dv2 
   out[9] = 2.0*f[18]*w1*w2+0.5773502691896258*f[27]*dv1*w2+0.5773502691896258*f[29]*dv2*w1+0.1666666666666667*f[39]*dv1*dv2 
   out[10] = 2.0*f[33]*w1*w2+0.5773502691896258*f[42]*dv1*w2+0.5773502691896258*f[44]*dv2*w1+0.1666666666666667*f[48]*dv1*dv2 
   out[11] = 2.0*f[1]*w2^2+1.154700538379252*f[4]*dv2*w2+0.149071198499986*f[10]*dv2^2+0.1666666666666667*f[1]*dv2^2 
   out[12] = 2.0*f[2]*w2^2+1.154700538379252*f[6]*dv2*w2+0.149071198499986*f[16]*dv2^2+0.1666666666666667*f[2]*dv2^2 
   out[13] = 2.0*f[8]*w2^2+1.154700538379252*f[14]*dv2*w2+0.149071198499986*f[25]*dv2^2+0.1666666666666667*f[8]*dv2^2 
   out[14] = 2.0*f[18]*w2^2+1.154700538379252*f[29]*dv2*w2+0.1666666666666667*f[18]*dv2^2 
   out[15] = 2.0*f[33]*w2^2+1.154700538379252*f[44]*dv2*w2+0.1666666666666667*f[33]*dv2^2 
end 
_M.energy[4] = function (f, out, dv, w) 
   local dv1, dv2 = dv[1], dv[2] 
   local w1, w2 = w[1], w[2] 
   out[1] = f[1]*w2^2+0.5773502691896258*f[4]*dv2*w2+f[1]*w1^2+0.5773502691896258*f[3]*dv1*w1+0.07453559924999298*f[10]*dv2^2+0.08333333333333333*f[1]*dv2^2+0.07453559924999298*f[9]*dv1^2+0.08333333333333333*f[1]*dv1^2 
   out[2] = f[2]*w2^2+0.5773502691896258*f[6]*dv2*w2+f[2]*w1^2+0.5773502691896258*f[5]*dv1*w1+0.07453559924999298*f[16]*dv2^2+0.08333333333333333*f[2]*dv2^2+0.07453559924999298*f[13]*dv1^2+0.08333333333333333*f[2]*dv1^2 
   out[3] = f[8]*w2^2+0.5773502691896258*f[14]*dv2*w2+f[8]*w1^2+0.5773502691896258*f[12]*dv1*w1+0.07453559924999298*f[25]*dv2^2+0.08333333333333333*f[8]*dv2^2+0.07453559924999298*f[24]*dv1^2+0.08333333333333333*f[8]*dv1^2 
   out[4] = f[18]*w2^2+0.5773502691896258*f[29]*dv2*w2+f[18]*w1^2+0.5773502691896258*f[27]*dv1*w1+0.08333333333333333*f[18]*dv2^2+0.08333333333333333*f[18]*dv1^2 
   out[5] = f[33]*w2^2+0.5773502691896258*f[44]*dv2*w2+f[33]*w1^2+0.5773502691896258*f[42]*dv1*w1+0.08333333333333333*f[33]*dv2^2+0.08333333333333333*f[33]*dv1^2 
end 
return _M 
