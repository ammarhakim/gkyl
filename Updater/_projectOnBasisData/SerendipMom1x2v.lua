local _M = { numDensity = {}, momentum = {} } 
 
_M.numDensity[1] = function (f, out, dv, w) 
   out[1] = 2.0*f[1] 
   out[2] = 2.0*f[2] 
end 
_M.numDensity[2] = function (f, out, dv, w) 
   out[1] = 2.0*f[1] 
   out[2] = 2.0*f[2] 
   out[3] = 2.0*f[8] 
end 
_M.numDensity[3] = function (f, out, dv, w) 
   out[1] = 2.0*f[1] 
   out[2] = 2.0*f[2] 
   out[3] = 2.0*f[8] 
   out[4] = 2.0*f[18] 
end 
_M.numDensity[4] = function (f, out, dv, w) 
   out[1] = 2.0*f[1] 
   out[2] = 2.0*f[2] 
   out[3] = 2.0*f[8] 
   out[4] = 2.0*f[18] 
   out[5] = 2.0*f[33] 
end 

 
_M.momentum[1] = function (f, out, dv, w) 
   local dv1, dv2 = dv[1], dv[2] 
   local w1, w2 = w[1], w[2] 
   out[1] = 2.0*f[1]*w1+0.5773502691896258*f[3]*dv1 
   out[2] = 2.0*f[2]*w1+0.5773502691896258*f[5]*dv1 
   out[3] = 2.0*f[1]*w2+0.5773502691896258*f[4]*dv2 
   out[4] = 2.0*f[2]*w2+0.5773502691896258*f[6]*dv2 
end 
_M.momentum[2] = function (f, out, dv, w) 
   local dv1, dv2 = dv[1], dv[2] 
   local w1, w2 = w[1], w[2] 
   out[1] = 2.0*f[1]*w1+0.5773502691896258*f[3]*dv1 
   out[2] = 2.0*f[2]*w1+0.5773502691896258*f[5]*dv1 
   out[3] = 2.0*f[8]*w1+0.5773502691896258*f[12]*dv1 
   out[4] = 2.0*f[1]*w2+0.5773502691896258*f[4]*dv2 
   out[5] = 2.0*f[2]*w2+0.5773502691896258*f[6]*dv2 
   out[6] = 2.0*f[8]*w2+0.5773502691896258*f[14]*dv2 
end 
_M.momentum[3] = function (f, out, dv, w) 
   local dv1, dv2 = dv[1], dv[2] 
   local w1, w2 = w[1], w[2] 
   out[1] = 2.0*f[1]*w1+0.5773502691896258*f[3]*dv1 
   out[2] = 2.0*f[2]*w1+0.5773502691896258*f[5]*dv1 
   out[3] = 2.0*f[8]*w1+0.5773502691896258*f[12]*dv1 
   out[4] = 2.0*f[18]*w1+0.5773502691896258*f[24]*dv1 
   out[5] = 2.0*f[1]*w2+0.5773502691896258*f[4]*dv2 
   out[6] = 2.0*f[2]*w2+0.5773502691896258*f[6]*dv2 
   out[7] = 2.0*f[8]*w2+0.5773502691896258*f[14]*dv2 
   out[8] = 2.0*f[18]*w2+0.5773502691896258*f[26]*dv2 
end 
_M.momentum[4] = function (f, out, dv, w) 
   local dv1, dv2 = dv[1], dv[2] 
   local w1, w2 = w[1], w[2] 
   out[1] = 2.0*f[1]*w1+0.5773502691896258*f[3]*dv1 
   out[2] = 2.0*f[2]*w1+0.5773502691896258*f[5]*dv1 
   out[3] = 2.0*f[8]*w1+0.5773502691896258*f[12]*dv1 
   out[4] = 2.0*f[18]*w1+0.5773502691896258*f[27]*dv1 
   out[5] = 2.0*f[33]*w1+0.5773502691896258*f[42]*dv1 
   out[6] = 2.0*f[1]*w2+0.5773502691896258*f[4]*dv2 
   out[7] = 2.0*f[2]*w2+0.5773502691896258*f[6]*dv2 
   out[8] = 2.0*f[8]*w2+0.5773502691896258*f[14]*dv2 
   out[9] = 2.0*f[18]*w2+0.5773502691896258*f[29]*dv2 
   out[10] = 2.0*f[33]*w2+0.5773502691896258*f[44]*dv2 
end 
return _M 
