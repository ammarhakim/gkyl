local _M = { pressureTensor = {}, energy = {} } 
 
_M.pressureTensor[1] = function (f, out, dv, w) 
   local dv1 = dv[1] 
   local w1 = w[1] 
   out[1] = 1.414213562373095*f[1]*w1^2+0.8164965809277261*f[3]*dv1*w1+0.1178511301977579*f[1]*dv1^2 
   out[2] = 1.414213562373095*f[2]*w1^2+0.1178511301977579*f[2]*dv1^2 
end 
_M.energy[1] = function (f, out, dv, w) 
   local dv1 = dv[1] 
   local w1 = w[1] 
   out[1] = 0.7071067811865475*f[1]*w1^2+0.408248290463863*f[3]*dv1*w1+0.05892556509887893*f[1]*dv1^2 
   out[2] = 0.7071067811865475*f[2]*w1^2+0.05892556509887893*f[2]*dv1^2 
end 
_M.pressureTensor[2] = function (f, out, dv, w) 
   local dv1 = dv[1] 
   local w1 = w[1] 
   out[1] = 1.414213562373095*f[1]*w1^2+0.8164965809277261*f[3]*dv1*w1+0.105409255338946*f[6]*dv1^2+0.1178511301977579*f[1]*dv1^2 
   out[2] = 1.414213562373095*f[2]*w1^2+0.8164965809277261*f[4]*dv1*w1+0.1178511301977579*f[2]*dv1^2 
   out[3] = 1.414213562373095*f[5]*w1^2+0.1178511301977579*f[5]*dv1^2 
end 
_M.energy[2] = function (f, out, dv, w) 
   local dv1 = dv[1] 
   local w1 = w[1] 
   out[1] = 0.7071067811865475*f[1]*w1^2+0.408248290463863*f[3]*dv1*w1+0.05270462766947297*f[6]*dv1^2+0.05892556509887893*f[1]*dv1^2 
   out[2] = 0.7071067811865475*f[2]*w1^2+0.408248290463863*f[4]*dv1*w1+0.05892556509887893*f[2]*dv1^2 
   out[3] = 0.7071067811865475*f[5]*w1^2+0.05892556509887893*f[5]*dv1^2 
end 
_M.pressureTensor[3] = function (f, out, dv, w) 
   local dv1 = dv[1] 
   local w1 = w[1] 
   out[1] = 1.414213562373095*f[1]*w1^2+0.8164965809277261*f[3]*dv1*w1+0.105409255338946*f[6]*dv1^2+0.1178511301977579*f[1]*dv1^2 
   out[2] = 1.414213562373095*f[2]*w1^2+0.8164965809277261*f[4]*dv1*w1+0.105409255338946*f[8]*dv1^2+0.1178511301977579*f[2]*dv1^2 
   out[3] = 1.414213562373095*f[5]*w1^2+0.8164965809277261*f[7]*dv1*w1+0.1178511301977579*f[5]*dv1^2 
   out[4] = 1.414213562373095*f[9]*w1^2+0.1178511301977579*f[9]*dv1^2 
end 
_M.energy[3] = function (f, out, dv, w) 
   local dv1 = dv[1] 
   local w1 = w[1] 
   out[1] = 0.7071067811865475*f[1]*w1^2+0.408248290463863*f[3]*dv1*w1+0.05270462766947297*f[6]*dv1^2+0.05892556509887893*f[1]*dv1^2 
   out[2] = 0.7071067811865475*f[2]*w1^2+0.408248290463863*f[4]*dv1*w1+0.05270462766947297*f[8]*dv1^2+0.05892556509887893*f[2]*dv1^2 
   out[3] = 0.7071067811865475*f[5]*w1^2+0.408248290463863*f[7]*dv1*w1+0.05892556509887893*f[5]*dv1^2 
   out[4] = 0.7071067811865475*f[9]*w1^2+0.05892556509887893*f[9]*dv1^2 
end 
_M.pressureTensor[4] = function (f, out, dv, w) 
   local dv1 = dv[1] 
   local w1 = w[1] 
   out[1] = 1.414213562373095*f[1]*w1^2+0.8164965809277261*f[3]*dv1*w1+0.105409255338946*f[6]*dv1^2+0.1178511301977579*f[1]*dv1^2 
   out[2] = 1.414213562373095*f[2]*w1^2+0.8164965809277261*f[4]*dv1*w1+0.105409255338946*f[8]*dv1^2+0.1178511301977579*f[2]*dv1^2 
   out[3] = 1.414213562373095*f[5]*w1^2+0.8164965809277261*f[7]*dv1*w1+0.105409255338946*f[11]*dv1^2+0.1178511301977579*f[5]*dv1^2 
   out[4] = 1.414213562373095*f[9]*w1^2+0.8164965809277261*f[12]*dv1*w1+0.1178511301977579*f[9]*dv1^2 
   out[5] = 1.414213562373095*f[14]*w1^2+0.1178511301977579*f[14]*dv1^2 
end 
_M.energy[4] = function (f, out, dv, w) 
   local dv1 = dv[1] 
   local w1 = w[1] 
   out[1] = 0.7071067811865475*f[1]*w1^2+0.408248290463863*f[3]*dv1*w1+0.05270462766947297*f[6]*dv1^2+0.05892556509887893*f[1]*dv1^2 
   out[2] = 0.7071067811865475*f[2]*w1^2+0.408248290463863*f[4]*dv1*w1+0.05270462766947297*f[8]*dv1^2+0.05892556509887893*f[2]*dv1^2 
   out[3] = 0.7071067811865475*f[5]*w1^2+0.408248290463863*f[7]*dv1*w1+0.05270462766947297*f[11]*dv1^2+0.05892556509887893*f[5]*dv1^2 
   out[4] = 0.7071067811865475*f[9]*w1^2+0.408248290463863*f[12]*dv1*w1+0.05892556509887893*f[9]*dv1^2 
   out[5] = 0.7071067811865475*f[14]*w1^2+0.05892556509887893*f[14]*dv1^2 
end 
return _M 
