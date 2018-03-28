local _M = {} 
_M[1] = function (dir, f, out) 
   if dir == 1  then 
   out[1] = f[1] 
   out[2] = -1.0*f[2] 
   out[3] = -1.0*f[3] 
   out[4] = f[4] 
   out[5] = f[5] 
   out[6] = -1.0*f[6] 
   out[7] = -1.0*f[7] 
   out[8] = f[8] 
   end 
end 
_M[2] = function (dir, f, out) 
   if dir == 1  then 
   out[1] = f[1] 
   out[2] = -1.0*f[2] 
   out[3] = -1.0*f[3] 
   out[4] = f[4] 
   out[5] = f[5] 
   out[6] = -1.0*f[6] 
   out[7] = -1.0*f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = f[11] 
   out[12] = -1.0*f[12] 
   out[13] = -1.0*f[13] 
   out[14] = f[14] 
   out[15] = f[15] 
   out[16] = -1.0*f[16] 
   out[17] = -1.0*f[17] 
   out[18] = -1.0*f[18] 
   out[19] = -1.0*f[19] 
   out[20] = f[20] 
   end 
end 
return _M 
