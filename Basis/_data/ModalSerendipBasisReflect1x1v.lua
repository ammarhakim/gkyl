local _M = {} 
_M[1] = function (dir, f, out) 
   if dir == 1  then 
   out[1] = f[1] 
   out[2] = -1.0*f[2] 
   out[3] = -1.0*f[3] 
   out[4] = f[4] 
   end 
end 
_M[2] = function (dir, f, out) 
   if dir == 1  then 
   out[1] = f[1] 
   out[2] = -1.0*f[2] 
   out[3] = -1.0*f[3] 
   out[4] = f[4] 
   out[5] = f[5] 
   out[6] = f[6] 
   out[7] = -1.0*f[7] 
   out[8] = -1.0*f[8] 
   end 
end 
return _M 
