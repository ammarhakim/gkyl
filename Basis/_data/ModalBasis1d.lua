local _M = {} 
_M[1] = function (z, b) 
   local z1 = z[1] 
   b[1] = 0.7071067811865475 
   b[2] = 1.224744871391589*z1 
end 
_M[2] = function (z, b) 
   local z1 = z[1] 
   b[1] = 0.7071067811865475 
   b[2] = 1.224744871391589*z1 
   b[3] = 2.371708245126284*z1^2-0.7905694150420947 
end 
_M[3] = function (z, b) 
   local z1 = z[1] 
   b[1] = 0.7071067811865475 
   b[2] = 1.224744871391589*z1 
   b[3] = 2.371708245126284*z1^2-0.7905694150420947 
   b[4] = 4.677071733467426*z1^3-2.806243040080456*z1 
end 
return _M 
