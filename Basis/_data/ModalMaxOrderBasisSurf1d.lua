local _M = {} 
-- polyOrder 1 
_M[1] = {} 
--    dir 1 
_M[1][1] = {} 
_M[1][1].upper = function (volIn, surfOut) 
   surfOut[1] = 1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
end 
_M[1][1].lower = function (volIn, surfOut) 
   surfOut[1] = 0.7071067811865475*volIn[1]-1.224744871391589*volIn[2] 
end 


-- polyOrder 2 
_M[2] = {} 
--    dir 1 
_M[2][1] = {} 
_M[2][1].upper = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[3]+1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
end 
_M[2][1].lower = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[3]-1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
end 


-- polyOrder 3 
_M[3] = {} 
--    dir 1 
_M[3][1] = {} 
_M[3][1].upper = function (volIn, surfOut) 
   surfOut[1] = 1.870828693386971*volIn[4]+1.58113883008419*volIn[3]+1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
end 
_M[3][1].lower = function (volIn, surfOut) 
   surfOut[1] = (-1.870828693386971*volIn[4])+1.58113883008419*volIn[3]-1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
end 


-- polyOrder 4 
_M[4] = {} 
--    dir 1 
_M[4][1] = {} 
_M[4][1].upper = function (volIn, surfOut) 
   surfOut[1] = 2.121320343559642*volIn[5]+1.870828693386971*volIn[4]+1.58113883008419*volIn[3]+1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
end 
_M[4][1].lower = function (volIn, surfOut) 
   surfOut[1] = 2.121320343559642*volIn[5]-1.870828693386971*volIn[4]+1.58113883008419*volIn[3]-1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
end 


return _M 
