local _M = {} 
-- polyOrder 1 
_M[1] = {} 
--    dir 1 
_M[1][1] = {} 
_M[1][1].upper = function (volIn, surfOut) 
   surfOut[1] = 1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[3] 
end 
_M[1][1].lower = function (volIn, surfOut) 
   surfOut[1] = 0.7071067811865475*volIn[1]-1.224744871391589*volIn[2] 
   surfOut[2] = 0.7071067811865475*volIn[3] 
end 

--    dir 2 
_M[1][2] = {} 
_M[1][2].upper = function (volIn, surfOut) 
   surfOut[1] = 1.224744871391589*volIn[3]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[2] 
end 
_M[1][2].lower = function (volIn, surfOut) 
   surfOut[1] = 0.7071067811865475*volIn[1]-1.224744871391589*volIn[3] 
   surfOut[2] = 0.7071067811865475*volIn[2] 
end 


-- polyOrder 2 
_M[2] = {} 
--    dir 1 
_M[2][1] = {} 
_M[2][1].upper = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[5]+1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.224744871391589*volIn[4]+0.7071067811865475*volIn[3] 
   surfOut[3] = 0.7071067811865475*volIn[6] 
end 
_M[2][1].lower = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[5]-1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[3]-1.224744871391589*volIn[4] 
   surfOut[3] = 0.7071067811865475*volIn[6] 
end 

--    dir 2 
_M[2][2] = {} 
_M[2][2].upper = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[6]+1.224744871391589*volIn[3]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.224744871391589*volIn[4]+0.7071067811865475*volIn[2] 
   surfOut[3] = 0.7071067811865475*volIn[5] 
end 
_M[2][2].lower = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[6]-1.224744871391589*volIn[3]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[2]-1.224744871391589*volIn[4] 
   surfOut[3] = 0.7071067811865475*volIn[5] 
end 


-- polyOrder 3 
_M[3] = {} 
--    dir 1 
_M[3][1] = {} 
_M[3][1].upper = function (volIn, surfOut) 
   surfOut[1] = 1.870828693386971*volIn[9]+1.58113883008419*volIn[5]+1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.58113883008419*volIn[7]+1.224744871391589*volIn[4]+0.7071067811865475*volIn[3] 
   surfOut[3] = 1.224744871391589*volIn[8]+0.7071067811865475*volIn[6] 
   surfOut[4] = 0.7071067811865475*volIn[10] 
end 
_M[3][1].lower = function (volIn, surfOut) 
   surfOut[1] = (-1.870828693386971*volIn[9])+1.58113883008419*volIn[5]-1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.58113883008419*volIn[7]-1.224744871391589*volIn[4]+0.7071067811865475*volIn[3] 
   surfOut[3] = 0.7071067811865475*volIn[6]-1.224744871391589*volIn[8] 
   surfOut[4] = 0.7071067811865475*volIn[10] 
end 

--    dir 2 
_M[3][2] = {} 
_M[3][2].upper = function (volIn, surfOut) 
   surfOut[1] = 1.870828693386971*volIn[10]+1.58113883008419*volIn[6]+1.224744871391589*volIn[3]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.58113883008419*volIn[8]+1.224744871391589*volIn[4]+0.7071067811865475*volIn[2] 
   surfOut[3] = 1.224744871391589*volIn[7]+0.7071067811865475*volIn[5] 
   surfOut[4] = 0.7071067811865475*volIn[9] 
end 
_M[3][2].lower = function (volIn, surfOut) 
   surfOut[1] = (-1.870828693386971*volIn[10])+1.58113883008419*volIn[6]-1.224744871391589*volIn[3]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.58113883008419*volIn[8]-1.224744871391589*volIn[4]+0.7071067811865475*volIn[2] 
   surfOut[3] = 0.7071067811865475*volIn[5]-1.224744871391589*volIn[7] 
   surfOut[4] = 0.7071067811865475*volIn[9] 
end 


-- polyOrder 4 
_M[4] = {} 
--    dir 1 
_M[4][1] = {} 
_M[4][1].upper = function (volIn, surfOut) 
   surfOut[1] = 2.121320343559642*volIn[14]+1.870828693386971*volIn[9]+1.58113883008419*volIn[5]+1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.870828693386971*volIn[12]+1.58113883008419*volIn[7]+1.224744871391589*volIn[4]+0.7071067811865475*volIn[3] 
   surfOut[3] = 1.58113883008419*volIn[11]+1.224744871391589*volIn[8]+0.7071067811865475*volIn[6] 
   surfOut[4] = 1.224744871391589*volIn[13]+0.7071067811865475*volIn[10] 
   surfOut[5] = 0.7071067811865475*volIn[15] 
end 
_M[4][1].lower = function (volIn, surfOut) 
   surfOut[1] = 2.121320343559642*volIn[14]-1.870828693386971*volIn[9]+1.58113883008419*volIn[5]-1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
   surfOut[2] = (-1.870828693386971*volIn[12])+1.58113883008419*volIn[7]-1.224744871391589*volIn[4]+0.7071067811865475*volIn[3] 
   surfOut[3] = 1.58113883008419*volIn[11]-1.224744871391589*volIn[8]+0.7071067811865475*volIn[6] 
   surfOut[4] = 0.7071067811865475*volIn[10]-1.224744871391589*volIn[13] 
   surfOut[5] = 0.7071067811865475*volIn[15] 
end 

--    dir 2 
_M[4][2] = {} 
_M[4][2].upper = function (volIn, surfOut) 
   surfOut[1] = 2.121320343559642*volIn[15]+1.870828693386971*volIn[10]+1.58113883008419*volIn[6]+1.224744871391589*volIn[3]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.870828693386971*volIn[13]+1.58113883008419*volIn[8]+1.224744871391589*volIn[4]+0.7071067811865475*volIn[2] 
   surfOut[3] = 1.58113883008419*volIn[11]+1.224744871391589*volIn[7]+0.7071067811865475*volIn[5] 
   surfOut[4] = 1.224744871391589*volIn[12]+0.7071067811865475*volIn[9] 
   surfOut[5] = 0.7071067811865475*volIn[14] 
end 
_M[4][2].lower = function (volIn, surfOut) 
   surfOut[1] = 2.121320343559642*volIn[15]-1.870828693386971*volIn[10]+1.58113883008419*volIn[6]-1.224744871391589*volIn[3]+0.7071067811865475*volIn[1] 
   surfOut[2] = (-1.870828693386971*volIn[13])+1.58113883008419*volIn[8]-1.224744871391589*volIn[4]+0.7071067811865475*volIn[2] 
   surfOut[3] = 1.58113883008419*volIn[11]-1.224744871391589*volIn[7]+0.7071067811865475*volIn[5] 
   surfOut[4] = 0.7071067811865475*volIn[9]-1.224744871391589*volIn[12] 
   surfOut[5] = 0.7071067811865475*volIn[14] 
end 


return _M 
