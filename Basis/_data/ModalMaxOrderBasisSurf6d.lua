local _M = {} 
-- polyOrder 1 
_M[1] = {} 
--    dir 1 
_M[1][1] = {} 
_M[1][1].upper = function (volIn, surfOut) 
   surfOut[1] = 1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[3] 
   surfOut[3] = 0.7071067811865475*volIn[4] 
   surfOut[4] = 0.7071067811865475*volIn[5] 
   surfOut[5] = 0.7071067811865475*volIn[6] 
   surfOut[6] = 0.7071067811865475*volIn[7] 
end 
_M[1][1].lower = function (volIn, surfOut) 
   surfOut[1] = 0.7071067811865475*volIn[1]-1.224744871391589*volIn[2] 
   surfOut[2] = 0.7071067811865475*volIn[3] 
   surfOut[3] = 0.7071067811865475*volIn[4] 
   surfOut[4] = 0.7071067811865475*volIn[5] 
   surfOut[5] = 0.7071067811865475*volIn[6] 
   surfOut[6] = 0.7071067811865475*volIn[7] 
end 

--    dir 2 
_M[1][2] = {} 
_M[1][2].upper = function (volIn, surfOut) 
   surfOut[1] = 1.224744871391589*volIn[3]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[2] 
   surfOut[3] = 0.7071067811865475*volIn[4] 
   surfOut[4] = 0.7071067811865475*volIn[5] 
   surfOut[5] = 0.7071067811865475*volIn[6] 
   surfOut[6] = 0.7071067811865475*volIn[7] 
end 
_M[1][2].lower = function (volIn, surfOut) 
   surfOut[1] = 0.7071067811865475*volIn[1]-1.224744871391589*volIn[3] 
   surfOut[2] = 0.7071067811865475*volIn[2] 
   surfOut[3] = 0.7071067811865475*volIn[4] 
   surfOut[4] = 0.7071067811865475*volIn[5] 
   surfOut[5] = 0.7071067811865475*volIn[6] 
   surfOut[6] = 0.7071067811865475*volIn[7] 
end 

--    dir 3 
_M[1][3] = {} 
_M[1][3].upper = function (volIn, surfOut) 
   surfOut[1] = 1.224744871391589*volIn[4]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[2] 
   surfOut[3] = 0.7071067811865475*volIn[3] 
   surfOut[4] = 0.7071067811865475*volIn[5] 
   surfOut[5] = 0.7071067811865475*volIn[6] 
   surfOut[6] = 0.7071067811865475*volIn[7] 
end 
_M[1][3].lower = function (volIn, surfOut) 
   surfOut[1] = 0.7071067811865475*volIn[1]-1.224744871391589*volIn[4] 
   surfOut[2] = 0.7071067811865475*volIn[2] 
   surfOut[3] = 0.7071067811865475*volIn[3] 
   surfOut[4] = 0.7071067811865475*volIn[5] 
   surfOut[5] = 0.7071067811865475*volIn[6] 
   surfOut[6] = 0.7071067811865475*volIn[7] 
end 

--    dir 4 
_M[1][4] = {} 
_M[1][4].upper = function (volIn, surfOut) 
   surfOut[1] = 1.224744871391589*volIn[5]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[2] 
   surfOut[3] = 0.7071067811865475*volIn[3] 
   surfOut[4] = 0.7071067811865475*volIn[4] 
   surfOut[5] = 0.7071067811865475*volIn[6] 
   surfOut[6] = 0.7071067811865475*volIn[7] 
end 
_M[1][4].lower = function (volIn, surfOut) 
   surfOut[1] = 0.7071067811865475*volIn[1]-1.224744871391589*volIn[5] 
   surfOut[2] = 0.7071067811865475*volIn[2] 
   surfOut[3] = 0.7071067811865475*volIn[3] 
   surfOut[4] = 0.7071067811865475*volIn[4] 
   surfOut[5] = 0.7071067811865475*volIn[6] 
   surfOut[6] = 0.7071067811865475*volIn[7] 
end 

--    dir 5 
_M[1][5] = {} 
_M[1][5].upper = function (volIn, surfOut) 
   surfOut[1] = 1.224744871391589*volIn[6]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[2] 
   surfOut[3] = 0.7071067811865475*volIn[3] 
   surfOut[4] = 0.7071067811865475*volIn[4] 
   surfOut[5] = 0.7071067811865475*volIn[5] 
   surfOut[6] = 0.7071067811865475*volIn[7] 
end 
_M[1][5].lower = function (volIn, surfOut) 
   surfOut[1] = 0.7071067811865475*volIn[1]-1.224744871391589*volIn[6] 
   surfOut[2] = 0.7071067811865475*volIn[2] 
   surfOut[3] = 0.7071067811865475*volIn[3] 
   surfOut[4] = 0.7071067811865475*volIn[4] 
   surfOut[5] = 0.7071067811865475*volIn[5] 
   surfOut[6] = 0.7071067811865475*volIn[7] 
end 

--    dir 6 
_M[1][6] = {} 
_M[1][6].upper = function (volIn, surfOut) 
   surfOut[1] = 1.224744871391589*volIn[7]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[2] 
   surfOut[3] = 0.7071067811865475*volIn[3] 
   surfOut[4] = 0.7071067811865475*volIn[4] 
   surfOut[5] = 0.7071067811865475*volIn[5] 
   surfOut[6] = 0.7071067811865475*volIn[6] 
end 
_M[1][6].lower = function (volIn, surfOut) 
   surfOut[1] = 0.7071067811865475*volIn[1]-1.224744871391589*volIn[7] 
   surfOut[2] = 0.7071067811865475*volIn[2] 
   surfOut[3] = 0.7071067811865475*volIn[3] 
   surfOut[4] = 0.7071067811865475*volIn[4] 
   surfOut[5] = 0.7071067811865475*volIn[5] 
   surfOut[6] = 0.7071067811865475*volIn[6] 
end 


-- polyOrder 2 
_M[2] = {} 
--    dir 1 
_M[2][1] = {} 
_M[2][1].upper = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[23]+1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.224744871391589*volIn[8]+0.7071067811865475*volIn[3] 
   surfOut[3] = 1.224744871391589*volIn[9]+0.7071067811865475*volIn[4] 
   surfOut[4] = 1.224744871391589*volIn[11]+0.7071067811865475*volIn[5] 
   surfOut[5] = 1.224744871391589*volIn[14]+0.7071067811865475*volIn[6] 
   surfOut[6] = 1.224744871391589*volIn[18]+0.7071067811865475*volIn[7] 
   surfOut[7] = 0.7071067811865475*volIn[10] 
   surfOut[8] = 0.7071067811865475*volIn[12] 
   surfOut[9] = 0.7071067811865475*volIn[13] 
   surfOut[10] = 0.7071067811865475*volIn[15] 
   surfOut[11] = 0.7071067811865475*volIn[16] 
   surfOut[12] = 0.7071067811865475*volIn[17] 
   surfOut[13] = 0.7071067811865475*volIn[19] 
   surfOut[14] = 0.7071067811865475*volIn[20] 
   surfOut[15] = 0.7071067811865475*volIn[21] 
   surfOut[16] = 0.7071067811865475*volIn[22] 
   surfOut[17] = 0.7071067811865475*volIn[24] 
   surfOut[18] = 0.7071067811865475*volIn[25] 
   surfOut[19] = 0.7071067811865475*volIn[26] 
   surfOut[20] = 0.7071067811865475*volIn[27] 
   surfOut[21] = 0.7071067811865475*volIn[28] 
end 
_M[2][1].lower = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[23]-1.224744871391589*volIn[2]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[3]-1.224744871391589*volIn[8] 
   surfOut[3] = 0.7071067811865475*volIn[4]-1.224744871391589*volIn[9] 
   surfOut[4] = 0.7071067811865475*volIn[5]-1.224744871391589*volIn[11] 
   surfOut[5] = 0.7071067811865475*volIn[6]-1.224744871391589*volIn[14] 
   surfOut[6] = 0.7071067811865475*volIn[7]-1.224744871391589*volIn[18] 
   surfOut[7] = 0.7071067811865475*volIn[10] 
   surfOut[8] = 0.7071067811865475*volIn[12] 
   surfOut[9] = 0.7071067811865475*volIn[13] 
   surfOut[10] = 0.7071067811865475*volIn[15] 
   surfOut[11] = 0.7071067811865475*volIn[16] 
   surfOut[12] = 0.7071067811865475*volIn[17] 
   surfOut[13] = 0.7071067811865475*volIn[19] 
   surfOut[14] = 0.7071067811865475*volIn[20] 
   surfOut[15] = 0.7071067811865475*volIn[21] 
   surfOut[16] = 0.7071067811865475*volIn[22] 
   surfOut[17] = 0.7071067811865475*volIn[24] 
   surfOut[18] = 0.7071067811865475*volIn[25] 
   surfOut[19] = 0.7071067811865475*volIn[26] 
   surfOut[20] = 0.7071067811865475*volIn[27] 
   surfOut[21] = 0.7071067811865475*volIn[28] 
end 

--    dir 2 
_M[2][2] = {} 
_M[2][2].upper = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[24]+1.224744871391589*volIn[3]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.224744871391589*volIn[8]+0.7071067811865475*volIn[2] 
   surfOut[3] = 1.224744871391589*volIn[10]+0.7071067811865475*volIn[4] 
   surfOut[4] = 1.224744871391589*volIn[12]+0.7071067811865475*volIn[5] 
   surfOut[5] = 1.224744871391589*volIn[15]+0.7071067811865475*volIn[6] 
   surfOut[6] = 1.224744871391589*volIn[19]+0.7071067811865475*volIn[7] 
   surfOut[7] = 0.7071067811865475*volIn[9] 
   surfOut[8] = 0.7071067811865475*volIn[11] 
   surfOut[9] = 0.7071067811865475*volIn[13] 
   surfOut[10] = 0.7071067811865475*volIn[14] 
   surfOut[11] = 0.7071067811865475*volIn[16] 
   surfOut[12] = 0.7071067811865475*volIn[17] 
   surfOut[13] = 0.7071067811865475*volIn[18] 
   surfOut[14] = 0.7071067811865475*volIn[20] 
   surfOut[15] = 0.7071067811865475*volIn[21] 
   surfOut[16] = 0.7071067811865475*volIn[22] 
   surfOut[17] = 0.7071067811865475*volIn[23] 
   surfOut[18] = 0.7071067811865475*volIn[25] 
   surfOut[19] = 0.7071067811865475*volIn[26] 
   surfOut[20] = 0.7071067811865475*volIn[27] 
   surfOut[21] = 0.7071067811865475*volIn[28] 
end 
_M[2][2].lower = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[24]-1.224744871391589*volIn[3]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[2]-1.224744871391589*volIn[8] 
   surfOut[3] = 0.7071067811865475*volIn[4]-1.224744871391589*volIn[10] 
   surfOut[4] = 0.7071067811865475*volIn[5]-1.224744871391589*volIn[12] 
   surfOut[5] = 0.7071067811865475*volIn[6]-1.224744871391589*volIn[15] 
   surfOut[6] = 0.7071067811865475*volIn[7]-1.224744871391589*volIn[19] 
   surfOut[7] = 0.7071067811865475*volIn[9] 
   surfOut[8] = 0.7071067811865475*volIn[11] 
   surfOut[9] = 0.7071067811865475*volIn[13] 
   surfOut[10] = 0.7071067811865475*volIn[14] 
   surfOut[11] = 0.7071067811865475*volIn[16] 
   surfOut[12] = 0.7071067811865475*volIn[17] 
   surfOut[13] = 0.7071067811865475*volIn[18] 
   surfOut[14] = 0.7071067811865475*volIn[20] 
   surfOut[15] = 0.7071067811865475*volIn[21] 
   surfOut[16] = 0.7071067811865475*volIn[22] 
   surfOut[17] = 0.7071067811865475*volIn[23] 
   surfOut[18] = 0.7071067811865475*volIn[25] 
   surfOut[19] = 0.7071067811865475*volIn[26] 
   surfOut[20] = 0.7071067811865475*volIn[27] 
   surfOut[21] = 0.7071067811865475*volIn[28] 
end 

--    dir 3 
_M[2][3] = {} 
_M[2][3].upper = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[25]+1.224744871391589*volIn[4]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.224744871391589*volIn[9]+0.7071067811865475*volIn[2] 
   surfOut[3] = 1.224744871391589*volIn[10]+0.7071067811865475*volIn[3] 
   surfOut[4] = 1.224744871391589*volIn[13]+0.7071067811865475*volIn[5] 
   surfOut[5] = 1.224744871391589*volIn[16]+0.7071067811865475*volIn[6] 
   surfOut[6] = 1.224744871391589*volIn[20]+0.7071067811865475*volIn[7] 
   surfOut[7] = 0.7071067811865475*volIn[8] 
   surfOut[8] = 0.7071067811865475*volIn[11] 
   surfOut[9] = 0.7071067811865475*volIn[12] 
   surfOut[10] = 0.7071067811865475*volIn[14] 
   surfOut[11] = 0.7071067811865475*volIn[15] 
   surfOut[12] = 0.7071067811865475*volIn[17] 
   surfOut[13] = 0.7071067811865475*volIn[18] 
   surfOut[14] = 0.7071067811865475*volIn[19] 
   surfOut[15] = 0.7071067811865475*volIn[21] 
   surfOut[16] = 0.7071067811865475*volIn[22] 
   surfOut[17] = 0.7071067811865475*volIn[23] 
   surfOut[18] = 0.7071067811865475*volIn[24] 
   surfOut[19] = 0.7071067811865475*volIn[26] 
   surfOut[20] = 0.7071067811865475*volIn[27] 
   surfOut[21] = 0.7071067811865475*volIn[28] 
end 
_M[2][3].lower = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[25]-1.224744871391589*volIn[4]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[2]-1.224744871391589*volIn[9] 
   surfOut[3] = 0.7071067811865475*volIn[3]-1.224744871391589*volIn[10] 
   surfOut[4] = 0.7071067811865475*volIn[5]-1.224744871391589*volIn[13] 
   surfOut[5] = 0.7071067811865475*volIn[6]-1.224744871391589*volIn[16] 
   surfOut[6] = 0.7071067811865475*volIn[7]-1.224744871391589*volIn[20] 
   surfOut[7] = 0.7071067811865475*volIn[8] 
   surfOut[8] = 0.7071067811865475*volIn[11] 
   surfOut[9] = 0.7071067811865475*volIn[12] 
   surfOut[10] = 0.7071067811865475*volIn[14] 
   surfOut[11] = 0.7071067811865475*volIn[15] 
   surfOut[12] = 0.7071067811865475*volIn[17] 
   surfOut[13] = 0.7071067811865475*volIn[18] 
   surfOut[14] = 0.7071067811865475*volIn[19] 
   surfOut[15] = 0.7071067811865475*volIn[21] 
   surfOut[16] = 0.7071067811865475*volIn[22] 
   surfOut[17] = 0.7071067811865475*volIn[23] 
   surfOut[18] = 0.7071067811865475*volIn[24] 
   surfOut[19] = 0.7071067811865475*volIn[26] 
   surfOut[20] = 0.7071067811865475*volIn[27] 
   surfOut[21] = 0.7071067811865475*volIn[28] 
end 

--    dir 4 
_M[2][4] = {} 
_M[2][4].upper = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[26]+1.224744871391589*volIn[5]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.224744871391589*volIn[11]+0.7071067811865475*volIn[2] 
   surfOut[3] = 1.224744871391589*volIn[12]+0.7071067811865475*volIn[3] 
   surfOut[4] = 1.224744871391589*volIn[13]+0.7071067811865475*volIn[4] 
   surfOut[5] = 1.224744871391589*volIn[17]+0.7071067811865475*volIn[6] 
   surfOut[6] = 1.224744871391589*volIn[21]+0.7071067811865475*volIn[7] 
   surfOut[7] = 0.7071067811865475*volIn[8] 
   surfOut[8] = 0.7071067811865475*volIn[9] 
   surfOut[9] = 0.7071067811865475*volIn[10] 
   surfOut[10] = 0.7071067811865475*volIn[14] 
   surfOut[11] = 0.7071067811865475*volIn[15] 
   surfOut[12] = 0.7071067811865475*volIn[16] 
   surfOut[13] = 0.7071067811865475*volIn[18] 
   surfOut[14] = 0.7071067811865475*volIn[19] 
   surfOut[15] = 0.7071067811865475*volIn[20] 
   surfOut[16] = 0.7071067811865475*volIn[22] 
   surfOut[17] = 0.7071067811865475*volIn[23] 
   surfOut[18] = 0.7071067811865475*volIn[24] 
   surfOut[19] = 0.7071067811865475*volIn[25] 
   surfOut[20] = 0.7071067811865475*volIn[27] 
   surfOut[21] = 0.7071067811865475*volIn[28] 
end 
_M[2][4].lower = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[26]-1.224744871391589*volIn[5]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[2]-1.224744871391589*volIn[11] 
   surfOut[3] = 0.7071067811865475*volIn[3]-1.224744871391589*volIn[12] 
   surfOut[4] = 0.7071067811865475*volIn[4]-1.224744871391589*volIn[13] 
   surfOut[5] = 0.7071067811865475*volIn[6]-1.224744871391589*volIn[17] 
   surfOut[6] = 0.7071067811865475*volIn[7]-1.224744871391589*volIn[21] 
   surfOut[7] = 0.7071067811865475*volIn[8] 
   surfOut[8] = 0.7071067811865475*volIn[9] 
   surfOut[9] = 0.7071067811865475*volIn[10] 
   surfOut[10] = 0.7071067811865475*volIn[14] 
   surfOut[11] = 0.7071067811865475*volIn[15] 
   surfOut[12] = 0.7071067811865475*volIn[16] 
   surfOut[13] = 0.7071067811865475*volIn[18] 
   surfOut[14] = 0.7071067811865475*volIn[19] 
   surfOut[15] = 0.7071067811865475*volIn[20] 
   surfOut[16] = 0.7071067811865475*volIn[22] 
   surfOut[17] = 0.7071067811865475*volIn[23] 
   surfOut[18] = 0.7071067811865475*volIn[24] 
   surfOut[19] = 0.7071067811865475*volIn[25] 
   surfOut[20] = 0.7071067811865475*volIn[27] 
   surfOut[21] = 0.7071067811865475*volIn[28] 
end 

--    dir 5 
_M[2][5] = {} 
_M[2][5].upper = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[27]+1.224744871391589*volIn[6]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.224744871391589*volIn[14]+0.7071067811865475*volIn[2] 
   surfOut[3] = 1.224744871391589*volIn[15]+0.7071067811865475*volIn[3] 
   surfOut[4] = 1.224744871391589*volIn[16]+0.7071067811865475*volIn[4] 
   surfOut[5] = 1.224744871391589*volIn[17]+0.7071067811865475*volIn[5] 
   surfOut[6] = 1.224744871391589*volIn[22]+0.7071067811865475*volIn[7] 
   surfOut[7] = 0.7071067811865475*volIn[8] 
   surfOut[8] = 0.7071067811865475*volIn[9] 
   surfOut[9] = 0.7071067811865475*volIn[10] 
   surfOut[10] = 0.7071067811865475*volIn[11] 
   surfOut[11] = 0.7071067811865475*volIn[12] 
   surfOut[12] = 0.7071067811865475*volIn[13] 
   surfOut[13] = 0.7071067811865475*volIn[18] 
   surfOut[14] = 0.7071067811865475*volIn[19] 
   surfOut[15] = 0.7071067811865475*volIn[20] 
   surfOut[16] = 0.7071067811865475*volIn[21] 
   surfOut[17] = 0.7071067811865475*volIn[23] 
   surfOut[18] = 0.7071067811865475*volIn[24] 
   surfOut[19] = 0.7071067811865475*volIn[25] 
   surfOut[20] = 0.7071067811865475*volIn[26] 
   surfOut[21] = 0.7071067811865475*volIn[28] 
end 
_M[2][5].lower = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[27]-1.224744871391589*volIn[6]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[2]-1.224744871391589*volIn[14] 
   surfOut[3] = 0.7071067811865475*volIn[3]-1.224744871391589*volIn[15] 
   surfOut[4] = 0.7071067811865475*volIn[4]-1.224744871391589*volIn[16] 
   surfOut[5] = 0.7071067811865475*volIn[5]-1.224744871391589*volIn[17] 
   surfOut[6] = 0.7071067811865475*volIn[7]-1.224744871391589*volIn[22] 
   surfOut[7] = 0.7071067811865475*volIn[8] 
   surfOut[8] = 0.7071067811865475*volIn[9] 
   surfOut[9] = 0.7071067811865475*volIn[10] 
   surfOut[10] = 0.7071067811865475*volIn[11] 
   surfOut[11] = 0.7071067811865475*volIn[12] 
   surfOut[12] = 0.7071067811865475*volIn[13] 
   surfOut[13] = 0.7071067811865475*volIn[18] 
   surfOut[14] = 0.7071067811865475*volIn[19] 
   surfOut[15] = 0.7071067811865475*volIn[20] 
   surfOut[16] = 0.7071067811865475*volIn[21] 
   surfOut[17] = 0.7071067811865475*volIn[23] 
   surfOut[18] = 0.7071067811865475*volIn[24] 
   surfOut[19] = 0.7071067811865475*volIn[25] 
   surfOut[20] = 0.7071067811865475*volIn[26] 
   surfOut[21] = 0.7071067811865475*volIn[28] 
end 

--    dir 6 
_M[2][6] = {} 
_M[2][6].upper = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[28]+1.224744871391589*volIn[7]+0.7071067811865475*volIn[1] 
   surfOut[2] = 1.224744871391589*volIn[18]+0.7071067811865475*volIn[2] 
   surfOut[3] = 1.224744871391589*volIn[19]+0.7071067811865475*volIn[3] 
   surfOut[4] = 1.224744871391589*volIn[20]+0.7071067811865475*volIn[4] 
   surfOut[5] = 1.224744871391589*volIn[21]+0.7071067811865475*volIn[5] 
   surfOut[6] = 1.224744871391589*volIn[22]+0.7071067811865475*volIn[6] 
   surfOut[7] = 0.7071067811865475*volIn[8] 
   surfOut[8] = 0.7071067811865475*volIn[9] 
   surfOut[9] = 0.7071067811865475*volIn[10] 
   surfOut[10] = 0.7071067811865475*volIn[11] 
   surfOut[11] = 0.7071067811865475*volIn[12] 
   surfOut[12] = 0.7071067811865475*volIn[13] 
   surfOut[13] = 0.7071067811865475*volIn[14] 
   surfOut[14] = 0.7071067811865475*volIn[15] 
   surfOut[15] = 0.7071067811865475*volIn[16] 
   surfOut[16] = 0.7071067811865475*volIn[17] 
   surfOut[17] = 0.7071067811865475*volIn[23] 
   surfOut[18] = 0.7071067811865475*volIn[24] 
   surfOut[19] = 0.7071067811865475*volIn[25] 
   surfOut[20] = 0.7071067811865475*volIn[26] 
   surfOut[21] = 0.7071067811865475*volIn[27] 
end 
_M[2][6].lower = function (volIn, surfOut) 
   surfOut[1] = 1.58113883008419*volIn[28]-1.224744871391589*volIn[7]+0.7071067811865475*volIn[1] 
   surfOut[2] = 0.7071067811865475*volIn[2]-1.224744871391589*volIn[18] 
   surfOut[3] = 0.7071067811865475*volIn[3]-1.224744871391589*volIn[19] 
   surfOut[4] = 0.7071067811865475*volIn[4]-1.224744871391589*volIn[20] 
   surfOut[5] = 0.7071067811865475*volIn[5]-1.224744871391589*volIn[21] 
   surfOut[6] = 0.7071067811865475*volIn[6]-1.224744871391589*volIn[22] 
   surfOut[7] = 0.7071067811865475*volIn[8] 
   surfOut[8] = 0.7071067811865475*volIn[9] 
   surfOut[9] = 0.7071067811865475*volIn[10] 
   surfOut[10] = 0.7071067811865475*volIn[11] 
   surfOut[11] = 0.7071067811865475*volIn[12] 
   surfOut[12] = 0.7071067811865475*volIn[13] 
   surfOut[13] = 0.7071067811865475*volIn[14] 
   surfOut[14] = 0.7071067811865475*volIn[15] 
   surfOut[15] = 0.7071067811865475*volIn[16] 
   surfOut[16] = 0.7071067811865475*volIn[17] 
   surfOut[17] = 0.7071067811865475*volIn[23] 
   surfOut[18] = 0.7071067811865475*volIn[24] 
   surfOut[19] = 0.7071067811865475*volIn[25] 
   surfOut[20] = 0.7071067811865475*volIn[26] 
   surfOut[21] = 0.7071067811865475*volIn[27] 
end 


return _M 
