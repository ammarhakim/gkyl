local _M = {} 
_M[1] = function (dir, f, out) 
   if dir == 1  then 
   out[1] = f[1] 
   out[2] = -1.0*f[2] 
   out[3] = f[3] 
   out[4] = f[4] 
   end 
   if dir == 2  then 
   out[1] = f[1] 
   out[2] = f[2] 
   out[3] = -1.0*f[3] 
   out[4] = f[4] 
   end 
   if dir == 3  then 
   out[1] = f[1] 
   out[2] = f[2] 
   out[3] = f[3] 
   out[4] = -1.0*f[4] 
   end 
end 
_M[2] = function (dir, f, out) 
   if dir == 1  then 
   out[1] = f[1] 
   out[2] = -1.0*f[2] 
   out[3] = f[3] 
   out[4] = f[4] 
   out[5] = -1.0*f[5] 
   out[6] = -1.0*f[6] 
   out[7] = f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   end 
   if dir == 2  then 
   out[1] = f[1] 
   out[2] = f[2] 
   out[3] = -1.0*f[3] 
   out[4] = f[4] 
   out[5] = -1.0*f[5] 
   out[6] = f[6] 
   out[7] = -1.0*f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   end 
   if dir == 3  then 
   out[1] = f[1] 
   out[2] = f[2] 
   out[3] = f[3] 
   out[4] = -1.0*f[4] 
   out[5] = f[5] 
   out[6] = -1.0*f[6] 
   out[7] = -1.0*f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   end 
end 
_M[3] = function (dir, f, out) 
   if dir == 1  then 
   out[1] = f[1] 
   out[2] = -1.0*f[2] 
   out[3] = f[3] 
   out[4] = f[4] 
   out[5] = -1.0*f[5] 
   out[6] = -1.0*f[6] 
   out[7] = f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = -1.0*f[11] 
   out[12] = f[12] 
   out[13] = -1.0*f[13] 
   out[14] = f[14] 
   out[15] = f[15] 
   out[16] = -1.0*f[16] 
   out[17] = f[17] 
   out[18] = -1.0*f[18] 
   out[19] = f[19] 
   out[20] = f[20] 
   end 
   if dir == 2  then 
   out[1] = f[1] 
   out[2] = f[2] 
   out[3] = -1.0*f[3] 
   out[4] = f[4] 
   out[5] = -1.0*f[5] 
   out[6] = f[6] 
   out[7] = -1.0*f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = -1.0*f[11] 
   out[12] = -1.0*f[12] 
   out[13] = f[13] 
   out[14] = f[14] 
   out[15] = f[15] 
   out[16] = f[16] 
   out[17] = -1.0*f[17] 
   out[18] = f[18] 
   out[19] = -1.0*f[19] 
   out[20] = f[20] 
   end 
   if dir == 3  then 
   out[1] = f[1] 
   out[2] = f[2] 
   out[3] = f[3] 
   out[4] = -1.0*f[4] 
   out[5] = f[5] 
   out[6] = -1.0*f[6] 
   out[7] = -1.0*f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = -1.0*f[11] 
   out[12] = f[12] 
   out[13] = f[13] 
   out[14] = -1.0*f[14] 
   out[15] = -1.0*f[15] 
   out[16] = f[16] 
   out[17] = f[17] 
   out[18] = f[18] 
   out[19] = f[19] 
   out[20] = -1.0*f[20] 
   end 
end 
_M[4] = function (dir, f, out) 
   if dir == 1  then 
   out[1] = f[1] 
   out[2] = -1.0*f[2] 
   out[3] = f[3] 
   out[4] = f[4] 
   out[5] = -1.0*f[5] 
   out[6] = -1.0*f[6] 
   out[7] = f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = -1.0*f[11] 
   out[12] = f[12] 
   out[13] = -1.0*f[13] 
   out[14] = f[14] 
   out[15] = f[15] 
   out[16] = -1.0*f[16] 
   out[17] = f[17] 
   out[18] = -1.0*f[18] 
   out[19] = f[19] 
   out[20] = f[20] 
   out[21] = f[21] 
   out[22] = -1.0*f[22] 
   out[23] = -1.0*f[23] 
   out[24] = f[24] 
   out[25] = f[25] 
   out[26] = f[26] 
   out[27] = -1.0*f[27] 
   out[28] = -1.0*f[28] 
   out[29] = -1.0*f[29] 
   out[30] = f[30] 
   out[31] = -1.0*f[31] 
   out[32] = f[32] 
   out[33] = f[33] 
   out[34] = f[34] 
   out[35] = f[35] 
   end 
   if dir == 2  then 
   out[1] = f[1] 
   out[2] = f[2] 
   out[3] = -1.0*f[3] 
   out[4] = f[4] 
   out[5] = -1.0*f[5] 
   out[6] = f[6] 
   out[7] = -1.0*f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = -1.0*f[11] 
   out[12] = -1.0*f[12] 
   out[13] = f[13] 
   out[14] = f[14] 
   out[15] = f[15] 
   out[16] = f[16] 
   out[17] = -1.0*f[17] 
   out[18] = f[18] 
   out[19] = -1.0*f[19] 
   out[20] = f[20] 
   out[21] = -1.0*f[21] 
   out[22] = f[22] 
   out[23] = -1.0*f[23] 
   out[24] = f[24] 
   out[25] = f[25] 
   out[26] = f[26] 
   out[27] = -1.0*f[27] 
   out[28] = -1.0*f[28] 
   out[29] = f[29] 
   out[30] = -1.0*f[30] 
   out[31] = f[31] 
   out[32] = -1.0*f[32] 
   out[33] = f[33] 
   out[34] = f[34] 
   out[35] = f[35] 
   end 
   if dir == 3  then 
   out[1] = f[1] 
   out[2] = f[2] 
   out[3] = f[3] 
   out[4] = -1.0*f[4] 
   out[5] = f[5] 
   out[6] = -1.0*f[6] 
   out[7] = -1.0*f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = -1.0*f[11] 
   out[12] = f[12] 
   out[13] = f[13] 
   out[14] = -1.0*f[14] 
   out[15] = -1.0*f[15] 
   out[16] = f[16] 
   out[17] = f[17] 
   out[18] = f[18] 
   out[19] = f[19] 
   out[20] = -1.0*f[20] 
   out[21] = -1.0*f[21] 
   out[22] = -1.0*f[22] 
   out[23] = f[23] 
   out[24] = f[24] 
   out[25] = f[25] 
   out[26] = f[26] 
   out[27] = f[27] 
   out[28] = f[28] 
   out[29] = -1.0*f[29] 
   out[30] = -1.0*f[30] 
   out[31] = -1.0*f[31] 
   out[32] = -1.0*f[32] 
   out[33] = f[33] 
   out[34] = f[34] 
   out[35] = f[35] 
   end 
end 
_M[5] = function (dir, f, out) 
   if dir == 1  then 
   out[1] = f[1] 
   out[2] = -1.0*f[2] 
   out[3] = f[3] 
   out[4] = f[4] 
   out[5] = -1.0*f[5] 
   out[6] = -1.0*f[6] 
   out[7] = f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = -1.0*f[11] 
   out[12] = f[12] 
   out[13] = -1.0*f[13] 
   out[14] = f[14] 
   out[15] = f[15] 
   out[16] = -1.0*f[16] 
   out[17] = f[17] 
   out[18] = -1.0*f[18] 
   out[19] = f[19] 
   out[20] = f[20] 
   out[21] = f[21] 
   out[22] = -1.0*f[22] 
   out[23] = -1.0*f[23] 
   out[24] = f[24] 
   out[25] = f[25] 
   out[26] = f[26] 
   out[27] = -1.0*f[27] 
   out[28] = -1.0*f[28] 
   out[29] = -1.0*f[29] 
   out[30] = f[30] 
   out[31] = -1.0*f[31] 
   out[32] = f[32] 
   out[33] = f[33] 
   out[34] = f[34] 
   out[35] = f[35] 
   out[36] = f[36] 
   out[37] = f[37] 
   out[38] = -1.0*f[38] 
   out[39] = -1.0*f[39] 
   out[40] = -1.0*f[40] 
   out[41] = -1.0*f[41] 
   out[42] = -1.0*f[42] 
   out[43] = f[43] 
   out[44] = -1.0*f[44] 
   out[45] = f[45] 
   out[46] = f[46] 
   out[47] = f[47] 
   out[48] = f[48] 
   out[49] = -1.0*f[49] 
   out[50] = f[50] 
   out[51] = f[51] 
   out[52] = -1.0*f[52] 
   out[53] = f[53] 
   out[54] = -1.0*f[54] 
   out[55] = f[55] 
   out[56] = f[56] 
   end 
   if dir == 2  then 
   out[1] = f[1] 
   out[2] = f[2] 
   out[3] = -1.0*f[3] 
   out[4] = f[4] 
   out[5] = -1.0*f[5] 
   out[6] = f[6] 
   out[7] = -1.0*f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = -1.0*f[11] 
   out[12] = -1.0*f[12] 
   out[13] = f[13] 
   out[14] = f[14] 
   out[15] = f[15] 
   out[16] = f[16] 
   out[17] = -1.0*f[17] 
   out[18] = f[18] 
   out[19] = -1.0*f[19] 
   out[20] = f[20] 
   out[21] = -1.0*f[21] 
   out[22] = f[22] 
   out[23] = -1.0*f[23] 
   out[24] = f[24] 
   out[25] = f[25] 
   out[26] = f[26] 
   out[27] = -1.0*f[27] 
   out[28] = -1.0*f[28] 
   out[29] = f[29] 
   out[30] = -1.0*f[30] 
   out[31] = f[31] 
   out[32] = -1.0*f[32] 
   out[33] = f[33] 
   out[34] = f[34] 
   out[35] = f[35] 
   out[36] = f[36] 
   out[37] = -1.0*f[37] 
   out[38] = f[38] 
   out[39] = -1.0*f[39] 
   out[40] = -1.0*f[40] 
   out[41] = -1.0*f[41] 
   out[42] = f[42] 
   out[43] = -1.0*f[43] 
   out[44] = f[44] 
   out[45] = -1.0*f[45] 
   out[46] = f[46] 
   out[47] = f[47] 
   out[48] = -1.0*f[48] 
   out[49] = f[49] 
   out[50] = f[50] 
   out[51] = f[51] 
   out[52] = f[52] 
   out[53] = -1.0*f[53] 
   out[54] = f[54] 
   out[55] = -1.0*f[55] 
   out[56] = f[56] 
   end 
   if dir == 3  then 
   out[1] = f[1] 
   out[2] = f[2] 
   out[3] = f[3] 
   out[4] = -1.0*f[4] 
   out[5] = f[5] 
   out[6] = -1.0*f[6] 
   out[7] = -1.0*f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = -1.0*f[11] 
   out[12] = f[12] 
   out[13] = f[13] 
   out[14] = -1.0*f[14] 
   out[15] = -1.0*f[15] 
   out[16] = f[16] 
   out[17] = f[17] 
   out[18] = f[18] 
   out[19] = f[19] 
   out[20] = -1.0*f[20] 
   out[21] = -1.0*f[21] 
   out[22] = -1.0*f[22] 
   out[23] = f[23] 
   out[24] = f[24] 
   out[25] = f[25] 
   out[26] = f[26] 
   out[27] = f[27] 
   out[28] = f[28] 
   out[29] = -1.0*f[29] 
   out[30] = -1.0*f[30] 
   out[31] = -1.0*f[31] 
   out[32] = -1.0*f[32] 
   out[33] = f[33] 
   out[34] = f[34] 
   out[35] = f[35] 
   out[36] = -1.0*f[36] 
   out[37] = f[37] 
   out[38] = f[38] 
   out[39] = -1.0*f[39] 
   out[40] = -1.0*f[40] 
   out[41] = -1.0*f[41] 
   out[42] = f[42] 
   out[43] = f[43] 
   out[44] = f[44] 
   out[45] = f[45] 
   out[46] = -1.0*f[46] 
   out[47] = -1.0*f[47] 
   out[48] = f[48] 
   out[49] = f[49] 
   out[50] = -1.0*f[50] 
   out[51] = -1.0*f[51] 
   out[52] = f[52] 
   out[53] = f[53] 
   out[54] = f[54] 
   out[55] = f[55] 
   out[56] = -1.0*f[56] 
   end 
end 
_M[6] = function (dir, f, out) 
   if dir == 1  then 
   out[1] = f[1] 
   out[2] = -1.0*f[2] 
   out[3] = f[3] 
   out[4] = f[4] 
   out[5] = -1.0*f[5] 
   out[6] = -1.0*f[6] 
   out[7] = f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = -1.0*f[11] 
   out[12] = f[12] 
   out[13] = -1.0*f[13] 
   out[14] = f[14] 
   out[15] = f[15] 
   out[16] = -1.0*f[16] 
   out[17] = f[17] 
   out[18] = -1.0*f[18] 
   out[19] = f[19] 
   out[20] = f[20] 
   out[21] = f[21] 
   out[22] = -1.0*f[22] 
   out[23] = -1.0*f[23] 
   out[24] = f[24] 
   out[25] = f[25] 
   out[26] = f[26] 
   out[27] = -1.0*f[27] 
   out[28] = -1.0*f[28] 
   out[29] = -1.0*f[29] 
   out[30] = f[30] 
   out[31] = -1.0*f[31] 
   out[32] = f[32] 
   out[33] = f[33] 
   out[34] = f[34] 
   out[35] = f[35] 
   out[36] = f[36] 
   out[37] = f[37] 
   out[38] = -1.0*f[38] 
   out[39] = -1.0*f[39] 
   out[40] = -1.0*f[40] 
   out[41] = -1.0*f[41] 
   out[42] = -1.0*f[42] 
   out[43] = f[43] 
   out[44] = -1.0*f[44] 
   out[45] = f[45] 
   out[46] = f[46] 
   out[47] = f[47] 
   out[48] = f[48] 
   out[49] = -1.0*f[49] 
   out[50] = f[50] 
   out[51] = f[51] 
   out[52] = -1.0*f[52] 
   out[53] = f[53] 
   out[54] = -1.0*f[54] 
   out[55] = f[55] 
   out[56] = f[56] 
   out[57] = f[57] 
   out[58] = -1.0*f[58] 
   out[59] = f[59] 
   out[60] = -1.0*f[60] 
   out[61] = -1.0*f[61] 
   out[62] = f[62] 
   out[63] = -1.0*f[63] 
   out[64] = -1.0*f[64] 
   out[65] = -1.0*f[65] 
   out[66] = f[66] 
   out[67] = f[67] 
   out[68] = -1.0*f[68] 
   out[69] = -1.0*f[69] 
   out[70] = f[70] 
   out[71] = f[71] 
   out[72] = f[72] 
   out[73] = f[73] 
   out[74] = f[74] 
   out[75] = f[75] 
   out[76] = -1.0*f[76] 
   out[77] = -1.0*f[77] 
   out[78] = -1.0*f[78] 
   out[79] = f[79] 
   out[80] = -1.0*f[80] 
   out[81] = f[81] 
   out[82] = f[82] 
   out[83] = f[83] 
   out[84] = f[84] 
   end 
   if dir == 2  then 
   out[1] = f[1] 
   out[2] = f[2] 
   out[3] = -1.0*f[3] 
   out[4] = f[4] 
   out[5] = -1.0*f[5] 
   out[6] = f[6] 
   out[7] = -1.0*f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = -1.0*f[11] 
   out[12] = -1.0*f[12] 
   out[13] = f[13] 
   out[14] = f[14] 
   out[15] = f[15] 
   out[16] = f[16] 
   out[17] = -1.0*f[17] 
   out[18] = f[18] 
   out[19] = -1.0*f[19] 
   out[20] = f[20] 
   out[21] = -1.0*f[21] 
   out[22] = f[22] 
   out[23] = -1.0*f[23] 
   out[24] = f[24] 
   out[25] = f[25] 
   out[26] = f[26] 
   out[27] = -1.0*f[27] 
   out[28] = -1.0*f[28] 
   out[29] = f[29] 
   out[30] = -1.0*f[30] 
   out[31] = f[31] 
   out[32] = -1.0*f[32] 
   out[33] = f[33] 
   out[34] = f[34] 
   out[35] = f[35] 
   out[36] = f[36] 
   out[37] = -1.0*f[37] 
   out[38] = f[38] 
   out[39] = -1.0*f[39] 
   out[40] = -1.0*f[40] 
   out[41] = -1.0*f[41] 
   out[42] = f[42] 
   out[43] = -1.0*f[43] 
   out[44] = f[44] 
   out[45] = -1.0*f[45] 
   out[46] = f[46] 
   out[47] = f[47] 
   out[48] = -1.0*f[48] 
   out[49] = f[49] 
   out[50] = f[50] 
   out[51] = f[51] 
   out[52] = f[52] 
   out[53] = -1.0*f[53] 
   out[54] = f[54] 
   out[55] = -1.0*f[55] 
   out[56] = f[56] 
   out[57] = f[57] 
   out[58] = f[58] 
   out[59] = -1.0*f[59] 
   out[60] = -1.0*f[60] 
   out[61] = -1.0*f[61] 
   out[62] = -1.0*f[62] 
   out[63] = f[63] 
   out[64] = -1.0*f[64] 
   out[65] = f[65] 
   out[66] = -1.0*f[66] 
   out[67] = -1.0*f[67] 
   out[68] = f[68] 
   out[69] = -1.0*f[69] 
   out[70] = f[70] 
   out[71] = f[71] 
   out[72] = f[72] 
   out[73] = f[73] 
   out[74] = f[74] 
   out[75] = f[75] 
   out[76] = -1.0*f[76] 
   out[77] = -1.0*f[77] 
   out[78] = f[78] 
   out[79] = -1.0*f[79] 
   out[80] = f[80] 
   out[81] = -1.0*f[81] 
   out[82] = f[82] 
   out[83] = f[83] 
   out[84] = f[84] 
   end 
   if dir == 3  then 
   out[1] = f[1] 
   out[2] = f[2] 
   out[3] = f[3] 
   out[4] = -1.0*f[4] 
   out[5] = f[5] 
   out[6] = -1.0*f[6] 
   out[7] = -1.0*f[7] 
   out[8] = f[8] 
   out[9] = f[9] 
   out[10] = f[10] 
   out[11] = -1.0*f[11] 
   out[12] = f[12] 
   out[13] = f[13] 
   out[14] = -1.0*f[14] 
   out[15] = -1.0*f[15] 
   out[16] = f[16] 
   out[17] = f[17] 
   out[18] = f[18] 
   out[19] = f[19] 
   out[20] = -1.0*f[20] 
   out[21] = -1.0*f[21] 
   out[22] = -1.0*f[22] 
   out[23] = f[23] 
   out[24] = f[24] 
   out[25] = f[25] 
   out[26] = f[26] 
   out[27] = f[27] 
   out[28] = f[28] 
   out[29] = -1.0*f[29] 
   out[30] = -1.0*f[30] 
   out[31] = -1.0*f[31] 
   out[32] = -1.0*f[32] 
   out[33] = f[33] 
   out[34] = f[34] 
   out[35] = f[35] 
   out[36] = -1.0*f[36] 
   out[37] = f[37] 
   out[38] = f[38] 
   out[39] = -1.0*f[39] 
   out[40] = -1.0*f[40] 
   out[41] = -1.0*f[41] 
   out[42] = f[42] 
   out[43] = f[43] 
   out[44] = f[44] 
   out[45] = f[45] 
   out[46] = -1.0*f[46] 
   out[47] = -1.0*f[47] 
   out[48] = f[48] 
   out[49] = f[49] 
   out[50] = -1.0*f[50] 
   out[51] = -1.0*f[51] 
   out[52] = f[52] 
   out[53] = f[53] 
   out[54] = f[54] 
   out[55] = f[55] 
   out[56] = -1.0*f[56] 
   out[57] = f[57] 
   out[58] = -1.0*f[58] 
   out[59] = -1.0*f[59] 
   out[60] = f[60] 
   out[61] = f[61] 
   out[62] = -1.0*f[62] 
   out[63] = -1.0*f[63] 
   out[64] = f[64] 
   out[65] = -1.0*f[65] 
   out[66] = -1.0*f[66] 
   out[67] = -1.0*f[67] 
   out[68] = -1.0*f[68] 
   out[69] = f[69] 
   out[70] = f[70] 
   out[71] = f[71] 
   out[72] = f[72] 
   out[73] = f[73] 
   out[74] = f[74] 
   out[75] = f[75] 
   out[76] = f[76] 
   out[77] = f[77] 
   out[78] = -1.0*f[78] 
   out[79] = -1.0*f[79] 
   out[80] = -1.0*f[80] 
   out[81] = -1.0*f[81] 
   out[82] = f[82] 
   out[83] = f[83] 
   out[84] = f[84] 
   end 
end 
return _M 
