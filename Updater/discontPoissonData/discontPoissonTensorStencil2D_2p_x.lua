local Lin = require "Lib.Linalg"
local function stencilFn(dx)
  local _M = {}

  _M[1] = Lin.Mat(9,9)
  _M[1][1][1] = 3.75/dx[1]^2
  _M[1][1][2] = 4.763139720814412/dx[1]^2
  _M[1][1][3] = 0.0
  _M[1][1][4] = 0.0
  _M[1][1][5] = 2.683281572999748/dx[1]^2
  _M[1][1][6] = 0.0
  _M[1][1][7] = 0.0
  _M[1][1][8] = 0.0
  _M[1][1][9] = 0.0
  _M[1][2][1] = -4.763139720814412/dx[1]^2
  _M[1][2][2] = -5.8125/dx[1]^2
  _M[1][2][3] = 0.0
  _M[1][2][4] = 0.0
  _M[1][2][5] = -2.953149801483155/dx[1]^2
  _M[1][2][6] = 0.0
  _M[1][2][7] = 0.0
  _M[1][2][8] = 0.0
  _M[1][2][9] = 0.0
  _M[1][3][1] = 0.0
  _M[1][3][2] = 0.0
  _M[1][3][3] = 3.75/dx[1]^2
  _M[1][3][4] = 4.763139720814412/dx[1]^2
  _M[1][3][5] = 0.0
  _M[1][3][6] = 0.0
  _M[1][3][7] = 2.683281572999748/dx[1]^2
  _M[1][3][8] = 0.0
  _M[1][3][9] = 0.0
  _M[1][4][1] = 0.0
  _M[1][4][2] = 0.0
  _M[1][4][3] = -4.763139720814412/dx[1]^2
  _M[1][4][4] = -5.8125/dx[1]^2
  _M[1][4][5] = 0.0
  _M[1][4][6] = 0.0
  _M[1][4][7] = -2.953149801483155/dx[1]^2
  _M[1][4][8] = 0.0
  _M[1][4][9] = 0.0
  _M[1][5][1] = 1.677050983124842/dx[1]^2
  _M[1][5][2] = 1.210307295689818/dx[1]^2
  _M[1][5][3] = 0.0
  _M[1][5][4] = 0.0
  _M[1][5][5] = -0.5625/dx[1]^2
  _M[1][5][6] = 0.0
  _M[1][5][7] = 0.0
  _M[1][5][8] = 0.0
  _M[1][5][9] = 0.0
  _M[1][6][1] = 0.0
  _M[1][6][2] = 0.0
  _M[1][6][3] = 0.0
  _M[1][6][4] = 0.0
  _M[1][6][5] = 0.0
  _M[1][6][6] = 3.75/dx[1]^2
  _M[1][6][7] = 0.0
  _M[1][6][8] = 4.763139720814412/dx[1]^2
  _M[1][6][9] = 2.683281572999748/dx[1]^2
  _M[1][7][1] = 0.0
  _M[1][7][2] = 0.0
  _M[1][7][3] = 1.677050983124842/dx[1]^2
  _M[1][7][4] = 1.210307295689818/dx[1]^2
  _M[1][7][5] = 0.0
  _M[1][7][6] = 0.0
  _M[1][7][7] = -0.5625/dx[1]^2
  _M[1][7][8] = 0.0
  _M[1][7][9] = 0.0
  _M[1][8][1] = 0.0
  _M[1][8][2] = 0.0
  _M[1][8][3] = 0.0
  _M[1][8][4] = 0.0
  _M[1][8][5] = 0.0
  _M[1][8][6] = -4.763139720814412/dx[1]^2
  _M[1][8][7] = 0.0
  _M[1][8][8] = -5.8125/dx[1]^2
  _M[1][8][9] = -2.953149801483155/dx[1]^2
  _M[1][9][1] = 0.0
  _M[1][9][2] = 0.0
  _M[1][9][3] = 0.0
  _M[1][9][4] = 0.0
  _M[1][9][5] = 0.0
  _M[1][9][6] = 1.677050983124842/dx[1]^2
  _M[1][9][7] = 0.0
  _M[1][9][8] = 1.210307295689818/dx[1]^2
  _M[1][9][9] = -0.5625/dx[1]^2
  _M[2] = Lin.Mat(9,9)
  _M[2][1][1] = -7.5/dx[1]^2
  _M[2][1][2] = 0.0
  _M[2][1][3] = 0.0
  _M[2][1][4] = 0.0
  _M[2][1][5] = -5.366563145999495/dx[1]^2
  _M[2][1][6] = 0.0
  _M[2][1][7] = 0.0
  _M[2][1][8] = 0.0
  _M[2][1][9] = 0.0
  _M[2][2][1] = 0.0
  _M[2][2][2] = -21.375/dx[1]^2
  _M[2][2][3] = 0.0
  _M[2][2][4] = 0.0
  _M[2][2][5] = 0.0
  _M[2][2][6] = 0.0
  _M[2][2][7] = 0.0
  _M[2][2][8] = 0.0
  _M[2][2][9] = 0.0
  _M[2][3][1] = 0.0
  _M[2][3][2] = 0.0
  _M[2][3][3] = -7.5/dx[1]^2
  _M[2][3][4] = 0.0
  _M[2][3][5] = 0.0
  _M[2][3][6] = 0.0
  _M[2][3][7] = -5.366563145999495/dx[1]^2
  _M[2][3][8] = 0.0
  _M[2][3][9] = 0.0
  _M[2][4][1] = 0.0
  _M[2][4][2] = 0.0
  _M[2][4][3] = 0.0
  _M[2][4][4] = -21.375/dx[1]^2
  _M[2][4][5] = 0.0
  _M[2][4][6] = 0.0
  _M[2][4][7] = 0.0
  _M[2][4][8] = 0.0
  _M[2][4][9] = 0.0
  _M[2][5][1] = -3.354101966249685/dx[1]^2
  _M[2][5][2] = 0.0
  _M[2][5][3] = 0.0
  _M[2][5][4] = 0.0
  _M[2][5][5] = -25.125/dx[1]^2
  _M[2][5][6] = 0.0
  _M[2][5][7] = 0.0
  _M[2][5][8] = 0.0
  _M[2][5][9] = 0.0
  _M[2][6][1] = 0.0
  _M[2][6][2] = 0.0
  _M[2][6][3] = 0.0
  _M[2][6][4] = 0.0
  _M[2][6][5] = 0.0
  _M[2][6][6] = -7.5/dx[1]^2
  _M[2][6][7] = 0.0
  _M[2][6][8] = 0.0
  _M[2][6][9] = -5.366563145999495/dx[1]^2
  _M[2][7][1] = 0.0
  _M[2][7][2] = 0.0
  _M[2][7][3] = -3.354101966249684/dx[1]^2
  _M[2][7][4] = 0.0
  _M[2][7][5] = 0.0
  _M[2][7][6] = 0.0
  _M[2][7][7] = -25.125/dx[1]^2
  _M[2][7][8] = 0.0
  _M[2][7][9] = 0.0
  _M[2][8][1] = 0.0
  _M[2][8][2] = 0.0
  _M[2][8][3] = 0.0
  _M[2][8][4] = 0.0
  _M[2][8][5] = 0.0
  _M[2][8][6] = 0.0
  _M[2][8][7] = 0.0
  _M[2][8][8] = -21.375/dx[1]^2
  _M[2][8][9] = 0.0
  _M[2][9][1] = 0.0
  _M[2][9][2] = 0.0
  _M[2][9][3] = 0.0
  _M[2][9][4] = 0.0
  _M[2][9][5] = 0.0
  _M[2][9][6] = -3.354101966249685/dx[1]^2
  _M[2][9][7] = 0.0
  _M[2][9][8] = 0.0
  _M[2][9][9] = -25.125/dx[1]^2
  _M[3] = Lin.Mat(9,9)
  _M[3][1][1] = 3.75/dx[1]^2
  _M[3][1][2] = -4.763139720814412/dx[1]^2
  _M[3][1][3] = 0.0
  _M[3][1][4] = 0.0
  _M[3][1][5] = 2.683281572999748/dx[1]^2
  _M[3][1][6] = 0.0
  _M[3][1][7] = 0.0
  _M[3][1][8] = 0.0
  _M[3][1][9] = 0.0
  _M[3][2][1] = 4.763139720814412/dx[1]^2
  _M[3][2][2] = -5.8125/dx[1]^2
  _M[3][2][3] = 0.0
  _M[3][2][4] = 0.0
  _M[3][2][5] = 2.953149801483155/dx[1]^2
  _M[3][2][6] = 0.0
  _M[3][2][7] = 0.0
  _M[3][2][8] = 0.0
  _M[3][2][9] = 0.0
  _M[3][3][1] = 0.0
  _M[3][3][2] = 0.0
  _M[3][3][3] = 3.75/dx[1]^2
  _M[3][3][4] = -4.763139720814412/dx[1]^2
  _M[3][3][5] = 0.0
  _M[3][3][6] = 0.0
  _M[3][3][7] = 2.683281572999748/dx[1]^2
  _M[3][3][8] = 0.0
  _M[3][3][9] = 0.0
  _M[3][4][1] = 0.0
  _M[3][4][2] = 0.0
  _M[3][4][3] = 4.763139720814412/dx[1]^2
  _M[3][4][4] = -5.8125/dx[1]^2
  _M[3][4][5] = 0.0
  _M[3][4][6] = 0.0
  _M[3][4][7] = 2.953149801483155/dx[1]^2
  _M[3][4][8] = 0.0
  _M[3][4][9] = 0.0
  _M[3][5][1] = 1.677050983124842/dx[1]^2
  _M[3][5][2] = -1.210307295689818/dx[1]^2
  _M[3][5][3] = 0.0
  _M[3][5][4] = 0.0
  _M[3][5][5] = -0.5625/dx[1]^2
  _M[3][5][6] = 0.0
  _M[3][5][7] = 0.0
  _M[3][5][8] = 0.0
  _M[3][5][9] = 0.0
  _M[3][6][1] = 0.0
  _M[3][6][2] = 0.0
  _M[3][6][3] = 0.0
  _M[3][6][4] = 0.0
  _M[3][6][5] = 0.0
  _M[3][6][6] = 3.75/dx[1]^2
  _M[3][6][7] = 0.0
  _M[3][6][8] = -4.763139720814412/dx[1]^2
  _M[3][6][9] = 2.683281572999748/dx[1]^2
  _M[3][7][1] = 0.0
  _M[3][7][2] = 0.0
  _M[3][7][3] = 1.677050983124842/dx[1]^2
  _M[3][7][4] = -1.210307295689818/dx[1]^2
  _M[3][7][5] = 0.0
  _M[3][7][6] = 0.0
  _M[3][7][7] = -0.5625/dx[1]^2
  _M[3][7][8] = 0.0
  _M[3][7][9] = 0.0
  _M[3][8][1] = 0.0
  _M[3][8][2] = 0.0
  _M[3][8][3] = 0.0
  _M[3][8][4] = 0.0
  _M[3][8][5] = 0.0
  _M[3][8][6] = 4.763139720814412/dx[1]^2
  _M[3][8][7] = 0.0
  _M[3][8][8] = -5.8125/dx[1]^2
  _M[3][8][9] = 2.953149801483155/dx[1]^2
  _M[3][9][1] = 0.0
  _M[3][9][2] = 0.0
  _M[3][9][3] = 0.0
  _M[3][9][4] = 0.0
  _M[3][9][5] = 0.0
  _M[3][9][6] = 1.677050983124842/dx[1]^2
  _M[3][9][7] = 0.0
  _M[3][9][8] = -1.210307295689818/dx[1]^2
  _M[3][9][9] = -0.5625/dx[1]^2
  return(_M)
end

return(stencilFn)