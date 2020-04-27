local Lin = require "Lib.Linalg"
local function stencilFn(dx)
  local _M = {}

  _M[1] = Lin.Mat(4,4)
  _M[1][1][1] = 2.25/dx[2]^2
  _M[1][1][2] = 0.0
  _M[1][1][3] = 2.165063509461096/dx[2]^2
  _M[1][1][4] = 0.0
  _M[1][2][1] = 0.0
  _M[1][2][2] = 2.25/dx[2]^2
  _M[1][2][3] = 0.0
  _M[1][2][4] = 2.165063509461096/dx[2]^2
  _M[1][3][1] = -2.165063509461096/dx[2]^2
  _M[1][3][2] = 0.0
  _M[1][3][3] = -1.75/dx[2]^2
  _M[1][3][4] = 0.0
  _M[1][4][1] = 0.0
  _M[1][4][2] = -2.165063509461096/dx[2]^2
  _M[1][4][3] = 0.0
  _M[1][4][4] = -1.75/dx[2]^2
  _M[2] = Lin.Mat(4,4)
  _M[2][1][1] = -4.5/dx[2]^2
  _M[2][1][2] = 0.0
  _M[2][1][3] = 0.0
  _M[2][1][4] = 0.0
  _M[2][2][1] = 0.0
  _M[2][2][2] = -4.5/dx[2]^2
  _M[2][2][3] = 0.0
  _M[2][2][4] = 0.0
  _M[2][3][1] = 0.0
  _M[2][3][2] = 0.0
  _M[2][3][3] = -11.5/dx[2]^2
  _M[2][3][4] = 0.0
  _M[2][4][1] = 0.0
  _M[2][4][2] = 0.0
  _M[2][4][3] = 0.0
  _M[2][4][4] = -11.5/dx[2]^2
  _M[3] = Lin.Mat(4,4)
  _M[3][1][1] = 2.25/dx[2]^2
  _M[3][1][2] = 0.0
  _M[3][1][3] = -2.165063509461096/dx[2]^2
  _M[3][1][4] = 0.0
  _M[3][2][1] = 0.0
  _M[3][2][2] = 2.25/dx[2]^2
  _M[3][2][3] = 0.0
  _M[3][2][4] = -2.165063509461096/dx[2]^2
  _M[3][3][1] = 2.165063509461096/dx[2]^2
  _M[3][3][2] = 0.0
  _M[3][3][3] = -1.75/dx[2]^2
  _M[3][3][4] = 0.0
  _M[3][4][1] = 0.0
  _M[3][4][2] = 2.165063509461096/dx[2]^2
  _M[3][4][3] = 0.0
  _M[3][4][4] = -1.75/dx[2]^2
  return(_M)
end

return(stencilFn)