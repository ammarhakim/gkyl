local Lin = require "Lib.Linalg"
local function stencilFn(dx,
                         DxxC, DyyC, DxyC,
                         DxxL, DyyL, DxyL,
                         DxxR, DyyR, DxyR,
                         DxxB, DyyB, DxyB,
                         DxxT, DyyT, DxyT,
                         bcDx, bcNx, bcValx,
                         bcDy, bcNy, bcValy)
  local _M = {}

  local surfDyyT = {}
  surfDyyT[1] = (-0.408248290463863*DyyT[3])+0.408248290463863*DyyC[3]+0.3535533905932737*DyyT[1]+0.3535533905932737*DyyC[1]
  surfDyyT[2] = (-0.408248290463863*DyyT[4])+0.408248290463863*DyyC[4]+0.3535533905932737*DyyT[2]+0.3535533905932737*DyyC[2]

  local surfDyyB = {}
  surfDyyB[1] = (-0.408248290463863*DyyC[3])+0.408248290463863*DyyB[3]+0.3535533905932737*DyyC[1]+0.3535533905932737*DyyB[1]
  surfDyyB[2] = (-0.408248290463863*DyyC[4])+0.408248290463863*DyyB[4]+0.3535533905932737*DyyC[2]+0.3535533905932737*DyyB[2]

  local surfDxxL = {}
  surfDxxL[1] = 0.408248290463863*DxxL[2]-0.408248290463863*DxxC[2]+0.3535533905932737*DxxL[1]+0.3535533905932737*DxxC[1]
  surfDxxL[2] = 0.408248290463863*DxxL[4]-0.408248290463863*DxxC[4]+0.3535533905932737*DxxL[3]+0.3535533905932737*DxxC[3]

  local surfDxxR = {}
  surfDxxR[1] = (-0.408248290463863*DxxR[2])+0.408248290463863*DxxC[2]+0.3535533905932737*DxxR[1]+0.3535533905932737*DxxC[1]
  surfDxxR[2] = (-0.408248290463863*DxxR[4])+0.408248290463863*DxxC[4]+0.3535533905932737*DxxR[3]+0.3535533905932737*DxxC[3]

  local surfDxyT = {}
  surfDxyT[1] = (-0.408248290463863*DxyT[3])+0.408248290463863*DxyC[3]+0.3535533905932737*DxyT[1]+0.3535533905932737*DxyC[1]
  surfDxyT[2] = (-0.408248290463863*DxyT[4])+0.408248290463863*DxyC[4]+0.3535533905932737*DxyT[2]+0.3535533905932737*DxyC[2]

  local surfDxyB = {}
  surfDxyB[1] = (-0.408248290463863*DxyC[3])+0.408248290463863*DxyB[3]+0.3535533905932737*DxyC[1]+0.3535533905932737*DxyB[1]
  surfDxyB[2] = (-0.408248290463863*DxyC[4])+0.408248290463863*DxyB[4]+0.3535533905932737*DxyC[2]+0.3535533905932737*DxyB[2]

  local surfDxyL = {}
  surfDxyL[1] = 0.408248290463863*DxyL[2]-0.408248290463863*DxyC[2]+0.3535533905932737*DxyL[1]+0.3535533905932737*DxyC[1]
  surfDxyL[2] = 0.408248290463863*DxyL[4]-0.408248290463863*DxyC[4]+0.3535533905932737*DxyL[3]+0.3535533905932737*DxyC[3]

  surfDxyR = {}
  surfDxyR[1] = (-0.408248290463863*DxyR[2])+0.408248290463863*DxyC[2]+0.3535533905932737*DxyR[1]+0.3535533905932737*DxyC[1]
  surfDxyR[2] = (-0.408248290463863*DxyR[4])+0.408248290463863*DxyC[4]+0.3535533905932737*DxyR[3]+0.3535533905932737*DxyC[3]

  _M[1] = Lin.Mat(4,4)
  _M[1][1][1] = 0.0
  _M[1][1][2] = 0.0
  _M[1][1][3] = 0.0
  _M[1][1][4] = 0.0
  _M[1][2][1] = 0.0
  _M[1][2][2] = 0.0
  _M[1][2][3] = 0.0
  _M[1][2][4] = 0.0
  _M[1][3][1] = 0.0
  _M[1][3][2] = 0.0
  _M[1][3][3] = 0.0
  _M[1][3][4] = 0.0
  _M[1][4][1] = 0.0
  _M[1][4][2] = 0.0
  _M[1][4][3] = 0.0
  _M[1][4][4] = 0.0
  _M[2] = Lin.Mat(4,4)
  _M[2][1][1] = (1.590990257669731*surfDxxL[1])/dx[1]^2
  _M[2][1][2] = (1.530931089239486*surfDxxL[1])/dx[1]^2
  _M[2][1][3] = (1.590990257669731*surfDxxL[2])/dx[1]^2-(1.224744871391589*surfDxyL[1])/(dx[1]*dx[2])
  _M[2][1][4] = (1.530931089239486*surfDxxL[2])/dx[1]^2-(1.414213562373095*surfDxyL[1])/(dx[1]*dx[2])
  _M[2][2][1] = -(2.755675960631073*surfDxxL[1])/dx[1]^2
  _M[2][2][2] = -(2.651650429449552*surfDxxL[1])/dx[1]^2
  _M[2][2][3] = (2.121320343559642*surfDxyL[1])/(dx[1]*dx[2])-(2.755675960631073*surfDxxL[2])/dx[1]^2
  _M[2][2][4] = (2.449489742783178*surfDxyL[1])/(dx[1]*dx[2])-(2.651650429449552*surfDxxL[2])/dx[1]^2
  _M[2][3][1] = (1.590990257669731*surfDxxL[2])/dx[1]^2
  _M[2][3][2] = (1.530931089239486*surfDxxL[2])/dx[1]^2
  _M[2][3][3] = (1.590990257669731*surfDxxL[1])/dx[1]^2-(1.224744871391589*surfDxyL[2])/(dx[1]*dx[2])
  _M[2][3][4] = (1.530931089239486*surfDxxL[1])/dx[1]^2-(1.414213562373095*surfDxyL[2])/(dx[1]*dx[2])
  _M[2][4][1] = -(2.755675960631073*surfDxxL[2])/dx[1]^2
  _M[2][4][2] = -(2.651650429449552*surfDxxL[2])/dx[1]^2
  _M[2][4][3] = (2.121320343559642*surfDxyL[2])/(dx[1]*dx[2])-(2.755675960631073*surfDxxL[1])/dx[1]^2
  _M[2][4][4] = (2.449489742783178*surfDxyL[2])/(dx[1]*dx[2])-(2.651650429449552*surfDxxL[1])/dx[1]^2
  _M[3] = Lin.Mat(4,4)
  _M[3][1][1] = 0.0
  _M[3][1][2] = 0.0
  _M[3][1][3] = 0.0
  _M[3][1][4] = 0.0
  _M[3][2][1] = 0.0
  _M[3][2][2] = 0.0
  _M[3][2][3] = 0.0
  _M[3][2][4] = 0.0
  _M[3][3][1] = 0.0
  _M[3][3][2] = 0.0
  _M[3][3][3] = 0.0
  _M[3][3][4] = 0.0
  _M[3][4][1] = 0.0
  _M[3][4][2] = 0.0
  _M[3][4][3] = 0.0
  _M[3][4][4] = 0.0
  _M[4] = Lin.Mat(4,4)
  _M[4][1][1] = 0.0
  _M[4][1][2] = 0.0
  _M[4][1][3] = 0.0
  _M[4][1][4] = 0.0
  _M[4][2][1] = 0.0
  _M[4][2][2] = 0.0
  _M[4][2][3] = 0.0
  _M[4][2][4] = 0.0
  _M[4][3][1] = 0.0
  _M[4][3][2] = 0.0
  _M[4][3][3] = 0.0
  _M[4][3][4] = 0.0
  _M[4][4][1] = 0.0
  _M[4][4][2] = 0.0
  _M[4][4][3] = 0.0
  _M[4][4][4] = 0.0
  _M[5] = Lin.Mat(4,4)
  _M[5][1][1] = (-(114.5512985522207*surfDxxR[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))-(114.5512985522207*surfDxxL[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(114.5512985522207*dx[1]^2*surfDyyT[1]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(25.45584412271572*surfDxxR[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(25.45584412271572*surfDxxL[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(25.45584412271572*dx[1]^2*surfDyyT[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(84.85281374238573*dx[1]^2*surfDyyB[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][1][2] = (-(114.5512985522207*dx[1]^2*surfDyyT[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))-(110.227038425243*surfDxxR[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(110.227038425243*surfDxxL[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(88.18163074019442*dx[1]*surfDxyT[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(146.9693845669907*dx[1]*surfDxyB[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(25.45584412271572*dx[1]^2*surfDyyT[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(84.85281374238573*dx[1]^2*surfDyyB[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(24.49489742783179*surfDxxR[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(24.49489742783179*surfDxxL[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(19.59591794226543*dx[1]*surfDxyT[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][1][3] = (-(114.5512985522207*dx[2]^2*surfDxxR[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))-(114.5512985522207*dx[2]^2*surfDxxL[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(88.18163074019442*dx[1]*surfDxyR[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(88.18163074019442*dx[1]*surfDxyL[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(110.227038425243*dx[1]^2*surfDyyT[1]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(25.45584412271572*dx[2]^2*surfDxxR[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(25.45584412271572*dx[2]^2*surfDxxL[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(19.59591794226543*dx[1]*surfDxyR[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(19.59591794226543*dx[1]*surfDxyL[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(24.49489742783179*dx[1]^2*surfDyyT[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(176.3632614803889*dx[1]^2*surfDyyB[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][1][4] = (-(110.227038425243*dx[1]^2*surfDyyT[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))-(110.227038425243*dx[2]^2*surfDxxR[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(110.227038425243*dx[2]^2*surfDxxL[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(101.8233764908629*dx[1]*surfDxyT[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(101.8233764908629*dx[1]*surfDxyR[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(101.8233764908629*dx[1]*surfDxyL[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(305.4701294725887*dx[1]*surfDxyB[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(24.49489742783179*dx[1]^2*surfDyyT[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(176.3632614803889*dx[1]^2*surfDyyB[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(24.49489742783179*dx[2]^2*surfDxxR[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(24.49489742783179*dx[2]^2*surfDxxL[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(22.62741699796953*dx[1]*surfDxyT[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(22.62741699796953*dx[1]*surfDxyR[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(22.62741699796953*dx[1]*surfDxyL[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][2][1] = (-(114.5512985522207*dx[1]^2*surfDyyT[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))-(198.4086691654373*surfDxxR[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(198.4086691654373*surfDxxL[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(25.45584412271572*dx[1]^2*surfDyyT[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(84.85281374238573*dx[1]^2*surfDyyB[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(44.0908153700972*surfDxxR[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(44.0908153700972*surfDxxL[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][2][2] = (88.18163074019442*dx[1]*dx[2]*surfDxyT[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(146.9693845669907*dx[1]*dx[2]*surfDxyB[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(190.9188309203678*surfDxxR[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(190.9188309203678*surfDxxL[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(432.0*DxxC[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(114.5512985522207*dx[1]^2*surfDyyT[1]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(19.59591794226543*dx[1]*dx[2]*surfDxyT[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(42.42640687119286*surfDxxR[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(42.42640687119286*surfDxxL[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*DxxC[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(25.45584412271572*dx[1]^2*surfDyyT[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(84.85281374238573*dx[1]^2*surfDyyB[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][2][3] = (-(110.227038425243*dx[1]^2*surfDyyT[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))-(198.4086691654373*dx[2]^2*surfDxxR[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(198.4086691654373*dx[2]^2*surfDxxL[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(152.7350647362943*dx[1]*surfDxyR[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(152.7350647362943*dx[1]*surfDxyL[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(432.0*DxyC[1]*dx[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(24.49489742783179*dx[1]^2*surfDyyT[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(176.3632614803889*dx[1]^2*surfDyyB[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(44.0908153700972*dx[2]^2*surfDxxR[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(44.0908153700972*dx[2]^2*surfDxxL[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(33.9411254969543*dx[1]*surfDxyR[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(33.9411254969543*dx[1]*surfDxyL[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*DxyC[1]*dx[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][2][4] = (-(432.0*dx[2]^2*DxxC[3]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))+(101.8233764908629*dx[1]*dx[2]*surfDxyT[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(305.4701294725887*dx[1]*dx[2]*surfDxyB[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(190.9188309203678*dx[2]^2*surfDxxR[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(190.9188309203678*dx[2]^2*surfDxxL[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(432.0*dx[1]*DxyC[2]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(176.3632614803889*dx[1]*surfDxyR[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(176.3632614803889*dx[1]*surfDxyL[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(110.227038425243*dx[1]^2*surfDyyT[1]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*dx[2]^2*DxxC[3]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(22.62741699796953*dx[1]*dx[2]*surfDxyT[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(42.42640687119286*dx[2]^2*surfDxxR[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(42.42640687119286*dx[2]^2*surfDxxL[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*dx[1]*DxyC[2]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(39.19183588453087*dx[1]*surfDxyR[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(39.19183588453087*dx[1]*surfDxyL[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(24.49489742783179*dx[1]^2*surfDyyT[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(176.3632614803889*dx[1]^2*surfDyyB[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][3][1] = (-(114.5512985522207*dx[2]^2*surfDxxR[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))-(114.5512985522207*dx[2]^2*surfDxxL[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(198.4086691654373*dx[1]^2*surfDyyT[1]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(25.45584412271572*dx[2]^2*surfDxxR[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(25.45584412271572*dx[2]^2*surfDxxL[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(44.0908153700972*dx[1]^2*surfDyyT[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(146.9693845669907*dx[1]^2*surfDyyB[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][3][2] = (-(198.4086691654373*dx[1]^2*surfDyyT[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))-(110.227038425243*dx[2]^2*surfDxxR[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(110.227038425243*dx[2]^2*surfDxxL[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(152.7350647362943*dx[1]*surfDxyT[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(254.5584412271572*dx[1]*surfDxyB[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(432.0*DxyC[1]*dx[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(44.0908153700972*dx[1]^2*surfDyyT[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(146.9693845669907*dx[1]^2*surfDyyB[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(24.49489742783179*dx[2]^2*surfDxxR[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(24.49489742783179*dx[2]^2*surfDxxL[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(33.9411254969543*dx[1]*surfDxyT[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*DxyC[1]*dx[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][3][3] = (88.18163074019442*dx[1]*dx[2]*surfDxyR[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(88.18163074019442*dx[1]*dx[2]*surfDxyL[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(114.5512985522207*surfDxxR[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(114.5512985522207*surfDxxL[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(190.9188309203678*dx[1]^2*surfDyyT[1]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(432.0*DyyC[1]*dx[1]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(19.59591794226543*dx[1]*dx[2]*surfDxyR[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(19.59591794226543*dx[1]*dx[2]*surfDxyL[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(25.45584412271572*surfDxxR[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(25.45584412271572*surfDxxL[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(42.42640687119286*dx[1]^2*surfDyyT[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(305.4701294725887*dx[1]^2*surfDyyB[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*DyyC[1]*dx[1]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][3][4] = (-(432.0*dx[1]*dx[2]*DxyC[3]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))-(190.9188309203678*dx[1]^2*surfDyyT[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(101.8233764908629*dx[1]*dx[2]*surfDxyR[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(101.8233764908629*dx[1]*dx[2]*surfDxyL[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(110.227038425243*surfDxxR[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(110.227038425243*surfDxxL[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(176.3632614803889*dx[1]*surfDxyT[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(529.0897844411666*dx[1]*surfDxyB[1]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(432.0*dx[1]^2*DyyC[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*dx[1]*dx[2]*DxyC[3]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(42.42640687119286*dx[1]^2*surfDyyT[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(305.4701294725887*dx[1]^2*surfDyyB[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(22.62741699796953*dx[1]*dx[2]*surfDxyR[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(22.62741699796953*dx[1]*dx[2]*surfDxyL[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(24.49489742783179*surfDxxR[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(24.49489742783179*surfDxxL[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(39.19183588453087*dx[1]*surfDxyT[1]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*dx[1]^2*DyyC[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][4][1] = (-(198.4086691654373*dx[1]^2*surfDyyT[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))-(198.4086691654373*dx[2]^2*surfDxxR[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(198.4086691654373*dx[2]^2*surfDxxL[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(44.0908153700972*dx[1]^2*surfDyyT[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(146.9693845669907*dx[1]^2*surfDyyB[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(44.0908153700972*dx[2]^2*surfDxxR[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(44.0908153700972*dx[2]^2*surfDxxL[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][4][2] = (-(432.0*dx[2]^2*DxxC[3]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))+(152.7350647362943*dx[1]*dx[2]*surfDxyT[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(254.5584412271572*dx[1]*dx[2]*surfDxyB[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(190.9188309203678*dx[2]^2*surfDxxR[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(190.9188309203678*dx[2]^2*surfDxxL[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(432.0*dx[1]*DxyC[2]*dx[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(198.4086691654373*dx[1]^2*surfDyyT[1]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*dx[2]^2*DxxC[3]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(33.9411254969543*dx[1]*dx[2]*surfDxyT[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(42.42640687119286*dx[2]^2*surfDxxR[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(42.42640687119286*dx[2]^2*surfDxxL[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*dx[1]*DxyC[2]*dx[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(44.0908153700972*dx[1]^2*surfDyyT[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(146.9693845669907*dx[1]^2*surfDyyB[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][4][3] = (-(432.0*dx[1]*dx[2]*DxyC[3]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))-(190.9188309203678*dx[1]^2*surfDyyT[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(152.7350647362943*dx[1]*dx[2]*surfDxyR[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(152.7350647362943*dx[1]*dx[2]*surfDxyL[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(198.4086691654373*surfDxxR[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(198.4086691654373*surfDxxL[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(432.0*dx[1]^2*DyyC[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*dx[1]*dx[2]*DxyC[3]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(42.42640687119286*dx[1]^2*surfDyyT[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(305.4701294725887*dx[1]^2*surfDyyB[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(33.9411254969543*dx[1]*dx[2]*surfDxyR[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(33.9411254969543*dx[1]*dx[2]*surfDxyL[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(44.0908153700972*surfDxxR[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(44.0908153700972*surfDxxL[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*dx[1]^2*DyyC[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[5][4][4] = (-(864.0*dx[1]*dx[2]*DxyC[4]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy))+(176.3632614803889*dx[1]*dx[2]*surfDxyT[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(176.3632614803889*dx[1]*dx[2]*surfDxyR[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(176.3632614803889*dx[1]*dx[2]*surfDxyL[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(529.0897844411666*dx[1]*dx[2]*surfDxyB[2]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(190.9188309203678*surfDxxR[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(190.9188309203678*surfDxxL[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(432.0*DxxC[1]*dx[2]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(190.9188309203678*dx[1]^2*surfDyyT[1]*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(432.0*DyyC[1]*dx[1]^2*bcNy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(192.0*dx[1]*dx[2]*DxyC[4]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(39.19183588453087*dx[1]*dx[2]*surfDxyT[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)-(39.19183588453087*dx[1]*dx[2]*surfDxyR[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(39.19183588453087*dx[1]*dx[2]*surfDxyL[2]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(42.42640687119286*surfDxxR[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(42.42640687119286*surfDxxL[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*DxxC[1]*dx[2]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(42.42640687119286*dx[1]^2*surfDyyT[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(305.4701294725887*dx[1]^2*surfDyyB[1]*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)+(96.0*DyyC[1]*dx[1]^2*bcDy)/(72.0*dx[1]^2*dx[2]^2*bcNy-16.0*dx[1]^2*dx[2]^2*bcDy)
  _M[6] = Lin.Mat(4,4)
  _M[6][1][1] = (114.5512985522207*surfDyyT[1]*bcNy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)-(25.45584412271572*surfDyyT[1]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)+(16.97056274847715*surfDyyB[1]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)
  _M[6][1][2] = (114.5512985522207*dx[1]*surfDyyT[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(88.18163074019442*surfDxyT[1]*dx[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(29.39387691339815*surfDxyB[1]*dx[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(25.45584412271572*dx[1]*surfDyyT[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(16.97056274847715*dx[1]*surfDyyB[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(19.59591794226543*surfDxyT[1]*dx[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)
  _M[6][1][3] = (-(110.227038425243*surfDyyT[1]*bcNy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy))+(24.49489742783179*surfDyyT[1]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)-(19.59591794226543*surfDyyB[1]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)
  _M[6][1][4] = (-(110.227038425243*dx[1]*surfDyyT[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy))-(101.8233764908629*surfDxyT[1]*dx[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(33.9411254969543*surfDxyB[1]*dx[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(24.49489742783179*dx[1]*surfDyyT[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(19.59591794226543*dx[1]*surfDyyB[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(22.62741699796953*surfDxyT[1]*dx[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)
  _M[6][2][1] = (114.5512985522207*surfDyyT[2]*bcNy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)-(25.45584412271572*surfDyyT[2]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)+(16.97056274847715*surfDyyB[2]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)
  _M[6][2][2] = (88.18163074019442*dx[2]*surfDxyT[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(29.39387691339815*dx[2]*surfDxyB[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(114.5512985522207*dx[1]*surfDyyT[1]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(19.59591794226543*dx[2]*surfDxyT[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(25.45584412271572*dx[1]*surfDyyT[1]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(16.97056274847715*dx[1]*surfDyyB[1]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)
  _M[6][2][3] = (-(110.227038425243*surfDyyT[2]*bcNy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy))+(24.49489742783179*surfDyyT[2]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)-(19.59591794226543*surfDyyB[2]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)
  _M[6][2][4] = (-(101.8233764908629*dx[2]*surfDxyT[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy))+(33.9411254969543*dx[2]*surfDxyB[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(110.227038425243*dx[1]*surfDyyT[1]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(22.62741699796953*dx[2]*surfDxyT[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(24.49489742783179*dx[1]*surfDyyT[1]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(19.59591794226543*dx[1]*surfDyyB[1]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)
  _M[6][3][1] = (198.4086691654373*surfDyyT[1]*bcNy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)-(44.0908153700972*surfDyyT[1]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)-(29.39387691339815*surfDyyB[1]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)
  _M[6][3][2] = (198.4086691654373*dx[1]*surfDyyT[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(152.7350647362943*surfDxyT[1]*dx[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(50.91168824543144*surfDxyB[1]*dx[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(44.0908153700972*dx[1]*surfDyyT[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(29.39387691339815*dx[1]*surfDyyB[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(33.9411254969543*surfDxyT[1]*dx[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)
  _M[6][3][3] = (-(190.9188309203678*surfDyyT[1]*bcNy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy))+(42.42640687119286*surfDyyT[1]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)+(33.9411254969543*surfDyyB[1]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)
  _M[6][3][4] = (-(190.9188309203678*dx[1]*surfDyyT[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy))-(176.3632614803889*surfDxyT[1]*dx[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(58.7877538267963*surfDxyB[1]*dx[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(42.42640687119286*dx[1]*surfDyyT[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(33.9411254969543*dx[1]*surfDyyB[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(39.19183588453087*surfDxyT[1]*dx[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)
  _M[6][4][1] = (198.4086691654373*surfDyyT[2]*bcNy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)-(44.0908153700972*surfDyyT[2]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)-(29.39387691339815*surfDyyB[2]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)
  _M[6][4][2] = (152.7350647362943*dx[2]*surfDxyT[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(50.91168824543144*dx[2]*surfDxyB[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(198.4086691654373*dx[1]*surfDyyT[1]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(33.9411254969543*dx[2]*surfDxyT[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(44.0908153700972*dx[1]*surfDyyT[1]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(29.39387691339815*dx[1]*surfDyyB[1]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)
  _M[6][4][3] = (-(190.9188309203678*surfDyyT[2]*bcNy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy))+(42.42640687119286*surfDyyT[2]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)+(33.9411254969543*surfDyyB[2]*bcDy)/(72.0*dx[2]^2*bcNy-16.0*dx[2]^2*bcDy)
  _M[6][4][4] = (-(176.3632614803889*dx[2]*surfDxyT[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy))-(58.7877538267963*dx[2]*surfDxyB[2]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)-(190.9188309203678*dx[1]*surfDyyT[1]*bcNy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(39.19183588453087*dx[2]*surfDxyT[2]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(42.42640687119286*dx[1]*surfDyyT[1]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)+(33.9411254969543*dx[1]*surfDyyB[1]*bcDy)/(72.0*dx[1]*dx[2]^2*bcNy-16.0*dx[1]*dx[2]^2*bcDy)
  _M[7] = Lin.Mat(4,4)
  _M[7][1][1] = 0.0
  _M[7][1][2] = 0.0
  _M[7][1][3] = 0.0
  _M[7][1][4] = 0.0
  _M[7][2][1] = 0.0
  _M[7][2][2] = 0.0
  _M[7][2][3] = 0.0
  _M[7][2][4] = 0.0
  _M[7][3][1] = 0.0
  _M[7][3][2] = 0.0
  _M[7][3][3] = 0.0
  _M[7][3][4] = 0.0
  _M[7][4][1] = 0.0
  _M[7][4][2] = 0.0
  _M[7][4][3] = 0.0
  _M[7][4][4] = 0.0
  _M[8] = Lin.Mat(4,4)
  _M[8][1][1] = (1.590990257669731*surfDxxR[1])/dx[1]^2
  _M[8][1][2] = -(1.530931089239486*surfDxxR[1])/dx[1]^2
  _M[8][1][3] = (1.590990257669731*surfDxxR[2])/dx[1]^2+(1.224744871391589*surfDxyR[1])/(dx[1]*dx[2])
  _M[8][1][4] = (-(1.530931089239486*surfDxxR[2])/dx[1]^2)-(1.414213562373095*surfDxyR[1])/(dx[1]*dx[2])
  _M[8][2][1] = (2.755675960631073*surfDxxR[1])/dx[1]^2
  _M[8][2][2] = -(2.651650429449552*surfDxxR[1])/dx[1]^2
  _M[8][2][3] = (2.755675960631073*surfDxxR[2])/dx[1]^2+(2.121320343559642*surfDxyR[1])/(dx[1]*dx[2])
  _M[8][2][4] = (-(2.651650429449552*surfDxxR[2])/dx[1]^2)-(2.449489742783178*surfDxyR[1])/(dx[1]*dx[2])
  _M[8][3][1] = (1.590990257669731*surfDxxR[2])/dx[1]^2
  _M[8][3][2] = -(1.530931089239486*surfDxxR[2])/dx[1]^2
  _M[8][3][3] = (1.224744871391589*surfDxyR[2])/(dx[1]*dx[2])+(1.590990257669731*surfDxxR[1])/dx[1]^2
  _M[8][3][4] = (-(1.414213562373095*surfDxyR[2])/(dx[1]*dx[2]))-(1.530931089239486*surfDxxR[1])/dx[1]^2
  _M[8][4][1] = (2.755675960631073*surfDxxR[2])/dx[1]^2
  _M[8][4][2] = -(2.651650429449552*surfDxxR[2])/dx[1]^2
  _M[8][4][3] = (2.121320343559642*surfDxyR[2])/(dx[1]*dx[2])+(2.755675960631073*surfDxxR[1])/dx[1]^2
  _M[8][4][4] = (-(2.449489742783178*surfDxyR[2])/(dx[1]*dx[2]))-(2.651650429449552*surfDxxR[1])/dx[1]^2
  _M[9] = Lin.Mat(4,4)
  _M[9][1][1] = 0.0
  _M[9][1][2] = 0.0
  _M[9][1][3] = 0.0
  _M[9][1][4] = 0.0
  _M[9][2][1] = 0.0
  _M[9][2][2] = 0.0
  _M[9][2][3] = 0.0
  _M[9][2][4] = 0.0
  _M[9][3][1] = 0.0
  _M[9][3][2] = 0.0
  _M[9][3][3] = 0.0
  _M[9][3][4] = 0.0
  _M[9][4][1] = 0.0
  _M[9][4][2] = 0.0
  _M[9][4][3] = 0.0
  _M[9][4][4] = 0.0
  _M[10] = Lin.Vec(4)
  _M[10][1] = 0.0 + -(50.91168824543144*surfDyyB[1]*bcValy)/(18.0*dx[2]^2*bcNy-4.0*dx[2]^2*bcDy)
  _M[10][2] = 0.0 + -(50.91168824543144*surfDyyB[2]*bcValy)/(18.0*dx[2]^2*bcNy-4.0*dx[2]^2*bcDy)
  _M[10][3] = 0.0 + (88.18163074019442*surfDyyB[1]*bcValy)/(18.0*dx[2]^2*bcNy-4.0*dx[2]^2*bcDy)
  _M[10][4] = 0.0 + (88.18163074019442*surfDyyB[2]*bcValy)/(18.0*dx[2]^2*bcNy-4.0*dx[2]^2*bcDy)
  return(_M)
end

return(stencilFn)