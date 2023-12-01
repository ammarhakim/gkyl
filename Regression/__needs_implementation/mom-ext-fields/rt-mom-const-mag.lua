-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments()

local gasGamma = 5.0/3.0
local lower = -1.0
local upper = 1.0
local cells = 64

local epsilon0 = 1.0
local mu0 = 1.0
local chargeElc = -1.0
local massElc = 1.0
local n0 = 1.0
local elcTemp = 1.0/2048.0 -- vte/c = 1/32
-- using plasma beta definition beta = vt^2/vA^2
-- plasma beta = 1/1024
local B0 = 1.0

-- NOTE: with these parameters, rhoe = vte/OmegaCe = 1/32

local cfl = 1.0
local tEnd = 4.0*math.pi
local nFrame = 1

-- Gaussian pulse
local function pulse2D(n, x, y, x0, y0, sigma)
   local v2 = (x - x0)^2 + (y - y0)^2
   return n/(2*math.pi*sigma^2)*math.exp(-v2/(2*sigma^2))
end

momentApp = Moments.App {
   tEnd   = tEnd,
   nFrame = nFrame,
   lower  = {lower, lower},
   upper  = {upper, upper},
   cells  = {cells, cells},
   cfl    = cfl,

   -- Decomposition for configuration space.
   decompCuts = {1, 1},    -- Cuts in each configuration direction.

   -- Electrons.
   elc = Moments.Species {
      charge = chargeElc, mass = massElc,

      equation    = Moments.Euler { gasGamma = gasGamma },
      -- Initial conditions.
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local rhoe = pulse2D(massElc, x, y, 0.0, 0.0, 0.25)
         local rhovx_e = 0.0
         local rhovy_e = 0.0
         local rhovz_e = 0.0
         local u_e = rhoe*elcTemp/(gasGamma - 1) + 0.5*(rhovx_e*rhovx_e + rhovy_e*rhovy_e + rhovz_e*rhovz_e)/rhoe
         return rhoe, rhovx_e, rhovy_e, rhovz_e, u_e
      end,
      bcx = { Moments.Species.bcWall, Moments.Species.bcWall },
      bcy = { Moments.Species.bcWall, Moments.Species.bcWall },
   },


   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      bcx = { Moments.Field.bcReflect, Moments.Field.bcReflect },
      bcy = { Moments.Field.bcReflect, Moments.Field.bcReflect },

      is_ext_em_static = true,
      ext_em_func = function (t, xn)
         local x, y = xn[1], xn[2]
         return 0.0, 0.0, 0.0, 0.0, 0.0, B0
      end
   },

}

-- Run application.
momentApp:run()
