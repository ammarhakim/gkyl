-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments()

local gasGamma = 5.0/3.0
local lower = -10.0
local upper = 10.0
local cells = 16

local epsilon0 = 1.0
local mu0 = 1.0
local chargeElc = -1.0
local massElc = 1.0
local chargeIon = 1.0
local massIon = 1836.153
local n0 = 1.0
local elcTemp = 1.0/2048.0 -- vte/c = 1/32
local ionTemp = elcTemp
-- using plasma beta definition beta = vt^2/vA^2
-- plasma beta = 1/1024
local B0 = 1.0

-- NOTE: with these parameters, rhoe = vte/OmegaCe = 1/32
--                              rhoi = vti/OmegaCi ~ 1.33
--                              Lx = 2*320 rhoe ~ 2*7.5 rhoi

-- constant electric field for uniform E x B drift
local E0 = 0.1

local cfl = 1.0
local tEnd = 200.0*math.pi
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

   -- boundary conditions for configuration space
   periodicDirs = {1, 2}, -- periodic directions

   -- Electrons.
   elc = Moments.Species {
      charge = chargeElc, mass = massElc,

      equation    = Moments.Euler { gasGamma = gasGamma },
      -- Initial conditions.
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local rhoe = massElc
         local rhovx_e = massElc*E0/B0
         local rhovy_e = -massElc*E0/B0
         local rhovz_e = 0.0
         local u_e = rhoe*elcTemp/(gasGamma - 1) + 0.5*(rhovx_e*rhovx_e + rhovy_e*rhovy_e + rhovz_e*rhovz_e)/rhoe
         return rhoe, rhovx_e, rhovy_e, rhovz_e, u_e
      end,
   },

   -- Ions.
   ion = Moments.Species {
      charge = chargeIon, mass = massIon,

      equation    = Moments.Euler { gasGamma = gasGamma },
      -- Initial conditions.
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local rhoi = massIon
         local rhovx_i = massIon*E0/B0
         local rhovy_i = -massIon*E0/B0
         local rhovz_i = 0.0
         local u_i = rhoi*ionTemp/(gasGamma - 1) + 0.5*(rhovx_i*rhovx_i + rhovy_i*rhovy_i + rhovz_i*rhovz_i)/rhoi
         return rhoi, rhovx_i, rhovy_i, rhovz_i, u_i
      end,
   },

   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      is_ext_em_static = true,
      ext_em_func = function (t, xn)
         local x, y = xn[1], xn[2]
         return E0, E0, 0.0, 0.0, 0.0, B0
      end,
   },
}

-- Run application.
momentApp:run()
