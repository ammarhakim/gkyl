-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments()

gasGamma = 1.4 -- gas adiabatic constant

-- inflow conditions
rhoIn = 1.0
prIn = 1.0
csIn = math.sqrt(gasGamma*prIn/rhoIn)
uIn = 8.0*csIn -- supersonic inflow
erIn = 0.5*rhoIn*uIn^2 + prIn/(gasGamma-1)

function wedge(x, y, x0, y0, angle)
   return (
      x>x0 and (y<y0+(x-x0)*math.tan(0.5*angle)) and (y>y0-(x-x0)*math.tan(0.5*angle))
   ) and true or false
end

print("rhoIn", rhoIn)

function bc_lower_func(t, nc, skin, ctx)
   return rhoIn, rhoIn*uIn, 0.0, 0.0, erIn
end

-- create app
eulerApp = Moments.App {
   tEnd = 5*2.0/uIn, -- end time
   nFrame = 10, -- number of output frame
   lower = {-0.1, 0.0}, -- lower left corner
   upper = {1.0, 0.55}, -- upper right corner
   cells = {200, 100}, -- number of cells
   cflFrac = 0.9, -- CFL fraction
   
   -- decomposition stuff
   decompCuts = {1, 1}, -- cuts in each direction

   -- electrons
   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,

      equation = Moments.Euler { gasGamma = gasGamma },
      -- initial conditions
      init = function (t, xn)
         return 1e-6*rhoIn, 0.0, 0.0, 0.0, 1e-6*erIn
      end,
      evolve = true, -- evolve species?

      -- outer boundary conditions
      bcx = { Moments.Species.bcFunc, Moments.Species.bcCopy },
      bcy = { Moments.Species.bcWall, Moments.Species.bcCopy },
      bcx_lower_func = bc_lower_func,

      -- has interior (embedded) boundary conditions?
      hasSsBnd = true,
      -- mask defining the embedded boundary
      inOutFunc = function (t, xn)
         local x, y = xn[1], xn[2]
         local xc, yc = 0.0, 0.0
         local rad = 0.5
         return  wedge(x,y, 0.0, 0.0, 30*math.pi/180) and -1.0 or 1.0
      end,
      -- boundary conditions to be applied on the embedded boundary
      ssBc = { Moments.Species.bcWall },
   },   
}

-- run application
eulerApp:run()
