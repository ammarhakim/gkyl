-- Gkyl ------------------------------------------------------------------------
--
-- Fix M0 by rescaling, and fix M1 and M2 by iteration.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"

-- .............. USER INPUTS ................ --

-- Basic parameters
local n0   = 1.0                             -- Density.
local u0   = 0.0                             -- Flow speed.
local vt   = 1.0/3.0                         -- Thermal speed..
local B0   = 1.0                             -- Background magnetic field
-- Large bump on tail of Maxwellian:
local ab   = math.sqrt(0.1)                  -- Amplitude of bump.
local ub   = 4*math.sqrt(((3*vt/2)^2)/3)     -- Location of bump.
local sb   = 0.12                            -- Softening factor to avoid divergence.
local vtb  = 1.0                             -- Thermal speed of Maxwellian in bump.  

-- Initial distribution function: maxwellian with bump on tail
local function bumpMaxwell(x,vx,n,u,vth,bA,bU,bS,bVth)
   local Pi   = math.pi
   local vSq  = ((vx-u)/(math.sqrt(2.0)*vth))^2
   local vbSq = ((vx-u)/(math.sqrt(2.0)*bVth))^2

   return (n/math.sqrt(2.0*Pi*vth))*math.exp(-vSq)
         +(n/math.sqrt(2.0*Pi*bVth))*math.exp(-vbSq)*(bA^2)/((vx-bU)^2+bS^2)
end
local initialF = function (t, xn)
   local x, v = xn[1], xn[2]
   return bumpMaxwell(x,v,n0,u0,vt,ab,ub,sb,vtb)
end

-- Phase space mesh.
local lower     = {0.0, -8.0*vt}
local upper     = {1.0,  8.0*vt}
local numCells  = {2, 32}

local polyOrder = 1 -- Polynomial order of the basis.
local cDim = 1  -- Number of configuration space dimensions.

-- Evolve and output time
local tEnd   = 100  -- End time
local nFrame = 50   -- Number of frames to write

-- .............. END OF USER INPUTS (maybe) ................ --

local function createGrid(lo, up, nCells, pDirs)
   pDirs = pDirs or {}
   local gridOut = Grid.RectCart {
      lower = lo,  cells        = nCells,
      upper = up,  periodicDirs = pDirs,
   }
   return gridOut
end

local function createConfBasis(dim, pOrder, bKind)
   bKind = bKind or "Ser"
   local basis
   if (bKind=="Ser") then
      basis = Basis.CartModalSerendipity { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Tensor") then
      basis = Basis.CartModalTensor { ndim = dim, polyOrder = pOrder }
   else
      assert(false,"Invalid basis")
   end
   return basis
end

local function createPhaseBasis(cdim, vdim, pOrder, bKind)
   bKind = bKind or "Ser"
   local pdim = cdim+vdim
   local basis
   if (bKind=="Ser") then
      basis = Basis.CartModalSerendipity { ndim = pdim, polyOrder = pOrder }
   elseif (bKind=="Tensor") then
      basis = Basis.CartModalTensor { ndim = pdim, polyOrder = pOrder }
   else
      assert(false,"Invalid basis")
   end
   return basis
end

local function createField(grid, basis, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = {1, 1},
      metaData      = {polyOrder  = basis:polyOrder(),
                       basisType  = basis:id(),},
   }
   return fld
end

local pDim = #lower   -- Length of lower
local vDim = pDim-cDim

local lowerC, upperC, numCellsC = {}, {}, {}
for d=1,cDim do
   lowerC[d] = lower[d]
   upperC[d] = upper[d]
   numCellsC[d] = numCells[d]
end

-- Create phase and conf space grids and bases.
local phaseGrid  = createGrid(lower, upper, numCells)
local phaseBasis = createPhaseBasis(cDim, vDim, polyOrder)
local confGrid   = createGrid(lowerC, upperC, numCellsC)
local confBasis  = createConfBasis(confGrid:ndim(), polyOrder)

-- Create various fields needed.
local distf  = createField(phaseGrid, phaseBasis)   -- create a field
local fM     = createField(phaseGrid, phaseBasis)
local m0     = createField(confGrid, confBasis) 
local m0scl  = createField(confGrid, confBasis) 
local dm0    = createField(confGrid, confBasis) 
local ddm0   = createField(confGrid, confBasis) 
local m0_new = createField(confGrid, confBasis)
local m1     = createField(confGrid, confBasis, vDim)
local dm1    = createField(confGrid, confBasis, vDim) 
local ddm1   = createField(confGrid, confBasis, vDim) 
local m1_new = createField(confGrid, confBasis, vDim)
local m2     = createField(confGrid, confBasis)
local dm2    = createField(confGrid, confBasis)
local ddm2   = createField(confGrid, confBasis)
local m2_new = createField(confGrid, confBasis)
local m2flow = createField(confGrid, confBasis)
local uDrift = createField(confGrid, confBasis, vDim)
local vtSq   = createField(confGrid, confBasis)

-- Create projection updaters.
local projPhaseScalar = Updater.ProjectOnBasis {
   onGrid   = phaseGrid,
   basis    = phaseBasis,
   evaluate = function (t, xn) return 1.0 end   -- Set later.
}
local projConfScalar = Updater.ProjectOnBasis {
   onGrid   = confGrid,
   basis    = confBasis,
   evaluate = function (t, xn) return 1.0 end   -- Set later.
}
local maxwellian = Updater.MaxwellianOnBasis {
   onGrid     = phaseGrid,   confGrid  = confGrid,
   phaseBasis = phaseBasis,  confBasis = confBasis,
}

-- Moment updaters.
local calcNumDensity = Updater.DistFuncMomentCalc {
   onGrid     = phaseGrid,   confBasis = confBasis,
   phaseBasis = phaseBasis,  moment    = "M0",
}
local calcMomDensity = Updater.DistFuncMomentCalc {
   onGrid     = phaseGrid,   confBasis = confBasis,
   phaseBasis = phaseBasis,  moment    = "M1i",
}
local calcKinEnergyDensity = Updater.DistFuncMomentCalc {
   onGrid     = phaseGrid,   confBasis = confBasis,
   phaseBasis = phaseBasis,  moment    = "M2",
}

-- Integrated moments updaters.
local intm0 = DataStruct.DynVector {
      numComponents = 1,
}
local intm1= DataStruct.DynVector {
      numComponents = 1,
}
local intm2 = DataStruct.DynVector {
      numComponents = 1,
}
local intQuantC = Updater.CartFieldIntegratedQuantCalc {
   onGrid        = confGrid,
   basis         = confBasis,
   numComponents = 1,
   quantity      = "V",
}

-- Binary operation updaters.
local weakDivide = Updater.CartFieldBinOp {
   onGrid     = confGrid,
   weakBasis  = confBasis,
   operation  = "Divide",
}
local weakMultiply = Updater.CartFieldBinOp {
   onGrid     = confGrid,
   weakBasis  = confBasis,
   operation  = "Multiply",
}
local weakMultiplyConfPhase = Updater.CartFieldBinOp {
   onGrid     = phaseGrid,
   weakBasis  = phaseBasis,
   operation  = "Multiply",
   fieldBasis = confBasis,
}

-- Project the initial distribution function.
projPhaseScalar:setFunc(initialF)
projPhaseScalar:advance(0.0, {}, {distf})

-- Compute the initial M0, M1i and M2, and write them out.
calcNumDensity:advance(0.0, {distf}, {m0})
calcMomDensity:advance(0.0, {distf}, {m1})
calcKinEnergyDensity:advance(0.0, {distf}, {m2})
intQuantC:advance(0.0, {m0}, {intm0})
intQuantC:advance(0.0, {m1}, {intm1})
intQuantC:advance(0.0, {m2}, {intm2})
m0:write("m0_in.bp")
intm0:write("intm0_in.bp")
m1:write("m1_in.bp")
intm1:write("intm1_in.bp")
m2:write("m2_in.bp")
intm2:write("intm2_in.bp")

-- Compute u and vt^2.
weakDivide:advance(0., {m0, m1}, {uDrift})
weakMultiply:advance(0., {uDrift, m1}, {m2flow})
vtSq:combine(1., m2, -1., m2flow)
weakDivide:advance(0., {m0, vtSq}, {vtSq})

-- Project the Maxwellian onto the basis.
maxwellian:advance(0.0, {m0,uDrift,vtSq}, {fM})
fM:write("0.bp")

-- Rescale the Maxwellian.
calcNumDensity:advance(0.0, {fM}, {m0_new})
weakDivide:advance(0., {m0_new, m0}, {m0scl})
weakMultiplyConfPhase:advance(0., {m0scl, fM}, {fM})
-- Compute the moments.
calcNumDensity:advance(0.0, {fM}, {m0_new})
calcMomDensity:advance(0.0, {fM}, {m1_new})
calcKinEnergyDensity:advance(0.0, {fM}, {m2_new})
intQuantC:advance(0.0, {m0_new}, {intm0})
intQuantC:advance(0.0, {m1_new}, {intm1})
intQuantC:advance(0.0, {m2_new}, {intm2}) 

m0_new:write("m0_0.bp")
intm0:write("intm0_0.bp")
m1_new:write("m1_0.bp")
intm1:write("intm1_0.bp")
m2_new:write("m2_0.bp")
intm2:write("intm2_0.bp")
  
-- Fix the moments with rescaling and iteration
for n=1,tEnd do
   -- Fix m1.
   ddm1:combine(1., m1, -1., m1_new)
   dm1:combine(1., dm1, 1., ddm1)
   m1_new:combine(1., m1, 1., dm1)
   -- Fix m2.
   ddm2:combine(1., m2, -1., m2_new)
   dm2:combine(1., dm2, 1., ddm2) 
   m2_new:combine(1., m2, 1., dm2)
   
   -- Compute u and vt^2 with the moments of the rescaled Maxwellian.
   weakDivide:advance(0., {m0_new, m1_new}, {uDrift})
   weakMultiply:advance(0., {uDrift, m1_new}, {m2flow})
   vtSq:combine(1., m2_new, -1., m2flow)
   weakDivide:advance(0., {m0_new, vtSq}, {vtSq})
   
   -- Project the Maxwellian onto the basis.
   maxwellian:advance(0.0, {m0_new,uDrift,vtSq}, {fM})
   -- Rescale the Maxwellian.
   calcNumDensity:advance(0.0, {fM}, {m0_new})
   weakDivide:advance(0., {m0_new, m0}, {m0scl})
   weakMultiplyConfPhase:advance(0., {m0scl, fM}, {fM})
   -- Compute the moments.
   calcNumDensity:advance(0.0, {fM}, {m0_new})
   calcMomDensity:advance(0.0, {fM}, {m1_new})
   calcKinEnergyDensity:advance(0.0, {fM}, {m2_new})

   -- Output the results
   if n%nFrame==0 then
      intQuantC:advance(0.0, {m0_new}, {intm0})
      intQuantC:advance(0.0, {m1_new}, {intm1})
      intQuantC:advance(0.0, {m2_new}, {intm2}) 
      filename = string.format("%s.bp",n/nFrame) 
      m0_new:write("m0_"..filename)
      intm0:write("intm0_"..filename)
      m1_new:write("m1_"..filename)
      intm1:write("intm1_"..filename)
      m2_new:write("m2_"..filename)
      intm2:write("intm2_"..filename)
      fM:write(filename)
   end
end
