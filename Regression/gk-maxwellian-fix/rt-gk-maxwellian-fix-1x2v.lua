-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to fix GK maxwellian with iteration method, 1x2v.
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"
local Time = require "Lib.Time"

-- Reference density, drift velocity and thermal speed.
local n0 = 1.0
local u0 = 0.50
local vt0 = 1.0
local B0 = 0.5
local mass = 1.0

-- Functions defining the density, drift velocity and thermal speed squared.
local eps_n = 0.0
local eps_u = 0.0
local eps_vt = 0.0
local m0Func = function (t, xn) return n0*(1.0+eps_n*math.cos(2.*math.pi*xn[1])) end
local uDriftFunc = function (t, xn) return u0*(1.0+eps_u*math.cos(2.*math.pi*xn[1])) end
local vtSqFunc = function (t, xn) return (vt0*(1.0+eps_vt*math.cos(2.*math.pi*xn[1])))^2 end
local bmagFunc = function (t, xn) return B0 end


-- Initial distribution function
local initialF = function(t, xn)
   local x, vpar, mu = xn[1], xn[2], xn[3]
   local fOut = (m0Func(t,xn)/(math.sqrt(2.*math.pi*vtSqFunc(t,xn))^3))*math.exp(-((vpar-uDriftFunc(t,xn))^2+2*mu*bmagFunc(t,xn)/mass)/(2.*(vtSqFunc(t,xn))))
   return fOut
end

-- Phase space mesh.
local lower = {-0.50, -6.0*vt0, 0.0}
local upper = { 0.50,  6.0*vt0, mass*(6.0*(vt0)^2)/(2*B0)}
local numCells = {8, 8, 8}
-- Dimensions.
local pDim = #numCells
local vDim = 3 

local polyOrder = 1 -- Polynomial order of the basis.

local function createGrid(lo, up, nCells, pDirs)
   pDirs = pDirs or {}
   local gridOut = Grid.RectCart {
      lower = lo,
      cells = nCells,
      upper = up,
      periodicDirs = pDirs,
   }
   return gridOut
end

local function createBasis(dim, pOrder, bKind)
   bKind = bKind or "Ser"
   local basis
   if (bKind=="Ser") then
      basis = Basis.CartModalSerendipity { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Max") then
      basis = Basis.CartModalOrder { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Tensor") then
      basis = Basis.CartModalTensor { ndim = dim, polyOrder = pOrder }
   else
      assert(false,"Invalid basis")
   end
   return basis
end

local function createField(grid, basis, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid = grid,
      numComponents = basis:numBasis()*vComp,
      ghost = {1, 1},
      metaData = {
	 polyOrder  = basis:polyOrder(),
	 basisType  = basis:id(),
      },
   }
   return fld
end

-- Create phase and conf space grids and bases.
local phaseGrid = createGrid(lower, upper, numCells)
local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
local confGrid = createGrid({lower[1]}, {upper[1]}, {numCells[1]})
local confBasis = createBasis(confGrid:ndim(), polyOrder)

-- Create various fields needed.
local distf = createField(phaseGrid, phaseBasis)
local fM = createField(phaseGrid, phaseBasis)
local fOut = createField(phaseGrid, phaseBasis)

local m0 = createField(confGrid, confBasis) 
local m1 = createField(confGrid, confBasis)
local m2 = createField(confGrid, confBasis)
local m2flow = createField(confGrid, confBasis)
local uDrift = createField(confGrid, confBasis)
local vtSq = createField(confGrid, confBasis)
local bmag = createField(confGrid, confBasis)
local m0_ic = createField(confGrid, confBasis) 
local m1_ic = createField(confGrid, confBasis)
local m2_ic = createField(confGrid, confBasis)

-- Create projection updaters.
local projPhaseScalar = Updater.ProjectOnBasis {
   onGrid = phaseGrid,
   basis = phaseBasis,
   evaluate = function (t, xn) return 1.0 end   -- Set later.
}
local projConfScalar = Updater.ProjectOnBasis {
   onGrid = confGrid,
   basis = confBasis,
   evaluate = function (t, xn) return 1.0 end   -- Set later.
}
local maxwellian = Updater.MaxwellianOnBasis {
   onGrid = phaseGrid,
   confGrid = confGrid,
   phaseBasis = phaseBasis,
   confBasis = confBasis,
   mass = mass,
}

-- Moment updaters.
local calcNumDensity = Updater.DistFuncMomentCalc {
   onGrid = phaseGrid,
   confBasis = confBasis,
   phaseBasis = phaseBasis,
   moment = "GkM0",
   gkfacs = {mass, bmag},
}
local calcMomDensity = Updater.DistFuncMomentCalc {
   onGrid = phaseGrid,
   confBasis = confBasis,
   phaseBasis = phaseBasis,
   moment = "GkM1",
   gkfacs = {mass, bmag},
}
local calcKinEnergyDensity = Updater.DistFuncMomentCalc {
   onGrid = phaseGrid,
   confBasis = confBasis,
   phaseBasis = phaseBasis,
   moment = "GkM2",
   gkfacs = {mass, bmag},
}

-- Binary operation updaters.
local weakDivide = Updater.CartFieldBinOp {
   onGrid = confGrid,
   weakBasis = confBasis,
   operation = "Divide",
}
local weakMultiply = Updater.CartFieldBinOp {
   onGrid = confGrid,
   weakBasis = confBasis,
   operation = "Multiply",
}
local weakMultiplyConfPhase = Updater.CartFieldBinOp {
   onGrid = phaseGrid,
   weakBasis = phaseBasis,
   operation = "Multiply",
   fieldBasis = confBasis,
}

-- Project the initial distribution function.
projPhaseScalar:setFunc(initialF)
projPhaseScalar:advance(0.0, {}, {distf})
projConfScalar:setFunc(function(t,xn) return bmagFunc(t,xn) end)
projConfScalar:advance(0.0, {}, {bmag})

-- Compute the initial M0, M1i and M2
calcNumDensity:advance(0.0, {distf}, {m0})
calcMomDensity:advance(0.0, {distf}, {m1})
calcKinEnergyDensity:advance(0.0, {distf}, {m2})

-- Compute u and vt^2.
weakDivide:advance(0., {m0, m1}, {uDrift})
weakMultiply:advance(0., {uDrift, m1}, {m2flow})
vtSq:combine(1., m2, -1., m2flow)
weakDivide:advance(0., {m0, vtSq}, {vtSq})
vtSq:scale(1./vDim)

-- Project the Maxwellian onto the basis.
maxwellian:advance(0.0, {m0, uDrift, vtSq, bmag}, {fM})

-- write out initial data
distf:write("distf.bp")
m0:write("m0.bp")
m1:write("m1.bp")
m2:write("m2.bp")
fM:write("fM-uncorrected.bp")

-- Compute the M0, M1i and M2 from projected Maxwellian
calcNumDensity:advance(0.0, {fM}, {m0_ic})
calcMomDensity:advance(0.0, {fM}, {m1_ic})
calcKinEnergyDensity:advance(0.0, {fM}, {m2_ic})

m0_ic:write("m0-uncorrected.bp")
m1_ic:write("m1-uncorrected.bp")
m2_ic:write("m2-uncorrected.bp")

-- Create updater to fix Maxwellian
local iterFix = Updater.IterGkMaxwellianFix {
   onGrid = phaseGrid,
   phaseBasis = phaseBasis,
   
   confGrid = confGrid,
   confBasis = confBasis,
  
   mass = mass,
   bmag = bmag,

   maxIter = 50,
   relEps = 1e-14,
}

iterFix:advance(0.0, { fM, m0, m1, m2, bmag }, { fOut } )

-- Compute the final M0, M1i and M2
calcNumDensity:advance(0.0, {fOut}, {m0})
calcMomDensity:advance(0.0, {fOut}, {m1})
calcKinEnergyDensity:advance(0.0, {fOut}, {m2})

fOut:write("fM-corrected.bp")
m0:write("m0-corrected.bp")
m1:write("m1-corrected.bp")
m2:write("m2-corrected.bp")
