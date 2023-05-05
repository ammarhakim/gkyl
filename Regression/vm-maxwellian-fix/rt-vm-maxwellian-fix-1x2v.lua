-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to fix maxwellian with iteration method, 1x2v. 
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"
local Time = require "Lib.Time"

-- Reference density, drift velocity and thermal speed.
local n0 = 1.0
local u0 = {0.0, 0.75}
local vt0 = 1.0   

-- Phase space mesh.
local lower = {-0.50, -6.0*vt0, -6.0*vt0}
local upper = { 0.50,  6.0*vt0,  6.0*vt0}
local numCells = {8, 8, 8}

local polyOrder = 1 -- Polynomial order of the basis.
local pDim = #numCells
local vDim = #u0
local cDim = pDim - vDim

local lowerC, upperC, numCellsC = {}, {}, {}
for d=1, cDim do
   lowerC[d] = lower[d]
   upperC[d] = upper[d]
   numCellsC[d] = numCells[d]
end

-- Functions defining the density, drift velocity and thermal speed squared.
local eps_n = 0.1
local eps_u = 0.1
local eps_vt = 0.1
local m0Func = function (t, xn) return n0*(1.0+eps_n*math.cos(2.*math.pi*xn[1])) end
local uDriftFunc = function (t, xn) return u0[1]*(1.0+eps_u*math.cos(2.*math.pi*xn[1])), u0[2]*(1.0+eps_u*math.cos(2.*math.pi*xn[1]))  end
local vtSqFunc = function (t, xn) return (vt0*(1.0+eps_vt*math.cos(2.*math.pi*xn[1])))^2 end

-- Initial distribution function
local initialF = function(t, xn)
   local x, vx, vy = xn[1], xn[2], xn[3]

   local den = m0Func(t,xn)
   local ux, uy = uDriftFunc(t,xn)
   local vtsq = vtSqFunc(t,xn)

   return (den/math.sqrt((2.*math.pi*vtsq)^vDim))*math.exp(-((vx-ux)^2+(vy-uy)^2)/(2.*vtsq))
end

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
local phaseBasis = createPhaseBasis(cDim, vDim, polyOrder)
local confGrid = createGrid(lowerC, upperC, numCellsC)
local confBasis = createConfBasis(confGrid:ndim(), polyOrder)

-- Create various fields needed.
local distf = createField(phaseGrid, phaseBasis)
local fM = createField(phaseGrid, phaseBasis)
local fOut = createField(phaseGrid, phaseBasis)

local m0 = createField(confGrid, confBasis) 
local m1i = createField(confGrid, confBasis, vDim)
local m2 = createField(confGrid, confBasis)
local m2flow = createField(confGrid, confBasis)
local uDrift = createField(confGrid, confBasis, vDim)
local vtSq = createField(confGrid, confBasis)

local m0_ic = createField(confGrid, confBasis) 
local m1i_ic = createField(confGrid, confBasis, vDim)
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
}

-- Moment updaters.
local calcNumDensity = Updater.DistFuncMomentCalc {
   onGrid = phaseGrid,
   confBasis = confBasis,
   phaseBasis = phaseBasis,
   moment = "M0",
}
local calcMomDensity = Updater.DistFuncMomentCalc {
   onGrid = phaseGrid,
   confBasis = confBasis,
   phaseBasis = phaseBasis,
   moment = "M1i",
}
local calcKinEnergyDensity = Updater.DistFuncMomentCalc {
   onGrid = phaseGrid,
   confBasis = confBasis,
   phaseBasis = phaseBasis,
   moment = "M2",
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

-- Compute the initial M0, M1i and M2
calcNumDensity:advance(0.0, {distf}, {m0})
calcMomDensity:advance(0.0, {distf}, {m1i})
calcKinEnergyDensity:advance(0.0, {distf}, {m2})

-- Compute u and vt^2.
weakDivide:advance(0., {m0, m1i}, {uDrift})
weakMultiply:advance(0., {uDrift, m1i}, {m2flow})
vtSq:combine(1., m2, -1., m2flow)
weakDivide:advance(0., {m0, vtSq}, {vtSq})
vtSq:scale(1./vDim)

-- Project the Maxwellian onto the basis.
maxwellian:advance(0.0, {m0, uDrift, vtSq}, {fM})

-- write out initial data
distf:write("distf.bp")
m0:write("m0.bp")
m1i:write("m1i.bp")
m2:write("m2.bp")
fM:write("fM-uncorrected.bp")

-- Compute the M0, M1i and M2 from projected Maxwellian
calcNumDensity:advance(0.0, {fM}, {m0_ic})
calcMomDensity:advance(0.0, {fM}, {m1i_ic})
calcKinEnergyDensity:advance(0.0, {fM}, {m2_ic})

m0_ic:write("m0-uncorrected.bp")
m1i_ic:write("m1i-uncorrected.bp")
m2_ic:write("m2-uncorrected.bp")

-- Create updater to fix Maxwellian
local iterFix = Updater.IterMaxwellianFix {
   onGrid = phaseGrid,
   phaseBasis = phaseBasis,
   
   confGrid = confGrid,
   confBasis = confBasis,

   maxIter = 50,
   relEps = 1e-14,

}

iterFix:advance(0.0, { fM, m0, m1i, m2 }, { fOut } )

-- Compute the final M0, M1i and M2
calcNumDensity:advance(0.0, {fOut}, {m0})
calcMomDensity:advance(0.0, {fOut}, {m1i})
calcKinEnergyDensity:advance(0.0, {fOut}, {m2})

fOut:write("fM-corrected.bp")
m0:write("m0-corrected.bp")
m1i:write("m1i-corrected.bp")
m2:write("m2-corrected.bp")
