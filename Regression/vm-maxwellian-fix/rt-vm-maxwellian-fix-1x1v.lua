-- Gkyl ------------------------------------------------------------------------
--
--
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"
local Time = require "Lib.Time"

-- Reference density, drift velocity and thermal speed.
local n0 = 1.0
local u0 = 1.0
local vt0 = 1.0   

-- Functions defining the density, drift velocity and thermal speed squared.
local m0Func     = function (t, xn) return n0 end
local uDriftFunc = function (t, xn) return u0 end
local vtSqFunc   = function (t, xn) return vt0^2 end

-- Initial distribution function
local initialF = function(t, xn)
   local x, vx = xn[1], xn[2]

   local den = m0Func(t,xn)
   local ux = uDriftFunc(t,xn)
   local vtsq = vtSqFunc(t,xn)

   return (den/math.sqrt(2.*math.pi*vtsq))*math.exp(-((vx-ux)^2)/(2.*vtsq))
end

-- Phase space mesh.
local lower = {-0.50, -6.0*vt0}
local upper = { 0.50,  6.0*vt0}
local numCells = {8, 8}

local polyOrder = 1 -- Polynomial order of the basis.
local cDim = 1
local vDim = 1
local pDim = cDim + vDim

local lowerC, upperC, numCellsC = {}, {}, {}
for d=1, cDim do
   lowerC[d] = lower[d]
   upperC[d] = upper[d]
   numCellsC[d] = numCells[d]
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
local m0 = createField(confGrid, confBasis) 
local m1i = createField(confGrid, confBasis, vDim)
local m2 = createField(confGrid, confBasis)
local uDrift = createField(confGrid, confBasis, vDim)
local vtSq = createField(confGrid, confBasis)
