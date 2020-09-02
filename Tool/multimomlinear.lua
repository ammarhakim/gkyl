-- Gkyl ------------------------------------------------------------------------
--
-- Linear dispersion solver for multi-moment multifluid equations.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct = require "DataStruct"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local argparse = require "Lib.argparse"
local complex = require "sci.complex"
local ffi = require "ffi"
local lfs = require "Lib.lfs"

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("multimomlinear")
   :description [[
Solves the linear dispersion relation for multi-moment multifluid
equations. Arbitrary number of species are supported. Each plasma
species can be either an isothermal fluid (including cold fluid),
five-moment or a ten-moment fluid. Some support for Hammett-Perkin
type collisionless closures is provided. The field equations can
either be full Maxwell equations or, for electrostatic problems,
Poisson equations.

For full details on how the dispersion solver works see Developer
Documentation on Gkeyll RTD website.
]]

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS)

-- clear out contents of matrix
local function matrixClear(m, val)
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do
	 m[i][j] = val
      end
   end
end


-- Base species type
local Species = Proto()

function Species:init(tbl)
   self.mass, self.charge = tbl.mass, tbl.charge
end

-- Isothermal species type
local IsoThermal = Proto(Species)

function IsoThermal:init(tbl)
   IsoThermal.super.init(self, tbl)
   
   self.density = tbl.density -- number density
   
   assert(#tbl.velocity == 3, "Must provide all 3 components of velocity")
   self.velocity = tbl.velocity -- velocity vector
   
   self.temperature = tbl.temperature -- species temperature (0.0 is okay)
end

-- Number of equations in system
function IsoThermal:numEquations()
   return 4
end

-- Euler (5M) species type
local Euler = Proto(Species)

function Euler:init(tbl)
   Euler.super.init(self, tbl)

   self.density = tbl.density -- number density
   
   assert(#tbl.velocity == 3, "Must provide all 3 components of velocity")
   self.velocity = tbl.velocity -- velocity vector
   self.pressure = tbl.pressure  -- pressure
   
   self.gasGamma = tbl.gasGamma and tbl.gasGamma or 5.0/3.0
end

-- Number of equations in system
function Euler:numEquations()
   return 5
end

-- Construct matrix with contributions to moment-equation part of
-- dispersion matrix: this is a numEquations X numEquations sized
-- block
function Euler:calcMomDispMat(k, E, B)
   local D = Lin.ComplexMat(5, 5)
   matrixClear(D, 0.0)

   local gamma = self.gasGamma
   local n, p = self.density, self.pressure
   local u = self.velocity
   local m, qbym = self.mass, self.charge/self.mass

   D[1][1] = k[3]*u[3]+k[2]*u[2]+k[1]*u[1] + 0*1i 
   D[1][2] = k[1]*n + 0*1i 
   D[1][3] = k[2]*n + 0*1i 
   D[1][4] = k[3]*n + 0*1i 
   D[1][5] = 0 + 0*1i 
   D[2][1] = 0 + 0*1i 
   D[2][2] = k[3]*u[3]+k[2]*u[2]+k[1]*u[1] + 0*1i 
   D[2][3] = 0 + B[3]*qbym*1i 
   D[2][4] = 0 + -B[2]*qbym*1i 
   D[2][5] = k[1]/(m*n) + 0*1i 
   D[3][1] = 0 + 0*1i 
   D[3][2] = 0 + -B[3]*qbym*1i 
   D[3][3] = k[3]*u[3]+k[2]*u[2]+k[1]*u[1] + 0*1i 
   D[3][4] = 0 + B[1]*qbym*1i 
   D[3][5] = k[2]/(m*n) + 0*1i 
   D[4][1] = 0 + 0*1i 
   D[4][2] = 0 + B[2]*qbym*1i 
   D[4][3] = 0 + -B[1]*qbym*1i 
   D[4][4] = k[3]*u[3]+k[2]*u[2]+k[1]*u[1] + 0*1i 
   D[4][5] = k[3]/(m*n) + 0*1i 
   D[5][1] = 0 + 0*1i 
   D[5][2] = k[1]*p*gamma + 0*1i 
   D[5][3] = k[2]*p*gamma + 0*1i 
   D[5][4] = k[3]*p*gamma + 0*1i 
   D[5][5] = k[3]*u[3]+k[2]*u[2]+k[1]*u[1] + 0*1i 
   
   return D
end

-- Construct matrix with contributions to field part of dispersion
-- matrix: this is a numEquations X 6 sized block (there are 6 maxwell
-- equations)
function Euler:calcFieldDispMat(k, E, B)
   local D = Lin.ComplexMat(5, 6)
   matrixClear(D, 0.0)

   local u = self.velocity
   local qbym = self.charge/self.mass

   D[2][1] = qbym*1i 
   D[2][5] = -u[3]*qbym*1i 
   D[2][6] = u[2]*qbym*1i 
   D[3][2] = qbym*1i 
   D[3][4] = u[3]*qbym*1i 
   D[3][6] = -u[1]*qbym*1i 
   D[4][3] = qbym*1i 
   D[4][4] = -u[2]*qbym*1i 
   D[4][5] = u[1]*qbym*1i
 
   return D
end

-- Ten-moment species type
local TenMoment = Proto(Species)

function TenMoment:init(tbl)
   TenMoment.super.init(self, tbl)

   self.density = tbl.density -- number density
   
   assert(#tbl.velocity == 3, "Must provide all 3 components of velocity")
   self.velocity = tbl.velocity -- velocity vector

   assert(#tbl.pressure == 6, "All six components of pressure (Pxx, Pxy, Pxz, Pyy, Pyz, Pzz) must be provided")
   self.pressureTensor = tbl.pressueTensor -- pressure
end

-- Number of equations in system
function TenMoment:numEquations()
   return 10
end

-- Field equations
local Maxwell = Proto()

function Maxwell:init(tbl)
   self.epsilon0, self.mu0 = tbl.epsilon0, tbl.mu0

   assert(#tbl.electricField == 3, "Must all 3 components of equilibrium electric field")
   self.electricFld = tbl.electricField

   assert(#tbl.magneticField == 3, "Must all 3 components of equilibrium magnetic field")
   self.magneticFld = tbl.magneticField
end

-- Number of equations in system
function Maxwell:numEquations()
   return 6
end

-- Calculate contribution to field part of dispersion matrix
function Maxwell:calcFieldDispMat(k)
   local D = Lin.ComplexMat(6, 6)
   matrixClear(D, 0.0)
   
   local c = self.epsilon0*self.mu0

   D[1][5] = k[3]*c^2
   D[1][6] = -k[2]*c^2
   D[2][4] = -k[3]*c^2 
   D[2][6] = k[1]*c^2 
   D[3][4] = k[2]*c^2 
   D[3][5] = -k[1]*c^2
   D[4][2] = -k[3]
   D[4][3] = k[2]
   D[5][1] = k[3]
   D[5][3] = -k[1] 
   D[6][1] = -k[2]
   D[6][2] = k[1]
end

-- Calculate contribution to moment part of dispersion matrix
function Maxwell:calcMomDispMat(k, q, n, u)
   -- 3 components of current, contributing to n and u[1..3]
   local D = Lin.ComplexMat(3, 4)
   matrixClear(D, 0.0)
   
   local epsilon0 = self.epsilon0

   D[1][1] = -(u[1]*q)/epsilon0*1i 
   D[1][2] = -(n*q)/epsilon0*1i 
   D[1][3] = 0*1i 
   D[1][4] = 0*1i 
   D[2][1] = -(u[2]*q)/epsilon0*1i 
   D[2][2] = 0*1i 
   D[2][3] = -(n*q)/epsilon0*1i 
   D[2][4] = 0*1i 
   D[3][1] = -(u[3]*q)/epsilon0*1i 
   D[3][2] = 0*1i 
   D[3][3] = 0*1i 
   D[3][4] = -(n*q)/epsilon0*1i    
end

-- Just example: will replace with user-specified stuff
local elc = Euler {
   mass = 1.0, 
   charge = -1.0,
   density = 1.0,
   velocity = {0.0, 0.0, 0.0},
   pressure = 0.25,
}

local ion = Euler {
   mass = 25.0, 
   charge = 1.0,
   density = 1.0,
   velocity = {0.0, 0.0, 0.0},
   pressure = 0.25,
}

local field = Maxwell {
   epsilon0 = 1.0, mu0 = 1.0,

   electricField = {0.0, 0.0, 0.0},
   magneticField = {0.0, 0.0, 1.0},
}


-- construct list of species
local speciesList = {elc, ion}
local kvec = {0.5, 0.0, 0.0} -- wavevector

-- total number of equations 
local numEqn = field:numEquations()
for _, s in ipairs(speciesList) do
   numEqn = numEqn + s:numEquations()
end

-- Allocate space for dispersion matrix
local dispMat = Lin.ComplexMat(numEqn, numEqn)
matrixClear(dispMat, 0.0)

-- compute starting index into each species sub-matrices in dispMat
local speciesIdx = {}
local currStart = 1
for sidx, s in ipairs(speciesList) do
   speciesIdx[sidx] = currStart
   currStart = currStart+s:numEquations()
end
local fieldIdx = currStart -- starting index into EM field

-- loop over each species, inserting contribution to dispersion matrix
for sidx, s in ipairs(speciesList) do
   -- get this species' contribution to moment equation part of
   -- dispersion matrix
   local Dmom = s:calcMomDispMat(kvec, field.electricFld, field.magneticFld)

   -- accumulate contribution to dispMat
   local i0, j0 = speciesIdx[sidx]-1, speciesIdx[sidx]-1
   for i = 1, s:numEquations() do
      for j = 1, s:numEquations() do
	 dispMat[i+i0][j+j0] = dispMat[i+i0][j+j0] + Dmom[i][j]
      end
   end

   -- get species' contribution to field part of dispersion matrix
   local Dfld = s:calcFieldDispMat(kvec, field.electricFld, field.magneticFld)
   local i0, j0 = speciesIdx[sidx]-1, fieldIdx-1
   for i = 1, s:numEquations() do
      for j = 1, field:numEquations() do
	 dispMat[i+i0][j+j0] = dispMat[i+i0][j+j0] + Dfld[i][j]
      end
   end
   
end
