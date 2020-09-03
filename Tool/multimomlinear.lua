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

-- add flags and options
parser:flag("-e --example", "Fetch example input file", false)
parser:option("-i --input", "Input file")

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS)

if args.example then
   -- write out example to terminal
   io.write([[
-- Input file for multimomlinear Tool.

local Species = require "Tool.LinearSpecies"

-- Electrons
elc = Species.Euler {
   mass = 1.0, -- mass
   charge = -1.0, -- charge
   density = 1.0, -- number density
   velocity = {1.0, 0.0, 0.0}, -- velocity vector
   pressure = 0.1, -- pressure
}

-- Ions
ion = Species.Euler {
   mass = 25.0, -- mass
   charge = -.0, -- charge
   density = 1.0, -- number density
   velocity = {0.0, 0.0, 0.0}, -- velocity vector
   pressure = 0.1, -- pressure
}

-- EM field
field = Species.Maxwell {
   epsilon0 = 1.0, mu0 = 1.0,

   electricField = {0.0, 0.0, 0.0}, -- background electric field
   magneticField = {0.0, 0.0, 1.0}, -- background magnetic field
}

-- list of species to include in dispersion relation
speciesList = { elc, ion }

-- Wave vector
kvector = {1.0, 0.0, 0.0}
]]
   )
   
   return
end



-- FFI call prototypes 
ffi.cdef [[
  typedef struct EigenEigen EigenEigen;

  EigenEigen* new_EigenEigen(int N, double *mre, double *mim);
  void delete_EigenEigen(EigenEigen *ee);
  void compute_EigenEigen(EigenEigen* ee, double *ere, double *eim);
]]

-- clear out contents of matrix
local function matrixClear(m, val)
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do
	 m[i][j] = val
      end
   end
end

-- Increment src matrix in sub-region by dest matrix starting from
-- row,col location
local function matrixSubIncr(dest, src, row, col)
   for r = 1, src:numRows() do
      for c = 1, src:numCols() do
	 dest[r+row-1][c+col-1] = dest[r+row-1][c+col-1] + src[r][c]
      end
   end
end

-- Creates eigenvalue solver from given complex matrix
local function createEigenSolver(D)
   local N = D:numRows() -- matrix is square
   local Dre, Dim = Lin.Mat(N,N), Lin.Mat(N,N)

   for i = 1, D:numRows() do
      for j = 1, D:numCols() do
	 Dre[i][j], Dim[i][j] = D[i][j].re, D[i][j].im
      end
   end

   return ffi.C.new_EigenEigen(D:numRows(), Dre:data(), Dim:data())
end

-- Computes eigensystem from created solver
local function computeEigenSystem(N, solver)
   local ere, eim = Lin.Vec(N), Lin.Vec(N)

   ffi.C.compute_EigenEigen(solver, ere:data(), eim:data())

   return ere, eim
end

local count = 1 -- just a counter for number of wave-vectors

-- Solves the dispersion relation in the EM case. Inputs are the wave
-- vector, the list of species objects and the field
-- object. Eigenvalues corresponding to kvec are appended to the
-- eigValues DynVector.
local function solveDispEM(kvec, speciesList, field, eigValues)
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
      -- species' contribution to moment equations
      local Dmom = s:calcMomDispMat(kvec, field.electricFld, field.magneticFld)
      matrixSubIncr(dispMat, Dmom, speciesIdx[sidx], speciesIdx[sidx])

      -- species' contribution to field equations
      local Dfld = s:calcFieldDispMat(kvec, field.electricFld, field.magneticFld)
      matrixSubIncr(dispMat, Dfld, speciesIdx[sidx], fieldIdx)
   end

   -- get field contribution to field equation
   local Dfld = field:calcFieldDispMat(kvec)
   matrixSubIncr(dispMat, Dfld, fieldIdx, fieldIdx)

   -- for each species compute contribution to moment equations (this is
   -- from current sources.)
   for sidx, s in ipairs(speciesList) do
      local Dmom = field:calcMomDispMat(kvec, s.charge, s.density, s.velocity)
      matrixSubIncr(dispMat, Dmom, fieldIdx, speciesIdx[sidx])
   end

   -- create eigenvalue solver
   local eigSolver = createEigenSolver(dispMat)
   -- compute eigensystem
   local evr, evi = computeEigenSystem(dispMat:numRows(), eigSolver)

   -- append to output field
   local kw = { kvec[1], kvec[2], kvec[3] }
   for i = 1, #evr do
      kw[3+2*i-1], kw[3+2*i] = evr[i], evi[i]
   end
   eigValues:appendData(count, kw )
   count = count+1
end

-- open input file and read contents: this loads 'speciesList' and
-- 'field' into this module's namespace
if args.input then
   local inpFile = assert(loadfile(args.input))
   inpFile() -- load contents of input file
else
   print("Must specify an input file to run!")
   return
end

assert(speciesList, "Input file must have a 'speciesList' table!")
assert(field, "Input file must have a 'field'!")

-- modify the output prefix to ensure files are written to the proper
-- place
GKYL_OUT_PREFIX = lfs.currentdir() .. "/" .. args.input:sub(1, args.input:len()-4)

-- total number of equations 
local numEqn = field:numEquations()
for _, s in ipairs(speciesList) do
   numEqn = numEqn + s:numEquations()
end
-- allocate dynvector for output. Arranged as:
-- kx, ky, kz, w1r, w1i, w2r, w2i, ...
local eigValues = DataStruct.DynVector { numComponents = 3 + 2*numEqn }

-- solve dispersion relation
solveDispEM(kvector, speciesList, field, eigValues)

-- write out eigenvaules
eigValues:write("frequencies.bp", 0.0, 0)



