-- Gkyl ------------------------------------------------------------------------
--
-- Finds all the restart files for simulation with given input file, detects
-- which of those contain Cartesian grid fields, and interpolates
-- those fields onto the grid with resolution given by the user.
-- Common use: to restart a simulation using a higher resolution.
--
-- Use:
--   gkyl refinerestart -i input_file.lua -r output_resolution -o output_directory
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local lfs         = require "lfs"
local AdiosReader = require "Io.AdiosReader"
local Logger      = require "Lib.Logger"
local Grid        = require "Grid"
local Basis       = require "Basis"
local DataStruct  = require "DataStruct"
local Mpi         = require "Comm.Mpi"
local ffi         = require "ffi"
local Adios       = require "Io.Adios"
local Updater     = require "Updater"
local argparse    = require "Lib.argparse"
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"

local log = Logger { logToFile = true }
local verboseLog = function (msg) end -- default no messages are written
local verboseLogger = function (msg) log(msg) end

local function replaceHyphen(str) return str:gsub("%-", "---") end

local getFilePath = function(file)
   -- Extract directory from full filepath-name 'file' (may or may not be in it).
   local fileRev  = string.reverse(file)
   local rsS, rsE = string.find(fileRev,"/")
   local filePath
   if rsS and rsE then
      filePath = string.reverse(string.sub(fileRev, rsS, string.len(file)))
   else
      filePath = lfs.currentdir()
   end
   return filePath
end

local getFileName = function(file)
   -- Extract file name from full filepath-name 'file'.
   local fileRev  = string.reverse(file)
   local rsS, rsE = string.find(fileRev,"/")
   local fileName
   if rsS and rsE then
      fileName = string.reverse(string.sub(fileRev, 1, rsS-1))
   else
      fileName = file
   end
   return fileName
end

local findRestart = function(inputFile)
   -- Given the input file for a simulation (inputFile)
   -- create a list of all restart files with a CartGridField in it.
   local restartFiles = {}

   local path = getFilePath(inputFile)
   
   -- Loop over files and make a list of restart files.
   for file in lfs.dir(path) do
      if file ~= "." and file ~= ".." and string.sub(file,1,2) ~= ".!" then
         local fullNm = path .. "/" .. file
         local sS, sE = string.find(fullNm, '_restart.bp')
         if sS and sE then 
            -- Don't include files that are not defined on a mesh.
            local r1 = AdiosReader.Reader(fullNm)
            if r1:hasAttr("lowerBounds") and r1:hasAttr("upperBounds") and
               r1:hasAttr("numCells") and r1:hasAttr("polyOrder") then
               table.insert(restartFiles, fullNm)
            end
         end
      end
   end

   return restartFiles
end

local function createGrid(lo,up,nCells)
   local gridOut = Grid.RectCart {
      lower = lo,
      upper = up,
      cells = nCells,
   }
   return gridOut
end

local function createBasis(dim, pOrder, bKind)
   local basis
   if bKind=="Ser" or bKind=="serendipity" then
      basis = Basis.CartModalSerendipity { ndim = dim, polyOrder = pOrder }
   elseif bKind=="Max" or bKind=="maximal-order" then
      basis = Basis.CartModalMaxOrder { ndim = dim, polyOrder = pOrder }
   elseif bKind=="Tensor" or bKind=="tensor" then
      basis = Basis.CartModalTensor { ndim = dim, polyOrder = pOrder }
   else
      assert(false,"Invalid basis")
   end
   return basis
end

local function createField(grid, basis, comp)
   comp = comp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = comp,
      ghost         = {0, 0},
   }
   return fld
end

local interpFile = function(file, simName, nCellsOut, outpath)
   -- Interpolate the CartGridField(s) in 'file' to a grid with
   -- 'nCellsOut' cells and write the new file to 'outpath'.
   -- 'file' is produced with a simulation whose input file
   -- without .lua is 'simName'

   -- Create donor and target grids:
   local frh = AdiosReader.Reader(file)
   local numCells_do    = frh:getAttr('numCells')._values
   local lowerBounds_do = frh:getAttr('lowerBounds')._values
   local upperBounds_do = frh:getAttr('upperBounds')._values
   local polyOrder_do   = frh:getAttr('polyOrder')._values[1]
   local basisType_do   = frh:getAttr('basisType')._values[1]
   frh:close()

   local numCells_tar = {}
   for d =1, #numCells_do do numCells_tar[d] = nCellsOut[d] end

   local grid = { _do  = createGrid(lowerBounds_do, upperBounds_do, numCells_do),
                  _tar = createGrid(lowerBounds_do, upperBounds_do, numCells_tar) }
                 
   local basis = createBasis(grid["_do"]:ndim(), polyOrder_do, basisType_do)

   -- Updater to interpolate a CartField from one grid to another.
   local interpUpd = Updater.CartFieldInterpolate {
      onGrid   = grid["_tar"],  onBasis   = basis,
      fromGrid = grid["_do"],   fromBasis = basis,
   }

   -- Make a list of fields to interpolate and obtain the shape of each
   -- (mainly to know the number of components) without reading the field.
   local fd = Adios.read_open_file(file, Mpi.COMM_WORLD)
   assert(not (fd == nil), string.format("Cannot open ADIOS file %s for reading", file))
   local varShapes, time, frame = {}, nil, nil
   for i = 0, fd.nvars-1 do
      local nm    = ffi.string(fd.var_namelist[i])
      local obj   = Adios.inq_var_byid(fd, i)
      if obj.ndim > 0 then
         varShapes[nm] = {}
         for j = 0, obj.ndim-1 do varShapes[nm][j+1] = tonumber(obj.dims[j]) end
      else
         if nm == "time"  then time  = ffi.cast("double *",obj.value)[0] end
         if nm == "frame" then frame = ffi.cast("int *",obj.value)[0] end
      end
   end
   Adios.read_close(fd)

   -- Create data fields.
   local fld_do, fld_tar = {}, {}
   local elemType
   for varNm, varSh in pairs(varShapes) do
      fld_do[varNm]  = createField(grid["_do"] , basis, varSh[#varSh]) 
      fld_tar[varNm] = createField(grid["_tar"], basis, varSh[#varSh]) 
      elemType = elemType and elemType or fld_do[varNm]:elemType()
   end

   -- Read in donor field.
   local fieldIo = { _do  = AdiosCartFieldIo{elemType   = elemType,  method = "MPI",
                                             writeGhost = true,
                                             metaData   = { polyOrder = basis:polyOrder(),
                                                            basisType = basis:id()},},
                     _tar = AdiosCartFieldIo{elemType   = elemType,  method = "MPI",
                                             writeGhost = true,
                                             metaData   = { polyOrder = basis:polyOrder(),
                                                            basisType = basis:id()},}, }
   -- Change GKYL_OUT_PREFIX and craft outSuffix so that they
   -- work with the way they are used in AdiosCartFieldIo.
   GKYL_OUT_PREFIX = getFilePath(file) .. simName
   local fileSuffix, _ = string.gsub(replaceHyphen(file), replaceHyphen(GKYL_OUT_PREFIX).."_", "")
   fieldIo["_do"]:read(fld_do, fileSuffix, true)

   -- Interpolate.
   for varNm, varSh in pairs(varShapes) do
      interpUpd:advance(0.0,{fld_do[varNm]},{fld_tar[varNm]})
   end

   -- Write target field to file:
   -- Change GKYL_OUT_PREFIX and craft outSuffix so that they
   -- work with the way they are used in AdiosCartFieldIo.
   -- Also: have to change the name of the file (added _refined) because otherwise
   --       AdiosCartFieldIo uses the same group and therefore the donor grid.
   GKYL_OUT_PREFIX = outpath .. simName
   fieldIo["_tar"]:write(fld_tar, string.sub(fileSuffix,1,-4).."_refined.bp", time, frame, true)
     
end

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("help")
   :description [[Interpolate restart files onto another mesh.]]
parser:option("-i --input", "Input file name, e.g. vlasov-shock.lua")
parser:option("-r --resolution", "Resolution of grid to interpolate to, e.g. 32,64")
parser:option("-o --outdir", "Output directory, e.g. ./newres/")

-- Parse command line parameters.
local args = parser:parse(GKYL_COMMANDS)

if (not args.input) or (not args.resolution) or (not args.outdir) then
   print("Must specify simulation input file, output resolution and output directory with -i, -r and -o")
else
   inputFile, outres, outdir = args.input, args.resolution, args.outdir
   local resFiles = findRestart(inputFile)
   numCells_tar = {}
   for n in string.gmatch( outres, "(%d+)" ) do table.insert(numCells_tar, tonumber(n)) end
   verboseLogger("  Will output new restart files with "..outres.." cells to "..outdir.."\n")

   -- Get input file name w/o .lua.
   local simName = string.sub(getFileName(inputFile),1,-5)

   -- Interpolate each field to new mesh and writ it out.
   for i, f in ipairs(resFiles) do interpFile(f, simName, numCells_tar, outdir) end
end
