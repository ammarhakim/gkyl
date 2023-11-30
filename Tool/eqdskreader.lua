-- Gkyl ------------------------------------------------------------------------
--
-- Reads eqdsk files and converts to BP files
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Alloc = require "Lib.Alloc"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local argparse = require "Lib.argparse"

local Mpi = require "Comm.Mpi"
local Adios = require "Io.Adios"

GKYL_ADIOS2_MPI = GKYL_ADIOS2_MPI or Adios.init_mpi(Mpi.COMM_WORLD)

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("eqdskreader")
   :description [[
Reads eqdsk files and converts the data into ADIOS-BP so they can be
read in gkyel and postgkyl. The eqdsk file format is somewhat strange
and it is not clear if this tool can read every eqdsk file out
there. Some caution is needed.
]]

-- add flags and options
parser:option("-i --input", "Input eqdsk file to read")

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS)

-- modify prefix to write files based on eqdsk file name
GKYL_OUT_PREFIX = lfs.currentdir() .. "/" .. args.input

-- Reads and returns header information
local function readHeader(fh)
   local ln1 = fh:read()
   local n = ln1:len()
   
   local s,e = ln1:find("%d+")
   local nx = tonumber(ln1:sub(s,e))

   ln1 = ln1:sub(e+1, n)
   local s,e = ln1:find("%d+")
   local ny = tonumber(ln1:sub(s,e))

   return nx, ny
end

-- reads eqdsk file all into memory at once. Further splitting into
-- pieces is not done in this function.
--
-- fh : Handle to open file
-- returns nx, ny, data
local function readEqdsk(fh)
   local nx, ny = readHeader(fh)
   local data = Alloc.Double()

   numMatch = "[+-]?%d+%.%d+[eE][+-]%d+"
   for ln in fh:lines() do
      local n = ln:len()
      local ln1 = ln:sub(1, n) -- copy so we can mod
   
      while true do
	 local s,e = ln1:find(numMatch)
	 if s ~= nil and e ~= nil then
	    data:push( tonumber(ln1:sub(s,e)) )
	    ln1 = ln1:sub(e+1, n)
	 else
	    break
	 end
      end
   end

   return nx, ny, data
end

-- open file and read data
fh = io.open(args.input, "r")
local nx, ny, data = readEqdsk(fh)
fh:close()

-- function to allocate fields
local function makeField(meta, grid)
   return DataStruct.Field {
      onGrid = grid,
      numComponents = 1,
      metaData = meta,
   }
end


-- extract various pieces of from read data
local rdim, zdim = data[1], data[2]
local rcenter, rleft = data[3], data[4]
local zmid = data[5]

local rmaxis, zmaxis = data[6], data[7]
local simag, sibry = data[8], data[9]
local bcenter = data[10]

local current = data[11]

-- store scalars in a table to put into meta-data file
local metaData = {
   rdim = rdim, zdim = zdim,
   rcenter = rcenter, rleft = rleft,
   zmid = zmid,

   rmaxis = rmaxis, zmaxis = zmaxis,
   simag = simag, sibry = sibry,
   bcenter = bcenter,

   current = current,
}

-- create grids
local grid1d = Grid.RectCart {
   lower = { simag },
   upper = { sibry },
   cells = { nx }
}
local grid2d = Grid.RectCart {
   lower = { rleft, zmid-zdim/2},
   upper = { rleft+rdim, zmid+zdim/2 },
   cells = { nx, ny }
}

-- allocate various fields
fpol = makeField(meta, grid1d)
pres = makeField(meta, grid1d)
ffprime = makeField(meta, grid1d)
pprime = makeField(meta, grid1d)
psi = makeField(meta, grid2d)
qpsi = makeField(meta, grid1d)

-- construct 1D arrays
for i = 1, nx do
   fpol:get(i)[1] = data[i+20]
   pres:get(i)[1] = data[i+20 + nx]
   ffprime:get(i)[1] = data[i+20 + 2*nx]
   pprime:get(i)[1] = data[i+20 + 3*nx]
   qpsi:get(i)[1] = data[i+20 + 4*nx + nx*ny]
end

local indexer = psi:indexer()
for i = 1, nx do
   for j = 1, ny do
      psi:get(indexer(i,j))[1] = data[(i-1)+20 + (j-1)*nx + 4*nx + 1]
   end
end

-- write data to BP file
fpol:write("fpol.bp")
pres:write("pres.bp")
ffprime:write("ffprime.bp")
pprime:write("pprime.bp")
psi:write("psi.bp")
qpsi:write("qpsi.bp")

