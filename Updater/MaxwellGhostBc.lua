-- System libraries
local xsys = require "xsys"

-- Gkyl libraries.
local CartDecomp     = require "Lib.CartDecomp"
local CartFieldBinOp = require "Updater.CartFieldBinOp"
local DataStruct     = require "DataStruct"
local Grid           = require "Grid"
local Lin            = require "Lib.Linalg"
local LinearDecomp   = require "Lib.LinearDecomp"
local Mpi            = require "Comm.Mpi"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local Proto          = require "Lib.Proto"
local Range          = require "Lib.Range"
local UpdaterBase    = require "Updater.Base"
local MomDecl        = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"
local CartFieldBinOp= require "Updater.CartFieldBinOp"
local DistFuncMomentCalc  = require "Updater.DistFuncMomentCalc"

-- Boundary condition updater.
local MaxwellGhostBc = Proto(UpdaterBase)
local dirlabel = {"X", "Y", "Z"}

function MaxwellGhostBc:init(tbl)
   MaxwellGhostBc.super.init(self, tbl) -- Setup base object.
   self._boundaryGrid = tbl.boundaryGrid
   self._boundaryGrid = tbl.confBouundaryGrid
   self._dir = tbl.dir
   self._edge = tbl.edge
   --self.localGhostRangeWOCorners = tbl.localGhostRangeWOCorners
   self.myGLobalGhostRange = tbl.myGLobalGhostRange
   self.ghostFld = tbl.ghostFld
   self._dirlabel = dirlabel[self._dir]
end

-- The advance method
function MaxwellGhostBc:advance(tCurr, inFld, outFld)
   local fIn = inFld
   local fOut = outfld[1]
   fOut:copyRangeToRange(fIn,self.localGhostRangeWOCorners, self.localGhostRangeWOCorners )
end


function MaxwellGhostBc:getDir() return self._dir end

function MaxwellGhostBc:getEdge() return self._edge end

function MaxwellGhostBc:label() return "Flux"..self._dirlabel..self._edge end

function MaxwellGhostBc:getBoundaryGrid() return self._boundaryGrid end

return MaxwellGhostBc
