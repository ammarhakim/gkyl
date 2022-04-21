-- Gkyl ------------------------------------------------------------------------
-- 
-- Compute the contribution to div{ P } from one species for a 1D field
-- equation (GkNoCurrentEpar) to be used in GK, assuming vanishing currents
-- and quasineutrality.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase  = require "Updater.Base"
local DataStruct   = require "DataStruct"
local LinearDecomp = require "Lib.LinearDecomp"
local Lin          = require "Lib.Linalg"
local Proto        = require "Lib.Proto"
local xsys         = require "xsys"
local ModDecl      = require "Updater.zeroCurrentGkEparDivPData.zeroCurrentGkEpar_divP_mod_decl"

-- Inherit the base Updater from UpdaterBase updater object.
local ZeroJGkEparDivP = Proto(UpdaterBase)

function ZeroJGkEparDivP:init(tbl)
   ZeroJGkEparDivP.super.init(self, tbl)   -- Setup base object.

   -- Get grid and basis.
   self._grid  = assert(tbl.onGrid, "ZeroCurrentGkEparDivP: must specify a grid.")
   self._basis = assert(tbl.basis, "ZeroCurrentGkEparDivP: must specify a basis.")

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self.dim = self._basis:ndim()
   local nm, p = self._basis:id(), self._basis:polyOrder()

   self._divPKer = ModDecl.selectDivPKer(self._basis:id(), self.dim, self._basis:polyOrder())

   self.dx = Lin.Vec(self.dim)

end

function ZeroJGkEparDivP:_advance(tCurr, inFlds, outFlds)

   local charge = assert(inFlds[1], "ZeroJGkEparDivP: Must specify the charge in 'inFld[1]'.")
   local delparFac = assert(inFlds[2], "ZeroJGkEparDivP: Must specify delparFac=cmag/(J*B) in 'inFld[2]'.")
   local dlnbmagdz = assert(inFlds[3], "ZeroJGkEparDivP: Must specify dlnbmagdz in 'inFld[3]'.")
   local m2par = assert(inFlds[4], "ZeroJGkEparDivP: Must specify m2par in 'inFld[4]'.")
   local m2perp = assert(inFlds[5], "ZeroJGkEparDivP: Must specify m2perp in 'inFld[5]'.")
   local divPFld = assert(outFlds[1], "ZeroJGkEparDivP: Must specify an output field 'outFld[1]'.")

   local indexer = divPFld:genIndexer()

   local grid = self._grid
   local divPRangeDecomp
   if self.onGhosts then
      divPRangeDecomp = LinearDecomp.LinearDecompRange {
         range = divPFld:localExtRange(), numSplit = grid:numSharedProcs() }
   else
      divPRangeDecomp = LinearDecomp.LinearDecompRange {
         range = divPFld:localRange(), numSplit = grid:numSharedProcs() }
   end

   local delparFacItr, dlnbmagdzItr = delparFac:get(1), dlnbmagdz:get(1)
   local m2parItr, m2perpItr = m2par:get(1), m2perp:get(1)
   local divPItr = divPFld:get(1)

   local tId = grid:subGridSharedId()   -- Local thread ID.
   for idx in divPRangeDecomp:rowMajorIter(tId) do
      grid:setIndex(idx)
      grid:getDx(self.dx)

      delparFac:fill(indexer(idx), delparFacItr)
      dlnbmagdz:fill(indexer(idx), dlnbmagdzItr)
      m2par:fill(indexer(idx), m2parItr)
      m2perp:fill(indexer(idx), m2perpItr)
      divPFld:fill(indexer(idx), divPItr)

      self._divPKer(self.dx[1], charge, delparFacItr:data(), dlnbmagdzItr:data(), m2parItr:data(), m2perpItr:data(), divPItr:data())
   end

end

return ZeroJGkEparDivP
