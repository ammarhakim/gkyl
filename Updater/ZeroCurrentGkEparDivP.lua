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

   self.dim = self.basis:ndim()
   local nm, p = self._basis:id(), self._basis:polyOrder()

   self._divPKer = ModDecl.selectDivPKer(self.basis:id(), self.dim, self.basis:polyOrder())

   self.dx = Lin.Vec(self.dim)

end

function ZeroJGkEparDivP:_advance(tCurr, inFlds, outFlds)

   local charge = assert(inFld[1], "ZeroJGkEparDivP: Must specify the charge in 'inFld[1]'.")
   local cmag = assert(inFld[2], "ZeroJGkEparDivP: Must specify cmag in 'inFld[2]'.")
   local jacobTotInv = assert(inFld[3], "ZeroJGkEparDivP: Must specify jacobTotInv in 'inFld[3]'.")
   local jacobGeoInv = assert(inFld[4], "ZeroJGkEparDivP: Must specify jacobGeoInv in 'inFld[4]'.")
   local bmagInv = assert(inFld[5], "ZeroJGkEparDivP: Must specify bmagInv in 'inFld[5]'.")
   local dbmagdz = assert(inFld[6], "ZeroJGkEparDivP: Must specify dbmagdz in 'inFld[6]'.")
   local m2par = assert(inFld[7], "ZeroJGkEparDivP: Must specify m2par in 'inFld[7]'.")
   local m2perp = assert(inFld[8], "ZeroJGkEparDivP: Must specify m2perp in 'inFld[8]'.")
   local divPFld = assert(outFld[1], "ZeroJGkEparDivP: Must specify an output field 'outFld[1]'.")

   local grid = self.grid
   if self.onGhosts then
      divPRangeDecomp = LinearDecomp.LinearDecompRange {
         range = divPFld:localExtRange(), numSplit = grid:numSharedProcs() }
   else
      divPRangeDecomp = LinearDecomp.LinearDecompRange {
         range = divPFld:localRange(), numSplit = grid:numSharedProcs() }
   end

   local cmagItr, jacobTotInvItr, jacobGeoInvItr = cmag:get(1), jacobTotInv:get(1), jacobGeoInv:get(1)
   local bmagInvItr, dbmagdzItr = bmagInv:get(1), dbmagdz:get(1)
   local m2parItr, m2perpItr = m2par:get(1), m2perp:get(1)
   local divPFldItr = divPFld:get(1)

   local tId = grid:subGridSharedId()   -- Local thread ID.
   for idx in divPRangeDecomp:rowMajorIter(tId) do
      grid:setIndex(idx)
      grid:getDx(self.dx)

      cmag:fill(indexer(idx), cmagItr)
      jacobTotInv:fill(indexer(idx), jacobTotInvItr)
      jacobGeoInv:fill(indexer(idx), jacobGeoInvItr)
      bmagInv:fill(indexer(idx), bmagInvItr)
      dbmagdz:fill(indexer(idx), dbmagdzItr)
      m2par:fill(indexer(idx), m2parItr)
      m2perp:fill(indexer(idx), m2perpItr)
      divPFld:fill(indexer(idx), divPItr)

      self._divPKer(self.dx[1], charge, cmagItr:data(), jacobTotInvItr:data(), jacobGeoInvItr:data(), bmagInvItr:data(), dbmagdzItr:data(), m2parItr:data(), m2perpItr:data(), divPItr:data())
   end

end

return ZeroJGkEparDivP
