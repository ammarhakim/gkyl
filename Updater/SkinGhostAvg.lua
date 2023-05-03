-- Gkyl ------------------------------------------------------------------------
--
-- Updater to Average the skin and ghost cells.
--
--
-- The Skin cells in the specified bcDir will be modified so that the boundary
-- vals become the average of the input field's skin and ghost cell boundary vals
-- Intended for 3x.
--
--
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Range       = require "Lib.Range"
local Proto       = require "Lib.Proto"
local skinGhostAvgDecl = require "Updater.skinGhostAvgData.skinGhostAvgDecl"



local SkinGhostAvg= Proto(UpdaterBase)


function SkinGhostAvg:init(tbl)
   self.grid = assert(tbl.grid, "Updater.SkinGhostAvg: Must provide grid object using 'grid'")
   self.basis = assert(tbl.basis, "Updater.SkinGhostAvg: Must provide basis using 'basis'")
   self.edge = tbl.edge
   self.bcDir = assert(tbl.bcDir, "Updater.skinGhostAvg: Must provide direction using 'bcDir'")
   self.lowerGhost = assert(tbl.lowerGhost, "Updater.SkinGhostAvg: Must provide number of lower ghost cells using 'lowerGhost'")
   self.upperGhost = assert(tbl.upperGhost, "Updater.SkinGhostAvg: Must provide number of upper ghost cells using 'upperGhost'")
   local advArgs = assert(tbl.advanceArgs, "Updater.SkinGhostAvg: Must provide a sample field in 'advanceArgs'")
   self.sampleFld = advArgs[1]

   assert(self.edge=='upper' or self.edge=='lower', "Updader.SkinGhostAvg: 'edge' must be 'upper' or 'lower'")
   assert(self.basis:id()=="serendipity", "Updater.SkinGhostAvg: only implemented for modal serendipity basis")
   assert(self.grid:ndim()==3, "Updater.SkinGhostAvg: only implemented for 3D")
   assert(self.grid:ndim() == self.basis:ndim(), "Updater.SkinGhostAvg: dimensions of basis and grid must match")
   assert(self.basis:polyOrder()==1, "Updater.SkinGhostAvg: only implemented for polyOrder = 1")

   --Define Necessary Ranges
   self.globalRange = self.grid:globalRange()
   if self.edge=='lower' then
      self.globalSkinInRange = self.grid:globalRange():lowerSkin(self.bcDir, self.lowerGhost)
      self.skinGhostAvgFunc = skinGhostAvgDecl.selectAvg('lower', self.basis:id(),self.grid:ndim(), self.basis:polyOrder())
      self.ghostModifier = -self.lowerGhost
   end

   if self.edge=='upper' then
      self.globalSkinInRange = self.grid:globalRange():upperSkin(self.bcDir, self.upperGhost)
      self.skinGhostAvgFunc = skinGhostAvgDecl.selectAvg('upper', self.basis:id(),self.grid:ndim(), self.basis:polyOrder())
      self.ghostModifier = self.upperGhost
   end
   self.localSkinInRange = self.sampleFld:localRange():intersect(self.globalSkinInRange)
end

function SkinGhostAvg:advance(inFld)
   local indexer=inFld:genIndexer()
   local fitrSkin = inFld:get(1)
   local fitrGhost = inFld:get(1)

   for idx in self.localSkinInRange:rowMajorIter() do
      local idxGhost = idx:copy()
      idxGhost[self.bcDir] = idx[self.bcDir] + self.ghostModifier
      inFld:fill(indexer(idx), fitrSkin)
      inFld:fill(indexer(idxGhost),fitrGhost)
      self.skinGhostAvgFunc(fitrSkin:data(),fitrGhost:data())
   end

end

return SkinGhostAvg
