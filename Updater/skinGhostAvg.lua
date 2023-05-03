local Range       = require "Lib.Range"
local Proto       = require "Lib.Proto"
local skinGhostAvgDecl = require "Updater.skinGhostAvgData.skinGhostAvgDecl"



local skinGhostAvg= Proto(UpdaterBase)


function skinGhostAvg:init(tbl)
   self.grid = tbl.grid
   self.basis = tbl.basis
   self.edge = tbl.edge
   self.bcDir = tbl.bcDir
   self.lowerGhost = tbl.lowerGhost
   self.upperGhost = tbl.upperGhost

   --Define Necessary Ranges
   self.globalRange = self.grid:globalRange()
   if self.edge=='lower' then
      self.globalGhostInRange = self.grid:globalRange():lowerGhost(self.bcDir,self.lowerGhost)
      self.globalSkinInRange = self.grid:globalRange():lowerSkin(self.bcDir, self.lowerGhost)
      self.skinGhostAvgFunc = skinGhostAvgDecl.selectBcModifier('lower', self.basis:id(),self.grid:ndim(), self.basis:polyOrder())
      self.ghostModifier = -self.lowerGhost
   end

   if self.edge=='upper' then
      self.globalGhostInRange = self.grid:globalRange():upperGhost(self.bcDir,self.upperGhost)
      self.globalSkinInRange = self.grid:globalRange():upperSkin(self.bcDir, self.upperGhost)
      self.skinGhostAvgFunc = skinGhostAvgDecl.selectBcModifier('upper', self.basis:id(),self.grid:ndim(), self.basis:polyOrder())
      self.ghostModifier = self.upperGhost
   end
end

function skinGhostAvg:advance(inFld)
   local localGhostInRange = inFld:localRange():intersect(self.globalGhostInRange)
   local localSkinInRange = inFld:localRange():intersect(self.globalSkinInRange)
   lv,uv = localSkinInRange:lowerAsVec(), localSkinInRange:upperAsVec()

   -- Note about the way this loop is done (this comment can be deleted later after review)
   -- I want to access the ghost cell and the skin cell from the same rank so that I can fill both pointers,
   -- then pass them to the kernel. The way I do this is by creating idxGhost which is incremented/decremented by one
   -- depending on whether it is an upper or lower boundary

   --Maybe there is a better way to do this loop, but I have not seen it yet.

   --Also, I think this loop will only work in serial as is
   -- The rank containing the skin cell will not have access to the ghost cell.
   local indexer=inFld:genIndexer()
   local fitrSkin = inFld:get(1)
   local fitrGhost = inFld:get(1)

   for idx in localSkinInRange:rowMajorIter() do
      local idxGhost = idx:copy()
      idxGhost[self.bcDir] = idx[self.bcDir] + self.ghostModifier
      inFld:fill(indexer(idx), fitrSkin)
      inFld:fill(indexer(idxGhost),fitrGhost)
      self.skinGhostAvgFunc(fitrSkin:data(),fitrGhost:data())
   end

end

return skinGhostAvg
