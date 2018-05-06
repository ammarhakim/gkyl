-- Gkyl ------------------------------------------------------------------------
--
-- Tests for weak multiplication and division binary operations.
-- Here we use the calculation of primitive moments as a test case.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Basis = require "Basis"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"
local CalcDiagnostic = require "Updater.CalcDiagnostic"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats = Unit.stats

------------------------------------------------------------------------------------
---                    Test in 1x
------------------------------------------------------------------------------------
function test_binOp1x(nx, p, writeMatrix)
   writeMatrix = writeMatrix or false
   -- Phase-space and config-space grids.
   local confGrid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {nx},
   }
   -- Basis functions.
   local basis = Basis.CartModalSerendipity { ndim = 1, polyOrder = p }
   io.write("nx=",nx," polyOrder=", p, "\n")

   -- Zeroth and first (parallel) moments.
   local numDens = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1},
   }
   local mom1A = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1},
   }
   local mom1 = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1},
   }
   -- Calculated and analytica parallel flow speed.
   local Ui = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1},
   }
   local UiA = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1},
   }

   -- Initialize number density and first (parallel) moment.
   local initDens = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = basis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x = xn[1]
--                    return math.exp(-((x-mu)/(math.sqrt(2)*sig))^2)
--                    return 0.90-x
--                    return math.exp(-x)-0.5
                    return math.exp(-x)
                 end
   }
   initDens:advance(0.,0.,{},{numDens})
   local initMom1i = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = basis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x = xn[1]
--                      return (math.exp(-((x-mu)/(math.sqrt(2)*sig))^2))*math.tanh(x-0.5)
--                      return (0.90-x)*math.tanh(x-0.5)
--                    return (0.90-x)*math.sin(4.0*2.0*math.pi*x)
--                    return (math.exp(-x)-0.5)*math.sin(4.0*2.0*math.pi*x)
                    return (math.exp(-x))*math.sin(4.0*2.0*math.pi*x)
                 end
   }
   initMom1i:advance(0.,0.,{},{mom1A})
   -- Analytic parallel flow speed.
   local initUi = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = basis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x = xn[1]
--                    return math.exp(-((x-mu)/(math.sqrt(2)*sig))^2)
--                    return math.tanh(x-0.5)
                    return math.sin(4.0*2.0*math.pi*x)
                 end
   }
   initUi:advance(0.,0.,{},{UiA})

   -- Compute parallel frow speed.
   local calcUi = Updater.CartFieldBinOp {
      onGrid     = confGrid,
      weakBasis  = basis,
      operation  = "Divide",
   }

   print("Computing parallel flow speed...")
   local t1 = os.clock()
   calcUi:advance(0.,0.,{numDens,mom1A},{Ui})
   local t2 = os.clock()
   io.write("Ui computation took total of ", t2-t1, " s\n")

   local calcMom1 = Updater.CartFieldBinOp {
      onGrid     = confGrid,
      weakBasis  = basis,
      operation  = "Multiply",
   }
   print("Computing parallel momentum...")
   local t1 = os.clock()
   calcMom1:advance(0.,0.,{numDens,UiA},{mom1})
   local t2 = os.clock()
   io.write("Mom1 computation took total of ", t2-t1, " s\n")

   local errU = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost = {1, 1}
   }
   local errMom1 = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost = {1, 1}
   }

   errU:combine(1.0, UiA, -1.0, Ui)
   errMom1:combine(1.0, mom1A, -1.0, mom1)

   numDens:write("numDens-1d.bp", 0.0)
   mom1A:write("mom1A-1d.bp", 0.0)
   mom1:write("mom1-1d.bp", 0.0)
   Ui:write("Ui-1d.bp", 0.0)
   UiA:write("UiA-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = confGrid,
      basis         = basis,
      numComponents = 1,
      quantity      = "V2"
   }
   local dynVecUi = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {errU}, {dynVecUi})
   local tmUi, lvUi = dynVecUi:lastData()
   io.write("Average RMS in Ui error = ", math.sqrt(lvUi[1]), "\n")

   local dynVecMom1 = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {errMom1}, {dynVecMom1})
   local tmMom1, lvMom1 = dynVecMom1:lastData()
   io.write("Average RMS in Mom1 error = ", math.sqrt(lvMom1[1]), "\n")
   print()
   return math.sqrt(lvUi[1]), math.sqrt(lvMom1[1])
end

------------------------------------------------------------------------------------
---                    Test in 2x
------------------------------------------------------------------------------------
function test_binOp2x(nx, ny, p, writeMatrix)
   writeMatrix = writeMatrix or false
   -- Phase-space and config-space grids.
   local confGrid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {nx, ny},
   }
   -- Basis functions.
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = p }
   io.write("nx=",nx," ny=",ny," polyOrder=", p, "\n")

   -- Zeroth and first (parallel) moments.
   local numDens = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1},
   }
   local mom1A = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1},
   }
   local mom1 = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1},
   }
   -- Calculated and analytica parallel flow speed.
   local Ui = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1},
   }
   local UiA = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1},
   }

   -- Initialize number density and first (parallel) moment.
   local initDens = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = basis,
      evaluate = function (t,xn)
                    local mu   = 0.5
                    local sig  = 0.04
                    local x, y = xn[1], xn[2]
--                    return math.exp(-((x-mu)/(math.sqrt(2)*sig))^2)
--                    return 0.90-x
--                    return math.exp(-x)-0.5
                    return math.exp(-x)*math.exp(-((y-mu)/(math.sqrt(2)*sig))^2)
                 end
   }
   initDens:advance(0.,0.,{},{numDens})
   local initMom1i = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = basis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x, y = xn[1], xn[2]
--                    return (math.exp(-((x-mu)/(math.sqrt(2)*sig))^2))*math.tanh(x-0.5)
--                    return (0.90-x)*math.tanh(x-0.5)
--                    return (0.90-x)*math.sin(4.0*2.0*math.pi*x)
--                    return (math.exp(-x)-0.5)*math.sin(4.0*2.0*math.pi*x)
--                    return (math.exp(-x))*math.sin(4.0*2.0*math.pi*x)
                    return math.exp(-x)*math.exp(-((y-mu)/(math.sqrt(2)*sig))^2)*
                           math.sin(4.0*2.0*math.pi*x)  
                 end
   }
   initMom1i:advance(0.,0.,{},{mom1A})
   -- Analytic parallel flow speed.
   local initUi = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = basis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x, y = xn[1], xn[2]
--                    return math.exp(-((x-mu)/(math.sqrt(2)*sig))^2)
--                    return math.tanh(x-0.5)
                    return math.sin(4.0*2.0*math.pi*x)
                 end
   }
   initUi:advance(0.,0.,{},{UiA})

   -- Compute parallel frow speed.
   local calcUi = Updater.CartFieldBinOp {
      onGrid     = confGrid,
      weakBasis  = basis,
      operation  = "Divide",
   }

   print("Computing parallel flow speed...")
   local t1 = os.clock()
   calcUi:advance(0.,0.,{numDens,mom1A},{Ui})
   local t2 = os.clock()
   io.write("Ui computation took total of ", t2-t1, " s\n")

   local calcMom1 = Updater.CartFieldBinOp {
      onGrid     = confGrid,
      weakBasis  = basis,
      operation  = "Multiply",
   }
   print("Computing parallel momentum...")
   local t1 = os.clock()
   calcMom1:advance(0.,0.,{numDens,UiA},{mom1})
   local t2 = os.clock()
   io.write("Mom1 computation took total of ", t2-t1, " s\n")

   local errU = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost = {1, 1}
   }
   local errMom1 = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = basis:numBasis(),
      ghost = {1, 1}
   }

   errU:combine(1.0, UiA, -1.0, Ui)
   errMom1:combine(1.0, mom1A, -1.0, mom1)

   numDens:write("numDens-1d.bp", 0.0)
   mom1A:write("mom1A-1d.bp", 0.0)
   mom1:write("mom1-1d.bp", 0.0)
   Ui:write("Ui-1d.bp", 0.0)
   UiA:write("UiA-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = confGrid,
      basis         = basis,
      numComponents = 1,
      quantity      = "V2"
   }
   local dynVecUi = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {errU}, {dynVecUi})
   local tmUi, lvUi = dynVecUi:lastData()
   io.write("Average RMS in Ui error = ", math.sqrt(lvUi[1]), "\n")

   local dynVecMom1 = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {errMom1}, {dynVecMom1})
   local tmMom1, lvMom1 = dynVecMom1:lastData()
   io.write("Average RMS in Mom1 error = ", math.sqrt(lvMom1[1]), "\n")
   print()
   return math.sqrt(lvUi[1]), math.sqrt(lvMom1[1])
end

------------------------------------------------------------------------------------
---                    Test in 1x1V
------------------------------------------------------------------------------------
function test_binOp1x1v(nx, nv, p, writeMatrix)
   writeMatrix = writeMatrix or false
   -- Phase-space and config-space grids.
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0, 6.0},
      cells = {nx, nv},
   }
   local confGrid = Grid.RectCart {
      lower = { phaseGrid:lower(1) },
      upper = { phaseGrid:upper(1) },
      cells = { phaseGrid:numCells(1) },
   }
   -- Basis functions.
   local phaseBasis = Basis.CartModalSerendipity { ndim = 2, polyOrder = p }
   local confBasis  = Basis.CartModalSerendipity { ndim = 1, polyOrder = p }
   io.write("nx=",nx," nv=",nv," polyOrder=", p, "\n")

   -- Phase space distribution function.
   local distF = DataStruct.Field {
      onGrid        = phaseGrid,
      numComponents = phaseBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Analytic and calculated numDens*distF.
   local ndistF = DataStruct.Field {
      onGrid        = phaseGrid,
      numComponents = phaseBasis:numBasis(),
      ghost         = {1, 1},
   }
   local ndistFA = DataStruct.Field {
      onGrid        = phaseGrid,
      numComponents = phaseBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Zeroth and first (parallel) moments.
   local numDens = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }
   local mom1A = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }
   local mom1 = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Calculated and analytica parallel flow speed.
   local Ui = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }
   local UiA = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }

   -- Initialize distribution function. We will (weak) multiply this by numDens.
   local initDistF = Updater.ProjectOnBasis {
      onGrid   = phaseGrid,
      basis    = phaseBasis,
      evaluate = function (t,xn)
                    local muV  = 0.0
                    local sigV = 0.8
                    local x, v = xn[1], xn[2]
                    return math.exp(-((v-muV)/(math.sqrt(2)*sigV))^2)
--                    return 0.90-x
--                    return math.exp(-x)-0.5
--                    return math.exp(-x)
                 end
   }
   initDistF:advance(0.,0.,{},{distF})
   -- Initialize number density and first (parallel) moment.
   local initDens = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = confBasis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x = xn[1]
--                    return math.exp(-((x-mu)/(math.sqrt(2)*sig))^2)
--                    return 0.90-x
--                    return math.exp(-x)-0.5
                    return math.exp(-x)
                 end
   }
   initDens:advance(0.,0.,{},{numDens})
   local initMom1i = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = confBasis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x = xn[1]
--                      return (math.exp(-((x-mu)/(math.sqrt(2)*sig))^2))*math.tanh(x-0.5)
--                      return (0.90-x)*math.tanh(x-0.5)
--                    return (0.90-x)*math.sin(4.0*2.0*math.pi*x)
--                    return (math.exp(-x)-0.5)*math.sin(4.0*2.0*math.pi*x)
                    return (math.exp(-x))*math.sin(4.0*2.0*math.pi*x)
                 end
   }
   initMom1i:advance(0.,0.,{},{mom1A})
   -- Analytic numDens*distF.
   local initnDistFA = Updater.ProjectOnBasis {
      onGrid   = phaseGrid,
      basis    = phaseBasis,
      evaluate = function (t,xn)
                    local muV  = 0.0
                    local sigV = 0.8
                    local x, v = xn[1], xn[2]
                    return math.exp(-x)*math.exp(-((v-muV)/(math.sqrt(2)*sigV))^2)
--                    return 0.90-x
--                    return math.exp(-x)-0.5
--                    return math.exp(-x)
                 end
   }
   initnDistFA:advance(0.,0.,{},{ndistFA})
   -- Analytic parallel flow speed.
   local initUi = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = confBasis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x = xn[1]
--                    return math.exp(-((x-mu)/(math.sqrt(2)*sig))^2)
--                    return math.tanh(x-0.5)
                    return math.sin(4.0*2.0*math.pi*x)
                 end
   }
   initUi:advance(0.,0.,{},{UiA})

   -- Compute parallel flow speed.
   local fldDiv = Updater.CartFieldBinOp {
      onGrid     = confGrid,
      weakBasis  = phaseBasis,
      fieldBasis = confBasis,
      operation  = "Divide",
   }

   print("Computing parallel flow speed...")
   local t1 = os.clock()
   fldDiv:advance(0.,0.,{numDens,mom1A},{Ui})
   local t2 = os.clock()
   io.write("Ui computation took total of ", t2-t1, " s\n")

   local fldMult = Updater.CartFieldBinOp {
      onGrid     = confGrid,
      weakBasis  = phaseBasis,
      fieldBasis = confBasis,
      operation  = "Multiply",
   }
   print("Computing parallel momentum...")
   local t1 = os.clock()
   fldMult:advance(0.,0.,{numDens,UiA},{mom1})
   local t2 = os.clock()
   io.write("Mom1 computation took total of ", t2-t1, " s\n")

   print("Multiply numDens by distF...")
   local t1 = os.clock()
   fldMult:advance(0.,0.,{numDens,distF},{ndistF})
   local t2 = os.clock()
   io.write("Mom1 computation took total of ", t2-t1, " s\n")

   local errF = DataStruct.Field {
      onGrid        = phaseGrid,
      numComponents = phaseBasis:numBasis(),
      ghost = {1, 1}
   }
   local errU = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost = {1, 1}
   }
   local errMom1 = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost = {1, 1}
   }

   errF:combine(1.0, ndistFA, -1.0, ndistF)
   errU:combine(1.0, UiA, -1.0, Ui)
   errMom1:combine(1.0, mom1A, -1.0, mom1)

   ndistF:write("ndistF.bp", 0.0)
   ndistFA:write("ndistFA.bp", 0.0)
   numDens:write("numDens.bp", 0.0)
   mom1A:write("mom1A.bp", 0.0)
   mom1:write("mom1.bp", 0.0)
   Ui:write("Ui.bp", 0.0)
   UiA:write("UiA.bp", 0.0)

   local calcIntF = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = phaseGrid,
      basis         = phaseBasis,
      numComponents = 1,
      quantity      = "V2"
   }
   local dynVecF = DataStruct.DynVector { numComponents = 1 }
   calcIntF:advance(0.0, 0.0, {errF}, {dynVecF})
   local tmF, lvF = dynVecF:lastData()
   io.write("Average RMS in F error = ", math.sqrt(lvF[1]), "\n")

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = confGrid,
      basis         = confBasis,
      numComponents = 1,
      quantity      = "V2"
   }
   local dynVecUi = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {errU}, {dynVecUi})
   local tmUi, lvUi = dynVecUi:lastData()
   io.write("Average RMS in Ui error = ", math.sqrt(lvUi[1]), "\n")

   local dynVecMom1 = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, 0.0, {errMom1}, {dynVecMom1})
   local tmMom1, lvMom1 = dynVecMom1:lastData()
   io.write("Average RMS in Mom1 error = ", math.sqrt(lvMom1[1]), "\n")
   print()
   return math.sqrt(lvF[1]), math.sqrt(lvUi[1]), math.sqrt(lvMom1[1])
end

------------------------------------------------------------------------------------
---                    Test in 1x2v
------------------------------------------------------------------------------------
function test_binOp1x2v(nx, nvx, nvy, p, writeMatrix)
   writeMatrix = writeMatrix or false
   -- Phase-space and config-space grids.
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -1.0, -1.0},
      upper = {1.0, 1.0, 1.0},
      cells = {nx, nvx, nvy},
   }
   local confGrid = Grid.RectCart {
      lower = { phaseGrid:lower(1) },
      upper = { phaseGrid:upper(1) },
      cells = { phaseGrid:numCells(1) },
   }
   -- Basis functions.
   local phaseBasis = Basis.CartModalSerendipity { ndim = 3, polyOrder = p }
   local confBasis  = Basis.CartModalSerendipity { ndim = 1, polyOrder = p }
   io.write("nx=",nx," nvx=",nvx," nvy=",nvy," polyOrder=", p, "\n")

   -- Zeroth and first (parallel) moments.
   local numDens = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }
   local mom1A = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis()*2,
      ghost         = {1, 1},
   }
   local mom1 = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis()*2,
      ghost         = {1, 1},
   }
   -- Calculated and analytic parallel flow speed.
   local Ui = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis()*2,
      ghost         = {1, 1},
   }
   local UiA = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis()*2,
      ghost         = {1, 1},
   }

   -- Initialize number density and first (parallel) moment.
   local initDens = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = confBasis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x = xn[1]
--                    return math.exp(-((x-mu)/(math.sqrt(2)*sig))^2)
--                    return 0.90-x
--                    return math.exp(-x)-0.5
                    return math.exp(-x)
                 end
   }
   initDens:advance(0.,0.,{},{numDens})
   local initMom1i = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = confBasis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x = xn[1]
--                      return (math.exp(-((x-mu)/(math.sqrt(2)*sig))^2))*math.tanh(x-0.5)
--                      return (0.90-x)*math.tanh(x-0.5)
--                    return (0.90-x)*math.sin(4.0*2.0*math.pi*x)
--                    return (math.exp(-x)-0.5)*math.sin(4.0*2.0*math.pi*x)
                    return (math.exp(-x))*math.sin(4.0*2.0*math.pi*x), 
                           (math.exp(-x))*math.cos(4.0*2.0*math.pi*x)

                 end
   }
   initMom1i:advance(0.,0.,{},{mom1A})
   -- Analytic parallel flow speed.
   local initUi = Updater.ProjectOnBasis {
      onGrid   = confGrid,
      basis    = confBasis,
      evaluate = function (t,xn)
                    local mu  = 0.5
                    local sig = 0.04
                    local x = xn[1]
--                    return math.exp(-((x-mu)/(math.sqrt(2)*sig))^2)
--                    return math.tanh(x-0.5)
                    return math.sin(4.0*2.0*math.pi*x), 
                           math.cos(4.0*2.0*math.pi*x)
                 end
   }
   initUi:advance(0.,0.,{},{UiA})

   -- Compute parallel frow speed.
   local calcUi = Updater.CartFieldBinOp {
      onGrid     = confGrid,
      weakBasis  = phaseBasis,
      fieldBasis = confBasis,
      operation  = "Divide",
   }

   print("Computing parallel flow speed...")
   local t1 = os.clock()
   calcUi:advance(0.,0.,{numDens,mom1A},{Ui})
   local t2 = os.clock()
   io.write("Ui computation took total of ", t2-t1, " s\n")

   local calcMom1 = Updater.CartFieldBinOp {
      onGrid     = confGrid,
      weakBasis  = phaseBasis,
      fieldBasis = confBasis,
      operation  = "Multiply",
   }
   print("Computing parallel momentum...")
   local t1 = os.clock()
   calcMom1:advance(0.,0.,{numDens,UiA},{mom1})
   local t2 = os.clock()
   io.write("Mom1 computation took total of ", t2-t1, " s\n")

   local errU = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis()*2,
      ghost = {1, 1}
   }
   local errMom1 = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis()*2,
      ghost = {1, 1}
   }

   errU:combine(1.0, UiA, -1.0, Ui)
   errMom1:combine(1.0, mom1A, -1.0, mom1)

   numDens:write("numDens-1d.bp", 0.0)
   mom1A:write("mom1A-1d.bp", 0.0)
   mom1:write("mom1-1d.bp", 0.0)
   Ui:write("Ui-1d.bp", 0.0)
   UiA:write("UiA-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = confGrid,
      basis         = confBasis,
      numComponents = 2,
      quantity      = "V2"
   }
   local dynVecUi = DataStruct.DynVector { numComponents = 2 }
   calcInt:advance(0.0, 0.0, {errU}, {dynVecUi})
   local tmUi, lvUi = dynVecUi:lastData()
   io.write("Average RMS in Ui error = (", math.sqrt(lvUi[1]),",", math.sqrt(lvUi[2]),") \n")

   local dynVecMom1 = DataStruct.DynVector { numComponents = 2 }
   calcInt:advance(0.0, 0.0, {errMom1}, {dynVecMom1})
   local tmMom1, lvMom1 = dynVecMom1:lastData()
   io.write("Average RMS in Mom1 error = (", math.sqrt(lvMom1[1]),",", math.sqrt(lvMom1[2]),") \n")
   print()
   return math.sqrt(lvUi[1]), math.sqrt(lvUi[2]), math.sqrt(lvMom1[1]), math.sqrt(lvMom1[2])
end

-------------------------------------------------------
-- Call tests in 1x for various resolutions
-------------------------------------------------------
function binOp1x_conv(p)
  print(" ")
  print("--- Testing convergence of 1x BinOp updater with p=",p," ---")
  errUi1, errMom11 = test_binOp1x(32, p)
  errUi2, errMom12 = test_binOp1x(64, p)
  errUi3, errMom13 = test_binOp1x(128, p)
  print("Division order:", math.log10(errUi1/errUi2)/math.log10(2.0), 
                           math.log10(errUi2/errUi3)/math.log10(2.0))
  print("Multiply order:", math.log10(errMom11/errMom12)/math.log10(2.0), 
                           math.log10(errMom12/errMom13)/math.log10(2.0))
--  assert_close(1.0, err1/err2/4.0, .01)
--  assert_close(1.0, err2/err3/4.0, .01)
  print()
end

-------------------------------------------------------
-- Call tests in 2x for various resolutions
-------------------------------------------------------
function binOp2x_conv(p)
  print(" ")
  print("--- Testing convergence of 2x BinOp updater with p=",p," ---")
  errUi1, errMom11 = test_binOp2x(32, 32, p)
  errUi2, errMom12 = test_binOp2x(64, 64, p)
  errUi3, errMom13 = test_binOp2x(128, 128, p)
  print("Division order:", math.log10(errUi1/errUi2)/math.log10(2.0), 
                           math.log10(errUi2/errUi3)/math.log10(2.0))
  print("Multiply order:", math.log10(errMom11/errMom12)/math.log10(2.0), 
                           math.log10(errMom12/errMom13)/math.log10(2.0))
--  assert_close(1.0, err1/err2/4.0, .01)
--  assert_close(1.0, err2/err3/4.0, .01)
  print()
end

-------------------------------------------------------
-- Call tests in 1x1V for various resolutions
-------------------------------------------------------
function binOp1x1v_conv(p)
  print(" ")
  print("--- Testing convergence of BinOp updater with p=",p," ---")
  errF1, errUi1, errMom11 = test_binOp1x1v(32, 8, p)
  errF2, errUi2, errMom12 = test_binOp1x1v(64, 8, p)
  errF3, errUi3, errMom13 = test_binOp1x1v(128, 8, p)
  print("phaseMultiply order:", math.log10(errF1/errF2)/math.log10(2.0), 
                                math.log10(errF2/errF3)/math.log10(2.0))
  print("Division order:", math.log10(errUi1/errUi2)/math.log10(2.0), 
                           math.log10(errUi2/errUi3)/math.log10(2.0))
  print("confMultiply order:", math.log10(errMom11/errMom12)/math.log10(2.0), 
                               math.log10(errMom12/errMom13)/math.log10(2.0))
--  assert_close(1.0, err1/err2/4.0, .01)
--  assert_close(1.0, err2/err3/4.0, .01)
  print()
end

-------------------------------------------------------
-- Call tests in 1x2v for various resolutions
-------------------------------------------------------
function binOp1x2v_conv(p)
  print("--- Testing convergence of BinOp updater with p=",p," ---")
  errUi1x, errUi1y, errMom11x, errMom11y = test_binOp1x2v(32, 32, 32, p)
  errUi2x, errUi2y, errMom12x, errMom12y = test_binOp1x2v(64, 64, 64, p)
  errUi3x, errUi3y, errMom13x, errMom13y = test_binOp1x2v(128, 128, 128, p)
  print("Division order: (", math.log10(errUi1x/errUi2x)/math.log10(2.0),",", 
    math.log10(errUi1y/errUi2y)/math.log10(2.0),")  (", 
    math.log10(errUi2x/errUi3x)/math.log10(2.0),",", math.log10(errUi2x/errUi3x)/math.log10(2.0))
  print("Multiply order: (", math.log10(errMom11x/errMom12x)/math.log10(2.0),",", 
    math.log10(errMom11y/errMom12y)/math.log10(2.0),")  (", 
    math.log10(errMom12x/errMom13x)/math.log10(2.0),",", math.log10(errMom12x/errMom13x)/math.log10(2.0))
--  assert_close(1.0, err1/err2/4.0, .01)
--  assert_close(1.0, err2/err3/4.0, .01)
  print()
end

-- run tests
local t1 = os.clock()
binOp1x_conv(1)
-- binOp1x_conv(2)
-- binOp2x_conv(1)
-- binOp1x1v_conv(1)
-- binOp1x1v_conv(2)
-- binOp1x2v_conv(1)
-- binOp1x2v_conv(2)
local t2 = os.clock()

print()
if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
io.write("Total test time: ", t2-t1, " s\n")
