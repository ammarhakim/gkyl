#....................................................#
#.
#.pgkylUtil.py
#.Manaure Francisquez.
#.September 2018.
#.
#.This file contains functions and operators used by 
#.other scripts for post-processing Gkeyll data.
#.
#.
#....................................................#

#.Import libraries.
import postgkyl as pg
import numpy as np
import adios as ad
#.These are used for creating directories.
import os
from os import path
import errno
import shutil

sqrt2    = np.sqrt(2.0)
rsqrt2   = np.sqrt(2.0)/2.0
rsqrt2Cu = 1.0/np.sqrt(2.0**3)
sqrt3    = np.sqrt(3.0)
sqrt3d2  = np.sqrt(3.0/2.0)
sqrt5    = np.sqrt(5.0)
sqrt7    = np.sqrt(7.0)

#.Function to check existence of file/directory.......#
def checkDir(dirIn):
  if os.path.exists(os.path.dirname(dirIn)):
     return True
  else:
     return False

#.Function to check existence and/or make directory.......#
def checkMkdir(dirIn):
  if not os.path.exists(os.path.dirname(dirIn)):
    try:
      os.makedirs(os.path.dirname(dirIn))
    except OSError as exc: # Guard against race condition
      if exc.errno != errno.EEXIST:
        raise

#.Establish last frame outputted (useful for unfinished runs):
def findLastFrame(absFileName,fileExt):
  #.Input is the file name with its absolute address attached to it.
  #.Indicate file type: fileExt='bp' for ADIOS files, ='h5' for HDF5.
  cF         = 0
  moreFrames = os.path.isfile(absFileName+str(cF+1)+'.'+fileExt)
  while moreFrames:
    cF = cF+1
    moreFrames = os.path.isfile(absFileName+str(cF+1)+'.'+fileExt)
  return cF

#.Obtain the true grid (not for interpolated data..........#
def getRawGrid(dataFile,**opKey):
  pgData = pg.GData(dataFile)     #.Read data with pgkyl.
  dimOut = pgData.getNumDims()
  xNodal = pgData.getGrid()

  #.If desired, output cell center values of grid coordinates instead of nodal coordinates.
  if 'location' in opKey:
    if opKey['location']=='center':
      xOut = [[] for i in range(dimOut)]
      for i in range(dimOut):
        nNodes  = np.shape(xNodal[i])[0]
        xOut[i] = np.zeros(nNodes-1)
        xOut[i] = np.multiply(0.5,xNodal[i][0:nNodes-1]+xNodal[i][1:nNodes])
    else:
      xOut = xNodal
  else:
    xOut = xNodal

  nxOut = np.zeros(dimOut,dtype='int')
  lxOut = np.zeros(dimOut,dtype='double')
  dxOut = np.zeros(dimOut,dtype='double')
  for i in range(dimOut):
    nxOut[i] = np.size(xOut[i])
    lxOut[i] = xOut[i][-1]-xOut[i][0]
    dxOut[i] = xOut[i][ 1]-xOut[i][0]

  return xOut, dimOut, nxOut, lxOut, dxOut

#.Establish the grid......................................#
def getGrid(dataFile,p,basisType,**opKey):
  pgData         = pg.GData(dataFile)                       #.Read data with pgkyl.
  pgInterp       = pg.GInterpModal(pgData, p, basisType)    #.Interpolate data.
  xNodal, dataInterp = pgInterp.interpolate()
  dimOut         = np.shape(xNodal)[0]			    #.Number of dimensions in data.

  #.If desired, output cell center values of grid coordinates instead of nodal coordinates.
  if 'location' in opKey:
    if opKey['location']=='center':
      xOut = [[] for i in range(dimOut)]
      for i in range(dimOut): 
        nNodes  = np.shape(xNodal[i])[0]
        xOut[i] = np.zeros(nNodes-1)
        xOut[i] = np.multiply(0.5,xNodal[i][0:nNodes-1]+xNodal[i][1:nNodes])
    else:
      xOut = xNodal
  else:
    xOut = xNodal

  nxOut = np.zeros(dimOut,dtype='int')
  lxOut = np.zeros(dimOut,dtype='double')
  dxOut = np.zeros(dimOut,dtype='double')
  for i in range(dimOut):
    nxOut[i] = np.size(xOut[i])
    lxOut[i] = xOut[i][-1]-xOut[i][0]
    dxOut[i] = xOut[i][ 1]-xOut[i][0]
  return xOut, dimOut, nxOut, lxOut, dxOut

#.Obtain raw DG data.....................................#
def getRawData(dataFile):
  pgData  = pg.GData(dataFile)     #.Read data with pgkyl.
  dataOut = pgData.popValues()
  return dataOut

#.Interpolate DG data.....................................#
def getInterpData(dataFile,p,basisType,**opKey):
  pgData        = pg.GData(dataFile)                     #.Read data with pgkyl.
  pgInterp      = pg.GInterpModal(pgData, p, basisType)    #.Interpolate data.
  if 'comp' in opKey:
    xOut, dataOut = pgInterp.interpolate(opKey['comp'])
  else:
    xOut, dataOut = pgInterp.interpolate()
  return dataOut

#.Read the time variable in file..........................#
def getTime(dataFile):
#.Extract the time from file.
  hF      = ad.file(dataFile)
  timeOut = hF['time'].read()
  hF.close()
  return timeOut

#.........................................................#
#.This function finds the index of the grid point nearest to a given fix value.
def findNearestIndex(array,value):
  return (np.abs(array-value)).argmin()
#...end of findNearestIndex function...#

#.Evaluate function expanded in 1x Serendipity basis at certain points.
#.Limited to 0<p<4.
def evalF1x_e(fIn,xE,xcIn,dxIn,pOrderIn):
  NxE = np.size(xE)
  fEs = np.zeros(NxE)
  if pOrderIn == 1:
    for i in range(NxE):
      fEs[i] = rsqrt2*fIn[0] + sqrt3d2*fIn[1]*(xE[i]-xcIn)/(0.5*dxIn)
  elif pOrderIn == 2:
    for i in range(NxE):
      xi = (xE[i]-xcIn)/(0.5*dxIn)
      fEs[i] = rsqrt2*fIn[0] + sqrt3d2*xi*fIn[1] + 3.0*sqrt5*rsqrt2Cu*(xi**2-1.0/3.0)*fIn[2]
  elif pOrderIn == 3:
    for i in range(NxE):
      fEs[i] = rsqrt2*fIn[0] + sqrt3d2*xi*fIn[1] + 3.0*sqrt5*rsqrt2Cu*(xi**2-1.0/3.0)*fIn[2] + 5.0*sqrt7*rsqrt2Cu*(xi**3-3.0*xi/5.0)*fIn[3]
  return fEs

#.........................................................#
#.Plot cell-wise polynomial for 0<p<4 in 1D plot.
def plotLocalPoly(axisIn,xNodal1D,fIn,pIn,**opKey):
  hpOut = [0]*(np.size(xNodal1D)-1)
  for i in range(np.size(xNodal1D)-1):
    dxLoc = xNodal1D[i+1]-xNodal1D[i]
    xcLoc = 0.5*(xNodal1D[i+1]+xNodal1D[i])
    xLoc  = [xNodal1D[i]]
    for p in range(pIn):
      xLoc.append(xNodal1D[i]+float(p+1)*dxLoc/float(pIn))
    yLoc  = evalF1x_e(fIn[i],xLoc,xcLoc,dxLoc,pIn)

    if 'lines' in opKey:
        opKey['lines'][i].set_data(xLoc,yLoc)
    else:
      if 'color' in opKey:
        hpOut[i], = axisIn.plot(xLoc, yLoc, color=opKey['color'])
      else:
        hpOut[i], = axisIn.plot(xLoc, yLoc)

  if 'lines' not in opKey:
    print(np.shape(hpOut))
    return hpOut

#.Return a variable name to put on the figure.............# 
def commonVarName(fileVarName,**opKey):
  fVnameSplit = fileVarName.split("_")
  if len(fVnameSplit)>1:
    species = fVnameSplit[0]
    var     = fVnameSplit[1]
  else:
    var     = fVnameSplit[0]

  if var == 'GkM0':
    if species == 'electron':
      varStrOut = 'n_e'
    else:
      varStrOut = 'n_i'
    unitStrOut  = ' (m$^{-3}$)'
  if var == 'GkM1':
    if species == 'electron':
      varStrOut = 'n_eu_{\parallel,e}'
    else:
      varStrOut = 'n_eu_{\parallel,i}'
    unitStrOut  = ' (m$^{-2}$/s)'
  if var == 'uPar':
    if species == 'electron':
      varStrOut = 'u_{\parallel,e}'
    else:
      varStrOut = 'u_{\parallel,i}'
    unitStrOut  = ' (m/s)'
  elif var == 'vthSq':
    if species == 'electron':
      varStrOut = 'T_e'
    else:
      varStrOut = 'T_i'
    unitStrOut  = ' (eV)'
  elif var == 'phi':
    varStrOut  = 'phi'
    unitStrOut = ' (V)'
  elif var == 'field':
    if opKey['comp']==0:
      varStrOut  = 'E_x'
      unitStrOut = ''
    elif opKey['comp']==1:
      varStrOut  = 'E_y'
      unitStrOut = ''
    elif opKey['comp']==2:
      varStrOut  = 'E_z'
      unitStrOut = ''
    elif opKey['comp']==3:
      varStrOut  = 'B_x'
      unitStrOut = ''
    elif opKey['comp']==4:
      varStrOut  = 'B_y'
      unitStrOut = ''
    elif opKey['comp']==5:
      varStrOut  = 'B_z'
      unitStrOut = ''

  return varStrOut, unitStrOut

#.........................................................#
#.This function reads the time average if it is already computed
#.and stored in a file, or computes a new one (and stores it in
#.a file if saveAv=True).
def getTimeAv(dataDir,simName,varName,iFrame,fFrame,p,b,saveAv,tAvDir):
  #.Check or create post data directory.
  checkMkdir(tAvDir)
  #.Check if time average file already exists.
  tAvFile = tAvDir+simName+'_'+varName+'_TimeAv'+str(iFrame)+'-'+str(fFrame)+'.bp' 
  if not os.path.isfile(tAvFile):
    #.Compute time average and store it in new file.
    fileName = dataDir+simName+'_'+varName+'_%d.bp'
    x, gridDim, nx, lx, dx = getGrid(fileName % iFrame,p,b,location='center')

    q0AvT = np.zeros(nx)
    for nFr in range(iFrame,fFrame+1):
      #.Read 3D data into q0.
      q0AvT = np.add(q0AvT,np.squeeze(getInterpData(fileName % nFr,p,b)))

    q0AvT = np.divide(q0AvT,float(fFrame-iFrame+1))

    if saveAv:
      #.Save time average to a file for reuse.
      print(" ")
      print("Saving time average in "+tAvFile+" ...")
      #.Function to write DG coefficients to Gkeyll-style ADIOS file.
      sNumCells  = ""
      sOffsets   = ""
      for i in range(np.size(nx)):
        sNumCells += "{:d},".format(int(nx[i]))
        sOffsets  += "0,"
      #.ADIOS init.
      ad.init_noxml()
      ad.set_max_buffer_size(1000)
      groupId = ad.declare_group("CartFieldInterp", "")
      ad.select_method(groupId, "POSIX1", "", "")
      #.Define variables and attributes.
      ad.define_attribute_byvalue(groupId, "numCells", "", nx)
      lo = np.zeros(np.size(nx), dtype='double')
      up = np.zeros(np.size(nx), dtype='double')
      for i in range(np.size(nx)):
        lo[i], up[i] = x[i][0], x[i][-1]
      ad.define_attribute_byvalue(groupId, "lowerBounds", "", lo)
      ad.define_attribute_byvalue(groupId, "upperBounds", "", up)
      ad.define_var(groupId, "CartGridFieldInterpTimeAv", "",
              ad.DATATYPE.double,
              sNumCells, sNumCells, sOffsets)
      fh = ad.open("CartFieldInterp", tAvFile, 'w')
      ad.write(fh, "CartGridFieldInterpTimeAv", q0AvT)
      ad.close(fh)
      ad.finalize()
      #.Deal with weird file output where a '.bp.0' file is created.
      if len(tAvFile.split('/')) > 1:
          nm = tAvFile.split('/')[-1]
      else:
          nm = tAvFile
      shutil.move(tAvFile + '.dir/' + nm + '.0', tAvFile)
      shutil.rmtree(tAvFile + '.dir')
  else:
    #.Read time average from existent file.
    print(" ")
    print("Reading time average in "+tAvFile+" ...")
    hF    = ad.file(tAvFile)
    q0AvT = hF['CartGridFieldInterpTimeAv'].read()
    hF.close()

  return q0AvT
  
#.Set minimum and maximum values in an array..............#
def setMinMax(aIn,minIn,maxIn):
  if np.amin(aIn)<minIn:
    minOut = np.amin(aIn)
  else:
    minOut = minIn
  if np.amax(aIn)>maxIn:
    maxOut = np.amax(aIn)
  else:
    maxOut = maxIn
  return minOut, maxOut

#.Derivative along X......................................#
def derX(aIn,dx,xBC,acc):
#.aIn: 2D field.
#.xBC: integer indicating boundary condition along X.
#.acc: accuracy, 2 for 2nd order, 4 for 4th order.
  s0 = aIn.shape[0]
  s1 = aIn.shape[1]
  ax = np.zeros((s0, s1))
  rdxd2  = 1.0/(2.0*dx)
  rdx2d3 = 2.0/(3.0*dx)
  rdxd12 = 1.0/(12.0*dx)

  aBu = np.zeros((s0+acc, s1+acc))
  aBu[2:s0+2,2:s1+2] = aIn
  if xBC == 1:
    #.Even symmetry.
    aBu[2:s0+2,0]    = aIn[:,3]
    aBu[2:s0+2,1]    = aIn[:,2]
    aBu[2:s0+2,s1+2] = aIn[:,s1-2]
    aBu[2:s0+2,s1+3] = aIn[:,s1-3]
  elif xBC == -1:
    #.Odd symmetry.
    aBu[2:s0+2,0]    = -aIn[:,3]
    aBu[2:s0+2,1]    = -aIn[:,2]
    aBu[2:s0+2,s1+2] = -aIn[:,s1-2]
    aBu[2:s0+2,s1+3] = -aIn[:,s1-3]

  for j in range(0,s0):
    for i in range(0,s1):
      ax[j,i]=(aBu[j+2,i+3]-aBu[j+2,i+1])*rdx2d3-(aBu[j+2,i+4]-aBu[j+2,i])*rdxd12
  return ax

#.Derivative along Y......................................#
def derY(aIn,dy,yBC,acc):
#.aIn: 2D field.
#.yBC: integer indicating boundary condition along Y.
#.acc: accuracy, 2 for 2nd order, 4 for 4th order.
  s0   = aIn.shape[0]
  s1   = aIn.shape[1]
  ay   = np.zeros([s0, s1])
  rdyd2  = 1.0/(2.0*dy)
  rdy2d3 = 2.0/(3.0*dy)
  rdyd12 = 1.0/(12.0*dy)

  aBu = np.zeros([s0+acc, s1+acc])
  aBu[2:s0+2,2:s1+2] = aIn
  if yBC == 9:
    #.Periodic.
    aBu[0,2:s1+2]    = aIn[s0-3,:]
    aBu[1,2:s1+2]    = aIn[s0-1,:]
    aBu[s0+2,2:s1+2] = aIn[0,:]
    aBu[s0+3,2:s1+2] = aIn[1,:]

  for j in range(0,s0):
    for i in range(0,s1):
      ay[j,i]=(aBu[j+3,i+2]-aBu[j+1,i+2])*rdy2d3-(aBu[j+4,i+2]-aBu[j,i+2])*rdyd12
  return ay

#.Average all of the dimIn dimension......................#
def avAllDim(fIn,dimIn):
#...fIn: 2D field.
#...dimIn: dimension to be averaged.
    fAv = np.mean(fIn, axis=dimIn)
    return fAv

