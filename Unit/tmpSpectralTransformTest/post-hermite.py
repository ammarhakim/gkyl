from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplColors
#style.use('/ihome/mana/documents/gkeyll/postgkyl.mplstyle')
style.use('/Users/mana/Documents/Research/Gkeyll/code/postProcessingScripts/postgkyl.mplstyle')
import postgkyl as pg
import sys
#sys.path.insert(0, '/ihome/mana/documents/gkeyll/postProcessingScripts/')
sys.path.insert(0, '/Users/mana/Documents/Research/Gkeyll/code/postProcessingScripts/')
import pgkylUtil as pgu
import matplotlib.colors as colors

#.
#.Make plots from 1x1v Hermite spectrum tests.
#.Manaure Francisquez.
#.June 2019.
#.

outDir    = './'                      #.Output directory where figures are saved.
fileName  = 'hermite-test-1x1v-p1'    #.Root name of files to process.
polyOrder = 1

outFigureFile    = False        #.Output a figure file?.
figureFileFormat = '.png'       #[ Can be .png, .pdf, .ps, .eps, .svg.
#[ For data that changes by orders of magnitude in time one
#[ can request log color scale (default=False).
logColorBar      = True

plotxmFrames     = False         #.Plot x-m space in time (m is the Hermite index).
plotmSpectrum    = True         #.Plot m-spectrum sampled at x=samplePoint.

samplePoint = [0.5]

#.Tests to plot and their collisionality.
useTests = [4]

xLabel = r'x'
yLabel = r'Hermite index, $m$'

xLabel = r'Hermite index, $m$'
yLabel = r'$|f_m|^2$'

#........................... END OF USER INPUTS (maybe) ...........................#

basis = 'ms'

#.Some RGB colors. These are MATLAB-like.
defaultBlue   = [0, 0.4470, 0.7410]
defaultOrange = [0.8500, 0.3250, 0.0980]
defaultGreen  = [0.4660, 0.6740, 0.1880]
defaultPurple = [0.4940, 0.1840, 0.5560]
#.Colors in a single array.
defaultColors = [defaultBlue,defaultOrange,defaultGreen,defaultPurple]

#.LineStyles in a single array.
lineStyles = ['-','--',':','-.']

#.Font sizes used in the rest of the file.
xyLabelFontSize = 16
xyTickFontSize  = 12
colorBarLabelFontSize = 16

# ~~~~~~~~~~~ Plot a figure for each x-m frame ~~~~~~~~~~~~~~ #
if (plotxmFrames):
  #.Prepare figure.
  figProp1a = [7,6]
  ax1aPos   = [0.15, 0.15, 0.68, 0.82]
  cax1aPos  = [0.85, 0.13, 0.03, 0.835]
  fig1      = plt.figure(figsize=(figProp1a[0],figProp1a[1]))
  ax1a      = fig1.add_axes(ax1aPos)
  cbar_ax1a = fig1.add_axes(cax1aPos)
  ax1a.set_xlabel(xLabel, fontsize=xyLabelFontSize)
  ax1a.set_ylabel(yLabel, fontsize=xyLabelFontSize, labelpad=-1)
  for tick in ax1a.xaxis.get_major_ticks():
    tick.label.set_fontsize(xyTickFontSize)
  for tick in ax1a.yaxis.get_major_ticks():
    tick.label.set_fontsize(xyTickFontSize)
  #.Colormaps recommended: inferno, magma, plasma, viridis, Greys (NOT jet!)
  colormap  = 'inferno'
  plt.set_cmap(colormap)

  nT = useTests[0]
  prefix = './s'+str(nT)+'/'

  #.Hermite spectrum in each x-cell.
  fName = prefix+fileName+'_neut_Hermite2_'    #.Complete file name.
  
  nFrames                      = pgu.findLastFrame(fName,'bp')                #.Establish the number of frames.
  x, dim, nx, lx, dx           = pgu.getRawGrid(fName+'0.bp')                 #.Get the grid from file.
  xIn, dimIn, nxIn, lxIn, dxIn = pgu.getGrid(fName+'0.bp',polyOrder,basis)    #.Get the interpolated grid from file.

  #.Cell centers.
  xc   = 0.5*(x[0][0:-1]+x[0][1:])
  xcIn = 0.5*(xIn[0][0:-1]+xIn[0][1:])

  #.Will re-organize data into an object containing the cell in the first
  #.dimension, the hermite number in the second, and the expansion
  #.coefficients along x in the third.
  mModes  = (polyOrder+1)*(nx[1]-1)
  nSurfB  = polyOrder+1
  xmFld   = np.zeros((nx[0]-1,mModes,nSurfB))
  xmFldIn = np.zeros((nxIn[0]-1,mModes))

  #.Spatial array used for pcolormesh, which expects X, M to be one larger
  #.than xmFldIn in each direction.
  mSpace = np.linspace(0,mModes,mModes+1)
  X, M   = np.meshgrid(xIn[0],mSpace,indexing='ij')

  for nF in range(1):
    pgData = pgu.getRawData(fName+str(nF)+'.bp')    #.Read data with pgkyl.

    #.Reorganize data.
    for i in range(nx[0]-1):
      for j in range(mModes//(polyOrder+1)):
        for k in range(nSurfB):
#          xmFld[i,2*j,k]   = pgData[i,j,k]
#          xmFld[i,2*j+1,k] = pgData[i,j,nSurfB+k]
          xmFld[i,2*j,k]   = pgData[i,j,nSurfB*k]
          xmFld[i,2*j+1,k] = pgData[i,j,nSurfB*k+1]

    #.Calculate data interpolated along x.
    for i in range(nx[0]-1):
      for j in range(mModes):
        xmFldIn[(polyOrder+1)*i:(polyOrder+1)*(i+1),j] = pgu.evalF1x_e(np.squeeze(xmFld[i,j,:]),xcIn[(polyOrder+1)*i:(polyOrder+1)*(i+1)],xc[i],dx[0],polyOrder)


    if logColorBar:
      xmFldInSq = xmFldIn*xmFldIn
      valMin = np.amin(xmFldInSq)
      valMax = np.amax(xmFldInSq)
      hpl1a = ax1a.pcolormesh(X, M, xmFldInSq,
                              norm=colors.SymLogNorm(linthresh=0.03,linscale=0.03,
                              vmin=valMin,vmax=valMax),cmap='RdBu_r')
      cbara = plt.colorbar(hpl1a,ax=ax1a,cax=cbar_ax1a,extend='both')
      cbara.set_label(r'$\log |F(x,m)|^2$', rotation=270, labelpad=18, fontsize=colorBarLabelFontSize)
    else:
      valMin = np.amin(xmFldIn)
      valMax = np.amax(xmFldIn)
      hpl1a = ax1a.pcolormesh(X, M, xmFldIn)
      cbara = plt.colorbar(hpl1a,ax=ax1a,cax=cbar_ax1a)
      cbara.set_label(r'$F(x,m)$', rotation=270, labelpad=18, fontsize=colorBarLabelFontSize)
  
    if outFigureFile:
      plt.savefig(outDir+fileName+'xmFldSq_'+str(nF)+figureFileFormat)
    else:
      plt.show()


# ~~~~~~~~~~~ Plot a the Hermite spectrum at x=samplePoint ~~~~~~~~~~~~~~ #
if (plotmSpectrum):
  #.Prepare figure.
  figProp1a = [6,4]
  ax1aPos   = [0.15, 0.15, 0.78, 0.82]
#  cax1aPos  = [0.85, 0.13, 0.03, 0.835]
  fig1      = plt.figure(figsize=(figProp1a[0],figProp1a[1]))
  ax1a      = fig1.add_axes(ax1aPos)
#  cbar_ax1a = fig1.add_axes(cax1aPos)
  ax1a.set_xlabel(xLabel, fontsize=xyLabelFontSize)
  ax1a.set_ylabel(yLabel, fontsize=xyLabelFontSize, labelpad=-1)
  for tick in ax1a.xaxis.get_major_ticks():
    tick.label.set_fontsize(xyTickFontSize)
  for tick in ax1a.yaxis.get_major_ticks():
    tick.label.set_fontsize(xyTickFontSize)
  #.Colormaps recommended: inferno, magma, plasma, viridis, Greys (NOT jet!)
  colormap  = 'inferno'
  plt.set_cmap(colormap)

  nT = useTests[0]
  prefix = './s'+str(nT)+'/'

  #.Hermite spectrum in each x-cell.
  fName = prefix+fileName+'_neut_Hermite2_'    #.Complete file name.
  
  nFrames                      = pgu.findLastFrame(fName,'bp')                #.Establish the number of frames.
  x, dim, nx, lx, dx           = pgu.getRawGrid(fName+'0.bp')                 #.Get the grid from file.
  xIn, dimIn, nxIn, lxIn, dxIn = pgu.getGrid(fName+'0.bp',polyOrder,basis)    #.Get the interpolated grid from file.

  #.Cell centers.
  xc   = 0.5*(x[0][0:-1]+x[0][1:])
  xcIn = 0.5*(xIn[0][0:-1]+xIn[0][1:])

  #.Will re-organize data into an object containing the cell in the first
  #.dimension, the hermite number in the second, and the expansion
  #.coefficients along x in the third.
  mModes  = (polyOrder+1)*(nx[1]-1)
  nSurfB  = polyOrder+1
  xmFld   = np.zeros((nx[0]-1,mModes,nSurfB))
  xmFldIn = np.zeros((nxIn[0]-1,mModes))

  #.Spatial array used for pcolormesh, which expects X, M to be one larger
  #.than xmFldIn in each direction.
  mSpace = np.linspace(0,mModes,mModes+1)
  X, M   = np.meshgrid(xIn[0],mSpace,indexing='ij')

  for nF in range(1):
    pgData = pgu.getRawData(fName+str(nF)+'.bp')    #.Read data with pgkyl.

    #.Reorganize data. Notice that the order assumes that the map in getSolution
    #.is column-major: elements are written from Eigen matrix to Gkeyll CartField
    #.C++ array column-by-column.
    for i in range(nx[0]-1):
      for j in range(mModes//(polyOrder+1)):
        for k in range(nSurfB):
          xmFld[i,2*j,k]   = pgData[i,j,nSurfB*k]
          xmFld[i,2*j+1,k] = pgData[i,j,nSurfB*k+1]

    #.Calculate data interpolated along x.
    for i in range(nx[0]-1):
      for j in range(mModes):
        xmFldIn[(polyOrder+1)*i:(polyOrder+1)*(i+1),j] = pgu.evalF1x_e(np.squeeze(xmFld[i,j,:]),xcIn[(polyOrder+1)*i:(polyOrder+1)*(i+1)],xc[i],dx[0],polyOrder)

    iX = pgu.findNearestIndex(xcIn,samplePoint[0])

    xmFldInSq = xmFldIn*xmFldIn
    valMin = np.amin(xmFldInSq)
    valMax = np.amax(xmFldInSq)
    hpl1a = ax1a.semilogy(mSpace[0:-1], xmFldInSq[iX,:])
    plt.axis((mSpace[0],mSpace[-2],valMin*0.9,valMax*1.1))
  
    if outFigureFile:
      plt.savefig(outDir+fileName+'xmFldSq_'+str(nF)+figureFileFormat)
    else:
      plt.show()
