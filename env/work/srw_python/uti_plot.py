#############################################################################
# Simple 1D & 2D plotting utilities package (/ "wrapper" of featured scientific plotting libraries).
# Currently uses matplotlib, but may support other libraries in the future.
# v 0.01
#############################################################################

#numerical modules
from array import *
from math import *
import numpy as np

#plotting modules and plotting setup
import matplotlib as mpl

mpl.use('TkAgg')
import matplotlib.pyplot as pl
pl.ioff()

#system  modules
import time, os, sys

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

def __plot_1D__(ar1d, x_range, labels, fig, typ=111):
  lenAr1d = len(ar1d)
  if lenAr1d > x_range[2]: ar1d = np.array(ar1d[0:x_range[2]])
  else:
    if lenAr1d < x_range[2]:
      auxAr = array('d', [0]*x_range[2])
      for i in range(lenAr1d): auxAr[i] = ar1d[i]
      ar1d = np.array(auxAr)

  if isinstance(ar1d,(list,array)): ar1d = np.array(ar1d)

  x = np.linspace(x_range[0],x_range[1],x_range[2])
  ax = fig.add_subplot(typ)
  ax.plot(x,ar1d)
  ax.grid()
  ax.set_xlim(x[0],x[-1])
  ax.set_xlabel(labels[0])
  ax.set_ylabel(labels[1])
  if(len(labels) > 2):
    ax.set_title(labels[2])

#def plot1D(ar1d,x_range,labels=('energy, KeV','Ph/s/0.1%BW')):
def uti_plot1d(ar1d, x_range, labels=('energy [eV]', 'ph/s/0.1%bw')):
  #fig = pl.figure(figsize=(12,8))
  fig = pl.figure()
  __plot_1D__(ar1d,x_range,labels,fig)
  return fig

def __plot_2D__(ar2d, x_range, y_range, labels, fig, typ=111):
  totLen = int(x_range[2]*y_range[2])
  lenAr2d = len(ar2d)
  if lenAr2d > totLen: ar2d = np.array(ar2d[0:totLen])
  else:
    if lenAr2d < totLen:
      auxAr = array('d', [0]*lenAr2d)
      for i in range(lenAr2d): auxAr[i] = ar2d[i]
      ar2d = np.array(auxAr)
   
  #if isinstance(ar2d,(list,array)): ar2d = np.array(ar2d).reshape(x_range[2],y_range[2])
  if isinstance(ar2d,(list,array)): ar2d = np.array(ar2d)
  ar2d = ar2d.reshape(y_range[2],x_range[2])
  
  x = np.linspace(x_range[0],x_range[1],x_range[2])
  y = np.linspace(y_range[0],y_range[1],y_range[2])
  ax = fig.add_subplot(typ)
  #ax.pcolormesh(x,y,ar2d.T,cmap=pl.cm.Greys_r)
  ax.pcolormesh(x,y,ar2d,cmap=pl.cm.Greys_r)
  ax.set_xlim(x[0],x[-1])
  ax.set_ylim(y[0],y[-1])
  ax.set_xlabel(labels[0])
  ax.set_ylabel(labels[1])
  if(len(labels) > 2):
    ax.set_title(labels[2])

#def plot2D(ar2d,x_range,y_range,labels=('Horizontal position, mm','Vertical position, mm')):
def uti_plot2d(ar2d, x_range, y_range, labels=('Horizontal position [m]','Vertical position [m]')):
  #fig = pl.figure(figsize=(12,8))
  fig = pl.figure()
  __plot_2D__(ar2d,x_range,y_range,labels,fig)
  return fig

#def srw_ascii_load(fname):
def uti_data_file_load(_fname, _read_labels=1):
  nLinesHead = 11
  hlp = []
  #with open(_fname,'r') as f: hlp = f.readlines(nLinesHead)
  with open(_fname,'r') as f:
    for i in range(nLinesHead):
      hlp.append(f.readline())
  
  #ne, nx, ny, ns    = [   int(hlp[i].replace('#','').split()[0]) for i in [3,6,9,10] ]
  ne, nx, ny = [int(hlp[i].replace('#','').split()[0]) for i in [3,6,9]]
  ns = 1
  testStr = hlp[nLinesHead - 1]
  if testStr[0] == '#':
    ns = int(testStr.replace('#','').split()[0])

  e0,e1,x0,x1,y0,y1 = [float(hlp[i].replace('#','').split()[0]) for i in [1,2,4,5,7,8]]

  #data = np.squeeze( np.loadtxt(_fname, dtype=np.float64).reshape(ns,ny,nx,ne)[0] ) #get data from file.
  data = np.squeeze(np.loadtxt(_fname, dtype=np.float64)) #get data from file (C-aligned flat)
  
  allrange = e0, e1, ne, x0, x1, nx, y0, y1, ny

  arLabels = ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Intensity']
  arUnits = ['eV', 'm', 'm', 'ph/s/.1%bw/mm^2']
  
  if _read_labels:
    
    arTokens = hlp[0].split(' [')
    arLabels[3] = arTokens[0].replace('#','')
    arUnits[3] = '';
    if len(arTokens) > 1:
      arUnits[3] = arTokens[1].split('] ')[0]
    #print(arLabels[3], arUnits[3])

    for i in range(3):
      arTokens = hlp[i*3 + 1].split()
      nTokens = len(arTokens)
      nTokensLabel = nTokens - 3
      nTokensLabel_mi_1 = nTokensLabel - 1
      strLabel = ''
      for j in range(nTokensLabel):
        strLabel += arTokens[j + 2]
        if j < nTokensLabel_mi_1: strLabel += ' '
      arLabels[i] = strLabel
      arUnits[i] = arTokens[nTokens - 1].replace('[','').replace(']','')
  
  m = enum('T','V','H','E','HV','EV','EH','EHV')
   
  if ne==1 and nx==1 and ny==1 : mode = m.T
  if ne==1 and nx==1 and ny>1  : mode = m.V
  if ne==1 and nx>1  and ny==1 : mode = m.H
  if ne>1  and nx==1 and ny==1 : mode = m.E
  if ne==1 and nx>1  and ny>1  : mode = m.HV
  if ne>1  and nx==1 and ny>1  : mode = m.EV
  if ne>1  and nx>1  and ny==1 : mode = m.EH
  if ne>1  and nx>1  and ny>1  : mode = m.EHV

  #print(mode)
  return data, mode, allrange, arLabels, arUnits

def rescale(maxabsval,strval):
  mult = 1
  if maxabsval>=1.0e2 and maxabsval<1.0e5:
    mult = 1.0e-3
    strval = 'k'+strval
  if maxabsval>=1.0e5 and maxabsval<1.0e8:
    mult = 1.0e-6
    strval = 'M'+strval
  if maxabsval>=1.0e8 and maxabsval<1.0e11:
    mult = 1.0e-6
    strval = 'G'+strval
  if maxabsval>=1.0e-4 and maxabsval<1.0e-1:
    mult = 1.0e3
    strval = 'm'+strval
  if maxabsval>=1.0e-7 and maxabsval<1.0e-4:
    mult = 1.0e6
    strval = 'u'+strval
  if maxabsval>=1.0e-10 and maxabsval<1.0e-7:
    mult = 1.0e9
    strval = 'n'+strval
  if maxabsval>1.0e-13 and maxabsval<1.0e-10:
    mult = 1.0e12
    strval = 'p'+strval
  return mult, strval

def rescale_range(allrange, _ar_units, _ec=0, _xc=0, _yc=0):
  e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
  em = abs(np.array([e0,e1])).max()
  xm = abs(np.array([x0,x1])).max()
  ym = abs(np.array([y0,y1])).max()
  #mult_e, str_e = rescale(em,"eV")
  #mult_x, str_x = rescale(xm,"m")
  #mult_y, str_y = rescale(ym,"m")
  mult_e, str_e = rescale(em, _ar_units[0])
  mult_x, str_x = rescale(xm, _ar_units[1])
  mult_y, str_y = rescale(ym, _ar_units[2])
  allnewrange = e0*mult_e, e1*mult_e, ne, x0*mult_x, x1*mult_x, nx, y0*mult_y, y1*mult_y, ny, _ec*mult_e, _xc*mult_x, _yc*mult_y
  strval = str_e, str_x, str_y
  return allnewrange, strval

def __mode_T__(data, allrange, _ar_labels, _ar_units, _ec=0, _xc=0, _yc=0):
  #allrange, units = rescale_range(allrange)
  allrange, units = rescale_range(allrange, _ar_units, _ec, _xc, _yc)
  #e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
  e0, e1, ne, x0, x1, nx, y0, y1, ny, ec, xc, yc = allrange
  #toprint = (e0,e1,units[0], x0,x1,units[1], y0,y1,units[2], data[0]) #squeeze have reduced it to an array with one element.
  toprint = (e0,units[0], x0,units[1], y0,units[2], data[0],units[3]) #squeeze have reduced it to an array with one element.
  #sys.stdout.write('Total Flux for \nE: %f -> %f %s\nX: %f -> %f %s\nY: %f -> %f %s\n is %f Ph/s/0.1%BW' % toprint) 
  sys.stdout.write(_ar_labels[3] + ' for \n' + _ar_labels[0] + ': %f %s\n' + _ar_labels[1] + ': %f %s\n' + _ar_labels[2] + ': %f %s\n is %f %s' % toprint) 
  return None

def __mode_V__(data, allrange, _ar_labels, _ar_units):
  #allrange, units = rescale_range(allrange)
  allrange, units = rescale_range(allrange, _ar_units)
  #e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
  e0, e1, ne, x0, x1, nx, y0, y1, ny, ec, xc, yc = allrange
  range_y = y0, y1, ny
  #label = ("Vertical Position, ["+units[2]+"]","Ph/s/0.1%BW/mm^2")
  label = (_ar_labels[2] + ' [' + units[2] + ']', _ar_labels[3] + ' [' + _ar_units[3] + ']')
  fig = pl.figure(figsize=(4,4))
  __plot_1D__(data, range_y, label, fig)
  return fig

def __mode_H__(data, allrange, _ar_labels, _ar_units):
  #allrange, units = rescale_range(allrange)
  allrange, units = rescale_range(allrange, _ar_units)
  #e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
  e0, e1, ne, x0, x1, nx, y0, y1, ny, ec, xc, yc = allrange
  range_x = x0, x1, nx 
  #label = ("Horizontal Position, ["+units[1]+"]","Ph/s/0.1%BW/mm^2")
  label = (_ar_labels[1] + ' [' + units[1] + ']', _ar_labels[3] + ' [' + _ar_units[3] + ']')
  fig = pl.figure(figsize=(4,4))
  __plot_1D__(data, range_x, label, fig)
  return fig

def __mode_E__(data, allrange, _ar_labels, _ar_units):
  #allrange, units = rescale_range(allrange)
  allrange, units = rescale_range(allrange, _ar_units)
  e0, e1, ne, x0, x1, nx, y0, y1, ny, ec, xc, yc = allrange
  range_e = e0, e1, ne
  #label = ("Energy, ["+units[0]+"]","Ph/s/0.1%BW")
  label = (_ar_labels[0] + ' [' + units[0] + ']', _ar_labels[3] + ' [' + _ar_units[3] + ']')
  fig = pl.figure(figsize=(4,4))
  __plot_1D__(data,range_e,label,fig)
  return fig

def __mode_HV__(data, allrange, _ar_labels, _ar_units, _xc=0, _yc=0, _graphs_joined=1):
  #allrange, units = rescale_range(allrange)
  allrange, units = rescale_range(allrange, _ar_units, 0, _xc, _yc)
  e0, e1, ne, x0, x1, nx, y0, y1, ny, ec, xc, yc = allrange
  range_x = x0, x1, nx
  range_y = y0, y1, ny
  #print range_x, range_y
  #x = np.linspace(x0, x1, nx)
  #y = np.linspace(y0, y1, ny)
  #ix = np.where(abs(x)==abs(x).min())[0][0]
  #iy = np.where(abs(y)==abs(y).min())[0][0]
  #label2D = ("Horizontal Position, ["+units[1]+"]", "Vertical Position, ["+units[2]+"]")
  strTitle = _ar_labels[3]
  if (ne == 1) and (e0 > 0): strTitle += ' at ' + str(e0) + ' ' + units[0]
  label2D = (_ar_labels[1] + ' [' + units[1]+ ']', _ar_labels[2] + ' [' + units[2] + ']', strTitle)

  #label1H = ("Horizontal Position, ["+units[1]+"]","Ph/s/0.1%BW/mm^2")
  strTitle = 'At ' + _ar_labels[2] + ': ' + str(yc)
  if yc != 0: strTitle += ' ' + units[2]
  label1H = (_ar_labels[1] + ' [' + units[1] + ']', _ar_labels[3] + ' [' + _ar_units[3] + ']', strTitle)
  
  #label1V = ("Vertical Position, ["+units[2]+"]","Ph/s/0.1%BW/mm^2")
  strTitle = 'At ' + _ar_labels[1] + ': ' + str(xc)
  if xc != 0: strTitle += ' ' + units[1]
  label1V = (_ar_labels[2] + ' [' + units[2] + ']', _ar_labels[3] + ' [' + _ar_units[3] + ']', strTitle)

  fig = None
  if _graphs_joined:
    fig = pl.figure(figsize=(12,5))
    __plot_2D__(data, range_x, range_y, label2D, fig, 131) #showing graphs in one figure
  else: uti_plot2d(data, range_x, range_y, label2D)

  xStep = 0
  if nx > 1: xStep = (x1 - x0)/(nx - 1)
  yStep = 0
  if ny > 1: yStep = (y1 - y0)/(ny - 1)
  inperpOrd = 1 #interpolation order (1 to 3)
  
  #__plot_1D__(data[iy,:],range_x,label1H,fig,132)
  arCutX = array('d', [0]*nx)
  xx = x0
  for ix in range(nx):
    arCutX[ix] = uti_interp_2d(xx, yc, x0, xStep, nx, y0, yStep, ny, data, inperpOrd, 1, 0)
    xx += xStep
  if _graphs_joined: __plot_1D__(arCutX, range_x, label1H, fig, 132)
  else: uti_plot1d(arCutX, range_x, label1H)
  
  #__plot_1D__(data[:,ix],range_y,label1V,fig,133)
  arCutY = array('d', [0]*ny)
  yy = y0
  for iy in range(ny):
    arCutY[iy] = uti_interp_2d(xc, yy, x0, xStep, nx, y0, yStep, ny, data, inperpOrd, 1, 0)
    yy += yStep
  if _graphs_joined: __plot_1D__(arCutY, range_y, label1V, fig, 133)
  else: uti_plot1d(arCutY, range_y, label1V)
  
  return fig

#def __mode_EV__(data, allrange):
#  #to be updated
#  allrange, units = rescale_range(allrange)
#  e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
#  range_e = e0, e1, ne  
#  range_y = y0, y1, ny
#  e = np.linspace(e0, e1, ne)
#  y = np.linspace(y0, y1, ny)
#  ie = np.where(data.sum(axis=1)==data.sum(axis=1).max())[0][0]
#  iy = np.where(abs(y).min())[0][0]
#  label1E = ("Energy, ["+units[0]+"]","Ph/s/0.1%BW/mm^2")
#  label1V = ("Vertical Position, ["+units[2]+"]","Ph/s/0.1%BW/mm^2")
#  fig = pl.figure(figsize=(8,4))
#  __plot_1D__(data[iy,:],range_e,label1E,fig,121)
#  __plot_1D__(data[:,ie],range_y,label1V,fig,122)
#  return fig
  
#def __mode_EH__(data, allrange):
#  #to be updated
#  allrange, units = rescale_range(allrange)
#  e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
#  range_e = e0, e1, ne  
#  range_x = x0, x1, nx
#  e = np.linspace(e0, e1, ne)
#  x = np.linspace(x0, x1, nx)
#  ie = np.where(data.sum(axis=1)==data.sum(axis=1).max())[0][0]
#  ix = np.where(abs(y).min())[0][0]
#  label1E = ("Energy, ["+units[0]+"]","Ph/s/0.1%BW/mm^2")
#  label1H = ("Horizontal Position, ["+units[1]+"]","Ph/s/0.1%BW/mm^2")
#  fig = pl.figure(figsize=(8,4))
#  __plot_1D__(data[ix,:],range_e,label1E,fig,121)
#  __plot_1D__(data[:,ie],range_x,label1H,fig,122)
#  return fig

def __mode_EHV__(data, allrange, _ar_labels, _ar_units, _ec, _xc, _yc, _graphs_joined=1):
  #allrange, units = rescale_range(allrange)
  allrange, units = rescale_range(allrange, _ar_units, _ec, _xc, _yc)
  #e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
  e0, e1, ne, x0, x1, nx, y0, y1, ny, ec, xc, yc = allrange
  
  #e = np.linspace(e0, e1, ne)
  #x = np.linspace(x0, x1, nx)
  #y = np.linspace(y0, y1, ny)
  #ie = np.where(data.sum(axis=1)==data.sum(axis=1).max())[0][0]
  #ix = np.where(abs(x)==abs(x).min())[0][0]
  #iy = np.where(abs(y)==abs(y).min())[0][0]

  range_e = e0, e1, ne
  range_x = x0, x1, nx
  range_y = y0, y1, ny

  ie = 0
  if ne > 1:
    if ec > e1: ie = ne - 1
    elif ec > e0:
      eStep = (e1 - e0)/(ne - 1)
      if eStep > 0: ie = int(round((ec - e0)/eStep))
  ix = 0
  if nx > 1:
    if xc > x1: ix = nx - 1
    elif xc > x0:
      xStep = (x1 - x0)/(nx - 1)
      if xStep > 0: ix = int(round((xc - x0)/xStep))
  iy = 0
  if ny > 1:
    if yc > y1: iy = ny - 1
    elif yc > y0:
      yStep = (y1 - y0)/(ny - 1)
      if yStep > 0: iy = int(round((yc - y0)/yStep))
  
  #label2D = ("Horizontal Position, ["+units[1]+"]", "Vertical Position, ["+units[2]+"]")
  label2D = (_ar_labels[1] + ' [' + units[1] + ']', _ar_labels[2] + ' [' + units[2] + ']')
  
  #label1E = ("Energy, ["+units[0]+"]","Ph/s/0.1%BW/mm^2")
  label1E = (_ar_labels[0] + ' [' + units[0] + ']', _ar_labels[3] + ' [' + units[3] + ']')
  
  #label1H = ("Horizontal Position, ["+units[1]+"]","Ph/s/0.1%BW/mm^2")
  label1H = (_ar_labels[1] + ' [' + units[1] + ']', _ar_labels[3] + ' [' + units[3] + ']')
  
  #label1V = ("Vertical Position, ["+units[2]+"]","Ph/s/0.1%BW/mm^2")
  label1V = (_ar_labels[2] + ' [' + units[2] + ']', _ar_labels[3] + ' [' + units[3] + ']')

  arCutXY = array('d', [0]*nx*ny)
  perY = ne*nx
  i = 0
  for iiy in range(ny):
    perY_iiy = perY*iiy
    for iix in range(nx):
      arCutXY[i] = data[ie + ne*iix + perY_iiy]
      i += 1
      
  arCutE = array('d', [0]*ne)
  perX_ix = ne*ix
  perY_iy = perY*iy
  for iie in range(ne): arCutE[iie] = data[iie + perX_ix + perY_iy]

  arCutX = array('d', [0]*nx)
  for iix in range(nx): arCutX[iix] = data[ie + ne*iix + perY_iy]

  arCutY = array('d', [0]*ny)
  for iiy in range(ny): arCutY[iiy] = data[ie + perX_ix + perY*iiy]
  
  fig = pl.figure(figsize=(8,8))
  #__plot_2D__(data[:,:,ie], range_x, range_y, label2D, fig, 221)
  #__plot_1D__(data[ix,iy,:],range_e,label1E,fig,224)
  #__plot_1D__(data[ie,:,iy],range_x,label1H,fig,222)
  #__plot_1D__(data[:,ix,ie],range_y,label1V,fig,223)

  fig = None
  if _graphs_joined:
    fig = pl.figure(figsize=(12,5))
    __plot_2D__(arCutXY, range_x, range_y, label2D, fig, 221) #showing graphs in one figure
    __plot_1D__(arCutE, range_e, label1E, fig, 222)
    __plot_1D__(arCutX, range_x, label1X, fig, 223)
    __plot_1D__(arCutY, range_y, label1Y, fig, 224)
  else:
    uti_plot2d(arCutXY, range_x, range_y, label2D)
    uti_plot1d(arCutE, range_e, label1E)
    uti_plot1d(arCutX, range_x, label1X)
    uti_plot1d(arCutY, range_y, label1Y)
  return fig

#def srw_ascii_plot(fname):
def uti_data_file_plot(_fname, _read_labels=1, _e=0, _x=0, _y=0, _graphs_joined=1):
  #data, mode, allrange = srw_ascii_load(fname)
  data, mode, allrange, arLabels, arUnits = uti_data_file_load(_fname, _read_labels)
  #print(allrange)
  m = enum('T','V','H','E','HV','EV','EH','EHV')
  if mode==m.T:
    fig = __mode_T__(data, allrange, arLabels, arUnits, _e, _x, _y)
  elif mode==m.V:
    fig = __mode_V__(data, allrange)
  elif mode==m.H:
    fig = __mode_H__(data, allrange)
  elif mode==m.E:
    fig = __mode_E__(data, allrange, arLabels, arUnits)
  elif mode==m.HV:
    fig = __mode_HV__(data, allrange, arLabels, arUnits, _x, _y, _graphs_joined)
  #elif mode==m.EV:
  #  fig = __mode_EV__(data,allrange)
  #elif mode==m.EH:
  #  fig = __mode_EH__(data,allrange)
  elif mode==m.EHV:
    fig = __mode_EHV__(data, allrange, arLabels, arUnits, _e, _x, _y, _graphs_joined)
  return fig

def uti_plot_show():
    #required for matplotlib
    pl.show()

#****************************************************************************
def uti_interp_1d(_x, _x_min, _x_step, _nx, _ar_f, _ord=3, _ix_per=1, _ix_ofst=0):
    """
    Interpolate 1D function value tabulated on equidistant mesh, using polynomial interpolation
    :param _x: argument at which function value should be calculated
    :param _x_min: minimal argument value of the tabulated function
    :param _x_step: step of mesh at which function is tabulated
    :param _nx: number of points in mesh at which function is tabulated
    :param _ar_f: tabulated function list or array
    :param _ord: order of polynomial interpolation (1- linear, 2- quadratic, 3- cubic)
    :param _ix_per: argument index period of function data alignment (e.g. to interpolate one component of complex data, or in one dimension of multi-dimensional data)
    :param _ix_ofst: argument index offset of function data alignment
    :return: function value found by polynomial interpolation
    """
    if(_ord == 1):
        i0 = int(trunc((_x - _x_min)/_x_step + 1.e-09))
        if(i0 < 0):
            i0 = 0
        elif(i0 >= _nx - 1):
            i0 = _nx - 2
        i1 = i0 + 1
        f0 = _ar_f[i0*_ix_per + _ix_ofst]
        f1 = _ar_f[i1*_ix_per + _ix_ofst]
        t = (_x - (_x_min + _x_step*i0))/_x_step
        return f0 + (f1 - f0)*t
    elif(_ord == 2):
        i0 = int(round((_x - _x_min)/_x_step))
        if(i0 < 1):
            i0 = 1
        elif(i0 >= _nx - 1):
            i0 = _nx - 2
        im1 = i0 - 1
        i1 = i0 + 1
        t = (_x - (_x_min + _x_step*i0))/_x_step
        a0 = _ar_f[i0*_ix_per + _ix_ofst]
        fm1 = _ar_f[im1*_ix_per + _ix_ofst]
        f1 = _ar_f[i1*_ix_per + _ix_ofst]
        a1 = 0.5*(f1 - fm1)
        a2 = 0.5*(fm1 + f1 - 2*a0)
        return a0 + t*(a1 + t*a2)
    elif(_ord == 3):
        i0 = int(trunc((_x - _x_min)/_x_step + 1.e-09))
        if(i0 < 1):
            i0 = 1
        elif(i0 >= _nx - 2):
            i0 = _nx - 3
        im1 = i0 - 1
        i1 = i0 + 1
        i2 = i0 + 2
        t = (_x - (_x_min + _x_step*i0))/_x_step
        a0 = _ar_f[i0*_ix_per + _ix_ofst]
        fm1 = _ar_f[im1*_ix_per + _ix_ofst]
        f1 = _ar_f[i1*_ix_per + _ix_ofst]
        f2 = _ar_f[i2*_ix_per + _ix_ofst]
        a1 = -0.5*a0 + f1 - f2/6. - fm1/3.
        a2 = -a0 + 0.5*(f1 + fm1)
        a3 = 0.5*(a0 - f1) + (f2 - fm1)/6.
        return a0 + t*(a1 + t*(a2 + t*a3))
    return 0

#****************************************************************************
def uti_interp_2d(_x, _y, _x_min, _x_step, _nx, _y_min, _y_step, _ny, _ar_f, _ord=3, _ix_per=1, _ix_ofst=0):
    """
    Interpolate 2D function value tabulated on equidistant rectangular mesh and represented by C-aligned flat array, using polynomial interpolation
    :param _x: first argument at which function value should be calculated
    :param _y: second argument at which function value should be calculated
    :param _x_min: minimal value of the first argument of the tabulated function
    :param _x_step: step of the first argument at which function is tabulated
    :param _nx: number of points vs first argument at which function is tabulated
    :param _y_min: minimal value of the second argument of the tabulated function
    :param _y_step: step of the second argument at which function is tabulated
    :param _ny: number of points vs second argument at which function is tabulated
    :param _ar_f: function tabulated on 2D mesh, aligned as "flat" C-type list or array (first argument is changing most frequently)
    :param _ord: "order" of polynomial interpolation (1- bi-linear (on 4 points), 2- "bi-quadratic" (on 6 points), 3- "bi-cubic" (on 12 points))
    :param _ix_per: period of first argument index of the function data alignment (e.g. to interpolate one component of complex data, or in one dimension of multi-dimensional data)
    :param _ix_ofst: offset of the first argument index in function data alignment
    :return: function value found by 2D polynomial interpolation
    """
    if(_ord == 1): #bi-linear interpolation based on 4 points
        ix0 = int(trunc((_x - _x_min)/_x_step + 1.e-09))
        if(ix0 < 0):
            ix0 = 0
        elif(ix0 >= _nx - 1):
            ix0 = _nx - 2
        ix1 = ix0 + 1
        tx = (_x - (_x_min + _x_step*ix0))/_x_step
        
        iy0 = int(trunc((_y - _y_min)/_y_step + 1.e-09))
        if(iy0 < 0):
            iy0 = 0
        elif(iy0 >= _ny - 1):
            iy0 = _ny - 2
        iy1 = iy0 + 1
        ty = (_y - (_y_min + _y_step*iy0))/_y_step

        nx_ix_per = _nx*_ix_per
        iy0_nx_ix_per = iy0*nx_ix_per
        iy1_nx_ix_per = iy1*nx_ix_per
        ix0_ix_per_p_ix_ofst = ix0*_ix_per + _ix_ofst
        ix1_ix_per_p_ix_ofst = ix1*_ix_per + _ix_ofst
        a00 = _ar_f[iy0_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f10 = _ar_f[iy0_nx_ix_per + ix1_ix_per_p_ix_ofst]
        f01 = _ar_f[iy1_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f11 = _ar_f[iy1_nx_ix_per + ix1_ix_per_p_ix_ofst]
        a10 = f10 - a00
        a01 = f01 - a00
        a11 = a00 - f01 - f10 + f11
        return a00 + tx*(a10 + ty*a11) + ty*a01

    elif(_ord == 2): #bi-quadratic interpolation based on 6 points
        ix0 = int(round((_x - _x_min)/_x_step))
        if(ix0 < 1):
            ix0 = 1
        elif(ix0 >= _nx - 1):
            ix0 = _nx - 2
        ixm1 = ix0 - 1
        ix1 = ix0 + 1
        tx = (_x - (_x_min + _x_step*ix0))/_x_step

        iy0 = int(round((_y - _y_min)/_y_step))
        if(iy0 < 1):
            iy0 = 1
        elif(iy0 >= _ny - 1):
            iy0 = _ny - 2
        iym1 = iy0 - 1
        iy1 = iy0 + 1
        ty = (_y - (_y_min + _y_step*iy0))/_y_step

        nx_ix_per = _nx*_ix_per
        iym1_nx_ix_per = iym1*nx_ix_per
        iy0_nx_ix_per = iy0*nx_ix_per
        iy1_nx_ix_per = iy1*nx_ix_per
        ixm1_ix_per_p_ix_ofst = ixm1*_ix_per + _ix_ofst
        ix0_ix_per_p_ix_ofst = ix0*_ix_per + _ix_ofst
        ix1_ix_per_p_ix_ofst = ix1*_ix_per + _ix_ofst
        fm10 = _ar_f[iy0_nx_ix_per + ixm1_ix_per_p_ix_ofst]
        a00 = _ar_f[iy0_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f10 = _ar_f[iy0_nx_ix_per + ix1_ix_per_p_ix_ofst]
        f0m1 = _ar_f[iym1_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f01 = _ar_f[iy1_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f11 = _ar_f[iy1_nx_ix_per + ix1_ix_per_p_ix_ofst]
        a10 = 0.5*(f10 - fm10)
        a01 = 0.5*(f01 - f0m1)
        a11 = a00 - f01 - f10 + f11
        a20 = 0.5*(f10 + fm10) - a00
        a02 = 0.5*(f01 + f0m1) - a00
        return a00 + tx*(a10 + tx*a20 + ty*a11) + ty*(a01 + ty*a02)
    
    elif(_ord == 3): #bi-cubic interpolation based on 12 points
        ix0 = int(trunc((_x - _x_min)/_x_step + 1.e-09))
        if(ix0 < 1):
            ix0 = 1
        elif(ix0 >= _nx - 2):
            ix0 = _nx - 3
        ixm1 = ix0 - 1
        ix1 = ix0 + 1
        ix2 = ix0 + 2
        tx = (_x - (_x_min + _x_step*ix0))/_x_step

        iy0 = int(trunc((_y - _y_min)/_y_step + 1.e-09))
        if(iy0 < 1):
            iy0 = 1
        elif(iy0 >= _ny - 2):
            iy0 = _ny - 3
        iym1 = iy0 - 1
        iy1 = iy0 + 1
        iy2 = iy0 + 2
        ty = (_y - (_y_min + _y_step*iy0))/_y_step

        nx_ix_per = _nx*_ix_per
        iym1_nx_ix_per = iym1*nx_ix_per
        iy0_nx_ix_per = iy0*nx_ix_per
        iy1_nx_ix_per = iy1*nx_ix_per
        iy2_nx_ix_per = iy2*nx_ix_per
        ixm1_ix_per_p_ix_ofst = ixm1*_ix_per + _ix_ofst
        ix0_ix_per_p_ix_ofst = ix0*_ix_per + _ix_ofst
        ix1_ix_per_p_ix_ofst = ix1*_ix_per + _ix_ofst
        ix2_ix_per_p_ix_ofst = ix2*_ix_per + _ix_ofst
        f0m1 = _ar_f[iym1_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f1m1 = _ar_f[iym1_nx_ix_per + ix1_ix_per_p_ix_ofst]
        fm10 = _ar_f[iy0_nx_ix_per + ixm1_ix_per_p_ix_ofst]
        a00 = _ar_f[iy0_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f10 = _ar_f[iy0_nx_ix_per + ix1_ix_per_p_ix_ofst]
        f20 = _ar_f[iy0_nx_ix_per + ix2_ix_per_p_ix_ofst]
        fm11 = _ar_f[iy1_nx_ix_per + ixm1_ix_per_p_ix_ofst]
        f01 = _ar_f[iy1_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f11 = _ar_f[iy1_nx_ix_per + ix1_ix_per_p_ix_ofst]
        f21 = _ar_f[iy1_nx_ix_per + ix2_ix_per_p_ix_ofst]
        f02 = _ar_f[iy2_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f12 = _ar_f[iy2_nx_ix_per + ix1_ix_per_p_ix_ofst]
        a10 = -0.5*a00 + f10 - f20/6 - fm10/3
        a01 = -0.5*a00 + f01 - f02/6 - f0m1/3
        a11 = -0.5*(f01 + f10) + (f02 - f12 + f20 - f21)/6 + (f0m1 - f1m1 + fm10 - fm11)/3 + f11
        a20 = -a00 + 0.5*(f10 + fm10)
        a02 = -a00 + 0.5*(f01 + f0m1)
        a21 = a00 - f01 + 0.5*(f11 - f10 - fm10 + fm11)
        a12 = a00 - f10 + 0.5*(f11 - f01 - f0m1 + f1m1)
        a30 = 0.5*(a00 - f10) + (f20 - fm10)/6
        a03 = 0.5*(a00 - f01) + (f02 - f0m1)/6
        a31 = 0.5*(f01 + f10 - f11 - a00) + (f21 + fm10 - f20 - fm11)/6
        a13 = 0.5*(f10 - f11 - a00 + f01) + (f0m1 + f12 - f02 - f1m1)/6
        return a00 + tx*(a10 + tx*(a20 + tx*(a30 + ty*a31) + ty*a21) + ty*a11) + ty*(a01 + ty*(a02 + ty*(a03 + tx*a13) + tx*a12))
    return 0
