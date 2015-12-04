"""uti_plot_com module containing plot utilities not specific to a particular backend

.. moduleauthor:: 
"""
from copy import *
from array import *
import traceback
#import sys
#import numpy as np

import uti_math
import uti_io

#****************************************************************************
#def srw_ascii_load(fname):
def file_load(_fname, _read_labels=1):
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
    else: nLinesHead -= 1

    e0,e1,x0,x1,y0,y1 = [float(hlp[i].replace('#','').split()[0]) for i in [1,2,4,5,7,8]]

    #data = np.squeeze(np.loadtxt(_fname, dtype=np.float64)) #get data from file (C-aligned flat)
    data = uti_io.read_ascii_data_cols(_fname, '\t', _i_col_start=0, _i_col_end=0, _n_line_skip=nLinesHead)[0]

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

    m = _enum('T','V','H','E','HV','EV','EH','EHV')

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

#****************************************************************************
def rescale(maxabsval, strval):
    """Force labels to 1.0e-3 boundary which contains maxabsval
    :param double maxabsval: absolute value on axis
    :param str strval: units
    :return (multiplier, strval): axis multiplier, axis label
    """
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

#****************************************************************************
def rescale_range(allrange, _ar_units, _ec=0, _xc=0, _yc=0):
    """Adjust graph axis ranges and labels to be 1.0e-3 boundary

    :param tuple allrange: Order of ranges: e, x, y
    :param tuple _ar_units: units for ranges [e, x, y]
    :param 
    """
    e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange

    #em = abs(np.array([e0,e1])).max()
    #xm = abs(np.array([x0,x1])).max()
    #ym = abs(np.array([y0,y1])).max()

    abs_e0 = abs(e0); abs_e1 = abs(e1)
    em = abs_e0
    if(em < abs_e1): em = abs_e1

    abs_x0 = abs(x0); abs_x1 = abs(x1)
    xm = abs_x0
    if(xm < abs_x1): xm = abs_x1

    abs_y0 = abs(y0); abs_y1 = abs(y1)
    ym = abs_y0
    if(ym < abs_y1): ym = abs_y1
    
    #mult_e, str_e = _rescale(em,"eV")
    #mult_x, str_x = _rescale(xm,"m")
    #mult_y, str_y = _rescale(ym,"m")
    #mult_e, str_e = _rescale(em, _ar_units[0])
    #mult_x, str_x = _rescale(xm, _ar_units[1])
    #mult_y, str_y = _rescale(ym, _ar_units[2])
    mult_e, str_e = rescale(em, _ar_units[0])
    mult_x, str_x = rescale(xm, _ar_units[1])
    mult_y, str_y = rescale(ym, _ar_units[2])

    e0s = uti_math.num_round(e0*mult_e) #this is done to avoid in the labels numbers like "12.700000000001"
    e1s = uti_math.num_round(e1*mult_e)
    x0s = uti_math.num_round(x0*mult_x)
    x1s = uti_math.num_round(x1*mult_x)
    y0s = uti_math.num_round(y0*mult_y)
    y1s = uti_math.num_round(y1*mult_y)
    ecs = uti_math.num_round(_ec*mult_e)
    xcs = uti_math.num_round(_xc*mult_x)
    ycs = uti_math.num_round(_yc*mult_y)
    
    #allnewrange = e0*mult_e, e1*mult_e, ne, x0*mult_x, x1*mult_x, nx, y0*mult_y, y1*mult_y, ny, _ec*mult_e, _xc*mult_x, _yc*mult_y
    allnewrange = e0s, e1s, ne, x0s, x1s, nx, y0s, y1s, ny, ecs, xcs, ycs
    strval = str_e, str_x, str_y
    return allnewrange, strval
    
#****************************************************************************
def rescale_dim(_range, _base_unit=None):
    """Adjust range and units of a value ("dimension" of a plot) to be 1.0e-3 boundary

    :param list _range: min. and max. value of a range to be adjusted
    :param sting _base_unit: base unit (e.g. [m], [eV],...)
    :return: tuple containing new adjusted range and unit
    """

    #xm = abs(np.array([_range[0], _range[1]])).max()
    abs_x0 = abs(_range[0]); abs_x1 = abs(_range[1])
    xm = abs_x0
    if(xm < abs_x1): xm = abs_x1
    
    mult, unit = rescale(xm, _base_unit)
    newrange = deepcopy(_range)
    newrange[0] *= mult; newrange[1] *= mult
    return newrange, unit
 
#****************************************************************************
def _enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

