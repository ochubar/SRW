"""uti_plot backend for matplotlib

.. moduleauthor:: O. Chubar, N. Canestrari
"""
import os
import sys
import subprocess
import platform
from array import *
from math import *
import numpy as np
import uti_plot_com
import uti_math
import uti_plot

class Backend(object):
    def __init__(self, backend, fname_format):
        """Set matplotlib backend depending on IPython state, mpld3 availability and fname_format
        Tries `backend` based on what the execution context is (IPython, X11, etc.).
        Will throw exception if matplotlib or pyplot cannot be imported.

        :param str backend: attempted backend
        :param str fname_format: where to save figs
        """
        import matplotlib
        if self._running_in_ipython():
            backend = self._init_ipython(backend)
            import matplotlib.pyplot
            self._pl = matplotlib.pyplot
        else:
            if backend == uti_plot.DEFAULT_BACKEND:
                (backend, fname_format) = self._default_backend(fname_format)
            matplotlib.use(backend)
            import matplotlib.pyplot
            self._pl = matplotlib.pyplot
            matplotlib.pyplot.ioff()
            (backend, fname_format) = self._verify_pyplot(backend, fname_format)
        if fname_format is not None:
            print('Saving figures to ' + fname_format)
        self._backend = backend
        self._fname_format = fname_format
        self._figure_num = 0
  
    def uti_plot_show(self):
        'Render the plots generated thus far; Prompts user to do the generation'
        self._pyplot_show()

    #def plot1D(ar1d,x_range,labels=('energy, KeV','Ph/s/0.1%BW')):
    def uti_plot1d(self, ar1d, x_range, labels=('energy [eV]', 'ph/s/0.1%bw')):
        #fig = _pl.figure(figsize=(12,8))
        fig = self._pl.figure()
        self._plot_1D(ar1d,x_range,labels,fig)
        return self._maybe_savefig(fig)

    #def plot2D(ar2d,x_range,y_range,labels=('Horizontal position, mm','Vertical position, mm')):
    def uti_plot2d(self, ar2d, x_range, y_range, labels=('Horizontal Position [m]','Vertical Position [m]')):
        #fig = _pl.figure(figsize=(12,8))
        fig = self._pl.figure()
        self._plot_2D(ar2d,x_range,y_range,labels,fig)
        return self._maybe_savefig(fig)

    def uti_plot2d1d(self, data, x_range, y_range, xc, yc, labels, _graphs_joined=True):
        x0 = x_range[0]; x1 = x_range[1]; nx = x_range[2]
        #y0 = x_range[0]; y1 = y_range[1]; ny = y_range[2]
        y0 = y_range[0]; y1 = y_range[1]; ny = y_range[2] #OC090714

        label2D = labels[0]; label1H = labels[1]; label1V = labels[2]

        fig = None
        if _graphs_joined:
            #fig = _pl.figure(figsize=(12,5))
            fig = self._pl.figure(figsize=(15,5.3))
           
            self._plot_2D(data, x_range, y_range, label2D, fig, 131) #showing graphs in one panel
        else: self.uti_plot2d(data, x_range, y_range, label2D)

        xStep = 0
        if nx > 1: xStep = (x1 - x0)/(nx - 1)
        yStep = 0
        if ny > 1: yStep = (y1 - y0)/(ny - 1)
        inperpOrd = 1 #interpolation order to use (1 to 3)

        #_plot_1D(data[iy,:],range_x,label1H,fig,132)
        arCutX = array('d', [0]*nx)
        xx = x0
        for ix in range(nx):
            arCutX[ix] = uti_math.interp_2d(xx, yc, x0, xStep, nx, y0, yStep, ny, data, inperpOrd, 1, 0)
            xx += xStep
        if _graphs_joined: self._plot_1D(arCutX, x_range, label1H, fig, 132) #OC150814
        else: self.uti_plot1d(arCutX, x_range, label1H)

        #_plot_1D(data[:,ix],range_y,label1V,fig,133)
        arCutY = array('d', [0]*ny)
        yy = y0
        for iy in range(ny):
            #arCutY[iy] = _interp_2d(xc, yy, x0, xStep, nx, y0, yStep, ny, data, inperpOrd, 1, 0)
            arCutY[iy] = uti_math.interp_2d(xc, yy, x0, xStep, nx, y0, yStep, ny, data, inperpOrd, 1, 0)
            yy += yStep
        if _graphs_joined: self._plot_1D(arCutY, y_range, label1V, fig, 133)
        else: self.uti_plot1d(arCutY, y_range, label1V)

        if _graphs_joined: self._pl.tight_layout() #OC081115

        return self._maybe_savefig(fig)

    #def srw_ascii_plot(fname):
    def uti_data_file_plot(self, _fname, _read_labels=1, _e=0, _x=0, _y=0, _graphs_joined=1):
        #data, mode, allrange = srw_ascii_load(fname)
        #data, mode, allrange, arLabels, arUnits = _file_load(_fname, _read_labels)
        data, mode, allrange, arLabels, arUnits = uti_plot_com.file_load(_fname, _read_labels)

        #print(allrange)
        m = self._enum('T','V','H','E','HV','EV','EH','EHV')
        if mode==m.T:
            fig = self.__mode_T(data, allrange, arLabels, arUnits, _e, _x, _y)
        elif mode==m.V:
            fig = self.__mode_V(data, allrange)
        elif mode==m.H:
            fig = self.__mode_H(data, allrange)
        elif mode==m.E:
            fig = self.__mode_E(data, allrange, arLabels, arUnits)
        elif mode==m.HV:
            fig = self.__mode_HV(data, allrange, arLabels, arUnits, _x, _y, _graphs_joined)
        #elif mode==m.EV:
        #  fig = __mode_EV(data,allrange)
        #elif mode==m.EH:
        #  fig = __mode_EH(data,allrange)
        elif mode==m.EHV:
            fig = self.__mode_EHV(data, allrange, arLabels, arUnits, _e, _x, _y, _graphs_joined)
        return self._maybe_savefig(fig)

    def _enum(self, *sequential, **named): #Had to copy this in uti_plot_com
        enums = dict(zip(sequential, range(len(sequential))), **named)
        return type('Enum', (), enums)

    def _plot_1D(self, ar1d, x_range, labels, fig, typ=111):
        lenAr1d = len(ar1d)
        if lenAr1d > x_range[2]: ar1d = np.array(ar1d[0:x_range[2]])
        elif lenAr1d < x_range[2]:
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

    def _plot_2D(self, ar2d, x_range, y_range, labels, fig, typ=111):
        totLen = int(x_range[2]*y_range[2])
        lenAr2d = len(ar2d)
        if lenAr2d > totLen: ar2d = np.array(ar2d[0:totLen])
        elif lenAr2d < totLen:
            auxAr = array('d', [0]*lenAr2d)
            for i in range(lenAr2d): auxAr[i] = ar2d[i]
            ar2d = np.array(auxAr)

        #if isinstance(ar2d,(list,array)): ar2d = np.array(ar2d).reshape(x_range[2],y_range[2])
        if isinstance(ar2d,(list,array)): ar2d = np.array(ar2d)
        ar2d = ar2d.reshape(y_range[2],x_range[2])

        x = np.linspace(x_range[0],x_range[1],x_range[2])
        y = np.linspace(y_range[0],y_range[1],y_range[2])
        ax = fig.add_subplot(typ)
        #ax.pcolormesh(x,y,ar2d.T,cmap=_pl.cm.Greys_r)
        #ax.pcolormesh(x,y,ar2d,cmap=_pl.cm.Greys_r)
        ax.pcolormesh(x,y,ar2d,cmap=self._pl.cm.Greys_r) #OC150814
        
        ax.set_xlim(x[0],x[-1])
        ax.set_ylim(y[0],y[-1])
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        if(len(labels) > 2):
            ax.set_title(labels[2])

    def __mode_T(self, data, allrange, _ar_labels, _ar_units, _ec=0, _xc=0, _yc=0):
        #allrange, units = _rescale_range(allrange)
        #allrange, units = _rescale_range(allrange, _ar_units, _ec, _xc, _yc)
        allrange, units = uti_plot_com.rescale_range(allrange, _ar_units, _ec, _xc, _yc)

        #e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
        e0, e1, ne, x0, x1, nx, y0, y1, ny, ec, xc, yc = allrange
        #toprint = (e0,e1,units[0], x0,x1,units[1], y0,y1,units[2], data[0]) #squeeze have reduced it to an array with one element.
        toprint = (e0,units[0], x0,units[1], y0,units[2], data[0],units[3]) #squeeze have reduced it to an array with one element.
        #sys.stdout.write('Total Flux for \nE: %f -> %f %s\nX: %f -> %f %s\nY: %f -> %f %s\n is %f Ph/s/0.1%BW' % toprint) 
        sys.stdout.write(_ar_labels[3] + ' for \n' + _ar_labels[0] + ': %f %s\n' + _ar_labels[1] + ': %f %s\n' + _ar_labels[2] + ': %f %s\n is %f %s' % toprint) 
        return None

    def __mode_V(self, data, allrange, _ar_labels, _ar_units):
        #allrange, units = _rescale_range(allrange)
        #allrange, units = _rescale_range(allrange, _ar_units)
        allrange, units = uti_plot_com.rescale_range(allrange, _ar_units)

        #e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
        e0, e1, ne, x0, x1, nx, y0, y1, ny, ec, xc, yc = allrange
        range_y = y0, y1, ny
        #label = ("Vertical Position, ["+units[2]+"]","Ph/s/0.1%BW/mm^2")
        label = (_ar_labels[2] + ' [' + units[2] + ']', _ar_labels[3] + ' [' + _ar_units[3] + ']')
        fig = self._pl.figure(figsize=(4,4))
        self._plot_1D(data, range_y, label, fig)
        return fig

    def __mode_H(self, data, allrange, _ar_labels, _ar_units):
        #allrange, units = _rescale_range(allrange)
        #allrange, units = _rescale_range(allrange, _ar_units)
        allrange, units = uti_plot_com.rescale_range(allrange, _ar_units)

        #e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
        e0, e1, ne, x0, x1, nx, y0, y1, ny, ec, xc, yc = allrange
        range_x = x0, x1, nx 
        #label = ("Horizontal Position, ["+units[1]+"]","Ph/s/0.1%BW/mm^2")
        label = (_ar_labels[1] + ' [' + units[1] + ']', _ar_labels[3] + ' [' + _ar_units[3] + ']')
        fig = self._pl.figure(figsize=(4,4))
        self._plot_1D(data, range_x, label, fig)
        return fig

    def __mode_E(self, data, allrange, _ar_labels, _ar_units):
        #allrange, units = _rescale_range(allrange)
        #allrange, units = _rescale_range(allrange, _ar_units)
        allrange, units = uti_plot_com.rescale_range(allrange, _ar_units)

        e0, e1, ne, x0, x1, nx, y0, y1, ny, ec, xc, yc = allrange
        range_e = e0, e1, ne
        #label = ("Energy, ["+units[0]+"]","Ph/s/0.1%BW")
        label = (_ar_labels[0] + ' [' + units[0] + ']', _ar_labels[3] + ' [' + _ar_units[3] + ']')
        fig = self._pl.figure(figsize=(4,4))
        self._plot_1D(data,range_e,label,fig)
        return fig

    def __mode_HV(self, data, allrange, _ar_labels, _ar_units, _xc=0, _yc=0, _graphs_joined=True):
        #Could be moved to uti_plot_com.py (since there is not Matplotlib-specific content)
        #allrange, units = _rescale_range(allrange)
        #allrange, units = _rescale_range(allrange, _ar_units, 0, _xc, _yc)
        allrange, units = uti_plot_com.rescale_range(allrange, _ar_units, 0, _xc, _yc)

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
        #print(label2D)

        #label1H = ("Horizontal Position, ["+units[1]+"]","Ph/s/0.1%BW/mm^2")
        strTitle = 'At ' + _ar_labels[2] + ': ' + str(yc)
        if yc != 0: strTitle += ' ' + units[2]
        label1H = (_ar_labels[1] + ' [' + units[1] + ']', _ar_labels[3] + ' [' + _ar_units[3] + ']', strTitle)

        #label1V = ("Vertical Position, ["+units[2]+"]","Ph/s/0.1%BW/mm^2")
        strTitle = 'At ' + _ar_labels[1] + ': ' + str(xc)
        if xc != 0: strTitle += ' ' + units[1]
        label1V = (_ar_labels[2] + ' [' + units[2] + ']', _ar_labels[3] + ' [' + _ar_units[3] + ']', strTitle)

        #return plot_2D_1D(data, range_x, range_y, [label2D, label1H, label1V], _graphs_joined)
        return self.uti_plot2d1d(data, range_x, range_y, xc, yc, [label2D, label1H, label1V], _graphs_joined)

    ##    fig = None
    ##    if _graphs_joined:
    ##        fig = _pl.figure(figsize=(12,5))
    ##        _plot_2D(data, range_x, range_y, label2D, fig, 131) #showing graphs in one figure
    ##    else: uti_plot2d(data, range_x, range_y, label2D)
    ##
    ##    xStep = 0
    ##    if nx > 1: xStep = (x1 - x0)/(nx - 1)
    ##    yStep = 0
    ##    if ny > 1: yStep = (y1 - y0)/(ny - 1)
    ##    inperpOrd = 1 #interpolation order (1 to 3)
    ##
    ##    #_plot_1D(data[iy,:],range_x,label1H,fig,132)
    ##    arCutX = array('d', [0]*nx)
    ##    xx = x0
    ##    for ix in range(nx):
    ##        #arCutX[ix] = _interp_2d(xx, yc, x0, xStep, nx, y0, yStep, ny, data, inperpOrd, 1, 0)
    ##        arCutX[ix] = uti_math.interp_2d(xx, yc, x0, xStep, nx, y0, yStep, ny, data, inperpOrd, 1, 0)
    ##        xx += xStep
    ##    if _graphs_joined: _plot_1D(arCutX, range_x, label1H, fig, 132)
    ##    else: uti_plot1d(arCutX, range_x, label1H)
    ##
    ##    #_plot_1D(data[:,ix],range_y,label1V,fig,133)
    ##    arCutY = array('d', [0]*ny)
    ##    yy = y0
    ##    for iy in range(ny):
    ##        #arCutY[iy] = _interp_2d(xc, yy, x0, xStep, nx, y0, yStep, ny, data, inperpOrd, 1, 0)
    ##        arCutY[iy] = uti_math.interp_2d(xc, yy, x0, xStep, nx, y0, yStep, ny, data, inperpOrd, 1, 0)
    ##        yy += yStep
    ##    if _graphs_joined: _plot_1D(arCutY, range_y, label1V, fig, 133)
    ##    else: uti_plot1d(arCutY, range_y, label1V)
    ##
    ##    return fig

    #def __mode_EV(data, allrange):
    #  #to be updated
    #  allrange, units = _rescale_range(allrange)
    #  e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
    #  range_e = e0, e1, ne  
    #  range_y = y0, y1, ny
    #  e = np.linspace(e0, e1, ne)
    #  y = np.linspace(y0, y1, ny)
    #  ie = np.where(data.sum(axis=1)==data.sum(axis=1).max())[0][0]
    #  iy = np.where(abs(y).min())[0][0]
    #  label1E = ("Energy, ["+units[0]+"]","Ph/s/0.1%BW/mm^2")
    #  label1V = ("Vertical Position, ["+units[2]+"]","Ph/s/0.1%BW/mm^2")
    #  fig = _pl.figure(figsize=(8,4))
    #  _plot_1D(data[iy,:],range_e,label1E,fig,121)
    #  _plot_1D(data[:,ie],range_y,label1V,fig,122)
    #  return fig

    #def __mode_EH(data, allrange):
    #  #to be updated
    #  allrange, units = _rescale_range(allrange)
    #  e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
    #  range_e = e0, e1, ne  
    #  range_x = x0, x1, nx
    #  e = np.linspace(e0, e1, ne)
    #  x = np.linspace(x0, x1, nx)
    #  ie = np.where(data.sum(axis=1)==data.sum(axis=1).max())[0][0]
    #  ix = np.where(abs(y).min())[0][0]
    #  label1E = ("Energy, ["+units[0]+"]","Ph/s/0.1%BW/mm^2")
    #  label1H = ("Horizontal Position, ["+units[1]+"]","Ph/s/0.1%BW/mm^2")
    #  fig = _pl.figure(figsize=(8,4))
    #  _plot_1D(data[ix,:],range_e,label1E,fig,121)
    #  _plot_1D(data[:,ie],range_x,label1H,fig,122)
    #  return fig

    def __mode_EHV(self, data, allrange, _ar_labels, _ar_units, _ec, _xc, _yc, _graphs_joined=1):
        #allrange, units = _rescale_range(allrange)
        #allrange, units = _rescale_range(allrange, _ar_units, _ec, _xc, _yc)
        allrange, units = uti_plot_com.rescale_range(allrange, _ar_units, _ec, _xc, _yc)

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

        #fig = _pl.figure(figsize=(8,8))
        #_plot_2D(data[:,:,ie], range_x, range_y, label2D, fig, 221)
        #_plot_1D(data[ix,iy,:],range_e,label1E,fig,224)
        #_plot_1D(data[ie,:,iy],range_x,label1H,fig,222)
        #_plot_1D(data[:,ix,ie],range_y,label1V,fig,223)

        fig = None
        if _graphs_joined:
            fig = self._pl.figure(figsize=(12,5))
            self._plot_2D(arCutXY, range_x, range_y, label2D, fig, 221) #showing graphs in one figure
            self._plot_1D(arCutE, range_e, label1E, fig, 222)
            self._plot_1D(arCutX, range_x, label1X, fig, 223)
            self._plot_1D(arCutY, range_y, label1Y, fig, 224)
        else:
            self.uti_plot2d(arCutXY, range_x, range_y, label2D)
            self.uti_plot1d(arCutE, range_e, label1E)
            self.uti_plot1d(arCutX, range_x, label1X)
            self.uti_plot1d(arCutY, range_y, label1Y)
        return self._maybe_savefig(fig)

    def _maybe_savefig(self, figure):
        """If there is an ``fname_format``, then save the pyplot `figure`
        (if not None) to the computed file name
        """
        if figure is None:
            return figure
        if self._fname_format is None:
            return figure
        figure.savefig(self._fname_format.format(self._figure_num))
        self._figure_num += 1
        return figure

    def _pyplot_show(self):
        """'Render the plots, if the backend is interactive.
        Plots saved to files (``fname_format``), have already been written.
        """
        with self._HideGUIErrorOutput():
            self._pl.show()

    def _default_fname_format(self):
        """Use main program's name (or uti_plot if none) for default file format name.
        The default file type is png
        """
        try:
            base = os.path.splitext(
                os.path.basename(sys.modules['__main__']['__file__']))
        except:
            base = 'uti_plot'
        fname_format = base + '-{}.png'
        return fname_format

    def _default_backend(self, fname_format):
        """Selects backend based on execution context (X11, macosx, etc.).  If an
        interactive backend is not possible, return _default_file_backend()
        """
        if self._running_in_x11() or self._running_in_windows():
            backend = 'TkAgg'
        elif self._running_in_macosx():
            backend = 'macosx'
        else:
            (backend, fname_format) = self._default_file_backend(fname_format)
        return (backend, fname_format)

    def _default_file_backend(self, fname_format):
        """Returns "agg" for backend with fname_format.  If fname_format is none,
        calls _default_fname_format() to compute fname_format.
        """
        if fname_format is None:
            fname_format = self._default_fname_format()
        return ('agg', fname_format)

    def _init_ipython(self, backend):
        """Tries `backend`.  Returns ``inline`` backend if default requested or
        unable to load requested backend.
        """
        if backend == uti_plot.DEFAULT_BACKEND:
            backend = 'inline'
        is_mpld3 = backend == 'mpld3'
        b = 'inline' if is_mpld3 else backend
        get_ipython().magic('matplotlib ' + b)
        if is_mpld3:
            try:
                import mpld3
                mpld3.enable_notebook()
            except:
                print('mpld3: backend unavailable; plots will be inline')
                backend = 'inline'
        return backend

    def _verify_pyplot(self, backend, fname_format):
        """Create a pyplot figure and close it to be sure backend functions properly.
        If backend fails, configures _default_file_backend().
        """
        import matplotlib.pyplot
        pl = matplotlib.pyplot
        try:
            pl.figure(figsize=(0,0))
            pl.close('all')
        except:
            old = backend
            (backend, fname_format) = self._default_file_backend(fname_format)
            pl.switch_backend(backend)
            # If this raises an exception, uti_plot will set backend to None 
            pl.figure()
            pl.close('all')
            print(old + ': backend unavailable; plots will be saved to files')
        return (backend, fname_format)

    def _running_in_ipython(self):
        'Is the current interpreter IPython?'
        try:
            __IPYTHON__
            return True
        except NameError:
            return False

    def _running_in_macosx(self):
        'Is Mac OS X running?'
        return platform.system() == 'Darwin'

    def _running_in_windows(self):
        'Is Darwin running?'
        return platform.system() == 'Windows'

    def _running_in_x11(self):
        'Is X11 running?'
        d = os.environ.get('DISPLAY')
        if d is None or d == '':
            return False
        devnull = open(os.devnull, 'wb')
        tries = [
            # xset is not in core install on some OSes
            ['xset', '-q'],
            # This takes a couple of seconds
            ['xterm', '/bin/true']]
        for t in tries:
            try:
                subprocess.call(t, stdout=devnull, stderr=devnull)
                return True
            except:
                pass
        return false

    class _HideGUIErrorOutput(object):
        """Redirect low level stderr (fd #2) to os.devnull as context manager

        TkAgg has a bug::

            can't invoke "event" command:  application has been destroyed
                while executing
            "event generate $w <<ThemeChanged>>"
                (procedure "ttk::ThemeChanged" line 6)
                invoked from within
            "ttk::ThemeChanged"

        The output comes from C so it has to be redirected at the file descriptor
        level.  This class does the low level redirection so that message doesn't
        come out.  Example::

            with HideGUIErrorOutput():
                pl.show()
        """
        def __init__(self):
            pass

        def __enter__(self):
            'Redirects stderr to devnull, saving _prev_stderr'
            if sys.stderr is None or sys.__stderr__ is None:
                return
            sys.stderr.flush()
            fno = sys.__stderr__.fileno() #OC150814
            if(fno >= 0):
                null = os.open(os.devnull, os.O_WRONLY)
                #self._prev_stderr = os.dup(sys.__stderr__.fileno())
                self._prev_stderr = os.dup(fno)
                os.dup2(null, sys.__stderr__.fileno())
                os.close(null)

        def __exit__(self, exc_type, exc_value, traceback):
            'Redirects to _prev_stderr'
            if sys.stderr is None or sys.__stderr__ is None:
                return
            sys.__stderr__.flush()
            fno = sys.__stderr__.fileno() #OC150814
            if(fno >= 0):
                #os.dup2(self._prev_stderr, sys.__stderr__.fileno())
                os.dup2(self._prev_stderr, fno)
                os.close(self._prev_stderr)
