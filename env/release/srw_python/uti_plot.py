"""Simple 1D & 2D plotting utilities package for "Synchrotron Radiation Workshop" (SRW).

``uti_plot`` currently wraps ``matplotlib``, but other backends are
planned.  If no suitable backend is available, ``uti_plot_init`` sets
the backend to ``uti_plot_none`` so that the calling program is still
functional.  This is useful for systems where there is no graphing
library available, but you still want to see the results of the
SRW program.

Usage:

    import uti_plot as up

    up.uti_plot_init()
    uti_plot1d(...)
    uti_plot_show()

Modules:

    uti_plot
        This module, which loads all other modules dynamically

    uti_plot_matplotlib
        Does the actually plotting using matplotlib.pyplot.  Currently, loaded in all cases except
        when ``backend`` is ``None``

    test_uti_plot
        Simple tests for uti_plot

.. moduleauthor:: Rob Nagler <nagler@radiasoft.net>
"""
import sys
import uti_plot_com
import traceback

_backend = None

DEFAULT_BACKEND = '<default>'

def uti_plot_init(backend=DEFAULT_BACKEND, fname_format=None):
    """Initializes plotting engine with backend and, optionally, save plots to fname_format

    Tries to initialize `backend` as the plotting engine.  If not found, an
    error will be printed, and this module's functions will be no-ops.  If
    DEFAULT_BACKEND provided, an appropriate backend will be chosen and printed.
    Plots may also be saved if fname_format is supplied.

    You may call ``uti_plot_init(None)`` explicitly so that no plotting occurs.

    :param str backend: a matplot backend (TkAgg, etc.) or ``inline`` in IPython
    :param str fname_format: where to save plots. format field is a sequential plot number, starting at 0.
    """
    global _backend
    if backend is not None:
        try:
            import uti_plot_matplotlib
            _backend = uti_plot_matplotlib.Backend(backend, fname_format)
            return
        except:
            traceback.print_exc()
            print(backend + ': unable to import specified backend (or its dependency); no plots')
    elif fname_format is not None:
        #raise Value(fname_format + ': fname_format must be null if backend is None')
        raise ValueError(fname_format + ': fname_format must be null if backend is None')
    _backend = _BackendNone()

def uti_plot_show():
    """Display the plots"""
    #if '_backend' not in locals(): uti_plot_init() #?
    _backend.uti_plot_show()

def uti_plot1d(ar1d, x_range, labels=('Photon Energy [eV]', 'ph/s/0.1%bw'), units=None):
    """Generate one-dimensional line plot from given array

    :param array ar1d: data points
    :param list x_range: Passed to numpy.linspace(start sequence, stop sequnce, num samples)
    :param tuple labels: [x-axis, y-axis]
    """
    #if '_backend' not in locals(): uti_plot_init() #?
    
    if(units is not None):
        x_range, x_unit = uti_plot_com.rescale_dim(x_range, units[0])
        units = [x_unit, units[1]]
        strTitle = '' if(len(labels) < 3) else labels[2]
        labels = (labels[0] + ' [' + units[0] + ']', labels[1] + ' [' + units[1] + ']', strTitle)

    _backend.uti_plot1d(ar1d, x_range, labels)

def uti_plot1d_ir(ary, arx, labels=('Longitudinal Position [m]', 'Horizontal Position [m]'), units=None): #OC15112017
    """Generate one-dimensional line plot from given array

    :param array arx: abscissa array
    :param array ary: ordinate array
    :param tuple labels: [x-axis, y-axis]
    """
    #if '_backend' not in locals(): uti_plot_init() #?
    
    if(units is not None):
        #x_range = [min(arx), max(arx), len(arx)]
        #x_range, x_unit = uti_plot_com.rescale_dim(x_range, units[0])
        #units = [x_unit, units[1]]
        strTitle = '' if(len(labels) < 3) else labels[2]
        labels = (labels[0] + ' [' + units[0] + ']', labels[1] + ' [' + units[1] + ']', strTitle)

    _backend.uti_plot1d_ir(ary, arx, labels)

def uti_plot2d(ar2d, x_range, y_range, labels=('Horizontal Position [m]','Vertical Position [m]'), units=None):
    """Generate quad mesh plot from given "flattened" array

    :param array ar2d: data points
    :param list x_range: Passed to numpy.linspace(start sequence, stop sequnce, num samples)
    :param list y_range: y axis (same structure as x_range)
    :param tuple labels: [x-axis, y-axis]
    """
    #if '_backend' not in locals(): uti_plot_init() #?
    if(units is not None):
        x_range, x_unit = uti_plot_com.rescale_dim(x_range, units[0])
        y_range, y_unit = uti_plot_com.rescale_dim(y_range, units[1])
        units = [x_unit, y_unit,  units[2]]
        strTitle = '' if(len(labels) < 3) else labels[2]
        labels = (labels[0] + ' [' + units[0]+ ']', labels[1] + ' [' + units[1] + ']', strTitle)

    _backend.uti_plot2d(ar2d, x_range, y_range, labels)

def uti_plot2d1d(ar2d, x_range, y_range, x=0, y=0, labels=('Horizontal Position', 'Vertical Position', 'Intensity'), units=None, graphs_joined=True):
    """Generate 2d quad mesh plot from given "flattened" array, and 1d cuts passing through (x, y)

    :param array ar2d: data points
    :param list x_range: Passed to numpy.linspace(start sequence, stop sequnce, num samples)
    :param list y_range: y axis (same structure as x_range)
    :param x: x value for 1d cut
    :param y: y value for 1d cut
    :param tuple labels: [x-axis, y-axis, z-axis]
    :param tuple units: [x-axis, y-axis, z-axis]
    :param graphs_joined: switch specifying whether the 2d plot and 1d cuts have to be displayed in one panel or separately
    """
    #if '_backend' not in locals(): uti_plot_init() #?
    if(units is not None): #checking / re-scaling x, y
        x_range, x_unit = uti_plot_com.rescale_dim(x_range, units[0])
        y_range, y_unit = uti_plot_com.rescale_dim(y_range, units[1])
        units = [x_unit, y_unit,  units[2]]

        strTitle = labels[2]
        label2D = (labels[0] + ' [' + units[0]+ ']', labels[1] + ' [' + units[1] + ']', strTitle)

        strTitle = 'At ' + labels[1] + ': ' + str(y)
        if y != 0: strTitle += ' ' + units[1]
        label1X = (labels[0] + ' [' + units[0] + ']', labels[2] + ' [' + units[2] + ']', strTitle)

        strTitle = 'At ' + labels[0] + ': ' + str(x)
        if x != 0: strTitle += ' ' + units[0]
        label1Y = (labels[1] + ' [' + units[1] + ']', labels[2] + ' [' + units[2] + ']', strTitle)
        
    else: #OC081115
        strTitle = labels[2]
        label2D = (labels[0], labels[1], strTitle)

        strTitle = 'At ' + labels[1] + ': ' + str(y)
        label1X = (labels[0], labels[2], strTitle)

        strTitle = 'At ' + labels[0] + ': ' + str(x)
        label1Y = (labels[1], labels[2], strTitle)

    labels = [label2D, label1X, label1Y]

    _backend.uti_plot2d1d(ar2d, x_range, y_range, x, y, labels, graphs_joined)

def uti_plot_data_file(_fname, _read_labels=1, _e=0, _x=0, _y=0, _graphs_joined=True, #Same as uti_data_file_plot, but better fits function name decoration rules in this module (uti_plot*)
                       _multicolumn_data=False, _column_x=None, _column_y=None, #MR31102017
                       _scale='linear', _width_pixels=None):
    """Generate plot from configuration in _fname

    :param str _fname: config loaded from here
    :param bool _read_labels: whether to read labels from _fname
    :param float _e: photon energy adjustment
    :param float _x: horizonal position adjustment
    :param float _y: vertical position adjustment
    :param bool _graphs_joined: if true, all plots in a single figure
    :param bool _multicolumn_data: if true, visualize multicolumn data data
    :param str _column_x: column for horizontal axis
    :param str _column_x: column for vertical axis
    :param str _scale: the scale to use for plotting data (linear by default, but could use log, log2, log10)  
    :param int _width_pixels: the width of the final plot in pixels  
    """
    #if '_backend' not in locals(): uti_plot_init() #?
    _backend.uti_plot_data_file(_fname, _read_labels, _e, _x, _y, _graphs_joined,
                                _multicolumn_data, _column_x, _column_y, #MR31102017
                                _scale, _width_pixels)

#def uti_data_file_plot(_fname, _read_labels=1, _e=0, _x=0, _y=0, _graphs_joined=True):
#def uti_data_file_plot(_fname, _read_labels=1, _e=0, _x=0, _y=0, _graphs_joined=True, _traj_report=False, _traj_axis='x'): #MR29072016
#def uti_data_file_plot(_fname, _read_labels=1, _e=0, _x=0, _y=0, _graphs_joined=True, _traj_report=False, _traj_axis='x', _scale='linear', _width_pixels=None): #MR20012017  
def uti_data_file_plot(_fname, _read_labels=1, _e=0, _x=0, _y=0, _graphs_joined=True,
                       _multicolumn_data=False, _column_x=None, _column_y=None, #MR31102017
                       _scale='linear', _width_pixels=None):
    """Generate plot from configuration in _fname

    :param str _fname: config loaded from here
    :param bool _read_labels: whether to read labels from _fname
    :param float _e: photon energy adjustment
    :param float _x: horizonal position adjustment
    :param float _y: vertical position adjustment
    :param bool _graphs_joined: if true, all plots in a single figure
    :param bool _multicolumn_data: if true, visualize multicolumn data data
    :param str _column_x: column for horizontal axis
    :param str _column_x: column for vertical axis
    :param str _scale: the scale to use for plotting data (linear by default, but could use log, log2, log10)  
    :param int _width_pixels: the width of the final plot in pixels  
    """
    #if '_backend' not in locals(): uti_plot_init() #?
    #_backend.uti_data_file_plot(_fname, _read_labels, _e, _x, _y, _graphs_joined)
    #_backend.uti_data_file_plot(_fname, _read_labels, _e, _x, _y, _graphs_joined, _traj_report, _traj_axis) #MR29072016
    #_backend.uti_data_file_plot(_fname, _read_labels, _e, _x, _y, _graphs_joined, _traj_report, _traj_axis, _scale, _width_pixels) #MR20012017  
    #_backend.uti_data_file_plot(_fname, _read_labels, _e, _x, _y, _graphs_joined,
    #                            _multicolumn_data, _column_x, _column_y, #MR31102017
    #                            _scale, _width_pixels)
    uti_plot_data_file(_fname, _read_labels, _e, _x, _y, _graphs_joined, _multicolumn_data, _column_x, _column_y, _scale, _width_pixels) #OC16112017

class _BackendBase(object):
    def __getattr__(self, attr):
        return self._backend_call

class _BackendMissing(_BackendBase):
    def _backend_call(self, *args, **kwargs):
        uti_plot_init()
        method_name = sys._getframe(1).f_code.co_name
        func = getattr(_backend, method_name)
        return func(*args)

class _BackendNone(_BackendBase):
    def _backend_call(*args, **kwargs):
        pass

_backend = _BackendMissing()
