"""Material Utilities Module


Modules:

    calc_refl
    calc_refl_coated
    calc_refl_multilayer
    calc_delta_atten_len
    srwl_opt_setup_transm_from_material
    add_mat


.. moduleauthor:: Nathan Whittington
.. moduleauthor:: Roman Chernikov
"""
# updated 22-06-2026

from __future__ import print_function  # Python 2.7 compatibility

import inspect
from array import array

import numpy as np

try:
    from srwlib import SRWLOptT
except Exception:
    from .srwlib import SRWLOptT


xraydb = None
xraydb_found = False
xraydb_version = ''

try:
    import xraydb
    xraydb_found = True
    xraydb_version = getattr(xraydb, '__version__', '')
except Exception as exc:
    print("Warning: 'xraydb' can not be loaded: {}. Reflectivity is set to 1 and refractive optical elements use vacuum properties.".format(exc))


def _xraydb_function(_name):
    if not xraydb_found:
        return None
    return getattr(xraydb, _name, None)


mirror_reflectivity = _xraydb_function('mirror_reflectivity')
coated_reflectivity = _xraydb_function('coated_reflectivity')
multilayer_reflectivity = _xraydb_function('multilayer_reflectivity')
xray_delta_beta = _xraydb_function('xray_delta_beta')
get_material = _xraydb_function('get_material')
add_material = _xraydb_function('add_material')


def _warn_refl_fallback(_message):
    print("Warning: {} Reflectivity is set to 1.".format(_message))
    return 1


def _check_xraydb_function(_func, _func_name):
    if not xraydb_found:
        return _warn_refl_fallback("'xraydb' is not available.")
    if _func is None:
        version_text = " version {}".format(xraydb_version) if xraydb_version else ''
        return _warn_refl_fallback(
            "'{}' is not available in xraydb{}.".format(_func_name, version_text)
        )
    try:
        if 'output' not in inspect.signature(_func).parameters:
            version_text = " version {}".format(xraydb_version) if xraydb_version else ''
            return _warn_refl_fallback(
                "'{}' in xraydb{} does not support complex amplitude output.".format(
                    _func_name, version_text
                )
            )
    except (TypeError, ValueError):
        pass
    return None


def _check_refl_parameters(_n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin,
                           _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type):
    try:
        n_ph_en = int(_n_ph_en)
        n_ang = int(_n_ang)
        n_comp = int(_n_comp)
    except (TypeError, ValueError):
        return _warn_refl_fallback(
            "Numbers of photon energy, grazing angle, and component points must be integers."
        )
    if n_ph_en < 2:
        return _warn_refl_fallback(
            "At least two photon energy points are required by SRW reflectivity interpolation."
        )
    if n_ang < 2:
        return _warn_refl_fallback(
            "At least two grazing angle points are required by SRW reflectivity interpolation."
        )
    if n_comp not in (1, 2):
        return _warn_refl_fallback("Number of reflectivity components must be 1 or 2.")
    if _ph_en_scale_type not in ('lin', 'log'):
        return _warn_refl_fallback(
            "Invalid photon energy scale type '{}'; use 'lin' or 'log'.".format(_ph_en_scale_type)
        )
    if _ang_scale_type not in ('lin', 'log'):
        return _warn_refl_fallback(
            "Invalid grazing angle scale type '{}'; use 'lin' or 'log'.".format(_ang_scale_type)
        )
    if (_ph_en_scale_type == 'log') or (_ang_scale_type == 'log'):
        return _warn_refl_fallback(
            "Logarithmic sampling is not supported by SRW reflectivity interpolation."
        )
    try:
        if (_ph_en_start <= 0) or (_ph_en_fin <= 0):
            return _warn_refl_fallback("Photon energy limits must be positive.")
        if (_ang_start < 0) or (_ang_fin < 0):
            return _warn_refl_fallback("Grazing angle limits can not be negative.")
        if _ph_en_start == _ph_en_fin:
            return _warn_refl_fallback("Photon energy limits must define a non-zero range.")
        if _ang_start == _ang_fin:
            return _warn_refl_fallback("Grazing angle limits must define a non-zero range.")
    except TypeError:
        return _warn_refl_fallback("Photon energy and grazing angle limits must be numbers.")
    return None


def _sampling(_start, _fin, _n, _scale_type):
    if _scale_type == 'lin':
        return np.linspace(_start, _fin, _n)
    print("Warning: logarithmic scale is not supported by SRW reflectivity interpolation yet.")
    return np.logspace(np.log10(_start), np.log10(_fin), _n)


def _material_is_defined(_material, _density):
    if _density is not None:
        return True
    if (get_material is None) or (_material is None):
        return False
    try:
        return get_material(_material) is not None
    except Exception:
        return False


def _check_material(_material, _density, _description='material'):
    if _material_is_defined(_material, _density):
        return None
    return _warn_refl_fallback(
        "Density of {} '{}' was not found; specify the density or register the "
        "material with add_mat().".format(_description, _material)
    )


def _refractive_fallback(_ph_en, _message):
    print("Warning: {} Refractive index decrement is set to 0 and attenuation is disabled.".format(_message))
    try:
        energy = np.asarray(_ph_en, dtype=float)
    except (TypeError, ValueError):
        return 0.0, 1.e+23
    if energy.ndim == 0:
        return 0.0, 1.e+23
    return array('d', [0.0]*energy.size), array('d', [1.e+23]*energy.size)


def _resolve_material(_material, _density):
    if _density is not None:
        try:
            density = float(_density)
        except (TypeError, ValueError):
            return None
        if (not np.isfinite(density)) or (density <= 0):
            return None
        return _material, density
    if get_material is None:
        return None
    try:
        return get_material(_material)
    except Exception:
        return None


def _refl_to_srw_array(_refl_s, _refl_p, _n_tot):
    try:
        refl = np.ravel(_refl_s).view(np.float64)
        if _refl_p is not None:
            refl = np.concatenate((refl, np.ravel(_refl_p).view(np.float64)))
    except (TypeError, ValueError) as exc:
        return _warn_refl_fallback(
            "Calculated reflectivity can not be converted to an SRW array: {}.".format(exc)
        )
    if _n_tot != len(refl):
        return _warn_refl_fallback(
            "Calculated reflectivity array length {} does not match expected length {}.".format(
                len(refl), _n_tot
            )
        )
    if not np.all(np.isfinite(refl)):
        return _warn_refl_fallback("Calculated reflectivity contains invalid values.")
    return array('d', refl)


def calc_refl(
    _mat,
    _n_ph_en=1,
    _n_ang=1,
    _n_comp=1,
    _ph_en_start=0,
    _ph_en_fin=0,
    _ph_en_scale_type='lin',
    _ang_start=0,
    _ang_fin=0,
    _ang_scale_type='lin',
    _dens=None,
    _roughness=0.0
):
    """Calculate complex mirror reflectivity vs photon energy and grazing angle.

    The result is a C-aligned flat array ordered by polarization, grazing angle,
    photon energy, and real / imaginary component. If xraydb or required material
    data are unavailable, a warning is printed and perfect reflectivity is used.
    """
    fallback = _check_xraydb_function(mirror_reflectivity, 'mirror_reflectivity')
    if fallback is not None:
        return fallback
    fallback = _check_refl_parameters(
        _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin,
        _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type
    )
    if fallback is not None:
        return fallback
    fallback = _check_material(_mat, _dens)
    if fallback is not None:
        return fallback

    n_ph_en = int(_n_ph_en)
    n_ang = int(_n_ang)
    n_comp = int(_n_comp)
    ph_en = _sampling(_ph_en_start, _ph_en_fin, n_ph_en, _ph_en_scale_type)
    ang = _sampling(_ang_start, _ang_fin, n_ang, _ang_scale_type)

    try:
        refl_s = mirror_reflectivity(
            _mat, ang, ph_en, density=_dens, roughness=_roughness,
            polarization='s', output='amplitude'
        )
        refl_p = None
        if n_comp == 2:
            refl_p = mirror_reflectivity(
                _mat, ang, ph_en, density=_dens, roughness=_roughness,
                polarization='p', output='amplitude'
            )
    except Exception as exc:
        return _warn_refl_fallback(
            "xraydb mirror reflectivity calculation for '{}' failed: {}.".format(_mat, exc)
        )

    return _refl_to_srw_array(refl_s, refl_p, n_ph_en*n_ang*n_comp*2)


def calc_refl_coated(
    _coating,
    _coating_thick,
    _substr,
    _n_ph_en=1,
    _n_ang=1,
    _n_comp=1,
    _ph_en_start=0,
    _ph_en_fin=0,
    _ph_en_scale_type='lin',
    _ang_start=0,
    _ang_fin=0,
    _ang_scale_type='lin',
    _coating_dens=None,
    _substr_dens=None,
    _binder=None,
    _binder_dens=None,
    _binder_thick=0.0,
    _surf_roughness=0.0,
    _substr_roughness=0.0,
):
    """Calculate complex reflectivity of a coated mirror."""
    fallback = _check_xraydb_function(coated_reflectivity, 'coated_reflectivity')
    if fallback is not None:
        return fallback
    fallback = _check_refl_parameters(
        _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin,
        _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type
    )
    if fallback is not None:
        return fallback
    if _coating_thick is None:
        return _warn_refl_fallback("Coating thickness is not specified.")
    fallback = _check_material(_coating, _coating_dens, 'coating material')
    if fallback is not None:
        return fallback
    fallback = _check_material(_substr, _substr_dens, 'substrate material')
    if fallback is not None:
        return fallback
    if _binder is not None:
        fallback = _check_material(_binder, _binder_dens, 'binder material')
        if fallback is not None:
            return fallback

    n_ph_en = int(_n_ph_en)
    n_ang = int(_n_ang)
    n_comp = int(_n_comp)
    ph_en = _sampling(_ph_en_start, _ph_en_fin, n_ph_en, _ph_en_scale_type)
    ang = _sampling(_ang_start, _ang_fin, n_ang, _ang_scale_type)
    kwargs = dict(
        coating_dens=_coating_dens,
        surface_roughness=_surf_roughness,
        substrate_dens=_substr_dens,
        substrate_roughness=_substr_roughness,
        binder=_binder,
        binder_thick=_binder_thick,
        binder_dens=_binder_dens,
        output='amplitude',
    )

    try:
        refl_s = coated_reflectivity(
            _coating, _coating_thick, _substr, ang, ph_en,
            polarization='s', **kwargs
        )
        refl_p = None
        if n_comp == 2:
            refl_p = coated_reflectivity(
                _coating, _coating_thick, _substr, ang, ph_en,
                polarization='p', **kwargs
            )
    except Exception as exc:
        return _warn_refl_fallback(
            "xraydb coated reflectivity calculation failed: {}.".format(exc)
        )

    return _refl_to_srw_array(refl_s, refl_p, n_ph_en*n_ang*n_comp*2)


def calc_refl_multilayer(
    _stackup,
    _thickness,
    _substr,
    _n_periods,
    _n_ph_en=1,
    _n_ang=1,
    _n_comp=1,
    _ph_en_start=0,
    _ph_en_fin=0,
    _ph_en_scale_type='lin',
    _ang_start=0,
    _ang_fin=0,
    _ang_scale_type='lin',
    _dens=None,
    _substr_dens=None,
    _substr_rough=0,
    _surf_rough=0,
):
    """Calculate complex reflectivity of a periodic multilayer mirror."""
    fallback = _check_xraydb_function(multilayer_reflectivity, 'multilayer_reflectivity')
    if fallback is not None:
        return fallback
    fallback = _check_refl_parameters(
        _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin,
        _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type
    )
    if fallback is not None:
        return fallback
    if (_stackup is None) or (_thickness is None):
        return _warn_refl_fallback("Multilayer materials and thicknesses must be specified.")
    try:
        n_layers = len(_stackup)
        if n_layers != len(_thickness):
            return _warn_refl_fallback(
                "Number of multilayer materials does not match number of thicknesses."
            )
    except TypeError:
        return _warn_refl_fallback("Multilayer materials and thicknesses must be sequences.")
    try:
        n_periods = int(_n_periods)
    except (TypeError, ValueError):
        return _warn_refl_fallback("Number of multilayer periods must be an integer.")
    if n_periods < 1:
        return _warn_refl_fallback("Number of multilayer periods must be positive.")
    if _dens is not None:
        try:
            if len(_dens) != n_layers:
                return _warn_refl_fallback(
                    "Number of multilayer densities does not match number of materials."
                )
        except TypeError:
            return _warn_refl_fallback("Multilayer densities must be a sequence.")

    densities = [None]*len(_stackup) if _dens is None else _dens
    for material, density in zip(_stackup, densities):
        fallback = _check_material(material, density, 'multilayer material')
        if fallback is not None:
            return fallback
    fallback = _check_material(_substr, _substr_dens, 'substrate material')
    if fallback is not None:
        return fallback

    n_ph_en = int(_n_ph_en)
    n_ang = int(_n_ang)
    n_comp = int(_n_comp)
    ph_en = _sampling(_ph_en_start, _ph_en_fin, n_ph_en, _ph_en_scale_type)
    ang = _sampling(_ang_start, _ang_fin, n_ang, _ang_scale_type)
    kwargs = dict(
        n_periods=n_periods,
        density=_dens,
        substrate_density=_substr_dens,
        substrate_rough=_substr_rough,
        surface_rough=_surf_rough,
        output='amplitude',
    )

    try:
        refl_s = multilayer_reflectivity(
            _stackup, _thickness, _substr, ang, ph_en,
            polarization='s', **kwargs
        )
        refl_p = None
        if n_comp == 2:
            refl_p = multilayer_reflectivity(
                _stackup, _thickness, _substr, ang, ph_en,
                polarization='p', **kwargs
            )
    except Exception as exc:
        return _warn_refl_fallback(
            "xraydb multilayer reflectivity calculation failed: {}.".format(exc)
        )

    return _refl_to_srw_array(refl_s, refl_p, n_ph_en*n_ang*n_comp*2)


def calc_delta_atten_len(_mat, _ph_en, _dens=None):
    """Calculate material properties used by SRW refractive optical elements.

    :param _mat: material name or chemical formula
    :param _ph_en: photon energy [eV], or a sequence of photon energies
    :param _dens: material density [g/cm^3]; if omitted, use the xraydb material database
    :return: refractive index decrement and attenuation length [m]

    Scalar photon energy returns two floats. A sequence returns two ``array('d')``
    objects suitable for spectrally-dependent SRW transmission elements and CRLs.
    """
    if (not xraydb_found) or (xray_delta_beta is None):
        return _refractive_fallback(_ph_en, "'xraydb.xray_delta_beta' is unavailable.")

    try:
        energy = np.asarray(_ph_en, dtype=float)
    except (TypeError, ValueError):
        return _refractive_fallback(_ph_en, "Photon energy must be numeric.")
    if (energy.size < 1) or (not np.all(np.isfinite(energy))) or np.any(energy <= 0):
        return _refractive_fallback(_ph_en, "Photon energy values must be finite and positive.")

    material = _resolve_material(_mat, _dens)
    if material is None:
        return _refractive_fallback(
            _ph_en,
            "Density of material '{}' was not found; specify _dens or register the material with add_mat().".format(_mat)
        )
    formula, density = material

    try:
        energy_arg = float(energy) if energy.ndim == 0 else energy
        delta, _beta, atten_len_cm = xray_delta_beta(formula, density, energy_arg)
        delta = np.asarray(delta, dtype=float)
        atten_len_m = 0.01*np.asarray(atten_len_cm, dtype=float)
    except Exception as exc:
        return _refractive_fallback(
            _ph_en, "xraydb refractive-property calculation for '{}' failed: {}.".format(_mat, exc)
        )

    if (not np.all(np.isfinite(delta))) or (not np.all(np.isfinite(atten_len_m))) or np.any(atten_len_m <= 0):
        return _refractive_fallback(
            _ph_en, "Calculated refractive properties for '{}' contain invalid values.".format(_mat)
        )

    if energy.ndim == 0:
        return float(delta), float(atten_len_m)
    return array('d', delta.ravel()), array('d', atten_len_m.ravel())


def srwl_opt_setup_transm_from_material(
    _mat,
    _thick,
    _ph_en,
    _dens=None,
    _rx=1.e-03,
    _ry=1.e-03,
    _nx=2,
    _ny=2,
    _x=0,
    _y=0,
    _ext_tr=1,
    _fx=1.e+23,
    _fy=1.e+23
):
    """Set up a uniform transmission element from material data.

    :param _mat: material name or chemical formula
    :param _thick: material thickness [m]
    :param _ph_en: photon energy [eV], or a sequence of photon energies
    :param _dens: material density [g/cm^3]; if omitted, use the xraydb material database
    :param _rx: horizontal coordinate range [m]
    :param _ry: vertical coordinate range [m]
    :param _nx: number of points vs horizontal position
    :param _ny: number of points vs vertical position
    :param _x: horizontal transverse coordinate of center [m]
    :param _y: vertical transverse coordinate of center [m]
    :param _ext_tr: transmission outside the grid/mesh is zero (0), or same as boundary (1)
    :param _fx: estimated focal length in the horizontal plane [m]
    :param _fy: estimated focal length in the vertical plane [m]
    :return: SRWLOptT by default
    """
    delta, atten_len = calc_delta_atten_len(_mat, _ph_en, _dens)

    try:
        thick = float(_thick)
        if (not np.isfinite(thick)) or (thick < 0):
            raise ValueError
    except (TypeError, ValueError):
        print("Warning: material thickness is invalid. Transmission element is set to vacuum.")
        thick = 0.0

    try:
        energy = np.asarray(_ph_en, dtype=float).ravel()
        if energy.size < 1:
            raise ValueError
    except (TypeError, ValueError):
        print("Warning: photon energy is invalid. Transmission element mesh energy is set to 0.")
        energy = np.asarray([0.0])

    ne = int(energy.size)
    try:
        nx = int(_nx)
        ny = int(_ny)
        if (nx <= 0) or (ny <= 0):
            raise ValueError
    except (TypeError, ValueError):
        print("Warning: transmission mesh dimensions are invalid. A 2 x 2 mesh is used.")
        nx = 2
        ny = 2

    delta_arr = np.asarray(delta, dtype=float).ravel()
    atten_len_arr = np.asarray(atten_len, dtype=float).ravel()
    if (delta_arr.size != ne) or (atten_len_arr.size != ne):
        if (delta_arr.size == 1) and (atten_len_arr.size == 1):
            delta_arr = np.full(ne, delta_arr[0])
            atten_len_arr = np.full(ne, atten_len_arr[0])
        else:
            print("Warning: material data size is inconsistent with photon energy mesh. Transmission element is set to vacuum.")
            delta_arr = np.zeros(ne)
            atten_len_arr = np.full(ne, 1.e+23)

    if (not np.all(np.isfinite(delta_arr))) or (not np.all(np.isfinite(atten_len_arr))) or np.any(atten_len_arr <= 0):
        print("Warning: material data contain invalid values. Transmission element is set to vacuum.")
        delta_arr = np.zeros(ne)
        atten_len_arr = np.full(ne, 1.e+23)

    amp = np.exp(-0.5*thick/atten_len_arr)
    opd = -delta_arr*thick
    ar_tr_one_point = np.empty(2*ne, dtype=float)
    ar_tr_one_point[0::2] = amp
    ar_tr_one_point[1::2] = opd
    ar_tr = array('d', np.tile(ar_tr_one_point, nx*ny))

    op_t = SRWLOptT(
        _nx=nx, _ny=ny, _rx=_rx, _ry=_ry, _arTr=ar_tr, _extTr=_ext_tr,
        _Fx=_fx, _Fy=_fy, _x=_x, _y=_y, _ne=ne,
        _eStart=float(energy[0]), _eFin=float(energy[-1])
    )

    return op_t

def add_mat(name, formula, density, categories=None):
    """Add a material to the user-local xraydb material database."""
    if add_material is None:
        print("Warning: 'xraydb.add_material' is not available. Material was not added.")
        return
    if density is None:
        print("Warning: density of '{}' is not specified. Material was not added.".format(formula))
        return
    try:
        add_material(name, formula, density, categories)
    except Exception as exc:
        print("Warning: material '{}' was not added: {}".format(name, exc))
