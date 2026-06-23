"""Material Utilities Module


Modules:

    calc_refl_arr
    calc_coated_refl_arr
    calc_multilayer_refl_arr
    add_mat


.. moduleauthor:: Nathan Whittington
.. moduleauthor:: Roman Chernikov
"""
# updated 22-06-2026

from __future__ import print_function  # Python 2.7 compatibility

import inspect
from array import array

import numpy as np


xraydb = None
xraydb_found = False
xraydb_version = ''

try:
    import xraydb
    xraydb_found = True
    xraydb_version = getattr(xraydb, '__version__', '')
except Exception as exc:
    print("Warning: 'xraydb' can not be loaded: {}. Reflectivity is set to 1.".format(exc))


def _xraydb_function(_name):
    if not xraydb_found:
        return None
    return getattr(xraydb, _name, None)


mirror_reflectivity = _xraydb_function('mirror_reflectivity')
coated_reflectivity = _xraydb_function('coated_reflectivity')
multilayer_reflectivity = _xraydb_function('multilayer_reflectivity')
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


def calc_refl_arr(
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


def calc_coated_refl_arr(
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


def calc_multilayer_refl_arr(
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
