"""
utils: Utility function
------------------------------------------------------------------------
    
 Ver       Date      Who  Comments
-------  ----------  ---  ------------------------------------------
0.0.963  2016-03-14  JRM  Add pa_to_torr and torr_to_pa
0.0.964  2016-03-22  JRM  Added calc_adda_mm_per_px
"""
# -*- coding: utf-8 -*-
def pa_to_torr(pa):
    """pa_to_torr(pa)

    Convert pressure in Pascals to Torr

    Parameters
    ----------
    pa: number
        Pressure, in Pa

    Returns
    -------
    torr: number
        Pressure, in Torr
    """
    return(0.0075062*pa)

def torr_to_pa(torr):
    """torr_to_pa(torr)

    Convert pressure in Torr to Pascals

    Parameters
    ----------
    torr: number
        Pressure, in torr

    Returns
    -------
    pa: number
        Pressure, in Pascals
    """
    return(torr/0.0075062)

def calc_adda_mm_per_px(px_per_unit, x_px=1024, unit_fact=1.0e-9):
    """calc_adda_mm_per_px(px_per_unit, x_px=1024, unit_fact=1.0e-9)

    Compute the mm per pixel for the maximum pixel size (4096) of a Soft
    Imaging Systems ADDA-II slow scan interface. This is entered into
    the channel calibration in the Analysis.ini file and is corrected
    for subsampling (e.g. to 1024 pixels) for a given image within the
    AnalySIS software

    Parameters
    ----------
    px_per_unit: number
        The number of pixels per unit (e.g. nm) that are measured from
        a calibration standard using something like the KMAG Imaging-C
        module.
    x_px: number (1024)
        The image width set for the Channel in the analySIS software.
        The default is set to 1024.
    unit_factor: number (1.0e-9)
        This is the mutiplier in for the scale in SI units. The default
        is 1.0e-9 for nm.

    Returns
    -------
    px_per_mm: number
        The pixels/mm used by the analySIS inverse magnification
        calibration:
            px_per mm = slope * (1/magnification)
        These values are entered into arrays, one for the X-AXIS (width)
        and one for the Y-AXIS (height) for the channel in the
        Analysis.ini file. N.B. Correct for lens hysteresis when 
        measuring standards...

    Example for the Hitachi 4100:

    > from jmToolsPy3 import calc_adda_mm_per_px
    > print("800X: X= ", calc_adda_mm_per_px(150.92),
            " Y= ", calc_adda_mm_per_px(148.77))
    """
    adda_max = 4096 # Max pix size for ADDA-II
    scale_fact = adda_max/x_px
    un_per_px = px_per_unit /scale_fact
    mm_per_px = un_per_px*1000.*unit_fact
    px_per_mm = 1.0/mm_per_px
    return(round(px_per_mm, 3))
