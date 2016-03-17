"""
utils: Utility function
------------------------------------------------------------------------
    
 Ver       Date      Who  Comments
-------  ----------  ---  ------------------------------------------
0.0.963  2016-03-14  JRM  Add pa_to_torr and torr_to_pa
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
