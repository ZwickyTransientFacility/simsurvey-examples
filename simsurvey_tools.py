#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sncosmo

def load_ztf_fields(filename='data/ZTF_Fields.txt', mwebv=False, galactic=False):
    """Load the ZTF fields propose by Eran from the included file.

    Parameters
    ----------
    filename: [str]
        File name of the ASCII file containing the field definitions

    mwebv: [bool]
        Include the Milky Way E(B-V) from the file in the output

    Return
    ------
    Dictionary of np.arrays with the field coordinates, IDs (and exinction)
    """
    fields = np.genfromtxt(filename, comments='%')

    out = {'field_id': np.array(fields[:,0], dtype=int),
           'ra': fields[:,1],
           'dec': fields[:,2]}

    if mwebv:
        out['mwebv'] = fields[:,3]
    if galactic:
        out['l'] = fields[:,4]
        out['b'] = fields[:,5]

    return out

def load_ztf_ccds(filename='data/ZTF_corners.txt', num_segs=16):
    """
    """
    ccd_corners = np.genfromtxt(filename, skip_header=1)
    ccds = [ccd_corners[4*k:4*k+4, :2] for k in range(num_segs)]

    return ccds

def load_ztf_filters():
    """
    """
    bands = {
        'ztfi' : 'data/ztfi_eff.txt',
        'ztfr' : 'data/ztfr_eff.txt',
        'ztfg' : 'data/ztfg_eff.txt',
    }

    for bandname in bands.keys() :
        fname = bands[bandname]
        b = np.loadtxt(fname)
        band = sncosmo.Bandpass(b[:,0], b[:,1], name=bandname)
        sncosmo.registry.register(band)
