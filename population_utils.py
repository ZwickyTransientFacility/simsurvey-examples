"""
A set of functions to help in running simulations with simsurvey. (This is not a script.)
"""
__all__ = ['get_plan', 'get_tg_params', 'msip_jdrange']


import time
import os
import numpy as np
import pandas as pd
import simsurvey
import sncosmo
from astropy.cosmology import Planck15
import marshaltools
import simsurvey_tools as sst

msip_jdrange = (2458183.0, 2458484.0) # ZTF MSIP 

def get_plan(observation_library, jd_min=msip_jdrange[0], jd_max=msip_jdrange[1],
             pad=(-0., 0.), programid=(1), return_lib=False):
    """
    Obtain the `simsurvey.SurveyPlan` within a certain range set by the arguments of the
    method from an observation library. The survey plan will contain observations from
    `jd_min + pad[0]` to `jd_max + pad[1]`

    
    Parameters
    ----------
    observation_library : `pd.DataFrame`
        observation library as a pandas dataFrame. Must have the columns,
        `fieldId`, `diff_maglim`, `jd`, `rcid`, `fid`
    jd_min : float, defaults to 2458183.
        minimum time in survey
    jd_max : float, defaults to 2458484.
        maximum time in survey
    pad : tuple of floats, defaults to (0., 0.)
        survey is designed from `jd_min + pad[0]` to `jd_max + pad[1]`

    Returns
    -------
    survey_plan : instance of `simsurvey.SurveyPlan` 
    """
    observation_library.fid.replace({1:'p48g', 2:'p48r', 3:'p48i'}, inplace=True)
    query = 'jd > @jd_min and jd < @jd_max and programid == @programid'
    survey_lib = observation_library.query(query)
    survey_lib['skynoise'] = 0.2 * 10.0**(- 0.4 * survey_lib.diff_maglim) 

    print(f'survey_lib has {len(survey_lib)} observations\n')
    if return_lib:
        return survey_lib

    # Load fields and ccds in general
    ccds = sst.load_ztf_ccds(filename='data/ZTF_corners_rcid.txt', num_segs=64)
    fields = sst.load_ztf_fields(filename='data/ZTF_Fields.txt')

    print(ccds)
    print(fields)
    # Plan
    plan = simsurvey.SurveyPlan(time=survey_lib.jd.values,
                                band=survey_lib.fid.values,
                                obs_ccd=survey_lib.rcid.values,
                                obs_field=survey_lib.fieldId.values,
                                skynoise=survey_lib.skynoise.values,
                                zp=np.zeros(len(survey_lib)), # time),
                                band_dict={'1': 'p48g', '2': 'p48r', '3': 'p48i'},
                                ccds=ccds,
                                fields=fields)

    return plan

def get_tg_params(sfd98_dir, mjd_range, stored_params,template='salt2',
                  zmax=0.05, mwrv=3.1 ):  
    """
    Create a Transient generator with a set of previously selected parameters 

    Parameters
    ----------
    sfd98_dir :
    
    mjd_range :
    
    stored_params :
    
    template :
    
    zmax :
    """
    tr = simsurvey.get_transient_generator((0.0, zmax),
                                           transient='Ia',
                                           template=template,
                                           ra_range=(0,360),
                                           dec_range=(-30,90),
                                           mjd_range=mjd_range,
                                           sfd98_dir=sfd98_dir)
    
    tr._derived_properties['simul_parameters']['zcmb'] = stored_params['z']

    # Important that the zcmb be set before this step (that sets tr.ntransient)
    r_, d_ = simsurvey.utils.random.radec(tr.ntransient,
                                          ra_range=tr.ra_range,
                                          dec_range=tr.dec_range)
    tr._derived_properties['simul_parameters']['ra'] = r_
    tr._derived_properties['simul_parameters']['dec'] = d_
    mjd = np.random.uniform(mjd_range[0], mjd_range[1], len(r_))
    tr._derived_properties['simul_parameters']['mjd'] = mjd

    # MW E(B-V) can just be reset becaus it will be reread from the maps when needed
    tr._reset_mwebv_()

    # Now regenerate the lightcurve model parameters
    lc = tr.transient_coverage["lightcurve_prop"]
    tr._derived_properties["simul_parameters"]["lightcurve"] = \
        dict(mwr_v=np.ones_like(r_) * mwrv,
        x0=stored_params['x0'].values,
        x1=stored_params['x1'].values,
        c=stored_params['c'].values)

    return tr
