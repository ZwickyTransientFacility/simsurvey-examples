#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import simsurvey.cadence as simul

def load_ztf_fields(filename='data/ZTF_Fields.txt', mwebv=False):
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

    return out

def get_survey_plan_single_field(mjd_range=(58000, 58020),
                                 width=10.,
                                 height=10.,
                                 dt=3600.,
                                 bands=['r_r_gggg', 'r_r_ggii'],
                                 band_str='des%s',                       
                                 skynoise={'g': 800., 'r': 800., 'i': 1260.}):
    """Generates a very simple plan with just one field. This can be used to
    simulate some lightcurves without specifying a full survey.
    See `SNIa_early_lcs.py/ipynb` for examples of use
    
    Parameters
    ----------
    mjd_range: [tuple of floats]
        Start and end date (in units of days) for the plan

    width, height: [floats]
        Width and height of the single survey field

    dt: [float]
        Time interval between two observations a night

    bands: [list of strs]
        List of filters that should be used for the observations;
        The plan will alternate between the strings each night and
        each string specifies the bands (as single letters) to be used
        with separation `dt`. Use underscore ("_") to skip an observation.

    band_str: [str]
        Formatting string used to create the full name of the bandpass
        as registered in `sncosmo`.

    skynoise: [dict]
        Dictionary of skynoise values to be used for each band.

    Return
    ------
    SurveyPlan object
    """
    obs_days = np.arange(mjd_range[0], mjd_range[1]+1, 1)

    obs = {'time': [], 'field': [], 'band': [], 'skynoise': []}

    for l, d in enumerate(obs_days):
        for k, band in enumerate(bands[l % len(bands)]):
            if band != '_':
                obs['field'].append(0)
                obs['time'].append(d + k * dt / 86400.)    
                obs['band'].append(band_str%band)
                obs['skynoise'].append(skynoise[band])

    return simul.SurveyPlan(time=obs['time'], band=obs['band'],
                            obs_field=obs['field'],
                            skynoise=obs['skynoise'],
                            fields=dict(ra=[0.], dec=[0.]),
                            width=width, height=height)

def get_survey_plan_simple(mjd_range=(58000, 59095),
                           dec_range=(-30, 90),
                           ra_range=(0, 120),
                           dec_good=30.,
                           mwebv_max=0.4,
                           t_obs=45.,
                           t_night=28801.,
                           n_repeat=80,
                           bands=['r_r_gggg', 'r_r_ggii'],
                           band_str='des%s',                       
                           bands_alt=['rrrrrrrr', 'rrrrrrrr'],
                           fields=None,
                           skynoise={'g': 800., 'r': 800., 'i': 1260.}):
    """Generates a more sophisticated plan using a list of field.
    This can be used to simulate a full survey until Eric's code .
    
    See `SNIa_90d.ipynb` for an example of use.
    
    Parameters
    ----------
    mjd_range: [tuple of floats]
        Start and end date (in units of days) for the plan

    dec_range: [tuple of floats]
        Allowed declination range for the centers of the fields.
    
    ra_range: [tuple of floats]
        Allowed right ascension range for the centers of the fields *during the
        first night*. This is shifted by a little bit each night to account for
        Earth orbit the Sun.
    
    dec_good: [float]
        "Good" declination value. The fields for each night will be sorted
        according to how close they are to this values. Only `2 * n_repeat`
        fields will be used.  The first `n_repeat` fields are considered the
        "good" fields (e.g. near zenith) while the second set is considered
        the "alternate" fields.
    
    mwebv_max: [float]
        Maximum MW E(B-V) value for field to be used.
    
    t_obs: [float]
        Average time (in seconds) needed per observation (including slewing
        and other overhead)

    t_night: [float]
        Total available time (in seconds) per night.

    n_repeat: [int]
        Number of repeats for one band before moving to the next (or repeating
        the same one). 
    
    bands: [list of strs]
        List of filters that should be used for the observations;
        The plan will alternate between the strings each night and
        each string specifies the bands (as single letters) to be used
        with separation `n_repeat` times with seperation `t_obs`.
        Use underscore ("_") to use `bands_alt` and the seconds list of
        `n_repeat` fields.

    band_str: [str]
        Formatting string used to create the full name of the bandpass
        as registered in `sncosmo`.

    bands_alt: [list of strs]
        Similar to `bands` but only use if.

    fields: [dict similar to output of `load_ztf_fields()`]
        Dictionary of field definitions; should include MW E(B-V).
    
    skynoise: [dict]
        Dictionary of skynoise values to be used for each band.
    
    Return
    ------
    SurveyPlan object
    """
    if fields is None:
        fields = load_ztf_fields(mwebv=True)
        for k in fields.keys():
            fields[k] = fields[k][:906]

    ra_shift = 360. / 365.25
    obs_days = np.arange(mjd_range[0], mjd_range[1]+1, 1)

    obs = {'time': [], 'field': [], 'band': [], 'skynoise': []}

    for l, d in enumerate(obs_days):
        # Find fields that have their center within dec_range and in ra_range shifted
        # by ra_shift for each day
        field_mask = (((fields['ra'] - ra_shift*(d-obs_days[0]) + 180)%360 - 180
                       > ra_range[0]) &
                      ((fields['ra'] - ra_shift*(d-obs_days[0]) + 180)%360 - 180
                       < ra_range[1]) &
                      (fields['dec'] > dec_range[0]) &
                      (fields['dec'] < dec_range[1]) &
                      (fields['mwebv'] < mwebv_max))
        field_idx = np.where(field_mask)[0]

        # Sort first by ra then by dec
        field_idx = field_idx[np.argsort(fields['ra'][field_idx])]
        field_idx = field_idx[np.argsort(np.abs(dec_good - fields['dec'][field_idx]))]
        fields_good = fields['field_id'][field_idx[:n_repeat]]
        fields_alt = fields['field_id'][field_idx[n_repeat:2*n_repeat]]

        for k, t in enumerate(np.arange(0, t_night, t_obs)):
            b0, b1 = (l%len(bands), (k / len(fields_good)) % len(bands[l%len(bands)]))
            band = bands[b0][b1]
            if band != '_':
                obs['field'].append(fields_good[k % len(fields_good)])
            else:
                obs['field'].append(fields_alt[k % len(fields_alt)])
                band = bands_alt[b0][b1]

            obs['time'].append(d + t / 86400.)    
            obs['band'].append(band_str%band)
            obs['skynoise'].append(skynoise[band])

    return simul.SurveyPlan(time=obs['time'], band=obs['band'], obs_field=obs['field'],
                            skynoise=obs['skynoise'],
                            fields={k: v for k, v in fields.items()
                                    if k in ['ra', 'dec', 'field_id',
                                             'width', 'height']})
    
def get_survey_lcs_single_field(tr, progress_bar=True, notebook=False,
                                instprop=None, **kw):
    """
    """
    plan = get_survey_plan_single_field(**kw)
    if instprop is None:
        instprop = {"desg":{"gain":1.,"zp":30,"zpsys":'ab'},
                    "desr":{"gain":1.,"zp":30,"zpsys":'ab'},
                    "desi":{"gain":1.,"zp":30,"zpsys":'ab'}}

    survey = simul.SimulSurvey(
        generator=tr, 
        plan=plan, 
        instprop=instprop
    )
    
    lcs = survey.get_lightcurves(progress_bar=progress_bar, notebook=notebook)
    
    return survey, lcs

def get_survey_lcs_simple(tr, progress_bar=True, notebook=False, instprop=None, **kw):
    """
    """
    plan = get_survey_plan_simple(**kw)
    if instprop is None:
        instprop = {"desg":{"gain":1.,"zp":30,"zpsys":'ab'},
                    "desr":{"gain":1.,"zp":30,"zpsys":'ab'},
                    "desi":{"gain":1.,"zp":30,"zpsys":'ab'}}

    survey = simul.SimulSurvey(
        generator=tr, 
        plan=plan, 
        instprop=instprop
    )
    
    lcs = survey.get_lightcurves(progress_bar=progress_bar, notebook=notebook)
    
    return survey, lcs
