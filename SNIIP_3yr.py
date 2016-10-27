#! /usr/bin/env python
# -*- coding: utf-8 -*-

# See SNIIP_90d.ipynb for explanations 

import warnings
## No annoying warnings
warnings.filterwarnings('ignore')

import numpy as np
import simsurvey.cadence as simul
import sncosmo

import simsurvey_tools as sst

def random_parameters(redshifts, model,
                      mag=(-19.3, 0.1),
                      r_v=2., ebv_rate=0.11,
                      **kwargs):
    # Amplitude
    amp = []
    for z in redshifts:
        model.set(z=z)
        mabs = np.random.normal(mag[0], mag[1])
        model.set_source_peakabsmag(mabs, 'bessellb', 'vega')
        amp.append(model.get('amplitude'))

    return {
        'amplitude': np.array(amp),
        'hostr_v': r_v * np.ones(len(redshifts)),
        'hostebv': np.random.exponential(ebv_rate, len(redshifts))
    }

if __name__ == '__main__':
    mjd_range = (58000, 59095)

    dust = sncosmo.CCM89Dust()
    model = sncosmo.Model(source='s11-2005lc',
                          effects=[dust],
                          effect_names=['host'],
                          effect_frames=['rest'])

    transientprop = dict(lcmodel=model,
                         lcsimul_func=random_parameters,
                         lcsimul_prop=dict(mag=(-16.75, 0.98)))

    tr = simul.get_transient_generator([0.0, 0.1], ratekind='custom',
                                       ratefunc=lambda z: 1.5e-4,
                                       dec_range=[-40,90],
                                       mjd_range=[mjd_range[0] - 90,
                                                  mjd_range[1] + 25],
                                       transientprop=transientprop)


    survey, lcs = sst.get_survey_lcs_simple(tr, mjd_range=mjd_range)

    idx = lcs.meta['idx_orig']
    n_obs = np.zeros(survey.generator.ntransient)
    n_obs[idx] = np.array([len(a) for a in lcs])

    print 'SNe observed: %i out of %i'%(np.sum(n_obs > 0),
                                        survey.generator.ntransient)

    lcs.save('lcs_SNIIp_3yr.pkl')
