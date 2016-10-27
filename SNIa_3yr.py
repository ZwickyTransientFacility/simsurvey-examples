#! /usr/bin/env python
# -*- coding: utf-8 -*-

# See SNIa_90d.ipynb for explanations 

import warnings
## No annoying warnings
warnings.filterwarnings('ignore')

import numpy as np
import simsurvey.cadence as simul

import simsurvey_tools as sst

if __name__ == '__main__':
    mjd_range = (58000, 59095)

    tr = simul.get_sn_generator([0.,0.2], ratekind="basic", 
                                dec_range=[-40,90],
                                mjd_range=[mjd_range[0] - 60, mjd_range[-1] + 25])

    survey, lcs = sst.get_survey_lcs_simple(tr, mjd_range=mjd_range)

    idx = lcs.meta['idx_orig']
    n_obs = np.zeros(survey.generator.ntransient)
    n_obs[idx] = np.array([len(a) for a in lcs])

    print 'SNe observed: %i out of %i'%(np.sum(n_obs > 0), survey.generator.ntransient)

    lcs.save('lcs_SNIa_3yr.pkl')
