#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
## No annoying warnings
warnings.filterwarnings('ignore')

import numpy as np
import simsurvey.cadence as simul

import survey_tools as st

mjd_range = (58000, 58020)

tr = simul.get_sn_generator([0.,0.2], ratekind="basic", 
                            dec_range=[-1.,1.], ra_range=[-1.,1.],
                            mjd_range=[mjd_range[0] + 23, mjd_range[0] + 24],
                            ntransient=1000)

survey, lcs = st.get_survey_lcs_single_field(tr, mjd_range=mjd_range, notebook=True)

idx = lcs.meta['idx_orig']
n_obs = np.zeros(survey.generator.ntransient)
n_obs[idx] = np.array([len(a) for a in lcs])

print 'SNe observed: %i out of %i'%(np.sum(n_obs > 0), survey.generator.ntransient)

lcs.save('lcs_early.pkl')
