#!/usr/bin/env python
"""
nohup python run_sims.py --fname_params SALT2_params > test.log 2>&1 &
Batch mode:
$parallel --jobs 2 nohup python run_sims.py --fname_params < fnames.txt > mock_sims.log 2>&1 &

# Inputs
1. Observation Library : This is a csv file which records the observation strategy
details : set by script argument --obs_lib
- assumed to be the filename relative to the expected location `data_dir` 
- If data file is stored in a differnt location, either use links or change the hard coded location below
- Note that the Obs library should be a csv file (gzipped is fine), with the columns noted in `run_simulations.get_plan`
2. Filename containing the properties for the sequence of transients that will be simulated
details : set by script argument --fname_params, assumed to be relative to `data_dir`
- the argument should not include the `.csv` extension
- should be a csv file with headers, and parameters as in example
- The essential columns are x0,x1,c,z. It is good to have idx. Note that t0 will be ignored, and ra, dec is
  not supplied. `simsurvey` will supply these variables. You can check the light curves in meta with
  `idx_orig` reproducing the `idx` here. The other columns are not necessary, but having them (or other
  columns will not cause problems)
3. sfd98_dir : set this to what makes sense for your computer
"""
import os
import numpy as np
import pandas as pd
from population_utils import get_plan, get_tg_params, msip_jdrange
from astropy.cosmology import Planck15 as cosmo
import simsurvey
import sys
import argparse
import time


# HARD CODED INPUTS HERE
homedir = os.path.expanduser('~')
data_dir = os.path.join('data')
sfd98_dir = os.path.join(homedir, 'data/MWDUST/')


ts = time.time()
parser = argparse.ArgumentParser(description='simulation script')
parser.add_argument('--obs_lib',
                    help='observation library name, defaults to "df_eg_sim.csv.gz", see details',
                    default='df_eg_sim.csv.gz')
parser.add_argument('--fname_params', help='param filename, defaults to ',
                    default='SALT2_params')

args = parser.parse_args()
print(f'simsurvey version = {simsurvey.__version__}')
print(f'simsurvey location = {simsurvey.__file__}')


obs_lib_fname = os.path.join(data_dir, args.obs_lib)

print(f'Working with observation library file {obs_lib_fname}')


fname_params = os.path.join('data',  f'{args.fname_params}.csv')

sys.stdout.flush()
## Output params

# output file name
fname_lcs = f'data/lcs_{args.fname_params}.pkl'

sys.stdout.flush()


print(f'The dustmap directory is {sfd98_dir}')
sys.stdout.flush()
obs_lib = pd.read_csv(obs_lib_fname)
print(f'Original Obs Lib had {len(obs_lib)} entries\n')
obs_lib = obs_lib.query('@msip_jdrange[0] < jd < @msip_jdrange[1]')
print(f'For the given survey Obs Lib had {len(obs_lib)} entries\n')

obs_lib.fid.replace({1:'p48g', 2:'p48r', 3:'p48i'}, inplace=True)
p = get_plan(obs_lib, return_lib=False)

## Generate SN outside the observation window to accomodate tails and rises
## Not actually necessary if we are reading the stored params
mjd_range = (msip_jdrange[0] - 70., msip_jdrange[1] + 30)
print(f'The time range of observation is {msip_jdrange[0]}, {msip_jdrange[1]}')
print(f'While the time range of supernova simulations is {msip_jdrange[0] - 70}, {msip_jdrange[1] + 30}')
sys.stdout.flush()
stored_params = pd.read_csv(fname_params)

print(f'We will simulate {len(stored_params)} in this simulation\n')
sys.stdout.flush()
tr = get_tg_params(sfd98_dir, mjd_range, stored_params, template='salt2', zmax=0.05)
survey = simsurvey.SimulSurvey(generator=tr, plan=p)

messg = (f'fname_params = {fname_params}\n',
         f'fname_lcs = {fname_lcs}\n',
         f'obs_lib_fname = {obs_lib_fname}\n')
print(messg)
print('Now beginning to simulate light curves\n')
lcs = survey.get_lightcurves(progress_bar=True, notebook=False) # If you get an error because of the progress_bar, delete this line.
sys.stdout.flush()

print('Done generating LCs, have to save\n')
sys.stdout.flush()
lcs.save(fname_lcs)
te = time.time()
deltat = (te - ts)/60.
print(f'Took {deltat} mins\n')
print('DONE\n')
