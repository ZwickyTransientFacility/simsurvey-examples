# simsurvey-examples
Collection of examples of using the simsurvey package to simulate transient lightcurves

To run this code you will require `simsurvey` and its dependencies, which can be install using pip: ```pip install simsurvey```

The notebooks will instruct you to download additional files.

## Scripts

We include the script `run_sims.py` as an example script that uses `population_utils.py`. 
- Currently this works with the simsurvey in branch `rb` of simsurvey : https://github.com/ZwickyTransientFacility/simsurvey/tree/rb 
- Also download `df_eg_sim.csv.gz` from dropbox and replace the `data/df_eg_sim.csv.gz` by it
- You will have to set the `sfd98_dir` variable in `run_sims.py`
- You can work out the options from ```run_sims.py -h```, if you run this with default parameters after downloading the file from dropbox above, the program should run in a few minutes and produce a `data/lcs_SALT2_params.pkl` file. If you open the file using the usual `simsurvey` methods from this directory:
```
import simsurvey
lcs =  simsurvey.LightCurveCollection(load='data/lcs_SALT2_params.pkl')
print(f'The number of detected transients is {len(lcs.meta)}\n')
```
This will print 'The number of detected transients is n' where n is ~ 296.
