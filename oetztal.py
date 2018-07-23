from core import *
from plots import *

import os
from copy import deepcopy
from functools import partial

import matplotlib.pyplot as plt
import xarray as xr

import salem
from oggm import cfg, workflow, tasks, utils, graphics
from oggm.utils import get_demo_file
from oggm.core.flowline import FluxBasedModel, FileModel
FlowlineModel = partial(FluxBasedModel, inplace=False)

if __name__ == '__main__':
    cfg.initialize()
    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        cfg.PATHS['working_dir'] = os.environ.get("S_WORKDIR")
    else:
        #WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction'
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/find_initial_state/past_state_information'
        utils.mkdir(WORKING_DIR, reset=False)
        cfg.PATHS['working_dir'] = WORKING_DIR

    cfg.PATHS['plot_dir'] = os.path.join(cfg.PATHS['working_dir'], 'plots')
    cfg.PATHS['dem_file'] = get_demo_file('srtm_oetztal.tif')
    cfg.PATHS['climate_file'] = get_demo_file('HISTALP_oetztal.nc')

    # Use multiprocessing?
    cfg.PARAMS['use_multiprocessing'] = True

    # How many grid points around the glacier?
    cfg.PARAMS['border'] = 200

    cfg.PARAMS['run_mb_calibration'] = True
    cfg.PARAMS['optimize_inversion_params'] = False
    cfg.PARAMS['use_intersects'] = False

    # add to BASENAMES
    _doc = 'contains observed and searched glacier from synthetic experiment to find intial state'
    cfg.BASENAMES['synthetic_experiment'] = ('synthetic_experiment.pkl', _doc)

    plt.rcParams['figure.figsize'] = (8, 8)  # Default plot size

    # initialization
    rgi = get_demo_file('rgi_oetztal.shp')
    rgidf = salem.read_shapefile(rgi)
    gdirs = workflow.init_glacier_regions(rgidf)
    workflow.execute_entity_task(tasks.glacier_masks, gdirs)
    prepare_for_initializing(gdirs)
    synthetic_experiments_parallel(gdirs)
    volumes = pd.DataFrame()
    years = np.arange(1850, 1970, 170)
    for gdir in gdirs[:1]:

        if gdir.rgi_id != 'RGI50-11.00945':

            for yr in years:
                #find_possible_glaciers(gdir,gdir.read_pickle('synthetic_experiment'),yr)
                path = os.path.join(gdir.dir, 'result' + str(yr) + '.pkl')

                if os.path.isfile(path) and not pd.read_pickle(path).empty :
                    df = pd.read_pickle(path)

                else:
                    df = get_single_results(gdir, yr, gdir.read_pickle('synthetic_experiment'))
                    df.to_pickle(path)

                df['volume'] = df['model'].apply(lambda x: x.volume_m3)
                df['temp_bias'] = df['temp_bias'].apply(lambda x: float(x))
                max_model = deepcopy(df.loc[df.volume.idxmax(),'model'])

                max_obj = df.loc[df.volume.idxmax(),'objective']

                rp = gdir.get_filepath('model_run', filesuffix='experiment')
                ex_mod = FileModel(rp)
                v1850 = deepcopy(ex_mod.volume_m3)
                ex_mod.run_until(2000)
                v2000 = ex_mod.volume_m3
                volumes = volumes.append({'rgi':gdir.rgi_id,'ratio(experiment)': v1850/v2000,'ratio(max)':df.volume.max()/v2000,'objective(max)':max_obj, 'temp_bias':df.temp_bias.min()},ignore_index=True)
            plot_surface_col(gdir, df, ex_mod, 1850, '')
            fig, ax = plt.subplots(1, 2)
            max_model.reset_y0(1850)
            graphics.plot_modeloutput_map(gdirs=gdir, model=max_model, ax=ax[0],modelyr=1850)
            ex_mod.reset_y0(1850)
            graphics.plot_modeloutput_map(gdirs=gdir, model=ex_mod, ax=ax[1])

    volumes = volumes.set_index('rgi').sort_values('objective(max)')
    cols = ['ratio(experiment)','ratio(max)', 'temp_bias','objective(max)']

    with pd.option_context('display.max_rows', None, 'display.max_columns',None):
        print(volumes.reindex(columns=cols))


    plt.show()
