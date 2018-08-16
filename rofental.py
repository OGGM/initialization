from core import *
from plots import *

import os
from copy import deepcopy
from functools import partial

import matplotlib.pyplot as plt
import xarray as xr
import geopandas as gpd
import shapely.geometry as shpg

import salem
from oggm import cfg, workflow, tasks, utils, graphics
from oggm.utils import get_demo_file
from oggm.core.flowline import FluxBasedModel, FileModel
FlowlineModel = partial(FluxBasedModel, inplace=False)

if __name__ == '__main__':
    cfg.initialize()
    ON_CLUSTER = True

    # Local paths
    if ON_CLUSTER:
        cfg.PATHS['working_dir'] = os.environ.get("S_WORKDIR")
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/rofental'
        #WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction'
        #WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/find_initial_state/past_state_information'
        utils.mkdir(WORKING_DIR, reset=False)
        cfg.PATHS['working_dir'] = WORKING_DIR

    cfg.PATHS['plot_dir'] = os.path.join(cfg.PATHS['working_dir'], 'plots')
    utils.mkdir(cfg.PATHS['plot_dir'], reset=False)

    # Use multiprocessing?
    cfg.PARAMS['use_multiprocessing'] = True

    # How many grid points around the glacier?
    cfg.PARAMS['border'] = 200

    # Set to True for operational runs
    cfg.PARAMS['continue_on_error'] = True

    # Use HISTALP climate file
    cfg.PARAMS['baseline_climate'] = 'HISTALP'

    # Get RGI-glaciers
    rgi = utils.get_rgi_region_file('11', version='61',reset=False)
    rgidf = salem.read_shapefile(rgi)

    # We use intersects
    db = utils.get_rgi_intersects_region_file(version='61',region='11')
    cfg.set_intersects_db(db)

    # Get the Rofental Basin file
    path = utils.get_demo_file('rofental_hydrosheds.shp')
    basin = gpd.read_file(path)

    # Take all glaciers in the Rhone Basin
    in_bas = [basin.geometry.contains(shpg.Point(x, y))[0] for
              (x, y) in zip(rgidf.CenLon, rgidf.CenLat)]
    rgidf = rgidf.loc[in_bas]
    # Store them for later
    rgidf.to_file(os.path.join(WORKING_DIR, 'rgi_rofental.shp'))

    # Sort for more efficient parallel computing
    rgidf = rgidf.sort_values('Area', ascending=False)

    cfg.PARAMS['run_mb_calibration'] = True
    cfg.PARAMS['optimize_inversion_params'] = False

    # add to BASENAMES
    _doc = 'contains observed and searched glacier from synthetic experiment to find intial state'
    cfg.BASENAMES['synthetic_experiment'] = ('synthetic_experiment.pkl', _doc)

    plt.rcParams['figure.figsize'] = (8, 8)  # Default plot size

    # initialization
    gdirs = workflow.init_glacier_regions(rgidf,reset=False)
    prepare_for_initializing(gdirs)
    synthetic_experiments_parallel(gdirs)

    years = [1850]
    for gdir in gdirs:
        try:
            for yr in years:
                find_possible_glaciers(gdir,gdir.read_pickle('synthetic_experiment'),yr)
        except:
            print(gdir.rgi_id,' failed')