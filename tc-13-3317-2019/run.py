

import os
import sys
from copy import deepcopy
sys.path.append('../')
from initialization.core import *
from plots_paper import *

import matplotlib.pyplot as plt
import time
import geopandas as gpd
from oggm import cfg, workflow, utils
pd.options.mode.chained_assignment = None


def fitness_function(model1, model2):
    """
    calculates the objective value (difference in geometry)
    :param model1: oggm.flowline.FluxBasedModel
    :param model2: oggm.flowline.FluxBasedModel
    :return:       float
    """
    model1 = deepcopy(model1)
    model2 = deepcopy(model2)
    model2.run_until(2000)
    model1.run_until(2000)

    fls1 = model1.fls
    fls2 = model2.fls

    fitness = 0
    for i in range(len(model1.fls)):
        fitness = fitness + np.sum(
            abs(fls1[i].surface_h - fls2[i].surface_h)**2) + \
                    np.sum(abs(fls1[i].widths - fls2[i].widths)**2)

    return fitness


def find_median(df, epsilon):

    try:
        accept_df = df[df.fitness <= epsilon]
        quant_df = accept_df[accept_df.fitness <= accept_df.fitness.quantile(0.05)]
        # median state
        quant_df.loc[:, 'length'] = quant_df.model.apply(lambda x: x.length_m)
        quant_df = quant_df.sort_values('length', ascending=False)
        l = len(quant_df)
        if l % 2:
            index = int((l - 1) / 2)
        else:
            index = int(l / 2)
        return deepcopy(quant_df.iloc[index].model)

    except:

        return deepcopy(df.iloc[df.fitness.idxmin()].model)


if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        job_nr = int(os.environ.get('I'))
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/alps/'
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)

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

    # We use intersects
    db = utils.get_rgi_intersects_region_file(version='61', region='11')
    cfg.set_intersects_db(db)

    cfg.PARAMS['run_mb_calibration'] = True
    cfg.PARAMS['optimize_inversion_params'] = False

    # add to BASENAMES
    _doc = 'contains observed and searched glacier from synthetic experiment to find intial state'
    cfg.BASENAMES['synthetic_experiment'] = ('synthetic_experiment.pkl', _doc)

    # RGI file
    path = utils.get_rgi_region_file('11', version='61')
    rgidf = gpd.read_file(path)
    rgidf = rgidf[rgidf.RGIId == 'RGI60-11.00779']

    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=False)

    gdirs = workflow.init_glacier_regions(rgidf)

    if ON_CLUSTER:
        gdirs = gdirs[job_nr:len(gdirs):80]

    # preprocessing(gdirs)

    # experiments
    # synthetic_experiments_parallel(gdirs)

    t_0 = 1850
    t_e = 2000
    epsilon = 125

    # DataFrames for absolute and relative errors
    model_df = pd.DataFrame()
    time_df = pd.DataFrame()

    for gdir in gdirs:
        print(gdir.rgi_id)

        if os.path.isfile(os.path.join(gdir.dir, 'model_run_experiment.nc')):

            start = time.time()
            try:
                rp = gdir.get_filepath('model_run', filesuffix='_experiment')
                ex_mod = FileModel(rp)

                if ex_mod.area_km2_ts()[t_e] > 0.01:

                    df = find_possible_glaciers(gdir, t_0, t_e, 200)

                    df['volume'] = df.model.apply(lambda x: x.volume_km3)
                    df['length'] = df.model.apply(lambda x: x.length_m)

                    # median state
                    med_mod = find_median(df, epsilon)

                    # state with minimal fitness value (based on geometry)
                    min_mod = deepcopy(df.loc[df.fitness.idxmin(), 'model'])

                    # different fitness function values
                    df['fitness2'] = df.model.apply(lambda x: abs(x.area_km2_ts()[t_e] - ex_mod.area_km2_ts()[t_e]) ** 2)
                    df['fitness3'] = df.model.apply(lambda x: abs(x.length_m_ts()[t_e] - ex_mod.length_m_ts()[t_e]) ** 2)

                    # saves median state, minimum state and experiment model for error calculation
                    model_df.loc[gdir.rgi_id, 'median'] = deepcopy(med_mod)
                    model_df.loc[gdir.rgi_id, 'minimum'] = deepcopy(min_mod)
                    model_df.loc[gdir.rgi_id, 'experiment'] = deepcopy(ex_mod)
                    model_df.loc[gdir.rgi_id, 'fit2'] = deepcopy(df.loc[df.fitness2.idxmin(), 'model'])
                    model_df.loc[gdir.rgi_id, 'fit3'] = deepcopy(df.loc[df.fitness3.idxmin(), 'model'])

                    # saves time for each glacier state
                    time_df.loc[gdir.rgi_id, 'time'] = time.time() - start

                    try:
                        # plots for each single glacier
                        plot_experiment(gdir, ex_mod, t_0, cfg.PATHS['plot_dir'])
                        plot_candidates(gdir, df, t_0, 'step3', cfg.PATHS['plot_dir'])
                        plot_fitness_values(gdir, df, ex_mod, t_0, cfg.PATHS['plot_dir'])
                        plot_median(gdir, df, epsilon, ex_mod, t_0, t_e, cfg.PATHS['plot_dir'])
                        animation(gdir, df, ex_mod, med_mod, cfg.PATHS['plot_dir'])
                    except:
                        pass
                    plt.close()

            except:
                pass

        else:
            print(gdir.rgi_id, ' has no experiment')

    # saves model outputs and time measures for each job
    if ON_CLUSTER:
        model_df.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'models_' + str(job_nr) + '.pkl'), compression='gzip')
        time_df.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'time_' + str(job_nr) + '.pkl'), compression='gzip')