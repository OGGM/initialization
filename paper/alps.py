
import os
import sys
from copy import deepcopy
sys.path.append('../')
from reconstruction.core import *
from plots_paper import *

import matplotlib.pyplot as plt
import salem
import geopandas as gpd
from oggm import cfg, workflow, utils
pd.options.mode.chained_assignment = None
import time
from scipy import stats
from scipy.optimize import curve_fit

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

def func(x, a, b):
    return(x*a+b)

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
    '''
    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=False)

    gdirs = workflow.init_glacier_regions(rgidf)
    if ON_CLUSTER:
        gdirs = gdirs[job_nr:len(gdirs):80]

    preprocessing(gdirs)

    # experiments
    synthetic_experiments_parallel(gdirs)

    t_0 = 1850
    t_e = 2000

    model_df = pd.DataFrame()
    time_df = pd.DataFrame()


    for gdir in gdirs:

        if os.path.isfile(os.path.join(gdir.dir, 'model_run_experiment.nc')):
            start = time.time()
            try:
                rp = gdir.get_filepath('model_run', filesuffix='_experiment')
                ex_mod = FileModel(rp)

                if ex_mod.area_km2_ts()[2000] > 0.01:

                        df = find_possible_glaciers(gdir, t_0, t_e, 20)
                        med_mod = find_median(df, 125)
                        min_mod = deepcopy(df.loc[df.fitness.idxmin(), 'model'])

                        df['fitness2'] = df.model.apply(lambda x: abs(x.area_km2_ts()[2000] - ex_mod.area_km2_ts()[2000]) ** 2)
                        df['fitness3'] = df.model.apply(lambda x: abs(x.length_m_ts()[2000] - ex_mod.length_m_ts()[2000]) ** 2)

                        # saves median state, minimum state and experiment model
                        model_df.loc[gdir.rgi_id, 'median'] = deepcopy(med_mod)
                        model_df.loc[gdir.rgi_id, 'minimum'] = deepcopy(min_mod)
                        model_df.loc[gdir.rgi_id, 'experiment'] = deepcopy(ex_mod)
                        model_df.loc[gdir.rgi_id, 'fit2'] = deepcopy(df.loc[df.fitness2.idxmin(), 'model'])
                        model_df.loc[gdir.rgi_id, 'fit3'] = deepcopy(df.loc[df.fitness3.idxmin(), 'model'])

                        # time_df
                        time_df.loc[gdir.rgi_id, 'time'] = time.time()-start

                        try:
                            # plots
                            plot_candidates(gdir, df, t_0, 'step3', cfg.PATHS['plot_dir'])
                            plot_col_fitness_average(gdir, df, ex_mod, t_0, cfg.PATHS['plot_dir'])
                            plot_median_average(gdir, df, 125, ex_mod, t_0, t_e, cfg.PATHS['plot_dir'])
                            plot_compare_fitness(gdir, df, ex_mod, t_0, cfg.PATHS['plot_dir'])

                        except:
                            pass

            except:
                print(gdir.rgi_id+' failed')

    if ON_CLUSTER:
        model_df.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'models_'+str(job_nr)+'.pkl'), compression='gzip')
        time_df.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'time_'+str(job_nr)+'.pkl'), compression='gzip')


    model_df = pd.DataFrame()
    time_df = pd.DataFrame()
    for file in os.listdir(cfg.PATHS['working_dir']):
        print(file)
        if file.startswith('models'):
            p = os.path.join(cfg.PATHS['working_dir'], file)
            model_df = model_df.append(pd.read_pickle(p, compression='gzip'))
        if file.startswith('time'):
            p = os.path.join(cfg.PATHS['working_dir'], file)
            time_df = time_df.append(pd.read_pickle(p, compression='gzip'))

    model_df.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'models_merge.pkl'))
    time_df.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'time_merge.pkl'))



    model_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'models_merge.pkl'))
    median = model_df['median'].apply(lambda x: x.volume_km3_ts())
    minimum = model_df['minimum'].apply(lambda x: x.volume_km3_ts())
    experiment = model_df['experiment'].apply(lambda x: x.volume_km3_ts())
    fit2 = model_df['fit2'].apply(lambda x: x.volume_km3_ts())
    fit3 = model_df['fit3'].apply(lambda x: x.volume_km3_ts())

    # absolute errors
    error1 = median-experiment
    error2 = minimum-experiment

    # relative errors
    error3 = (median-experiment)/experiment
    error4 = (minimum-experiment)/experiment

    # logarithmic errors
    error5 = np.log(median/experiment)
    error6 = np.log(minimum/experiment)

    # errors for different fitness functions
    error7 = (fit2-experiment)/experiment
    error8 = (fit3-experiment)/experiment

    # plots
    plot_relative_error(error1, error2, 'abs', cfg.PATHS['plot_dir'], all=True)
    plot_relative_error(error3, error4, 'rel', cfg.PATHS['plot_dir'], all=True)
    plot_relative_error(error5, error6, 'log', cfg.PATHS['plot_dir'], all=True)

    plt.figure(figsize=(15, 10))
    plt.title('Relative error, n=' + str(len(error4.index)))
    error4.median().plot(label='geometry')
    plt.fill_between(error4.columns.values, error4.quantile(0.75), error4.quantile(0.25), alpha=0.5, zorder=3)
    error7.median().plot(label='area')
    plt.fill_between(error7.columns.values, error7.quantile(0.75), error7.quantile(0.25), alpha=0.5, zorder=2)
    error8.median().plot(label='length')
    plt.fill_between(error8.columns.values, error8.quantile(0.75), error8.quantile(0.25), alpha=0.5, zorder=1)

    plt.grid()
    plt.xlabel('Time')
    plt.ylabel(r'Median with interquartile range ')
    plt.legend(loc='best')
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'], 'errors', 'compare_median.png'),dpi=300)

    # mean
    plt.figure(figsize=(15, 10))
    plt.title('Relative error, n=' + str(len(error4.index)))
    error4.mean().plot(label='geometry')
    plt.fill_between(error4.columns.values, error4.mean()+error4.std(),
                     error4.mean() - error4.std(), alpha=0.5, zorder=3)
    error7.median().plot(label='area')
    plt.fill_between(error7.columns.values, error7.mean()+error7.std(),
                     error7.mean() - error7.std(), alpha=0.5, zorder=2)
    error8.median().plot(label='length')
    plt.fill_between(error8.columns.values, error8.mean()+error8.std(),
                     error8.mean() - error8.std(), alpha=0.5, zorder=1)

    plt.grid()
    plt.xlabel('Time')
    plt.ylabel(r'Mean with standard deviation ')
    plt.legend(loc='best')
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'], 'errors', 'compare_mean.png'), dpi=300)
    plt.show()

    '''

    time_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'time_merge.pkl'))
    time_df['time'] = time_df.time.apply(lambda x: x/60)
    time_df['area'] = rgidf.set_index('RGIId').loc[time_df.index].Area
    time_df.plot.scatter('time', 'area', alpha=0.7)
    plt.ylabel(r'Area ($km^2$)')
    plt.xlabel(r'Time (minutes)')

    slope, intercept, r_value, p_value, std_err = stats.linregress(x=time_df.time.values, y=time_df.area.values)
    print(slope, intercept)
    x = np.linspace(time_df.time.min(), 40, 100)
    y = slope*x + intercept
    plt.plot(x, y, color='C0')

    plt.show()

    popt, pcov = curve_fit(func, time_df.time.values, time_df.area.values)
    print(popt, pcov)