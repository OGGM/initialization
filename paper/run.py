

import os
import sys
from copy import deepcopy
sys.path.append('../')
from reconstruction.core import *
from plots_paper import *

import matplotlib.pyplot as plt
import salem
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




if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = False

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        cfg.PATHS['working_dir'] = WORKING_DIR
        rgidf = salem.read_shapefile(os.path.join('/home/users/julia/reconstruction','rgi', 'oetztal.shp'))
        job_nr = int(os.environ.get('I')) - 1
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/paper'
        cfg.PATHS['working_dir'] = WORKING_DIR
        utils.mkdir(WORKING_DIR, reset=False)
        rgidf = salem.read_shapefile('../rgi/oetztal.shp')

    cfg.PATHS['plot_dir'] = os.path.join(cfg.PATHS['working_dir'], 'plots2')
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

    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=False)

    gdirs = workflow.init_glacier_regions(rgidf)
    #preprocessing(gdirs)

    if ON_CLUSTER:
        gdirs = gdirs[job_nr:len(gdirs):20]

    # experiments
    #synthetic_experiments_parallel(gdirs)

    years = np.arange(1850, 1970, 10)
    #years = [1850]
    # DataFrames for absolute and relative errors
    #rel_df = pd.DataFrame()
    #abs_df = pd.DataFrame()
    #ids = [os.path.join(cfg.PATHS['plot_dir'], 'starting', 'starting_'+id+'.png') for id in ['RGI60-11.00897', 'RGI60-11.00779', 'RGI60-11.00739']]
    #union_plot(ids, cfg.PATHS['plot_dir'])

    abs_df = pd.read_pickle(os.path.join(WORKING_DIR, 'abs_error_new.pkl'))
    rel_df = pd.read_pickle(os.path.join(WORKING_DIR, 'rel_error_new.pkl'))
    '''
    for gdir in gdirs:
        print(gdir.rgi_id)
        df_list = {}

        if os.path.isfile(os.path.join(gdir.dir, 'model_run_experiment.nc')) :#and (gdir.rgi_id in ['RGI60-11.00897', 'RGI60-11.00779', 'RGI60-11.00739']):

            for yr in years:
                print(yr)
                df = find_possible_glaciers(gdir, yr, 2000, 200)
                #df = pd.read_pickle(os.path.join(gdir.dir,'result'+str(yr)+'.pkl'), compression='gzip')

                df['volume'] = df.model.apply(lambda x: x.volume_m3)
                df_list[str(yr)] = df

                rp = gdir.get_filepath('model_run', filesuffix='_experiment')
                ex_mod = FileModel(rp)

                abs_df.loc[gdir.rgi_id, 'exp_1850'] = ex_mod.volume_km3_ts()[yr]
                rel_df.loc[gdir.rgi_id, 'exp_1850'] = ex_mod.volume_km3_ts()[yr]


                x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * \
                    ex_mod.fls[-1].map_dx

                df['fitness'] = df.model.apply(lambda x: fitness_function(x, ex_mod))


                if yr == 1850:
                    #plot_compare_fitness(gdir, df, ex_mod, yr, cfg.PATHS['plot_dir'])

                    # errors
                    min_mod = deepcopy(df.loc[df.objective.idxmin(), 'model'])
                    v_min = min_mod.volume_km3_ts()[yr]
                    v_ex = ex_mod.volume_km3_ts()[yr]
                    abs_df.loc[gdir.rgi_id, '1850_min'] = v_ex - v_min
                    rel_df.loc[gdir.rgi_id, '1850_min'] = np.log(v_min / v_ex)
                try:
                    pass
                    #plot_experiment(gdir, ex_mod, yr, cfg.PATHS['plot_dir'])
                    #plot_candidates(gdir, df, yr, 'step3', cfg.PATHS['plot_dir'])
                    #plot_col_fitness(gdir, df, ex_mod, yr, cfg.PATHS['plot_dir'])
                except:
                    pass


                try:
                    m_mod = plot_median(gdir, df, ex_mod, yr,
                                        cfg.PATHS['plot_dir'])
                # if not, take min state
                except:
                    m_mod = deepcopy(df.loc[df.objective.idxmin(),
                                          'model'])

                # errors
                v_ex = ex_mod.volume_km3_ts()[yr]
                v_med = m_mod.volume_km3_ts()[yr]
                abs_df.loc[gdir.rgi_id,yr] = v_ex - v_med
                rel_df.loc[gdir.rgi_id, yr] = np.log(v_med/v_ex)



            plot_dir = os.path.join(cfg.PATHS['plot_dir'], 'starting')
            utils.mkdir(plot_dir, reset=False)
            #plot_fitness_over_time(gdir, df_list, ex_mod, plot_dir)
            plt.close()
        else:
            print(gdir.rgi_id, ' has no experiment')

    if ON_CLUSTER:
        abs_df.to_pickle(os.path.join(WORKING_DIR, 'abs_error'+str(job_nr)+'.pkl'))
        rel_df.to_pickle(os.path.join(WORKING_DIR, 'rel_error'+str(job_nr)+'.pkl'))
    else:
        abs_df.to_pickle(os.path.join(WORKING_DIR, 'abs_error_new.pkl'))
        rel_df.to_pickle(os.path.join(WORKING_DIR, 'rel_error_new.pkl'))

    abs_df = pd.read_pickle(os.path.join(WORKING_DIR, 'abs_error_new.pkl'))
    rel_df = pd.read_pickle(os.path.join(WORKING_DIR, 'rel_error_new.pkl'))
    '''
    #print(rel_df)
    rel_df = rel_df.replace([np.inf, -np.inf], np.nan).dropna().drop(columns=['1850_min', 'exp_1850'])
    print(rel_df.describe())
    print(abs_df.describe())
    print(abs_df.loc[abs_df[1850].idxmin(),:] )
    #plot_error_t02([rel_df, pd.DataFrame()], cfg.PATHS['plot_dir'],abs=False)

    #plot_error_t02([abs_df.drop(columns=['1850_min', 'exp_1850']), pd.DataFrame()], cfg.PATHS['plot_dir'], abs=True)
    #plt.show()

    #plot_min_vs_med([rel_df[[1850, '1850_min']], pd.DataFrame()],'', cfg.PATHS['plot_dir'])
    '''

    rel_df = pd.DataFrame()
    abs_df = pd.DataFrame()
    for error in os.listdir(os.path.join(WORKING_DIR, 'errors2')):
        if error.startswith('rel'):
            rel_df = rel_df.append(pd.read_pickle(os.path.join(WORKING_DIR, 'errors2',error)))
        elif error.startswith('abs'):
            abs_df = abs_df.append(pd.read_pickle(os.path.join(WORKING_DIR, 'errors2', error)))

    print(rel_df[[1850, '1850_min']])
    plot_min_vs_med([abs_df[[1850, '1850_min']].drop('RGI60-11.00746'), pd.DataFrame()], '', cfg.PATHS['plot_dir'])
    '''