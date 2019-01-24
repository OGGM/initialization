from core import *
from plots_paper import *

import os
from functools import partial

import matplotlib.pyplot as plt
import salem
from oggm import cfg, workflow,utils
pd.options.mode.chained_assignment = None
from matplotlib.colors import LinearSegmentedColormap
import scipy

def objective(model1, model2):
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
    objective=0
    for i in range(len(model1.fls)):
        objective = objective + np.sum(
            abs(fls1[i].surface_h - fls2[i].surface_h)**2) + \
                    np.sum(abs(fls1[i].widths - fls2[i].widths)**2)

    return objective


if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER = True

    # Local paths
    if ON_CLUSTER:
        WORKING_DIR = os.environ.get("S_WORKDIR")
        rgidf = salem.read_shapefile(os.path.join('/home/users/julia/reconstruction','rgi','oetztal.shp'))

    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/paper'
        utils.mkdir(WORKING_DIR, reset=False)
        rgidf = salem.read_shapefile(os.path.join(cfg.PATHS['working_dir'],'rgi','oetztal.shp'))


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

    workflow.execute_entity_task(tasks.glacier_masks, gdirs)
    prepare_for_initializing(gdirs)
    synthetic_experiments_parallel(gdirs[0:len(gdirs):4])

    years = np.arange(1850, 1885, 5)
    #years = [1850]

    rel_error_df = pd.DataFrame()
    abs_error_df = pd.DataFrame()

    for gdir in gdirs:

        df_list = {}

        if os.path.isfile(os.path.join(gdir.dir, 'model_run_experiment.nc')):

            for yr in years:
                print(gdir.rgi_id, yr)
                find_possible_glaciers(gdir,yr, 200)
                path = os.path.join(gdir.dir, 'result' + str(yr) + '.pkl')

                if os.path.isfile(path) and not pd.read_pickle(path).empty:
                    df = pd.read_pickle(path)

                else:
                    df = get_single_results(gdir, yr)
                    df.to_pickle(path)
                if df.empty:
                    continue

                rp = gdir.get_filepath('model_run', filesuffix='_experiment')
                ex_mod = FileModel(rp)


                df['objective'] = df.model.apply(objective, model2=ex_mod)
                df.to_pickle(path)

                df['volume'] = df['model'].apply(lambda x: x.volume_m3)
                df['temp_bias'] = df['temp_bias'].apply(lambda x: float(x))

                df_list[str(yr)]=df
                min_mod = deepcopy(df.loc[df.objective.idxmin(), 'model'])
                if yr==1850:
                    abs_error_df.loc[gdir.rgi_id, '1850_min'] = ex_mod.volume_km3_ts()[yr] - min_mod.volume_km3_ts()[yr]
                    rel_error_df.loc[gdir.rgi_id, '1850_min'] = np.log(min_mod.volume_km3_ts()[yr] / ex_mod.volume_km3_ts()[yr])
                    plot_compare_fitness(gdir, df, ex_mod, yr, cfg.PATHS['plot_dir'])

                try:
                    plot_experiment(gdir, ex_mod, yr, cfg.PATHS['plot_dir'])
                    plot_candidates(gdir, df, ex_mod, yr, 'step3', cfg.PATHS['plot_dir'])
                    plot_col_fitness(gdir, df, ex_mod, yr, cfg.PATHS['plot_dir'])
                    m_mod = plot_median(gdir, df, ex_mod, yr, cfg.PATHS['plot_dir'])

                except:
                    m_mod = df.loc[df['objective'].idxmin(), 'model']


                abs_error_df.loc[gdir.rgi_id,yr] = ex_mod.volume_km3_ts()[yr] -m_mod.volume_km3_ts()[yr]
                rel_error_df.loc[gdir.rgi_id, yr] = np.log( m_mod.volume_km3_ts()[yr]/ ex_mod.volume_km3_ts()[yr])

            plot_dir=os.path.join(cfg.PATHS['plot_dir'],'starting')
            utils.mkdir(plot_dir,reset=False)
            plot_fitness_over_time2(gdir,df_list,ex_mod,plot_dir)
            plt.close()


        else:
            print(gdir.rgi_id,' has no experiment')


    abs_error_df.to_pickle(os.path.join(WORKING_DIR, 'abs_error0.pkl'))
    rel_error_df.to_pickle(os.path.join(WORKING_DIR, 'rel_error0.pkl'))


# median_df = pd.read_pickle(os.path.join(WORKING_DIR, 'median.pkl'))
#plot_abs_error_t0([abs_error_df], r'Absoulte error in $t_0$',
#                  os.path.join(cfg.PATHS['plot_dir'], 'abs_error.png'))
#rel_error_df = pd.read_pickle(os.path.join(WORKING_DIR, 'rel_error.pkl'))
#plot_abs_error_t0([rel_error_df.dropna(axis=0)], r'Relative error in $t_0$',
#                  os.path.join(cfg.PATHS['plot_dir'], 'rel_error.png'))
# median_oe_df = pd.read_pickle('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/oetztal/median.pkl')

# median_df = pd.read_pickle('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/rofental/median.pkl')
# plot_median_vs_min([median_df,median_oe_df], cfg.PATHS['plot_dir'])
