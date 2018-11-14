from core import *
from plots_paper import *

import os
from copy import deepcopy
from functools import partial
import pickle

import matplotlib.pyplot as plt
import xarray as xr
import geopandas as gpd
import shapely.geometry as shpg
from scipy.optimize import curve_fit

import salem
from oggm import cfg, workflow, tasks, utils, graphics
from oggm.utils import get_demo_file
from oggm.core.flowline import FluxBasedModel, FileModel
FlowlineModel = partial(FluxBasedModel, inplace=False)




if __name__ == '__main__':
    cfg.initialize()

    ON_CLUSTER_NEW = False
    ON_CLUSTER_OLD = True

    # Local paths
    if ON_CLUSTER_NEW:
        WORKING_DIR = os.environ.get("S_WORKDIR")
    elif ON_CLUSTER_OLD:
        WORKING_DIR = 'out/rofental2'
    else:
        WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/rofental2'
        #WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction'
        #WORKING_DIR = '/home/juliaeis/Dokumente/OGGM/work_dir/find_initial_state/past_state_information'
        utils.mkdir(WORKING_DIR, reset=False)
    cfg.PATHS['working_dir'] = WORKING_DIR

    cfg.PATHS['plot_dir'] = os.path.join(cfg.PATHS['working_dir'], 'plots')
    #cfg.PATHS['plot_dir'] = '/home/juliaeis/Dokumente/eigene_paper/reconstruction/plots'
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
    #prepare_for_initializing(gdirs)
    #synthetic_experiments_parallel(gdir)


    #years = [1850,1875,1900,1925,1950]
    years = np.arange(1850, 1970,5)
    #years = [1850,1900]
    volumes = pd.DataFrame()
    median_df = pd.DataFrame()
    #ranged = pd.DataFrame(columns=years)

    for gdir in gdirs:
        df_list = {}
        to_range = []

        if os.path.isfile(os.path.join(gdir.dir,'synthetic_experiment.pkl')): #and gdir.rgi_id.endswith('897'):

            
            for yr in years:
                print(gdir.rgi_id,yr)

                #find_possible_glaciers(gdir,gdir.read_pickle('synthetic_experiment'),yr)
                path = os.path.join(gdir.dir, 'result' + str(yr) + '.pkl')

                if os.path.isfile(path) and not pd.read_pickle(path).empty:
                    df = pd.read_pickle(path)

                else:
                    df = get_single_results(gdir, yr, gdir.read_pickle(
                        'synthetic_experiment'))
                    df.to_pickle(path)

                rp = gdir.get_filepath('model_run', filesuffix='experiment')
                ex_mod = FileModel(rp)


                df['volume'] = df['model'].apply(lambda x: x.volume_m3)
                df['temp_bias'] = df['temp_bias'].apply(lambda x: float(x))

                df_list[str(yr)]=df

                #plot_experiment(gdir, df, ex_mod, yr, cfg.PATHS['plot_dir'])
                #plot_compare_fitness(gdir, df, ex_mod, yr, cfg.PATHS['plot_dir'])
                #plot_candidates(gdir, df, ex_mod, yr, 'step3',cfg.PATHS['plot_dir'])
                #plot_col_fitness(gdir, df, ex_mod, yr, cfg.PATHS['plot_dir'])

                try:
                    m_mod = plot_median(gdir, df, ex_mod, yr, cfg.PATHS['plot_dir'])

                except:
                    m_mod = df.loc[df['objective'].idxmin(), 'model']


                #median_df = median_df.append({'rgi': gdir.rgi_id, 'm_mod':m_mod,'ex_p':rp, 'min_mod':df.loc[df['objective'].idxmin(),'model']}, ignore_index=True)
                ex_mod.run_until(yr)
                m_mod.reset_y0(yr)

                median_df.loc[gdir.rgi_id,yr] = m_mod.volume_km3 - ex_mod.volume_km3
                plt.close()
                '''

                max_model = deepcopy(df.loc[df.volume.idxmax(), 'model'])

                max_obj = df.loc[df.volume.idxmax(), 'objective']

                rp = gdir.get_filepath('model_run', filesuffix='experiment')
                ex_mod = FileModel(rp)
                v1850 = deepcopy(ex_mod.volume_m3)
                print(gdir.rgi_id,v1850)
                ex_mod.run_until(2000)
                v2000 = ex_mod.volume_m3
                volumes = volumes.append(
                    {'rgi': gdir.rgi_id, 'ratio(experiment)': v1850 / v2000,
                     'ratio(max)': df.volume.max() / v2000,
                     'objective(max)': max_obj, 'temp_bias': df.temp_bias.min()},
                    ignore_index=True)
                r = df[df.objective<=100].volume.max()-df[df.objective<=100].volume.min()
                to_range.append(r)

                '''
            #plot_dir=os.path.join(cfg.PATHS['plot_dir'],'starting')
            #utils.mkdir(plot_dir,reset=False)
                
            #range.loc[gdir.rgi_id,:] = to_range

            #plot_fitness_over_time2(gdir,df_list,ex_mod,plot_dir)
            #plt.show()

        else:
            print(gdir.rgi_id,' has no experiment')


#range.to_pickle(os.path.join(WORKING_DIR,'range.pkl'))

#plot_range(rgidf,os.path.join(WORKING_DIR,'range.pkl'),cfg.PATHS['plot_dir'],False)

median_df.to_pickle(os.path.join(WORKING_DIR, 'median.pkl'))
median_df = pd.read_pickle(os.path.join(WORKING_DIR, 'median.pkl'))
print(median_df)
#median_oe_df = pd.read_pickle('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/oetztal/median.pkl')

#median_df = pd.read_pickle('/home/juliaeis/Dokumente/OGGM/work_dir/reconstruction/rofental/median.pkl')
#plot_median_vs_min([median_df,median_oe_df], cfg.PATHS['plot_dir'])
