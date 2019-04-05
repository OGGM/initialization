
import os
import sys
from copy import deepcopy
sys.path.append('../')
from reconstruction.core import *
from plots_paper import *
from reconstruction.animation import *

import matplotlib.pyplot as plt
import salem
import geopandas as gpd
from oggm import cfg, workflow, utils
pd.options.mode.chained_assignment = None
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import matplotlib as mpl
import time
from scipy import stats
from scipy.optimize import curve_fit
import random

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] =25
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 20 #30
mpl.rcParams['lines.linewidth'] = 3


def add_at(ax, t, loc=2):
    fp = dict(size=20)
    _at = AnchoredText(t, loc=loc, prop=fp ,borderpad=0)
    _at.patch.set_linewidth(3)
    ax.add_artist(_at)
    return _at


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
    rgidf = rgidf[rgidf.RGIId=='RGI60-11.00739']
    # sort for efficient using
    rgidf = rgidf.sort_values('Area', ascending=False)

    gdirs = workflow.init_glacier_regions(rgidf)

    if ON_CLUSTER:
        gdirs = gdirs[job_nr:len(gdirs):80]

    #preprocessing(gdirs)

    # experiments
    #synthetic_experiments_parallel(gdirs)


    t_0 = 1850
    t_e = 2000

    model_df = pd.DataFrame()
    time_df = pd.DataFrame()


    for gdir in gdirs:

        if os.path.isfile(os.path.join(gdir.dir, 'model_run_experiment.nc')) and gdir.rgi_id.endswith('00739'):

            start = time.time()
            #try:
            rp = gdir.get_filepath('model_run', filesuffix='_experiment')
            ex_mod = FileModel(rp)

            if ex_mod.area_km2_ts()[2000] > 0.01:

                df = find_possible_glaciers(gdir, t_0, t_e, 200)
                #df = df[df.fitness < 125]
                df.loc[:, 'length'] = df.model.apply(lambda x: x.length_m)
                df.loc[:, 'volume'] = df.volume/1e9

                med_mod = find_median(df, 125)
                min_mod = deepcopy(df.loc[df.fitness.idxmin(), 'model'])
                animation(gdir, df, ex_mod, med_mod, cfg.PATHS['plot_dir'])
                #plot_median_poster(gdir, df, 125, ex_mod, t_0, t_e, cfg.PATHS['plot_dir'])
                #plot_col_fitness_poster(gdir, df, ex_mod, t_0, cfg.PATHS['plot_dir'])
                #plot_candidates(gdir, df, t_0, 'step3', cfg.PATHS['plot_dir'])
                #plot_candidates(gdir, df, t_0, 'step2', cfg.PATHS['plot_dir'])
                #plot_candidates(gdir, df, t_0, 'step1', cfg.PATHS['plot_dir'])

            #except:
            #    pass

    '''

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

    model_df.loc[:, 'double'] = model_df.index.duplicated(keep='first')
    model_df = model_df[model_df['double'] == False]

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

    error3 = error3.replace([np.inf, -np.inf], np.nan).dropna()
    error4 = error4.replace([np.inf, -np.inf], np.nan).dropna()

    error3 = error3*100
    error4 = error4*100

    # logarithmic errors
    error5 = np.log(median/experiment)
    error6 = np.log(minimum/experiment)

    # errors for different fitness functions
    error7 = (fit2-experiment)/experiment
    error8 = (fit3-experiment)/experiment

    error7 = error7*100
    error8=error8*100



    print(error4.median()[[1850, 1900]])
    print(error7.median()[[1850, 1900]])
    print(error8.median()[[1850, 1900]])

    print(error4.quantile(0.95)[[1850, 1900]])
    print(error7.quantile(0.95)[[1850, 1900]])
    print(error8.quantile(0.95)[[1850, 1900]])


    df_1850 = pd.DataFrame()
    df_1850['Geometry'] = error3.loc[:, 1850]
    df_1850['Area'] = error7.loc[:, 1850]
    df_1850['Length'] = error8.loc[:, 1850]

    df_1850.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'compare_error.pkl'))


    df_1900 = pd.DataFrame()
    df_1900['Geometry'] = error3.loc[:, 1900]
    df_1900['Area'] = error7.loc[:, 1900]
    df_1900['Length'] = error8.loc[:, 1900]


    df_1900.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'compare_error1900.pkl'))

    import seaborn as sns


    df_1850 = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'compare_error.pkl'))
    df_1900 = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'compare_error1900.pkl'))


    df_1850 = (df_1850*100)
    df_1900 = (df_1900*100)


    df_1850 = df_1850[df_1850 < 600].dropna()
    df_1900 = df_1900[df_1900.index.isin(df_1850.index)]
    #df_1900 = df_1900[df_1900 < 600].dropna()


    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(25, 18))

    bins1 = np.arange(-100, 620, 14)  # The edges
    bins2 = np.arange(-100, 620, 14)

    #plt.figure()
    df_1850['Geometry'].hist(bins=bins1, alpha=0.5, color='C0', ax=ax2)
    df_1850['Area'].hist(bins=bins1, alpha=0.5, color='C1', ax=ax2)
    df_1850['Length'].hist(bins=bins1, alpha=0.5, color='C2', ax=ax2)


    sns.kdeplot(df_1850['Geometry'], shade=False, bw=0.1, ax=ax2, label='')
    sns.kdeplot(df_1850['Area'], shade=False, bw=0.1, ax=ax2, label='')
    sns.kdeplot(df_1850['Length'], shade=False, bw=0.1, ax=ax2, label='')



    df_1900['Geometry'].hist(bins=bins2, alpha=0.5, color='C0',
                             ax=ax1)
    df_1900['Area'].hist(bins=bins2, alpha=0.5, color='C1', ax=ax1)
    df_1900['Length'].hist(bins=bins2, alpha=0.5, color='C2',
                           ax=ax1)


    sns.kdeplot(df_1900['Geometry'], shade=False, bw=0.1, ax=ax1, label='')
    sns.kdeplot(df_1900['Area'], shade=False, bw=0.1, ax=ax1, label='')
    sns.kdeplot(df_1900['Length'], shade=False, bw=0.1, ax=ax1, label='')



    #ax1.set_ylim(None,2600)
    #ax2.set_ylim(None, 2600)

    #ax1.set_yticks([0,1000,2000])
    #ax2.set_yticks([0,500,1000])

    ax2.tick_params(axis='both', which='major', labelsize=50)
    ax1.tick_params(axis='both', which='major', labelsize=50)

    ax = f.add_subplot(111, frameon=False)
    ax.set_yticks([0, 500, 1000])
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off',
                    right='off', labelsize=60)
    plt.suptitle('Errors of reconstructed states \n (based on minimal fitness values), n=2317', y=1, fontsize=50)


    plt.xlabel(r'Volume error ($\%$)',fontsize=50)
    plt.ylabel('Frequency', fontsize=50)

    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)

    import matplotlib.patches as mpatches
    p2 = mpatches.Patch(color='C0', alpha=0.5, linewidth=0)
    p4 = mpatches.Patch(color='C1', alpha=0.4, linewidth=0)
    p6 = mpatches.Patch(color='C2', alpha=0.4, linewidth=0)

    leg = plt.legend(((p2), (p4), (p6)), ('Geometry', 'Area', 'Length'),
               fontsize='medium', framealpha=0)

    leg.set_title("Evaluation type", prop={'size': 'medium'})

    #ax.tick_params(axis='both', which='major', labelsize=100)

    plt.savefig(os.path.join(cfg.PATHS['plot_dir'], 'compare_median.png'),
                dpi=300, transparent=False)

    plt.show()


    plt.figure(figsize=(25, 18))

    sns.kdeplot(df_1850['Geometry'], shade=True, bw=0.1)
    sns.kdeplot(df_1850['Area'], shade=True, bw=0.1)
    sns.kdeplot(df_1850['Length'], shade=True, bw=0.1)

    plt.title('Bias of minimum states in 1850, n=2317')
    plt.xlabel(r'Volume error ($\%$)')
    plt.ylabel('Density')
    leg = plt.legend(fontsize='small', framealpha=0)
    leg.set_title("Evaluation type", prop={'size': 'small'})
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'], 'compare_median3.png'),
                dpi=300, transparent=True)

    plt.show()

    # plots
    plot_relative_error(error1, error2, 'abs', cfg.PATHS['plot_dir'], all=True)
    #plot_relative_error(error3, error4, 'rel', cfg.PATHS['plot_dir'], all=False)
    #plot_relative_error(error5, error6, 'log', cfg.PATHS['plot_dir'], all=True)


    import matplotlib.patches as mpatches

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(15, 10),
                                 gridspec_kw={'width_ratios': [1, 3],
                                              'wspace': 0.05})
    plt.subplots_adjust(hspace=.001)
    plt.suptitle(r'Relative error (Alps), n=' + str(len(error4.index)), y=0.99)

    p1, = ax2.plot(error4.columns.values, error4.median().values, zorder=5)
    ax2.fill_between(error4.columns.values, error4.quantile(0.95), error4.quantile(0.05), alpha=0.5, zorder=3)
    p3, = ax2.plot(error7.columns.values, error7.median().values, zorder=4)
    ax2.fill_between(error7.columns.values, error7.quantile(0.95), error7.quantile(0.05), alpha=0.4, zorder=2)
    p5, = ax2.plot(error8.columns.values, error8.median().values, zorder=4)
    ax2.fill_between(error8.columns.values, error8.quantile(0.95), error8.quantile(0.05), alpha=0.4, zorder=1)

    #plt.grid()
    ax2.set_xlabel('Time (years)')
    ax1.set_ylabel(r'Volume error ($\%$) ')
    p2 = mpatches.Patch(color='C0', alpha=0.5, linewidth=0)
    p4 = mpatches.Patch(color='C1', alpha=0.4, linewidth=0)
    p6 = mpatches.Patch(color='C2', alpha=0.4, linewidth=0)



    df1 = error4.loc[:, 1850]
    df1 = df1[(df1 > df1.quantile(0.05)) & (df1 < df1.quantile(0.95))]

    df2 = error7.loc[:, 1850]
    df2 = df2[(df2 > df2.quantile(0.05)) & (df2 < df2.quantile(0.95))]

    df3 = error8.loc[:, 1850]
    df3 = df3[(df3 > df3.quantile(0.05)) & (df3 < df3.quantile(0.95))]

    bins = np.arange(np.min([df1.min(), df2.min(),df3.min()]), np.max([df1.max(), df2.max(), df3.max()]), 30)

    df1.plot.hist(ax=ax1, bins=bins, color='C0', alpha=0.5, orientation='horizontal')
    df2.plot.hist(ax=ax1, bins=bins, color='C1', alpha=0.5,
                  orientation='horizontal')
    df3.plot.hist(ax=ax1, bins=bins, color='C2', alpha=0.5,
                  orientation='horizontal')
    ax1.invert_xaxis()
    ax1.set_xticks([1000,500,0])

    plt.suptitle(r'Relative error (Alps), n=' + str(len(error4.index)), y=0.99)
    ax1.set_title('1850', fontsize=20)
    ax2.set_title('1850-2000', fontsize=20)

    add_at(ax1, r"a", loc=1)
    add_at(ax2, r"b", loc=1)

    ax2.legend(((p1, p2), (p3, p4), (p5, p6)), ('Geometry', 'Area', 'Length'), loc='best', bbox_to_anchor=(0.95, 1))

    ax1.grid()
    ax2.grid()

    #plt.legend(loc='best')
    #ax1.tick_params(axis='both', which='major', labelsize=20)

    plt.savefig(os.path.join(cfg.PATHS['plot_dir'], 'compare_median.png'), dpi=300)
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'], 'compare_median.pdf'), dpi=300)
    plt.show()




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


    time_df = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'time_merge.pkl'))
    # time in hours
    time_df['time'] = time_df.time.apply(lambda x: (x/60/60))
    time_df.loc[:, 'double'] = time_df.index.duplicated(keep='first')

    time_df = time_df[time_df['double'] == False]
    print(time_df.time.min()*60*60)
    print(time_df.time.max()*60)
    print(time_df.time.sum()/24)

    # time per cpu hour
    time_df['time'] = time_df.time.apply(lambda x: x *28)
    print(time_df.time.min())
    print(time_df.time.max())
    print(time_df.time.sum())


    time_df['area'] = rgidf.set_index('RGIId').loc[time_df.index].Area

    fig = plt.figure(figsize=(15, 10))
    ax1 = fig.add_subplot(111)

    time_df.plot.scatter('time', 'area', alpha=0.7, ax=ax1)
    plt.ylabel(r'Glacier area ($km^2$)')
    plt.xlabel(r'Time (CPU hours)')
    plt.title('Reconstructions for the Alps, n=' + str(len(time_df)))

    slope, intercept, r_value, p_value, std_err = stats.linregress(x=time_df.time.values, y=time_df.area.values)

    x = np.linspace(time_df.time.min(), time_df.time.max(), 100)

    popt, pcov = curve_fit(func, time_df.time.values, time_df.area.values)
    slope, intercept, r_value, p_value, std_err = stats.linregress(time_df.time.values, time_df.area.values)

    print(r_value)

    y = func(x, popt[0], popt[1])
    y1 = func(x, popt[0] + pcov[0, 0] ** 0.5, popt[1] - pcov[1, 1] ** 0.5)
    y2 = func(x, popt[0] - pcov[0, 0] ** 0.5, popt[1] + pcov[1, 1] ** 0.5)

    ax1.plot(x, y, color='C0')
    ax1.fill_between(x, y1, y2, facecolor='C0', alpha=0.2)
    #plt.savefig(os.path.join(cfg.PATHS['plot_dir'], 'time', 'alps.png'), dpi=300)
    plt.close()


    popt2, pcov2 = curve_fit(func, time_df.area.values, time_df.time.values)

    df = pd.DataFrame()
    fig = plt.figure(figsize=(20, 10))
    ax2 = fig.add_subplot(111)

    for i in range(1, 20):
        path = utils.get_rgi_region_file(str(i).zfill(2), version='61')
        rgidf = gpd.read_file(path)
        region = path.split('.shp')[0].split('_')[-1]
        print(region)
        time = rgidf.Area.apply(func, a=popt2[0], b=popt2[1]).sum()
        df.loc[region, 'Time (days)'] = (time/24)

    print('Jahre:'+str(df['Time (days)'].sum()/365))

    #df =pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'estimated_time_per_region.pkl'))
    df['Time (days)'] = df['Time (days)'].apply(lambda x: (x*24)/1e4)
    print(df)
    print('10^4 stunden'+str(df['Time (days)'].sum()))

    df.plot.barh(legend=False, ax=ax2, fontsize=20)
    plt.xlabel(r'Time ($10^4$ CPU hours) ', fontsize=20)
    plt.title('Estimated running time per region', fontsize=25)
    plt.tight_layout()
    plt.savefig(os.path.join(cfg.PATHS['plot_dir'], 'time', 'per_region.png'), dpi=300)

    #df.to_pickle(os.path.join(cfg.PATHS['working_dir'],'estimated_time_per_region.pkl'))

    area = model_df['experiment'].apply(lambda x: x.area_km2_ts()[2000])
    rgidf = rgidf.set_index('RGIId')
    diff = pd.DataFrame()
    model_df.loc[:, 'double'] = model_df.index.duplicated(keep = False)
    area = area[model_df['double']==False]
    diff = pd.DataFrame()

    print(area.sum())
    print(rgidf[rgidf.index.isin(area.index)].Area.sum())


    for i in area.index:

        diff.loc[i, 'difference'] = area.loc[i]-rgidf.loc[i].Area
    diff.to_pickle(os.path.join(cfg.PATHS['working_dir'], 'difference.pkl'))

    diff = pd.read_pickle(os.path.join(cfg.PATHS['working_dir'], 'difference.pkl'))
    diff.plot.hist(bins=50)
    plt.xlabel(r'Difference in Area (km^2)')
    plt.show()

    print(diff.describe())

    '''

