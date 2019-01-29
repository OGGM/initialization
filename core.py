import os
import copy
from functools import partial

# External libs
import numpy as np
import pandas as pd
import multiprocessing as mp
from multiprocessing import Pool
from scipy.signal import argrelextrema

# locals
from oggm import workflow, tasks, utils
from oggm.core.inversion import mass_conservation_inversion
from oggm.core.flowline import FileModel


def _find_extrema(ts):
    """
    Needed to determine t_stag. Trajectories will be smoothed and
    extrema will be determined
    """
    # smooth to find maximum
    smooth_ts = ts.rolling(10).mean()
    # fill nan with 0, to avoid warning from np.greater
    smooth_ts = smooth_ts.fillna(0)
    # find extrema and take the first one
    extrema = argrelextrema(smooth_ts.values, np.greater, order=20)[0][0]
    return extrema


def _read_file_model(suffix, gdir):
    """
    Create FileModel from gdir and suffix
    """
    rp = gdir.get_filepath('model_run', filesuffix=suffix)
    fmod = FileModel(rp)
    return copy.deepcopy(fmod)


def _run_parallel_experiment(gdir):
    """
    Creates the synthetic experiment for one glacier. model_run_experiment.nc
    will be saved in working directory.
    """

    try:
        fls = gdir.read_pickle('model_flowlines')
        # try to run random climate with temperature bias -1
        try:
            model = tasks.run_random_climate(gdir, nyears=400, bias=0, seed=1,
                                             temperature_bias=-1,
                                             init_model_fls=fls)

        # perhaps temperature_bias -1 was to ambitious, try larger one
        except:
            model = tasks.run_random_climate(gdir, nyears=400, bias=0, seed=1,
                                             temperature_bias=-0.5,
                                             init_model_fls=fls)

        # construct observed glacier, previous glacier will be run forward from
        # 1850 - 2000 with past climate file
        fls = copy.deepcopy(model.fls)
        tasks.run_from_climate_data(gdir, ys=1850, ye=2000, init_model_fls=fls,
                                    output_filesuffix='_experiment')
    except:
        print('experiment failed : ' + str(gdir.rgi_id))


def _run_to_present(tupel, gdir, ys, ye):
    """
    Run glacier candidates forwards.
    """
    suffix = tupel[0]
    path = gdir.get_filepath('model_run', filesuffix=suffix)
    # does file already exists?
    if not os.path.exists(path):
        try:
            tasks.run_from_climate_data(gdir, ys=ys, ye=ye,
                                        output_filesuffix=suffix,
                                        init_model_fls=copy.deepcopy(
                                        tupel[1].fls))
            return suffix
        # oggm failed --> probaly "glacier exeeds boundaries"
        except:
            return None

    else:
        # does file contain a model?
        try:
            fmod = FileModel(path)
            return suffix
        except:
            return None


def _run_random_parallel(gdir, y0, list):
    """
    Paralleize the run_random_task.
    """
    pool = Pool()
    paths = pool.map(partial(_run_random_task, gdir=gdir, y0=y0), list)
    pool.close()
    pool.join()

    random_run_list = pd.DataFrame()
    for rp in paths:
        if rp != None:
            temp_bias = rp.split('.nc')[0].split('_')[-1]
            seed = rp.split('.nc')[0].split('_')[-2]
            suffix = str(y0)+'_random_' + str(seed)+'_' + str(temp_bias)
            v = pd.Series({'seed': seed, 'temp_bias': float(temp_bias),
                           'suffix': suffix})
            random_run_list = random_run_list.append(v, ignore_index=True)
    return random_run_list


def _run_random_task(tupel, gdir, y0):
    """
    Run random model to create lots of possible states
    """
    seed = tupel[0]
    temp_bias = tupel[1]
    fls = gdir.read_pickle('model_flowlines')
    suffix = str(y0) + '_random_'+str(seed) + '_' + str(temp_bias)

    # test if file already exist:
    path = gdir.get_filepath('model_run', filesuffix=suffix)

    # does file already exists?
    if not os.path.exists(path):
        try:
            tasks.run_random_climate(gdir, nyears=400, y0=y0, bias=0,
                                     seed=seed, temperature_bias=temp_bias,
                                     init_model_fls=copy.deepcopy(fls),
                                     output_filesuffix=suffix)
            return path
        # oggm failed --> probaly "glacier exeeds boundaries"
        except:
            return None

    else:
        # does file contain a model?
        try:
            fmod = FileModel(path)
            return path
        except:
            return None


def _run_file_model(suffix, gdir, ye):
    """
    Read FileModel and run it until ye
    """
    rp = gdir.get_filepath('model_run', filesuffix=suffix)
    fmod = FileModel(rp)
    fmod.run_until(ye)
    return copy.deepcopy(fmod)


def find_candidates(gdir, df, ys, ye, n):
    """
    Determine glacier candidates and run them to the date of observation
    :param gdir:    oggm.GlacierDirectories
    :param df:      pd.DataFrame (volume_m3_ts() from random climate runs)
    :param ys:      starting year
    :param ye:      year of observation
    :param n:       number of candidates
    :return:
    """
    indices = []
    # find nearest glacier state for each of the n volume classes (equidistant)
    for val in np.linspace(df.ts_section.min(), df.ts_section.max(), n):
        index = df.iloc[(df['ts_section'] - val).abs().argsort()][:1].index[0]
        if not index in indices:
            indices = np.append(indices, index)
    candidates = df.ix[indices]
    candidates = candidates.sort_values(['suffix', 'time'])
    candidates['fls_t0'] = None
    for suffix in candidates['suffix'].unique():
        rp = gdir.get_filepath('model_run', filesuffix=suffix)
        fmod = FileModel(rp)
        for i, t in candidates[candidates['suffix'] == suffix]['time'].iteritems():
            fmod.run_until(t)
            candidates.at[i, 'random_model_t0'] = copy.deepcopy(fmod)

    candidates = candidates.drop_duplicates()
    fls_list = []
    for i in candidates.index:
        s = candidates.loc[int(i), 'suffix'].split('_random')[-1]
        suffix = str(ys) + '_past' + s + '_'+str(int(candidates.loc[int(i), 'time']))
        fls = candidates.loc[int(i), 'random_model_t0']
        fls_list.append([suffix, fls])

    # run candidates until present
    pool = Pool()
    path_list = pool.map(partial(_run_to_present, gdir=gdir, ys=ys,
                                     ye=2000), fls_list)
    pool.close()
    pool.join()


def find_possible_glaciers(gdir, y0, n):
    # find good temp_bias_list
    random_df = find_temp_bias_range(gdir, y0)
    find_candidates(gdir, df=random_df, ys=y0, ye=2000,n=n)


def find_temp_bias_range(gdir, y0):
    """
    creates a pandas.DataFrame() with ALL created states. A subset of them will
    be tested later
    :param gdir:    oggm.GlacierDirectories
    :param y0:      int year of searched glaciers
    :return:        pandas.DataFrame()
    """
    t_eq = 0
    '''
    # try range (2,-2) first
    bias_list = [b.round(3) for b in np.arange(-2, 2, 0.05)]
    list = [(i**2, b) for i, b in enumerate(bias_list)]
    random_run_list = _run_random_parallel(gdir, y0, list)
    temp_b = -2
    while random_run_list['temp_bias'].min() == temp_b and temp_b >= -4:
        print(temp_b)
        n = len(random_run_list)
        list = [((i+n+1)**2, b.round(3)) for i, b in enumerate(np.arange(temp_b-1, temp_b-0.05, 0.05))]
        random_run_list = random_run_list.append(_run_random_parallel(gdir, y0, list), ignore_index=True)
        temp_b = temp_b-1
    '''
    # try range (2,-3) first  --> 100 runs
    bias_list = [b.round(3) for b in np.arange(-3, 2, 0.05)]
    list = [(i ** 2, b) for i, b in enumerate(bias_list)]
    random_run_list = _run_random_parallel(gdir, y0, list)

    # if temp bias = -3 does not create a glacier that exceeds boundary, we test further up to -5
    if random_run_list['temp_bias'].min() == -3:
        bias_list = [b.round(3) for b in np.arange(-5, -3, 0.05)]
        list = [(i ** 2, b) for i, b in enumerate(bias_list)]
        random_run_list = _run_random_parallel(gdir, y0, list)


    # check for zero glacier
    max_bias = random_run_list['temp_bias'].idxmax()
    p = gdir.get_filepath('model_run', filesuffix=random_run_list.loc[max_bias, 'suffix'])
    fmod = FileModel(p)

    if not fmod.volume_m3_ts().min() == 0:
        n = len(random_run_list)
        list = [((i + n + 1) ** 2, b.round(3)) for i, b in enumerate(np.arange(2.05, 3, 0.05))]
        random_run_list = random_run_list.append(_run_random_parallel(gdir, y0, list), ignore_index=True)
    i = 0
    # find t_eq
    for suffix in random_run_list['suffix'].values:
        if i < 10:
            try:
                rp = gdir.get_filepath('model_run', filesuffix=suffix)
                fmod = FileModel(rp)
                t = _find_extrema(fmod.volume_m3_ts())
                if t > t_eq:
                    t_eq = t
                i = i+1
            except:
                pass

    all = pd.DataFrame()
    for suffix in random_run_list['suffix']:
        try:
            rp = gdir.get_filepath('model_run', filesuffix=suffix)
            fmod = FileModel(rp)
            v = pd.DataFrame(fmod.volume_m3_ts()).reset_index()
            v = v[v['time'] >= t_eq]
            v = v.assign(suffix=lambda x: suffix)
            all = all.append(v, ignore_index=True)
        except:
            pass
    return all


def get_single_results(gdir,yr):
    """
    Creates a pd.DataFrame() containing all tested glaciers candidates in year
    yr. Read all "model_run+str(yr)+_past*.nc" files in gdir.dir
    :param gdir:    oggm.GlacierDirectory
    :param yr:      int, year of seachred glacier
    :return:
    """

    df = pd.DataFrame()
    prefix = 'model_run'+str(yr)+'_past'
    list = [f.split('model_run')[-1].split('.nc')[0] for f in os.listdir(gdir.dir) if
            f.startswith(prefix)]
    for f in list:
        try:
            rp = gdir.get_filepath('model_run',filesuffix=f)
            fmod = FileModel(rp)
            fmod_t = copy.deepcopy(fmod)
            fmod_t.run_until(2000)

            # read experiment
            ep = gdir.get_filepath('model_run',filesuffix='_experiment')
            emod = FileModel(ep)
            emod_t = copy.deepcopy(emod)
            emod_t.run_until(2000)

            obj = objective_value(fmod_t,emod_t)
            #fmod.reset_y0(yr)
            df = df.append({'model':copy.deepcopy(fmod),'objective':obj,'temp_bias':f.split('_')[-2],'time':f.split('_')[-1]},ignore_index=True)
        except:
            pass
    return df


def objective_value(model1, model2):
    """
    calculates the objective value (difference in geometry)
    :param model1: oggm.flowline.FluxBasedModel
    :param model2: oggm.flowline.FluxBasedModel
    :return:       float
    """

    model1 = copy.deepcopy(model1)
    model2 = copy.deepcopy(model2)
    model2.run_until(2000)
    model1.run_until(2000)

    fls1 = model1.fls
    fls2 = model2.fls
    objective=0
    for i in range(len(model1.fls)):
        objective = objective + np.sum(abs(fls1[i].surface_h-fls2[i].surface_h)**2)+ \
          np.sum(abs(fls1[i].widths-fls2[i].widths)**2)
    return objective


def prepare_for_initializing(gdirs):
    """
    oggm workflow for preparing initializing
    :param gdirs: list of oggm.GlacierDirectories
    :return None, but creates required files
    """
    list_tasks = [
        tasks.glacier_masks,
        tasks.compute_centerlines,
        tasks.initialize_flowlines,
        tasks.compute_downstream_line,
        tasks.compute_downstream_bedshape,
        tasks.catchment_area,
        tasks.catchment_intersections,
        tasks.catchment_width_geom,
        tasks.catchment_width_correction,
        tasks.process_histalp_data
    ]
    for task in list_tasks:
        workflow.execute_entity_task(task, gdirs)

    workflow.climate_tasks(gdirs)
    workflow.execute_entity_task(tasks.prepare_for_inversion, gdirs)

    #for gdir in gdirs:
    #    mass_conservation_inversion(gdir)

    workflow.execute_entity_task(mass_conservation_inversion, gdirs)
    #workflow.execute_entity_task(tasks.volume_inversion, gdirs)
    workflow.execute_entity_task(tasks.filter_inversion_output, gdirs)
    workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)


def synthetic_experiments_parallel(gdirs):
    """
    creates searched and observed glacier to test the method, need only to
    be run once

    :param gdirs: list of oggm.GlacierDirectories
    :return:
    """
    reset = True
    if os.path.isfile(gdirs[0].get_filepath('synthetic_experiment')):
        reset = utils.query_yes_no(
            'Running the function synthetic_experiments'
            ' will reset the previous results. Are you '
            ' sure you like to continue?')
    if not reset:
        return

    pool = mp.Pool()
    pool.map(_run_parallel_experiment, gdirs)
    pool.close()
    pool.join()
