import os
import copy
import pickle
from functools import partial

# External libs
import numpy as np
import pandas as pd
import multiprocessing as mp
from multiprocessing import Pool
from scipy.signal import argrelextrema

# locals
from oggm import workflow, tasks, utils, cfg
from oggm.core.inversion import mass_conservation_inversion
from oggm.core.flowline import FileModel, FluxBasedModel


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


def _single_calibration_run(gdir, mb_offset, ys,ye):
    """
    Creates the synthetic experiment for one glacier. model_geometry_experiment.nc
    will be saved in working directory.
    """

    # check, if this model_geometry already exists
    try:
        rp = gdir.get_filepath('model_geometry', filesuffix='_calibration_past_'+ str(mb_offset))
        model = FileModel(rp)

    # otherwise create calibration_run with mb_offset
    except:
        try:
            fls = gdir.read_pickle('model_flowlines')
            # run a 600 years random run with mb_offset
            tasks.run_random_climate(gdir, nyears=600, y0=ys, bias=mb_offset, seed=1, init_model_fls=fls,
                                     output_filesuffix='_calibration_random_'+str(mb_offset) )

            # construct s_OGGM --> previous glacier will be run forward from
            # ys - ye with past climate file

            tasks.run_from_climate_data(gdir, ys=ys, ye=ye, init_model_filesuffix='_calibration_random_'+str(mb_offset),
                                        bias=mb_offset, init_model_yr=600,
                                        output_filesuffix='_calibration_past_'+str(mb_offset))
            # return FileModel
            rp = gdir.get_filepath('model_geometry',filesuffix='_calibration_past_' + str(mb_offset))
            model = FileModel(rp)


        except:
            with open(os.path.join(gdir.dir, 'log.txt')) as log:

                error = list(log)[-1].split(';')[-1]

            return error

    return model



def _run_parallel_experiment(gdir, t0, te):
    """
    Creates the synthetic experiment for one glacier. model_geometry_synthetic_experiment.nc
    will be saved in working directory.
    """
    try:
        fls = gdir.read_pickle('model_flowlines')
        # try to run random climate with temperature bias -1
        tasks.run_random_climate(gdir, nyears=600, y0=t0, bias=0, seed=1, temperature_bias=-1, init_model_fls=fls,
                                 output_filesuffix='_random_synthetic_experiment')

        # construct observed glacier, previous glacier will be run forward from
        # t0 - te with past climate file
        tasks.run_from_climate_data(gdir, ys=t0, ye=te, init_model_filesuffix='_random_synthetic_experiment',
                                    init_model_yr=600, output_filesuffix='_synthetic_experiment')

        # remove output file from random run
        os.remove(gdir.get_filepath('model_geometry',filesuffix='_random_synthetic_experiment'))
        os.remove(gdir.get_filepath('model_diagnostics', filesuffix='_random_synthetic_experiment'))
    except:
        print('experiment failed : ' + str(gdir.rgi_id))


def _run_to_present(array, gdir, ys, ye, mb_offset):
    """
    Run glacier candidates forwards.
    """
    init_yr=array[0]
    init_filesuffix=array[1]

    s = init_filesuffix.split('_random')[-1]
    output_filesuffix = str(ys) + '_past' + s + '_'+str(int(init_yr))

    path = os.path.join(gdir.dir, str(ys), 'model_geometry' + output_filesuffix + '.nc')
    # does file already exists?
    if not os.path.exists(path):
        try:
            tasks.run_from_climate_data(gdir, ys=ys, ye=ye, bias=mb_offset, init_model_filesuffix=init_filesuffix,
                                        init_model_yr=init_yr, output_filesuffix=output_filesuffix)
            return output_filesuffix
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

def calibration_runs(gdirs, ys ):

    pool = Pool()
    pool.map(partial(find_mb_offset,ys=ys),gdirs)
    pool.close()
    pool.join()

def find_mb_offset(gdir, ys, a=-2000, b=2000):

    try:

        ye = gdir.rgi_date
        max_it = 15
        i = 0
        bounds = [a, b]

        df = pd.DataFrame()

        while i < max_it:
            mb_offset = round((bounds[0] + bounds[1]) / 2, 2)

            ex_mod2 = _single_calibration_run(gdir, mb_offset, ys, ye)
            if isinstance(ex_mod2, FileModel):
                diff = gdir.rgi_area_km2 - ex_mod2.area_km2_ts()[ye]
                error = ''
            # mb_offset needs to be set higher
            else:
                diff = -np.inf
                error = ex_mod2.split(':')[-1]


            df = df.append(pd.Series({'mb_offset': mb_offset, 'area_diff': diff, 'error':error}),
                           ignore_index=True)

            if (abs(diff) < 1e-4) or bounds[1] - bounds[0] <= 0.25:
                break

            elif diff<0:
                bounds[0] = mb_offset
            else:
                bounds[1] = mb_offset
            i += 1

        # best mb_offset found
        mb_offset = df.iloc[df.area_diff.abs().idxmin()].mb_offset

        df.to_csv(os.path.join(gdir.dir,'experiment_iteration.csv'))

        for file in os.listdir(gdir.dir):
            if file.startswith('model_geometry_calibration') and file.endswith('.nc') and not file.endswith('_'+str(mb_offset)+'.nc'):
                os.remove(os.path.join(gdir.dir,file))
            if file.startswith('model_diagnostics_calibration') and file.endswith('.nc') and not file.endswith('_'+str(mb_offset)+'.nc'):
                os.remove(os.path.join(gdir.dir,file))

    except:
        pass

def _run_random_parallel(gdir, y0, list, mb_offset):

    """
    Parallelize the run_random_task.
    """
    pool = Pool()
    paths = pool.map(partial(_run_random_task, gdir=gdir, y0=y0, mb_offset=mb_offset), list)
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


def _run_random_task(tupel, gdir, y0, mb_offset):
    """
    Run random model to create lots of possible states
    """
    seed = tupel[0]
    temp_bias = tupel[1]
    fls = gdir.read_pickle('model_flowlines')
    suffix = str(y0) + '_random_'+str(seed) + '_' + str(temp_bias)
    #path = gdir.get_filepath('model_rgeometry', filesuffix=suffix)
    path = os.path.join(gdir.dir, str(y0),'model_geometry'+suffix+'.nc')

    # does file already exists?
    if not os.path.exists(path):

        try:
            tasks.run_random_climate(gdir, nyears=600, y0=y0, bias=mb_offset,
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
    rp = gdir.get_filepath('model_geometry', filesuffix=suffix)
    fmod = FileModel(rp)
    fmod.run_until(ye)
    return copy.deepcopy(fmod)


def identification(gdir, list, ys, ye, n):
    """
    Determine glacier candidates and run them to the date of observation
    :param gdir:    oggm.GlacierDirectories
    :param df:      pd.DataFrame (volume_m3_ts() from random climate runs)
    :param ys:      starting year
    :param ye:      year of observation
    :param n:       number of candidates
    :return:
    """
    i = 0
    t_stag = 0
    # find t_stag
    for suffix in list['suffix'].values:
        if i < 10:
            try:
                rp = gdir.get_filepath('model_geometry', filesuffix=suffix)
                if not os.path.exists(rp):
                    rp = os.path.join(gdir.dir,str(ys),
                                      'model_geometry'+suffix+'.nc')
                fmod = FileModel(rp)
                t = _find_extrema(fmod.volume_m3_ts())
                if t > t_stag:
                    t_stag = t
                i = i+1
            except:
                pass

    # make sure that t_stag is not close to the end
    if t_stag>550:
        t_stag = 550

    df = pd.DataFrame()
    for suffix in list['suffix']:
        try:
            rp = gdir.get_filepath('model_geometry', filesuffix=suffix)
            if not os.path.exists(rp):
                rp = os.path.join(gdir.dir, str(ys),
                                  'model_geometry' + suffix + '.nc')
            fmod = FileModel(rp)
            v = pd.DataFrame(fmod.volume_m3_ts()).rename_axis('time').reset_index()
            v = v[v['time'] >= t_stag]
            v = v.assign(suffix=lambda x: suffix)
            df = df.append(v, ignore_index=True)
        except:
            pass

    indices = []
    # find nearest glacier state for each of the n volume classes (equidistant)
    for val in np.linspace(df.volume_m3.min(), df.volume_m3.max(), n):
        index = df.iloc[(df['volume_m3'] - val).abs().argsort()][:1].index[0]
        if not index in indices:
            indices = np.append(indices, index)
    candidates = df.loc[indices]
    candidates = candidates.sort_values(['suffix', 'time'])
    candidates = candidates.drop_duplicates()

    return candidates[['time', 'suffix']]

def find_possible_glaciers(gdir, y0, ye, n, ex_mod=None, mb_offset=0, delete=False):

    path = os.path.join(gdir.dir, 'result' + str(y0) + '.pkl')

    # if results are already there and number of candidates are the same, don't run it again
    if os.path.isfile(path):
        results = pd.read_pickle(path, compression='gzip')
        return results

    # 1. Generation of possible glacier states
    #    - Run random climate over 600 years with different temperature biases
    random_list = generation(gdir, y0, mb_offset)

    # 2. Identification of glacier candidates
    #    - Determine t_stag (begin of stagnation period)
    #    - Classification by volume (n equidistantly distributed classes)
    #    - Select one candidate by class
    #
    candidate_df = identification(gdir, random_list, ys=y0, ye=ye, n=n)

    # 3. Evaluation
    #    - Run each candidate forward from y0 to ye
    #    - Evaluate candidates based on the fitness function
    #    - Save all models in pd.DataFrame and write pickle
    #    - copy all model_geometry files to tarfile
    results = evaluation(gdir, candidate_df, y0, ye, ex_mod, mb_offset, delete)


    # find acceptable, 5th percentile and median
    if delete:
        # saves important outputs for evaluation based on experiment and based on fls
        save = {}

        # minimum
        save.update({'minimum_exp': results.loc[results.fitness.idxmin(), 'model'].split('/')[-1]})
        save.update({'minimum_fls': results.loc[results.fitness_fls.idxmin(), 'model'].split('/')[-1]})

        # acceptable
        results_exp = results[results.fitness <= 1]
        if len(results_exp) > 0:
            # acceptable
            save.update({'acc_min_exp': results_exp.loc[results_exp.length.idxmin(), 'model'].split('/')[-1]})
            save.update({'acc_max_exp': results_exp.loc[results_exp.length.idxmax(), 'model'].split('/')[-1]})

            # 5th percentile
            results_exp = results_exp[results_exp.fitness <= results_exp.fitness.quantile(0.05)]
            save.update({'perc_min_exp': results_exp.loc[results_exp.length.idxmin(), 'model'].split('/')[-1]})
            save.update({'perc_max_exp': results_exp.loc[results_exp.length.idxmax(), 'model'].split('/')[-1]})

            # median
            results_exp = results_exp.sort_values(by='length')
            l1 = len(results_exp)
            if l1 % 2:
                index_exp = int((l1 - 1) / 2)
            else:
                index_exp = int(l1 / 2)
            save.update({'median_exp': results_exp.iloc[index_exp].model.split('/')[-1]})

        results_fls = results[results.fitness_fls <= 1]
        if len(results_fls) > 0:
            save.update({'acc_min_fls': results_fls.loc[results_fls.length.idxmin(), 'model'].split('/')[-1]})
            save.update({'acc_max_fls': results_fls.loc[results_fls.length.idxmax(), 'model'].split('/')[-1]})

            results_fls = results_fls[results_fls.fitness_fls <= results_fls.fitness_fls.quantile(0.05)]
            save.update({'perc_min_fls': results_fls.loc[results_fls.length.idxmin(), 'model'].split('/')[-1]})
            save.update({'perc_max_fls': results_fls.loc[results_fls.length.idxmax(), 'model'].split('/')[-1]})

            results_fls = results_fls.sort_values(by='length')
            l2 = len(results_fls)
            if l2 % 2:
                index_fls = int((l2 - 1) / 2)
            else:
                index_fls = int(l2 / 2)

            save.update({'median_fls':results_exp.iloc[index_fls].model.split('/')[-1]})

        # save for later
        pickle.dump(save, open(os.path.join(gdir.dir,'initialization_output.pkl'), "wb"))


        # delete other files
        for file in os.listdir(gdir.dir):
            if file.startswith('model_geometry' + (str(y0))):
                if not file in set(save.values()):
                    os.remove(os.path.join(gdir.dir, file))

                    # remove diagnostic file, too
                    file = file.split('model_geometry')
                    file.insert(0, 'model_diagnostics')
                    file = os.path.join(gdir.dir,''.join(file))
                    try:
                        os.remove(file)
                    except:
                        pass

    else:

        # move all model_geometry* files from year y0 to a new directory --> avoids
        # that thousand of thousands files are created in gdir.dir

        utils.mkdir(os.path.join(gdir.dir, str(y0)), reset=False)
        for file in os.listdir(gdir.dir):
            if file.startswith('model_geometry' + (str(y0))):
                os.rename(os.path.join(gdir.dir, file),
                          os.path.join(gdir.dir, str(y0), file))
            elif file.startswith('model_diagnostics' + (str(y0))):
                os.remove(os.path.join(gdir.dir, file))

    return results


def generation(gdir, y0, mb_offset):

    """
    creates a pandas.DataFrame() with ALL created states. A subset of them will
    be tested later
    :param gdir:    oggm.GlacierDirectories
    :param y0:      int year of searched glaciers
    :return:        array
    """
    t_stag = 0
    # try range (2,-3) first  --> 100 runs
    bias_list = [b.round(3) for b in np.arange(-3, 2, 0.05)]
    list = [(i ** 2, b) for i, b in enumerate(bias_list)]
    random_run_list = _run_random_parallel(gdir, y0, list, mb_offset)

    # if temp bias = -3 does not create a glacier that exceeds boundary, we test further up to -5
    if random_run_list['temp_bias'].min() == -3:
        n = len(random_run_list)
        bias_list = [b.round(3) for b in np.arange(-5, -3, 0.05)]
        list = [((i+n+1) ** 2, b) for i, b in enumerate(bias_list)]
        random_run_list = random_run_list.append(_run_random_parallel(gdir, y0, list, mb_offset), ignore_index=True)

    # check for zero glacier
    max_bias = random_run_list['temp_bias'].idxmax()
    max_suffix = random_run_list.loc[max_bias, 'suffix']

    p = gdir.get_filepath('model_geometry',filesuffix=max_suffix)
    if not os.path.exists(p):
        p = os.path.join(gdir.dir, str(y0), 'model_geometry' + max_suffix + '.nc')
    fmod = FileModel(p)

    if not fmod.volume_m3_ts().min() == 0:
        n = len(random_run_list)
        list = [((i + n + 1) ** 2, b.round(3)) for i, b in enumerate(np.arange(2.05, 3, 0.05))]
        random_run_list = random_run_list.append(_run_random_parallel(gdir, y0, list, mb_offset), ignore_index=True)
    random_run_list = random_run_list.sort_values(by='temp_bias')
    return random_run_list


def evaluation(gdir, cand_df, y0, ye, emod, mb_offset, delete):

    """
    Creates a pd.DataFrame() containing all tested glaciers candidates in year
    yr. Read all "model_geometry+str(yr)+_past*.nc" files in gdir.dir
    :param gdir:        oggm.GlacierDirectory
    :param cand_df:     dataframe with init_model_filesuffix and init_time
    :param y0:          int, year of searched glacier
    :param ye:          int, year of observation
    :return:
    """

    # run candidates until present
    pool = Pool()
    suffix_list = pool.map(partial(_run_to_present, gdir=gdir, ys=y0,
                                 ye=ye, mb_offset=mb_offset), cand_df.to_numpy())
    pool.close()
    pool.join()

    df = pd.DataFrame()
    prefix = 'model_geomtry'+str(y0)+'_past'

    if emod is None:
        # read experiment
        ep = gdir.get_filepath('model_geometry', filesuffix='_synthetic_experiment')
        emod = FileModel(ep)
    emod_t = copy.deepcopy(emod)
    emod_t.run_until(ye)

    # get fls model
    fls = gdir.read_pickle('model_flowlines')
    fls_mod = FluxBasedModel(flowlines=fls)


    for f in suffix_list:

        try:
            # read past climate model runs and calculate objective
            rp = gdir.get_filepath('model_geometry', filesuffix=f)
            if not os.path.exists(rp):
                rp = os.path.join(gdir.dir, str(y0), 'model_geometry' + f + '.nc')

            fmod = FileModel(rp)
            fmod_t = copy.deepcopy(fmod)
            fmod_t.run_until(ye)
            fitness = fitness_value(fmod_t, emod_t, ye)
            fitness_fls = fitness_value_fls(fmod_t,fls_mod, ye)
            if not delete:
                df = df.append({'model': copy.deepcopy(fmod), 'fitness': fitness,
                                'fitness_fls':fitness_fls,'temp_bias': float(f.split('_')[-2]),
                                'time': f.split('_')[-1], 'volume': fmod.volume_km3},
                               ignore_index=True)
            else:
                df = df.append({'model':rp, 'fitness':fitness, 'temp_bias': float(f.split('_')[-2]),
                                'time': f.split('_')[-1], 'volume': fmod.volume_km3,'length': fmod.length_m,
                                'area': fmod.area_km2,'fitness_fls':fitness_fls,},
                               ignore_index=True)

        except:

            df = df.append({'model': None, 'fitness': None,
                            'temp_bias': float(f.split('_')[-2]),
                            'time': f.split('_')[-1], 'volume': None},ignore_index=True)


    if not delete:
        # save df with result models
        path = os.path.join(gdir.dir, 'result' + str(y0) + '.pkl')
        df.to_pickle(path, compression='gzip')

    return df

def fitness_value_fls(model1, model2, ye):
    """
    calculates the fitness value (difference in geometry)
    :param model1: oggm.flowline.FluxBasedModel
    :param model2: oggm.flowline.FluxBasedModel from fls (only year0)
    :param ye:     int, year of observation
    :return:       float
    """

    model1 = copy.deepcopy(model1)
    model2 = copy.deepcopy(model2)
    model2.run_until(0)
    model1.run_until(ye)

    fls1 = model1.fls
    fls2 = model2.fls
    fitness = 0
    m = 0
    for i in range(len(model1.fls)):
        fitness = fitness + np.sum(
            abs(fls1[i].surface_h - fls2[i].surface_h) ** 2) + \
                    np.sum(abs(fls1[i].widths - fls2[i].widths) ** 2)
        m = m + fls1[i].nx

    fitness = fitness / m
    fitness = fitness/125

    return fitness


def fitness_value(model1, model2, ye):
    """
    calculates the fitness value (difference in geometry)
    :param model1: oggm.flowline.FluxBasedModel
    :param model2: oggm.flowline.FluxBasedModel
    :param ye:     int, year of observation
    :return:       float
    """

    model1 = copy.deepcopy(model1)
    model2 = copy.deepcopy(model2)
    model2.run_until(ye)
    model1.run_until(ye)

    fls1 = model1.fls
    fls2 = model2.fls
    fitness = 0
    m = 0
    for i in range(len(model1.fls)):
        fitness = fitness + np.sum(
            abs(fls1[i].surface_h - fls2[i].surface_h) ** 2) + \
                    np.sum(abs(fls1[i].widths - fls2[i].widths) ** 2)
        m = m + fls1[i].nx

    fitness = fitness / m
    fitness = fitness/125

    return fitness


def preprocessing(gdirs):
    """
    oggm workflow for preparing initializing
    :param gdirs: list of oggm.GlacierDirectories from preprocessed level 2 onwards
    :return None, but creates required files
    """
    workflow.climate_tasks(gdirs)
    workflow.inversion_tasks(gdirs)
    workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)


def synthetic_experiments_parallel(gdirs, t0, te):
    """
    creates the synthetic experiments for all glaciers in gdirs in parallel, need only to
    be run once

    :param gdirs: list of oggm.GlacierDirectories
    :return:
    """
    reset = True
    if os.path.isfile(gdirs[0].get_filepath('model_geometry', filesuffix='_synthetic_experiment')):
        reset = utils.query_yes_no(
            'Running the function synthetic_experiments'
            ' will reset the previous results. Are you '
            ' sure you like to continue?')
    if not reset:
        return

    pool = mp.Pool()
    pool.map(partial(_run_parallel_experiment,t0=t0, te=te), gdirs)
    pool.close()
    pool.join()
