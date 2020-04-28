import os
import copy
import pickle
import shutil
import tarfile
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


def _run_experiment(gdir, temp_bias, bias, ys,ye):
    """
    Creates the synthetic experiment for one glacier. model_run_experiment.nc
    will be saved in working directory.
    """

    # check, if this experiment already exists
    try:
        rp = gdir.get_filepath('model_run', filesuffix='_advanced_experiment_' +str(temp_bias)+'_'+ str(bias))
        model = FileModel(rp)

    # otherwise create experiment
    except:
        try:
            fls = gdir.read_pickle('model_flowlines')
            model = tasks.run_random_climate(gdir, nyears=600, y0=ys, bias=bias, seed=1,
                                             temperature_bias=temp_bias,
                                             init_model_fls=fls,output_filesuffix='_random_experiment_'+str(temp_bias)+'_'+str(bias) )

            # construct observed glacier, previous glacier will be run forward from
            # 1917 - rgi_date with past climate file

            fls = copy.deepcopy(model.fls)
            tasks.run_from_climate_data(gdir, ys=ys, ye=ye, init_model_fls=fls,bias=bias,
                                        output_filesuffix='_advanced_experiment_'+str(temp_bias)+'_'+str(bias))
            # to return FileModel
            rp = gdir.get_filepath('model_run',filesuffix='_advanced_experiment_' + str(
                                       temp_bias) + '_' + str(bias))
            model = FileModel(rp)

        except:
            with open(os.path.join(gdir.dir,'log.txt')) as log:
                error=list(log)[-1].split(';')[-1]
            return error



    return model



def _run_parallel_experiment(gdir, t0, te):
    """
    Creates the synthetic experiment for one glacier. model_run_experiment.nc
    will be saved in working directory.
    """
    try:
        fls = gdir.read_pickle('model_flowlines')
        # try to run random climate with temperature bias -1

        model = tasks.run_random_climate(gdir, nyears=600, y0=t0, bias=0, seed=1,
                                         temperature_bias=-1,
                                         init_model_fls=fls)

        # construct observed glacier, previous glacier will be run forward from
        # 1850 - 2000 with past climate file
        b = fls[-1].bed_h

        fls = copy.deepcopy(model.fls)
        model = tasks.run_from_climate_data(gdir, ys=t0, ye=te, init_model_fls=fls,
                                    output_filesuffix='_experiment')
    except:
        print('experiment failed : ' + str(gdir.rgi_id))


def _run_to_present(tupel, gdir, ys, ye, bias):
    """
    Run glacier candidates forwards.
    """
    suffix = tupel[0]
    #path = gdir.get_filepath('model_run', filesuffix=suffix)
    path = os.path.join(gdir.dir, str(ys), 'model_run' + suffix + '.nc')
    # does file already exists?
    if not os.path.exists(path):
        try:
            tasks.run_from_climate_data(gdir, ys=ys, ye=ye, bias=bias,
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

def advanced_experiments(gdirs, temp_bias_list ,ys ):

    pool = Pool()
    pool.map(partial(find_residual,temp_bias_list=temp_bias_list,ys=ys),gdirs)
    pool.close()
    pool.join()

def find_residual(gdir, temp_bias_list, ys, a=-2000, b=2000):

    for temp_bias in temp_bias_list:
        try:

            fls = gdir.read_pickle('model_flowlines')
            mod = FluxBasedModel(flowlines=fls)

            ye = gdir.rgi_date
            max_it = 15
            i = 0
            bounds = [a, b]

            df = pd.DataFrame()

            while i < max_it:
                bias = round((bounds[0] + bounds[1]) / 2, 2)

                ex_mod2 = _run_experiment(gdir, temp_bias, bias, ys, ye)
                if isinstance(ex_mod2, FileModel):
                    diff = gdir.rgi_area_km2 - ex_mod2.area_km2_ts()[ye]
                    error = ''
                # bias needs to be set higher
                else:
                    diff = -np.inf
                    error = ex_mod2.split(':')[-1]


                df = df.append(pd.Series({'bias': bias, 'area_diff': diff, 'error':error}),
                               ignore_index=True)

                if (abs(diff) < 1e-4) or bounds[1] - bounds[0] <= 0.25:
                    break

                elif diff<0:
                    bounds[0] = bias
                else:
                    bounds[1] = bias
                i += 1

            # best bias found
            bias = df.iloc[df.area_diff.abs().idxmin()].bias
            df.to_csv(os.path.join(gdir.dir,'experiment_iteration.csv'))

            for file in os.listdir(gdir.dir):
                if file.startswith('model') and file.endswith('.nc') and not file.endswith('_'+str(bias)+'.nc'):
                    os.remove(os.path.join(gdir.dir,file))


        except:
            pass



def _run_random_parallel(gdir, y0, list, bias):

    """
    Parallelize the run_random_task.
    """
    pool = Pool()
    paths = pool.map(partial(_run_random_task, gdir=gdir, y0=y0, bias=bias), list)
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


def _run_random_task(tupel, gdir, y0, bias):
    """
    Run random model to create lots of possible states
    """
    seed = tupel[0]
    temp_bias = tupel[1]
    fls = gdir.read_pickle('model_flowlines')
    suffix = str(y0) + '_random_'+str(seed) + '_' + str(temp_bias)
    #path = gdir.get_filepath('model_run', filesuffix=suffix)
    path = os.path.join(gdir.dir, str(y0),'model_run'+suffix+'.nc')

    # does file already exists?
    if not os.path.exists(path):

        try:
            tasks.run_random_climate(gdir, nyears=600, y0=y0, bias=bias,
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
                rp = gdir.get_filepath('model_run', filesuffix=suffix)
                if not os.path.exists(rp):
                    rp = os.path.join(gdir.dir,str(ys),
                                      'model_run'+suffix+'.nc')
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
            rp = gdir.get_filepath('model_run', filesuffix=suffix)
            if not os.path.exists(rp):
                rp = os.path.join(gdir.dir, str(ys),
                                  'model_run' + suffix + '.nc')
            fmod = FileModel(rp)
            v = pd.DataFrame(fmod.volume_m3_ts()).reset_index()
            v = v[v['time'] >= t_stag]
            v = v.assign(suffix=lambda x: suffix)
            df = df.append(v, ignore_index=True)
        except:
            pass

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
        if not os.path.exists(rp):
            rp = os.path.join(gdir.dir, str(ys), 'model_run' + suffix + '.nc')
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
    return fls_list


def find_possible_glaciers(gdir, y0, ye, n, ex_mod=None, bias=0, delete=False):

    path = os.path.join(gdir.dir, 'result' + str(y0) + '.pkl')

    # if results are already there and number of candidates are the same, don't run it again
    if os.path.isfile(path):
        results = pd.read_pickle(path, compression='gzip')
        return results

    # 1. Generation of possible glacier states
    #    - Run random climate over 600 years with different temperature biases
    random_list = generation(gdir, y0, bias)

    # 2. Identification of glacier candidates
    #    - Determine t_stag(begin of stagnation period)
    #    - Classification by volume (n equidistantly distributed classes)
    #    - Select one candidate by class
    #
    candidate_list = identification(gdir, random_list, ys=y0, ye=ye, n=n)

    # 3. Evaluation
    #    - Run each candidate forward from y0 to ye
    #    - Evaluate candidates based on the fitnessfunction
    #    - Save all models in pd.DataFrame and write pickle
    #    - copy all model_run files to tarfile
    results = evaluation(gdir, candidate_list, y0, ye, ex_mod, bias, delete)


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
            if file.startswith('model_run' + (str(y0))):
                if not file in set(save.values()):
                    os.remove(os.path.join(gdir.dir, file))

                    # remove diagnostic file, too
                    file = file.split('model_run')
                    file.insert(0, 'model_diagnostics')
                    file = os.path.join(gdir.dir,''.join(file))
                    try:
                        os.remove(file)
                    except:
                        pass

    else:

        # move all model_run* files from year y0 to a new directory --> avoids
        # that thousand of thousands files are created in gdir.dir

        utils.mkdir(os.path.join(gdir.dir, str(y0)), reset=False)
        for file in os.listdir(gdir.dir):
            if file.startswith('model_run' + (str(y0))):
                os.rename(os.path.join(gdir.dir, file),
                          os.path.join(gdir.dir, str(y0), file))
            elif file.startswith('model_diagnostics' + (str(y0))):
                os.remove(os.path.join(gdir.dir, file))

    return results


def generation(gdir, y0, bias):

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
    random_run_list = _run_random_parallel(gdir, y0, list, bias)

    # if temp bias = -3 does not create a glacier that exceeds boundary, we test further up to -5
    if random_run_list['temp_bias'].min() == -3:
        n = len(random_run_list)
        bias_list = [b.round(3) for b in np.arange(-5, -3, 0.05)]
        list = [((i+n+1) ** 2, b) for i, b in enumerate(bias_list)]
        random_run_list = random_run_list.append(_run_random_parallel(gdir, y0, list, bias), ignore_index=True)

    # check for zero glacier
    max_bias = random_run_list['temp_bias'].idxmax()
    max_suffix = random_run_list.loc[max_bias, 'suffix']

    p = gdir.get_filepath('model_run',filesuffix=max_suffix)
    if not os.path.exists(p):
        p = os.path.join(gdir.dir, str(y0), 'model_run' + max_suffix + '.nc')
    fmod = FileModel(p)

    if not fmod.volume_m3_ts().min() == 0:
        n = len(random_run_list)
        list = [((i + n + 1) ** 2, b.round(3)) for i, b in enumerate(np.arange(2.05, 3, 0.05))]
        random_run_list = random_run_list.append(_run_random_parallel(gdir, y0, list, bias), ignore_index=True)
    random_run_list = random_run_list.sort_values(by='temp_bias')
    return random_run_list


def evaluation(gdir, fls_list, y0, ye, emod, bias, delete):

    """
    Creates a pd.DataFrame() containing all tested glaciers candidates in year
    yr. Read all "model_run+str(yr)+_past*.nc" files in gdir.dir
    :param gdir:        oggm.GlacierDirectory
    :param fls_list:    array with tupels [suffix, fls]
    :param y0:          int, year of searched glacier
    :param ye:          int, year of observation
    :return:
    """

    # run candidates until present
    pool = Pool()
    suffix_list = pool.map(partial(_run_to_present, gdir=gdir, ys=y0,
                                 ye=ye, bias=bias), fls_list)
    pool.close()
    pool.join()

    df = pd.DataFrame()
    prefix = 'model_run'+str(y0)+'_past'

    if emod is None:
        # read experiment
        ep = gdir.get_filepath('model_run', filesuffix='_experiment')
        emod = FileModel(ep)
    emod_t = copy.deepcopy(emod)
    emod_t.run_until(ye)

    # get fls model
    fls = gdir.read_pickle('model_flowlines')
    fls_mod = FluxBasedModel(flowlines=fls)


    for f in suffix_list:

        try:
            # read past climate model runs and calculate objective
            rp = gdir.get_filepath('model_run', filesuffix=f)
            if not os.path.exists(rp):
                rp = os.path.join(gdir.dir, str(y0), 'model_run' + f + '.nc')

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
    :param gdirs: list of oggm.GlacierDirectories
    :return None, but creates required files
    """
    workflow.gis_prepro_tasks(gdirs)
    workflow.climate_tasks(gdirs)
    workflow.inversion_tasks(gdirs)
    workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)


def synthetic_experiments_parallel(gdirs, t0, te):
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
    pool.map(partial(_run_parallel_experiment,t0=t0, te=te), gdirs)
    pool.close()
    pool.join()
