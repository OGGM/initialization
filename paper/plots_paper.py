import os
from functools import partial
from pylab import *
from oggm.core.flowline import FluxBasedModel, FileModel
import matplotlib.patches as mpatches
from oggm import graphics, tasks,utils
from matplotlib import cm
import xarray as xr
import pandas as pd
from multiprocessing import Pool
from copy import deepcopy
FlowlineModel = partial(FluxBasedModel, inplace=False)
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid.axes_grid import AxesGrid
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.legend_handler import HandlerLineCollection
import matplotlib.cm as cm
from matplotlib.path import Path
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
import matplotlib.cbook as cbook
#import ptitprince as pt
from matplotlib.collections import PolyCollection
import seaborn as sns
#from sklearn import preprocessing
import random
pd.options.mode.chained_assignment = None
warnings.filterwarnings("ignore")

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] = 30
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 25 #30
mpl.rcParams['lines.linewidth'] = 3


def add_at(ax, t, loc=2):
    fp = dict(size=30)
    _at = AnchoredText(t, loc=loc, prop=fp ,borderpad=0)
    _at.patch.set_linewidth(3)
    ax.add_artist(_at)
    return _at


def add_at2(ax, t, loc=2):
    fp = dict(size=10)
    _at = AnchoredText(t, loc=loc, prop=fp)
    _at.patch.set_linewidth(2)
    ax.add_artist(_at)
    return _at


class HandlerColorLineCollection(HandlerLineCollection):
    def create_artists(self, legend, artist, xdescent, ydescent,
                       width, height, fontsize, trans):
        x = np.linspace(0, width, self.get_numpoints(legend) + 1)
        y = np.zeros(
            self.get_numpoints(legend) + 1) + height / 2. - ydescent
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=artist.cmap,
                            transform=trans)
        lc.set_array(x)
        lc.set_linewidth(artist.get_linewidth())
        return [lc]


def plot_experiment(gdir, ex_mod, ys, plot_dir):

    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * \
        ex_mod.fls[-1].map_dx

    fig = plt.figure(figsize=(15, 14))
    grid = plt.GridSpec(2, 1, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[1, 0], sharex=ax1)

    if gdir.name != '':
        ax1.set_title(gdir.rgi_id+':'+gdir.name)
    else:
        ax1.set_title(gdir.rgi_id)

    # plot experiments, run until ys
    ex_mod = deepcopy(ex_mod)
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)
    i = np.where(ex_mod.fls[-1].thick > 0)[0][-1] + 10

    ax1.plot(x[:i], ex_mod.fls[-1].surface_h[:i], 'k:',
             label=r'$x_{'+str(ys)+'}^{exp}$', linewidth=3)
    ax1.plot(x[:i], ex_mod.fls[-1].bed_h[:i], 'k', label=r'$b$', linewidth=3)

    ex_mod.run_until(2000)

    ax2.plot(x[:i], ex_mod.fls[-1].surface_h[:i], 'k:',
             label=r'$x_{2000}^{exp = obs} $', linewidth=3)
    ax2.plot(x[:i], ex_mod.fls[-1].bed_h[:i], 'k', label=r'$b$', linewidth=3)

    # add figure names and legends
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)

    ax1.legend(loc=1)
    ax2.legend(loc=1)

    ax1.set_ylabel('Altitude (m)')
    ax1.set_xlabel('Distance along the main flowline (m)')
    ax2.set_ylabel('Altitude (m)')
    ax2.set_xlabel('Distance along the main flowline (m)')

    ax1.tick_params(axis='both', which='major')
    ax2.tick_params(axis='both', which='major')

    plot_dir = os.path.join(plot_dir, '00_experiment')
    utils.mkdir(plot_dir)
    fig_name = 'experiment_'+str(ys)+'_'+gdir.rgi_id
    plt.savefig(os.path.join(plot_dir, fig_name+'.pdf'), dpi=300)
    plt.savefig(os.path.join(plot_dir, fig_name+'.png'), dpi=300)
    plt.close()


def plot_candidates(gdir, df, yr, step, plot_dir):
    plot_dir = os.path.join(plot_dir, '06_candidates')
    utils.mkdir(plot_dir)
    fig, ax = plt.subplots(figsize=(10, 10))
    for file in os.listdir(os.path.join(gdir.dir, str(yr))):
        if file.startswith('model_run'+str(yr)+'_random'):
            suffix = file.split('model_run')[1].split('.nc')[0]
            rp = os.path.join(gdir.dir, str(yr), 'model_run'+suffix+'.nc')
            try:
                fmod = FileModel(rp)
                fmod.volume_m3_ts().plot(ax=ax, color='grey', label='',
                                         zorder=1)

            except:
                pass

    # last one again for labeling
    label = r'temperature bias $\in [$' + str(
        df['temp_bias'].min()) + ',' + str(df['temp_bias'].max()) + '$]$'
    df.time = df.time.apply(lambda x: int(x))
    t_eq = df['time'].sort_values().iloc[0]

    df['Fitness value'] = df.fitness
    plt.title(gdir.rgi_id)

    if step == 'step1':
        fmod.volume_m3_ts().plot(ax=ax, color='grey', label=label, zorder=1)
        plt.legend(loc=0, fontsize=23)
        plt.xlabel('Time (years)')
        plt.ylabel(r'Volume $(m^3)$')
        plt.savefig(os.path.join(plot_dir, 'candidates1_' + str(yr) + '_' +
                                 str(gdir.rgi_id) + '.png'), dpi=300)
    elif step == 'step2':
        ax.axvline(x=int(t_eq), color='k', zorder=1, label=r'$t_{stag}$')
        fmod.volume_m3_ts().plot(ax=ax, color='grey', label='', zorder=1)
        # black points
        df.plot.scatter(x='time', y='volume', ax=ax, color='k',
                        label='candidates', s=250, zorder=2)
        plt.legend(loc=0, fontsize=23)
        plt.xlabel('Time (years)')
        plt.ylabel(r'Volume $(m^3)$')
        plt.savefig(os.path.join(plot_dir, 'candidates2_' + str(yr) + '_' +
                                 str(gdir.rgi_id) + '.png'), dpi=300)
    elif step == 'step3':
        fmod.volume_m3_ts().plot(ax=ax, color='grey', label=None, zorder=1)
        ax.axvline(x=int(t_eq), color='k', zorder=1)

        cmap = matplotlib.cm.get_cmap('viridis')

        df.plot.scatter(x='time', y='volume', ax=ax, c='Fitness value',
                        colormap='viridis',
                        norm=mpl.colors.LogNorm(vmin=0.01, vmax=1000, clip=True),
                        s=250, edgecolors='k', zorder=2)
        # plot again points with objective == 0, without norm
        if len(df[df.fitness == 0]) > 0:
            df[df.fitness == 0].plot.scatter(x='time', y='volume', ax=ax,
                                               c=cmap(0), s=250,
                                               edgecolors='k', zorder=2)

        plt.xlabel('Time (years)')
        plt.ylabel(r'Volume $(m^3)$')
        plt.savefig(os.path.join(plot_dir, 'candidates3_' + str(yr) + '_' +
                                 str(gdir.rgi_id) + '.png'), dpi=300)
    plt.close()

    plt.figure(figsize=(15, 14))
    plt.hist(df.volume.values, bins=20)
    plt.xlabel(r'Volume $(m^3)$')
    plt.ylabel(r'Frequency')
    plt.title(gdir.rgi_id)
    plt.savefig(os.path.join(plot_dir, 'hist_candidates' + str(yr) + '_' +
                             str(gdir.rgi_id) + '.png'), dpi=300)
    plt.close()


def plot_compare_fitness(gdir, df, ex_mod, ys, plot_dir):

    plot_dir = os.path.join(plot_dir, '05_compare_fitness')
    utils.mkdir(plot_dir)

    fig = plt.figure(figsize=(15, 14))

    grid = plt.GridSpec(3, 1, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[1, 0], sharex=ax1)
    ax3 = plt.subplot(grid[2, 0], sharex=ax1)

    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)

    if gdir.name != '':
        ax1.set_title(gdir.rgi_id+': '+gdir.name, fontsize=30)
    elif gdir.rgi_id.endswith('779'):
        ax1.set_title(gdir.rgi_id + ': Guslarferner', fontsize=30)
    else:
        ax1.set_title(gdir.rgi_id, fontsize=30)

    df['fitness2'] = df.model.apply(
        lambda x: abs(x.area_m2_ts()[2000] - ex_mod.area_m2_ts()[2000])**2)
    df['fitness3'] = df.model.apply(
        lambda x: abs(x.length_m_ts()[2000] - ex_mod.length_m_ts()[2000])**2)

    norm = mpl.colors.LogNorm(vmin=df.fitness.min() + 0.01,
                              vmax=df.fitness.max(), clip=True)
    norm2 = mpl.colors.LogNorm(vmin=df.fitness2.min() + 0.01,
                               vmax=df.fitness2.max(), clip=True)
    norm3 = mpl.colors.LogNorm(vmin=df.fitness3.min() + 0.01,
                               vmax=df.fitness3.max(), clip=True)

    cmap = matplotlib.cm.get_cmap('viridis')

    # plot default objective
    df = df.sort_values('fitness', ascending=False)
    for i, model in df['model'].iteritems():
        if df.fitness.min() != df.fitness.max():
            color = cmap(norm(df.loc[i, 'fitness']))
        else:
            color = cmap(df.loc[i, 'fitness'])
        model.volume_m3_ts().plot(ax=ax1, color=[color],
                                  linewidth=3,
                                  label=r'$s_{1850-2000}^{exp}$')

    # plot objective based on area
    df = df.sort_values('fitness2', ascending=False)
    for i, model in df['model'].iteritems():
        if df.fitness2.min() != df.fitness2.max():
            color = cmap(norm2(int(df.loc[i, 'fitness2'])))
        else:
            color = cmap(df.loc[i, 'fitness2'])

        model.volume_m3_ts().plot(ax=ax2, color=[color],
                                  linewidth=3, label=r'$s_{1850-2000}^{exp}$')

    # plot objective based on length variations
    df = df.sort_values('fitness3', ascending=False)
    for i, model in df['model'].iteritems():
        if df.fitness3.min() != df.fitness3.max():
            color = cmap(norm3(int(df.loc[i, 'fitness3'])))
        else:
            color = cmap(df.loc[i, 'fitness3'])
        model.volume_m3_ts().plot(ax=ax3, color=[color],
                                  linewidth=3, label=r'$s_{1850-2000}^{exp}$')
    # plot experiments
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_m3_ts().plot(ax=ax1, color='red', linestyle=':', linewidth=3,
                               label=r'$s_{1850-2000}^{exp}$')
    ex_mod.volume_m3_ts().plot(ax=ax2, color='red', linestyle=':', linewidth=3,
                               label=r'$s_{1850-2000}^{exp}$')
    ex_mod.volume_m3_ts().plot(ax=ax3, color='red', linestyle=':', linewidth=3,
                               label=r'$s_{1850-2000}^{exp}$')
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax1, ax2, ax3])
    cbar = fig.colorbar(sm, cax=cax, **kw)
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label('Normalized fitness value', fontsize=30)

    # add figure names and legends
    add_at(ax1, r"a", loc=1)
    add_at(ax2, r"b", loc=1)
    add_at(ax3, r"c", loc=1)

    ax2.set_ylabel(r'Volume ($km^3$)', fontsize=30)
    ax3.set_xlabel(r'Time', fontsize=30)

    ax1.tick_params(axis='both', which='major', labelsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=30)

    plt.savefig(os.path.join(plot_dir, gdir.rgi_id+'.pdf'), dpi=300)
    #plt.show()
    plt.close()


def plot_fitness_over_time(gdir, df_list, ex_mod, plot_dir):

    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(111, yscale='linear')
    norm = mpl.colors.LogNorm(vmin=0.1, vmax=1e5)
    cmap = matplotlib.cm.get_cmap('viridis')

    volumes = np.linspace(df_list['1850'].volume.min(),
                          df_list['1850'].volume.max(), 100)
    # width of the patches
    yrs = [int(yr) for yr in list(df_list.keys())]
    yrs.sort()
    w = yrs[1] - yrs[0]

    for year, df in df_list.items():
        color = []
        for i in range(len(volumes)):
            if i != len(volumes) - 1:
                part = df[(df.volume >= volumes[i]) & (
                    df.volume <= volumes[i + 1])]
                color.append(part.objective.min())
            else:
                part = df[df.volume >= volumes[i]]
                color.append(part.objective.mean())  # or min

        # interpolate missing data
        missing = np.where(np.isnan(color))
        if len(missing[0]) != 0:
            xp = np.delete(range(len(volumes)), missing)
            fp = np.delete(color, missing)
            missing_y = np.interp(missing, xp, fp)
            for i, j in enumerate(missing[0]):
                color[j] = missing_y[0][i]

        year = np.zeros(len(volumes)) + int(year)
        for x, y, c in zip(volumes, year, color):
            if np.isnan(c):
                color = 'white'
            else:
                color = cmap(norm(c))
            ax.add_patch(
                Rectangle((y-(w/2), x), w, volumes[1] - volumes[0],
                          color=color))

    # add experiment in plot
    ex_mod.volume_m3_ts().plot(ax=ax, linestyle=':', color='red', linewidth=5)
    ax.set_xlim(yrs[0]-(w/2), yrs[-1]+(w/2))
    ax.set_ylim(volumes[0], volumes[-1])
    if gdir.rgi_id.endswith('779'):
        plt.title(gdir.rgi_id + ':  Guslarferner', fontsize=25)
        add_at(ax, r"b", loc=3)
    elif gdir.name == '':
        plt.title(gdir.rgi_id, fontsize=25)
        add_at(ax, r"c", loc=3)
    else:
        plt.title(gdir.rgi_id + ': ' + gdir.name, fontsize=25)
        add_at(ax, r"a", loc=3)

    plt.ylabel(r'Volume ($m^3$)', fontsize=25)
    plt.xlabel(r'Starting time', fontsize=25)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax])
    cbar = fig.colorbar(sm, cax=cax, **kw)
    cbar.ax.tick_params(labelsize=25)
    cbar.set_label('Fitness value', fontsize=25)
    plt.savefig(os.path.join(plot_dir, 'starting' '_' + gdir.rgi_id + '.png'))
    plt.close()


def plot_col_fitness(gdir, df, ex_mod, ys, plot_dir):

    plot_dir = os.path.join(plot_dir, '03_surface_by_fitness')
    utils.mkdir(plot_dir)
    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * \
        ex_mod.fls[-1].map_dx
    fig = plt.figure(figsize=(25, 18))
    grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1], sharey=ax1)
    ax3 = plt.subplot(grid[1, :])

    if gdir.name != '':
        plt.suptitle(gdir.rgi_id+': '+gdir.name, fontsize=30)
    elif gdir.rgi_id.endswith('779'):
        plt.suptitle(gdir.rgi_id+': Guslarferner', fontsize=30)
    else:
        plt.suptitle(gdir.rgi_id, fontsize=30)

    norm = mpl.colors.LogNorm(vmin=0.1, vmax=1e5)
    cmap = matplotlib.cm.get_cmap('viridis')

    #df = df.sort_values('objective', ascending=False)
    df = df.sort_values('fitness', ascending=False)

    for i, model in df['model'].iteritems():
        model = deepcopy(model)
        model.reset_y0(ys)
        #color = cmap(norm(df.loc[i, 'objective']))
        color = cmap(norm(df.loc[i, 'fitness']))
        ax1.plot(x, deepcopy(model.fls[-1].surface_h), color=color,
                 label='')
        model.volume_m3_ts().plot(ax=ax3, color=[color], label='')
        model.run_until(2000)

        ax2.plot(x, model.fls[-1].surface_h, color=color, label='')

    # plot experiments
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_m3_ts().plot(ax=ax3, color='red', linestyle=':', linewidth=3,
                               label='')
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)

    ax1.plot(x, ex_mod.fls[-1].surface_h, ':', color='red', label='',
             linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k', label='', linewidth=3)

    ex_mod.run_until(2000)

    ax2.plot(x, ex_mod.fls[-1].surface_h, ':', color='red', label='',
             linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', label='',
             linewidth=3)

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax1, ax2, ax3])
    cbar = fig.colorbar(sm, cax=cax, **kw)
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label('Fitness value', fontsize=30)

    # add figure names and x-/ylabels
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)
    add_at(ax3, r"c", loc=3)

    ax1.set_ylabel('Altitude (m)', fontsize=30)
    ax1.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax2.set_ylabel('Altitude (m)', fontsize=30)
    ax2.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax3.set_ylabel(r'Volume ($m^3$)', fontsize=30)
    ax3.set_xlabel('Time (years)', fontsize=30)

    ax1.tick_params(axis='both', which='major', labelsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=30)
    ax3.tick_params(axis='both', which='major', labelsize=30)
    ax3.yaxis.offsetText.set_fontsize(30)

    # add legend
    # legend

    t = np.linspace(0, 10, 200)
    x = np.cos(np.pi * t)
    y = np.sin(t)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap,
                        norm=plt.Normalize(0, 10), linewidth=3)
    lc2 = LineCollection(segments, color='k',
                         norm=plt.Normalize(0, 10), linewidth=3)
    lc3 = LineCollection(segments, color='r', linestyle=':',
                         norm=plt.Normalize(0, 10), linewidth=3)

    l1 = ax1.legend(handles=[lc3, lc, lc2], handler_map={
        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=[r'$x_{'+str(ys)+'}^{exp}$', r'$x_{'+str(ys)+'}$',
                            r'$b$'], loc=1)

    l2 = ax2.legend(handles=[lc3, lc, lc2],
                    handler_map={
                        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=[r'$x_{2000}^{obs}$', r'$x_{2000}$',
                            r'$b$'], loc=1)

    l3 = ax3.legend(handles=[lc3, lc],
                    handler_map={
                        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=[r'$s_{'+str(ys)+'-2000}^{exp}$', r'$s_{'+str(ys) +
                            '-2000}$'], loc=1)

    l1.set_zorder(0)
    l2.set_zorder(0)
    l3.set_zorder(0)

    ax3.set_xlim(xmin=1847, xmax=2003)
    fig_name = 'surface_'+str(ys)+'_'+gdir.rgi_id
    plt.savefig(os.path.join(plot_dir, fig_name+'.pdf'), dpi=300)
    plt.savefig(os.path.join(plot_dir, fig_name + '.png'), dpi=300)
    #plt.show()
    plt.close()

def plot_col_fitness_average(gdir, df, ex_mod, ys, plot_dir):

    plot_dir = os.path.join(plot_dir, '03_surface_by_fitness')
    utils.mkdir(plot_dir)
    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * \
        ex_mod.fls[-1].map_dx
    fig = plt.figure(figsize=(25, 18))
    grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1], sharey=ax1)
    ax3 = plt.subplot(grid[1, :])

    if gdir.name != '':
        plt.suptitle(gdir.rgi_id + ': ' + gdir.name, fontsize=30)
    elif gdir.rgi_id.endswith('779'):
        plt.suptitle(gdir.rgi_id + ': Guslarferner', fontsize=30)
    else:
        plt.suptitle(gdir.rgi_id, fontsize=30)

    norm = mpl.colors.LogNorm(vmin=0.01, vmax=1000)
    cmap = matplotlib.cm.get_cmap('viridis')

    # df = df.sort_values('objective', ascending=False)
    df = df.sort_values('fitness', ascending=False)

    for i, model in df['model'].iteritems():
        model = deepcopy(model)
        model.reset_y0(ys)
        # color = cmap(norm(df.loc[i, 'objective']))
        color = cmap(norm(df.loc[i, 'fitness']))
        ax1.plot(x, deepcopy(model.fls[-1].surface_h), color=color,
                 label='')
        model.volume_m3_ts().plot(ax=ax3, color=[color], label='')
        model.run_until(2000)

        ax2.plot(x, model.fls[-1].surface_h, color=color, label='')

    # plot experiments
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_m3_ts().plot(ax=ax3, color='red', linestyle=':',
                               linewidth=3,
                               label='')
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)

    ax1.plot(x, ex_mod.fls[-1].surface_h, ':', color='red', label='',
             linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k', label='', linewidth=3)

    ex_mod.run_until(2000)

    ax2.plot(x, ex_mod.fls[-1].surface_h, ':', color='red', label='',
             linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', label='',
             linewidth=3)

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax1, ax2, ax3])
    cbar = fig.colorbar(sm, cax=cax, **kw)
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label('Fitness value', fontsize=30)

    # add figure names and x-/ylabels
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)
    add_at(ax3, r"c", loc=3)

    ax1.set_ylabel('Altitude (m)', fontsize=30)
    ax1.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax2.set_ylabel('Altitude (m)', fontsize=30)
    ax2.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax3.set_ylabel(r'Volume ($m^3$)', fontsize=30)
    ax3.set_xlabel('Time (years)', fontsize=30)

    ax1.tick_params(axis='both', which='major', labelsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=30)
    ax3.tick_params(axis='both', which='major', labelsize=30)
    ax3.yaxis.offsetText.set_fontsize(30)

    # add legend
    # legend

    t = np.linspace(0, 10, 200)
    x = np.cos(np.pi * t)
    y = np.sin(t)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap,
                        norm=plt.Normalize(0, 10), linewidth=3)
    lc2 = LineCollection(segments, color='k',
                         norm=plt.Normalize(0, 10), linewidth=3)
    lc3 = LineCollection(segments, color='r', linestyle=':',
                         norm=plt.Normalize(0, 10), linewidth=3)

    l1 = ax1.legend(handles=[lc3, lc, lc2], handler_map={
        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=[r'$x_{' + str(ys) + '}^{exp}$',
                            r'$x_{' + str(ys) + '}$',
                            r'$b$'], loc=1)

    l2 = ax2.legend(handles=[lc3, lc, lc2],
                    handler_map={
                        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=[r'$x_{2000}^{obs}$', r'$x_{2000}$',
                            r'$b$'], loc=1)

    l3 = ax3.legend(handles=[lc3, lc],
                    handler_map={
                        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=[r'$s_{' + str(ys) + '-2000}^{exp}$',
                            r'$s_{' + str(ys) +
                            '-2000}$'], loc=1)

    l1.set_zorder(0)
    l2.set_zorder(0)
    l3.set_zorder(0)

    ax3.set_xlim(xmin=1847, xmax=2003)
    fig_name = 'surface_' + str(ys) + '_' + gdir.rgi_id
    plt.savefig(os.path.join(plot_dir, fig_name + '.pdf'), dpi=300)
    plt.savefig(os.path.join(plot_dir, fig_name + '.png'), dpi=300)
    #plt.show()
    plt.close()

def plot_median_average(gdir, df,eps, ex_mod, ys, ye, plot_dir):
    plot_dir = os.path.join(plot_dir, '04_median')
    utils.mkdir(plot_dir)
    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[
        -1].map_dx
    fig = plt.figure(figsize=(25, 18))
    grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1], sharey=ax1)
    ax3 = plt.subplot(grid[1, :])

    if gdir.name != '':
        plt.suptitle(gdir.rgi_id + ': ' + gdir.name, fontsize=30)
    elif gdir.rgi_id.endswith('779'):
        plt.suptitle(gdir.rgi_id + ': Guslarferner', fontsize=30)
    else:
        plt.suptitle(gdir.rgi_id, fontsize=30)

    df = df.sort_values('fitness', ascending=False)

    # acceptable glacier states
    df = df[df.fitness < eps]
    s_t0 = pd.DataFrame()
    s_te = pd.DataFrame()
    v = pd.DataFrame()
    for i in df.index:
        model = deepcopy(df.loc[i, 'model'])
        s_t0 = s_t0.append(pd.Series(model.fls[-1].surface_h),
                                 ignore_index=True)
        model.run_until(2000)
        s_te = s_te.append(pd.Series(model.fls[-1].surface_h),
                                 ignore_index=True)
        v = v.append(model.volume_m3_ts(),ignore_index=True)


    ax1.fill_between(x, s_t0.max().values, s_t0.min().values, alpha=0.3, color='grey',
                     label=r'$\mathcal{S}_{' + str(ys) + '}^{' + str(
                         eps) + '}$')
    ax2.fill_between(x, s_te.max().values, s_te.min().values, alpha=0.3, color='grey',
                     label=r'$\mathcal{S}_{' + str(ye) + '}^{' + str(
                         eps) + '}$')
    ax3.fill_between(model.volume_m3_ts().index,
                     v.max().values,
                     v.min().values, alpha=0.3, color='grey',
                     label=r'$\mathcal{S}_{' +str(ys)+'-'+str(ye) +'}^{' + str(eps) + '}$')

    # 5% quantile and median of 5% quantile

    df = df[df.fitness <= df.fitness.quantile(0.05)]
    s_t0 = pd.DataFrame()
    s_te = pd.DataFrame()
    v = pd.DataFrame()
    for i in df.index:
        model = deepcopy(df.loc[i, 'model'])
        s_t0 = s_t0.append(pd.Series(model.fls[-1].surface_h),
                           ignore_index=True)
        model.run_until(2000)
        s_te = s_te.append(pd.Series(model.fls[-1].surface_h),
                           ignore_index=True)
        v = v.append(model.volume_m3_ts(), ignore_index=True)

    ax1.fill_between(x, s_t0.max().values, s_t0.min().values, alpha=0.5,
                     label=r'$Q_{0.05}(\mathcal{S}_{'+str(ys)+'}^{'+str(eps)+'})$')

    ax2.fill_between(x, s_te.max().values, s_te.min().values, alpha=0.5,
                     label=r'$Q_{0.05}(\mathcal{S}_{2000}^{'+str(eps)+'})$')

    ax3.fill_between(v.columns, v.max().values, v.min().values, alpha=0.5, linewidth=3,
                     label=r'$Q_{0.05}(\mathcal{S}_{'+str(ys)+'-2000}^{'+str(eps)+'})$')

    # median of 5% quantile
    df.loc[:, 'length'] = df.model.apply(lambda x: x.length_m)
    df = df.sort_values('length', ascending=False)
    l = len(df)
    if l % 2:
        index = int((l - 1) / 2)
    else:
        index = int(l / 2)

    median_model = deepcopy(df.iloc[index].model)
    median_model.volume_m3_ts().plot(ax=ax3, linewidth=3, label='median')
    median_model.reset_y0(1850)
    median_model.run_until(ys)

    ax1.plot(x, median_model.fls[-1].surface_h, label='median',
             linewidth=3)
    median_model.run_until(2000)
    ax2.plot(x, median_model.fls[-1].surface_h, label='median',
             linewidth=3)

    # min model
    min_mod = deepcopy(df.loc[df.fitness.idxmin(), 'model'])
    min_mod.volume_m3_ts().plot(ax=ax3, color='r',
                                linewidth=3, label='minimum')
    min_mod.reset_y0(ys)
    min_mod.run_until(ys)

    ax1.plot(x, min_mod.fls[-1].surface_h, 'r', label='minimum',
             linewidth=3)

    min_mod.run_until(ye)

    ax2.plot(x, min_mod.fls[-1].surface_h, 'r', label='minimum',
             linewidth=3)

    # experiment
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_m3_ts().plot(ax=ax3, color='k', linestyle=':',
                               linewidth=3, label='')
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)

    ax1.plot(x, ex_mod.fls[-1].surface_h, 'k:', label='',
             linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k', label='', linewidth=3)

    ex_mod.run_until(2000)

    ax2.plot(x, ex_mod.fls[-1].surface_h, 'k:', label='',
             linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', label='', linewidth=3)

    # add figure names and x-/ylabels
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)
    add_at(ax3, r"c", loc=3)

    ax1.set_ylabel('Altitude (m)', fontsize=30)
    ax1.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax2.set_ylabel('Altitude (m)', fontsize=30)
    ax2.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax3.set_ylabel(r'Volume ($m^3$)', fontsize=30)
    ax3.set_xlabel('Time (years)', fontsize=30)

    ax1.tick_params(axis='both', which='major', labelsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=30)
    ax3.tick_params(axis='both', which='major', labelsize=30)
    ax3.yaxis.offsetText.set_fontsize(30)
    ax3.set_xlim(xmin=1847, xmax=2003)

    l1 = ax1.legend(loc=1, fontsize=30)
    l1.set_zorder(1)

    l2 = ax2.legend(loc=1, fontsize=30)
    l2.set_zorder(1)

    l3 = ax3.legend(loc=1, fontsize=30)
    l3.set_zorder(1)

    fig_name = 'median_'+str(ys)+'_'+gdir.rgi_id
    #plt.savefig(os.path.join(plot_dir, fig_name+'.pdf'), dpi=300)
    plt.savefig(os.path.join(plot_dir, fig_name+'.png'), dpi=300)

    #plt.show()
    plt.close()

    return median_model


def plot_median(gdir, df, ex_mod, ys, plot_dir):
    plot_dir = os.path.join(plot_dir, '04_median')
    utils.mkdir(plot_dir)
    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[
        -1].map_dx
    fig = plt.figure(figsize=(25, 18))
    grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1], sharey=ax1)
    ax3 = plt.subplot(grid[1, :])

    if gdir.name != '':
        plt.suptitle(gdir.rgi_id + ': ' + gdir.name, fontsize=30)
    elif gdir.rgi_id.endswith('779'):
        plt.suptitle(gdir.rgi_id + ': Guslarferner', fontsize=30)
    else:
        plt.suptitle(gdir.rgi_id, fontsize=30)

    # df = df.sort_values('objective', ascending=False)
    # df = df[df.objective < 100]
    df = df.sort_values('fitness', ascending=False)
    df = df[df.fitness < 100]

    min_id = df.volume.idxmin()
    max_id = df.volume.idxmax()

    min_model = deepcopy(df.loc[min_id, 'model'])
    max_model = deepcopy(df.loc[max_id, 'model'])

    min_model.reset_y0(ys)
    max_model.reset_y0(ys)
    ax1.fill_between(x, deepcopy(min_model.fls[-1].surface_h),
                     deepcopy(max_model.fls[-1].surface_h), alpha=0.3,
                     color='grey', label=r'$\mathcal{S}_{'+str(ys)+'}^{100}$')

    ax3.fill_between(min_model.volume_m3_ts().index, min_model.volume_m3_ts(),
                     max_model.volume_m3_ts(),alpha=0.3, color='grey',
                     label=r'$\mathcal{S}_{'+str(ys)+'-2000}^{100}$')
    min_model.run_until(2000)
    max_model.run_until(2000)
    ax2.fill_between(x, deepcopy(min_model.fls[-1].surface_h),
                     deepcopy(max_model.fls[-1].surface_h), alpha=0.3,
                     color='grey', label=r'$\mathcal{S}_{2000}^{100}$')
    # 5% quantile and median of 5% quantile

    #quant_df = df[df.objective <= df.objective.quantile(0.05)]
    quant_df = df[df.fitness <= df.fitness.quantile(0.05)]

    q_min_id = quant_df.volume.idxmin()
    q_max_id = quant_df.volume.idxmax()

    q_min_model = deepcopy(quant_df.loc[q_min_id, 'model'])
    q_max_model = deepcopy(quant_df.loc[q_max_id, 'model'])

    q_min_model.reset_y0(ys)
    q_max_model.reset_y0(ys)
    ax1.fill_between(x, deepcopy(q_min_model.fls[-1].surface_h),
                     deepcopy(q_max_model.fls[-1].surface_h), alpha=0.5,
                     label=r'$Q_{0.05}(\mathcal{S}_{'+str(ys)+'}^{100})$')

    ax3.fill_between(q_min_model.volume_m3_ts().index, q_min_model.volume_m3_ts(),
                     q_max_model.volume_m3_ts(), alpha=0.5, linewidth=3,
                     label=r'$Q_{0.05}(\mathcal{S}_{'+str(ys)+'-2000}^{100})$')
    q_min_model.run_until(2000)
    q_max_model.run_until(2000)
    ax2.fill_between(x, deepcopy(q_min_model.fls[-1].surface_h),
                     deepcopy(q_max_model.fls[-1].surface_h), alpha=0.5,
                     label=r'$Q_{0.05}(\mathcal{S}_{2000}^{100})$')

    # median of 5% quantile
    quant_df.loc[:, 'length'] = quant_df.model.apply(lambda x: x.length_m)
    quant_df = quant_df.sort_values('length', ascending=False)
    l = len(quant_df)
    if l % 2:
        index = int((l - 1) / 2)
    else:
        index = int(l / 2)

    median_model = deepcopy(quant_df.iloc[index].model)
    median_model.volume_m3_ts().plot(ax=ax3, linewidth=3, label='median')
    median_model.reset_y0(1850)
    median_model.run_until(ys)

    ax1.plot(x, median_model.fls[-1].surface_h, label='median',
             linewidth=3)
    median_model.run_until(2000)
    ax2.plot(x, median_model.fls[-1].surface_h, label='median',
             linewidth=3)




    # experiment
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_m3_ts().plot(ax=ax3, color='k', linestyle=':',
                               linewidth=3, label='')
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)

    ax1.plot(x, ex_mod.fls[-1].surface_h, 'k:', label='',
             linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k', label='', linewidth=3)

    ex_mod.run_until(2000)

    ax2.plot(x, ex_mod.fls[-1].surface_h, 'k:', label='',
             linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', label='', linewidth=3)

    # add figure names and x-/ylabels
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)
    add_at(ax3, r"c", loc=3)

    ax1.set_ylabel('Altitude (m)', fontsize=30)
    ax1.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax2.set_ylabel('Altitude (m)', fontsize=30)
    ax2.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax3.set_ylabel(r'Volume ($m^3$)', fontsize=30)
    ax3.set_xlabel('Time (years)', fontsize=30)

    ax1.tick_params(axis='both', which='major', labelsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=30)
    ax3.tick_params(axis='both', which='major', labelsize=30)
    ax3.yaxis.offsetText.set_fontsize(30)
    ax3.set_xlim(xmin=1847, xmax=2003)

    l1 = ax1.legend(loc=1, fontsize=30)
    l1.set_zorder(1)

    l2 = ax2.legend(loc=1, fontsize=30)
    l2.set_zorder(1)

    l3 = ax3.legend(loc=1, fontsize=30)
    l3.set_zorder(1)

    fig_name = 'median_'+str(ys)+'_'+gdir.rgi_id
    #plt.savefig(os.path.join(plot_dir, fig_name+'.pdf'), dpi=300)
    plt.savefig(os.path.join(plot_dir, fig_name+'.png'), dpi=300)
    #plt.show()
    plt.close()
    return median_model


def plot_median_vs_min(list, plot_dir):

    median = pd.concat(list, ignore_index=True)
    # median = median_df
    median['ex_mod'] = median['ex_p'].apply(lambda x: FileModel(x))
    median['diff_v'] = median.ex_mod.apply(
        lambda x: x.volume_km3_ts()[1850]) - median.m_mod.apply(
        lambda x: x.volume_km3_ts()[1850])
    median['diff2_v'] = median.ex_mod.apply(
        lambda x: x.volume_km3_ts()[1850]) - median.min_mod.apply(
        lambda x: x.volume_km3_ts()[1850])

    fig, ax = plt.subplots(figsize=(15, 10))
    boxprops = {'color': 'black', 'linewidth': 3}
    medianprops = {'color': 'C0', 'linewidth': 3}
    whiskerprops = {'linewidth': 3}
    capprops = {'linewidth': 3}
    flierprops = {'color': 'black', 'marker': '.', 'markerfacecolor': 'black',
                  'markersize': 10}

    median = median.drop(median.diff2_v.idxmin())
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                  alpha=0.7, lw=2)
    box = ax.boxplot(median[['diff_v', 'diff2_v']].values, boxprops=boxprops,
                     medianprops=medianprops,whiskerprops=whiskerprops,
                     capprops=capprops, flierprops=flierprops, widths=0.5,
                     patch_artist=True, labels=['median', 'minimum'])
    for i in range(len(box['boxes'])):
        box['boxes'][i].set_alpha(0.3)
        patch = patches.PathPatch(box['boxes'][i].get_path(), fill=False,
                                  edgecolor='black', lw=2)
        ax.add_patch(patch)
    ax.set_axisbelow(True)
    plt.ylabel(r'Differences in volume ($km^3$)')
    plt.tight_layout()
    #plt.savefig(os.path.join(plot_dir, 'boxplot.pdf'), dpi=300)
    plt.show()


def _adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def plot_abs_error_t0(list, ylabel, plot_dir):

    median = pd.concat(list,ignore_index=True)

    fig, ax = plt.subplots(figsize=(15, 10))
    boxprops = {'color': 'black', 'linewidth': 3}
    medianprops = {'color': 'C0', 'linewidth': 3}
    whiskerprops = {'linewidth': 3}
    capprops = {'linewidth': 3}
    flierprops = {'color': 'black', 'marker': '.', 'markerfacecolor': 'black',
                  'markersize': 5}

    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                  alpha=0.7, lw=2)

    ax.set_title('n = ' + str(len(median)))
    box = ax.boxplot(median.values, showfliers=False, boxprops=boxprops,
                     medianprops=medianprops,whiskerprops=whiskerprops,
                     capprops=capprops, flierprops=flierprops, widths=0.5,
                     patch_artist=True, labels=[1850, '', '', '', '', 1875, '',
                                                '', '', '', 1900, '', '', '',
                                                '', 1925, '', '', '', '', 1950,
                                                '', '', ''])
    for i in range(len(box['boxes'])):
        box['boxes'][i].set_alpha(0.3)
        patch = patches.PathPatch(box['boxes'][i].get_path(), fill=False,
                                  edgecolor='black', lw=2)
        ax.add_patch(patch)
    ax.set_axisbelow(True)
    # ax.violinplot(median.values, showextrema=False)
    plt.ylabel(ylabel)
    plt.xlabel(r'Starting time $t_0$')
    plt.tight_layout()

    # plt.savefig(plot_dir, dpi=300)
    plt.show()


def plot_error_t02(list, plot_dir, abs=True):

    merge = list[0].append(list[1])

    for rgi in merge.index:
        merge.loc[rgi, 'end'] = rgi.split('11.')[-1]
    to_delete = merge['end'].value_counts() > 1
    to_delete = to_delete[to_delete == True].index

    for rgi in list[1].index:
        if rgi.split('11.')[-1] in to_delete:
            list[1] = list[1].drop(rgi)
    median = list[0].append(list[1])

    fig, ax = plt.subplots(figsize=(10, 15))

    median = median.loc[:, median.columns.astype(float) % 10 == 0]

    pt.half_violinplot(data=median.values,split=True, label='', ax=ax,
                       linewidth=2, width=1.2, bw=.1, color='C0', alpha=0.1,
                       scale='area', inner=None, orient='h')

    sns.boxplot(data=median.values, width=0.2,
                boxprops={'facecolor': 'none', "zorder": 10, 'linewidth': 2,
                          'edgecolor': 'k'},
                medianprops={'color': 'C0', 'linewidth': 2},
                flierprops={'marker': '.', 'markerfacecolor': 'k',
                            'markeredgecolor': 'k', 'markersize': 8},
                capprops={'linewidth': 2, 'color': 'k'},
                whiskerprops={'linewidth': 2, 'color': 'k'}, orient='h',)

    ax = plt.gca()
    for art in ax.get_children():
        if isinstance(art, PolyCollection):

            art.set_facecolor(matplotlib.colors.colorConverter.to_rgba('C0', alpha=.5))
            art.set_edgecolor('k')

    ax.set_axisbelow(True)
    ax.set_yticklabels(median.columns)
    ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                  alpha=0.7, lw=2)

    plt.ylabel(r'Starting time $t_0$', fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    plt.tight_layout()

    if abs:
        plt.xlabel(r'Absolute error in $t_0$ ($km^3)$', fontsize=15)
        #ax.set_xlim((-0.25, 0.25))
    else:
        plt.xlabel(r'Logarithmic error in $t_0$', fontsize=15)
        ax.set_ylim(11.5, -1)
    plt.savefig(plot_dir, dpi=300)
    plt.show()


def plot_error_t03(list, plot_dir, gdirs, abs=True):

    merge = list[0].append(list[1])
    for rgi in merge.index:
        merge.loc[rgi, 'end'] = rgi.split('11.')[-1]
    to_delete = merge['end'].value_counts() > 1
    to_delete = to_delete[to_delete == True].index

    for rgi in list[1].index:
        if rgi.split('11.')[-1] in to_delete:
            list[1] = list[1].drop(rgi)
    median = list[0].append(list[1])

    fig, ax = plt.subplots(figsize=(10, 15))

    median = median.loc[:, median.columns.astype(float) % 10 == 0]
    for gdir in gdirs:
        if gdir.rgi_id in median.index:
            rp = gdir.get_filepath('model_run', filesuffix='experiment')
            ex_mod = FileModel(rp)
            print(ex_mod.volume_km3_ts()[1960])
    print(median)

    pt.half_violinplot(data=median[:30].values, split=True, label='', ax=ax,
                       linewidth=2, width=1.2, bw=.1, color='C0', alpha=0.1 ,
                       scale='area', inner=None, orient='h')

    sns.boxplot(data=median[:30].values, width=0.2,
                boxprops={'facecolor': 'none', "zorder": 10, 'linewidth': 2,
                          'edgecolor': 'k'},
                medianprops={'color': 'C0', 'linewidth': 2},
                flierprops={'marker': '.', 'markerfacecolor': 'k',
                            'markeredgecolor': 'k', 'markersize': 8},
                capprops={'linewidth': 2, 'color': 'k'},
                whiskerprops={'linewidth': 2, 'color': 'k'}, orient='h',)

    ax = plt.gca()
    for art in ax.get_children():
        if isinstance(art, PolyCollection):

            art.set_facecolor(matplotlib.colors.colorConverter.to_rgba('C0',alpha=.5))
            art.set_edgecolor('k')

    pt.half_violinplot(data=median[30:-1].values, split=True, label='', ax=ax,
                       linewidth=2, width=1.2,
                       bw=.1, color='C1', alpha=0.1, scale='area', inner=None,
                       orient='h')
    sns.boxplot(data=median[30:-1].values, width=0.2,
                boxprops={'facecolor': 'none', "zorder": 10, 'linewidth': 2,
                          'edgecolor': 'C1'},
                medianprops={'color': 'C1', 'linewidth': 2},
                flierprops={'marker': '.', 'markerfacecolor': 'C1',
                            'markeredgecolor': 'C1', 'markersize': 8},
                capprops={'linewidth': 2, 'color': 'C1'},
                whiskerprops={'linewidth': 2, 'color': 'C1'}, orient='h',
                )
    for art in ax.get_children():
        if isinstance(art, PolyCollection):
            art.set_facecolor(matplotlib.colors.colorConverter.to_rgba(art.get_facecolor()[0][:-1], alpha=.5))
            art.set_edgecolor('k')
    ax.set_axisbelow(True)

    ax.set_yticklabels(median.columns[:-1])
    ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                  alpha=0.7, lw=2)

    plt.ylabel(r'Starting time $t_0$', fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    plt.tight_layout()

    if abs:
        plt.xlabel(r'Absolute error in $t_0$ ($km^3)$', fontsize=15)
        ax.set_xlim((-0.25, 0.25))
    else:
        plt.xlabel(r'Logarithmic error in $t_0$', fontsize=15)
        ax.set_ylim(11.5, -1)
    #plt.savefig(plot_dir, dpi=300)
    plt.show()


def plot_min_vs_med(list, ylabel, plot_dir):

    merge = list[0].append(list[1])
    '''
    for rgi in merge.index:
        merge.loc[rgi, 'end'] = rgi.split('11.')[-1]
    to_delete = merge['end'].value_counts() > 1
    to_delete = to_delete[to_delete == True].index

    for rgi in list[1].index:
        if rgi.split('11.')[-1] in to_delete:
            list[1] = list[1].drop(rgi)
    '''
    median = list[0].append(list[1])
    print(median)


    fig, ax = plt.subplots(figsize=(10, 8))
    # median = median.loc[:,median.columns %10==0]

    pt.half_violinplot(data=median.values, split=True, label='', ax=ax,
                       linewidth=2, width=1, orient='h', bw=.075, color='C0',
                       alpha=0.1, scale='area', inner=None)

    sns.boxplot(data=median.values, width=0.2, orient='h',
                boxprops={'facecolor': 'none', "zorder": 10, 'linewidth': 2,
                          'edgecolor': 'k'},
                medianprops={'color': 'C0', 'linewidth': 2},
                flierprops={'marker': '.', 'markerfacecolor': 'k',
                            'markeredgecolor': 'k', 'markersize': 8},
                capprops={'linewidth': 2, 'color': 'k'},
                whiskerprops={'linewidth': 2, 'color': 'k'},
                )

    ax = plt.gca()
    for art in ax.get_children():
        if isinstance(art, PolyCollection):
            art.set_facecolor(
                matplotlib.colors.colorConverter.to_rgba('C0', alpha=.5))
            art.set_edgecolor('k')

    ax = sns.stripplot(data=median.values, color='C0', edgecolor="white",
                       size=3, jitter=1, zorder=0,alpha=0.5,orient='h')

    ax.set_axisbelow(True)

    ax.set_yticklabels(['median', 'minimum'], rotation=90,
                       verticalalignment="center")
    ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                  alpha=0.7, lw=2)

    plt.xlabel('Logarithmic error in 1850', fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=25)
    ax.set_ylim((1.25, -0.75))
    # plt.tight_layout()
    #plt.savefig(plot_dir, dpi=300)
    plt.show()


def union_plot(png_list, plot_dir):
    l = len(png_list)
    fig, ax = plt.subplots(l, 1, figsize=(20, 30))
    char = ['a', 'b', 'c']

    for i, ax in enumerate(fig.axes):
        image = matplotlib.image.imread(png_list[i])
        ax.imshow(image)
        ax.axis('off')

    fig.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(os.path.join(plot_dir, 'starting2.png'), dpi=400)
    #plt.show()


def plot_median_error(df, years, plot_dir,mean=False):
    plt.figure(figsize=(20, 10))


    if mean:
        for yr in years:
            r = [random.uniform(-0.5, 0.5) for x in range(len(df))]

            x = np.zeros(len(df.index)) + yr + r
            plt.plot(x, df.loc[:, yr].values,
                     'o', color='k', markersize=3, zorder=3)

        y1 = df.mean().values + df.std().values
        y2 = df.mean().values - df.std().values
        plt.plot(years, df.mean().values, color='grey', label='mean')
        plt.fill_between(years, y1, y2, alpha=0.5, color='grey', label='std')

    else:
        plt.ylim(-0.1,0.1)

    y3 = df.quantile(q=0.75).values
    y4 = df.quantile(q=0.25).values
    plt.plot(years, df.median().values, color='C0', label='median')
    plt.fill_between(years, y3, y4, alpha=0.5, color='C0', label='IQR')

    plt.xlabel(r'Starting time $t_0$')
    plt.ylabel(r'absolute error (km$^3$)')
    plt.grid()
    plt.legend(loc=4)
    plt.savefig(plot_dir, dpi=300)
    plt.show()

def plot_log_median_error(df, years, plot_dir,mean=False):
    fig = plt.figure(figsize=(20, 10))

    norm = mpl.colors.LogNorm(vmin=df.exp_1850.min(), vmax=df.exp_1850.max())
    cmap = matplotlib.cm.get_cmap('viridis')
    df = df.replace([np.inf, -np.inf], np.nan).dropna()

    if mean:
        for yr in years:

            for i in range(len(df)):
                x = yr + random.uniform(-0.5, 0.5)
                c = cmap(norm(df.iloc[i]['exp_1850']))
                plt.plot(x, df.iloc[i][yr], 'o', color=c, markersize=3, zorder=3)

        df = df.drop(columns=['1850_min', 'exp_1850'])

        y1 = df.mean().values + df.std().values
        y2 = df.mean().values - df.std().values
        plt.plot(years, df.mean().values, color='grey', label='mean')
        plt.fill_between(years, y1, y2, alpha=0.5, color='grey', label='std')

    else:
        df = df.drop(columns=['1850_min', 'exp_1850'])
        #plt.ylim(-0.1,0.1)

    y3 = df.quantile(q=0.75).values
    y4 = df.quantile(q=0.25).values
    plt.plot(years, df.median().values, color='C0', label='median')
    plt.fill_between(years, y3, y4, alpha=0.5, color='C0', label='IQR')

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm)
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label(r'volume (km$^3$)', fontsize=30)

    plt.xlabel(r'Starting time $t_0$')
    plt.ylabel(r'logarithmic error')
    plt.grid()
    plt.legend(loc=2)
    #plt.savefig(plot_dir, dpi=300)
    plt.show()

def plot_relative_error(df1, df2, type, plot_dir,all=False):

    plot_dir = os.path.join(plot_dir, 'errors')
    utils.mkdir(plot_dir)

    df1 = df1.replace([np.inf, -np.inf], np.nan).dropna()
    df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()

    if type == 'rel':
        title = 'Relative error'
        name = 'relative'
        unit = ''
    elif type == 'abs':
        title = 'Absolute error'
        name = 'absolute'
        unit = '$(km^3)$'
    elif type == 'log':
        title = 'Logarithmic error'
        name = 'logarithmic'
        unit = ''

    plt.figure(figsize=(18, 10))
    plt.title(title+', n=' + str(len(df1.index)))
    df1.median().plot(zorder=4, label='median state')
    plt.fill_between(df1.columns.values, df1.quantile(0.75),
                     df1.quantile(0.25),
                     facecolor='C0', alpha=0.5, interpolate=True)

    df2.median().plot(zorder=4, label='minimum state')
    plt.fill_between(df2.columns.values,
                     df2.quantile(0.75),
                     df2.quantile(0.25),
                     facecolor='C1', alpha=0.5, interpolate=True)

    plt.grid()
    plt.xlabel('Time')
    plt.ylabel(r'Median with interquartile range ' + unit)
    plt.legend(loc='best')
    plt.savefig(os.path.join(plot_dir, name+'_median.png'), dpi=300)
    plt.savefig(os.path.join(plot_dir, name + '_median.pdf'), dpi=300)

    df1 = df1.replace([np.inf, -np.inf], np.nan).dropna(axis=1)
    df2 = df2.replace([np.inf, -np.inf], np.nan).dropna(axis=1)

    plt.figure(figsize=(20, 10))
    plt.title(title+' , n=' + str(len(df1.index)))
    df1.mean().plot(zorder=4, label='median state')
    plt.fill_between(df1.columns.values,
                     df1.mean() + df1.std(),
                     df1.mean() - df1.std(),
                     facecolor='C0', alpha=0.5, interpolate=True)
    df2.mean().plot(zorder=4, label='minimum state')
    plt.fill_between(df2.columns.values,
                     df2.mean() + df2.std(),
                     df2.mean() - df2.std(),
                     facecolor='C1', alpha=0.5, interpolate=True)
    plt.grid()
    plt.xlabel('Time')
    plt.ylabel(r'Mean with standard deviation ' + unit)
    plt.legend(loc='best')
    plt.savefig(os.path.join(plot_dir, name+'_mean.png'), dpi=300)
    plt.savefig(os.path.join(plot_dir, name + '_mean.pdf'), dpi=300)
    #plt.show()

    if all:
        plt.figure(figsize=(18, 10))
        plt.title(title + ', n=' + str(len(df1.index)))
        df1.median().plot(zorder=4, label='median state')
        plt.fill_between(df1.columns.values, df1.quantile(0.75),
                         df1.quantile(0.25),
                         facecolor='C0', alpha=0.5, interpolate=True)

        plt.fill_between(df1.columns.values,
                         df1.max(),
                         df1.min(),
                         facecolor='C0', alpha=0.2, interpolate=True, zorder=1)

        plt.fill_between(df2.columns.values,
                         df2.max(),
                         df2.min(),
                         facecolor='C1', alpha=0.2, interpolate=True, zorder=1)

        for yr in df1.columns.values:
            if yr % 10 == 0:

                y1 = df1.loc[:, yr].values
                y2 = df2.loc[:, yr].values

                x1 = np.zeros(len(y1)) + yr - [random.uniform(0.5, 1.5) for x
                                               in range(len(y1))]
                x2 = np.zeros(len(y2)) + yr + [random.uniform(0.5, 1.5) for x
                                               in range(len(y2))]
                plt.plot(x1, y1, 'o', color='C0', markersize=2)
                plt.plot(x2, y2, 'o', color='C1', markersize=2)

        df2.median().plot(zorder=4, label='minimum state')
        plt.fill_between(df2.columns.values,
                         df2.quantile(0.75),
                         df2.quantile(0.25),
                         facecolor='C1', alpha=0.5, interpolate=True)

        plt.grid()
        plt.xlabel('Time')
        plt.ylabel(r'Median with interquartile range ' + unit)
        plt.legend(loc='best')
        plt.savefig(os.path.join(plot_dir, name + '_median_all.png'), dpi=300)
        plt.savefig(os.path.join(plot_dir, name + '_median_all.pdf'), dpi=300)

def plot_starting_error(df1, df2,df3, type, plot_dir):

    df1 = df1.replace([np.inf, -np.inf], np.nan).dropna(axis=1)
    df2 = df2.replace([np.inf, -np.inf], np.nan).dropna(axis=1)

    if type == 'rel':
        title = 'Relative error'
        name = 'relative'
        unit = ''
    elif type == 'abs':
        title = 'Absolute error'
        name = 'absolute'
        unit = '$(km^3)$'
    elif type == 'log':
        title = 'Logarithmic error'
        name = 'logarithmic'
        unit = ''

    norm = mpl.colors.LogNorm(vmin=0.001, vmax=3)
    cmap1 = matplotlib.cm.get_cmap('Blues_r')
    cmap2 = matplotlib.cm.get_cmap('Oranges_r')

    fig = plt.figure(figsize=(18, 10))
    plt.title(title + ', n=' + str(len(df1.columns)))
    df1.median(axis=1).plot(zorder=4, label='median state')
    plt.fill_between(df1.index.values, df1.quantile(0.75, axis=1),
                     df1.quantile(0.25, axis=1),
                     facecolor='C0', alpha=0.5, interpolate=True)

    plt.fill_between(df1.index.values,
                     df1.max(axis=1),
                     df1.min(axis=1),
                     facecolor='C0', alpha=0.2, interpolate=True, zorder=1)

    for yr in df1.index.values:
        if yr % 10 == 0:
            y1 = df1.loc[yr, :].values
            y2 = df2.loc[yr, :].values

            x1 = np.zeros(len(y1)) + yr - [random.uniform(0.5, 1.5) for x in range(len(y1))]
            x2 = np.zeros(len(y2)) + yr + [random.uniform(0.5, 1.5) for x in range(len(y2))]
            plt.plot(x1, y1, 'o', color='C0', markersize=2)
            plt.plot(x2, y2, 'o', color='C1', markersize=2)

    df2.median(axis=1).plot(zorder=4, label='minimum state')
    plt.fill_between(df2.index.values,
                     df2.quantile(0.75, axis=1),
                     df2.quantile(0.25, axis=1),
                     facecolor='C1', alpha=0.5, interpolate=True)

    plt.fill_between(df2.index.values,
                     df2.max(axis=1),
                     df2.min(axis=1),
                     facecolor='C1', alpha=0.2, interpolate=True, zorder=1)


    plt.grid()
    plt.xlabel('Time')
    plt.ylabel(r'Median with interquartile range ' + unit)
    plt.legend(loc='best')
    plt.savefig(os.path.join(plot_dir, name + '_median_all.png'), dpi=300)

    df1 = df1.replace([np.inf, -np.inf], np.nan).dropna(axis=1)
    df2 = df2.replace([np.inf, -np.inf], np.nan).dropna(axis=1)

    plt.figure(figsize=(15, 10))
    plt.title(title + ' , n=' + str(len(df1.columns)))
    df1.mean(axis=1).plot(zorder=4, label='median state')
    plt.fill_between(df1.index.values,
                     df1.mean(axis=1) + df1.std(axis=1),
                     df1.mean(axis=1) - df1.std(axis=1),
                     facecolor='C0', alpha=0.5, interpolate=True)
    df2.mean(axis=1).plot(zorder=4, label='minimum state')
    plt.fill_between(df2.index.values,
                     df2.mean(axis=1) + df2.std(axis=1),
                     df2.mean(axis=1) - df2.std(axis=1),
                     facecolor='C1', alpha=0.5, interpolate=True)
    plt.grid()
    plt.xlabel('Time')
    plt.ylabel(r'Mean with standard deviation ' + unit)
    plt.legend(loc='best')
    plt.savefig(os.path.join(plot_dir, name + '_mean.png'),
                dpi=300)
    plt.show()