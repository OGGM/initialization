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
import matplotlib.patches as patches
#from sklearn import preprocessing


import os


mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] = 30
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize']= 25 #30
mpl.rcParams['lines.linewidth']=3

def add_at(ax, t, loc=2):
    fp = dict(size=30)
    _at = AnchoredText(t, loc=loc, prop=fp,borderpad=0)
    _at.patch.set_linewidth(3)
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

def plot_experiment(gdir,df,ex_mod,ys, plot_dir):
    plot_dir = os.path.join(plot_dir,'00_experiment')
    utils.mkdir(plot_dir)

    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[-1].map_dx
    fig = plt.figure(figsize=(15,14))
    grid = plt.GridSpec(2, 1, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[1,0],sharex=ax1)

    if gdir.name != '':
        ax1.set_title(gdir.rgi_id+':'+gdir.name,fontsize=30)
    else:
        ax1.set_title(gdir.rgi_id+': Guslarferner', fontsize=30)

    # plot experiments
    ex_mod = deepcopy(ex_mod)
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)
    i = np.where(ex_mod.fls[-1].thick > 0)[0][-1] +10

    ax1.plot(x[:i], ex_mod.fls[-1].surface_h[:i], 'k:',label=r'$x_{'+str(ys)+'}^{exp}$',linewidth=3)
    ax1.plot(x[:i], ex_mod.fls[-1].bed_h[:i], 'k',label=r'$b$', linewidth=3)

    ex_mod.run_until(2000)

    ax2.plot(x[:i], ex_mod.fls[-1].surface_h[:i], 'k:',label=r'$x_{2000}^{exp = obs} $', linewidth=3)
    ax2.plot(x[:i], ex_mod.fls[-1].bed_h[:i], 'k', label=r'$b$',linewidth=3)

    # add figure names and legends
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)

    ax1.legend(loc=1)
    ax2.legend(loc=1)

    ax1.set_ylabel('Altitude (m)', fontsize=30)
    ax1.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax2.set_ylabel('Altitude (m)',fontsize=30)
    ax2.set_xlabel('Distance along the main flowline (m)',fontsize=30)

    ax1.tick_params(axis='both', which='major', labelsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=30)

    plt.savefig(os.path.join(plot_dir,'experiment_'+str(ys)+'_'+gdir.rgi_id+'.pdf'), dpi=300)
    plt.savefig(os.path.join(plot_dir, 'experiment_' + str(ys) + '_' + gdir.rgi_id + '.png'), dpi=300)
    plt.close()

def plot_candidates(gdir, df, ex_mod, yr, step,plot_dir):
    plot_dir = os.path.join(plot_dir,'06_candidates',gdir.rgi_id)
    utils.mkdir(plot_dir)
    fig, ax = plt.subplots(figsize=(10,10))
    for file in os.listdir(gdir.dir):
        if file.startswith('model_run'+str(yr)+'_random'):
            suffix = file.split('model_run')[1].split('.nc')[0]
            rp = gdir.get_filepath('model_run', filesuffix=suffix)
            try:
                fmod = FileModel(rp)
                fmod.volume_m3_ts().plot(ax=ax, color='grey', label='', zorder=1)

            except:
                pass

    # last one again for labeling
    label = r'temperature bias $\in [$' + str(
        df['temp_bias'].min()) + ',' + str(df['temp_bias'].max()) + '$]$'

    t_eq = df['time'].sort_values().iloc[0]

    df['Fitness value'] = df.objective
    df.time = df.time.apply(lambda x: int(x))

    plt.title(gdir.rgi_id)

    if step == 'step1':
        fmod.volume_m3_ts().plot(ax=ax, color='grey', label=label, zorder=1)
        plt.legend(loc=0, fontsize=23)
        plt.xlabel('Time (years)')
        plt.ylabel(r'Volume $(m^3)$')
        plt.savefig(os.path.join(plot_dir, 'candidates1_' + str(yr) + '_' + str(
                gdir.rgi_id) + '.png'), dpi=300)
    elif step == 'step2':
        ax.axvline(x=int(t_eq), color='k', zorder=1, label=r'$t_{stag}$')
        fmod.volume_m3_ts().plot(ax=ax, color='grey', label='', zorder=1)
        # black points
        df.plot.scatter(x='time', y='volume', ax=ax, color='k', label='candidates', s=250, zorder=2)
        plt.legend(loc=0, fontsize=23)
        plt.xlabel('Time (years)')
        plt.ylabel(r'Volume $(m^3)$')
        plt.savefig(os.path.join(plot_dir, 'candidates2_' + str(yr) + '_' + str(
            gdir.rgi_id) + '.png'), dpi=300)
    elif step == 'step3':
        fmod.volume_m3_ts().plot(ax=ax, color='grey', label=None, zorder=1)
        ax.axvline(x=int(t_eq), color='k', zorder=1)
        # colored points
        df.plot.scatter(x='time', y='volume', ax=ax, c='Fitness value', colormap='viridis', norm=mpl.colors.LogNorm(vmin=0.1, vmax=1e5), s=250,edgecolors='k', zorder=2)
        plt.xlabel('Time (years)')
        plt.ylabel(r'Volume $(m^3)$')
        plt.savefig(os.path.join(plot_dir, 'candidates3_' + str(yr) + '_' + str(gdir.rgi_id) + '.png'), dpi=300)

    plt.close()

    plt.figure(figsize=(15,14))
    plt.hist(df.volume.values, bins=20)
    plt.xlabel(r'Volume $(m^3)$')
    plt.ylabel(r'Frequency')
    plt.title(gdir.rgi_id)
    plt.savefig(os.path.join(plot_dir, 'hist_candidates' + str(yr) + '_' + str(gdir.rgi_id) + '.png'), dpi=300)
    plt.close()



def plot_compare_fitness(gdir,df,ex_mod,ys, plot_dir):

    plot_dir = os.path.join(plot_dir, '05_compare_fitness')
    utils.mkdir(plot_dir)
    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[-1].map_dx

    fig = plt.figure(figsize=(15,14))

    grid = plt.GridSpec(3, 1, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[1, 0],sharex=ax1)
    ax3 = plt.subplot(grid[2, 0],sharex=ax1)

    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)

    if gdir.name != '':
        ax1.set_title(gdir.rgi_id+':'+gdir.name,fontsize=30)
    else:
        ax1.set_title(gdir.rgi_id, fontsize=30)

    df['objective2'] = df.model.apply(lambda x: abs(x.area_m2_ts()[2000] - ex_mod.area_m2_ts()[2000])**2)
    df['objective3'] = df.model.apply(lambda x: abs(x.length_m_ts()[2000] - ex_mod.length_m_ts()[2000])**2)
    norm = mpl.colors.LogNorm(vmin=df.objective.min(), vmax=df.objective.max(),clip=True)
    norm2 = mpl.colors.LogNorm(vmin=df.objective2.min()+1,vmax=df.objective2.max(),clip=True)
    norm3 = mpl.colors.LogNorm(vmin=df.objective3.min() + 1, vmax=df.objective3.max(), clip=True)
    cmap = matplotlib.cm.get_cmap('viridis')

    df = df.sort_values('objective', ascending=False)
    for i, model in df['model'].iteritems():
        if df.objective.min() != df.objective.max():
            color = cmap(norm(df.loc[i, 'objective']))
        else:
            color = cmap(df.loc[i, 'objective'])
        model.volume_m3_ts().plot(ax=ax1, color=[color],
                                  linewidth=3,
                                  label=r'$s_{1850-2000}^{exp}$')
    df = df.sort_values('objective2', ascending=False)
    for i, model in df['model'].iteritems():
        if df.objective2.min() != df.objective2.max():
            color = cmap(norm2(int(df.loc[i, 'objective2'])))
        else:
            color = cmap(df.loc[i, 'objective2'])

        model.volume_m3_ts().plot(ax=ax2, color=[color],
                                   linewidth=3,
                                   label=r'$s_{1850-2000}^{exp}$')
    df = df.sort_values('objective3', ascending=False)
    for i, model in df['model'].iteritems():
        if df.objective3.min() != df.objective3.max():
            color = cmap(norm3(int(df.loc[i, 'objective3'])))
        else:
            color = cmap(df.loc[i, 'objective3'])
        model.volume_m3_ts().plot(ax=ax3, color=[color],
                                   linewidth=3,
                                   label=r'$s_{1850-2000}^{exp}$')
        #ax3.plot(x, deepcopy(model.fls[-1].surface_h), color=color,
        #         label='')

    # plot experiments
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_m3_ts().plot(ax=ax1, color='k', linestyle=':', linewidth=3,
                               label=r'$s_{1850-2000}^{exp}$')
    ex_mod.volume_m3_ts().plot(ax=ax2, color='k', linestyle=':', linewidth=3,
                               label=r'$s_{1850-2000}^{exp}$')
    ex_mod.volume_m3_ts().plot(ax=ax3, color='k', linestyle=':', linewidth=3, label=r'$s_{1850-2000}^{exp}$')
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)
    '''
    ax1.plot(x, ex_mod.fls[-1].surface_h, 'k:',label=r'$x_{1850}^{exp}$',linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k',label=r'$b$', linewidth=3)

    ax2.plot(x, ex_mod.fls[-1].surface_h, 'k:', label=r'$x_{1850}^{exp}$',
            linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', label=r'$b$', linewidth=3)

    ax3.plot(x, ex_mod.fls[-1].surface_h, 'k:', label=r'$x_{1850}^{exp}$',
             linewidth=3)
    ax3.plot(x, ex_mod.fls[-1].bed_h, 'k', label=r'$b$', linewidth=3)
    '''
    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax1,ax2,ax3])
    cbar = fig.colorbar(sm, cax=cax ,**kw)
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label('Normalized fitness value', fontsize=30)


    # add figure names and legends
    add_at(ax1, r"a", loc=1)
    add_at(ax2, r"b", loc=1)
    add_at(ax3, r"c", loc=1)

    #ax1.legend(loc=1)
    #ax2.legend(loc=1)
    #ax3.legend(loc=1)

    ax2.set_ylabel(r'Volume ($km^3$)', fontsize=30)
    #ax1.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    #ax2.set_ylabel('Altitude (m)',fontsize=30)
    #ax2.set_xlabel('Distance along the main flowline (m)',fontsize=30)
    #ax3.set_xlabel('Distance along the main flowline (m)', fontsize=30)

    ax1.tick_params(axis='both', which='major', labelsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=30)
    #ax3.tick_params(axis='both', which='major', labelsize=30)
    #ax3.yaxis.offsetText.set_fontsize(30)
    plt.show()
    #plt.savefig(os.path.join(plot_dir,gdir.rgi_id+'.pdf'), dpi=300)
    '''
    df['diff'] = df.objective.apply(lambda x: norm(x)) - df.objective2.apply(
        lambda x: norm2(x))
    plt.figure(figsize=(15, 14))
    m = max(abs(df['diff'].min()), df['diff'].max())
    norm3 =mpl.colors.Normalize(vmin=-m, vmax=m,clip=True)
    cmap2 = matplotlib.cm.get_cmap('RdBu')

    for i, model in df['model'].iteritems():
        model = deepcopy(model)
        model.reset_y0(ys)

        # calculate objective2
        color = cmap2(norm3(df.loc[i, 'diff']))
        plt.plot(x, deepcopy(model.fls[-1].surface_h), color=color,
                 label='')

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap2, norm=norm3)
    sm.set_array([])
    cbar = plt.colorbar(sm, **kw)
    cbar.set_label('Fitness value', fontsize=30)
    '''
    #plt.show()
    plt.close()

def plot_fitness_over_time2(gdir, df_list, ex_mod, plot_dir):
    from matplotlib.patches import Rectangle, Circle
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(111, yscale='linear')
    norm = mpl.colors.LogNorm(vmin=0.1, vmax=1e5)
    cmap = matplotlib.cm.get_cmap('viridis')

    volumes = np.linspace(df_list['1850'].volume.min(),
                          df_list['1850'].volume.max(), 100)

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
                Rectangle((y - 2.5, x), 5, volumes[1] - volumes[0],
                          color=color))

    # add experiment in plot
    ex_mod.volume_m3_ts().plot(ax=ax, linestyle=':', color='gold')
    ax.set_xlim(1847.5, 1967.5)
    ax.set_ylim(volumes[0], volumes[-1])
    plt.title(gdir.rgi_id + ': ' + gdir.name, fontsize=25)
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
    #plt.show()


def plot_col_fitness(gdir,df,ex_mod,ys, plot_dir):

    plot_dir = os.path.join(plot_dir,'03_surface_by_fitness', gdir.rgi_id)
    utils.mkdir(plot_dir)
    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[-1].map_dx
    fig = plt.figure(figsize=(25,18))
    grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1],sharey=ax1)
    ax3 = plt.subplot(grid[1, :])

    if gdir.name != '':
        plt.suptitle(gdir.rgi_id+': '+gdir.name,fontsize=30)
    else:
        plt.suptitle(gdir.rgi_id+': Guslarferner', fontsize=30)

    norm = mpl.colors.LogNorm(vmin=0.1, vmax=1e5)
    cmap = matplotlib.cm.get_cmap('viridis')

    df = df.sort_values('objective', ascending=False)

    for i, model in df['model'].iteritems():
        model = deepcopy(model)
        model.reset_y0(ys)
        color = cmap(norm(df.loc[i, 'objective']))

        ax1.plot(x, deepcopy(model.fls[-1].surface_h), color=color,
                 label='')
        model.volume_m3_ts().plot(ax=ax3, color=[color], label='')
        model.run_until(2000)

        ax2.plot(x, model.fls[-1].surface_h, color=color,label='')

    # plot experiments
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_m3_ts().plot(ax=ax3, color='gold', linestyle=':', linewidth=3, label='')
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)

    ax1.plot(x, ex_mod.fls[-1].surface_h, ':', color='gold',label='',linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k',label='', linewidth=3)

    ex_mod.run_until(2000)

    ax2.plot(x, ex_mod.fls[-1].surface_h, ':', color='gold',label='', linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', label='',linewidth=3)

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax1, ax2, ax3])
    cbar = fig.colorbar(sm, cax=cax ,**kw)
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label('Fitness value', fontsize=30)

    # add figure names and x-/ylabels
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)
    add_at(ax3, r"c", loc=3)

    ax1.set_ylabel('Altitude (m)', fontsize=30)
    ax1.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax2.set_ylabel('Altitude (m)',fontsize=30)
    ax2.set_xlabel('Distance along the main flowline (m)',fontsize=30)
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
    lc3 = LineCollection(segments, color='gold', linestyle=':',
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
                    labels=[r'$s_{'+str(ys)+'-2000}^{exp}$', r'$s_{'+str(ys)+'-2000}$'], loc=1)

    l1.set_zorder(0)
    l2.set_zorder(0)
    l3.set_zorder(0)


    ax3.set_xlim(xmin=1847, xmax=2003)
    plt.savefig(os.path.join(plot_dir,'surface_'+str(ys)+'_'+gdir.rgi_id+'.pdf'), dpi=300)
    plt.savefig(os.path.join(plot_dir, 'surface_' + str(ys) + '_' + gdir.rgi_id + '.png'), dpi=300)
    plt.close()


def plot_col_fitness2(gdir,df,ex_mod,ys, plot_dir):

    plot_dir = os.path.join(plot_dir,'03_surface_by_fitness')
    utils.mkdir(plot_dir)
    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[-1].map_dx
    fig = plt.figure(figsize=(25,18))
    grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1],sharey=ax1)
    ax3 = plt.subplot(grid[1, :])

    if gdir.name != '':
        plt.suptitle(gdir.rgi_id+': '+gdir.name,fontsize=30)
    else:
        plt.suptitle(gdir.rgi_id+': Guslarferner', fontsize=30)

    norm = mpl.colors.LogNorm(vmin=0.1, vmax=1e5)
    cmap = matplotlib.cm.get_cmap('viridis')

    df = df.sort_values('objective', ascending=False)

    for i, model in df['model'].iteritems():
        model = deepcopy(model)
        model.reset_y0(ys)
        color = cmap(norm(df.loc[i, 'objective']))

        ax1.plot(x, deepcopy(model.fls[-1].surface_h), color=color,
                 label='')
        model.volume_m3_ts().plot(ax=ax3, color=[color], label='')
        model.run_until(2000)

        ax2.plot(x, model.fls[-1].surface_h, color=color,label='')

    # plot experiments
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_m3_ts().plot(ax=ax3, color='gold', linestyle=':', linewidth=3, label='')
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)

    ax1.plot(x, ex_mod.fls[-1].surface_h, ':', color='gold',label='',linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k',label='', linewidth=3)

    ex_mod.run_until(2000)

    ax2.plot(x, ex_mod.fls[-1].surface_h, ':', color='gold',label='', linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', label='',linewidth=3)

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax1, ax2, ax3])
    cbar = fig.colorbar(sm, cax=cax ,**kw)
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label('Fitness value', fontsize=30)

    # add figure names and x-/ylabels
    add_at(ax1, r"a", loc=3)
    add_at(ax2, r"b", loc=3)
    add_at(ax3, r"c", loc=3)

    ax1.set_ylabel('Altitude (m)', fontsize=30)
    ax1.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax2.set_ylabel('Altitude (m)',fontsize=30)
    ax2.set_xlabel('Distance along the main flowline (m)',fontsize=30)
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
    lc3 = LineCollection(segments, color='gold', linestyle=':',
                         norm=plt.Normalize(0, 10), linewidth=3)

    l1 = ax1.legend(handles=[lc3, lc, lc2], handler_map={
        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=['true solution', 'reconstructed ice surface',
                            'bed topography'], loc=1,fontsize=22)

    l2 = ax2.legend(handles=[lc3, lc, lc2],
                    handler_map={
                        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=['true solution', 'reconstructed ice surface',
                            'bed topography'], loc=1, fontsize=22)

    l3 = ax3.legend(handles=[lc3],
                    handler_map={
                        lc: HandlerColorLineCollection(numpoints=100)},
                    labels=['true solution'], loc=1, fontsize=22)
    #ax3.legend(handles=[lc3], labels=['true solution'], loc=0)

    l1.set_zorder(0)
    l2.set_zorder(0)
    l3.set_zorder(0)


    ax3.set_xlim(xmin=ys-3, xmax=2003)
    plt.savefig(os.path.join(plot_dir,'01_Kesselwandferner.pdf'), dpi=300)
    plt.show()

def plot_median(gdir,df,ex_mod,ys, plot_dir):
    plot_dir = os.path.join(plot_dir,'04_median', gdir.rgi_id)
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
    else:
        plt.suptitle(gdir.rgi_id, fontsize=30)

    df = df.sort_values('objective', ascending=False)
    df = df[df['objective']<100]
    min_id = df.volume.idxmin()
    max_id = df.volume.idxmax()

    min_model = deepcopy(df.loc[min_id,'model'])
    max_model = deepcopy(df.loc[max_id, 'model'])

    min_model.reset_y0(ys)
    max_model.reset_y0(ys)
    ax1.fill_between(x, deepcopy(min_model.fls[-1].surface_h),
                     deepcopy(max_model.fls[-1].surface_h),alpha=0.3,
                     color='grey',label=r'$\mathcal{S}_{'+str(ys)+'}^{100}$')

    ax3.fill_between(min_model.volume_m3_ts().index,min_model.volume_m3_ts(),
                     max_model.volume_m3_ts(),alpha=0.3, color='grey',
                     label=r'$\mathcal{S}_{'+str(ys)+'-2000}^{100}$')
    min_model.run_until(2000)
    max_model.run_until(2000)
    ax2.fill_between(x, deepcopy(min_model.fls[-1].surface_h),
                     deepcopy(max_model.fls[-1].surface_h), alpha=0.3,
                     color='grey', label=r'$\mathcal{S}_{2000}^{100}$')
    # 5% quantile and median of 5% quantile
    quant_df = df[df.objective <= df.objective.quantile(0.05)]

    for model in quant_df.model:
        ax1.plot(x,model.fls[-1].surface_h,color='grey')
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
                     q_max_model.volume_m3_ts(), alpha=0.5,linewidth=3,
                     label=r'$Q_{0.05}(\mathcal{S}_{'+str(ys)+'-2000}^{100})$')
    q_min_model.run_until(2000)
    q_max_model.run_until(2000)
    ax2.fill_between(x, deepcopy(q_min_model.fls[-1].surface_h),
                     deepcopy(q_max_model.fls[-1].surface_h), alpha=0.5,
                     label=r'$Q_{0.05}(\mathcal{S}_{2000}^{100})$')

    # median of 5% quantile
    quant_df['length'] = quant_df.model.apply(lambda x: x.length_m)
    quant_df = quant_df.sort_values('length', ascending=False)
    l = len(quant_df)
    if l % 2:
        index = int((l - 1) / 2)
    else:
        index = int(l / 2)

    median_model = deepcopy(quant_df.iloc[index].model)
    median_model.volume_m3_ts().plot(ax=ax3,
                               linewidth=3, label='median')
    median_model.reset_y0(1850)
    median_model.run_until(ys)

    ax1.plot(x, median_model.fls[-1].surface_h, label='median',
             linewidth=3)
    median_model.run_until(2000)
    ax2.plot(x, median_model.fls[-1].surface_h, label='median',
             linewidth=3)

    #experiment
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

    l1 = ax1.legend(loc=1,fontsize=30)
    l1.set_zorder(1)

    l2 = ax2.legend(loc=1, fontsize=30)
    l2.set_zorder(1)

    l3 = ax3.legend(loc=1, fontsize=30)
    l3.set_zorder(1)

    plt.savefig(os.path.join(plot_dir,'median_'+str(ys)+'_'+gdir.rgi_id+'.pdf'), dpi=300)
    plt.savefig(os.path.join(plot_dir, 'median_' + str(ys)+'_'+gdir.rgi_id + '.png'), dpi=300)
    plt.close()
    return median_model


def plot_median_vs_min(list,plot_dir):

    median = pd.concat(list,ignore_index=True)
    # median = median_df
    median['ex_mod'] = median['ex_p'].apply(lambda x: FileModel(x))
    median['diff_v'] = median.ex_mod.apply(
        lambda x: x.volume_km3_ts()[1850]) - median.m_mod.apply(
        lambda x: x.volume_km3_ts()[1850])
    median['diff2_v'] = median.ex_mod.apply(
        lambda x: x.volume_km3_ts()[1850]) - median.min_mod.apply(
        lambda x: x.volume_km3_ts()[1850])

    fig, ax = plt.subplots(figsize=(15,10))
    boxprops = {'color': 'black', 'linewidth':3}
    medianprops = {'color': 'C0', 'linewidth': 3}
    whiskerprops = {'linewidth':3}
    capprops = {'linewidth': 3}
    flierprops = {'color': 'black','marker':'.', 'markerfacecolor':'black',
                  'markersize':10}

    median = median.drop(median.diff2_v.idxmin())
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                  alpha=0.7, lw=2)
    box = ax.boxplot(median[['diff_v','diff2_v']].values,boxprops=boxprops,
               medianprops=medianprops,whiskerprops=whiskerprops,
               capprops=capprops,flierprops=flierprops,widths=0.5,
               patch_artist=True, labels=['median','minimum'])
    for i in range(len(box['boxes'])):
        box['boxes'][i].set_alpha(0.3)
        patch = patches.PathPatch(box['boxes'][i].get_path(), fill=False,edgecolor='black', lw=2)
        ax.add_patch(patch)
    ax.set_axisbelow(True)
    plt.ylabel(r'Differences in volume ($km^3$)')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'boxplot.pdf'), dpi=300)
    plt.show()
