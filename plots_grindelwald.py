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
from core import objective_value

import os
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.legend_handler import HandlerLineCollection

mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] = 30
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize']=11 #25
mpl.rcParams['lines.linewidth']=3

def plot_surface_sep(gdir,df,ex_mod,ys, plot_dir):


    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[
        -1].map_dx

    fig1,ax1 = plt.subplots(figsize=(14,11))
    fig2,ax2 = plt.subplots(figsize=(14,11))
    fig3,ax3 = plt.subplots(figsize=(14,11))

    p2 = ax2.get_position()

    if gdir.name != '':
        fig1.suptitle(gdir.rgi_id + ':' + gdir.name, x=p2.x1 / 2, fontsize=30)
    else:
        fig1.suptitle(gdir.rgi_id, x=p2.x1 / 2, fontsize=30)

    import matplotlib as mpl
    import matplotlib.cm as cm

    norm = mpl.colors.LogNorm(vmin=0.1, vmax=1e5)
    # cmap = matplotlib.cm.get_cmap('RdYlGn_r')
    cmap = matplotlib.cm.get_cmap('RdBu_r')
    df = df.sort_values('objective', ascending=False)
    for i, model in df['model'].iteritems():
        model = deepcopy(model)
        model.reset_y0(ys)
        color = cmap(norm(df.loc[i, 'objective']))
        ax1.plot(x, deepcopy(model.fls[-1].surface_h), color=color,
                 linewidth=2)
        model.volume_m3_ts().plot(ax=ax3, color=color, linewidth=2)

        model.run_until(2000)

        ax2.plot(x, model.fls[-1].surface_h, color=color, linewidth=2)

    # plot experiments
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_m3_ts().plot(ax=ax3, color='k', linestyle=':', linewidth=3)
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)

    ax1.plot(x, ex_mod.fls[-1].surface_h, 'k:', linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k', linewidth=3)

    ex_mod.run_until(2000)

    ax2.plot(x, ex_mod.fls[-1].surface_h, 'k:', linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', linewidth=3)

    ax1.annotate(r'$t =  ' + str(ys) + '$', xy=(0.8, 0.9),
                 xycoords='axes fraction',
                 fontsize=30)
    ax2.annotate(r'$t =  2000$', xy=(0.8, 0.9), xycoords='axes fraction',
                 fontsize=30)

    ax1.set_ylabel('Altitude (m)', fontsize=30)
    ax1.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax2.set_ylabel('Altitude (m)', fontsize=30)
    ax2.set_xlabel('Distance along the main flowline (m)', fontsize=30)
    ax3.set_ylabel(r'Volume ($m^3$)', fontsize=30)
    ax3.set_xlabel('Time (years)', fontsize=30)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    #cax, kw = mpl.colorbar.make_axes([ax1, ax2, ax3])
    cbar = fig3.colorbar(sm)
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label('objective', fontsize=30)

    #ax1.tick_params(axis='both', which='major', labelsize=30)
    #ax2.tick_params(axis='both', which='major', labelsize=30)
    #ax3.tick_params(axis='both', which='major', labelsize=30)
    ax3.yaxis.offsetText.set_fontsize(20)

    fig1.savefig(os.path.join(plot_dir, 'surface_1850.png'), transparent=True,
                dpi=300)
    fig2.savefig(os.path.join(plot_dir, 'surface_2000.png'), transparent=True,
                 dpi=300)
    fig3.savefig(os.path.join(plot_dir, 'volume.png'), transparent=True,
                 dpi=300)
    # plt.savefig(os.path.join(plot_dir, 'surface_' +str(ys)+'_'+ gdir.rgi_id + '.png'), dpi=200)

    # plt.show()
    # plt.close()
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

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def plot_surface_col(gdir,df,ex_mod,ys, plot_dir):
    #df = df[df['objective']<=1000]


    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[-1].map_dx
    fig = plt.figure(figsize=cm2inch(40,24))
    grid = plt.GridSpec(2, 2, hspace=0.3, wspace=0.3)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1])
    ax3 = plt.subplot(grid[1, :])

    p2 = ax2.get_position()
    '''
    if gdir.name != '':
        plt.suptitle(gdir.rgi_id+':'+gdir.name,x=p2.x1/2,y=0.94,fontsize=20)
    else:
        plt.suptitle(gdir.rgi_id, x=p2.x1/2,y=0.94, fontsize=20)
    '''
    plt.suptitle(gdir.rgi_id+': Hintereisferner' , x=p2.x1 / 2, y=0.94,fontsize=25)
    import matplotlib as mpl
    import matplotlib.cm as cm

    norm = mpl.colors.LogNorm(vmin=0.1, vmax=1e5)
    norm = mpl.colors.LogNorm(vmin=df['objective'].min(),vmax=df['objective'].max())
    #cmap = matplotlib.cm.get_cmap('RdYlGn_r')
    cmap = matplotlib.cm.get_cmap('viridis')

    #cmap = matplotlib.cm.get_cmap('viridis')
    from colormap import Colormap
    c = Colormap()
    #cmap = c.cmap_linear('#0a3b70','navajowhite','#730421')

    df = df.sort_values('objective',ascending=False)
    for i,model in df['model'].iteritems():
        model = deepcopy(model)
        model.reset_y0(ys)
        color = cmap(norm(df.loc[i,'objective']))
        ax1.plot(x,deepcopy(model.fls[-1].surface_h),color=color, linewidth=2)
        model.volume_m3_ts().plot(ax=ax3,color=color, linewidth=2)

        model.run_until(2000)

        ax2.plot(x, model.fls[-1].surface_h, color=color, linewidth=2)

    '''
    # read diagnostics from ex_mod
    dp = gdir.get_filepath('model_diagnostics', filesuffix='experiment')
    ds = xr.open_dataset(dp)
    df = ds.to_dataframe()
    ax1.axhline(y=df.ela_m[df.index == ys].values)
    ax2.axhline(y=df.ela_m[df.index == 2000].values)
    '''
    # plot experiments
    ex_mod = deepcopy(ex_mod)
    ex_mod.volume_m3_ts().plot(ax=ax3, color='gold', linestyle=':', linewidth=3)
    ex_mod.reset_y0(1850)
    ex_mod.run_until(ys)

    ax1.plot(x, ex_mod.fls[-1].surface_h, ':',color='gold', linewidth=3)
    ax1.plot(x, ex_mod.fls[-1].bed_h, 'k', linewidth=3,label='bed topography')

    ex_mod.run_until(2000)

    ax2.plot(x, ex_mod.fls[-1].surface_h, ':',color='gold', linewidth=3)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', linewidth=3)

    ax1.annotate(r'$t =  '+str(ys)+'$', xy=(0.7,0.9), xycoords='axes fraction',
                 fontsize=20)
    ax2.annotate(r'$t =  2000$', xy=(0.7, 0.9), xycoords='axes fraction',
                 fontsize=20)

    ax1.set_ylabel('Altitude (m)', fontsize=20)
    ax1.set_xlabel('Distance along the main flowline (m)', fontsize=20)
    ax2.set_ylabel('Altitude (m)',fontsize=20)
    ax2.set_xlabel('Distance along the main flowline (m)',fontsize=20)
    ax3.set_ylabel(r'Volume ($m^3$)', fontsize=20)
    ax3.set_xlabel('Time (years)', fontsize=20)

    sm = plt.cm.ScalarMappable(cmap=cmap,norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax1,ax2,ax3])
    cbar = fig.colorbar(sm,cax=cax,**kw)
    cbar.ax.tick_params(labelsize=18)
    cbar.set_label('Fitness value', fontsize=20)

    ax1.tick_params(axis='both', which='major', labelsize=18)
    ax2.tick_params(axis='both', which='major', labelsize=18)
    ax3.tick_params(axis='both', which='major', labelsize=18)
    ax3.yaxis.offsetText.set_fontsize(18)

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

    l1=ax1.legend(handles=[lc3,lc,lc2], handler_map={lc: HandlerColorLineCollection(numpoints=100)},labels=['true solution', 'reconstruced ice surface', 'bed topography'],loc=3)

    l2=ax2.legend(handles=[lc3, lc, lc2],
               handler_map={lc: HandlerColorLineCollection(numpoints=100)},
               labels=['true solution', 'reconstruced ice surface', 'bed topography'],
               loc=3)
    ax3.legend(handles=[lc3],labels=['true solution'],loc=0)

    l1.set_zorder(0)
    l2.set_zorder(0)
    #plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'surface_1850_HEF_viridis.png'), transparent=True,dpi=300)
    #plt.savefig(os.path.join(plot_dir, 'surface_' +str(ys)+'_'+ gdir.rgi_id + '.png'), dpi=200)

    #plt.show()
    #plt.close()


def plot_surface_mean(gdir,df,ex_mod,t0,te,plot_dir):

    fig,ax1 = plt.subplots(figsize=(20, 15))
    ax2 = fig.add_axes([0.55, 0.66, 0.3, 0.2])

    ax1.set_title(gdir.rgi_id, fontsize=25)
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.95, box.height])
    ax1.annotate(r'$t = t_0 = ' + str(t0) + '$', xy=(0.1, 0.95),
                 xycoords='axes fraction',
                 fontsize=20)
    ax2.annotate(r'$t = ' + str(te) + '$', xy=(0.15, 0.85),
                 xycoords='axes fraction',
                 fontsize=15)

    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[
        -1].map_dx

    # add column for 2000 model
    df['model_te'] = df['model'].apply(lambda x: deepcopy(x))
    for model in df['model_te']:
        model.run_until(te)

    # 5% quantile and median of 5% quantile
    quant_df = df[df.objective <= df.objective.quantile(0.05)]

    # median of 5% quantile
    l = len(quant_df)
    if l%2:
        index = int((l-1)/2)
    else:
        index = int(l/2)
    model = deepcopy(quant_df.sort_values('volume').iloc[index]['model'])
    p1 = ax1.plot(x, deepcopy(model).fls[-1].surface_h, linewidth=2)
    model.run_until(te)
    ax2.plot(x, deepcopy(model).fls[-1].surface_h, linewidth=2)

    # 5%quantile at t0
    quant_list = quant_df['model'].apply(lambda x: deepcopy(x).fls[-1].surface_h).values
    quant_min = np.vstack(quant_list).max(axis=0)
    quant_max = np.vstack(quant_list).min(axis=0)
    ax1.fill_between(x, quant_min, quant_max, alpha=0.3)

    # same for te
    quant_list = quant_df['model_te'].apply(lambda x: deepcopy(x).fls[-1].surface_h).values
    quant_min = np.vstack(quant_list).max(axis=0)
    quant_max = np.vstack(quant_list).min(axis=0)
    ax2.fill_between(x, quant_min, quant_max, alpha=0.3)

    # all solutions with objective <100
    df100 = df[df['objective'] <= 100]
    if len(df100)>0:
        surface = df100['model'].apply(lambda x: deepcopy(x).fls[-1].surface_h).values
        max_t0 = np.vstack(surface).max(axis=0)
        min_t0 = np.vstack(surface).min(axis=0)
        ax1.fill_between(x, min_t0, max_t0, alpha=0.3,color='grey')

        surface = df100['model_te'].apply(lambda x: deepcopy(x).fls[-1].surface_h).values
        max_te = np.vstack(surface).max(axis=0)
        min_te = np.vstack(surface).min(axis=0)
        ax2.fill_between(x, min_te, max_te, alpha=0.3,color='grey')

    p2 = ax1.fill(np.NaN, np.NaN, alpha=0.3, color=p1[0].get_color(),
                  linewidth=0)
    p3 = ax1.fill(np.NaN, np.NaN, alpha=0.3, color='grey', linewidth=0)

    # print best objective value
    best = df.iloc[df.objective.idxmin()]
    ax1.plot(x,best['model'].fls[-1].surface_h,label='best')

    ex_mod =deepcopy(ex_mod)
    ex_mod.reset_y0(1850)
    ex_mod.run_until(1850)

    ax1.plot(x,deepcopy(ex_mod).fls[-1].surface_h,'k:',linewidth=2)
    ax1.plot(x,ex_mod.fls[-1].bed_h,'k', linewidth=2)
    ex_mod.run_until(2000)
    ax2.plot(x, deepcopy(ex_mod).fls[-1].surface_h, 'k:', linewidth=2)
    ax2.plot(x, ex_mod.fls[-1].bed_h, 'k', linewidth=2)

    try:
        ax1.legend([p3[0],(p2[0], p1[0])], [r'J(s) $\leq 100$' +'\n (N='+str(len(df100))+')',r'$\widetilde{J}(s)_{0.05}$'],loc='center left', bbox_to_anchor=(1, 0.5),fontsize=15)
    except:
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=15)
    ax1.set_xlabel('Distance along the Flowline',fontsize=25)
    ax1.set_ylabel('Altitude (m)',fontsize=25)

    ax2.set_xlabel('Distance along the Flowline (m)',fontsize=20)
    ax2.set_ylabel('Altitude (m)',fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=15)

    plot_dir = os.path.join(plot_dir,'best100')
    utils.mkdir(plot_dir,reset=False)
    #plt.savefig(os.path.join(plot_dir, 'best100_'+str(gdir.rgi_id)+'.png'))
    plt.show()
    plt.close()


def plot_candidates(gdir,df,yr,plot_dir):

    fig,ax = plt.subplots(figsize=(11,11))
    experiment = gdir.read_pickle('synthetic_experiment')
    prefix = 'model_run' + str(yr) + '_random'
    list = [f.split('model_run')[-1].split('.nc')[0] for f in
            os.listdir(gdir.dir) if
            f.startswith(prefix)]

    for f in list:
        try:
            rp = gdir.get_filepath('model_run', filesuffix=f)
            fmod = FileModel(rp)
            fmod.volume_m3_ts().plot(ax=ax,color='grey',label='',zorder=1)
        except:
            pass
    df['time'] = df['time'].apply(lambda x: float(x))
    t_eq = df['time'].sort_values().iloc[0]
    ax.axvline(x=t_eq,color='k',zorder=1, label=r'$t_{stag}$')

    df.plot.scatter(x='time',y='volume',ax=ax,c='black',s=250, zorder=2, label='candidates')
    df['Fitness value']=df['objective']
    from colormap import Colormap
    c = Colormap()
    cmap = c.cmap_linear('#0a3b70', 'navajowhite', '#730421')

    #df.plot.scatter(x='time',y='volume',ax=ax,c='Fitness value',
    #               colormap='viridis',norm=mpl.colors.LogNorm(vmin=0.1, vmax=1e5),s=250,
    #                zorder=2,edgecolors='k')
    label = r'temperature bias $\in [$' + str(
        df['temp_bias'].min()) + ',' + str(df['temp_bias'].max()) + '$]$'
    fmod.volume_m3_ts().plot(ax=ax, color='grey', label='', zorder=1)
    plt.ylabel(r'Volume $(m^3)$')
    plt.xlabel('Time (years)')
    plt.title(gdir.rgi_id)
    plt.legend()
    plt.savefig(os.path.join(plot_dir, 'method2.png'),
                dpi=300)
    #plt.close()
    #plt.show()