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


import os


mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] = 30
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize']=11 #25
mpl.rcParams['lines.linewidth']=3



def plot_surface_col(gdir,df,ex_mod,ys, plot_dir):
    #df = df[df['objective']<=1000]

    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[-1].map_dx
    fig = plt.figure(figsize=(20,15))
    grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1])
    ax3 = plt.subplot(grid[1, :])

    p2 = ax2.get_position()

    if gdir.name != '':
        plt.suptitle(gdir.rgi_id+':'+gdir.name,x=p2.x1/2,fontsize=20)
    else:
        plt.suptitle(gdir.rgi_id, x=p2.x1/2, fontsize=20)

    import matplotlib as mpl
    import matplotlib.cm as cm

    norm = mpl.colors.LogNorm(vmin=0.1, vmax=1e5)
    #norm = mpl.colors.LogNorm(vmin=df['objective'].min(),vmax=df['objective'].max())
    #cmap = matplotlib.cm.get_cmap('RdYlGn_r')
    cmap = matplotlib.cm.get_cmap('viridis')
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
    #plot real flowlines
    fls = gdir.read_pickle('model_flowlines')
    ax2.plot(x,fls[-1].surface_h,'k')
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

    ax1.annotate(r'$t =  '+str(ys)+'$', xy=(0.8,0.9), xycoords='axes fraction',
                 fontsize=15)
    ax2.annotate(r'$t =  2000$', xy=(0.8, 0.9), xycoords='axes fraction',
                 fontsize=15)

    ax1.set_ylabel('Altitude (m)', fontsize=15)
    ax1.set_xlabel('Distance along the main flowline (m)', fontsize=15)
    ax2.set_ylabel('Altitude (m)',fontsize=15)
    ax2.set_xlabel('Distance along the main flowline (m)',fontsize=15)
    ax3.set_ylabel(r'Volume ($m^3$)', fontsize=15)
    ax3.set_xlabel('Time (years)', fontsize=15)

    sm = plt.cm.ScalarMappable(cmap=cmap,norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax1,ax2,ax3])
    cbar = fig.colorbar(sm,cax=cax,**kw)
    cbar.ax.tick_params(labelsize=15)
    cbar.set_label('objective', fontsize=15)

    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax3.tick_params(axis='both', which='major', labelsize=15)
    ax3.yaxis.offsetText.set_fontsize(15)

    #plt.savefig(os.path.join(plot_dir,'surface_by_fitness','surface_'+str(ys)+'_'+gdir.rgi_id+'.pdf'), dpi=200)
    plt.savefig(os.path.join(plot_dir, 'surface_' +str(ys)+'_'+ gdir.rgi_id + '.png'), dpi=200)

    #plt.show()
    plt.close()


def plot_surface_mean(gdir,df,ex_mod,t0,te,plot_dir):

    fig,ax1 = plt.subplots(figsize=(20, 15))
    ax2 = fig.add_axes([0.55, 0.66, 0.3, 0.2])
    if gdir.name != "":
        ax1.set_title(gdir.rgi_id + ': ' + gdir.name, fontsize=30)
    else:
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

    ex_mod =deepcopy(ex_mod)
    ex_mod.reset_y0(1850)
    ex_mod.run_until(t0)

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

    plt.savefig(os.path.join(plot_dir, 'median_' +str(t0)+'_'+ gdir.rgi_id + '.png'))
    #plt.show()
    plt.close()


def plot_climate(gdir, plot_dir) :
    #plot_dir = os.path.join(plot_dir, 'surface')
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    d = xr.open_dataset(gdir.get_filepath('climate_monthly'))
    print(d)
    temp = d.temp.resample(freq='12MS', dim='time', how=np.mean).to_series()
    temp = temp[temp.index.year >= 1850]

    fig, ax1 = plt.subplots(figsize=(15, 10))
    del temp.index.name
    temp.plot( linewidth=3, label='Annual temp')
    temp.rolling(31, center=True, min_periods=15).mean().plot(linewidth=3,
        label='31-yr avg')
    ax1.legend(loc='best', fontsize=25)
    plt.title('HISTALP annual temperature:'+gdir.name, fontsize=25)
    plt.ylabel(r'degC', fontsize=25)
    plt.xlabel(r'Time', fontsize=25)

    ax1.tick_params(axis='both', which='major', labelsize=20)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir,'climate', gdir.rgi_id+'_Histalps.png'), dpi=300)
    #plt.show()

def plot_fitness_over_time(gdir,df_list,plot_dir):
    plt.figure()
    norm = mpl.colors.LogNorm(vmin=0.1, vmax=1e5)
    cmap = matplotlib.cm.get_cmap('viridis')



    for year in df_list.keys():
        df_list[year] = df_list[year].sort_values('objective',ascending=False)
        volume = df_list[year].volume
        yr = np.zeros(len(volume))+int(year)
        plt.scatter(x=yr,y=volume,c=df_list[year].objective,cmap=cmap,norm=norm,marker='s')
        plt.savefig(os.path.join(plot_dir, 'starting' '_' + gdir.rgi_id + '.png'))
    plt.show()
    plt.close()

def plot_fitness_over_time2(gdir,df_list,ex_mod,plot_dir):
    from matplotlib.patches import Rectangle,Circle
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111,yscale='linear')
    norm = mpl.colors.LogNorm(vmin=0.1, vmax=1e5)
    cmap = matplotlib.cm.get_cmap('viridis')

    volumes = np.linspace(df_list['1850'].volume.min(),df_list['1850'].volume.max(),100)

    for year,df in df_list.items():
        color=[]
        for i in range(len(volumes)):
            if i!= len(volumes)-1:
                part = df[(df.volume >= volumes[i]) & (df.volume<=volumes[i+1])]
                color.append(part.objective.min())
            else:
                part= df[df.volume >= volumes[i]]
                color.append(part.objective.mean()) # or min

        # interpolate missing data
        missing = np.where(np.isnan(color))
        if len(missing[0]) != 0:
            xp = np.delete(range(len(volumes)),missing)
            fp = np.delete(color,missing)
            missing_y = np.interp(missing, xp, fp)
            for i,j in enumerate(missing[0]):
                color[j] = missing_y[0][i]

        year = np.zeros(len(volumes))+int(year)
        for x,y,c in zip(volumes,year,color):
            if np.isnan(c):
                color='white'
            else:
                color=cmap(norm(c))
            ax.add_patch(Rectangle((y-2.5, x), 5, volumes[1]-volumes[0], color=color))

    # add experiment in plot
    ex_mod.volume_m3_ts().plot(ax=ax,linestyle=':',color='gold')
    ax.set_xlim(1847.5,1967.5)
    ax.set_ylim(volumes[0],volumes[-1])
    plt.title(gdir.rgi_id + ': '+ gdir.name, fontsize=25)
    plt.ylabel(r'Volume ($m^3$)', fontsize=25)
    plt.xlabel(r'Starting time', fontsize=25)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax, kw = mpl.colorbar.make_axes([ax])
    cbar = fig.colorbar(sm, cax=cax, **kw)
    cbar.ax.tick_params(labelsize=25)
    cbar.set_label('Fitness value', fontsize=25)
    #plt.savefig(os.path.join(plot_dir, 'starting' '_' + gdir.rgi_id + '.png'))

    #plt.close()


