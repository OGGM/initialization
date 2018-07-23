from functools import partial
from pylab import *
from oggm.core.flowline import FluxBasedModel, FileModel
from oggm import graphics, tasks
from matplotlib import cm
import xarray as xr
import pandas as pd
from multiprocessing import Pool
from copy import deepcopy
FlowlineModel = partial(FluxBasedModel, inplace=False)
from scipy.optimize import curve_fit

import os




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
    cmap = matplotlib.cm.get_cmap('RdYlGn_r')
    df = df.sort_values('objective',ascending=False)
    for i,model in df['model'].iteritems():
        model.reset_y0(ys)
        color = cmap(norm(df.loc[i,'objective']))
        ax1.plot(x,deepcopy(model.fls[-1].surface_h),color=color, linewidth=2)
        model.volume_m3_ts().plot(ax=ax3,color=color, linewidth=2)

        model.run_until(2000)

        ax2.plot(x, model.fls[-1].surface_h, color=color, linewidth=2)

    # read diagnostics from ex_mod
    dp = gdir.get_filepath('model_diagnostics', filesuffix='experiment')
    ds = xr.open_dataset(dp)
    df = ds.to_dataframe()
    ax1.axhline(y=df.ela_m[df.index == ys].values)
    ax2.axhline(y=df.ela_m[df.index == 2000].values)

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

    #plt.savefig(os.path.join(plot_dir, 'surface_'+str(ys)+'_'+gdir.rgi_id+'.pdf'), dpi=200)
    #plt.savefig(os.path.join(plot_dir, 'surface_' +str(ys)+'_'+ gdir.rgi_id + '.png'), dpi=200)

    #plt.show()
    #plt.close()