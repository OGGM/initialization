
import matplotlib as mpl
import os
import numpy as np
import copy

mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['font.size'] =30
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.labelweight'] = 'medium'
mpl.rcParams['legend.fontsize'] = 30 #30
mpl.rcParams['lines.linewidth'] = 3



def animation(gdir, df, ex_mod, med_mod, plot_dir):

    import matplotlib.pyplot as plt
    from matplotlib import animation
    mpl.style.use('default')

    plot_dir = os.path.join(plot_dir, 'animation')
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig = plt.figure(figsize=(20, 15))
    ax1 = plt.axes()
    fill = ax1.fill_between([], [], color='grey', alpha=0.1, label='accepted')
    fill2 = ax1.fill_between([], [], color='C0', alpha=0.5, label='best')
    time_text = ax1.text(0.5, 0.95, '', transform=ax1.transAxes, size=30)


    plotlays, plotcols, label, linestyle  = [3], ["C0","C1","k","k"] , ['median state', 'minimum state','synth. experiment','bed rock'],['-','-',':','-']
    lines = []
    for index in range(4):
        lobj = ax1.plot([], [], lw=4, color=plotcols[index], label=label[index], linestyle=linestyle[index])[0]
        lines.append(lobj)

    # experiment model
    ex_mod = copy.deepcopy(ex_mod)

    # minimum model
    min_mod = copy.deepcopy(df.loc[df.fitness.idxmin(), 'model'])

    # median_model
    med_mod = copy.deepcopy(med_mod)

    # acceptables
    df = df[df.fitness < 125]
    acc_max = copy.deepcopy(df.loc[df.volume.idxmax(), 'model'])
    acc_min = copy.deepcopy(df.loc[df.volume.idxmin(), 'model'])

    # 5th percentile

    df = df[df.fitness < df.fitness.quantile(0.05)]
    quant_max = copy.deepcopy(df.loc[df.volume.idxmax(), 'model'])
    quant_min = copy.deepcopy(df.loc[df.volume.idxmin(), 'model'])

    x = np.arange(ex_mod.fls[-1].nx) * ex_mod.fls[-1].dx * ex_mod.fls[-1].map_dx


    def init():
        ax1.plot(x, ex_mod.fls[-1].bed_h, 'k', label='bed rock')
        time_text.set_text('')
        for line in lines:
            line.set_data([],[])
        return lines



    def animate(t):
        if t ==1850:

            ex_mod.reset_y0(1850)
            min_mod.reset_y0(1850)
            med_mod.reset_y0(1850)
            acc_max.reset_y0(1850)
            acc_min.reset_y0(1850)
            quant_max.reset_y0(1850)
            quant_min.reset_y0(1850)

        else:

            ex_mod.run_until(t)
            min_mod.run_until(t)
            med_mod.run_until(t)
            acc_max.run_until(t)
            acc_min.run_until(t)
            quant_max.run_until(t)
            quant_min.run_until(t)


        time_text.set_text('time = %.1f' % t)

        y1 = ex_mod.fls[-1].bed_h
        y2 = min_mod.fls[-1].surface_h
        y3 = med_mod.fls[-1].surface_h
        y4 = acc_max.fls[-1].surface_h
        y5 = acc_min.fls[-1].surface_h
        y6 = quant_max.fls[-1].surface_h
        y7 = quant_min.fls[-1].surface_h
        y8 = ex_mod.fls[-1].surface_h

        xlist = [x, x, x, x]
        ylist = [y3, y2, y8, y1]
        ax1.collections.clear()
        fill = ax1.fill_between(x, y4, y5, color='grey', alpha=0.2, label='accepted')
        fill2 = ax1.fill_between(x, y6, y7, color='C0', alpha=0.3, label = 'best')

        #for index in range(0,1):
        for lnum, line in enumerate(lines):
            line.set_data(xlist[lnum], ylist[lnum]) # set data for each line separately.

        return (fill2,)+tuple(lines)+(fill,) + (time_text,)

    # call the animator.  blit=True means only re-draw the parts that have changed.
    ani = animation.FuncAnimation(fig, animate, frames=range(1850, 2001, 1),
                        init_func=init, blit=True, repeat=False)

    plt.plot(x, ex_mod.fls[-1].bed_h, lw=4, color='k')
    plt.legend(loc='best', fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=25, width=3)
    plt.xlabel('Distance along the Flowline (m)', fontsize=30)
    plt.ylabel('Altitude (m)', fontsize=30)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(3)


    if gdir.name != "":
        plt.title(gdir.rgi_id)
    else:
        plt.title(gdir.rgi_id , fontsize=30)
    ani.save(os.path.join(plot_dir, 'surface_animation_'+gdir.rgi_id+'.mp4'), dpi=100)

