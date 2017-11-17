import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def _color_losses(particles, colorcodes):
    colored_losses = np.zeros_like(particles, dtype=int)

    for c_index, colorcode in enumerate(colorcodes):
        if colorcode[0] != 'Other':
            for p_index, particle in np.ndenumerate(particles):
                contains_checks = all(sub in particle.lost
                                      for sub in colorcode[1])
                exclude_checks = all(sub not in particle.lost
                                     for sub in colorcode[2])
                options_empty = (len(colorcode[3])==0)
                options_checks = any(sub in particle.lost
                                      for sub in colorcode[3])
                if (contains_checks and exclude_checks
                    and (options_empty or options_checks)):
                    colored_losses[p_index] = c_index+1

    if colorcodes[-1][0] == 'Other':
        for index, entry in np.ndenumerate(colored_losses):
            if entry == 0:
                colored_losses[index] = len(colorcodes)

    return colored_losses

def acceptanceplot(trackgrid, colorcodes, colormap, show=True,
                   filename=None, beam=None, aperture=None):
    colored_losses = _color_losses(trackgrid.particles, colorcodes)

    xscale = 1E3
    xlabel = "$x$  [mm]"
    xpscale = 1E3
    xplabel = "$x'$  [mrad]"

    xmin = trackgrid.xmin*xscale
    xmax = trackgrid.xmax*xscale
    xpmin = trackgrid.xpmin*xpscale
    xpmax = trackgrid.xpmax*xpscale
    cols = len(colorcodes)
    colorlabels = [colorcode[0] for colorcode in colorcodes]

    fig, ax = plt.subplots()
    cax = ax.imshow(colored_losses.transpose(), interpolation='none',
                    extent=[xmin, xmax, xpmin, xpmax],
                    aspect=(xmax-xmin)/(xpmax-xpmin), cmap=colormap,
                    vmin=0.5, vmax=cols+0.5)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(xplabel)

    if beam is not None:
        beam = [[point[0]*xscale, point[1]*xpscale] for point in beam]
        ax.add_patch(Polygon(beam, closed=True, fill=True, edgecolor='none',
                     facecolor=(0.5, 0.5, 0.5, 0.5)))
    cbar = fig.colorbar(cax, ticks=[i+1 for i in range(cols)])
    cbar.ax.set_yticklabels(colorlabels)
    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename)
    if show is True:
        plt.show()
    return

def trajectoryplot(tracks, colorcodes, colormap, aperture=None,
                   show=True, filename=None, reducepoints=[1,1], linewidth=0.1):
    colored_losses = _color_losses(tracks.particles, colorcodes)

    xlabel = "$s$  [m]"
    xplabel = "$x$  [m]"

    cols = len(colorcodes)
    colorlabels = [colorcode[0] for colorcode in colorcodes]

    fig, ax = plt.subplots()

    for index, particle in np.ndenumerate(tracks.particles):
            if index[0]%reducepoints[0] == 0 and index[1]%reducepoints[1] == 0:
                data_s = [coordinate[0] for coordinate in particle.history]
                data_x = [coordinate[1] for coordinate in particle.history]
                cax = ax.plot(data_s, data_x, linewidth=linewidth,
                              color=colormap(colored_losses[index]-1))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(xplabel)

    if aperture is not None:
        for patch in aperture:
            ax.add_patch(Polygon(patch, closed=True, fill=True, edgecolor='none',
                         facecolor=(0, 0, 0, 0.65)))
    #cbar = fig.colorbar(cax, ticks=[i+1 for i in range(cols)])
    #cbar.ax.set_yticklabels(colorlabels)
    if filename is not None:
        plt.savefig(filename)
    if show is True:
        plt.show()
    return
