import numpy as np
from .particle import Particle

def track(particle, line):
    for element in line:
        element.track(particle)
        particle.update_history()
        if particle.lost != 'CIRCULATING':
            return
    return

class TrackGrid:
    """Array of particles tracked through line starting from gridpoints"""
    def __init__(self, line, xmin, xmax, xres, xpmin, xpmax, xpres):
        self.xmin = xmin
        self.xmax = xmax
        self.xres = xres
        self.xpmin = xpmin
        self.xpmax = xpmax
        self.xpres = xpres

        self.nx = round((xmax-xmin)/xres)
        self.npx = round((xpmax-xpmin)/xpres)

        self.particles = np.zeros((self.nx, self.npx), dtype=object)

        for index, _ in np.ndenumerate(self.particles):
            # due to matrix indexing we start at xmin,xpmax
            self.particles[index] = Particle(xmin+index[0]*xres,
                                             xpmax-index[1]*xpres)
            track(self.particles[index], line)

class TrackList:
    """List of particles tracked through line starting from initial conditions"""
    def __init__(self, line, inits):
        self.particles = np.zeros((len(inits), 1), dtype=object)

        for index, _ in np.ndenumerate(self.particles):
            self.particles[index] = Particle(inits[index[0]][0],
                                             inits[index[0]][1])
            track(self.particles[index], line)
            self.particles[index].lost = str(index[0])
