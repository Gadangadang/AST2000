import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts
#import part_4
#import Traj_test



utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

planet_times, planet_positions = np.load("Planet_orbits.npy",allow_pickle=True)
orbits = planet_positions


planet0 = (orbits[:][0][0],orbits[:][1][0])
planet1 = (orbits[:][1][0],orbits[:][1][1])


r2 = np.array([(planet1[0]-planet0[0])**2,(planet1[1]-planet0[1])**2])


r = np.sqrt(r2[0]+r2[1])

r_min = np.min(r)
a = int(np.where(r == r_min)[0])



N = 1400000
time = 7.6*40
dt = time/N

tid = a*dt


def interpol(a,b):

    f = interpolate.interp1d(a,b,axis=0)
    return f

fx0 = interpol(planet_times,orbits[:][0][0])
fy0 = interpol(planet_times,orbits[:][1][0])
fx1 = interpol(planet_times,orbits[:][0][1])
fy1 = interpol(planet_times,orbits[:][1][1])
fx2 = interpol(planet_times,orbits[:][0][2])
fy2 = interpol(planet_times,orbits[:][1][2])
fx3 = interpol(planet_times,orbits[:][0][3])
fy3 = interpol(planet_times,orbits[:][1][3])
fx4 = interpol(planet_times,orbits[:][0][4])
fy4 = interpol(planet_times,orbits[:][1][4])
fx5 = interpol(planet_times,orbits[:][0][5])
fy5 = interpol(planet_times,orbits[:][1][5])
fx6 = interpol(planet_times,orbits[:][0][6])
fy6 = interpol(planet_times,orbits[:][1][6])
fx7 = interpol(planet_times,orbits[:][0][7])
fy7 = interpol(planet_times,orbits[:][1][7])

f_plan = [(fx0,fy0),(fx1,fy1),(fx2,fy2),(fx3,fy3),(fx4,fy4),(fx5,fy5),(fx6,fy6),(fx7,fy7)]


def f_planet(i,t):
    fx = f_plan[i][0]
    fy = f_plan[i][1]
    return tuple([fx(t),fy(t)])

"""
xrocket = Traj_test.x[:,0]
yrocket = Traj_test.x[:,1]
time = np.zeros(len(xrocket))
dt = (7.6*40+3)/1400000


time[0] = 0
for i in range(len(time)-1):
    time[i+1] = time[i] +  dt


xrock = interpol(time,xrocket)
yrock = interpol(time,yrocket)
"""
fig = plt.figure()
ax = plt.axes(xlim=(-9, 9), ylim=(-9, 9))
ax.set_xlabel("AU")
ax.set_ylabel("AU")

line, = ax.plot([], [],"-", lw=2)

plotlays, plotcols = [2], ["black","red"]
lines = []
for index in range(2):
    lobj = ax.plot([],[],lw=2,label="Planet {}".format(index),color=plotcols[index])[0]
    lines.append(lobj)


x1,y1 = [],[]
x2,y2 = [],[]
x3,y3 = [],[]

# initialization function
def init():
    for line in lines:
        line.set_data([],[])
    return lines


# animation function
def animate(i):
    t = 0.01*i
    tid = 15.75
    x = f_planet(0,tid + t)[0]
    y = f_planet(0,tid + t)[1]
    x1.append(x)
    y1.append(y)

    x = f_planet(1,tid + t)[0]
    y = f_planet(1,tid + t)[1]
    x2.append(x)
    y2.append(y)
    """
    x =f_planet(1,tid + t)[0]-f_planet(0,tid + t)[0]
    y =f_planet(1,tid + t)[1]-f_planet(0,tid + t)[1]
    x3.append(x)
    y3.append(y)"""

    xlist = [x1, x2]
    ylist = [y1, y2]

    #for index in range(0,1):
    for lnum,line in enumerate(lines):
        if lnum is not 2:
            line.set_data(xlist[lnum], ylist[lnum]) # set data for each line separately.
            line.set_label('Planet ' + str(lnum))
        else:
            line.set_data(xlist[lnum], ylist[lnum]) # set data for each line separately.
            line.set_label('Distance')

        legend = plt.legend()
    ax.set_title("Time (yr) {}".format(tid + t))
    return lines



ani = animation.FuncAnimation(fig, animate,init_func=init,
     frames=1000, interval=10, blit= True)

ani.save('orbit_with_rocket.mp4', writer="ffmpeg")
fig.show()
