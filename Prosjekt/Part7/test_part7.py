import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts
#import part_4



#utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

planet_times, planet_positions = np.load("Planet_orbits_part7.npy",allow_pickle=True)
orbits = planet_positions

planet0 = (orbits[:][0][0],orbits[:][1][0])
planet1 = (orbits[:][1][0],orbits[:][1][1])


r2 = np.array([(planet1[0]-planet0[0])**2,(planet1[1]-planet0[1])**2])


r = np.sqrt(r2[0]+r2[1])

r_min = np.min(r)
a = int(np.where(r == r_min)[0])

N = 1400000
time = 7.6*40+3
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
    return np.array([fx(t),fy(t),0])



def plot_rute(tid,simulation_time):
    t = np.linspace(tid,tid+simulation_time,101)
    t2 = np.linspace(tid,tid+simulation_time,101)
    plt.plot(fx0(t),fy0(t),label="Orbit 0")
    plt.plot(fx1(t2),fy1(t2),label="Orbit 1")
    plt.plot(fx0(tid),fy0(tid),"go",label="Planet 0")
    plt.plot(fx1(tid),fy1(tid),"gx",label="Planet 1")
    plt.axis("equal")
    plt.xlabel("X in AU")
    plt.ylabel("Y in AU")
    plt.title("Planet Position T = {}".format(tid))
    plt.legend()
    plt.savefig("interpolate_position.jpeg")
    plt.show()
