import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts

from ast2000tools.relativity import RelativityExperiments

utils.check_for_newer_version()
seed = utils.get_seed('Sgfrette')
experiments = RelativityExperiments(seed)
system = SolarSystem(seed)

planet_idx = 1 # I want to perform the experiment near planet 4
#experiments.black_hole_descent(planet_idx, number_of_light_signals=30, consider_light_travel=True, text_file_dir='.', filename_1='black_hole_descent_frame_1_True.xml', filename_2='black_hole_descent_frame_2_True.xml')


x,y = np.loadtxt("black_hole_descent_frame_1.txt")
x1,y1 = np.loadtxt("black_hole_descent_frame_1_True_with_light_travel.txt")

x_without = []
y_without = []

x_with = []
y_with = []

for i in range(1,len(x)):

    x_without.append(x[i]-x[i-1]+i)
    y_without.append(y[i]-y[i-1])
    x_with.append(x1[i]-x1[i-1]+i)
    y_with.append(y1[i]-y1[i-1])



plt.plot(x_without,y_without,label="Without light travel time")
plt.plot(x_with,y_with,label="With light travel time")
plt.legend()
plt.xlabel("Signal count")
plt.ylabel("Time [s]")
plt.title("Shell Frame")
plt.savefig("Frame1.jpeg")
plt.show()

x,y = np.loadtxt("black_hole_descent_frame_2.txt")
x1,y1 = np.loadtxt("black_hole_descent_frame_2_True_with_light_travel.txt")

x_without = []
y_without = []

x_with = []
y_with = []

for i in range(1,len(x)):

    x_without.append(x[i]-x[i-1]+i)
    y_without.append(y[i]-y[i-1])
    x_with.append(x1[i]-x1[i-1]+i)
    y_with.append(y1[i]-y1[i-1])



plt.plot(x_without,y_without,label="Without light travel time")
plt.plot(x_with,y_with,label="With light travel time")
plt.legend()
plt.xlabel("Signal count")
plt.ylabel("Time [s]")
plt.title("Falling Frame")
plt.savefig("Frame2.jpeg")
plt.show()
