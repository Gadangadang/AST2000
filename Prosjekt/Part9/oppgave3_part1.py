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
#experiments.black_hole_descent(planet_idx, number_of_light_signals=30, consider_light_travel=False, text_file_dir='.', filename_1='black_hole_descent_frame_1.xml', filename_2='black_hole_descent_frame_2.xml')

# Finn E/m
m_kg = 2.2E7*m_sun
Mm = G/(c**2)*m_kg
EoverM = np.sqrt(1-2*Mm/(1*AU))*1/(np.sqrt(1-0.26**2))
print(EoverM)

frame1 = []
frame2 = []
#Finn avstand fra rakett til sort hull:

infile = open("black_hole_descent_frame_1.txt", "r")
for line in infile: # for loop to split each row into seperate strings
    num = line.split() #creating a nested list wrd which contains all 5 lists from .dat file
    for i in num:
        frame1.append(float(i))

infile = open("black_hole_descent_frame_2.txt", "r")
for line in infile: # for loop to split each row into seperate strings
    num = line.split() #creating a nested list wrd which contains all 5 lists from .dat file
    for i in num:
        frame2.append(float(i))

print(frame1)
dtstart= frame1[1]-frame1[0]
dtaustart= frame2[1]-frame2[0]

dtslutt= frame1[-1]-frame1[-2]
dtauslutt= frame2[-1]-frame2[-2]

def r(dt,dtau):
    return 2*Mm/(1-(EoverM*dtau/dt)**2)
Schwarz = 2*m_kg*G/c**2

print("Schwartchild radii is {:g}km".format(Schwarz/1000))

print("r_start = {}AU, r_start = {} Schwarz".format(r(dtstart,dtaustart)/AU,r(dtstart,dtaustart)/Schwarz))
print("r_slutt = {}AU, r_slutt = {} Schwarz".format(r(dtslutt,dtauslutt)/AU,r(dtslutt,dtauslutt)/Schwarz))
