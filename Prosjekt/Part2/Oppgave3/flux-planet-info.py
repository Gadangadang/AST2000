import numpy as np
import matplotlib.pyplot as plt
from numba import jit

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
"""
This is the code for c3
"""
#utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)


#Reading code, separates values in order (time,flux) into to lists
infile = open("light_curve_data_carl.txt", "r")
infile.readline()

time = []
flux = []
for line in infile:
    word = line.split()
    time.append(float(word[0]))
    flux.append(float(word[1]))


#Setting values based on observations on plot
t0 = 0.0005#time right before planet enters line of sight
t1 = 0.0006#time right after whole planet is within line of sight
t2 = 0.0084# right before planet leaves line of sight

#Calculate values for mass, velocity, radii and density of planet,
#also radii of star
v_star = 5*10**(-4) #AU/yr
m_planet = 0.000668*m_sun #kg     m_star**(1/3)*v_star_radial*P**(1/3)/(2*np.pi*G)**(1/3)
v_planet_cm = 2.17 #AU/yr  v_star*(m_star/m_planet)

Radius_planet = (v_star+v_planet_cm)*(t1-t0)/2*AU*1000 # meters
density_planet = m_planet/((4/3)*np.pi*Radius_planet**3) #kg/m^3
Radius_star = (v_star+v_planet_cm)*(t2-t0)/2*AU #meters

print("Planet radius = {}".format(Radius_planet))
print("Planet density ={}".format(density_planet))
print("Star radius = {}".format(Radius_star))
