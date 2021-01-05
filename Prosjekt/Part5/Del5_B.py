import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts
import test

#utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)
codes = [8208,389] # The codes recieved from the group teacher (in this case just one)
shortcuts = SpaceMissionShortcuts(mission, codes)

fuel_mass_after_lanuch= 10000
Gr = G_sol
m_ship = (mission.spacecraft_mass + fuel_mass_after_lanuch)/m_sun

def grav_sun(r):
    return -Gr*m_ship*system.star_mass*r/np.linalg.norm(r)**3
def grav_planet(r,t):
    a = np.zeros(2)
    for i in range(8):
        a -= Gr*m_ship*system.masses[i]*\
        (r-np.asarray([test.f_planet(i,t)[0],test.f_planet(i,t)[1]]))\
        /np.linalg.norm((r-test.f_planet(i,t)))**3
    return a

def ship_traj(initial_time,initial_pos,initial_vel,simulation_time,dt_length,mass_ship):
    dt = dt_length
    N = int((initial_time+simulation_time)/dt)
    x = np.zeros((len(range(int(initial_time/dt),N))+1,2),float)
    x[0,:] = initial_pos
    v = initial_vel
    t = initial_time
    a_i = np.asarray(grav_sun(x[0,:])+ grav_planet(x[0,:],t))/m_ship
    j = 0
    for i in range(int(initial_time/dt),N):
        t+=dt
        x[j+1,:] = x[j,:]+(v*dt + 0.5*a_i*dt**2)
        j +=1
        a_iplus1 = np.asarray(grav_sun(x[j,:])+ grav_planet(x[j,:],t))/m_ship
        v = v + 0.5*( a_i + a_iplus1 )*dt
        a_i = a_iplus1
    return x, v, t
