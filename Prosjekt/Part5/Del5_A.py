import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts
import test

utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)
codes = [8208,389] # The codes recieved from the group teacher (in this case just one)
shortcuts = SpaceMissionShortcuts(mission, codes)



shortcuts.place_spacecraft_on_escape_trajectory(17.5,160000,pi/3, 10000)
fuel_mass_after_lanuch, time_after_launch, rocket_position_after_launch,\
 rocket_velocity_after_launch = shortcuts.get_launch_results()

print("Fuel: {} Time: {} Pos: {} Vel: {}".format(fuel_mass_after_lanuch, time_after_launch, rocket_position_after_launch,\
 rocket_velocity_after_launch))

print(np.linalg.norm(rocket_position_after_launch))

fuel_mass_after_lanuch= 10000

"""planet_times, planet_positions = np.load("Planet_orbits.npy",allow_pickle=True)
Orbits = planet_positions"""


Gr = G_sol
m_ship = (mission.spacecraft_mass + fuel_mass_after_lanuch)/m_sun
# Awesome code
def grav_sun(r):
    return -Gr*m_ship*system.star_mass*r/np.linalg.norm(r)**3
def grav_planet(r,f_planet,t):
    a = np.zeros(2)

    for i in range(8):
        a -= Gr*m_ship*system.masses[i]*\
        (r-np.sqrt(test.f_planet(i,t)[0]**2+test.f_planet(i,t)[1]**2))\
        /np.linalg.norm((r-test.f_planet(i,t)))**3
    return a

def ship_traj(initial_time,initial_pos,initial_vel,simulation_time,dt_length):
    dt = dt_length
    N = int((initial_time+simulation_time)/dt_length)




    x = np.zeros((N,2),float)

    x[0,:] = initial_pos
    v = initial_vel
    t = initial_time
    a_i = np.asarray(grav_sun(x[0,:])+ grav_planet(x[0,:],test.f_planet,initial_time))

    for i in range(int(initial_time/dt),N-1):
        t+=dt
        x[i+1,:] = x[i,:]+(v*dt + 0.5*a_i*dt**2)
        print(x[i+1,:])
        a_iplus1 = np.asarray(grav_sun(x[i,:])+ grav_planet(x[i,:],test.f_planet,t))
        v = v+ 0.5*( a_i + a_iplus1 )*dt
        a_i = a_iplus1
    return x, v, t




if __name__ == "__main__":
    x,v,t = ship_traj(time_after_launch,rocket_position_after_launch,\
    rocket_velocity_after_launch,1.5,(7.6*40+3)/1400000)
    print("Final position: {} Final Velocity: {} Final Time: {}".\
    format(x[-1,0],v,t))
"""
    plt.plot(x[0,:],x[1,:])
    plt.axis("equal")
    plt.show()"""
