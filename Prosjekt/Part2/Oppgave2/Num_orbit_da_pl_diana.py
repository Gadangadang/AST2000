import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem

utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed(6421)
mission = SpaceMission(seed)
system = SolarSystem(seed)

"""
This is the code for task B2
"""

"""
All the code under is the same as for num_orbit
"""
Planets = np.zeros((7,5))
for i in range(len(Planets)):
    Planet0_x0 = [system.initial_positions[0,i],system.initial_positions[1,i]]
    Planet0_v0 = [system.initial_velocities[0][i],system.initial_velocities[1][i]]
    Planets[i,0] = Planet0_x0[0]
    Planets[i,1] = Planet0_x0[1]
    Planets[i,2] = Planet0_v0[0]
    Planets[i,3] = Planet0_v0[1]
    Planets[i,4] = i

Gr = G_sol
def integrator(Planet):
    def gravity(x,y):
        grav = -(Gr*system.star_mass)/(np.sqrt(x**2+y**2))**2
        ax = grav*x /(np.sqrt(x**2+y**2))
        ay = grav*y/(np.sqrt(x**2+y**2))
        return np.asarray([ax,ay])

    x0x,x0y,v0x,v0y,planet_index = Planet

    N = 1400000
    time = 7.6*40 #time for 20 rotations for our home planet
    #print(time)
    dt = time/N
    t = dt
    x = np.zeros((N,2),float)
    x[0,:] = [x0x,x0y]
    v = np.asarray([v0x,v0y])
    a_i = gravity(x[0,0],x[0,1])
    A_1 = 0
    A_2 = 0
    for i in range(N-1):
        t+=dt
        x[i+1,:] = x[i,:]+(v*dt + 0.5*a_i*dt**2)
        a_iplus1 = gravity(x[i+1,0],x[i+1,1])
        v += 0.5*( a_i + a_iplus1 )*dt
        a_i = a_iplus1
        """
        The parameters for the if test were randomly choosen
        """
        if t>=2 and t<=3:
            A_1 += np.linalg.norm(v)*np.linalg.norm(x[i+1])*0.5*dt
        if t>=5 and t<=6:
            A_2 += np.linalg.norm(v)*np.linalg.norm(x[i+1])*0.5*dt
    return x, planet_index,A_1,A_2



def plot_func(Planet):
    planet_orbit,planet_index,A_1,A_2 = integrator(Planet)
    plt.plot(planet_orbit[:,0], planet_orbit[:,1],label="{}".format(planet_index))
    print("Areal 1 for planet {} = {}".format(planet_index,A_1))
    print("Areal 2 for planet {} = {}".format(planet_index,A_2))
for Planet in Planets:
    plot_func(Planet)


plt.show()

#system.generate_system_snapshot()
