import numpy as np
import matplotlib.pyplot as plt
from numba import jit

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem

#utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed(6421)
mission = SpaceMission(seed)
system = SolarSystem(seed)
"""
This is the code for task A1 and A2.
"""

"""a,e lists with eccentricities and semi_major_axes for each planet"""
e = [i for i in system.eccentricities]
a = [i for i in system.semi_major_axes]

name = [i for i in range(8)]
x= zip(e,a,name)

def orbit_analytical(e,a,f):
    r = (a*(1-e**2))/(1+e*np.cos(f))
    return r

angle = np.linspace(0,2*np.pi,1000)

"""Loop that calculates the orbit vector, and generates x,y coordinates for each point"""
for i in x:
    r = orbit_analytical(i[0],i[1],angle)
    x = r*np.cos(angle)
    y = r*np.sin(angle)
    plt.plot(x, y, "--",label="{}".format(i[2]))

"""Loops through each planet and makes array with neccesary info"""
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
sunmass = system.star_mass

"""Numerical integrator that calculates orbit with leapfrog method"""
def integrator(Planet):
    def gravity(x,y):
        grav = -(Gr*sunmass)/(np.sqrt(x**2+y**2))**2
        ax = grav*x /(np.sqrt(x**2+y**2))
        ay = grav*y/(np.sqrt(x**2+y**2))
        return np.asarray([ax,ay])
    x0x,x0y,v0x,v0y,planet_index = Planet

    N = 1400000
    time = 7.6*40 #time for 20 rotations for our home planet
    dt = time/N
    t = dt
    x = np.zeros((N,2),float)
    x[0,:] = [x0x,x0y]
    v = np.asarray([v0x,v0y])
    a_i = gravity(x[0,0],x[0,1])
    for i in range(N-1):
        t+=dt
        x[i+1,:] = x[i,:]+(v*dt + 0.5*a_i*dt**2)
        a_iplus1 = gravity(x[i+1,0],x[i+1,1])
        v += 0.5*( a_i + a_iplus1 )*dt
        a_i = a_iplus1
    return x, planet_index

verifu = np.zeros(2, 3, N=1400000)
#Plotting
def plot_func(Planet):
    planet_orbit,planet_index = integrator(Planet)
    plt.plot(planet_orbit[:,0], planet_orbit[:,1],label="{}".format(planet_index))

for Planet in Planets:
    plot_func(Planet)


plt.legend()
plt.xlabel("x-position (AU)")
plt.ylabel("y-position (AU)")
plt.savefig("Numeric-analy-moresteps.jpeg")

plt.show()
