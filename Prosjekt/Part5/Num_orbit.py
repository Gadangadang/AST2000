import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem

utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)



Planets = np.zeros((8,5))
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
    def gravity(x,y,i):
        grav = -(Gr*system.star_mass)/(np.sqrt(x**2+y**2))**2
        ax = grav*x /(np.sqrt(x**2+y**2))
        ay = grav*y/(np.sqrt(x**2+y**2))
        return np.asarray([ax,ay])
    x0x,x0y,v0x,v0y,planet_index = Planet
    N = 1400000
    time = 7.6*20+3 #time for 20 rotations for our home planet
    dt = time/N
    t = dt
    x = np.zeros((N,2),float)
    x[0,:] = [x0x,x0y]
    v = np.asarray([v0x,v0y])
    a_i = gravity(x[0,0],x[0,1],planet_index)
    for i in range(N-1):
        t+=dt
        x[i+1,:] = x[i,:]+(v*dt + 0.5*a_i*dt**2)
        a_iplus1 = gravity(x[i+1,0],x[i+1,1],planet_index)
        v += 0.5*( a_i + a_iplus1 )*dt
        a_i = a_iplus1
    return x

Orbitsx = []
Orbitsy = []
for planet in Planets:
    Orbitsx.append(integrator(planet)[:,0])
    Orbitsy.append(integrator(planet)[:,1])
Orbits = np.asarray([Orbitsx,Orbitsy])
"""Orbits = np.asarray([np.dstack((Orbitsx)),np.dstack(Orbitsy)])"""
system.verify_planet_positions(7.6*20+3,Orbits, filename = 'Planet_orbits.npy')
