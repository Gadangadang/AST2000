import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import matplotlib.path as mpath
import matplotlib.patches as mpatches

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem

utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

sigma = sigma


"""Loops through each planet and makes array with neccesary info"""
Planets = np.zeros((8,5))
for i in range(len(Planets)):
    Planet0_x0 = [system.initial_positions[0,i],system.initial_positions[1,i]]
    Planet0_v0 = [system.initial_velocities[0][i],system.initial_velocities[1][i]]
    Planets[i,0] = Planet0_x0[0]
    Planets[i,1] = Planet0_x0[1]
    Planets[i,2] = Planet0_v0[0]
    Planets[i,3] = Planet0_v0[1]
    Planets[i,4] = i

def temp_planet(i,avg_distanse):
    temp = (system.star_temperature**4/4*((system.star_radius*1000/AU)/(avg_distanse))**2)**(1/4)
    return temp

def habitable_zone(temp):
    r = np.sqrt((system.star_temperature**4)/(4*temp**4))*system.star_radius*1000/AU
    return r

for planet in Planets:
    temp = temp_planet(int(planet[4]),np.sqrt(planet[0]**2+planet[1]**2))

    print("Planet {} at distance {:g} has a temp of {:g}K, {:g}C".format\
    (planet[4],np.sqrt(planet[0]**2+planet[1]**2),temp,temp-273.15))

def run():
    Gr = G_sol
    sunmass = system.star_mass

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


    theta = np.linspace(0, 2*np.pi, 100)
    temperature = np.linspace(260,400,140)
    for i in temperature:
        r = habitable_zone(i)
        x1 = r*np.cos(theta)
        x2 = r*np.sin(theta)
        plt.plot(x1,x2,color="#98FB98")


    def plot_func(Planet):
        planet_orbit,planet_index = integrator(Planet)
        plt.plot(planet_orbit[:,0], planet_orbit[:,1],"--",label="{}".format(planet_index))
        plt.plot(planet_orbit[-1,0], planet_orbit[-1,1],"ko")

    for Planet in Planets:
        plot_func(Planet)



    plt.legend(loc="lower left")
    plt.axis("equal")
    plt.grid()
    plt.xlabel("AU")
    plt.ylabel("AU")
    plt.savefig("habitable_zone.jpeg")
    plt.show()
run()


#Solar panel calculation
A = (40/0.12)/(sigma*system.star_temperature**4*((system.star_radius*1000)\
/(np.sqrt(system.initial_positions[0,1]**2 + system.initial_positions[1,1]**2)*AU))**2)
print("Area of solar panel is {:g}m^2, with length and width {:g}m on planet 1".format(A, np.sqrt(A)), )
