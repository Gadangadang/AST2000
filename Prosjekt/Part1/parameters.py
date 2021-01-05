
import numpy as np
from numba import jit
import matplotlib.pyplot as plt
from ast2000tools.constants import *
import seaborn as sns

np.random.seed(123)

#CONSTANTS
k = k_B          #Boltzmann's Constant

#PHYSICAL VARIABLES
T = 3E3               #Temperature in Kelvin
L = 10E-7        #Box Length in meters
N = 10E5                 #Number of particles
m = m_H2          #Mass of individual particle in kg
short_side = L/2 - np.sqrt((1/16)*L**2)/2
long_side = L/2 + np.sqrt((1/16)*L**2)/2

#INITIAL CONDITIONS
sigma = np.sqrt((k*T)/m)     #The standard deviation of our particle velocities

x =  np.random.uniform(0,L, size = (int(N), 3))  # Position vector, fill in uniform distribution in all 3 dimensions
v =  np.random.normal(0,sigma, size = (int(N), 3)) # Velocity vector, fill in correct distribution

#SIMULATION VARIABLES
time = 10E-9                       #Simulation Runtime in Seconds
steps = 1000                 #Number of Steps Taken in Simulation
dt = time/steps                 #Simulation Step Length in Seconds

#RUNNING THE CONDITIONAL INTEGRATION LOOP
#@jit(nopython=True)
def all_part(x,v):
    exiting = 0         #The total number of particles that have exited the gas box
    f = 0               #Used to calculate Force/Second later in the simulation
    col_wall = 0        #Collisions with wall
    p_zdirection = 0 #momentum in z direction each time step
    refill = np.array([0, 0, L])
    condition_exit = np.zeros_like(x, dtype = bool)

    for i in range(int(steps)):
        x+= dt*v # Change in position
        condition1 = np.logical_and(np.greater_equal(x,L), np.greater(v,0))
        condition2 = np.logical_and(np.less_equal(x,0), np.less(v,0))
        condition_collision = np.logical_or(condition1,condition2)

        condition3 = np.logical_and(np.less_equal(x[:,2], 0), np.less(v[:,2], 0))
        condition4 = np.logical_and(np.greater_equal(x[:,0], short_side), np.less_equal(x[:,0],long_side))
        condition5 = np.logical_and(np.greater_equal(x[:,1], short_side), np.less_equal(x[:,1],long_side))

        temp = np.logical_and(condition3, condition4, condition5)
        condition_exit[:,0] = temp
        condition_exit[:,1] = temp
        condition_exit[:,2] = temp
        exiting += np.sum(temp)
        p_zdirection += np.sum(temp*abs(v[:,2]))

        refill_mat = condition_exit.astype(np.int8)
        flip_bool = np.logical_and(condition_collision, np.logical_not(condition_exit)).astype(np.int8)

        x[:,2] = x[:,2] + refill_mat[:,0]*L

        col_wall = col_wall + np.sum(flip_bool)
        flip = flip_bool.astype(np.int8)*2-1
        v = -v*flip


    f+= (2*m*p_zdirection)/time #force on top wall

    return exiting, col_wall, f


exiting, col_wall, f = all_part(x,v)

particles_per_second = exiting/time  #The number of particles exiting per second
mean_force = f/steps                 #The box force averaged over all time steps
box_mass = particles_per_second*m




import numpy as np
from numba import jit
from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem

utils.check_for_newer_version()
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

exiting_count, collision_count, force = all_part(x,v)


home_planet_idx = 0 # The home planet always has index 0
home_planet_mass = mission.system.masses[home_planet_idx]*m_sun
home_planet_radii = mission.system.radii[home_planet_idx]*1000 #km
home_planet_grav = home_planet_mass*G/home_planet_radii**2

"""Tot_box = 100*(G*home_planet_mass)/(mean_force*(home_planet_radii)**2)"""
Tot_box = mission.spacecraft_area/(2*L**2)
Tot_force = (Tot_box)*f
init_rocket_mass = mission.spacecraft_mass


f= open("Parameter.py","w+")
f.write("mean_force = {:g}\n".format(mean_force))
f.write("box_mass = {:g}\n".format(box_mass))
f.write("Tot_box = {:g}\n".format(Tot_box))
f.write("Tot_force = {:g}\n".format(Tot_force))
f.write("particles_per_second = {:g}\n".format(particles_per_second))
f.write("home_planet_mass = {:g}\n".format(home_planet_mass))
f.write("home_planet_radii = {:g}\n".format(home_planet_radii))
f.write("init_rocket_mass = {:g}\n".format(init_rocket_mass))
f.write("home_planet_idx = {}\n".format(home_planet_idx))
