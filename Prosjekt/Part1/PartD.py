import numpy as np
from numba import jit
from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from PartBC import *
from ast2000tools.solar_system import SolarSystem


#Prerequisite for program
utils.check_for_newer_version()
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

#Constants
exiting_count, collision_count, force = all_part(x,v)

steps = 1000
mean_force = force/steps
home_planet_idx = 0 # The home planet always has index 0
home_planet_mass = mission.system.masses[home_planet_idx]*m_sun
home_planet_radii = mission.system.radii[home_planet_idx]
home_planet_grav = home_planet_mass*G/home_planet_radii**2

particles_per_second = exiting/time  #The number of particles exiting per second
box_mass = particles_per_second*m          #The total fuel loss per second
Tot_box = (G*home_planet_mass)/(mean_force*(mission.system.radii[home_planet_idx])**2)
Tot_force = Tot_box*mean_force
init_rocket_mass = mission.spacecraft_mass
speed_boost = np.sqrt((G*home_planet_mass)/home_planet_radii)#Choose whatever in m/s

#Calculate total fuel consumption for given speed boost
def fuel_convert(tot_force, fuel_consumption, initial_rocket_mass, speed_boost):
    delta_t = (speed_boost*initial_rocket_mass)/tot_force
    fuel_consumed = delta_t*fuel_consumption

    return delta_t, fuel_consumed

delta_t, tot_fuel = fuel_convert(Tot_force,box_mass,\
init_rocket_mass,speed_boost )

if __name__ == '__main__':
    """print('The gas box exerts a thrust of {:g} N.'.format(mean_force))
    """
    print("Total amound of boxes is at least {:g}".format(Tot_box))
    print("Total force from rocket is {:g}".format(Tot_box*mean_force))
    print("Total fuel used is {:g} kg during {:g}s".format(delta_t, tot_fuel))
    print("{:g}".format(box_mass*Tot_box))
