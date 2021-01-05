import numpy as np
from numba import jit
from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
import Parameter as par
import matplotlib.pyplot as plt



#Prerequisite for program
utils.check_for_newer_version()
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

#Constants
home_planet_idx = par.home_planet_idx # The home planet always has index 0
home_planet_mass = par.home_planet_mass
home_planet_radii = par.home_planet_radii #in km
init_rocket_mass = par.init_rocket_mass
rocket_fuel = 10000 #kg


def home_planet_grav(r):
    grav = -home_planet_mass*G/(home_planet_radii+r)**2
    return grav


mass_fuel_Sec = par.box_mass*par.Tot_box #kg/s
rocket_mass = init_rocket_mass + rocket_fuel
v_esc = np.sqrt((G*2*home_planet_mass)/(home_planet_radii))
Tot_rocket_force = par.Tot_force



t = 0
v = 0
x = 0 #Definerer nullpunkt på planetoverflaten

Max_time = 1200 #sec
steps = 10E6
dt = Max_time/steps

vel = []
pos = []
time = []

while v <= v_esc:

    F = Tot_rocket_force + home_planet_grav(x)

    a = F/(rocket_mass)
    v = v + a*dt
    x = x + v*dt
    t = t + dt
    vel.append(v)
    pos.append(x)
    time.append(t)

    rocket_mass -= mass_fuel_Sec*dt

print(AU)
def convert_to_solar_system_frame_simple(final_height_above_surface,final_upward_speed,launch_duration):

    """
    Here you can implement the coordinate system conversion for
    challenge F of Part 1 of the project.
    """
    tangential_speed = 2*pi/(system.rotational_periods[0]*day)/(AU/yr) #AU/yr
    rocket_position_before_launch = [system.initial_positions[0][0],system.initial_position[1][0]]#AU
    rocket_velocity_after_launch = [final_upward_speed/(AU/yr),system.initial_velocities[1][0]+tangential_speed]#AU/yr
    rocket_position_after_launch = [system.initial_positions[0][0]+final_height_above_surface/AU,rocket_velocity_after_launch[1]*launch_duration/yr]
    time_after_launch = launch_duration/yr

    return rocket_position_before_launch, \
           rocket_position_after_launch,  \
           rocket_velocity_after_launch,  \
           time_after_launch
print(convert_to_solar_system_frame_simple(x,v,t))

"""
print("Time: {}\n Position over planet: {:g}\nVelocity: {:g}\n \
 Rocket mass: {:g}\nIf pos, then esc {}"\
 .format(t,x,v,rocket_mass, v-v_esc))
print("Total fuel consumtion is {} kg".format(mass_fuel_Sec*t))

plt.plot(time,vel,Label="Velocity")

plt.legend()
plt.show()
plt.plot(time,pos,Label="Posisjon above ground")

plt.legend()
plt.show()
"""
