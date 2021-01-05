import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from ast2000tools.shortcuts import SpaceMissionShortcuts

utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)
codes = [8208,389] # The codes recieved from the group teacher (in this case just one)
shortcuts = SpaceMissionShortcuts(mission, codes)

shortcuts.place_spacecraft_on_escape_trajectory(9E-10,40000,pi/3, 20000)
fuel_mass_after_launch, time_after_launch, rocket_position_after_launch,\
 rocket_velocity_after_launch= shortcuts.get_launch_results()
fuel_mass_after_launch=20000

# 1 Formula for radial velocity
def dp_shift(lamda):
    r = c*lamda/mission.reference_wavelength
    return r

theta1,theta2 = mission.star_direction_angles[0],mission.star_direction_angles[1]

lamda1, lamda2 = mission.star_doppler_shifts_at_sun[0],mission.star_doppler_shifts_at_sun[1]

#Radial velocity for star with respect to reference star
r1 = dp_shift(lamda1)
r2 = dp_shift(lamda2)
d_vel = [r1,r2]
d_vel = tuple(d_vel)


#Coordinate change
def unit_convert(phi1,phi2,vel):
    A = np.array([[np.sin(phi2),-np.sin(phi1)],[-np.cos(phi2),np.cos(phi1)]])
    vel_vec = np.array([vel[0],vel[1]])
    new_vel = 1/np.sin(phi2-phi1)*np.dot(A,vel_vec)
    return new_vel

#c.1.3
def dp_shift_3(lamda,measured):
    r = c*(lamda-measured)/mission.reference_wavelength
    return r

me_la1,me_la2 = mission.measure_star_doppler_shifts()[0],\
mission.measure_star_doppler_shifts()[1]

r1 = dp_shift_3(lamda1,me_la1)
r2 = dp_shift_3(lamda2,me_la2)
d_vel = [r1,r2]
d_vel = tuple(d_vel)



#Generalized code

def velocity_rocket(m_lamda1,m_lamda2,phi1,phi2):
    lamda1, lamda2 = mission.star_doppler_shifts_at_sun[0],mission.star_doppler_shifts_at_sun[1]

    def dp_shift_3(lamda,m_lamda):
        r = c_AU_pr_s*yr*(lamda-m_lamda)/mission.reference_wavelength
        return r

    def unit_convert(phi1,phi2,vel):
        A = np.array([[np.sin(phi2),-np.sin(phi1)],[-np.cos(phi2),np.cos(phi1)]])
        vel_vec = np.array([vel[0],vel[1]])
        new_vel = 1/np.sin(phi2-phi1)*np.dot(A,vel_vec)
        return new_vel

    r1 = dp_shift_3(lamda1,m_lamda1)
    r2 = dp_shift_3(lamda2,m_lamda2)
    d_vel = [r1,r2]
    d_vel = tuple(d_vel)
    vx,vy = unit_convert(theta1,theta2,d_vel)[0],unit_convert(theta1,theta2,d_vel)[1]
    return vx,vy

phi1,phi2 = mission.star_direction_angles[0],mission.star_direction_angles[1]
me_la1,me_la2 = mission.measure_star_doppler_shifts()[0],\
mission.measure_star_doppler_shifts()[1]

vx,vy = velocity_rocket(me_la1,me_la2,phi1,phi2)
print((vx,vy))
