import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts
import test
from Del5_B import ship_traj, grav_sun,grav_planet

#utils.check_for_newer_version()
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)
codes = [8208,389] # The codes recieved from the group teacher (in this case just one)
shortcuts = SpaceMissionShortcuts(mission, codes)

tid = 15.75

simulation_time = 20
x = np.array([])

shortcuts.place_spacecraft_on_escape_trajectory(150000,15,tid,20000000,0, 80000)
fuel_mass_after_lanuch, time_after_launch, rocket_position_after_launch,\
 rocket_velocity_after_launch = shortcuts.get_launch_results()
fuel_mass_after_lanuch= 80000

Gr = G_sol
m_ship = (mission.spacecraft_mass + fuel_mass_after_lanuch)/m_sun

x0,v0,t0  = rocket_position_after_launch,\
rocket_velocity_after_launch,time_after_launch


time = np.array([0,2.5]) #simulation_time
#siste delta_v = [0.3,0.2]
delta_v = np.array([[0.70672,0.66795]])
#delta_v = np.array([[0.745,0.772],[0.1,-0.425],[0.4,0.2]])
def fuel_convert(rocket_thrust, rocket_mass_loss_rate, initial_rocket_mass, speed_boost):
    delta_t = (speed_boost*yr/AU*initial_rocket_mass)/rocket_thrust
    fuel_consumed = delta_t*rocket_mass_loss_rate/m_sun
    return fuel_consumed


def boost(x0,v0,t0,t_float,delta_v,mass_lost):
    x1,v1,t1 = ship_traj(t0,x0,\
    v0+delta_v,t_float,(7.6*40+3)/1400000,m_ship-mass_lost)
    return x1,v1,t1

for i in range(len(time)):
    if i ==0:
        mass_lost = fuel_convert(150000,15,m_ship,delta_v[i])
        x2,v0,t0 = boost(x0,v0,tid,time[1],delta_v[i],mass_lost)
    if i == 1:
        hyp = test.f_planet(1,time[i]+tid) - x0
        print(hyp)
        dt_planet = 1E-8
        v_planet = (np.asarray(test.f_planet(1,time[i]+tid+dt_planet))-\
        np.asarray(test.f_planet(1,time[i]+tid-dt_planet)))/(2*dt_planet)
        tan_vec = np.asarray([-hyp[1]/np.linalg.norm(hyp),\
        hyp[0]/np.linalg.norm(hyp)])
        v_stable = np.sqrt(Gr*system.masses[1]/np.linalg.norm(hyp))
        print(v_stable)
        v_stable *= tan_vec
        mass_lost -= fuel_convert(150000,15,m_ship,-v_stable-v0)
        x2,v0,t0  = boost(x0,v0,tid + time[i],simulation_time-time[i],\
        -v_stable+v_planet - v0,mass_lost)
        print(-v_stable-v_planet - v0)
    if i ==0:
        x = x2
    else:
        x = np.concatenate((x,x2))
    x0 = x2[-1]

def interpol(a,b):

    f = interpolate.interp1d(a,b,axis=0)
    return f

times = np.zeros(len(x))
times[0] = tid
dt = (simulation_time)/len(x)
N = len(x)
for i in range(N-1):
    times[i+1] = times[i] + dt

fx = interpol(times,x[:,0])
fy = interpol(times,x[:,1])
v0_with_boost = rocket_velocity_after_launch + delta_v[0]

if __name__ == "__main__":

    l = np.linalg.norm(x[-1]*np.sqrt(system.masses[1]/(10*system.star_mass)))
    distance_from_planet = np.linalg.norm(test.f_planet(1,tid+simulation_time)-x[-1])
    print(distance_from_planet,l)
    print(distance_from_planet<l)
    plt.plot(x[:,0],x[:,1],label="Rocket")
    plt.axis("equal")
    test.plot_rute(tid,simulation_time)
