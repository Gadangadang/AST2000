import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from ast2000tools.shortcuts import SpaceMissionShortcuts
import Falling
#utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

#Simulering av landing
def sim_landing(landing_sequence,thruster_boost,time_siq_begin,time_of_launch,\
                distance_for_boost, size_of_shute,height_og_shute):

    t,pos_act, vel = landing_sequence.orient()
    dir_planet = -pos_act/np.linalg.norm(pos_act)
    """
    Launch lander with lander vel
    """
    dt_length = 1E-1
    x,v,t = Falling.falling(pos_act,vel*37/40,dt_length,\
    mission.lander_area,mission.lander_mass,distance_for_boost,height_og_shute,thruster_boost)
    time1 = time_siq_begin + t
    ax.scatter(x[-1,0],x[-1,1],label='Time of shute')
    """
    Launch parachute
    """
    #Ny simulering av at raketten "faller".
    dt_length = 0.05
    x1,v,t = Falling.falling(x[-1],np.asarray(v),dt_length,\
    size_of_shute,mission.lander_mass,distance_for_boost,0,thruster_boost)
    x_tot = np.concatenate((x,x1))
    time1 += t
    if np.linalg.norm(v)<3 and np.linalg.norm(x_tot[-1])<10:
        print('You have crash landed!')
    if np.linalg.norm(v)>3 and np.linalg.norm(x_tot[-1])<10:
        print('You have successfully landed!')

    #Finn radiell hastighet i siste tidssteg
    r = x_tot[-1]
    v_r = v -(np.linalg.norm(r)*2*np.pi/(system.rotational_periods[1]*day))\
    * np.asarray([-r[1]/np.linalg.norm(r),r[0]/np.linalg.norm(r),0])
    print(np.dot(v_r,-r/np.linalg.norm(r)))
    return x_tot, time1

if __name__ == '__main__':

    # Print a message if a newer version of ast2000tools is available
    utils.check_for_newer_version()
    seed = utils.get_seed('Sgfrette')
    mission = SpaceMission(seed)
    system = SolarSystem(seed)

    # Load the mission instance saved after completing Part 5
    mission = SpaceMission.load('mission_after_part_6.pickle')
    landing_sequence = mission.ongoing_landing_sequence

    #Sett paramtre og liknende
    t,pos_act, vel = landing_sequence.orient()
    time_of_launch = 0
    thruster_boost = 550
    distance_for_boost = 50
    height_og_shute = 2000
    size_of_shute = 16 #m^2
    time_siq_begin= 0
    fig, ax = plt.subplots()

    #Simuler landing
    x,t = sim_landing(landing_sequence,thruster_boost,time_siq_begin,time_of_launch,\
                    distance_for_boost, size_of_shute, height_og_shute)

    print("t = {}h".format(t/60/60))

    #Plotting
    circle = plt.Circle((0,0),system.radii[1]*1000, color='peru')
    ax.add_artist(circle)
    ax.plot(x[:,0],x[:,1],label='lander')
    ax.plot(x[-1,0],x[-1,1],"ro",label='Landing Finish')
    ax.set(xlim=(-np.linalg.norm(x[0]),np.linalg.norm(x[0])), \
    ylim=(-np.linalg.norm(x[0]),np.linalg.norm(x[0])))
    plt.legend(loc="lower left")
    plt.xlabel("Distanse [m]")
    plt.ylabel("Distanse [m]")

    plt.savefig("Sim_landing_m4.jpeg")
    plt.show()
