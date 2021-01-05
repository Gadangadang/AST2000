import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts

import Part4
"""
Confirms the manual orientation
"""
import Traj_test2
"""
Simulates orbit of rocket
"""
import test
"""
Interpolates planet orbits
"""

def fuel_convert(rocket_thrust, fuel_consumtion, initial_rocket_mass, speed_boost):
    delta_t = (speed_boost*yr/AU*initial_rocket_mass)/rocket_thrust
    fuel_consumed = delta_t*fuel_consumtion/m_sun
    return fuel_consumed

def send_spacecraft(interplanetary_travel,
                    destination_planet_idx,
                    boost_times,fuel_mass_after_launch,
                    fuel_consumtion,thrust):
    # Insert awesome code here..
    mass_ship = (mission.spacecraft_mass + fuel_mass_after_launch)/m_sun
    dt_check = 0.3
    interplanetary_travel.start_video()
    interplanetary_travel.look_in_direction_of_planet(0)

    rocket_travel = np.zeros((len(boost_times),2))
    #First check for position, time and velocity
    time, pos_act, vel = interplanetary_travel.orient()
    i = 0
    p = 0
    for time_check in boost_times: #Checks for every 0.1 yrs if boost is necessary
        if time_check == 0:

            #First boost to send rocket towards destination
            interplanetary_travel.boost(Traj_test2.v0_with_boost-vel)
            mass_ship = mass_ship - fuel_convert(thrust,fuel_consumtion, mass_ship,\
            Traj_test2.v0_with_boost-vel) #Calculate change in mass

        interplanetary_travel.coast_until_time(15.75+time_check) #Coast untill that change in time

        #Takes picture every 0.1 yrs after 3.5 yrs of travel

        if time_check > 3.5:
            interplanetary_travel.take_picture("Pick_number_{}_4.xml".format(i))
            i+=1


        #After 0.55 yrs, change direction of camera towards destination planet
        if time_check > 0.25:
            interplanetary_travel.look_in_direction_of_planet(1)

        #Find orientation for rocket
        time, pos_act, vel = interplanetary_travel.orient()

        rocket_travel[p,0] = pos_act[0]
        rocket_travel[p,1] = pos_act[1]
        #Find simulated values to compare with real travel
        exs_pos_without_boost,v,t = Traj_test2.ship_traj(15.75+time_check,pos_act,\
        vel,dt_check,(7.6*40+3)/1400000,mass_ship)


        pos_sim = np.asarray([Traj_test2.fx(15.75+time_check\
        +dt_check),Traj_test2.fy(15.75+time_check+dt_check)])

        #Find eventual rate of change for velocity
        delta_v = (pos_sim - exs_pos_without_boost[-1])/dt_check

        #Check if boost is necessary for correction, if abs(boost) required to align
        #real and simulated rocket trajectory is below 0.1 AU/yr, we dont bother
        if np.linalg.norm(delta_v)>0.1:
            interplanetary_travel.boost(delta_v)
            mass_ship = mass_ship - fuel_convert(thrust,fuel_consumtion, mass_ship,delta_v)
            time, pos_act, vel = interplanetary_travel.orient()

        p += 1

    #Injection boosts
    hyp = test.f_planet(1,time) - pos_act
    dt_planet = 1E-8

    v_planet = (np.asarray(test.f_planet(1,time+dt_planet))-\
    np.asarray(test.f_planet(1,time-dt_planet)))/(2*dt_planet)

    tan_vec = np.asarray([-hyp[1]/np.linalg.norm(hyp),\
    hyp[0]/np.linalg.norm(hyp)])
    v_stable = np.sqrt(G_sol*system.masses[1]/np.linalg.norm(hyp))
    v_stable *= tan_vec
    interplanetary_travel.boost(-v_stable+v_planet - vel)

    rock_pos =np.zeros((9000,2))
    distanse = []
    j = 0
    t = np.linspace(19.75,20.25,1000)
    t_before = np.linspace(15.75,19.75,1000)
    for i in t:
        interplanetary_travel.coast_until_time(i)
        time, pos_act, vel = interplanetary_travel.orient()
        rock_pos[j,0] = pos_act[0]
        rock_pos[j,1] = pos_act[1]
        j += 1

        hyp = test.f_planet(1,i) - pos_act
        r = np.linalg.norm(hyp)

        if j == 50:
            time, pos_act, vel = interplanetary_travel.orient()
            p_retning = (test.f_planet(1,i)-pos_act)/np.linalg.norm((test.f_planet(1,i)-pos_act))
            v_retning = np.array([p_retning[1],-p_retning[0]])
            boost_low = np.array([0.05,0.05])*(-v_retning)
            interplanetary_travel.boost(boost_low)
            #plt.plot(rock_pos[j,0]-test.fx1(i),rock_pos[j,1]-test.fy1(i),"go",label='First boost')



        if j == 265:
            time, pos_act, vel = interplanetary_travel.orient()
            p_retning = (test.f_planet(1,i)-pos_act)/np.linalg.norm((test.f_planet(1,i)-pos_act))
            v_retning = np.array([p_retning[1],-p_retning[0]])
            boost_low = np.array([0.1,0.1])*(-v_retning)
            interplanetary_travel.boost(boost_low)
            #plt.plot(rock_pos[j,0]-test.fx1(i),rock_pos[j,1]-test.fy1(i),"go",label='Second boost')

        if j >= 265:
            distanse.append(r)


        target = 0.00176189005
        if np.linalg.norm(test.f_planet(1,i)-pos_act) < target:
            print("Distanse = {:g}AU, Inside by : {:g}AU".format(r,r-target))
        else:
            print("Distanse = {:g}AU, Miss by : {:g}AU".format(r,r-target))


    #Find if rocket follows planet for 7 yrs
    t = np.linspace(20.25,24.25,8000)




    for i in t:
        interplanetary_travel.coast_until_time(i)
        time, pos_act, vel = interplanetary_travel.orient()
        rock_pos[j,0] = pos_act[0]
        rock_pos[j,1] = pos_act[1]

        hyp = test.f_planet(1,i) - pos_act
        r = np.linalg.norm(hyp)
        sintheta = hyp[1]/r
        costheta = hyp[0]/r


        j +=1

        print("r: {}".format(r))



    t_ny =  np.linspace(15.75,24.25,9000)
    t = np.linspace(19.75,24.25,9000)

    rocket_orbit = np.array([rock_pos[:,0]-test.fx1(t),rock_pos[:,1]-test.fy1(t)])
    semimajor_axis = np.max(abs(rocket_orbit[0]))
    semiminor_axis = np.max(abs(rocket_orbit[1]))

    print(semiminor_axis/semimajor_axis)
    if semimajor_axis > semiminor_axis:
        e = np.sqrt(1-(semimajor_axis/semiminor_axis)**2)
    else:
        e = np.sqrt(1-(semiminor_axis/semimajor_axis)**2)
    print("a = {}, b = {}, e = {}".format(semimajor_axis,semiminor_axis,e))

    #t = np.concatenate((t,t_float))
    distanse = np.asarray(distanse)
    R = np.sum(distanse)
    snitt_r = R/len(distanse)

    T = 2*np.pi*np.sqrt(snitt_r**3/(G_sol*system.masses[1]))
    print("Average orbit time in lower orbit is {:g} yrs, or {:g}days ".format(T,T*365))
    p

    #Plotting
    plt.plot(test.fx1(t_ny),test.fy1(t_ny),label="Orbit planet 1")
    plt.plot(test.fx0(t_ny),test.fy0(t_ny),label="Orbit planet 0")
    plt.plot(rocket_travel[:,0],rocket_travel[:,1],label='Travel rocket')
    plt.plot(rock_pos[:,0],rock_pos[:,1],label='Orbit stable')
    plt.legend()
    plt.title('Rocket and planet')
    plt.xlabel('AU')
    plt.ylabel('AU')
    plt.savefig('Actual_Launch_orbit.png')
    plt.show()

    plt.plot(rock_pos[:,0]-test.fx1(t),rock_pos[:,1]-test.fy1(t),label='Orbit stable')
    plt.plot(rock_pos[-1,0]-test.fx1(t[-1]),rock_pos[-1,1]-test.fy1(t[-1]),"ro",label='Last position of rocket')
    plt.plot(0,0,"bo")
    plt.legend()
    plt.title('Rocket and planet')
    plt.xlabel('AU')
    plt.ylabel('AU')
    plt.savefig('Rocket_w_boost_focus.jpeg')
    plt.show()



    #End of simulation

    interplanetary_travel.take_picture("Pick_at_end.xml")
    interplanetary_travel.finish_video(filename='travel_video_4_yr.xml', number_of_frames=2000)

    interplanetary_travel.record_destination(destination_planet_idx)


# Prevent the following code from executing when calling `import part_5`
if __name__ == '__main__':
    #utils.check_for_newer_version()
    # Construct SpaceMission instance for my mission
    seed = utils.get_seed('Sgfrette')
    mission = SpaceMission.load("Part4_godkjent.pickle")
    system = SolarSystem(seed)
    codes = [8208,389] # The codes recieved from the group teacher (in this case just one)
    shortcuts = SpaceMissionShortcuts(mission, codes)
    fuel_consumtion = 15
    thrust = 150000
    shortcuts.place_spacecraft_on_escape_trajectory(thrust,fuel_consumtion,\
    15.75,20000000,0, 80000)
    #Part4.launch_and_orient_spacecraft(mission,0,15.75, fuel_consumtion,thrust)



    # Initiate the real interplanetary travel
    interplanetary_travel = mission.begin_interplanetary_travel()

    fuel_mass_after_launch, time_after_launch, rocket_position_after_launch,\
    rocket_velocity_after_launch = shortcuts.get_launch_results()
    fuel_mass_after_launch = 80000
    print(rocket_velocity_after_launch)
    boost_times = np.linspace(0,4,41) #simulation_time

    destination_planet_idx = 1

    # Send the spacecraft on its way to the destination planet
    send_spacecraft(interplanetary_travel,
                    destination_planet_idx,
                    boost_times,fuel_mass_after_launch,
                    fuel_consumtion,thrust)

    # Save the mission instance after entering orbit
    SpaceMission.save('mission_after_part_5.pickle', mission)
