import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts
import atmos_part7



def land_spacecraft(landing_sequence,
                    initial_time,
                    duration,
                    time_of_lander_launch,
                    lander_launch_delta_v,
                    parachute_area,
                    parachute_deployment_height,
                    thruster_force,
                    thruster_activation_height):
    """
    Here you can implement the commanding of the spacecraft for
    challenge C of Part 7 of the project.
    """
    landing_sequence.adjust_parachute_area(parachute_area)          # Set the parachute area to 42 m^2
    landing_sequence.adjust_landing_thruster(force=thruster_force) # Set the thruster force to 500 Newtons



     # Generate a still picture
    landing_sequence.start_video()
    landing_sequence.look_in_direction_of_motion()
    t = np.linspace(initial_time,initial_time + duration,50001)
    j = 0
    p = True
    s = True
    f = True
    go = True

    #Simuler faktisk landing
    for i in t:
        landing_sequence.fall_until_time(i)
        time, position, vel = landing_sequence.orient()
        distanse = np.linalg.norm(position)-system.radii[1]*1000
        vr = np.dot(vel,-position/np.linalg.norm(position))
        print("x = {:.5e} moh vr = {:.5f} m/s" .format(distanse,vr))

        #Launch Lander
        if j == 100:
            landing_sequence.launch_lander(lander_launch_delta_v)
            print(np.linalg.norm(position)-system.radii[1]*1000,"meters above ground")

        #Deploy parachute
        if np.linalg.norm(position) <= parachute_deployment_height and s:
            landing_sequence.deploy_parachute()
            s = False

        #Activate thrusters
        if np.linalg.norm(position) <= thruster_activation_height and f:
            landing_sequence.activate_landing_thruster()
            f = False

        #Hvis landa, stopp loop
        if np.linalg.norm(position) < (system.radii[1]*1000 +0.0001):
            time_of_stop = i
            break

        #Endre retning på kamera
        if np.linalg.norm(position) < (system.radii[1]*1000 + 4000) and go:
            landing_sequence.look_in_fixed_direction(polar_angle=0, azimuth_angle=pi)
            go = False

        j += 1

    landing_sequence.take_picture(filename='extreme_case2.xml')
    landing_sequence.finish_video(filename="extreme_case2_vid.xml",radial_camera_offset=0)






if __name__ == '__main__':

    # Print a message if a newer version of ast2000tools is available
    utils.check_for_newer_version()
    seed = utils.get_seed('Sgfrette')
    mission = SpaceMission(seed)
    system = SolarSystem(seed)

    # Load the mission instance saved after completing Part 6
    mission = SpaceMission.load('mission_after_part_6.pickle')

    # Obtain the ongoing landing sequence
    landing_sequence = mission.ongoing_landing_sequence
    SpaceMission.save("mission_before_landing.pickle",mission)
    # Obtain the initial conditions for the landing simulation
    initial_time, initial_position, initial_velocity = landing_sequence.orient()

    # Run planning simulations to determine landing parameters here...
    r_low = 2500 + system.radii[1]*1000#m, 2500m is a rough estimate for a good altitude
    vr = 3 # m/s
    C_d = 1

    #Kan bruke likning, men satt til å være kostant lik 16m^2
    parachute_area = G*(system.masses[1]*m_sun*mission.lander_mass)/r_low**2*\
    2/(C_d*atmos_part7.density_func(r_low)*vr**2)
    print("Parachute size is {}m^2".format(parachute_area))
    duration = 5*60*60
    time_of_lander_launch = initial_time
    lander_launch_delta_v = np.array([-3*initial_velocity[0]/40,-3*initial_velocity[1]/40,0])
    parachute_deployment_height = system.radii[1]*1000 + 1000
    thruster_force = 550
    thruster_activation_height = system.radii[1]*1000 + 50


    # Land the lander safely on the surface of the planet
    land_spacecraft(landing_sequence,
                    initial_time,
                    duration,
                    time_of_lander_launch,
                    lander_launch_delta_v,
                    parachute_area,
                    parachute_deployment_height,
                    thruster_force,
                    thruster_activation_height)
