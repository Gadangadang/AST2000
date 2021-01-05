import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts
#import part_3
import test_part6

#Få raketten i lavere orbit
def enter_low_orbit(landing_sequence):
    """
    Here you can implement the commanding of the spacecraft for
    challenge A of Part 6 of the project.
    """
    distanse = []

    rock_pos =np.zeros((1000,2))
    j = 0
    p = 0
    t = np.linspace(0,6*24*60*60,1000) # Run for 1 days
    landing_sequence.look_in_direction_of_planet(planet_idx=1)
    landing_sequence.start_video()
    time, pos_act, vel = landing_sequence.orient()

    #Injeksjonsmanøver
    hyp = -pos_act

    dt_planet = 1E-8
    v_planet = (np.asarray(test_part6.f_planet(1,dt_planet))-\
    np.asarray(test_part6.f_planet(1,dt_planet)))/(2*dt_planet)

    tan_vec = np.asarray([-hyp[1]/np.linalg.norm(hyp),\
    hyp[0]/np.linalg.norm(hyp),0])

    v_stable = np.sqrt(G*system.masses[1]*m_sun/np.linalg.norm(hyp))
    v_stable *= tan_vec
    landing_sequence.boost(-v_stable+v_planet - vel)

    #Loop for å få raketten nærmere planeten
    for i in t:
        landing_sequence.fall_until_time(i)
        time, pos_act, vel = landing_sequence.orient()
        rock_pos[j,0] = pos_act[0]
        rock_pos[j,1] = pos_act[1]


        hyp = pos_act
        r = np.linalg.norm(hyp) - system.radii[1]*1000
        print("Distanse is {}km".format(r/1000))
        distanse.append(r)


        #Send raketten nærmere, med andre ord brems veldig opp, og fall så nærme
        #planeten vi klarer
        if j == 5:
            time, pos_act, vel = landing_sequence.orient()
            p_retning = -pos_act/np.linalg.norm(-pos_act)

            v_retning = np.array([p_retning[1],-p_retning[0],p_retning[2]])
            boost_low = np.array([900,900,0])*(-v_retning)
            landing_sequence.boost(boost_low)

        #Sett inn ny injeksjonsmanøver slik at vi holder banen i punktet som er nærmest
        #[900,900]m/s, j = 199 step for 6500kmoh, [800,800]m/s j = 210 for 11800kmoh
        if j == 199:
            hyp = -pos_act

            dt_planet = 1E-8
            v_planet = (np.asarray(test_part6.f_planet(1,dt_planet))-\
            np.asarray(test_part6.f_planet(1,dt_planet)))/(2*dt_planet)

            tan_vec = np.asarray([-hyp[1]/np.linalg.norm(hyp),\
            hyp[0]/np.linalg.norm(hyp),0])
            v_stable = np.sqrt(G*system.masses[1]*m_sun/np.linalg.norm(hyp))
            v_stable *= tan_vec
            landing_sequence.boost(-v_stable+v_planet - vel)

        #Bilder for eventuell valg av landingsposisjon
        if j == 850:
            landing_sequence.take_picture("pic_{}.xml".format(int(j/10)))
        if j == 880:
            landing_sequence.take_picture("pic_{}.xml".format(int(j/10)))
        if j == 910:
            landing_sequence.take_picture("pic_{}.xml".format(int(j/10)))
        if j == 940:
            landing_sequence.take_picture("pic_{}.xml".format(int(j/10)))
        if j == 970:
            landing_sequence.take_picture("pic_{}.xml".format(int(j/10)))
            """break, her stopper vi fordi vi ønsker senere å starte raketten i
             denne posisjonen, men er kommentert ut for å se hvordan rakettbanen
              utvikler seg over hele tiden vi simulerer"""
        if j == 1000:
            landing_sequence.take_picture("pic_{}.xml".format(int(j/10)))
        if j == 820:
            landing_sequence.take_picture("pic_{}.xml".format(int(j/10)))

        j += 1

    landing_sequence.finish_video(filename="Lower_orbit.xml")

    #Plotting sekvens

    plt.plot(t,distanse)
    plt.xlabel("Time")
    plt.ylabel("Distanse [m]")
    plt.savefig("Distanse.jpeg")
    plt.show()

    fig, ax = plt.subplots()
    circle = plt.Circle((0,0),system.radii[1]*1000, color='peru',label="Planet 1")
    ax.add_artist(circle)

    area = (system.radii[1]*1000)**2
    ax.plot(rock_pos[:,0],rock_pos[:,1],"b--",label="Rocket Path")
    ax.plot(rock_pos[-1,0],rock_pos[-1,1],"ko",label="Rocket")


    plt.legend(loc="lower left")
    plt.xlabel("Distanse [m]")
    plt.ylabel("Distanse [m]")
    plt.axis("equal")

    plt.savefig("Rocket_orbit.jpeg")
    plt.show()


# Prevent the following code from executing when calling `import part_6`
if __name__ == '__main__':

    # Print a message if a newer version of ast2000tools is available
    #utils.check_for_newer_version()
    seed = utils.get_seed('Sgfrette')
    mission = SpaceMission(seed)
    system = SolarSystem(seed)


    # Load the mission instance saved after completing Part 5
    mission = SpaceMission.load('mission_after_part_5.pickle')

    # Initiate the landing sequence
    landing_sequence = mission.begin_landing_sequence()

    # Send the spacecraft into a low circular orbit
    enter_low_orbit(landing_sequence)


    # Perform spectral line fitting to estimate mean_molecular_mass here...


    # Save the mission instance after scouting
    SpaceMission.save('mission_after_part_6.pickle', mission)
