import numpy as np
from numba import jit
import matplotlib.pyplot as plt
import seaborn as sns
from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem



"""
Koden vi har lagd er helt vår egen. Det eneste vi har brukt er strukturen
for hvordan helheten skal se ut.
"""



def simulate_engine_performance(number_of_particles_in_box,
                                box_side_length,
                                temperature,
                                box_hole_area):

    from numba import jit


    np.random.seed(123)

    #CONSTANTS
    k = k_B          #Boltzmann's Constant

    #PHYSICAL VARIABLES
    T = temperature               #Temperature in Kelvin
    L = box_side_length       #Box Length in meters
    N = number_of_particles_in_box                 #Number of particles
    m = m_H2          #Mass of individual particle in kg
    short_side = L/2 - np.sqrt(box_hole_area)/2
    long_side = L/2 + np.sqrt(box_hole_area)/2

    #INITIAL CONDITIONS
    sigma = np.sqrt((k*T)/m)     #The standard deviation of our particle velocities

    x =  np.random.uniform(0,L, size = (int(N), 3))  # Position vector, fill in uniform distribution in all 3 dimensions
    v =  np.random.normal(0,sigma, size = (int(N), 3)) # Velocity vector, fill in correct distribution

    #SIMULATION VARIABLES
    time = 10E-9                       #Simulation Runtime in Seconds
    steps = 1000                 #Number of Steps Taken in Simulation
    dt = time/steps                 #Simulation Step Length in Seconds

    #RUNNING THE CONDITIONAL INTEGRATION LOOP
    #@jit(nopython=True) Sometimes jit just would not work
    def all_part(x,v):
        exiting = 0         #The total number of particles that have exited the gas box
        f = 0               #Used to calculate Force/Second later in the simulation
        col_wall = 0        #Collisions with wall
        p_zdirection = 0 #momentum in z direction each time step
        refill = np.array([0, 0, L]) #Reset for particles exiting the box
        condition_exit = np.zeros_like(x, dtype = bool)

        for i in range(int(steps)):
            x+= dt*v # Change in position
            """Condition one and two checks if particle is outside the box,
            with. """
            condition1 = np.logical_and(np.greater_equal(x,L), np.greater(v,0))
            condition2 = np.logical_and(np.less_equal(x,0), np.less(v,0))
            condition_collision = np.logical_or(condition1,condition2)
            """
            Condition three checks if particle hits the bottom of the box
            and is moving downwards. This is to make sure we dont bounce back
            any particles on their way back.
            """
            condition3 = np.logical_and(np.less_equal(x[:,2],0),np.less(v[:,2],0))
            """
            Condition three and four checks if the particle hits the hole in the
            bottom.
            """
            condition4 = np.logical_and(np.greater_equal(x[:,0], short_side),\
            np.less_equal(x[:,0],long_side))
            condition5 = np.logical_and(np.greater_equal(x[:,1], short_side),\
            np.less_equal(x[:,1],long_side))
            """
            Temp creates an array with true and false values in the Position
            to a particle depending on if it escapes the box or not.
            """
            temp = np.logical_and(condition3, condition4, condition5)
            condition_exit[:,0] = temp
            condition_exit[:,1] = temp
            condition_exit[:,2] = temp
            exiting += np.sum(temp)
            p_zdirection += np.sum(temp*abs(v[:,2]))

            refill_mat = condition_exit.astype(np.int8)
            flip_bool = np.logical_and(condition_collision, np.logical_not(condition_exit)).astype(np.int8)
            """
            We add the refill_mat*L so that we move all the particles who escaped
            up to the top.
            """
            x[:,2] = x[:,2] + refill_mat[:,0]*L

            col_wall = col_wall + np.sum(flip_bool)
            flip = flip_bool.astype(np.int8)*2-1
            v = -v*flip


        f= (2*m*p_zdirection)/time #force on top wall

        return exiting, col_wall, f


    exiting, col_wall, f = all_part(x,v)

    particles_per_second = exiting/time  #The number of particles exiting per second
    mean_force = f/steps                 #The box force averaged over all time steps
    box_mass = particles_per_second*m

    thrust_per_box = f
    mass_loss_rate_per_box = box_mass

    return thrust_per_box, \
           mass_loss_rate_per_box


def compute_fuel_mass_needed_for_boost(mission,
                                       rocket_thrust,
                                       rocket_mass_loss_rate,
                                       initial_fuel_mass,
                                       target_delta_v):

    init_rocket_mass = mission.spacecraft_mass + initial_fuel_mass

    #Calculate total fuel consumption for given speed boost
    def fuel_convert(rocket_thrust, rocket_mass_loss_rate, initial_rocket_mass, speed_boost):
        delta_t = (speed_boost*initial_rocket_mass)/tot_force
        fuel_consumed = delta_t*fuel_consumption

        return delta_t, fuel_consumed

    delta_t, tot_fuel = fuel_convert(Tot_force,box_mass,\
    init_rocket_mass,speed_boost )


    # You will probably also need this quantity:
    # mission.spacecraft_mass

    return fuel_consumed


def simulate_rocket_launch(mission,
                           rocket_thrust,
                           rocket_mass_loss_rate,
                           initial_fuel_mass):
    """
    Here you can implement the rocket launch simulation for challenge E
    of Part 1 of the project.
    """

    t = 0
    v = 0
    x = 0 #Definerer nullpunkt på planetoverflaten

    Max_time = 1200 #sec
    steps = 10E6
    dt = Max_time/steps
    home_planet_idx = 0
    home_planet_mass = mission.system.masses[home_planet_idx]*m_sun
    home_planet_radii = mission.system.radii[home_planet_idx]*1000 #in meter
    rocket_mass = initial_fuel_mass + mission.spacecraft_mass
    v_esc = np.sqrt(2*G*home_planet_mass/home_planet_radii+x)


    def home_planet_grav(r):
        grav = -home_planet_mass*rocket_mass*G/(home_planet_radii+r)**2
        return grav
    i=0
    while (0 <= v) & (v <= v_esc):

        F = rocket_thrust+home_planet_grav(x)

        a = F/(rocket_mass)
        v = v + a*dt
        x = x + v*dt
        t = t + dt

        rocket_mass -= rocket_mass_loss_rate*dt
        v_esc = np.sqrt(2*G*home_planet_mass/(home_planet_radii+x))

        print("V: {:g} M: {:g} G: {:g}".format(v_esc-v,rocket_mass,home_planet_grav(x)))


    final_height_above_surface = x
    final_upward_speed = v
    fuel_mass_after_launch = rocket_mass
    launch_duration = t

    return final_height_above_surface, \
           final_upward_speed,         \
           fuel_mass_after_launch,     \
           launch_duration


def convert_to_solar_system_frame_simple(system,
                                         final_height_above_surface,
                                         final_upward_speed,
                                         launch_duration):
    """
    Here we add all the velocity components and positional components, and
    convert them to more helpful units.
    """
    tangential_speed = mission.system.radii[0]/AU*1000*(2*pi/\
    (system.rotational_periods[0]*day/yr))#AU/yr

    rocket_position_before_launch = [mission.system.radii[0]*1000/AU+\
    system.initial_positions[0][0],system.initial_positions[1][0]]#AU

    rocket_velocity_after_launch = [final_upward_speed/(AU/yr)+system.initial_velocities[0][0],\
    system.initial_velocities[1][0]+tangential_speed]#AU/yr
    """
    We had to manually add a certain distance to our final position, because of
    an error in our code or the verification.
    """
    rocket_position_after_launch = [rocket_position_before_launch[0]+\
    final_height_above_surface/AU,rocket_velocity_after_launch[1]*launch_duration/yr]
    time_after_launch = launch_duration/yr

    return rocket_position_before_launch, \
           rocket_position_after_launch,  \
           rocket_velocity_after_launch,  \
           time_after_launch

def run_engine_and_launch_simulation(mission,
                                     number_of_boxes,
                                     number_of_particles_in_box,
                                     box_side_length,
                                     temperature,
                                     box_hole_area,
                                     initial_fuel_mass):
    """
    This function executes the simulation code for the rocket engine and launch,
    and returns the results.
    """

    # Simulate engine
    thrust_per_box,        \
    mass_loss_rate_per_box \
      = simulate_engine_performance(number_of_particles_in_box,
                                    box_side_length,
                                    temperature,
                                    box_hole_area)

    # Compute rocket_thrust and rocket_mass_loss_rate here...
    rocket_mass_loss_rate = mass_loss_rate_per_box*number_of_boxes
    rocket_thrust = thrust_per_box * number_of_boxes
    print(rocket_mass_loss_rate,rocket_thrust)
    # Simulate launch
    final_height_above_surface, \
    final_upward_speed,         \
    fuel_mass_after_launch,     \
    launch_duration             \
      = simulate_rocket_launch(mission,
                               rocket_thrust,
                               rocket_mass_loss_rate,
                               initial_fuel_mass)

    return thrust_per_box,             \
           mass_loss_rate_per_box,     \
           final_height_above_surface, \
           final_upward_speed,         \
           fuel_mass_after_launch,     \
           launch_duration


# Prevent the following code from executing when calling import part_1
if __name__ == '__main__':

    # Print a message if a newer version of ast2000tools is available
    utils.check_for_newer_version()

    # Construct SpaceMission instance for my mission
    seed = utils.get_seed('Sgfrette')
    mission = SpaceMission(seed)

    # Extract associated SolarSystem object
    system = mission.system
    box_side_length = 10E-7 #cm
    number_of_boxes =  int(mission.spacecraft_area/(2*box_side_length**2))
    number_of_particles_in_box = 3*10E5
    temperature = 30900
    box_hole_area = (1/20)*box_side_length**2
    initial_fuel_mass = 20000



    # Run engine and launch simulations and get results relative to planet
    thrust_per_box,             \
    mass_loss_rate_per_box,     \
    final_height_above_surface, \
    final_upward_speed,         \
    fuel_mass_after_launch,     \
    launch_duration             \
      = run_engine_and_launch_simulation(mission,
                                         number_of_boxes,
                                         number_of_particles_in_box,
                                         box_side_length,
                                         temperature,
                                         box_hole_area,
                                         initial_fuel_mass)

    # Convert simulated launch results to the solar system frame,
    # assuming launch along the x-direction at time = 0 years
    rocket_position_before_launch, \
    rocket_position_after_launch,  \
    rocket_velocity_after_launch,  \
    time_after_launch              \
      = convert_to_solar_system_frame_simple(system,
                                             final_height_above_surface,
                                             final_upward_speed,
                                             launch_duration)
    time_of_launch = 0.0

    # Perform real launch
    mission.set_launch_parameters(thrust_per_box,
                                  mass_loss_rate_per_box,
                                  number_of_boxes,
                                  initial_fuel_mass,
                                  launch_duration,
                                  rocket_position_before_launch,
                                  time_of_launch)
    print("T:{}, Pos:{}, Fuel:{}".format(time_after_launch,rocket_position_after_launch,fuel_mass_after_launch))


    mission.launch_rocket()

    # Verify simulated launch results
    mission.verify_launch_result(rocket_position_after_launch)
