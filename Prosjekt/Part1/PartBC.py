import numpy as np
from numba import jit
import matplotlib.pyplot as plt
import seaborn as sns

np.random.seed(123)

#CONSTANTS
k = 1.38064852e-23          #Boltzmann's Constant

#PHYSICAL VARIABLES
T = 3E3               #Temperature in Kelvin
L = 10E-6         #Box Length in meters
N = 10E5                 #Number of particles
m = 3.3474472e-27           #Mass of individual particle in kg
short_side = L/2 - np.sqrt((1/8)*L**2)
long_side = L - np.sqrt((1/8)*L**2)
#INITIAL CONDITIONS
sigma = np.sqrt((k*T)/m)     #The standard deviation of our particle velocities

x =  np.random.uniform(0,L, size = (int(N), 3))  # Position vector, fill in uniform distribution in all 3 dimensions
v =  np.random.normal(0,sigma, size = (int(N), 3)) # Velocity vector, fill in correct distribution

#SIMULATION VARIABLES
time = 10E-11                       #Simulation Runtime in Seconds
steps = 100                 #Number of Steps Taken in Simulation
dt = time/steps                 #Simulation Step Length in Seconds
"""
sns.set(color_codes=True)
sns.distplot(v[:,0], bins=1000)
plt.xlabel("Velocity")
plt.ylabel("Percent of particles")
plt.show()
sns.set(color_codes=True)
sns.distplot([np.linalg.norm(v[i]) for i in range(int(N))], bins=1000)
plt.xlabel("a")
plt.ylabel("b")
plt.show()
"""
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


    f+= (2*m*p_zdirection*exiting)/dt #force on top wall

    return exiting, col_wall, f


exiting, collision, f = all_part(x,v)



if __name__=="__main__":
    particles_per_second = exiting/time  #The number of particles exiting per second
    mean_force = f/steps                 #The box force averaged over all time steps
    box_mass = particles_per_second*m          #The total fuel loss per second

    print('There are {:g} particles exiting the gas box per second.'\
    .format(particles_per_second))
    print('The gas box exerts a thrust of {:g} N.'.format(mean_force))
    print('The box has lost a mass of {:g} kg/s.'.format(box_mass))

    """
    sns.set(color_codes=True)
    sns.distplot(v[:,0], bins=1000)
    plt.xlabel("Velocity")
    plt.ylabel("Percent of particles")
    plt.show()
    sns.set(color_codes=True)
    sns.distplot([np.linalg.norm(v[i]) for i in range(int(N))], bins=1000)
    plt.xlabel("a")
    plt.ylabel("b")
    plt.show()
    """
