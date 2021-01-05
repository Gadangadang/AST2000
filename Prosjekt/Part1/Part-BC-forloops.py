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

#INITIAL CONDITIONS
sigma = np.sqrt((k*T)/m)     #The standard deviation of our particle velocities

x =  np.random.uniform(0,L, size = (int(N), 3))  # Position vector, fill in uniform distribution in all 3 dimensions
v =  np.random.normal(0,sigma, size = (int(N), 3)) # Velocity vector, fill in correct distribution

#SIMULATION VARIABLES
time = 10E-11                       #Simulation Runtime in Seconds
steps = 1000                 #Number of Steps Taken in Simulation
dt = time/steps                  #Simulation Step Length in Seconds
"""
#plot velocity up against probability of particles
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
@jit(nopython=True)
def all_part(x,v):
    exiting = 0         #The total number of particles that have exited the gas box
    f = 0               #Used to calculate Force/Second later in the simulation
    col_wall = 0        #Collisions with wall

    for i in range(int(steps)):
        """
        Each loop must begin by allowing each particle to move depending on its velocity
        """

        x+= dt*v # Change in position
        p_zdirection = 0 #momentum in z direction each time step

        for j in range(int(N)): #Check each particle
            """Check if collision with top area of box,
             then calculate momentum in that direction"""


            for k in range(3):  #Check if x,y,z- components are outside box

                """Check if a parcticle's posistion comp is outside box,
                and check direction of same comp in velocity to correct path"""
                if (x[j,k] <= 0 and v[j,k] < 0) or (x[j,k] >= L and v[j,k] > 0):

                    """Change velocity comp to correct path of particle"""
                    v[j,k] = -v[j,k]
                    col_wall = col_wall + 1.0 #Add for each collision with wall

                    """Check if exited bottom floor of box"""
                    if x[j,2] <= 0:
                        p_zdirection += m*abs(v[j,2])
                        exiting += 1
                        x[j,:] = np.random.uniform(0,L,3)


    f+= (2*p_zdirection)/T #force on top wall

    return exiting, col_wall, f


exiting, collision, f = all_part(x,v)

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
print(exiting)
print(collision)

particles_per_second = exiting/time  #The number of particles exiting per second
mean_force = f/steps                 #The box force averaged over all time steps
box_mass = particles_per_second*m          #The total fuel loss per second

print('There are {:g} particles exiting the gas box per second.'\
.format(particles_per_second))
print('The gas box exerts a thrust of {:g} N.'.format(mean_force))
print('The box has lost a mass of {:g} kg/s.'.format(box_mass))
