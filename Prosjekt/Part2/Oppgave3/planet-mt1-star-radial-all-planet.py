import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem

utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)
"""
This is the code for c4
"""

sunmass = system.star_mass

sigma = 0.0000126*1/5
N = 1400000
V_rad = []
Vofsystem = 46

Gr = G_sol
def integrator(planets,sun,planets_index):

    def gravity_on_planet(x_sun,y_sun,x_planet,y_planet,index):

        r = np.asarray([x_sun-x_planet, y_sun-y_planet])
        grav = (Gr*sunmass*masses[index])/(np.linalg.norm(r))**2
        ax = grav*r[0] /(np.linalg.norm(r))
        ay = grav*r[1]/(np.linalg.norm(r))

        return np.asarray([ax,ay],float)

    number_of_planets = int(len(planets))
    cm = np.array([0,0])
    masses = [system.masses[i] for i in planets_index]
    x0xsun,x0ysun,v0xsun,v0ysun = sun

    N = 1400000
    time = 7.6*20

    dt = time/N
    t = dt
    x_planets = np.asarray([np.zeros((N,2),float) for i in range(number_of_planets)])
    v_planets = [np.asarray([planets[i,2],planets[i,3]]) for i in range(number_of_planets)]
    for i in range(number_of_planets):
        x_planets[i,0,:] = [planets[i,0],planets[i,1]]
    a_i_planets = np.zeros((number_of_planets,2))
    a_iplus1_planets = np.zeros((number_of_planets,2))

    for i in range(number_of_planets):
        a_i_planets[i] = gravity_on_planet(x0xsun,x0ysun,planets[i,0],planets[i,1],i)

    x_sun = np.zeros((N,2),float)
    x_sun[0,:] = [x0xsun,x0ysun]
    v_sun = np.asarray([v0xsun,v0ysun],float)
    a_i_sun = -sum(a_i_planets)/sunmass

    for i in range(N-1):
        t+=dt
        """Compute centermass for xy coordinates"""
        cm_x = [x_planets[g,i,0]*masses[g] for g in range(number_of_planets)]
        cm_y = [x_planets[g,i,1]*masses[g] for g in range(number_of_planets)]
        cm = [(sum(cm_x)+sunmass*x_sun[i,0])/(sum(masses)+sunmass), \
        (sum(cm_y)+sunmass*x_sun[i,1])/(sum(masses)+sunmass)]

        for j in range(number_of_planets):
            x_planets[j,i+1,:] = x_planets[j,i,:]+(v_planets[j]*dt + 0.5*a_i_planets[j]/masses[j]*dt**2)-cm
        x_sun[i+1,:] = x_sun[i,:]+(v_sun*dt + 0.5*a_i_sun*dt**2)-cm
        for k in range(number_of_planets):
            a_iplus1_planets[k] = gravity_on_planet(x_sun[i+1,0],x_sun[i+1,1],x_planets[k,i+1,0],\
            x_planets[k,i+1,1],k)
        a_iplus1_sun = -sum(a_i_planets)/sunmass
        for l in range(number_of_planets):
            v_planets[l] += 0.5*(a_i_planets[l] + a_iplus1_planets[l])/masses[l]*dt
        a_i_planets = a_iplus1_planets
        v_sun += 0.5*( a_i_sun + a_iplus1_sun )*dt
        a_i_sun = a_iplus1_sun

        #radial velocity is equal to abs of sun velocity vector with y coordinate vector direction
        v_inp = np.linalg.norm(v_sun)*x_sun[i+1,1]/np.linalg.norm(x_sun[i+1])
        V_rad.append(v_inp)
    return x_planets,x_sun,V_rad, dt

"""Input index of planets wanted in simulation"""
planets_index = [0,1,2,3,4,5,6,7]
planets = []

"""Generate list with info for numerical integration"""
for i in planets_index:
    planets.append([system.initial_positions[0,i],system.initial_positions[1,i],\
    system.initial_velocities[0][i],system.initial_velocities[1][i]])

planets= np.asarray(planets)
Sun = [0,0,0,0]

planet_orbit, sun_orbit,V_rad, dt = integrator(planets,Sun,planets_index)

#Writing sequence
V_rad = np.asarray(V_rad) + Vofsystem
v_rad = V_rad[0::500]
time = []

time = np.asarray(time)

plt.plot(time,v_rad)
plt.xlabel("Years")
plt.ylabel("AU/Years")
plt.savefig("radial-velocity-allplanets.jpeg")
plt.show()

"""
write code
"""
"""infile = open("Radial_vel_allplanets_reduced.txt","w")
for i in range(len(v_rad)):
    infile.write("{} {}\n".format(i*dt,v_rad[i]))
    time.append(i)
infile.close()"""
