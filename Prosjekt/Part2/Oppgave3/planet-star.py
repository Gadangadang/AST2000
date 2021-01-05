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

sunmass = system.star_mass
planetmass = system.masses[0]

"""
This code is for c1
"""

Gr = G_sol
def integrator(planet,sun):

    def gravity_on_sun(x_sun,y_sun,x_planet,y_planet):

        r = np.asarray([x_sun-x_planet, y_sun-y_planet])
        grav = -(Gr*sunmass*planetmass)/(np.linalg.norm(r))**2
        ax = grav*r[0] /(np.linalg.norm(r))
        ay = grav*r[1]/(np.linalg.norm(r))

        return np.asarray([ax,ay],float)

    x0x,x0y,v0x,v0y = planet
    x0xsun,x0ysun,v0xsun,v0ysun = sun

    N = 76000#choosen for 10 000 per year
    time = 7.6 #time for 20 rotations for our home planet
    #print(time)
    dt = time/N
    t = dt

    x_planet = np.zeros((N,2),float)
    x_sun = np.zeros((N,2),float)

    cm = np.array([0,0])

    x_planet[0,:] = [x0x,x0y]
    x_sun[0,:] = [x0xsun,x0ysun]

    v_planet = np.asarray([v0x,v0y],float)
    v_sun = np.asarray([v0xsun,v0ysun],float)

    F_sun0 = gravity_on_sun(x0xsun,x0ysun,x0x,x0y)
    a_i_sun = F_sun0/sunmass
    a_i_planet = -F_sun0/planetmass

    for i in range(N-1):
        t+=dt
        cm = np.asarray([(sunmass*x_sun[i,0] + planetmass*\
        x_planet[i,0])/(sunmass+planetmass),\
        (sunmass*x_sun[i,1] + planetmass*\
        x_planet[i,1])/(sunmass+planetmass)])

        x_planet[i+1,:] = x_planet[i,:]+(v_planet*dt + 0.5*a_i_planet*dt**2)-cm
        x_sun[i+1,:] = x_sun[i,:]+(v_sun*dt + 0.5*a_i_sun*dt**2)-cm

        F_sun = gravity_on_sun(x_sun[i+1,0],x_sun[i+1,1],x_planet[i+1,0],\
        x_planet[i+1,1])
        a_iplus1_sun = F_sun/sunmass
        a_iplus1_pl = -F_sun/planetmass
        v_planet += 0.5*( a_i_planet + a_iplus1_pl )*dt
        a_i_planet = a_iplus1_pl
        v_sun += 0.5*( a_i_sun + a_iplus1_sun )*dt
        a_i_sun = a_iplus1_sun
        """
        For testing E.
        """
        """if i in var:
            E_1 = 0.5*my*(np.linalg.norm(v_sun-v_planet))**2- Gr*masstot*my/(np.linalg.norm(x_planet-x_sun))**2
            print(E_1)"""


    return x_planet,x_sun


Planet0 = [system.initial_positions[0,0],system.initial_positions[1,0],\
system.initial_velocities[0][0],system.initial_velocities[1][0]]

Sun = [0,0,0,0]

planet_orbit, sun_orbit = integrator(Planet0,Sun)

plt.plot(planet_orbit[:,0], planet_orbit[:,1], label="Planet 0")
plt.plot(sun_orbit[:,0], sun_orbit[:,1], label="Sun")

plt.xlabel("x-position (AU)")
plt.ylabel("y-position (AU)")
plt.axis("equal")
plt.legend()
plt.savefig("Sun-planet.jpeg")

plt.show()
